using Dates, Printf
import TOML

# ── Shared live state ──────────────────────────────────────────────────────────
mutable struct LiveState
    rt         :: Any
    input_dict :: Any
    last_nc    :: String
    output_dir :: String
end
const LIVE = LiveState(nothing, nothing, "", "")

# ── Plain-variable model state (no Observables) ────────────────────────────────
const STATE = Dict{String,Any}(
    "running"  => false,
    "step"     => 0,
    "time"     => "--:--:--",
    "progress" => 0.0,
    "nactive"  => 0,
)

const LOG_LINES = String[]
const LOG_LOCK  = ReentrantLock()

function log!(msg::String)
    lock(LOG_LOCK) do
        push!(LOG_LINES, msg)
        length(LOG_LINES) > 300 && popfirst!(LOG_LINES)
    end
    @info msg
end

get_log_snapshot() = lock(LOG_LOCK) do; copy(LOG_LINES); end

# ── Control channel ────────────────────────────────────────────────────────────
const CONTROL_CH = Channel{Symbol}(8)

pause_model!()  = STATE["running"] && put!(CONTROL_CH, :pause)
resume_model!() = put!(CONTROL_CH, :resume)
function stop_model!()
    put!(CONTROL_CH, :stop)
    STATE["running"] = false
end

# ── Schema-driven configuration (new path) ────────────────────────────────────
# The browser is the source of truth for the RunSpec as a TOML-shaped dict;
# the server stores the latest POST here and run_model_spec! consumes it.
const SPEC_DICT = Ref{Any}(nothing)

# Default config shown in the GUI before the user loads/edits anything: the
# bundled era5_conv.toml if present, else a bare RunSpec scaffold.  TOML only —
# no legacy FLEXPART option files (see CLAUDE.md hard rule).
function default_spec_dict()
    for cfg in ("era5_conv.toml", "example_dust.toml")
        p = joinpath(@__DIR__, "..", "configs", cfg)
        if isfile(p)
            try
                return runspec_to_dict(read_runspec(p))
            catch e
                @warn "default_spec_dict: could not load $cfg" exception=e
            end
        end
    end
    return runspec_to_dict(RunSpec())
end

# Shared stepping loop (used by both the legacy and the RunSpec paths).
function _run_loop!(rt, kw::NamedTuple, outdir::AbstractString)
    ctx = rt.ctx; st = rt.state; rc = ctx.runconfig
    bench = get(kw, :bench, nothing)
    dir = Int(rc.direction); ts = rc.timestep
    t0   = dir == -1 ? rc.enddate   : rc.startdate
    tend = dir == -1 ? rc.startdate : rc.enddate
    nsteps = round(Int, abs(Dates.value(tend - t0)) / abs(Dates.value(ts)))
    lstep  = max(1, round(Int, abs(Dates.value(rc.outstep)) / abs(Dates.value(ts))))
    log!("nsteps=$nsteps, output every $lstep steps")
    read_met_step!(rt, t0); log!("Met loaded. Running…")
    step_nc = 0; itime = 0; paused = false
    while itime < nsteps
        while isready(CONTROL_CH)
            cmd = take!(CONTROL_CH)
            cmd == :stop && (log!("Stopped at step $itime."); @goto done)
            cmd == :pause && (paused = true; log!("Paused at step $itime."); break)
        end
        while paused
            if isready(CONTROL_CH)
                cmd = take!(CONTROL_CH)
                cmd == :resume && (paused = false; log!("Resumed at step $itime."))
                cmd == :stop && (log!("Stopped while paused."); @goto done)
            else sleep(0.05) end
        end
        itime += 1; cur = t0 + itime * ts
        EASLEE.timeit!(bench, :metread) do; read_met_step!(rt, cur); end
        do_out = (itime % lstep == 0 || itime == nsteps)
        do_out && EASLEE.Types.reset_accumulators!(ctx.output)
        EASLEE.step_model!(st, ctx, cur; do_output=do_out, kw...)
        State.sync_dtime_ns!(st)
        if do_out
            step_nc += 1
            nc_path = joinpath(outdir, @sprintf("easlee_step%05d.nc", step_nc))
            initnetcdf!(ctx.output, nc_path); writenetcdf(ctx.output, nc_path, 1)
            LIVE.last_nc = nc_path
            STATE["nactive"]  = count(==(0x01), st.status)
            STATE["time"]     = Dates.format(cur, "yyyy-mm-dd HH:MM:SS")
            STATE["progress"] = Float64(itime) / Float64(nsteps)
            STATE["step"]     = step_nc
            try
                load_column_from_nc!(nc_path); load_crosssection_from_nc!(nc_path); load_particles_from_live!()
            catch e
                log!("figure update failed: $e")
            end
            log!(@sprintf("step %d/%d  t=%s  active=%d", itime, nsteps,
                 Dates.format(cur, "HH:MM"), STATE["nactive"]))
        end
    end
    log!("Run complete. $(step_nc) output files in $outdir")
    @label done
    if bench !== nothing
        io = IOBuffer(); EASLEE.bench_report(bench; io=io)
        for ln in split(String(take!(io)), '\n'); isempty(ln) || log!(ln); end
    end
    return nothing
end

# Run from a RunSpec dict (schema form / TOML).  Assembles the physics objects
# via build_run, then steps through the shared loop.
function run_model_spec!(cfg::Dict)
    STATE["running"] = true
    task = lowercase(string(get(cfg, "task", "transport")))
    try
        if task == "validation"
            # Validation runs go through the SAME entry point as the CLI
            # (scripts/easlee.jl), so the GUI and command line can never drift.
            _run_task_via_easlee!(cfg, task)
        else
            log!("Building runtime from RunSpec…")
            spec = validate_runspec(dict_to_runspec(cfg))
            log!("Run \"$(spec.name)\": $(spec.direction)  $(spec.start) → $(spec.stop)")
            rt, obj = build_run(spec; prefetch=false)
            LIVE.rt = rt
            @eval Main rt = $rt; @eval Main ctx = $(rt.ctx); @eval Main state = $(rt.state)
            outdir = spec.output.path; mkpath(outdir); LIVE.output_dir = rstrip(outdir, '/')
            _run_loop!(rt, step_kwargs(obj), outdir)
        end
    catch e
        log!("ERROR: $e")
        for l in first(stacktrace(catch_backtrace()), 8); log!("  $l"); end
    end
    STATE["running"] = false
end

# Dispatch a non-transport task to scripts/easlee.jl as a subprocess, streaming its
# stdout/stderr into the GUI log.  The form/loaded config is written to a temp TOML.
function _run_task_via_easlee!(cfg::Dict, task::String)
    repo   = abspath(joinpath(@__DIR__, ".."))
    tmp    = joinpath(tempdir(), "easlee_gui_$(getpid())_$(task).toml")
    open(tmp, "w") do io; TOML.print(io, cfg); end
    easlee = joinpath(repo, "scripts", "easlee.jl")
    log!("task=$task → dispatching via scripts/easlee.jl  ($(basename(tmp)))")
    cmd = `$(Base.julia_cmd()) --project=$(joinpath(repo,"easlee")) $easlee $tmp`
    open(pipeline(cmd; stderr=stdout), "r") do io
        for ln in eachline(io); log!(rstrip(ln)); end
    end
    # Subprocesses write NetCDF but can't touch our FIG_CACHE — load the newest
    # output file so the Figures tab works for validation runs too.
    outdir = get(get(cfg, "output", Dict()), "path", "")
    if isdir(outdir)
        ncs = filter(f -> endswith(f, ".nc"), readdir(outdir; join=true))
        if !isempty(ncs)
            latest = ncs[argmax(mtime.(ncs))]; LIVE.last_nc = latest
            try
                load_column_from_nc!(latest); load_crosssection_from_nc!(latest)
                log!("figures loaded from $(basename(latest))")
            catch e
                log!("figure load failed: $e")
            end
        end
    end
    log!("task=$task finished.")
end

