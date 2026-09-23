#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# easlee.jl — THE single entry point.  Everything is driven by ONE TOML config:
#
#     julia --project=easlee scripts/easlee.jl  <config.toml>
#
# The config's top-level `task` field selects what EASLEE does; the user never
# edits a script.  Tasks:
#
#   task = "transport"   (default)  plain LPDM run  → grid_<name>.nc
#                                    [run]/[meteorology]/[[release]]/[output]/[physics]
#
#   task = "validation"              DATEM tracer-experiment validation (CAPTEX/ETEX):
#                                    run releases, sample stations, ATMES-II stats.
#                                    Driven by [validation] (obs_dir/mode/direction).
#
# Optional [run].threads sets the worker thread count (default 8).
#
# This launcher reads only the `task` (cheap TOML parse) and dispatches to the
# matching handler as a subprocess, passing the SAME config through.  No ENV vars,
# no per-task scripts for the user to touch.
# ─────────────────────────────────────────────────────────────────────────────
using TOML

const REPO = abspath(joinpath(@__DIR__, ".."))

function usage()
    println("""
    usage:  julia --project=easlee scripts/easlee.jl <config.toml>
            the config's  task = "transport"|"validation"  selects the run.
    """)
    exit(1)
end

isempty(ARGS) && usage()
cfg = abspath(ARGS[1])
isfile(cfg) || (println("config not found: $cfg"); exit(1))

d        = TOML.parsefile(cfg)
task     = lowercase(get(d, "task", "transport"))
nthreads = string(get(get(d, "run", Dict()), "threads", 8))
jl       = `$(Base.julia_cmd()) --project=$(joinpath(REPO,"easlee"))`
script(s) = joinpath(REPO, "scripts", s)

println("── EASLEE ──  task=$task  config=$(basename(cfg))  threads=$nthreads")

if task == "transport"
    # The engine: read_runspec(TOML) → build_run → step loop → grid_<name>.nc
    run(`$jl -t $nthreads $(script("run_config.jl")) $cfg`)

elseif task == "validation"
    # Map the TOML ([validation] + standard blocks) → the knobs validate_captex.jl
    # reads.  The user only edits the TOML; the env vars are an internal detail.
    v    = get(d, "validation", Dict())
    run0 = get(d, "run", Dict())
    met0 = get(get(d, "meteorology", [Dict()]), 1, Dict())
    rel0 = get(get(d, "release", [Dict()]), 1, Dict())
    phys = get(d, "physics", Dict())
    turb = get(get(phys, "turbulence", Dict()), "scheme", "hanna")
    conv = get(get(phys, "convection", Dict()), "scheme", "massflux")
    relz = get(v, "relz", [0.0, 20.0])
    env = Dict{String,String}(
        "CAPTEX_OBS" => string(get(v, "obs_dir", "")),
        "CAPTEX_MET" => string(get(met0, "path", "")),
        "MODE"       => string(get(v, "mode", "run")),
        "DIRECTION"  => string(get(v, "direction", "forward")),
        "RELEASE"    => string(get(v, "release", "all")),
        "MAXBACK"    => string(get(v, "maxback", 8)),
        "RUNDAYS"    => string(get(v, "rundays", 4)),
        "DXOUT"      => string(get(v, "dxout", 0.25)),
        "TS"         => string(get(run0, "timestep", "15min")),
        "TURB"       => string(turb),
        "CONV"       => string(conv),
        "NP"         => string(get(rel0, "particles", 200000)),
        "CTL"        => string(get(v, "ctl", 5.0)),
        "IFINE"      => string(get(v, "ifine", 4)),
        "RELZ0"      => string(relz[1]),
        "RELZ1"      => string(relz[end]),
    )
    withenv(env...) do
        run(`$jl -t $nthreads $(joinpath(REPO, "validation", "validate_captex.jl"))`)
    end

else
    println("unknown task=\"$task\" — use transport | validation")
    exit(1)
end
