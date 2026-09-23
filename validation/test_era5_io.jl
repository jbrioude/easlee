# test_era5_io.jl ── isolate and time the ERA5 met I/O.
#
# Builds the Pacific ERA5 runtime and times: (1) the bootstrap (first 2-file load),
# (2) each per-step read_met_step! — only the hourly *rotations* actually load a file,
# so you see the true cost per global-ERA5 file, and (3) optionally advection, so you
# can observe how the async prefetch interacts with the threaded advance!.
#
# Toggles (environment variables):
#   PREFETCH=0|1   async met prefetch on/off            (default 0 = synchronous)
#   NSTEPS=N       number of 15-min steps to march       (default 8 ≈ 2 h ≈ 2 rotations)
#   NP=N           if >0, place N parcels and time advance each step (default 0)
#   TURB=none|hanna  turbulence for the advance timing    (default none)
#   METDIR=path    override the met directory             (default from pathnames file)
#
# Examples:
#   julia -t 6   --project=easlee scripts/test_era5_io.jl                 # sync I/O, no advance
#   PREFETCH=1 julia -t 6,1 --project=easlee scripts/test_era5_io.jl      # async prefetch on :interactive
#   NP=40000 julia -t 6 --project=easlee scripts/test_era5_io.jl          # + time advection each step
#
# Where the time goes (to instrument further): read_met_step! → readmetfiles[_prefetch]!
#   → _load_met_file! → loadgribmetfiles! (eccodes GRIB decode)
#   → metgribfieldcalculation! (derived fields: rho, theta, levdot, zlev, drhodz …)
# Both run over the FULL global grid (Ni×Nj×Nk); cropping to the domain is the main lever.

using Dates, Printf, Random
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl"))
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))   # TOML → runtime (no legacy files)
using .EASLEE, .EASLEE.State, .EASLEE.MetAdapter, .EASLEE.Operators
const ll2idx = Operators.latlontoindyindx

const PREFETCH = get(ENV, "PREFETCH", "0") == "1"
const NSTEPS   = parse(Int, get(ENV, "NSTEPS", "8"))
const NP       = parse(Int, get(ENV, "NP", "0"))
const TURB     = get(ENV, "TURB", "none")
RNG = Random.default_rng()

@printf("=== ERA5 I/O test ===  threads=%d (interactive=%d)  PREFETCH=%s  NSTEPS=%d  NP=%d  TURB=%s\n",
        Threads.nthreads(), Threads.nthreads(:interactive), PREFETCH, NSTEPS, NP, TURB)
flush(stdout)

# ── Build runtime from TOML (the bootstrap @time includes the first 2-file met load) ─
const CFG = get(ENV, "CFG", joinpath(REPO, "configs", "era5_pacific.toml"))
spec = read_runspec(CFG)
haskey(ENV, "METDIR") && (spec.meteorology[1].path = ENV["METDIR"])
for r in spec.release; r.particles = max(NP, 1); end

println("\n-- build_run (incl. first 2 met files) --"); flush(stdout)
rt = nothing
@time (rt, _obj) = build_run(spec; extra_slots=max(NP,1), prefetch=PREFETCH)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
gc  = ctx.nests[1].gridconfig
ts  = rc.timestep; loopstart = rc.startdate
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = TURB == "hanna" ? EASLEE.Types.Hanna(false) : EASLEE.Types.NoTurb()
    try; ctx.nests[n].nestrunconfig.ctl = 2.0f0; catch; end
end

# ── Report the met grid (why each file is heavy) ─────────────────────────────────
let m = ctx.nests[1].metfields3D
    nb = sizeof(m.Ux)
    @printf("\nmet grid: Ni=%d Nj=%d Nk=%d   Ux=%s %s  (%.0f MB/field)\n",
            Int(gc.Ni), Int(gc.Nj), Int(ctx.nests[1].levelconfig.Nk),
            typeof(m.Ux), size(m.Ux), nb/1e6)
end

# ── Optionally place parcels so we can time advance each step ────────────────────
if NP > 0
    fill!(state.status, 0x00)
    for s in 1:NP
        lon = 150.0 + rand(RNG)*80; lat = rand(RNG)*50
        indy, indx = ll2idx(gc, lat, lon)
        state.status[s]=0x01; state.nest[s]=Int8(1); state.relbox[s]=Int32(1)
        state.indx[s]=Float32(indx); state.indy[s]=Float32(indy); state.indl[s]=Float32(2+rand(RNG)*30)
        state.Dtime[s]=loopstart; state.Dtimeinit[s]=loopstart; state.Dtime_ns[s]=Int64(Dates.value(loopstart))
        state.is_new[s]=false; state.dt[s]=1f0; state.mass[s]=1f0; state.rng_state[s]=UInt64(0)
    end
end

# ── March the clock; time read_met_step! (loads only on hourly rotation) ─────────
function run!(rt, ctx, state, nsteps, loopstart, ts, NP)
    println("\nstep   currentdate   read_met[s]   rotated?   advance[s]"); flush(stdout)
    total_io = 0.0; nrot = 0
    for it in 1:nsteps
        cd = loopstart + (it-1)*ts
        t_io = @elapsed read_met_step!(rt, cd)
        rotated = t_io > 0.5                       # a real file load takes seconds
        total_io += t_io; rotated && (nrot += 1)
        t_adv = 0.0
        if NP > 0
            indt = MetAdapter.met_time_index_array(cd, ctx.nests)
            t_adv = @elapsed Operators.advance_active_parcels!(state, ctx, cd, indt)
            Operators.reassign_escaped_parcels!(state, ctx, cd); State.sync_dtime_ns!(state)
            @inbounds for i in 1:NP; state.status[i]=0x01; state.Dtime[i]=cd+ts; end
        end
        @printf("%3d   %s   %9.2f   %-8s   %9.4f\n",
                it, Dates.format(cd, "mm-dd HH:MM"), t_io, rotated ? "LOAD" : "cached", t_adv)
        flush(stdout)
    end
    @printf("\nsummary: %d file load(s) in %d steps, total read_met=%.1f s  (~%.1f s/file)\n",
            nrot, nsteps, total_io, nrot > 0 ? total_io/nrot : 0.0)
end

run!(rt, ctx, state, NSTEPS, loopstart, ts, NP)
