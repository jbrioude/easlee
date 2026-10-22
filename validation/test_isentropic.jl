# test_isentropic.jl — hybrid kinematic/isentropic vertical transport.
#
#   1. Radiative-tendency fields: loaded, stratospheric magnitudes ~ K/day.
#   2. θ conservation (adiabatic limit, diabatic=false): stratospheric parcels
#      must hold their isentrope; compare against the kinematic-only drift.
#   3. Diabatic drift: with real avg_ttlwr/ttswr the θ drift must be ~ the
#      integrated heating rate (K/day scale), not kinematic noise.
#   4. Blend zone: parcels spanning σ 0.5–0.05 — no pile-up at the seam.
#
# Usage:  julia --threads=8 --project=../easlee test_isentropic.jl

using Dates, Printf, Statistics, Random
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")

global gpu = false
npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-52s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

spec = read_runspec(joinpath(@__DIR__,"..","configs","era5_conv.toml")); spec.physics.vertical.mode = "isentropic"
rt, _obj = build_run(spec; prefetch=false, extra_slots=30_000)
ctx   = rt.ctx
state = rt.state
rc    = ctx.runconfig
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end
read_met_step!(rt, rc.startdate)

m3 = rt.metfields3D[1]
nz = size(m3["zlev"], 3)

# ── 1. tendency fields ────────────────────────────────────────────────────────
println("── radiative tendencies ──────────────────────────────────────────")
@assert haskey(m3, "ttlwr") && haskey(m3, "ttswr") "tendencies not in namelist"
# stratospheric sample: level with θ ≈ 500 K in a mid-domain column
thq = [m3["theta"][720, 240, k, 1] for k in 1:nz]
k500 = findfirst(k -> thq[k] > 500.0f0, 1:nz)
lw = m3["ttlwr"][:, :, k500, 1]; sw = m3["ttswr"][:, :, k500, 1]
@printf("  θ≈500 K (k=%d): LW %.2f K/day (median), SW %.2f K/day (median, max %.2f)\n",
        k500, 86400 * median(lw), 86400 * median(sw), 86400 * maximum(sw))
check("LW cooling at θ≈500 K", -5.0 < 86400 * median(lw) < -0.05,
      @sprintf("%.2f K/day", 86400 * median(lw)))
check("SW heating present (day side)", 86400 * maximum(sw) > 0.2,
      @sprintf("max %.2f K/day", 86400 * maximum(sw)))

# ── parcel placement helper: N parcels at θ ∈ [450, 550] K over a region ────
const NSTRAT = 5000
const P0     = 20_001          # first manual slot (after the RELEASES block)
function place_parcels!(rngseed)
    rng = MersenneTwister(rngseed)
    for q in 0:NSTRAT-1
        i = P0 + q
        ix = 700.0f0 + 200.0f0 * rand(rng, Float32)     # ~lon 0–50E of domain
        iy = 160.0f0 + 160.0f0 * rand(rng, Float32)     # mid-latitudes
        th_t = 450.0f0 + 100.0f0 * rand(rng, Float32)
        thp = [m3["theta"][round(Int, ix), round(Int, iy), k, 1] for k in 1:nz]
        kk = findfirst(k -> thp[k] > th_t, 2:nz)
        kk === nothing && (kk = nz - 2)
        k = kk + 1
        f = (th_t - thp[k-1]) / max(thp[k] - thp[k-1], 1.0f-3)
        state.status[i] = 0x01; state.nest[i] = 1; state.is_new[i] = false
        state.indx[i] = ix; state.indy[i] = iy
        # θ is a full-level field: slot k sits at indl k − ½ (layout contract)
        state.indl[i] = clamp(Float32(k) - 1.5f0 + f, 2.0f0, Float32(nz - 1))
        state.mass[i] = 1.0f0
        state.Dtime[i] = rc.startdate
    end
end
# Grid θ at the parcel, TIME-INTERPOLATED at `date` — measuring against a
# fixed met slot would alias the local ∂θ/∂t (≈ K/hour at 20 km) into the
# apparent parcel drift.
function grid_theta(i, date)
    n  = Int(state.nest[i])
    th = rt.metfields3D[n]["theta"]
    indt = EASLEE.MetAdapter.met_time_index_array(date, ctx.nests)[n]
    EASLEE.Operators.interpolation4Dsingle(th, state.indx[i], state.indy[i],   # class F: stacked
        clamp(state.indl[i], 1.0f0, Float32(nz - 1)), indt)
end
function clear_parcels!()
    for q in 0:NSTRAT-1
        state.status[P0+q] = 0x00; state.mass[P0+q] = 0.0f0
    end
end
function fresh_met!()
    for n in 1:rt.nest
        rt.ctx.nests[n].nestrunconfig.t_index = (0, 0)
    end
    read_met_step!(rt, rc.startdate)
end
function run_steps!(; hyb)
    ts = rc.timestep
    nsteps = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(ts)))
    for itime in 0:nsteps
        d = rc.startdate + itime * ts
        read_met_step!(rt, d)
        EASLEE.step_model!(state, ctx, d; do_output=false, hybridvert=hyb)
        State.sync_dtime_ns!(state)
    end
end
alive() = [P0+q for q in 0:NSTRAT-1 if state.status[P0+q] == 0x01]

# ── 2. kinematic baseline drift ───────────────────────────────────────────────
println("\n── θ drift: kinematic vs adiabatic-isentropic vs diabatic ───────")
place_parcels!(11)
th0 = Dict(i => grid_theta(i, rc.startdate) for i in alive())
run_steps!(hyb=nothing)
ids = [i for i in alive() if haskey(th0, i)]
drift_kin = [abs(grid_theta(i, rc.enddate) - th0[i]) for i in ids]
@printf("  kinematic        : median |Δθ| = %.2f K, p95 = %.2f K  (n=%d)\n",
        median(drift_kin), quantile(drift_kin, 0.95), length(ids))
clear_parcels!()

# ── 3. adiabatic isentropic: θ must be conserved ─────────────────────────────
fresh_met!()
place_parcels!(11)
hyb0 = EASLEE.HybridVert(length(state.status), rt; diabatic=false)
th0 = Dict(i => grid_theta(i, rc.startdate) for i in alive())
run_steps!(hyb=hyb0)
ids = [i for i in alive() if haskey(th0, i)]
drift_adi = [abs(grid_theta(i, rc.enddate) - th0[i]) for i in ids]
@printf("  adiabatic hybrid : median |Δθ| = %.3f K, p95 = %.2f K\n",
        median(drift_adi), quantile(drift_adi, 0.95))
check("θ conserved on isentropes", median(drift_adi) < 0.2 &&
      median(drift_adi) < 0.1 * max(median(drift_kin), 0.5),
      @sprintf("median %.3f K vs kinematic %.2f K", median(drift_adi), median(drift_kin)))
clear_parcels!()

# ── 4. diabatic drift consistent with heating rates ───────────────────────────
fresh_met!()
place_parcels!(11)
hyb1 = EASLEE.HybridVert(length(state.status), rt; diabatic=true)
th0 = Dict(i => grid_theta(i, rc.startdate) for i in alive())
run_steps!(hyb=hyb1)
ids = [i for i in alive() if haskey(th0, i)]
sdrift = [grid_theta(i, rc.enddate) - th0[i] for i in ids]
@printf("  diabatic hybrid  : median Δθ = %+.3f K, median |Δθ| = %.3f K over the run\n",
        median(sdrift), median(abs.(sdrift)))
check("diabatic drift is K/day-scale, not kinematic noise",
      median(abs.(sdrift)) < 0.2 * max(median(drift_kin), 0.5) &&
      median(abs.(sdrift)) > 0.0005,
      @sprintf("median |Δθ| %.3f K (kinematic %.2f K)", median(abs.(sdrift)), median(drift_kin)))
clear_parcels!()

# ── PT-level vs model-level winds (the along-isentrope sampling error) ───────
println("\n── isentropic vs model-level winds ───────────────────────────────")
fresh_met!()
place_parcels!(11)
if hyb1.pt !== nothing
    EASLEE.Operators.read_pt_step!(hyb1.pt, rc.startdate)
    du = Float32[]
    for i in alive()
        gc = ctx.nests[1].gridconfig
        lat, lon = EASLEE.Operators.indyindxtolatlon(gc, state.indy[i], state.indx[i])
        thp = grid_theta(i, rc.startdate)
        uth, vth = EASLEE.Operators._pt_uv(hyb1.pt, Float32(lat), Float32(lon),
                                           Float32(thp), rc.startdate)
        uml = EASLEE.Operators.interpolateU(gc.Uedge, ctx.nests[1].metfields3D.Ux,
                                            state.indx[i], state.indy[i],
                                            clamp(state.indl[i], 1.0f0, Float32(nz-1)), 1.0f0)
        vml = EASLEE.Operators.interpolateV(gc.Vedge, ctx.nests[1].metfields3D.Vy,
                                            state.indx[i], state.indy[i],
                                            clamp(state.indl[i], 1.0f0, Float32(nz-1)), 1.0f0)
        (isfinite(uth) && isfinite(uml)) || continue
        push!(du, sqrt((uth - uml)^2 + (vth - vml)^2))
    end
    @printf("  |u_θ − u_ml| at θ∈[450,550] K: median %.2f m/s, p95 %.2f m/s (n=%d)\n",
            median(du), quantile(du, 0.95), length(du))
    check("PT winds loaded & consistent with model levels",
          length(du) > 1000 && 0.0 < median(du) < 5.0,
          @sprintf("median %.2f m/s", median(du)))
else
    check("PT winds loaded", false, "no IFS_PT_MET files found")
end
clear_parcels!()

# ── 5. blend zone: no pile-up at the kinematic/isentropic seam ───────────────
println("\n── blend zone σ ∈ [0.05, 0.5] ────────────────────────────────────")
fresh_met!()
rng = MersenneTwister(7)
for q in 0:NSTRAT-1
    i = P0 + q
    ix = 700.0f0 + 200.0f0 * rand(rng, Float32)
    iy = 160.0f0 + 160.0f0 * rand(rng, Float32)
    sg = 0.05f0 + 0.45f0 * rand(rng, Float32)
    pp = [m3["puv"][round(Int, ix), round(Int, iy), k, 1] for k in 1:nz]
    tgt = sg * pp[1]
    kk = findfirst(k -> pp[k] < tgt, 2:nz); kk === nothing && (kk = nz - 2)
    k = kk + 1
    f = (pp[k-1] - tgt) / max(pp[k-1] - pp[k], 1.0f0)
    state.status[i] = 0x01; state.nest[i] = 1; state.is_new[i] = false
    state.indx[i] = ix; state.indy[i] = iy
    state.indl[i] = clamp(Float32(k) - 1.5f0 + f, 2.0f0, Float32(nz - 1))   # puv: class F
    state.mass[i] = 1.0f0; state.Dtime[i] = rc.startdate
end
sigma_of(i) = begin
    n = Int(state.nest[i]); pv = rt.metfields3D[n]["puv"]
    ii = clamp(round(Int, state.indx[i]), 1, size(pv, 1) - 1)
    jj = clamp(round(Int, state.indy[i]), 1, size(pv, 2) - 1)
    EASLEE.Operators.interpolation4Dsingle(pv, Float32(ii), Float32(jj),
        clamp(state.indl[i], 1.0f0, Float32(nz - 1)), 1.0f0) / pv[ii, jj, 1, 1]
end
edges = range(0.05, 0.5, length=10)
h0 = [count(i -> edges[b] <= sigma_of(i) < edges[b+1], alive()) for b in 1:9]
hyb2 = EASLEE.HybridVert(length(state.status), rt; diabatic=true)
run_steps!(hyb=hyb2)
h1 = [count(i -> edges[b] <= sigma_of(i) < edges[b+1], alive()) for b in 1:9]
println("  σ bin     before   after")
for b in 1:9
    @printf("  %.2f–%.2f  %6d  %6d\n", edges[b], edges[b+1], h0[b], h1[b])
end
ratio = maximum(h1[b] / max(h0[b], 1) for b in 2:8)
check("no pile-up at the blend seam", ratio < 1.5,
      @sprintf("max bin growth ×%.2f", ratio))

println()
@printf("=== %d passed, %d failed ===\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
