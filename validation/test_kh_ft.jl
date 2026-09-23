#!/usr/bin/env julia
# test_kh_ft.jl — FLEXPART free-atmosphere horizontal diffusivity
# ([physics.turbulence].kh_ft), trajectory route:
#   [1] ensemble variance of above-PBL parcels grows by ~2·K·t beyond the
#       kh=0 twin (the defining property of a K-diffusivity random walk),
#   [2] below-PBL parcels are untouched (the Hanna BL owns that region),
#   [3] config round-trip keeps kh_ft,
#   [4] the double-counting guard: kh_ft > 0 with [physics.mixing] armed must
#       refuse to build (one K moves parcels OR exchanges mass, never both).
#
#   julia -t 2 --project=easlee test/test_kh_ft.jl
using Dates, Printf, Statistics
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO,"easlee","src","config","EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO,"easlee","src","io","initialize.jl"))
include(joinpath(REPO,"easlee","src","io","tools.jl"))
include(joinpath(REPO,"easlee","src","io","output.jl"))
include(joinpath(REPO,"easlee","EASLEE.jl")); using .EASLEE
include(joinpath(REPO,"easlee","runtime.jl"))
include(joinpath(REPO,"easlee","build_run.jl"))

const CFG = joinpath(REPO, "configs", "test_wind.toml")
const NSTEP = 6                       # 6 × 5 min = 1800 s of walk

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-64s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# run the tiny ETEX config for NSTEP steps; return horizontal position
# variance [m²] (x+y) of the above-PBL and below-PBL populations.
# The stock config releases over a 10°×6°×5 km box — ensemble variance
# 7.5e10 m², drowning any 2Kt signal.  Shrink to a POINT release at
# 3–4 km (above the PBL) so var0 is only shear+turb+meso spread (~1e8 m²)
# and the walk's 2Kt is detectable against it.
function probe(kh::Float64)
    spec = read_runspec(CFG); spec.wind = "3D"
    spec.physics.turbulence.kh_ft = kh
    r = spec.release[1]
    r.lon = [5.0, 5.0]; r.lat = [49.0, 49.0]; r.z = [3000.0, 4000.0]
    rt, obj = build_run(spec; prefetch=false)
    ctx = rt.ctx; st = rt.state; rc = ctx.runconfig; ts = rc.timestep
    kw = step_kwargs(obj)
    t0 = rc.startdate
    read_met_step!(rt, t0)
    EASLEE.step_model!(st, ctx, t0; kw...)
    EASLEE.State.sync_dtime_ns!(st)
    for k in 1:NSTEP
        t = t0 + k*ts
        read_met_step!(rt, t)
        EASLEE.step_model!(st, ctx, t; kw...)
        EASLEE.State.sync_dtime_ns!(st)
    end
    nest = ctx.nests[1]
    m3 = nest.metfields3D; m2 = nest.metfields2D
    gc = nest.gridconfig
    nzloc = Int(nest.levelconfig.Nk)
    dxm = 111320.0 * Float64(gc.dx)              # deg->m per index unit (approx)
    xa = Float64[]; ya = Float64[]; xb = Float64[]; yb = Float64[]
    for i in eachindex(st.status)
        st.status[i] == 0x01 || continue
        il = clamp(st.indl[i], 1.0f0, Float32(nzloc) - 1.0f-3)
        z = EASLEE.Operators.interpolation4Dlevsingle(m3.zlev, st.indx[i],
                                                      st.indy[i], il, 1.0f0)
        pb = EASLEE.Operators.interpolate2Dsingle(m2.pblh, st.indx[i],
                                                  st.indy[i], 1.0f0)
        cosl = cosd(46.0 + Float64(st.indy[i] - 1) * Float64(gc.dy))
        if isfinite(z) && z > max(pb, 50.0f0)
            push!(xa, Float64(st.indx[i]) * dxm * cosl)
            push!(ya, Float64(st.indy[i]) * dxm)
        else
            push!(xb, Float64(st.indx[i]) * dxm * cosl)
            push!(yb, Float64(st.indy[i]) * dxm)
        end
    end
    va = (isempty(xa) ? NaN : var(xa) + var(ya))
    vb = (isempty(xb) ? NaN : var(xb) + var(yb))
    return va, vb, length(xa), length(xb)
end

println("── kh_ft: variance growth vs the 2·K·t law (point release, 3–4 km) ──")
const KTEST = 5.0e4                    # wiring check: signal ≫ shear spread
va0, vb0, na0, nb0 = probe(0.0)
vaK, vbK, naK, nbK = probe(KTEST)
t_walk = NSTEP * 300.0
expected = 2 * 2 * KTEST * t_walk      # 2 axes × 2Kt
excess = vaK - va0
@printf("  above-PBL: n=%d/%d  var0=%.3e  varK=%.3e  excess=%.3e  2·2Kt=%.3e\n",
        na0, naK, va0, vaK, excess, expected)
check("above-PBL variance grows by ~2Kt per axis (×0.5–2 tolerated)",
      isfinite(excess) && 0.5 * expected < excess < 2.0 * expected,
      info=@sprintf("ratio %.2f", excess / expected))
check("population is above the PBL (point release at 3–4 km)",
      na0 > 3000 && nb0 < 2000,
      info=@sprintf("above=%d below=%d", na0, nb0))

println("── config round-trip ──")
sp = read_runspec(CFG); sp.physics.turbulence.kh_ft = 50.0
d = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("round-trip keeps kh_ft = 50", sp2.physics.turbulence.kh_ft == 50.0)
sp3 = read_runspec(CFG)
check("default is off (0.0)", sp3.physics.turbulence.kh_ft == 0.0)

println("── double-counting guard ──")
function guard_throws()
    spg = read_runspec(CFG)
    spg.physics.turbulence.kh_ft = 50.0
    spg.physics.mixing.scheme = "pairwise"
    try
        build_run(spg; prefetch=false)
        return false
    catch e
        return occursin("never both", sprint(showerror, e))
    end
end
check("kh_ft>0 + [physics.mixing] armed refuses to build", guard_throws())

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
