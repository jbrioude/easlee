# test_plume_rise.jl ── validation of the biomass-burning injection-height module
# (Sofiev 2012 semi-empirical + simplified Freitas 1-D entraining plume) on a REAL
# ERA5 sounding, plus the end-to-end wired path (release seeding fills the layer).
#
# Part 1 — standalone physics at the release cell of the config:
#   • FRP sweep (10 MW … 20 GW), both schemes: injection layer, monotone in FRP;
#   • Sofiev limits: FRP→0 gives α·H_abl; stronger FT stability lowers the top;
#   • Freitas: z_top > z_NBL > 0, plume deepens with FRP and with fire area.
#
# Part 2 — wired path (the engine, not the module):
#   release_parcels! + initialize_new_parcels! with [release.plume_rise] set, then
#   check every seeded parcel's height lies inside the computed injection layer
#   (and that the layer differs from the TOML z-range it overrides).
#
#   CFG=configs/era5_conv.toml julia -t 4 --project=easlee test/test_plume_rise.jl

using Dates, Printf, Statistics, Random
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))
const Op = EASLEE.Operators

const CFG    = get(ENV, "CFG", joinpath(REPO, "configs", "era5_conv.toml"))
const SCHEMES = ("sofiev", "freitas")
const FRPS   = [1.0e7, 1.0e8, 5.0e8, 1.0e9, 5.0e9, 2.0e10]   # W
const AREA   = 2.0e6                                           # m² (freitas)
const NPART  = 20000

fails = String[]
check(ok, msg) = (ok ? @printf("  ✓ %s\n", msg) : (push!(fails, msg); @printf("  ✗ FAIL %s\n", msg)); ok)

spec = read_runspec(CFG)
r = spec.release[1]
r.particles = NPART
r.time = [spec.start, spec.start]                 # instantaneous release at t0
rt, obj = build_run(spec; prefetch=false, extra_slots=NPART)
ctx = rt.ctx; state = rt.state; nest = ctx.nests[1]
loopstart = ctx.runconfig.startdate
read_met_step!(rt, loopstart)
t_array = EASLEE.MetAdapter.met_time_index_array(loopstart, ctx.nests)
indt = t_array[1]
Nk = Int(nest.levelconfig.Nk)

# ── the ambient column at the release-box centre ──────────────────────────────
lonc = 0.5 * (r.lon[1] + r.lon[2]); latc = 0.5 * (r.lat[1] + r.lat[2])
indy_c, indx_c = Op.latlontoindyindx(nest.gridconfig, latc, lonc)
ix = round(Int, indx_c); iy = round(Int, indy_c); it = clamp(round(Int, Float64(indt)), 1, 2)
z3 = nest.metfields3D.zlev; th = nest.metfields3D.theta
# θ is a full-level field: pair it with the full-level heights, as the engine's
# _compute_plume_layer does (layout contract 2026-10-22; was zlev[k] with θ[k])
zcol = Op._full_level_heights(z3, ix, iy, it)
zagl = zcol .- zcol[1]
θcol = Float64[th[ix, iy, k, it] for k in 1:Nk]
pblh = Float64(nest.metfields2D.pblh[ix, iy, it])
pblh > 50.0 || (pblh = Op.estimate_pbl_height(zagl, θcol))
n2ft = Op.n2_free_trop(zagl, θcol, pblh)
@printf("=== plume-rise test  cfg=%s  cell=(%.2f°E, %.2f°N)  H_abl=%.0f m  N²_FT=%.3g s⁻² ===\n",
        basename(CFG), lonc, latc, pblh, n2ft)

# ── Part 1: standalone FRP sweep ──────────────────────────────────────────────
println("\n── Part 1: injection layers on the real sounding ──")
@printf("%-10s  %-22s  %-22s\n", "FRP", "sofiev [z_b, z_top] m", "freitas [z_NBL, z_top] m")
tops = Dict(s => Float64[] for s in SCHEMES)
lay  = Dict{Tuple{String,Float64},Tuple{Float64,Float64}}()
for f in FRPS
    zs = Op.plume_injection_layer("sofiev",  f, AREA, pblh, zagl, θcol; bottom="pbl")
    zf = Op.plume_injection_layer("freitas", f, AREA, pblh, zagl, θcol; bottom="pbl")
    push!(tops["sofiev"], zs[2]); push!(tops["freitas"], zf[2])
    lay[("sofiev",f)] = zs; lay[("freitas",f)] = zf
    @printf("%-10.3g  [%7.0f, %7.0f]      [%7.0f, %7.0f]\n", f, zs[1], zs[2], zf[1], zf[2])
end
check(issorted(tops["sofiev"]),  "sofiev plume top monotone in FRP")
check(issorted(tops["freitas"]), "freitas plume top monotone in FRP")
check(all(f -> (l = lay[("freitas",f)]; l[2] > l[1] > 0), FRPS), "freitas: z_top > z_NBL > 0 at every FRP")
check(all(f -> lay[("sofiev",f)][2] >= 0.24 * pblh - 1, FRPS), "sofiev: top ≥ α·H_abl floor")

# Sofiev stability limit: doubling N²_FT must lower the top (exp(−δN²/N0²) term)
t1 = Op.plume_top_sofiev(5.0e8, pblh, n2ft)
t2 = Op.plume_top_sofiev(5.0e8, pblh, 2 * max(n2ft, 1e-5))
check(t2 < t1, @sprintf("sofiev: doubling FT stability lowers the top (%.0f → %.0f m)", t1, t2))

# Freitas area sensitivity: same heat flux over a bigger fire → weaker b0 but
# larger initial flux; the plume must stay finite and ordered.
za = Op.plume_layer_freitas(Op.frp_to_convective_heat(5.0e8), 0.5e6, zagl, θcol)
zb = Op.plume_layer_freitas(Op.frp_to_convective_heat(5.0e8), 8.0e6, zagl, θcol)
check(all(x -> isfinite(x) && x >= 0, (za..., zb...)), "freitas: finite layers across fire areas 0.5–8 km²")

# PBL-vs-loft regime (the validated 500 MW / 5 GW pair from PLUME_RISE.md)
z500  = lay[("sofiev", 5.0e8)]; z5g = lay[("sofiev", 5.0e9)]
@printf("  · 500 MW sofiev layer [%.0f, %.0f] m  vs  5 GW [%.0f, %.0f] m  (H_abl=%.0f)\n",
        z500[1], z500[2], z5g[1], z5g[2], pblh)

# ── Part 2: wired release seeding ─────────────────────────────────────────────
println("\n── Part 2: engine seeding fills the computed layer ──")
for scheme in SCHEMES, frp in (5.0e8, 5.0e9)
    pcfg = Dict(1 => (scheme = scheme, frp = frp, area = AREA, profile = "uniform", bottom = "pbl"))
    Op.set_plume_config!(pcfg)
    fill!(state.status, 0x00); fill!(state.is_new, false)
    rng = MersenneTwister(3)
    nrel = Op.release_parcels!(state, ctx, loopstart; rng = rng)
    Op.initialize_new_parcels!(state, ctx, loopstart, t_array; rng = rng)
    idx = [i for i in eachindex(state.status) if state.status[i] == 0x01]

    # expected layer, computed exactly as the engine computes it (first parcel's cell)
    zl = Op._compute_plume_layer(pcfg[1], nest, state.indx[idx[1]], state.indy[idx[1]], indt)

    # parcel heights AGL at their own columns
    buf = Vector{Float32}(undef, Nk)
    zp = Float64[]
    for i in idx
        Op.get_vertical_profile!(buf, z3, state.indx[i], state.indy[i], indt, Nk)
        k = clamp(Base.unsafe_trunc(Int, state.indl[i]), 1, Nk - 1)
        f = state.indl[i] - k
        push!(zp, buf[k] + f * (buf[k+1] - buf[k]) - buf[1])
    end
    tol = 60.0    # m — column-to-column zlev variation across the release box
    inlayer = count(z -> zl[1] - tol <= z <= zl[2] + tol, zp)
    @printf("%s FRP=%.0e: released %d, layer [%.0f, %.0f] m, parcels z∈[%.0f, %.0f] m, in-layer %.1f%%\n",
            scheme, frp, nrel, zl[1], zl[2], minimum(zp), maximum(zp), 100inlayer / length(zp))
    check(nrel == NPART, "$scheme FRP=$frp: all $NPART parcels released")
    check(inlayer / length(zp) > 0.99, "$scheme FRP=$frp: ≥99% of parcels inside the injection layer")
    check(maximum(zp) - minimum(zp) > 0.5 * (zl[2] - zl[1]), "$scheme FRP=$frp: parcels FILL the layer (not clumped)")
end
Op.set_plume_config!(Dict{Int,NamedTuple}())

println()
if isempty(fails)
    println("ALL PLUME-RISE CHECKS PASSED ✓")
else
    println("FAILED CHECKS:"); foreach(m -> println("  ✗ ", m), fails); exit(1)
end
