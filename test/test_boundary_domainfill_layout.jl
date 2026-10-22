#!/usr/bin/env julia
# test_boundary_domainfill_layout.jl — the stacked vertical layout (contract
# 2026-10-22) at the consumers that pair layers or heights with full-level
# (class F) fields:
#   [1] compute_domain_mass: layer [zlev[k], zlev[k+1]) × ρ of full level k+1,
#   [2] inject_boundary_parcels!: inflow of a wind confined to ONE full level
#       lands in the ONE layer that holds it,
#   [3] build_domainfill / domainfill_init!: per-layer weights ρ[k+1]·Δz, kept
#       layers by full-level height, parcels placed in the weighted layers,
#   [4] domainfill_inject!: same single-level inflow check on the box faces,
#   [5] plume-rise / volcanic columns: class-F values paired with the full-level
#       heights zf(k) = ½(zlev[k-1] + zlev[k]), zf(1) = ground,
#   [6] the Main-scope interpolation4Dsingle (io/tools.jl) equals the engine's.
# Synthetic met: half levels every 1000 m, a distinct ρ per slot, so a read
# half a level off is visible in every number.
#
#   julia --project=easlee test/test_boundary_domainfill_layout.jl
using Printf, Dates, Random

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
const Op = EASLEE.Operators

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic met ─────────────────────────────────────────────────────────────
const NI, NJ, NZ = 12, 10, 8
const DZ = 1000.0
dims = (NI, NJ, NZ, 2)
zlev = zeros(Float32, dims)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * DZ); end       # half levels, zlev[1] = 0
rhos = [1.0 + 0.1 * s for s in 1:NZ]                                 # slot s (full level s)
rho  = zeros(Float32, dims); for s in 1:NZ; rho[:, :, s, :] .= Float32(rhos[s]); end
Ux = zeros(Float32, dims); Vy = zeros(Float32, dims)
θs = [s <= 3 ? 300.0 : 310.0 for s in 1:NZ]                          # jump between full levels 3 and 4
Ts = [300.0 - s for s in 1:NZ]
ps = [1.0e5 - 1.0e4 * (s - 1) for s in 1:NZ]
theta = zeros(Float32, dims); T = zeros(Float32, dims); puv = zeros(Float32, dims)
for s in 1:NZ
    theta[:, :, s, :] .= Float32(θs[s]); T[:, :, s, :] .= Float32(Ts[s]); puv[:, :, s, :] .= Float32(ps[s])
end
m3 = (Ux = Ux, Vy = Vy, rho = rho, zlev = zlev, theta = theta, T = T, puv = puv)
m2 = (pblh = zeros(Float32, NI, NJ, 2),)                             # 0 → PBL top from θ
const DLON = 0.25f0; const DLAT = 0.25f0
GC = EASLEE.Types.Rectilinear(Int32(NI), Int32(NJ), Int32(NZ), DLON, DLAT, 0.0f0, -1.0f0,
        Float32(1 / (DLON * 111320.0)), Float32(1 / (DLAT * 111320.0)),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, true, DLON, DLAT, false)
nest = (id = 1, gridconfig = GC, metfields3D = m3, metfields2D = m2,
        levelconfig = (Nk = Int32(NZ),))
const DT_MS = 600_000
ctx = (runconfig = (timestep = Millisecond(DT_MS),),
       output = (Niout = Int32(NI - 1), Njout = Int32(NJ - 1)),
       nests = [nest])
zf = [k == 1 ? 0.0 : (k - 1.5) * DZ for k in 1:NZ]                   # full-level heights

println("="^78); println("[1] compute_domain_mass"); println("="^78)
area = 0.0
for j in 1:(NJ - 1), i in 1:(NI - 1)
    global area += Float64(Op._face_dx_m(GC, Float32(i) + 0.5f0, Float32(j) + 0.5f0)) *
                   Float64(Op._face_dy_m(GC, Float32(i) + 0.5f0, Float32(j) + 0.5f0))
end
col_new = sum(rhos[k + 1] * DZ for k in 1:(NZ - 1))                  # ρ of the full level each layer holds
col_old = sum(0.5 * (rhos[k] + rhos[k + 1]) * DZ for k in 1:(NZ - 1))
M = Float64(EASLEE.compute_domain_mass(nest, 1.0f0))
@printf("  column mass: model %.2f kg/m²   contract %.2f   old co-located read %.2f\n",
        M / area, col_new, col_old)
check("domain air mass = Σ ρ[k+1]·Δz per column", isapprox(M / area, col_new; rtol = 1e-4))

println("="^78); println("[2] inject_boundary_parcels!: wind on full level 5 only"); println("="^78)
Ux[:, :, 5, :] .= 10.0f0                         # full level 5 = layer [zlev[4], zlev[5]) = indl ∈ [4, 5)
dy_m = Float64(Op._face_dy_m(GC, 1.0f0, 1.5f0))
flux_cell = rhos[5] * 10.0 * dy_m * DZ * DT_MS / 1000
mpp = Float32(flux_cell / 50)                    # 50 parcels per face cell
st = EASLEE.State.ParcelState(2000)
n = EASLEE.inject_boundary_parcels!(st, ctx, DateTime(2025, 1, 1), [1.0f0], mpp;
                                    rng = MersenneTwister(1))
act = [i for i in eachindex(st.status) if st.status[i] == 0x01]
nl4 = count(i -> 4.0f0 <= st.indl[i] < 5.0f0, act)
@printf("  injected %d (expected %d ± %d), in layer [4,5): %d, elsewhere: %d\n",
        n, 50 * (NJ - 1), NJ - 1, nl4, length(act) - nl4)
check("injected count = ρ[5]·u[5]·dy·Δz·dt / m per face cell", abs(n - 50 * (NJ - 1)) <= NJ - 1)
check("every injected parcel is in the layer of full level 5 (indl ∈ [4,5))",
      n > 0 && nl4 == length(act) == n)
check("injected on the west face (u > 0 inflow)", all(i -> st.indx[i] < 1.5f0, act))

println("="^78); println("[3] build_domainfill / domainfill_init!"); println("="^78)
Ux .= 0.0f0
lonb = [0.5, 2.0]; latb = [-0.5, 0.5]; ztop = 3200.0; NP = 3900
df = EASLEE.build_domainfill(nest, 1.0f0, lonb, latb, ztop, NP)
@printf("  keep %s  weights %s  (contract ρ[k+1]·Δz = %s)\n", string(df.keep_lev),
        string(round.(df.lev_w; digits = 1)), string([rhos[k + 1] * DZ for k in df.keep_lev]))
check("kept layers = those whose full level (mid) ≤ ztop: [1, 2, 3]", df.keep_lev == [1, 2, 3])
check("layer weights = ρ[k+1]·(zlev[k+1] − zlev[k])",
      all(isapprox(df.lev_w[i], rhos[k + 1] * DZ; rtol = 1e-5) for (i, k) in enumerate(df.keep_lev)))
check("lev_h = half-level heights zlev[k]", df.lev_h ≈ [(k - 1) * DZ for k in 1:NZ])
# box air mass with the same horizontal sweep as build_domainfill
abox = 0.0; latq = latb[1] + DLAT / 2
while latq < latb[2]
    lonq = lonb[1] + DLON / 2
    while lonq < lonb[2]
        global abox += (DLAT * 111_320.0) * (DLON * 111_320.0 * cosd(latq)); lonq += DLON
    end
    global latq += DLAT
end
mpp_exp = abox * sum(rhos[k + 1] * DZ for k in 1:3) / NP
@printf("  mass_per_parcel %.6g kg (contract %.6g)\n", df.mass_per_parcel, mpp_exp)
check("mass_per_parcel = box air mass of the kept layers / nparts",
      isapprox(Float64(df.mass_per_parcel), mpp_exp; rtol = 1e-5))
st = EASLEE.State.ParcelState(NP + 2000)
np_ = EASLEE.domainfill_init!(st, ctx, df, DateTime(2025, 1, 1), 1.0f0; rng = MersenneTwister(2))
act = [i for i in eachindex(st.status) if st.status[i] == 0x01]
cnt = [count(i -> k <= st.indl[i] < k + 1, act) for k in 1:3]
@printf("  placed %d: per layer %s (weights ∝ %s)\n", np_, string(cnt),
        string(round.(NP .* df.lev_w ./ sum(df.lev_w); digits = 1)))
check("all parcels placed, inside the kept layers", np_ == NP && sum(cnt) == NP)
check("parcels per layer ∝ ρ[k+1]·Δz", cnt == [1200, 1300, 1400])

println("="^78); println("[4] domainfill_inject!: wind on full level 3 only"); println("="^78)
Ux[:, :, 3, :] .= 10.0f0                         # full level 3 = layer 2 = indl ∈ [2, 3)
nbefore = length(act)
ninj = EASLEE.domainfill_inject!(st, ctx, df, DateTime(2025, 1, 1), 1.0f0, DT_MS / 1000;
                                 rng = MersenneTwister(3))
new = [i for i in eachindex(st.status) if st.status[i] == 0x01 && !(i in act)]
nexp = 0.0; latq = latb[1] + DLAT / 2
while latq < latb[2]
    global nexp += rhos[3] * 10.0 * (DLAT * 111_320.0) * DZ * (DT_MS / 1000) / df.mass_per_parcel
    global latq += DLAT
end
nl2 = count(i -> 2.0f0 <= st.indl[i] < 3.0f0, new)
@printf("  injected %d (expected %.1f), in layer [2,3): %d, elsewhere: %d\n",
        ninj, nexp, nl2, length(new) - nl2)
check("injected count = Σ ρ[3]·u[3]·dy·Δz·dt / m", abs(ninj - nexp) <= 4 + 3 * sqrt(nexp))
check("every injected parcel is in the layer of full level 3 (indl ∈ [2,3))",
      ninj > 0 && nl2 == length(new) == ninj)
Ux .= 0.0f0

println("="^78); println("[5] plume-rise / volcanic columns: full-level heights"); println("="^78)
zc = Op._full_level_heights(zlev, 4, 4, 1)
check("_full_level_heights = [ground, ½(zlev[k-1]+zlev[k])...]", zc ≈ zf; info = string(zc[1:4]))
zv, Tv, pv, qv = Op._volcanic_column_arrays(nest, 4.0f0, 4.0f0, 1.0f0)
check("volcanic column: z = zf(k) paired with T, p of slot k",
      zv ≈ zf && Tv ≈ Ts && pv ≈ ps && length(qv) == NZ)
lay = Op._compute_plume_layer((scheme = "sofiev", frp = 2.0e10, area = 1.0e6, bottom = "pbl"),
                              nest, 4.0f0, 4.0f0, 1.0f0)
@printf("  sofiev layer [%.0f, %.0f] m — bottom = PBL top from θ: contract zf(4) = %.0f m, old zlev[4] = %.0f m\n",
        lay[1], lay[2], zf[4], (4 - 1) * DZ)
check("plume layer bottom = θ-derived PBL top at the full-level height zf(4)",
      isapprox(lay[1], zf[4]; atol = 1e-6) && lay[2] > lay[1])

println("="^78); println("[6] Main-scope interpolation4Dsingle (io/tools.jl)"); println("="^78)
f4 = rand(MersenneTwister(4), Float32, 5, 5, NZ, 2)
maxd = 0.0f0
for z in (1.0f0, 1.2f0, 1.5f0, 2.3f0, 5.0f0, Float32(NZ) - 0.5f0, Float32(NZ) - 0.2f0),
    t in (1.0f0, 1.3f0, 2.0f0)
    global maxd = max(maxd, abs(interpolation4Dsingle(f4, 2.3f0, 3.6f0, z, t) -
                                Op.interpolation4Dsingle(f4, 2.3f0, 3.6f0, z, t)))
end
check("Main copy = engine read (surface ramp 2z−1, top clamp, t = 2)", maxd <= 1.0f-6;
      info = @sprintf("(max |Δ| = %.2g)", maxd))
check("indl = 1 reads the surface slot 1", interpolation4Dsingle(f4, 2.0f0, 2.0f0, 1.0f0, 1.0f0) == f4[2, 2, 1, 1])

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
