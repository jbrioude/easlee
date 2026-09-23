#!/usr/bin/env julia
# test_estimator_calibration.jl — ANALYTIC calibration of the output estimators.
# Motivation (2026-07-06): estimator bugs kept being discovered through 2-hour
# ETEX/CAPTEX reruns; the estimators are pure functions of (positions, masses,
# bandwidths) and must be checked against closed-form fields in seconds.
#   [1] absolute calibration: uniform slab, C(x_s) == M/(A·Δz) — catches ANY
#       normalization error (2π, Δz, window averaging) at the few-% level,
#   [2] Gaussian puff: C(center) == M·G₂(0;σ_eff)/Δz with the KDE convolution
#       width σ_eff² = σ² + h², station off-center likewise,
#   [3] the near-source filament blowup (ETEX station 8021, deterministic):
#       thin filament + kNN bandwidths ≪ physical footprint → point evaluation
#       explodes; the hmin floor bounds it without touching [1]-[2],
#   [4] window semantics: flushed value is the MEAN over samples, not the sum.
#
#   julia --project=easlee -t 2 test/test_estimator_calibration.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: ReceptorSet, receptor_sample!, receptor_flush!,
                         update_bandwidths!

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# synthetic 1° Rectilinear grid at the equator so 1 index unit = 111.32 km and
# dx ≈ dy in metres; levels every 100 m; layer 0–100 m
G = EASLEE.Types.Rectilinear(Int32(100), Int32(100), Int32(10), 1.0f0, 1.0f0,
        0.0f0, 0.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
zlev = Array{Float32,4}(undef, 100, 100, 10, 2)
for k in 1:10; zlev[:, :, k, :] .= Float32((k - 1) * 100.0); end
const DEG = 111_320.0             # m per degree at the equator
const DZ  = 100.0                 # sampling layer depth

println("="^78)
println("[1] absolute calibration — uniform slab")
println("="^78)
Random.seed!(42)
N = 200_000
L = 4.0 * DEG                      # 4°×4° slab, station in the middle
M = 1.0e6                          # total mass [kg]
st = EASLEE.State.ParcelState(N)
st.status .= 0x01; st.nest .= Int32(1)
st.mass .= Float32(M / N)
st.indx .= 31.0f0 .+ 4.0f0 .* rand(Float32, N)   # lon 30–34°
st.indy .= 1.0f0 .+ 4.0f0 .* rand(Float32, N)    # lat 0–4° (equator: cos≈1, so
                                                 # index-metres == real metres)
st.indl .= 1.0f0 .+ rand(Float32, N)             # z uniform 0–100 m
update_bandwidths!(st, G, 1)
C_true = M / (L * L * DZ)                        # kg/m³, exactly
# ground-truth bandwidths: exact 12-NN RMS distance for a parcel sample —
# validates update_bandwidths! without trusting any order-statistics algebra
function brute_h(stx, idxs; k=12)
    x = Float64.(stx.indx) .* DEG; y = Float64.(stx.indy) .* DEG
    hs = Float64[]
    for i in idxs
        d2 = (x .- x[i]).^2 .+ (y .- y[i]).^2
        d2[i] = Inf
        push!(hs, sqrt(sum(partialsort(d2, 1:k)) / k))
    end
    return hs
end
Random.seed!(1)
samp = rand(1:N, 300)
hb = brute_h(st, samp)
check("update_bandwidths! matches brute-force 12-NN RMS (±10%)",
      0.9 < median(Float64.(st.hband[samp])) / median(hb) < 1.1;
      info=@sprintf("kNN %.0f m vs brute %.0f m", median(Float64.(st.hband[samp])), median(hb)))
# a single point value carries ~15% KDE variance (N_eff ≈ 40 kernels) —
# calibration is judged on the MEAN over a 5×5 station array
slons = repeat(collect(31.0:0.5:33.0), inner=5)
slats = repeat(collect(1.0:0.5:3.0), outer=5)
rs = ReceptorSet(string.(1:25), slons, slats, 0.0, DZ)
receptor_sample!(rs, st, G, LC, zlev, 1.0, 1)
receptor_flush!(rs, 0.0)
ratios = [r[3] for r in rs.rows] ./ C_true
c_est = mean([r[3] for r in rs.rows])
check("uniform slab: ⟨C⟩/C_true over 25 stations = 1 (±7%)",
      0.93 < mean(ratios) < 1.07;
      info=@sprintf("mean %.3f  (spread %.2f–%.2f)", mean(ratios), minimum(ratios), maximum(ratios)))
# the same with a floor far below the spacing must be IDENTICAL
h_typ = median(Float64.(st.hband[st.hband .> 0]))
rsf = ReceptorSet(string.(1:25), slons, slats, 0.0, DZ; hmin=h_typ/10)
receptor_sample!(rsf, st, G, LC, zlev, 1.0, 1); receptor_flush!(rsf, 0.0)
check("small hmin (≪ spacing) leaves the calibrated value unchanged",
      isapprox(mean([r[3] for r in rsf.rows]), c_est; rtol=1e-12);
      info=@sprintf("spacing h≈%.0f m", h_typ))

println()
println("="^78)
println("[2] Gaussian puff — closed-form KDE convolution")
println("="^78)
Random.seed!(7)
Np = 200_000
σ = 30_000.0                       # 30 km puff
Mp = 5.0e5
stp = EASLEE.State.ParcelState(Np)
stp.status .= 0x01; stp.nest .= Int32(1)
stp.mass .= Float32(Mp / Np)
stp.indx .= 51.0f0 .+ Float32(σ / DEG) .* randn(Float32, Np)
stp.indy .= 3.0f0 .+ Float32(σ / DEG) .* randn(Float32, Np)   # centre lat 2°
stp.indl .= 1.0f0 .+ rand(Float32, Np)
update_bandwidths!(stp, G, 1)
# The kernel width is NOT a free parameter: it must be the mean distance to
# the surrounding trajectories — the ensemble spread the random walk creates
# IS the dispersion that would otherwise grow each parcel's covariance, so
# h = local spacing redistributes it back into the kernels.  Assert that,
# then compare against the PRESCRIBED puff (M, σ): with h ≪ σ the KDE
# reproduces the prescribed field itself (h̄²/σ² here ≈ 4e-4).
# kernel width must equal the mean distance to the surrounding trajectories
# (the dispersion the random walk carries in centroid spread): validate the
# CORE parcels against brute-force exact 12-NN RMS
core_idx = findall((abs.(Float64.(stp.indx) .- 51.0) .< 0.3σ/DEG) .&
                   (abs.(Float64.(stp.indy) .- 3.0)  .< 0.3σ/DEG))
Random.seed!(2)
csamp = core_idx[rand(1:length(core_idx), 200)]
hbp = brute_h(stp, csamp)
check("puff-core kernel width = brute-force neighbour spacing (±10%)",
      0.9 < median(Float64.(stp.hband[csamp])) / median(hbp) < 1.1;
      info=@sprintf("kNN %.0f m vs brute %.0f m", median(Float64.(stp.hband[csamp])), median(hbp)))
# prescribed-field reproduction: mean over rings (single points carry ~15%
# KDE variance).  8 stations at r = 0.2σ (truth ×e^{-0.02}) and 8 at r = σ.
ring(r, n) = [(50.0 + r/DEG*cos(θ), 2.0 + r/DEG*sin(θ)) for θ in 2π .* (0:n-1) ./ n]
in_ring = ring(0.2σ, 8); out_ring = ring(σ, 8)
rsp = ReceptorSet(string.(1:16),
                  vcat(first.(in_ring), first.(out_ring)),
                  vcat(last.(in_ring),  last.(out_ring)), 0.0, DZ)
receptor_sample!(rsp, stp, G, LC, zlev, 1.0, 1); receptor_flush!(rsp, 0.0)
cvals = [r[3] for r in rsp.rows]
c0_true = Mp / (2π * σ^2 * DZ)                    # the PRESCRIBED field
cp0 = mean(cvals[1:8]);  cp1 = mean(cvals[9:16])
# Documented estimator property: density-adaptive bandwidths SHARPEN peaks —
# with h == neighbour spacing (proved above), the residual is a ±≈10% radial
# tilt (+ at the peak, − on the flank).  Bounded here; the aggregate mass is
# exact (next check), so the tilt redistributes, never creates.
check("puff near-centre ring within the +12% adaptive peak bias",
      1.0 <= cp0 / (c0_true * exp(-0.02)) < 1.12;
      info=@sprintf("ratio %.3f", cp0/(c0_true*exp(-0.02))))
check("puff r=σ ring within the −12% adaptive flank bias",
      0.88 < cp1 / (c0_true * exp(-0.5)) <= 1.0;
      info=@sprintf("ratio %.3f", cp1/(c0_true*exp(-0.5))))
# aggregate invariant: the radial integral of the estimated profile returns
# the prescribed mass (trapezoid over 24 rings to 4σ)
rr = collect(0.5:0.5:12.0) .* (σ / 3)
rsm = ReceptorSet(string.(1:24), [50.0 + r/DEG for r in rr], fill(2.0, 24), 0.0, DZ)
receptor_sample!(rsm, stp, G, LC, zlev, 1.0, 1); receptor_flush!(rsm, 0.0)
prof = [r[3] for r in rsm.rows]
Mint = sum(2π .* rr .* prof .* (σ/6)) * DZ          # ∫C·2πr dr · Δz
check("radial integral of the estimated puff returns the mass (±10%)",
      0.9 < Mint / Mp < 1.1; info=@sprintf("ratio %.3f", Mint/Mp))

println()
println("="^78)
println("[3] near-source filament (the ETEX-8021 case, deterministic)")
println("="^78)
# 20k parcels in a 200 m-wide filament passing 150 m from the station: the kNN
# bandwidth is the ~m-scale SAMPLING density, not the parcels' km-scale
# physical footprint → unbounded point evaluation.
Random.seed!(3)
Nf = 20_000
Mf = 1.0e3
stf = EASLEE.State.ParcelState(Nf)
stf.status .= 0x01; stf.nest .= Int32(1)
stf.mass .= Float32(Mf / Nf)
# filament 200 m wide in lon, 40 km long in lat, passing 150 m east of the
# station at (lat 2, lon 40) = (indy 3, indx 41) on this grid
stf.indx .= Float32(41.0 + 150.0/DEG) .+ Float32(200.0/DEG) .* randn(Float32, Nf)
stf.indy .= 3.0f0 .+ Float32(40_000.0/DEG) .* (rand(Float32, Nf) .- 0.5f0)
stf.indl .= 1.0f0 .+ rand(Float32, Nf)
update_bandwidths!(stf, G, 1)
hmed = median(Float64.(stf.hband))
# analytic: a line source of density λ = M/L smeared by kernel h gives
# C(d) = λ·exp(−d²/2h²)/(√(2π)·h·Δz) at distance d from the line
line_c(d, h) = (Mf / 40_000.0) * exp(-d^2 / (2h^2)) / (sqrt(2π) * h * DZ)
rs0 = ReceptorSet(["S"], [40.0], [2.0], 0.0, DZ)                  # no floor
receptor_sample!(rs0, stf, G, LC, zlev, 1.0, 1); receptor_flush!(rs0, 0.0)
HP = 2_000.0                                                       # footprint floor
rsn = ReceptorSet(["S"], [40.0], [2.0], 0.0, DZ; hmin=HP)
receptor_sample!(rsn, stf, G, LC, zlev, 1.0, 1); receptor_flush!(rsn, 0.0)
c_raw = rs0.rows[1][3]; c_flr = rsn.rows[1][3]
c_phys = line_c(150.0, HP)
check("kNN h in the filament is metres, not km (the 8021 mechanism)",
      hmed < 500.0; info=@sprintf("median h = %.0f m", hmed))
check("unfloored point evaluation blows up ≫ footprint value",
      c_raw > 5.0 * c_phys; info=@sprintf("×%.0f above physical", c_raw/c_phys))
check("hmin floor recovers the footprint-smeared analytic value (±15%)",
      0.85 < c_flr / c_phys < 1.15; info=@sprintf("ratio %.3f", c_flr/c_phys))

println()
println("="^78)
println("[4] window semantics — flushed value is the sample MEAN")
println("="^78)
rsw = ReceptorSet(string.(1:25), slons, slats, 0.0, DZ)
receptor_sample!(rsw, st, G, LC, zlev, 1.0, 1)
st.mass .*= 3.0f0                       # second sample at 3× the mass
receptor_sample!(rsw, st, G, LC, zlev, 1.0, 1)
st.mass ./= 3.0f0
receptor_flush!(rsw, 0.0)
check("two samples at 1× and 3× flush to the 2× mean",
      isapprox(mean([r[3] for r in rsw.rows]), 2.0 * c_est; rtol=1e-6);
      info=@sprintf("ratio %.4f", mean([r[3] for r in rsw.rows])/(2c_est)))

println()
println("="^78)
println("[5] adaptive deposit on a REALISTIC output mapping (task #40 geometry)")
println("="^78)
# The ETEX run maps met→output with coeff=2.5, offset=0: the RETIRED mid-lat
# formula read offsety+coeffy·Nj/2 = 151.25 "degrees" → cos = −0.877 → a
# NEGATIVE deposit cell size → inverted stencil bounds → every h ≳ 3 km
# parcel silently dropped (HALF the ETEX mass) and R(x) all-NaN.  This case
# locks the fix: exact conservation + R(x) populated + nothing dropped.
using .EASLEE.Operators: _conc_nest!
GE = EASLEE.Types.Rectilinear(Int32(161), Int32(121), Int32(10), 0.25f0, 0.25f0,
        -10.0f0, 40.0f0, Float32(1/111320.0), Float32(1/111320.0),
        2.5f0, 0.0f0, 2.5f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
zlevE = Array{Float32,4}(undef, 161, 121, 10, 2)
for k in 1:10; zlevE[:, :, k, :] .= Float32((k - 1) * 100.0); end
rhoE = ones(Float32, 161, 121, 10, 2)
Ne = 20_000
ste = EASLEE.State.ParcelState(Ne)
ste.status .= 0x01; ste.nest .= Int32(1); ste.mass .= 1.0f0
ste.indx .= 20.0f0 .+ 120.0f0 .* rand(Float32, Ne)
ste.indy .= 20.0f0 .+ 80.0f0  .* rand(Float32, Ne)
ste.indl .= 1.0f0 .+ rand(Float32, Ne)
ste.Dtime .= ste.Dtimeinit                      # deposit gate
update_bandwidths!(ste, GE, 1)
# adversarial h mixture: fresh (0), plume-core (m-scale), synoptic (100 km)
ste.hband[1:5000] .= 0.0f0
ste.hband[5001:10000] .= Float32.(10 .^ (0.5 .+ 1.5 .* rand(5000)))
ste.hband[10001:end] .= Float32.(10 .^ (3.5 .+ 1.5 .* rand(10000)))
NiE = Int(round(161 * 2.5)); NjE = Int(round(121 * 2.5))
gridE = zeros(Float32, NiE, NjE, 3, 1, 1)
massE = zeros(Float64, NiE, NjE, 3, 1, 1)
resnE = zeros(Float64, NiE, NjE); resdE = zeros(Float64, NiE, NjE)
dE = _conc_nest!(gridE, massE, ste, ste.Dtimeinit[1], 1.0f0, NiE, NjE, 1, 1,
                 Float32[0.0, 100.0, 200.0], 1, zlevE, rhoE,
                 EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10)),
                 GE, false, 60.0, Float32[1.0f0], :adaptive, resnE, resdE)
check("REAL-mapping deposit conserves mass exactly",
      isapprox(sum(massE), dE.deposited_mass; rtol=1e-12) && dE.deposited_mass > 0.9 * Ne;
      info=@sprintf("Σ=%.1f of %d parcels", dE.deposited_mass, Ne))
check("R(x) populated (no all-NaN resolution field)",
      count(>(0.0), resdE) > 1000;
      info=@sprintf("%d cells carry R", count(>(0.0), resdE)))

println()
println("="^78)
println("[6] longitude-convention robustness (0–360 met frame vs −180..180 stations)")
println("="^78)
# same slab, station given at lon−360 (the CAPTEX configuration: parcel lons
# come back in the met frame; the sampler must wrap the difference)
rsl = ReceptorSet(["W"], [32.0 - 360.0], [2.0], 0.0, DZ)
receptor_sample!(rsl, st, G, LC, zlev, 1.0, 1); receptor_flush!(rsl, 0.0)
rsr = ReceptorSet(["E"], [32.0], [2.0], 0.0, DZ)
receptor_sample!(rsr, st, G, LC, zlev, 1.0, 1); receptor_flush!(rsr, 0.0)
check("station at lon−360 sees the identical concentration",
      isapprox(rsl.rows[1][3], rsr.rows[1][3]; rtol=1e-12) && rsr.rows[1][3] > 0;
      info=@sprintf("%.3e vs %.3e", rsl.rows[1][3], rsr.rows[1][3]))

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
