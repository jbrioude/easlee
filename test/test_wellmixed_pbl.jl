# test_wellmixed_pbl.jl — well-mixed criterion for the Hanna turbulence scheme
# with the blended PBL entrainment zone (option b, on-the-fly).
#
# Parcels are initialised uniformly in AIR MASS over a 0–4 km column spanning
# the boundary layer and the free troposphere, then advanced 24 h with
# turbulence only.  A correct scheme keeps the distribution uniform in air
# mass — in particular NO pile-up or void at the PBL top, the classic failure
# mode of a hard σ_w discontinuity at h.  Also reports the PBL↔FT exchange
# fraction (must be finite — closed-box bouncing would give 0 — but modest).
#
# Layout contract 2026-10-22: the synthetic column is given as the met holds it —
# zlev on the half levels, rho class F (slot 1 = surface, slot k at
# zf(k) = ½(zlev[k-1] + zlev[k])) — and the ρ / ∂ρ/∂z profiles the turbulence
# reads are built by the production half-level builders (met/interp.jl).  Until
# then rho was given co-located with zlev (rho[k] at height[k]).
#
# Usage:  julia --threads=8 --project=../easlee test_wellmixed_pbl.jl

using Printf, Statistics, Random
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.Operators: turbulence_hanna_v1, hanna_wstar, interpolation1Dsingleirreg,
                         full_to_half_levels!, half_level_gradient!
using .EASLEE.Types: Hanna

const NZ   = 81                 # 50 m spacing, 0–4000 m
const ZTOP = 4000.0f0
const N    = 100_000
const DT   = 900.0f0
const NSTP = 96                 # 24 h

height  = Float32.(range(0, ZTOP, length=NZ))                 # zlev (half levels)
rho_of(z) = 1.2f0 * exp(-z / 8500.0f0)
zf      = Float32[k == 1 ? height[1] : 0.5f0 * (height[k-1] + height[k]) for k in 1:NZ]
rhoF    = rho_of.(zf)                                         # class F, as archived
rhograd = half_level_gradient!(similar(rhoF), rhoF, height, NZ)
rho     = full_to_half_levels!(copy(rhoF), NZ)               # ρ at the half levels

# Uniform-in-air-mass sampling over the column (analytic density)
rhoz   = rho_of.(height)
layerm = [0.5f0*(rhoz[k]+rhoz[k+1]) * (height[k+1]-height[k]) for k in 1:NZ-1]
cm = cumsum(layerm); cm ./= cm[end]

function run_case(name, ustar, pblh, ol)
    wstar = hanna_wstar(ustar, pblh, ol)
    hanna = Hanna(true)
    rng = MersenneTwister(7)

    z  = Vector{Float32}(undef, N)
    wp = randn(rng, Float32, N)
    for p in 1:N
        u = rand(rng, Float32)
        k = searchsortedfirst(cm, u)
        f = k == 1 ? u / cm[1] : (u - cm[k-1]) / (cm[k] - cm[k-1])
        z[p] = height[k] + f * (height[k+1] - height[k])
    end
    below0 = count(p -> z[p] < pblh, 1:N)

    Threads.@threads :static for p in 1:N
        indl = interpolation1Dsingleirreg(height, z[p])
        w    = wp[p]
        ldt  = 10.0f0
        for s in 1:NSTP
            zp, _, w, ldt, _, _ = turbulence_hanna_v1(hanna, indl, DT, ldt, 5.0f0, w,
                                                      height, rho, rhograd,
                                                      ustar, pblh, wstar, ol)
            indl = interpolation1Dsingleirreg(height, zp)
            z[p] = zp
        end
        wp[p] = w
    end

    # Histogram in 200 m bins, compare against the air-mass distribution
    nb = 20
    edges = range(0, ZTOP, length=nb+1)
    hcount = zeros(Int, nb)
    for p in 1:N
        b = clamp(ceil(Int, z[p] / (ZTOP/nb)), 1, nb)
        hcount[b] += 1
    end
    expectm = [begin
        zb1, zb2 = edges[b], edges[b+1]
        1.2f0 * 8500.0f0 * (exp(-zb1/8500.0f0) - exp(-zb2/8500.0f0))
    end for b in 1:nb]
    expectn = expectm ./ sum(expectm) .* N

    devs = abs.(hcount .- expectn) ./ expectn
    below1 = count(p -> z[p] < pblh, 1:N)
    @printf("\n── %s  (u*=%.2f  h=%.0f  L=%.0f  w*=%.2f) ──\n", name, ustar, pblh, ol, wstar)
    println("   z bin [m]      expected   final    dev")
    for b in 1:nb
        flag = edges[b] < pblh < edges[b+1] ? "  ← PBL top" : ""
        @printf("  %5.0f–%-5.0f  %8.0f  %7d  %5.1f%%%s\n",
                edges[b], edges[b+1], expectn[b], hcount[b], 100devs[b], flag)
    end
    maxdev = maximum(devs)
    # σ of counting noise ≈ sqrt(N_bin); 5σ tolerance plus 3% systematic
    tol = 5 / sqrt(minimum(expectn)) + 0.03
    ok  = maxdev < tol
    @printf("  max deviation %.1f%% (tolerance %.1f%%)  —  %s\n",
            100maxdev, 100tol, ok ? "PASS" : "FAIL")
    @printf("  parcels below h: %.1f%% → %.1f%% (air-mass fraction %.1f%%)\n",
            100below0/N, 100below1/N,
            100 * (1 - exp(-pblh/8500.0f0)) / (1 - exp(-ZTOP/8500.0f0)))
    return ok
end

println("=== Well-mixed criterion across the PBL entrainment zone ===")
ok1 = run_case("convective", 0.45f0, 1200.0f0, -60.0f0)
ok2 = run_case("stable",     0.25f0,  400.0f0, 150.0f0)
ok3 = run_case("neutral",    0.40f0,  900.0f0, 1.0f6)

println()
println(ok1 && ok2 && ok3 ? "=== ALL PASS ===" : "=== FAILURES — see above ===")
exit(ok1 && ok2 && ok3 ? 0 : 1)
