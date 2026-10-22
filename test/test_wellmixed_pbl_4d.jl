# test_wellmixed_pbl_4d.jl — 2026-09-22.  The well-mixed criterion applied to the
# DEVICE-SIDE Hanna (`turbulence_hanna_4d`), the physics test that matters for a
# Langevin scheme: parcels started uniform in air mass over a 0-4 km column must
# stay uniform after 24 h of turbulence only — no pile-up or void at the PBL top.
# Same column, same cases and the same tolerance as test_wellmixed_pbl.jl (CPU),
# with the column presented as 4-D met arrays and the per-parcel xorshift RNG.
# Layout contract 2026-10-22: rho4 is class F (slot 1 = surface, slot k at
# zf(k) = ½(zlev[k-1] + zlev[k])), as the readers store it; the kernel forms ρ and
# ∂ρ/∂z at the half levels itself (no drhodz array).  Until then rho4 was
# co-located with zlev and a precomputed gradient array was passed.
using Printf, Statistics, Random
include("../easlee/EASLEE.jl"); using .EASLEE
using .EASLEE.Operators: turbulence_hanna_4d, hanna_wstar, interpolation4D_irregzsingle, interpolation4Dlevsingle
const NZ = 81; const ZTOP = 4000.0f0; const N = 100_000; const DT = 900.0f0; const NSTP = 96
height  = Float32.(range(0, ZTOP, length=NZ))
rho     = 1.2f0 .* exp.(-height ./ 8500.0f0)                 # analytic, at the half levels
zf      = Float32[k == 1 ? height[1] : 0.5f0 * (height[k-1] + height[k]) for k in 1:NZ]
rhoF    = 1.2f0 .* exp.(-zf ./ 8500.0f0)                     # class F, as archived
mk4(v) = (a = Array{Float32,4}(undef, 2, 2, NZ, 2); for i in 1:2, j in 1:2, k in 1:NZ, t in 1:2; a[i,j,k,t] = v[k]; end; a)
zlev4, rho4 = mk4(height), mk4(rhoF)
layerm = [0.5f0*(rho[k]+rho[k+1]) * (height[k+1]-height[k]) for k in 1:NZ-1]
cm = cumsum(layerm); cm ./= cm[end]
function run_case(name, ustar, pblh, ol)
    wstar = hanna_wstar(ustar, pblh, ol)
    rng = MersenneTwister(7)
    z = Vector{Float32}(undef, N); wp = randn(rng, Float32, N)
    for p in 1:N
        u = rand(rng, Float32); k = searchsortedfirst(cm, u)
        f = k == 1 ? u / cm[1] : (u - cm[k-1]) / (cm[k] - cm[k-1])
        z[p] = height[k] + f * (height[k+1] - height[k])
    end
    below0 = count(p -> z[p] < pblh, 1:N)
    Threads.@threads :static for p in 1:N
        indl = interpolation4D_irregzsingle(zlev4, 1.5f0, 1.5f0, z[p], 1.0f0, Int32(1), Int32(NZ))
        w = wp[p]; ldt = 10.0f0; st = UInt64(0x9E3779B97F4A7C15) * UInt64(p) + UInt64(11)
        for s in 1:NSTP
            zp, _, w, ldt, _, _, st = turbulence_hanna_4d(indl, 1.5f0, 1.5f0, DT, ldt, 5.0f0, w,
                                                          1.0f0, Int32(NZ), st, zlev4, rho4,
                                                          ustar, pblh, wstar, ol, 0.1f0)
            indl = interpolation4D_irregzsingle(zlev4, 1.5f0, 1.5f0, zp, 1.0f0, Int32(1), Int32(NZ))
            z[p] = zp
        end
        wp[p] = w
    end
    nb = 20; edges = range(0, ZTOP, length=nb+1); hcount = zeros(Int, nb)
    for p in 1:N; hcount[clamp(ceil(Int, z[p] / (ZTOP/nb)), 1, nb)] += 1; end
    expectm = [1.2f0 * 8500.0f0 * (exp(-edges[b]/8500.0f0) - exp(-edges[b+1]/8500.0f0)) for b in 1:nb]
    expectn = expectm ./ sum(expectm) .* N
    devs = abs.(hcount .- expectn) ./ expectn
    below1 = count(p -> z[p] < pblh, 1:N)
    @printf("\n── %s (device-side Hanna)  u*=%.2f h=%.0f L=%.0f w*=%.2f ──\n", name, ustar, pblh, ol, wstar)
    for b in 1:nb
        flag = edges[b] < pblh < edges[b+1] ? "  ← PBL top" : ""
        @printf("  %5.0f–%-5.0f  %8.0f  %7d  %5.1f%%%s\n", edges[b], edges[b+1], expectn[b], hcount[b], 100devs[b], flag)
    end
    maxdev = maximum(devs); tol = 5 / sqrt(minimum(expectn)) + 0.03; ok = maxdev < tol
    @printf("  max deviation %.1f%% (tolerance %.1f%%)  —  %s\n", 100maxdev, 100tol, ok ? "PASS" : "FAIL")
    @printf("  parcels below h: %.1f%% → %.1f%% (air-mass fraction %.1f%%)\n", 100below0/N, 100below1/N,
            100 * (1 - exp(-pblh/8500.0f0)) / (1 - exp(-ZTOP/8500.0f0)))
    return ok
end
ok1 = run_case("convective", 0.35f0, 1200.0f0, -50.0f0)
ok2 = run_case("stable",     0.15f0,  300.0f0,  80.0f0)
ok3 = run_case("neutral",    0.30f0,  800.0f0, Inf32)
println("\n", (ok1 && ok2 && ok3) ? "ALL PASS" : "FAIL")
