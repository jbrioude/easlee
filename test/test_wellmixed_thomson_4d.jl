# test_wellmixed_thomson_4d.jl — 2026-09-22.  Well-mixed criterion for the KA
# THOMSON body (`turbulence_v3_4d`) against the CPU body (`turbulence_v3`) on the
# same Hanna-filled TKE profile (tke = 0.5·σ_w², as hanna_fill_tke! writes it).
#
# Why this test exists: the KA body formed the Langevin drift term as a forward
# difference of the CLAMPED σ_w at the two bracketing levels, while the CPU body
# uses (wshare/σ_w)·∂tke/∂z — consistent with the value σ_w² = 2·wshare·tke.  The
# inconsistent form produces exactly the failure this criterion detects: spurious
# drift in the TKE-floor region above the PBL top (pile-up or void there), and it
# showed up as a +1.7 % mass offset between the KA and threads advances in an
# all-switches-off control.  After the 2026-09-22 alignment both bodies should
# pass, with the same PBL air-mass fraction.
# Layout contract 2026-10-22: θ and ρ are given class F (slot 1 = surface, slot k
# at zf(k) = ½(zlev[k-1] + zlev[k])); the CPU body gets them through the
# production half-level builders, the KA body forms the same half-level values
# from the 4-D arrays (no drhodz array).  Until then both were co-located with zlev.
using Printf, Statistics, Random
include("../easlee/EASLEE.jl"); using .EASLEE
using .EASLEE.Operators: turbulence_v3, turbulence_v3_4d, hanna_wstar, hanna_sigw_tlw,
                         interpolation1Dsingleirreg, interpolation4D_irregzsingle, VMixClosure,
                         full_to_half_levels!, half_level_gradient!
using .EASLEE.Types: Thomson
const NZ = 81; const ZTOP = 4000.0f0; const N = 100_000; const DT = 900.0f0; const NSTP = 96
height  = Float32.(range(0, ZTOP, length=NZ))
zf      = Float32[k == 1 ? height[1] : 0.5f0 * (height[k-1] + height[k]) for k in 1:NZ]
rhoF    = 1.2f0 .* exp.(-zf ./ 8500.0f0)                          # class F, as archived
pttF    = Float32[300.0f0 + 0.004f0 * z for z in zf]               # stable-ish θ profile for the length scale
rhograd = half_level_gradient!(similar(rhoF), rhoF, height, NZ)   # CPU profiles: the production builders
rho_h   = full_to_half_levels!(copy(rhoF), NZ)
ptt_h   = full_to_half_levels!(copy(pttF), NZ)
mk4(v) = (a = Array{Float32,4}(undef, 2, 2, NZ, 2); for i in 1:2, j in 1:2, k in 1:NZ, t in 1:2; a[i,j,k,t] = v[k]; end; a)
zlev4, rho4, ptt4 = mk4(height), mk4(rhoF), mk4(pttF)
vm = VMixClosure(false, NaN32, NaN32)                  # the Hanna-fill closure (wshare = 1)
rho    = 1.2f0 .* exp.(-height ./ 8500.0f0)            # analytic, for the air-mass sampling
layerm = [0.5f0*(rho[k]+rho[k+1]) * (height[k+1]-height[k]) for k in 1:NZ-1]
cm = cumsum(layerm); cm ./= cm[end]
expect_frac(pblh) = (1 - exp(-pblh/8500.0f0)) / (1 - exp(-ZTOP/8500.0f0))
function seed!(z, rng)
    for p in 1:N
        u = rand(rng, Float32); k = searchsortedfirst(cm, u)
        f = k == 1 ? u / cm[1] : (u - cm[k-1]) / (cm[k] - cm[k-1])
        z[p] = height[k] + f * (height[k+1] - height[k])
    end
end
function histo(z)
    nb = 20; hcount = zeros(Int, nb)
    for p in 1:N; hcount[clamp(ceil(Int, z[p] / (ZTOP/nb)), 1, nb)] += 1; end
    edges = range(0, ZTOP, length=nb+1)
    expectm = [1.2f0 * 8500.0f0 * (exp(-edges[b]/8500.0f0) - exp(-edges[b+1]/8500.0f0)) for b in 1:nb]
    expectn = expectm ./ sum(expectm) .* N
    maximum(abs.(hcount .- expectn) ./ expectn), 5 / sqrt(minimum(expectn)) + 0.03
end
function run_case(name, ustar, pblh, ol)
    wstar = hanna_wstar(ustar, pblh, ol)
    tke = Float32[0.5f0 * hanna_sigw_tlw(z, pblh, ustar, wstar, ol, 0.1f0)[1]^2 for z in height]
    tke4 = mk4(tke)
    res = Dict{String,Tuple{Float64,Float64,Float64}}()
    for body in ("CPU turbulence_v3", "KA turbulence_v3_4d")
        rng = MersenneTwister(7); z = Vector{Float32}(undef, N); seed!(z, rng); wp = randn(rng, Float32, N)
        Threads.@threads :static for p in 1:N
            w = wp[p]; ldt = 10.0f0
            if body == "CPU turbulence_v3"
                indl = interpolation1Dsingleirreg(height, z[p])
                for s in 1:NSTP
                    zp, _, w, ldt, _, _ = turbulence_v3(Thomson(true, false), indl, DT, ldt, 5.0f0, w,
                                                       height, ptt_h, tke, rho_h, rhograd; vm=vm)
                    indl = interpolation1Dsingleirreg(height, zp); z[p] = zp
                end
            else
                indl = interpolation4D_irregzsingle(zlev4, 1.5f0, 1.5f0, z[p], 1.0f0, Int32(1), Int32(NZ))
                st = UInt64(0x9E3779B97F4A7C15) * UInt64(p) + UInt64(3)
                for s in 1:NSTP
                    zp, _, w, ldt, _, _, st = turbulence_v3_4d(Thomson(true, false), indl, 1.5f0, 1.5f0, DT, ldt, 5.0f0, w,
                                                               1.0f0, Int32(NZ), st, zlev4, ptt4, tke4, rho4, vm)
                    indl = interpolation4D_irregzsingle(zlev4, 1.5f0, 1.5f0, zp, 1.0f0, Int32(1), Int32(NZ)); z[p] = zp
                end
            end
        end
        maxdev, tol = histo(z)
        res[body] = (maxdev, tol, count(p -> z[p] < pblh, 1:N) / N)
    end
    @printf("\n── %s  u*=%.2f h=%.0f L=%.0f  (air-mass fraction below h: %.1f%%) ──\n", name, ustar, pblh, ol, 100expect_frac(pblh))
    ok = true
    for body in ("CPU turbulence_v3", "KA turbulence_v3_4d")
        md, tol, fr = res[body]; pass = md < tol; ok &= pass
        @printf("  %-22s max deviation %5.1f%% (tol %.1f%%)  below h %.1f%%   %s\n", body, 100md, 100tol, 100fr, pass ? "PASS" : "FAIL")
    end
    return ok
end
ok = run_case("convective", 0.35f0, 1200.0f0, -50.0f0) & run_case("stable", 0.15f0, 300.0f0, 80.0f0) & run_case("neutral", 0.30f0, 800.0f0, Inf32)
println("\n", ok ? "ALL PASS" : "FAIL")
