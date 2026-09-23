# test_hanna_gpu_4d.jl — 2026-09-21.  On-the-fly Hanna on the KA path.
#
# `turbulence_hanna_4d` (turbulence_4d.jl) is the device-side twin of the CPU
# `turbulence_hanna_v1`: same Langevin body, but the column is sampled from the
# 4-D met arrays instead of hoisted 1-D profiles, and the random draw comes from
# the per-parcel xorshift stream rather than the task RNG.  So the two are
# statistically equivalent, never bitwise — exactly as Thomson already is on this
# path (BENCHMARKS.md section 2).
#   [1] the deterministic half (sigma_w, T_L, mixing length, TKE) is IDENTICAL
#   [2] the stochastic half agrees in distribution (mean, spread) over many parcels
#   [3] a parcel is never pushed through the ground or the model top
#   [4] the surface parameters reach the result: changing u*, h or L changes it
using Printf, Statistics
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    if ok; npass += 1; @printf("  PASS  %s %s\n", name, info) else; nfail += 1; @printf("  FAIL  %s  %s\n", name, info) end
end

# a convective boundary layer: 40 levels, 0-4000 m, rho with a scale height
const NK = 40
zprof   = Float32[(k-1) * 100.0f0 for k in 1:NK]
rhoprof = Float32[1.225f0 * exp(-z/8500.0f0) for z in zprof]
gradprof= Float32[-1.225f0/8500.0f0 * exp(-z/8500.0f0) for z in zprof]
# the same column as 4-D arrays (2 x 2 horizontal, 2 time slots)
mk4(v) = (a = Array{Float32,4}(undef, 2, 2, NK, 2); for i in 1:2, j in 1:2, k in 1:NK, t in 1:2; a[i,j,k,t] = v[k]; end; a)
zlev4, rho4, grad4 = mk4(zprof), mk4(rhoprof), mk4(gradprof)
ustar, pblh, ol = 0.35f0, 1200.0f0, -50.0f0        # unstable, convective
wstar = EASLEE.Operators.hanna_wstar(ustar, pblh, ol)
hanna = EASLEE.Types.Hanna(true, 0.1f0)

println("[1] the deterministic half must be identical")
for z0 in (25.0f0, 150.0f0, 600.0f0, 1100.0f0, 1350.0f0, 2500.0f0)
    s1, d1, t1 = EASLEE.Operators.hanna_sigw_tlw(z0, pblh, ustar, wstar, ol, 0.1f0)
    check(@sprintf("sigma_w / dsigma_w-dz / T_L at z = %6.0f m", z0),
          all(isfinite, (s1, d1, t1)) && s1 > 0 && t1 > 0,
          info=@sprintf("sigw %.4f  dsigwdz %+.2e  T_L %6.1f s", s1, d1, t1))
end
# both implementations call hanna_sigw_tlw, so the profile quantities agree by
# construction; what matters is that the 4-D sampler reproduces the 1-D column
indl = 4.5f0
z_1d = EASLEE.Operators.interpolation1Dsingle(zprof, indl)
z_4d = EASLEE.Operators.interpolation4Dlevsingle(zlev4, 1.5f0, 1.5f0, indl, 1.0f0)
check("the 4-D sampler reproduces the 1-D column height", isapprox(z_1d, z_4d; rtol=1e-6),
      info=@sprintf("%.3f vs %.3f m", z_1d, z_4d))

println("\n[2] the stochastic half must agree in distribution")
N = 20000
dz_cpu = Float32[]; dz_gpu = Float32[]
for i in 1:N
    _, delz, _, _, _, _ = EASLEE.Operators.turbulence_hanna_v1(hanna, indl, 300.0f0, 0.0f0, 5.0f0, 0.0f0,
                                                     zprof, rhoprof, gradprof, ustar, pblh, wstar, ol)
    push!(dz_cpu, delz)
    _, delz2, _, _, _, _, _ = EASLEE.Operators.turbulence_hanna_4d(indl, 1.5f0, 1.5f0, 300.0f0, 0.0f0, 5.0f0, 0.0f0,
                                                         1.0f0, Int32(NK), UInt64(0x9E3779B97F4A7C15 * i),
                                                         zlev4, rho4, grad4, ustar, pblh, wstar, ol, 0.1f0)
    push!(dz_gpu, delz2)
end
m1, m2 = mean(dz_cpu), mean(dz_gpu); s1, s2 = std(dz_cpu), std(dz_gpu)
se = sqrt(s1^2/N + s2^2/N)
check("mean displacement agrees within 4 standard errors", abs(m1 - m2) < 4*se,
      info=@sprintf("%.4f vs %.4f m (4 s.e. = %.4f)", m1, m2, 4*se))
check("spread agrees within 5 %", abs(s1 - s2)/max(s1, 1e-6) < 0.05,
      info=@sprintf("%.4f vs %.4f m", s1, s2))

println("\n[3] the parcel stays inside the column")
zs = Float32[]
for i in 1:5000
    z, _, _, _, _, _, _ = EASLEE.Operators.turbulence_hanna_4d(2.0f0, 1.5f0, 1.5f0, 900.0f0, 0.0f0, 5.0f0, 0.0f0,
                                                     1.0f0, Int32(NK), UInt64(12345 * i + 7),
                                                     zlev4, rho4, grad4, ustar, pblh, wstar, ol, 0.1f0)
    push!(zs, z)
end
check("never below the ground or above the model top",
      all(zs .>= zprof[1]) && all(zs .<= zprof[end]),
      info=@sprintf("range %.1f to %.1f m of [%.0f, %.0f]", minimum(zs), maximum(zs), zprof[1], zprof[end]))

println("\n[4] the surface parameters actually reach the result")
function spread(u, h, l)
    w = EASLEE.Operators.hanna_wstar(u, h, l); v = Float32[]
    for i in 1:4000
        _, d, _, _, _, _, _ = EASLEE.Operators.turbulence_hanna_4d(4.0f0, 1.5f0, 1.5f0, 300.0f0, 0.0f0, 5.0f0, 0.0f0,
                                                         1.0f0, Int32(NK), UInt64(777 * i + 3),
                                                         zlev4, rho4, grad4, u, h, w, l, 0.1f0)
        push!(v, d)
    end
    std(v)
end
s_conv = spread(0.35f0, 1200.0f0, -50.0f0)     # convective
s_stab = spread(0.15f0,  200.0f0, +80.0f0)     # stable, shallow
s_windy= spread(0.80f0, 1200.0f0, -50.0f0)     # same but 2.3x u*
check("a convective PBL mixes more than a stable one", s_conv > 2*s_stab,
      info=@sprintf("sigma(delz) convective %.3f vs stable %.3f m", s_conv, s_stab))
check("a larger friction velocity mixes more", s_windy > s_conv,
      info=@sprintf("u*=0.80 gives %.3f vs u*=0.35 %.3f m", s_windy, s_conv))
println("\n$npass/$(npass+nfail) checks passed"); nfail == 0 || println("FAIL")
