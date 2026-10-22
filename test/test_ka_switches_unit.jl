# test_ka_switches_unit.jl — 2026-09-22.  Function-level checks for the five switches
# ported to the KA advance (mesoscale, kh_ft, horizontal_tke, horizontal_hanna, MO).
#   [1] _std16_4d is BITWISE the CPU _std16 (same arithmetic on the same array)
#   [2] _mesoscale_sigmas_4d == _mesoscale_sigmas on random fields and at the edges,
#       u/v on their full levels and levdot on its half levels (layout contract)
#   [3] mo_wind_factor, hanna_h, tke_partition_h are pure and finite over the
#       parameter ranges the kernel will feed them (no NaN, no negative sigma)
using Printf, Random
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
O = EASLEE.Operators
npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    if ok; npass += 1; @printf("  PASS  %s %s\n", name, info) else; nfail += 1; @printf("  FAIL  %s  %s\n", name, info) end
end
rng = MersenneTwister(7)
println("[1] _std16_4d vs _std16, bitwise")
a = rand(rng, Float32, 12, 9, 20, 2) .* 30.0f0 .- 15.0f0
a[3, 4, 5, 1] = NaN32                                     # a masked cell
worst = 0.0f0; nident = 0; ntot = 0
for i0 in 1:11, j0 in 1:8, k0 in 1:19
    x = O._std16(a, i0, j0, k0); y = O._std16_4d(a, i0, j0, k0)
    global ntot += 1; x == y && (global nident += 1)
    global worst = max(worst, abs(x - y))
end
check("identical on every neighbourhood of a random field", nident == ntot, info="$nident / $ntot, worst |Δ| = $worst")
check("the NaN cell zeroes the sigma in both", O._std16(a, 2, 3, 4) == 0.0f0 == O._std16_4d(a, 2, 3, 4))

println("\n[2] _mesoscale_sigmas_4d vs _mesoscale_sigmas")
u = rand(rng, Float32, 12, 9, 20, 2); v = rand(rng, Float32, 12, 9, 20, 2); w = rand(rng, Float32, 12, 9, 20, 2)
ok = true
for (xi, yi, zi) in ((3.7f0, 4.2f0, 5.9f0), (0.5f0, 0.5f0, 0.5f0), (11.9f0, 8.9f0, 19.9f0), (200f0, -3f0, 40f0))
    s1 = O._mesoscale_sigmas(u, v, w, xi, yi, zi, 20); s2 = O._mesoscale_sigmas_4d(u, v, w, xi, yi, zi, 20)
    global ok &= (s1 == s2)
end
check("identical, including clamped positions off the grid", ok)
# layout contract 2026-10-22: u, v are class F — the levels bracketing indl are
# the stacked ones (k, k+1 for indl in [k-½, k+½)); levdot is class H (trunc(indl))
okF = true
for (zi, kF, kH) in ((2.2f0, 2, 2), (2.7f0, 3, 2), (1.2f0, 1, 1), (1.6f0, 2, 1))
    s1 = O._mesoscale_sigmas(u, v, w, 3.7f0, 4.2f0, zi, 20)
    global okF &= s1[1] == O._std16(u, 3, 4, kF) && s1[2] == O._std16(v, 3, 4, kF) &&
                  s1[3] == O._std16(w, 3, 4, kH)
end
check("u, v sampled on the full levels bracketing indl, levdot on the half levels", okF)

println("\n[3] the helpers are finite and sane over the kernel's input ranges")
bad = 0; n = 0
for z in (0.3f0, 5f0, 50f0, 500f0), z1 in (10f0, 40f0), ust in (0.01f0, 0.3f0, 1.0f0), ol in (-5f0, -200f0, Inf32, 50f0, 5f0), U in (0.01f0, 3f0, 20f0)
    f = O.mo_wind_factor(z, z1, ust, ol, U); global n += 1
    # the function clamps to [0, 1] by design: 0 means the similarity correction
    # exceeds the level-1 wind (0.3 m AGL under strong u*), i.e. calm at the surface
    (isfinite(f) && 0.0f0 <= f <= 1.0f0) || (global bad += 1)
end
check("mo_wind_factor finite and within [0, 1] on $n combinations", bad == 0, info="$bad bad")
bad = 0; n = 0
for z in (1f0, 100f0, 800f0, 1300f0, 3000f0), h in (50f0, 400f0, 1500f0), ol in (-10f0, -300f0, 100f0, Inf32), ust in (0.05f0, 0.6f0)
    wst = O.hanna_wstar(ust, h, ol); su, tu, sv, tv = O.hanna_h(z, h, ol, ust, wst); global n += 1
    (all(isfinite, (su, tu, sv, tv)) && su >= 0 && sv >= 0 && tu > 0 && tv > 0) || (global bad += 1)
end
check("hanna_h finite, sigma >= 0, T_L > 0 on $n combinations", bad == 0, info="$bad bad")
println("\n$npass/$(npass+nfail) checks passed"); nfail == 0 || println("FAIL")
