# test_convection_adjoint.jl — the backward mass-flux convection kernel is the adjoint
# (mass-weighted transpose) of the forward one, and both keep a uniform ensemble uniform.
#   julia --project=easlee test/test_convection_adjoint.jl
using Printf, Random, Statistics
include(joinpath(@__DIR__, "..", "easlee", "src", "physics", "convection.jl"))

npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

# ── synthetic column: 30 levels, 300 m apart, exponential density, an updraught that
#    entrains between levels 3 and 10 and detrains up to 22, a weaker downdraught 12 → 2
nz = 30
z   = Float32[300.0 * (k - 1) for k in 1:nz]
rho = Float32[1.2 * exp(-z[k] / 8000.0) for k in 1:nz]
mu  = zeros(Float32, nz); md = zeros(Float32, nz)
for k in 3:22
    mu[k] = k <= 10 ? 0.08f0 * (k - 2) / 8 : 0.08f0 * (22 - k) / 12
end
for k in 2:12
    md[k] = -0.03f0 * (k <= 6 ? (k - 1) / 5 : (12 - k) / 6)     # archived downdraught is NEGATIVE
end
dt = 900.0f0
dm = zeros(Float32, nz); m1 = zeros(Float32, nz); S = zeros(Float32, nz)
mdmag = Float32[max(-md[k], 0f0) for k in 1:nz]
_build_conv_column!(dm, m1, S, z, rho, mu, mdmag, dt, nz)
eu = [dm[k] > 0 ? max(_uflux(mu, k) - _uflux(mu, k - 1), 0f0) * dt / dm[k] : 0f0 for k in 1:nz-1]
@printf("column: entrained fraction per step max %.3f at layer %d; Σdm %.1f kg/m², Σm1 %.1f\n", maximum(eu), argmax(eu), sum(dm), sum(m1))

# ── Monte-Carlo transition matrices
N = 40_000
rng = MersenneTwister(7)
rfn() = rand(rng, Float32)
P = zeros(nz - 1, nz - 1); Q = zeros(nz - 1, nz - 1)
for k in 1:nz-1, _ in 1:N
    P[k, _convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, false)] += 1 / N
    Q[k, _convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, true)]  += 1 / N
end
# adjoint identity: Q(l→k)·dm[l] == P(k→l)·dm[k]   (mass moved from k to l forward = mass traced back from l to k)
A = [P[k, l] * dm[k] for k in 1:nz-1, l in 1:nz-1]        # forward mass flow k→l
B = [Q[l, k] * dm[l] for k in 1:nz-1, l in 1:nz-1]        # backward mass flow l→k, same indexing
offdiag = [k != l for k in 1:nz-1, l in 1:nz-1]
moved = sum(A[offdiag])
err = sum(abs.(A[offdiag] .- B[offdiag])) / moved
check("Σ|A−B| over off-diagonal mass flow / Σ moved mass < 0.06 (MC noise)", err < 0.06, @sprintf("%.4f (moved %.2f of %.1f kg/m²)", err, moved, sum(dm)))
# the largest single flows agree individually
idx = sortperm(vec(A[offdiag]); rev=true)[1:8]
Ao = vec(A[offdiag])[idx]; Bo = vec(B[offdiag])[idx]
check("8 largest flows agree to 10 %", all(abs.(Ao .- Bo) ./ Ao .< 0.10), join([@sprintf("%.3f/%.3f", a, b) for (a, b) in zip(Ao, Bo)], " "))
# the transpose identity would FAIL for the forward kernel applied backward (the old code):
Bold = [P[l, k] * dm[l] for k in 1:nz-1, l in 1:nz-1]
errold = sum(abs.(A[offdiag] .- Bold[offdiag])) / moved
check("(control) forward kernel used backward is NOT adjoint: error > 0.5", errold > 0.5, @sprintf("%.3f", errold))
# well-mixedness: a uniform-in-air-mass ensemble stays uniform under both kernels
M = 400_000
w = dm[1:nz-1] ./ sum(dm[1:nz-1])
starts = [rand(rng) for _ in 1:M]; cw = cumsum(w)
outF = zeros(nz - 1); outB = zeros(nz - 1)
for s in starts
    k = clamp(searchsortedfirst(cw, s), 1, nz - 1)
    outF[_convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, false)] += 1 / M
    outB[_convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, true)]  += 1 / M
end
devF = maximum(abs.(outF ./ w .- 1)); devB = maximum(abs.(outB ./ w .- 1))
check("uniform ensemble stays uniform, forward (max |χ−1| < 0.05)", devF < 0.05, @sprintf("%.3f", devF))
check("uniform ensemble stays uniform, backward (max |χ−1| < 0.05)", devB < 0.05, @sprintf("%.3f", devB))
println("\n$(npass) passed, $(nfail) failed"); nfail == 0 || exit(1)
