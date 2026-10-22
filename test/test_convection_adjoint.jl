# test_convection_adjoint.jl — the backward mass-flux convection kernel is the adjoint
# (mass-weighted transpose) of the forward one, and both keep a uniform ensemble uniform.
# Columns follow the stacked layout of the readers: zlev on the half levels (zlev[1] = 0),
# rho/umf/dmf on the full levels, slot k ≥ 2 between zlev[k-1] and zlev[k], slot 1 the
# surface (fluxes 0).  The scheme's layers are the full levels k = 2 … nz.
#   julia --project=easlee test/test_convection_adjoint.jl
using Printf, Random, Statistics
include(joinpath(@__DIR__, "..", "easlee", "src", "met", "interp.jl"))     # get_vertical_profile!
include(joinpath(@__DIR__, "..", "easlee", "src", "gpu", "rng.jl"))        # parcel xorshift streams
include(joinpath(@__DIR__, "..", "easlee", "src", "physics", "convection.jl"))

npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

# ── synthetic column: 30 half levels 300 m apart, exponential density on the full levels,
#    an updraught that entrains between full levels 3 and 10 and detrains up to 22, a
#    weaker downdraught 12 → 2
nz = 30
z   = Float32[300.0 * (k - 1) for k in 1:nz]                                  # half levels
zf  = Float32[k == 1 ? 0.0 : 0.5 * (z[k-1] + z[k]) for k in 1:nz]            # full levels
rho = Float32[1.2 * exp(-zf[k] / 8000.0) for k in 1:nz]
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
L = 2:nz                                                                      # the layers
eu = [dm[k] > 0 ? max(_uflux(mu, k) - _uflux(mu, k - 1), 0f0) * dt / dm[k] : 0f0 for k in L]
@printf("column: entrained fraction per step max %.3f at layer %d; Σdm %.1f kg/m², Σm1 %.1f\n", maximum(eu), L[argmax(eu)], sum(dm), sum(m1))
check("layer air mass dm[k] = ρ[k]·(zlev[k] − zlev[k−1]), dm[1] = 0",
      dm[1] == 0f0 && all(abs(dm[k] - rho[k] * (z[k] - z[k-1])) <= 1f-5 * dm[k] for k in L),
      @sprintf("dm[2] %.2f  dm[%d] %.2f kg/m²", dm[2], nz, dm[nz]))

# ── Monte-Carlo transition matrices (row/column a ↔ layer L[a] = a + 1)
N = 40_000
rng = MersenneTwister(7)
rfn() = rand(rng, Float32)
nl = length(L)
P = zeros(nl, nl); Q = zeros(nl, nl)
for (a, k) in enumerate(L), _ in 1:N
    P[a, _convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, false) - 1] += 1 / N
    Q[a, _convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, true) - 1]  += 1 / N
end
dmL = dm[L]
# adjoint identity: Q(l→k)·dm[l] == P(k→l)·dm[k]   (mass moved from k to l forward = mass traced back from l to k)
A = [P[k, l] * dmL[k] for k in 1:nl, l in 1:nl]        # forward mass flow k→l
B = [Q[l, k] * dmL[l] for k in 1:nl, l in 1:nl]        # backward mass flow l→k, same indexing
offdiag = [k != l for k in 1:nl, l in 1:nl]
moved = sum(A[offdiag])
err = sum(abs.(A[offdiag] .- B[offdiag])) / moved
check("Σ|A−B| over off-diagonal mass flow / Σ moved mass < 0.06 (MC noise)", err < 0.06, @sprintf("%.4f (moved %.2f of %.1f kg/m²)", err, moved, sum(dm)))
# the largest single flows agree individually
idx = sortperm(vec(A[offdiag]); rev=true)[1:8]
Ao = vec(A[offdiag])[idx]; Bo = vec(B[offdiag])[idx]
check("8 largest flows agree to 10 %", all(abs.(Ao .- Bo) ./ Ao .< 0.10), join([@sprintf("%.3f/%.3f", a, b) for (a, b) in zip(Ao, Bo)], " "))
# the transpose identity would FAIL for the forward kernel applied backward (the old code):
Bold = [P[l, k] * dmL[l] for k in 1:nl, l in 1:nl]
errold = sum(abs.(A[offdiag] .- Bold[offdiag])) / moved
check("(control) forward kernel used backward is NOT adjoint: error > 0.5", errold > 0.5, @sprintf("%.3f", errold))
# well-mixedness: a uniform-in-air-mass ensemble stays uniform under both kernels
M = 400_000
w = dmL ./ sum(dmL)
starts = [rand(rng) for _ in 1:M]; cw = cumsum(w)
outF = zeros(nl); outB = zeros(nl)
for s in starts
    k = L[clamp(searchsortedfirst(cw, s), 1, nl)]
    outF[_convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, false) - 1] += 1 / M
    outB[_convect_column(rfn, k, dm, m1, S, mu, mdmag, dt, nz, true) - 1]  += 1 / M
end
devF = maximum(abs.(outF ./ w .- 1)); devB = maximum(abs.(outB ./ w .- 1))
check("uniform ensemble stays uniform, forward (max |χ−1| < 0.05)", devF < 0.05, @sprintf("%.3f", devF))
check("uniform ensemble stays uniform, backward (max |χ−1| < 0.05)", devB < 0.05, @sprintf("%.3f", devB))

# ── End to end through _convect_nest!: parcel → layer assignment (indl ∈ [k−1, k) is
#    layer k), landing heights and the flux placement on a stretched, ERA5-like column
#    (half-level spacing 20 m at the ground → 700 m aloft).  Fluxes on the full levels:
#    the updraught grows linearly from the ground to 3 km and stops ABRUPTLY above the
#    full level ktop, so no parcel may cross zlev[ktop+1] upward (the interface flux
#    there is ½(umf[ktop+1] + umf[ktop+2]) = 0); applied one layer too high, as the
#    collocated reading did, risers would cross it.
nz2 = 60
dz2 = Float32[20 + 680 * ((k - 1) / (nz2 - 2))^1.3 for k in 1:nz2-1]
z2  = Float32[0; cumsum(dz2)]
zf2 = Float32[k == 1 ? 0 : 0.5 * (z2[k-1] + z2[k]) for k in 1:nz2]
rho2 = Float32[1.2 * exp(-zf2[k] / 7500) for k in 1:nz2]
ktop = findlast(k -> zf2[k] < 12000, 1:nz2)
mu2 = Float32[(k == 1 || k > ktop) ? 0 : 0.06 * min(zf2[k] / 3000, 1) for k in 1:nz2]
md2 = Float32[(k == 1 || zf2[k] > 4000) ? 0 : -0.02 * sin(pi * zf2[k] / 4000) for k in 1:nz2]
@printf("\nERA5-like column: %d half levels, top %.0f m; updraught top: full level %d (zf %.0f m), zlev[%d] = %.0f m\n",
        nz2, z2[end], ktop, zf2[ktop], ktop + 1, z2[ktop+1])
col(v) = (a = zeros(Float32, 3, 3, nz2, 2); for i in 1:3, j in 1:3, t in 1:2; a[i, j, :, t] .= v; end; a)
nest = (metfields3D = (zlev = col(z2), rho = col(rho2)), gridconfig = (Ni = 3, Nj = 3))
cf   = (umf = col(mu2), dmf = col(md2))
dm2 = Float64[k == 1 ? 0 : rho2[k] * (z2[k] - z2[k-1]) for k in 1:nz2]     # the air mass, independent of the scheme
L2 = 2:nz2; w2 = dm2[L2] ./ sum(dm2[L2]); cw2 = cumsum(w2)
NP = 1_000_000
seed_state(np) = (indx = fill(2.0f0, np), indy = fill(2.0f0, np), indl = zeros(Float32, np),
                  rng_state = [rng_seed(i) for i in 1:np])
function seed_uniform!(st, rng)        # uniform in air mass: layer ∝ dm, uniform in height inside it
    for i in eachindex(st.indl)
        k = L2[clamp(searchsortedfirst(cw2, rand(rng)), 1, length(L2))]
        st.indl[i] = Float32(k - 1) + rand(rng, Float32)
    end
end
zof(st, i) = _z_at(z2, st.indl[i], nz2)
layer_of(st, i) = clamp(Base.unsafe_trunc(Int, st.indl[i]) + 1, 2, nz2)
chi_dev(st) = begin                     # max |χ−1| in units of the Poisson σ, over the layers
    cnt = zeros(Int, nz2); for i in eachindex(st.indl); cnt[layer_of(st, i)] += 1; end
    maximum(abs(cnt[k] - NP * w2[k-1]) / sqrt(NP * w2[k-1]) for k in L2)
end
rng2 = MersenneTwister(3)
st = seed_state(NP); seed_uniform!(st, rng2)
z0 = [zof(st, i) for i in 1:NP]
dts = 300f0
_convect_nest!(st, nest, cf, dts, 1.0f0, 1:NP)
z1 = [zof(st, i) for i in 1:NP]
mpp = sum(dm2[L2]) / NP
Fif(k) = _uflux(mu2, k)
up(k) = count(i -> z0[i] < z2[k] <= z1[i], 1:NP) * mpp / dts       # upward crossings of zlev[k]
nabove = count(i -> z0[i] < z2[ktop+1] && z1[i] > z2[ktop+1], 1:NP)
check("no parcel crosses zlev[ktop+1] upward (interface flux 0 there)", nabove == 0,
      @sprintf("%d crossings; through zlev[ktop] %.4f vs ½·umf[ktop] %.4f kg m⁻² s⁻¹", nabove, up(ktop), Fif(ktop)))
act = [k for k in 2:ktop if Fif(k) > 0.01f0]
Fm  = [up(k) for k in act]; Fa = [Fif(k) for k in act]; Fs = [Fif(k - 1) for k in act]
e_here  = sum(abs.(Fm .- Fa)) / sum(Fa)
e_shift = sum(abs.(Fm .- Fs)) / sum(Fa)
check("parcel up-flux through zlev[k] = ½(umf[k]+umf[k+1]) (L1 < 0.03)", e_here < 0.03,
      @sprintf("L1 %.3f  (vs the profile one layer lower: %.3f); Σ ratio %.3f", e_here, e_shift, sum(Fm) / sum(Fa)))
nground = count(i -> z0[i] < z2[2] && z1[i] > z2[2], 1:NP)
check("the ground layer [0, zlev[2]] feeds the updraught", nground > 0,
      @sprintf("%d risers, expected %.0f", nground, max(Fif(2) - Fif(1), 0f0) * dts / mpp))
check("landing heights inside [0, zlev[nz-1]], finite", all(i -> isfinite(z1[i]) && 0 <= z1[i] <= z2[nz2-1], 1:NP), "")
for _ in 1:19; _convect_nest!(st, nest, cf, dts, 1.0f0, 1:NP); end
dF = chi_dev(st)
st = seed_state(NP); seed_uniform!(st, rng2)
for _ in 1:20; _convect_nest!(st, nest, cf, dts, 1.0f0, 1:NP; backward = true); end
dB = chi_dev(st)
check("uniform ensemble stays uniform over 20 steps, forward (max dev < 5σ)", dF < 5, @sprintf("%.2f σ", dF))
check("uniform ensemble stays uniform over 20 steps, backward (max dev < 5σ)", dB < 5, @sprintf("%.2f σ", dB))
println("\n$(npass) passed, $(nfail) failed"); nfail == 0 || exit(1)
