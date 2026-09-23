#!/usr/bin/env julia
# test_vertical_modes.jl — SYNTHETIC-TRUTH tests of the vertical transport chains
# (the descent-anomaly investigation, tasks #27/#28):
#   [1] isentropic pinning (hybrid_vertical! at w = 1): a parcel given a fake
#       kinematic displacement must be pulled back EXACTLY onto its prognostic
#       isentrope (the kinematic vertical is erased, not stacked),
#   [2] diabatic drift: with a known θ̇ = (θ/T)·Q the parcel descends the θ
#       ladder by exactly θ̇·dt per step; θ_p stays prognostic (no grid re-seed),
#   [3] no drivers → no motion: θ̇ = 0 and no kinematic offset leaves indl and
#       θ_p unchanged (machine precision),
#   [4] blend gating: a parcel in the deep troposphere (σ > σ_r) is untouched
#       and its θ_p tracks the grid.
#
#   julia --project=easlee test/test_vertical_modes.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: HybridVert, hybrid_vertical!
using .EASLEE.Operators: IsentropicFields

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic stratospheric column ───────────────────────────────────────────
# 40 levels, z = (k−1)·500 m (0–19.5 km … extend: use 60 levels to 29.5 km),
# θ(k) linear in z above the tropopause: θ = 380 + 35·(z−16 km)/1 km above 16 km,
# and p from a scale-height profile so σ = p/ps puts the strat at w = 1.
const NZ = 60
dims = (6, 6, NZ, 2)
zlev  = zeros(Float32, dims); theta = zeros(Float32, dims)
puv   = zeros(Float32, dims); Tarr  = zeros(Float32, dims)
rho   = zeros(Float32, dims)
for k in 1:NZ
    z = (k - 1) * 500.0
    zlev[:, :, k, :]  .= Float32(z)
    p = 101325.0 * exp(-z / 7000.0)
    puv[:, :, k, :]   .= Float32(p)
    T = z < 16000.0 ? 288.0 - 0.0065 * z : 200.0
    Tarr[:, :, k, :]  .= Float32(T)
    rho[:, :, k, :]   .= Float32(p / (287.05 * T))
    θ = T * (1.0e5 / p)^0.2857
    theta[:, :, k, :] .= Float32(θ)
end
# uniform LW cooling Q [K/s] in tt1; no tt2
const Q_KDAY = -1.0
tt1 = fill(Float32(Q_KDAY / 86400.0), dims)
m3 = (zlev = zlev, theta = theta, puv = puv, T = Tarr, rho = rho,
      Ux = zeros(Float32, dims), Vy = zeros(Float32, dims))
LC = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
GC = (Ni = Int32(6), Nj = Int32(6), Uedge = true, Vedge = true, Wedge = true)
RC = (currentdate = DateTime(2025, 5, 22, 12), direction = Int8(1))
fctx = (nests = [(metfields3D = m3, levelconfig = LC, gridconfig = GC)],
        runconfig = RC)

# θ at a continuous index (mirror of the module's interpolation)
θat(indl) = begin
    k0 = clamp(floor(Int, indl), 1, NZ - 1)
    theta[3, 3, k0, 1] + (indl - k0) * (theta[3, 3, k0+1, 1] - theta[3, 3, k0, 1])
end

function mkhyb(; diabatic=true)
    HybridVert(Union{IsentropicFields,Nothing}[IsentropicFields{typeof(tt1)}(tt1, nothing)],
               nothing,                        # no PT winds in the fixture
               zeros(Float32, 50), 0.3f0, 0.1f0, diabatic)
end
function mkstate(indl0)
    st = EASLEE.State.ParcelState(50)
    st.status[1] = 0x01; st.nest[1] = Int32(1)
    st.indx[1] = 3.0f0; st.indy[1] = 3.0f0
    st.indl[1] = Float32(indl0)
    return st
end

# strat test level: 26 km → indl = 53 (z = 26,000 m); σ = p/ps ≈ 0.024 → w = 1
const IND26 = 53.0

println("="^80)
println("[1] isentropic pinning erases a fake kinematic displacement (w = 1)")
println("="^80)
st = mkstate(IND26); hyb = mkhyb(diabatic=false)
hybrid_vertical!(hyb, st, fctx, 900.0, [1.0])           # init θ_p at 26 km
θp0 = hyb.theta_p[1]
st.indl[1] += 2.0f0                                     # fake kinematic +1 km jump
hybrid_vertical!(hyb, st, fctx, 900.0, [1.0])
check("parcel pulled back to its isentrope (Δindl < 0.01 lvl)",
      abs(st.indl[1] - IND26) < 0.01,
      info=@sprintf("indl %.4f (target %.1f)", st.indl[1], IND26))
check("θ_p unchanged by the excursion", isapprox(hyb.theta_p[1], θp0; rtol=1e-5),
      info=@sprintf("θ_p %.2f K", hyb.theta_p[1]))
# repeat with a DOWNWARD excursion
st.indl[1] -= 3.0f0
hybrid_vertical!(hyb, st, fctx, 900.0, [1.0])
check("downward excursion erased too", abs(st.indl[1] - IND26) < 0.01)

println(); println("="^80)
println("[2] diabatic drift: descent = θ̇·dt through the ladder, θ_p prognostic")
println("="^80)
st2 = mkstate(IND26); hyb2 = mkhyb(diabatic=true)
hybrid_vertical!(hyb2, st2, fctx, 900.0, [1.0])         # init
θp_init = hyb2.theta_p[1]
nsteps = 96                                              # one day at 15 min
for _ in 1:nsteps
    hybrid_vertical!(hyb2, st2, fctx, 900.0, [1.0])
end
# expected: θ̇ = (θ/T)·Q with θ≈θp, T=200 K → Δθ_day = θp/200·(−1 K)
θp_exp = θp_init
for _ in 1:nsteps
    global θp_exp += θp_exp / 200.0 * Float32(Q_KDAY / 86400.0) * 900.0
end
check("θ_p follows the exact θ̇ integral over one day",
      isapprox(hyb2.theta_p[1], θp_exp; rtol=2e-4),
      info=@sprintf("θ_p %.3f vs %.3f K", hyb2.theta_p[1], θp_exp))
check("parcel indl sits on the θ_p ladder position",
      isapprox(θat(st2.indl[1]), hyb2.theta_p[1]; rtol=2e-4))
# dθ/dz here = 35·(500/1000)/500 m ... compute directly from the ladder:
dθdz = (theta[3,3,54,1] - theta[3,3,53,1]) / 500.0
dz_day = (hyb2.theta_p[1] - θp_init) / dθdz
check("implied descent ≈ (θ/T)·Q/(dθ/dz) (km/day scale sane)",
      -200.0 < dz_day < -50.0,   # fixture: (θ/T)·Q/(dθ/dz) = 2.9·(−1)/23.7 ≈ −121 m/day
      info=@sprintf("%.1f m/day (θ̇ %.2f K/day, dθ/dz %.1f K/km)",
                    dz_day, (hyb2.theta_p[1]-θp_init), dθdz*1000))

println(); println("="^80)
println("[3] no drivers → no motion")
println("="^80)
st3 = mkstate(IND26); hyb3 = mkhyb(diabatic=false)
hybrid_vertical!(hyb3, st3, fctx, 900.0, [1.0])
i0 = st3.indl[1]; t0 = hyb3.theta_p[1]
for _ in 1:96; hybrid_vertical!(hyb3, st3, fctx, 900.0, [1.0]); end
check("indl static over one day", abs(st3.indl[1] - i0) < 1e-3,
      info=@sprintf("drift %.2e lvl", st3.indl[1] - i0))
check("θ_p static", isapprox(hyb3.theta_p[1], t0; rtol=1e-6))

println(); println("="^80)
println("[4] tropospheric parcel untouched (σ > σ_r)")
println("="^80)
st4 = mkstate(5.0)                                       # 2 km: σ ≈ 0.75
hyb4 = mkhyb(diabatic=true)
hybrid_vertical!(hyb4, st4, fctx, 900.0, [1.0])
check("indl untouched in the kinematic region", st4.indl[1] == 5.0f0)
check("θ_p tracks the grid there",
      isapprox(hyb4.theta_p[1], θat(5.0); rtol=1e-5))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
