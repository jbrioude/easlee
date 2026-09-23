#!/usr/bin/env julia
# test_mesoscale.jl — mesoscale velocity fluctuations (FLEXPART usig/vsig/wsig):
# the sub-grid wind uncertainty term added to the advection.  Checks:
#   [1] _std16 — exact σ of the 2×2×2×2 grid-point neighbourhood (checkerboard
#       analytic, constant → 0, NaN → 0 guard),
#   [2] _mesoscale_update! OU statistics — stationary std = fturbmeso·σ and
#       lag-1 autocorrelation = r (the FLEXPART turbulence_mesoscale recursion),
#   [3] fresh-parcel full-variance draw (no spin-up bias),
#   [4] [physics.turbulence] mesoscale/fturbmeso TOML round-trip + validation.
#
#   julia --project=easlee scripts/test_mesoscale.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: _std16, _mesoscale_sigmas, _mesoscale_update!

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[1] _std16: neighbourhood σ (analytic)")
println("="^78)
# checkerboard ±A around 10: exactly 8 high / 8 low in any 2×2×2×2 block → σ = A
A = 3.0f0
cb = Array{Float32}(undef, 4, 4, 4, 2)
for l in 1:2, k in 1:4, j in 1:4, i in 1:4
    cb[i, j, k, l] = 10.0f0 + (isodd(i + j + k + l) ? A : -A)
end
check("checkerboard ±3 → σ = 3", isapprox(_std16(cb, 1, 1, 1), A; atol=1e-4))
check("constant field → σ = 0", _std16(fill(7.5f0, 4, 4, 4, 2), 2, 2, 2) == 0.0f0)
nanf = fill(1.0f0, 4, 4, 4, 2); nanf[2, 2, 2, 1] = NaN32
check("NaN in the neighbourhood → σ = 0 (masked-cell guard)",
      _std16(nanf, 1, 1, 1) == 0.0f0)
us, vs, ws = _mesoscale_sigmas(cb, fill(7.5f0, 4, 4, 4, 2), cb, 1.4f0, 1.6f0, 1.2f0, 4)
check("_mesoscale_sigmas routes the three fields", us == A && vs == 0.0f0 && ws == A)

println()
println("="^78)
println("[2] OU statistics: stationary std = fm·σ, lag-1 autocorrelation = r")
println("="^78)
Random.seed!(1234)
st = EASLEE.State.ParcelState(1)
usig, vsig, wsig = 5.0f0, 3.0f0, 0.02f0
fm = 0.16f0
r  = exp(-2.0f0 * 900.0f0 / 3600.0f0)        # Δt=15 min, T=1 h → r = e^-0.5
rs = sqrt(1.0f0 - r^2)
N = 200_000
u = Vector{Float32}(undef, N)
st.mesou[1] = 1.0f-6                          # arm past the fresh-draw branch
for n in 1:N
    _mesoscale_update!(st, 1, usig, vsig, wsig, r, rs, fm)
    u[n] = st.mesou[1]
end
burn = 100
σ_emp = std(u[burn:end])
ac1 = cor(u[burn:end-1], u[burn+1:end])
check("stationary std ≈ fturbmeso·usig", isapprox(σ_emp, fm * usig; rtol=0.02);
      info=@sprintf("(%.3f vs %.3f)", σ_emp, fm * usig))
check("lag-1 autocorrelation ≈ r = e^-0.5", isapprox(ac1, r; atol=0.01);
      info=@sprintf("(%.3f vs %.3f)", ac1, r))
check("mean ≈ 0", abs(mean(u[burn:end])) < 0.02 * fm * usig)

println()
println("="^78)
println("[3] fresh parcel: full-variance draw, w in level-index units")
println("="^78)
Random.seed!(7)
draws = Float32[]
for _ in 1:20_000
    s2 = EASLEE.State.ParcelState(1)
    _mesoscale_update!(s2, 1, usig, vsig, wsig, r, rs, fm)
    push!(draws, s2.mesou[1])
end
check("first-draw std ≈ fm·usig (no OU spin-up deficit)",
      isapprox(std(draws), fm * usig; rtol=0.03);
      info=@sprintf("(%.3f vs %.3f)", std(draws), fm * usig))

println()
println("="^78)
println("[4] TOML round-trip + validation")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "era5_pacific.toml"))
spec.physics.turbulence.mesoscale = true
spec.physics.turbulence.fturbmeso = 0.12
EasleeConfig.validate_runspec(spec)
d = runspec_to_dict(spec)
td = d["physics"]["turbulence"]
check("round-trip writes mesoscale + fturbmeso",
      td["mesoscale"] === true && td["fturbmeso"] == 0.12)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse-back preserves the fields",
      spec2.physics.turbulence.mesoscale && spec2.physics.turbulence.fturbmeso == 0.12)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
