#!/usr/bin/env julia
# test_massweighted_mixing.jl — reduced-mass DIRECTED exchange
# (MixingConfig exchange = :massweighted, ADAPTIVE_KERNEL §7):
#   dq_i = Δt·λ_ij·[2m_j/(m_i+m_j)]·(q_j−q_i)  (mirrored for j)
# Checks on the adjacency + gather kernel (the single changed point):
#   [1] species MASS Σ m·q conserved exactly for wild mass ratios (50–80×),
#   [2] uniform q is an exact fixed point for any mass distribution,
#   [3] equal-mass limit identical to :equal (regression guarantee),
#   [4] direction: small parcels move fast, large parcels barely (entrainment),
#   [5] directed factor bounded in (0,2) — no stiffness from mass disparity.
#
#   julia --project=easlee -t 2 test/test_massweighted_mixing.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: _pairs_to_adjacency, _gather_exchange_dc!, MixingConfig

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

Random.seed!(3)
N = 400
pairs_i = Int[]; pairs_j = Int[]; lambdas = Float64[]
for _ in 1:1500                                    # random pair graph
    i = rand(1:N); j = rand(1:N)
    i == j && continue
    push!(pairs_i, min(i, j)); push!(pairs_j, max(i, j)); push!(lambdas, 0.02 * rand())
end
m = 10.0 .^ (2.0 .* rand(N))                       # masses spanning ~100×
q0 = rand(N, 3); q0[:, 3] .= 7.5                   # species 3 uniform

function stepq(q, adj...; dt=1.0)
    dc = Matrix{Float64}(undef, size(q)...)
    _gather_exchange_dc!(dc, q, adj..., dt)
    return q .+ dc
end

println("="^78)
println("[1–2] conservation + uniform fixed point (mass-weighted, 100× masses)")
println("="^78)
adjw = _pairs_to_adjacency(N, pairs_i, pairs_j, lambdas; masses=m)
q = copy(q0)
tot0 = [sum(m .* q[:, s]) for s in 1:3]
for _ in 1:200; global q = stepq(q, adjw...); end
tot1 = [sum(m .* q[:, s]) for s in 1:3]
check("Σ m·q conserved to round-off (both species)",
      all(abs.(tot1 .- tot0) .< 1e-9 .* abs.(tot0));
      info=@sprintf("(rel err %.1e)", maximum(abs.(tot1 .- tot0) ./ abs.(tot0))))
check("uniform q is an exact fixed point", all(q[:, 3] .== 7.5))
check("mixing actually converges toward the mass-weighted mean",
      std(q[:, 1]) < 0.25 * std(q0[:, 1]);
      info=@sprintf("(std %.3f → %.3f)", std(q0[:, 1]), std(q[:, 1])))

println()
println("="^78)
println("[3] equal-mass limit ≡ current behaviour")
println("="^78)
adj_eq  = _pairs_to_adjacency(N, pairs_i, pairs_j, lambdas)
adj_eqm = _pairs_to_adjacency(N, pairs_i, pairs_j, lambdas; masses=fill(3.7, N))
check("adjacency weights identical for equal masses",
      adj_eq[3] == adj_eqm[3])

println()
println("="^78)
println("[4–5] directionality + boundedness")
println("="^78)
# one heavy background parcel (m=100) paired with one light plume parcel (m=1)
adj2 = _pairs_to_adjacency(2, [1], [2], [0.1]; masses=[100.0, 1.0])
q2 = [0.0 0.0; 100.0 0.0]                          # plume carries q=100
dc2 = Matrix{Float64}(undef, 2, 2)
_gather_exchange_dc!(dc2, q2, adj2..., 1.0)
check("light parcel moves ~2λ·Δq, heavy parcel ~λ·Δq/50",
      isapprox(dc2[2, 1], 0.1 * (2 * 100 / 101) * (0.0 - 100.0); rtol=1e-12) &&
      isapprox(dc2[1, 1], 0.1 * (2 * 1 / 101) * (100.0 - 0.0); rtol=1e-12);
      info=@sprintf("(heavy %+0.3f, light %+0.3f)", dc2[1, 1], dc2[2, 1]))
check("directed factors bounded in (0,2)",
      all(0.0 .< adj2[3] ./ 0.1 .< 2.0))
check("pair species mass conserved (m₁dq₁ = −m₂dq₂)",
      isapprox(100.0 * dc2[1, 1], -1.0 * dc2[2, 1]; rtol=1e-12))
check("MixingConfig accepts/validates the option",
      MixingConfig(exchange=:massweighted).exchange === :massweighted &&
      (try MixingConfig(exchange=:bogus); false catch; true end))

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
