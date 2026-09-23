#!/usr/bin/env julia
# test_transport_mixing.jl — [physics.mixing] inter-trajectory mixing config:
#   [1] defaults: TransportMixing() off; an off spec emits no mixing table,
#   [2] TOML parse: configs/volcanic_mixing.toml → pairwise, K_h,
#       sigma_v as a length string, background, every,
#   [3] round-trip: runspec_to_dict → dict_to_runspec preserves every field,
#   [4] sugar: mixing = "pairwise" (bare string) enables the scheme,
#   [5] validation: bad scheme / K_h ≤ 0 / negative widths / every < 1 rejected.
#
#   julia --project=easlee test/test_transport_mixing.jl
using Printf

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[1] defaults")
println("="^78)
mx = EasleeConfig.TransportMixing()
check("scheme defaults to off", mx.scheme == "off")
check("K_h / sigma0_h / sigma_v / nsig / max_neighbors / background / every",
      mx.K_h == 1000.0 && mx.sigma0_h == 0.0 && mx.sigma_v == 0.0 &&
      mx.nsig == 2.0 && mx.max_neighbors == 30 && mx.background && mx.every == 1)

println()
println("="^78)
println("[2] TOML parse (configs/volcanic_mixing.toml)")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "volcanic_mixing.toml"))
m = spec.physics.mixing
check("scheme = pairwise", m.scheme == "pairwise")
check("K_h = 1000 m²/s", m.K_h == 1000.0)
check("sigma_v = \"300m\" parses to 300 m", m.sigma_v == 300.0)
check("sigma0_h omitted → 0 (auto: met grid spacing)", m.sigma0_h == 0.0)
check("background = true, every = 1", m.background && m.every == 1)
check("nomix twin parses with mixing off",
      read_runspec(joinpath(REPO, "configs",
                            "volcanic_nomix.toml")).physics.mixing.scheme == "off")

println()
println("="^78)
println("[3] round-trip")
println("="^78)
spec.physics.mixing.sigma0_h = 25000.0
spec.physics.mixing.nsig = 2.5
spec.physics.mixing.max_neighbors = 20
spec.physics.mixing.every = 3
dct = runspec_to_dict(spec)
md = get(dct["physics"], "mixing", Dict())
check("pairwise spec emits the [physics.mixing] table",
      md isa AbstractDict && md["scheme"] == "pairwise" && md["K_h"] == 1000.0 &&
      md["sigma0_h"] == 25000.0 && md["sigma_v"] == 300.0 && md["nsig"] == 2.5 &&
      md["max_neighbors"] == 20 && md["background"] == true && md["every"] == 3)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(dct)
m2 = spec2.physics.mixing
check("parse-back preserves every field",
      m2.scheme == "pairwise" && m2.K_h == 1000.0 && m2.sigma0_h == 25000.0 &&
      m2.sigma_v == 300.0 && m2.nsig == 2.5 && m2.max_neighbors == 20 &&
      m2.background && m2.every == 3)
spec.physics.mixing = EasleeConfig.TransportMixing()      # off again
check("off spec emits NO mixing table",
      !haskey(runspec_to_dict(spec)["physics"], "mixing"))

println()
println("="^78)
println("[4] sugar: mixing = \"pairwise\"")
println("="^78)
dct2 = deepcopy(dct)
dct2["physics"]["mixing"] = "pairwise"
m3 = EasleeConfig.ConfigTOML.dict_to_runspec(dct2).physics.mixing
check("bare string enables pairwise with defaults",
      m3.scheme == "pairwise" && m3.K_h == 1000.0 && m3.sigma0_h == 0.0 &&
      m3.background && m3.every == 1)

println()
println("="^78)
println("[5] validation")
println("="^78)
function rejects(mutate!)
    s = EasleeConfig.ConfigTOML.dict_to_runspec(deepcopy(dct))
    mutate!(s.physics.mixing)
    try EasleeConfig.validate_runspec(s); false catch; true end
end
check("valid pairwise spec validates",
      (EasleeConfig.validate_runspec(
           EasleeConfig.ConfigTOML.dict_to_runspec(deepcopy(dct))); true))
check("scheme = \"bogus\" rejected",    rejects(m -> m.scheme = "bogus"))
check("K_h = 0 rejected",               rejects(m -> m.K_h = 0.0))
check("sigma0_h < 0 rejected",          rejects(m -> m.sigma0_h = -1.0))
check("sigma_v < 0 rejected",           rejects(m -> m.sigma_v = -1.0))
check("every = 0 rejected",             rejects(m -> m.every = 0))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
