# test_source_mass_quantum.jl — the emitted parcel mass equals the integrated flux, for the
# sea-spray (Ovadnevaite et al. 2014) and dust (Kok et al. 2014) sources: 12 steps on a synthetic
# grid with cell- and time-varying winds; the parcels' masses are summed from the state and
# compared with Σ_steps Σ_cells F·A·Δt evaluated independently from the pure flux functions.
#   julia --project=easlee test/test_source_mass_quantum.jl
using Printf, Dates, Random
global gpu = false; global chem = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "scripts", "_source_quantum_harness.jl"))
npass = 0; nfail = 0
function check(label, ok, detail="")
    global npass, nfail; ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end
r = source_quantum_harness()
for (name, x) in (("sea spray", r.spray), ("dust", r.dust))
    println("[$name] cells $(x.ncells), spawned $(x.n) parcels of $(x.q) kg")
    check("$name: parcels carry exactly n × quantum", x.spawned == x.n * x.q && x.nparcels == x.n, @sprintf("%.6e vs %.6e kg", x.spawned, x.n * x.q))
    check("$name: emitted counter = parcel mass", x.emitted == x.spawned)
    check("$name: spawned + unspawned = ∫F dA dt (independent, rtol 1e-9)", isapprox(x.spawned + x.unspawned, x.integral; rtol=1e-9), @sprintf("%.10e vs %.10e", x.spawned + x.unspawned, x.integral))
    check("$name: every cell residual < one quantum", x.maxacc < x.q, @sprintf("max %.3g of %.3g", x.maxacc, x.q))
    check("$name: spawned / integral within ncells·q of 1", abs(x.spawned / x.integral - 1) <= x.ncells * x.q / x.integral + 1e-12, @sprintf("%.5f", x.spawned / x.integral))
    check("$name: parcels at the surface layer (indl < 2)", x.indl < 2.0, @sprintf("max indl %.2f", x.indl))
end
println("\n$(npass) passed, $(nfail) failed"); nfail == 0 || exit(1)
