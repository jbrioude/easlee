#!/usr/bin/env julia
# test_species_convert.jl — [species.convert] stochastic in-flight conversion:
#   [1] conversion statistics: fraction converted over one step ≈ 1−exp(−dt/τ)
#       (binomial tolerance); mass rescaled by the ratio; expectation of the
#       total mass conserved × ratio,
#   [2] repeated steps: exponential decay of the parent population; converted
#       parcels never convert again (shadow box has no target); airm reset,
#   [3] budget accumulator = parent mass converted,
#   [4] config: [species.convert] parse / auto mass ratio / validation.
#
#   julia --project=easlee test/test_species_convert.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: SpeciesConvert, apply_convert!, any_convert

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

N = 200_000
τ = 3600.0; dt = 600.0
pexp = -expm1(-dt / τ)
mkcv() = SpeciesConvert(Float32[1/τ, 0], Int32[2, 0], Float32[1.5, 1.0], zeros(2))
function mkstate()
    st = EASLEE.State.ParcelState(N)
    st.status .= 0x01; st.relbox .= Int32(1); st.mass .= 2.0f0
    st.airm .= 7.0f0
    return st
end

println("="^78)
println("[1] one-step statistics")
println("="^78)
st = mkstate(); cv = mkcv()
tot = apply_convert!(st, cv, dt)
nconv = count(==(Int32(2)), st.relbox)
check("converted fraction ≈ 1−exp(−dt/τ)",
      abs(nconv / N - pexp) < 4 * sqrt(pexp * (1 - pexp) / N),
      info=@sprintf("%.4f vs %.4f", nconv / N, pexp))
i2 = findfirst(==(Int32(2)), st.relbox)
i1 = findfirst(==(Int32(1)), st.relbox)
check("converted parcel: mass × ratio, airm reset",
      st.mass[i2] == 3.0f0 && st.airm[i2] == 0.0f0)
check("unconverted parcel untouched",
      st.mass[i1] == 2.0f0 && st.airm[i1] == 7.0f0)
check("budget = parent mass converted",
      isapprox(cv.converted[1], nconv * 2.0; rtol=1e-12) &&
      isapprox(tot, cv.converted[1]; rtol=1e-12))
check("any_convert true / all-zero targets false",
      any_convert(cv) && !any_convert(SpeciesConvert(Float32[0], Int32[0], Float32[1], zeros(1))))

println()
println("="^78)
println("[2] repeated steps: exponential parent decay, daughters inert")
println("="^78)
st2 = mkstate(); cv2 = mkcv()
nsteps = 12
for _ in 1:nsteps
    apply_convert!(st2, cv2, dt)
end
frac_parent = count(==(Int32(1)), st2.relbox) / N
fexp = exp(-nsteps * dt / τ)
check("parent survival ≈ exp(−nsteps·dt/τ)",
      abs(frac_parent - fexp) < 5 * sqrt(fexp * (1 - fexp) / N),
      info=@sprintf("%.4f vs %.4f", frac_parent, fexp))
check("daughters converted exactly once (mass = 2·1.5, never 2·1.5²)",
      all(st2.mass[i] == 3.0f0 for i in eachindex(st2.relbox) if st2.relbox[i] == Int32(2)))

println()
println("="^78)
println("[3] config: [species.convert]")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-16T00:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[
    Dict{String,Any}("id"=>"SO2","kind"=>"gas","molar_mass"=>64.06,
        "convert"=>Dict{String,Any}("to"=>"SO4a","tau"=>"25d")),
    Dict{String,Any}("id"=>"SO4a","kind"=>"aerosol","aerosol_type"=>"sulfate",
        "molar_mass"=>96.06,"density"=>1770.0)],
  "release"=>[Dict{String,Any}("species"=>"SO2","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>100,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
EasleeConfig.validate_runspec(sp)
s1 = sp.species[1]
check("parse + validate", s1.convert_to == "SO4a" &&
      isapprox(s1.convert_tau_s, 25 * 86400.0) && s1.convert_mass_ratio == 0.0)
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
s2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt).species[1]
check("round-trip", s2.convert_to == "SO4a" && isapprox(s2.convert_tau_s, 25 * 86400.0))
d2 = deepcopy(d); d2["species"][1]["convert"]["to"] = "NOPE"
bad = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d2)); false
      catch; true end
check("convert.to must name a defined species", bad)
d3 = deepcopy(d); d3["species"][1]["convert"]["tau"] = "0s"
bad3 = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d3)); false
       catch; true end
check("tau = 0 rejected", bad3)

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
