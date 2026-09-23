#!/usr/bin/env julia
# test_output_grouping.jl — [output].group: what the accumulator's 4th
# dimension is (initoutgrid column grouping + accumulation mapping):
#   [1] "species": one column per species; outclass maps boxes to species
#       columns; grid memory independent of box count (the 3400-box/70 GB
#       OOM class of trap #5),
#   [2] "release": identity mapping — bit-identical legacy per-box columns,
#   [3] defaults: forward → species, backward → release; IOUTPUTFOREACHRELEASE=0
#       → single collapsed column,
#   [4] _conc_nest! deposits through outclass (two boxes, same species, one
#       column) and relmass-per-column backward normalisation,
#   [5] config: [output].group parse / validation / round-trip.
#
#   julia --project=easlee test/test_output_grouping.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: _conc_nest!

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic runconfig/input: 5 boxes carrying 2 species (1,1,2,1,2) ────────
function mkinputs(; foreach="1", outgroup=nothing, direction=1)
    runconfig = Dict{String,Any}(
        "OUTGRID" => Dict("NUMXGRID"=>["10"], "NUMYGRID"=>["8"],
                          "OUTHEIGHTS"=>["0.0, 100.0, 200.0"],
                          "DXOUT"=>["1.0"], "DYOUT"=>["1.0"],
                          "OUTLON0"=>["0.0"], "OUTLAT0"=>["0.0"]),
        "nbreleasebox" => 5,
        "direction"    => direction,
        "nest"         => 1,
    )
    cmd = Dict{String,Any}("IOUTPUTFOREACHRELEASE"=>[foreach])
    outgroup !== nothing && (cmd["OUTPUTGROUP"] = [outgroup])
    input = Dict{String,Any}(
        "COMMAND" => cmd,
        "RELEASES_CTRL" => Dict{String,Any}("SPECNUM_REL" =>
            ["1", "1", "2", "1", "2"]),
    )
    metconfig = [Dict{String,Any}("projection"=>"latlon",
                                  "longitudes"=>[0.0], "latitudes"=>[0.0])]
    return runconfig, input, metconfig
end

println("="^78)
println("[1] group = species: one column per species")
println("="^78)
rc1, in1, mc1 = mkinputs(outgroup="species")
oc1 = initoutgrid(rc1, in1, mc1)
check("Nc = 2 (species), not 5 (boxes)", rc1["Nc"] == 2)
check("concarray 4th dim = 2", size(oc1["concarray"], 4) == 2)
check("outclass maps (1,1,2,1,2) → columns (1,1,2,1,2)",
      oc1["outclass"] == Int32[1, 1, 2, 1, 2])
check("classspecies = [1, 2]", oc1["classspecies"] == Int32[1, 2])
check("outgroup recorded", oc1["outgroup"] == "species")

println()
println("="^78)
println("[2] group = release: bit-identical legacy per-box columns")
println("="^78)
rc2, in2, mc2 = mkinputs(outgroup="release")
oc2 = initoutgrid(rc2, in2, mc2)
check("Nc = 5 (boxes)", rc2["Nc"] == 5)
check("outclass = identity", oc2["outclass"] == Int32[1, 2, 3, 4, 5])
check("classspecies = per-box species", oc2["classspecies"] == Int32[1, 1, 2, 1, 2])
rc2b, in2b, mc2b = mkinputs(outgroup="box")        # alias
oc2b = initoutgrid(rc2b, in2b, mc2b)
check("\"box\" alias = release", rc2b["Nc"] == 5 && oc2b["outgroup"] == "release")

println()
println("="^78)
println("[3] defaults: forward → species, backward → release; collapse")
println("="^78)
rc3, in3, mc3 = mkinputs()                          # no OUTPUTGROUP, forward
oc3 = initoutgrid(rc3, in3, mc3)
check("forward default = species", rc3["Nc"] == 2 && oc3["outgroup"] == "species")
rc4, in4, mc4 = mkinputs(direction=-1)              # backward
oc4 = initoutgrid(rc4, in4, mc4)
check("backward default = release", rc4["Nc"] == 5 && oc4["outgroup"] == "release")
rc5, in5, mc5 = mkinputs(foreach="0")               # legacy collapse
oc5 = initoutgrid(rc5, in5, mc5)
check("IOUTPUTFOREACHRELEASE=0 → single column, species 0 (mixed)",
      rc5["Nc"] == 1 && oc5["classspecies"] == Int32[0])

println()
println("="^78)
println("[4] deposit through outclass + per-column backward norm")
println("="^78)
# 10 levels of synthetic met; two parcels from DIFFERENT boxes (1 and 4) of the
# SAME species must land in the SAME column under species grouping.
const NZ = 10
dims = (12, 10, NZ, 2)
zlev = zeros(Float32, dims); rho = fill(1.0f0, dims)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * 30.0); end
GEO = EASLEE.Types.Rectilinear(Int32(12), Int32(10), Int32(NZ), 1.0f0, 1.0f0,
        0.0f0, -50.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC  = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
t0  = DateTime(2025, 1, 1)
function mkst()
    st = EASLEE.State.ParcelState(10)
    for (i, rb) in ((1, 1), (2, 4))                 # boxes 1 and 4 → species 1
        st.status[i] = 0x01; st.nest[i] = Int32(1); st.relbox[i] = Int32(rb)
        st.indx[i] = 5.0f0; st.indy[i] = 5.0f0; st.indl[i] = 3.0f0
        st.mass[i] = 2.0f0
        st.Dtime[i] = t0
    end
    return st
end
st  = mkst()
zi  = Float32[0.0, 100.0, 200.0]
g5  = zeros(Float32, 12, 10, 3, 2, 1); m5 = zeros(Float64, size(g5))
ocl = Int32[1, 1, 2, 1, 2]
d = _conc_nest!(g5, m5, st, t0, 1.0f0, 12, 10, 2, 1, zi, 1,
                zlev, rho, LC, GEO, false, 60.0, Float32[1.0f0, 1.0f0], :none;
                outclass = ocl)
check("both parcels deposit (4 kg total)", isapprox(sum(m5), 4.0; rtol=1e-12))
check("… into the SAME species column", isapprox(sum(m5[:, :, :, 1, 1]), 4.0) &&
      sum(m5[:, :, :, 2, 1]) == 0.0)
# empty outclass = identity (legacy): boxes 1 and 4 in their own columns
g6 = zeros(Float32, 12, 10, 3, 5, 1); m6 = zeros(Float64, size(g6))
st2 = mkst()
_conc_nest!(g6, m6, st2, t0, 1.0f0, 12, 10, 5, 1, zi, 1,
            zlev, rho, LC, GEO, false, 60.0, ones(Float32, 5), :none)
check("empty outclass = identity (legacy per-box)",
      isapprox(sum(m6[:, :, :, 1, 1]), 2.0) && isapprox(sum(m6[:, :, :, 4, 1]), 2.0))
# backward: w = dt/M_col with the PER-COLUMN released mass
g7 = zeros(Float32, 12, 10, 3, 2, 1); m7 = zeros(Float64, size(g7))
st3 = mkst()
relm = Float32[8.0, 4.0]                            # column 1 released 8 kg
d7 = _conc_nest!(g7, m7, st3, t0, 1.0f0, 12, 10, 2, 1, zi, 1,
                 zlev, rho, LC, GEO, true, 60.0, relm, :none; outclass = ocl)
check("backward norm uses per-column relmass (2 × 2kg × 60s/8kg = 30)",
      isapprox(sum(m7[:, :, :, 1, 1]), 2 * 2.0 * 60.0 / 8.0; rtol=1e-6))

println()
println("="^78)
println("[5] config: [output].group parse / validation / round-trip")
println("="^78)
d5 = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-16T00:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "output"=>Dict{String,Any}("group"=>"release"),
  "species"=>[Dict{String,Any}("id"=>"X","kind"=>"gas")],
  "release"=>[Dict{String,Any}("species"=>"X","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>10,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d5)
check("parse group=release", sp.output.group == "release")
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
check("round-trip", sp2.output.group == "release")
d5["output"]["group"] = "banana"
bad = try; EasleeConfig.ConfigTOML.dict_to_runspec(d5); false; catch; true; end
check("invalid group value rejected", bad)
d5["output"]["group"] = ""
spd = EasleeConfig.ConfigTOML.dict_to_runspec(d5)
check("empty group = auto (writer omits it)", spd.output.group == "" &&
      !haskey(get(EasleeConfig.ConfigTOML.runspec_to_dict(spd), "output", Dict()), "group"))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
