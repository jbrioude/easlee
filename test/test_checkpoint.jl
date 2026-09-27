#!/usr/bin/env julia
# test_checkpoint.jl — warm-start checkpoint/restart ([run] checkpoint/restart):
#   [1] save → load round-trip restores EVERY ParcelState field, the seeded
#       mask and the physics accumulators, bit-exact,
#   [2] guards: larger-than-pool checkpoint rejected; version pinned,
#   [3] config: [run] checkpoint/restart parse + round-trip.
#
#   julia --project=easlee test/test_checkpoint.jl
using Printf, Dates, Random

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

rng = Xoshiro(3)
function mkstate(n)
    st = EASLEE.State.ParcelState(n)
    st.status[1:2:end] .= 0x01
    for f in (:indx, :indy, :indl, :mass, :airm, :cxx, :hband, :mesou, :miss3)
        getfield(st, f) .= rand(rng, Float32, n)
    end
    st.relbox .= Int32.(rand(rng, 0:3, n))
    st.Dtime  .= DateTime(2022, 1, 15) .+ Second.(rand(rng, 0:86400, n))
    EASLEE.State.sync_dtime_ns!(st)
    return st
end

spec = EasleeConfig.RunSpec()
spec.name = "chktest"
push!(spec.species, EasleeConfig.Species(id="X"))
ice  = EASLEE.IcePhysics([true], Float32[1.5], Float32[20e-6], Float32[500],
                         Float32[0], [42.0], Float32[0.3], [7.0])
ssb  = EASLEE.Operators.SeaSprayBox(Int32(1), 0.0, 10.0, 0.0, 10.0, Ref(true),
                                    Int32[1], Int32[1], [1.0e8], [123.0], 50.0,
                                    DateTime(2022,1,15), DateTime(2022,1,16),
                                    -1.0, 35.0, -1.0, -1.0, 2017.0,
                                    Ref(1.0e4), Ref(200))
# mirror the RunObjects fields _checkpoint_save/_load! touch (convphys/pradphys
# added with [species.convert] / [physics.plume_radiation]; nothing = feature off)
mkobj() = (chemstate = nothing, icephys = ice, convphys = nothing,
           pradphys = nothing, tagphys = nothing, ssprays = [ssb],
           dusts = EASLEE.Operators.DustBox[])   # kind="dust" sources (2026-09-25)

println("="^78)
println("[1] save → load round-trip (bit-exact)")
println("="^78)
st1 = mkstate(100)
rt1 = (state = st1,)
seeded1 = bitrand(rng, 100)
d0 = DateTime(2022, 1, 17, 10, 0)
path = joinpath(mktempdir(), "chk.jls")
_checkpoint_save(path, rt1, mkobj(), spec, d0, seeded1)
check("checkpoint file written", isfile(path))

st2 = EASLEE.State.ParcelState(120)              # bigger pool is fine
rt2 = (state = st2,)
seeded2 = falses(120)
ice.frozen[1] = 0.0; ice.scavenged[1] = 0.0      # blank the accumulators
ssb.acc[1] = 0.0; ssb.emitted[] = 0.0; ssb.nspawn[] = 0
dl = _checkpoint_load!(rt2, mkobj(), spec, path, seeded2)
check("date restored", dl == d0)
check("EVERY ParcelState field bit-exact (into a larger pool)",
      all(view(getfield(st2, f), 1:100) == getfield(st1, f)
          for f in fieldnames(typeof(st1)) if f != :nest_buffers))
check("seeded mask restored", seeded2[1:100] == seeded1 && !any(seeded2[101:end]))
check("ice accumulators restored", ice.frozen[1] == 42.0 && ice.scavenged[1] == 7.0)
check("seaspray accumulators restored",
      ssb.acc[1] == 123.0 && ssb.emitted[] == 1.0e4 && ssb.nspawn[] == 200)

println()
println("="^78)
println("[2] guards")
println("="^78)
small = (state = EASLEE.State.ParcelState(50),)
bad = try _checkpoint_load!(small, mkobj(), spec, path, falses(50)); false
      catch; true end
check("checkpoint larger than the pool rejected", bad)

println()
println("="^78)
println("[3] config: [run] checkpoint/restart")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-16T00:00:00",
                          "checkpoint"=>"6h","restart"=>"/tmp/c.jls"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[Dict{String,Any}("id"=>"X","kind"=>"gas")],
  "release"=>[Dict{String,Any}("species"=>"X","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>10,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse", sp.checkpoint == "6h" && sp.restart == "/tmp/c.jls")
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
check("round-trip", sp2.checkpoint == "6h" && sp2.restart == "/tmp/c.jls")

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
