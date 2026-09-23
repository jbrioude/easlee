#!/usr/bin/env julia
# test_met_gap_guard.jl — the met-coverage guard ([meteorology].max_gap):
#   _check_met_window refuses to run when the AVAILABLE timeline has a coverage
#   gap (missing files) larger than max_gap WITHIN the run window — a big gap
#   means a MISSING met file that would be silently interpolated/held across.
#   [1] continuous hourly met over the window → OK,
#   [2] a >6 h internal gap straddling the window → error,
#   [3] a >6 h trailing overrun (stop past last met) → error; a ≤6 h one → OK,
#   [4] max_gap = 0 disables the guard,
#   [5] start before the first met time → error (the original check),
#   [6] config: [meteorology].max_gap parse + round-trip.
#
#   julia --project=easlee test/test_met_gap_guard.jl
using Printf, Dates

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
    @printf("  %-60s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# write a synthetic AVAILABLE file (2 header lines + "YYYYMMDD HHMMSS  file  ON DISC")
function mkavail(times::Vector{DateTime})
    p = joinpath(mktempdir(), "AVAILABLE")
    open(p, "w") do io
        println(io, "ERA5_CAPTEX"); println(io, "YYYYMMDD HHMISS"); println(io, "___ ___ ___")
        for t in times
            println(io, Dates.format(t, "yyyymmdd"), " ", Dates.format(t, "HHMMSS"),
                    "   f.grb  ON DISC")
        end
    end
    return p
end
errs(f) = try f(); false catch; true end

d0 = DateTime(1983, 10, 14, 7)
hourly(a, b) = [a + Hour(k) for k in 0:Int(Dates.value(b - a) ÷ 3_600_000)]

println("="^78); println("[1] continuous hourly met over the window → OK"); println("="^78)
av = mkavail(hourly(d0, DateTime(1983, 10, 17, 19)))
check("no gap, window inside archive → no error",
      !errs(() -> _check_met_window(av, DateTime(1983,10,14,16), DateTime(1983,10,17,16), "ERA5")))

println(); println("="^78); println("[2] a >6 h internal gap straddling the window → error"); println("="^78)
# hourly to Oct-16 07:00, then jump to Oct-17 07:00 (24 h gap: the missing Oct-17 day)
gapped = vcat(hourly(d0, DateTime(1983,10,17,7)),  # NOTE this is continuous…
              )
# build a real internal hole: drop everything on Oct-16 after 07:00 up to Oct-17 07:00
holed = filter(t -> !(t > DateTime(1983,10,16,7) && t < DateTime(1983,10,17,7)),
               hourly(d0, DateTime(1983,10,18,7)))
avh = mkavail(holed)
check("24 h internal hole inside window → error",
      errs(() -> _check_met_window(avh, DateTime(1983,10,14,16), DateTime(1983,10,17,16), "ERA5")))
check("same hole but window ends BEFORE it → no error",
      !errs(() -> _check_met_window(avh, DateTime(1983,10,14,16), DateTime(1983,10,16,6), "ERA5")))

println(); println("="^78); println("[3] trailing overrun"); println("="^78)
avt = mkavail(hourly(d0, DateTime(1983,10,17,7)))   # last met 07:00
check("stop 9 h past last met (>6 h) → error",
      errs(() -> _check_met_window(avt, DateTime(1983,10,14,16), DateTime(1983,10,17,16), "ERA5")))
check("stop 3 h past last met (≤6 h) → OK (held frame)",
      !errs(() -> _check_met_window(avt, DateTime(1983,10,14,16), DateTime(1983,10,17,10), "ERA5")))

println(); println("="^78); println("[4] max_gap = 0 disables the guard"); println("="^78)
check("disabled: 9 h overrun → no error",
      !errs(() -> _check_met_window(avt, DateTime(1983,10,14,16), DateTime(1983,10,17,16), "ERA5"; max_gap_s=0.0)))
check("disabled: 24 h internal hole → no error",
      !errs(() -> _check_met_window(avh, DateTime(1983,10,14,16), DateTime(1983,10,17,16), "ERA5"; max_gap_s=0.0)))

println(); println("="^78); println("[5] start before first met → error (original check)"); println("="^78)
check("start before the archive → error",
      errs(() -> _check_met_window(av, DateTime(1983,10,14,0), DateTime(1983,10,17,16), "ERA5")))

println(); println("="^78); println("[6] config: [meteorology].max_gap parse + round-trip"); println("="^78)
base = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"1983-10-14T16:00:00","stop"=>"1983-10-17T16:00:00"),
  "species"=>[Dict{String,Any}("id"=>"X","kind"=>"gas")],
  "release"=>[Dict{String,Any}("species"=>"X","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["1983-10-14T16:00:00","1983-10-14T17:00:00"],
      "particles"=>10,"mass"=>"1kg")])
d = deepcopy(base); d["meteorology"]=[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E","max_gap"=>"3h")]
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("max_gap = 3h parses to 10800 s", sp.meteorology[1].max_met_gap_s == 10800.0)
d2 = deepcopy(base); d2["meteorology"]=[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")]
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(d2)
check("default max_gap = 6 h when unset", sp2.meteorology[1].max_met_gap_s == 6*3600.0)
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp3 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
check("round-trip preserves max_gap = 3h", sp3.meteorology[1].max_met_gap_s == 10800.0)
check("default is NOT serialized (kept implicit)",
      !haskey(EasleeConfig.ConfigTOML.runspec_to_dict(sp2)["meteorology"][1], "max_gap"))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
