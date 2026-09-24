# test_available_cache.jl — the pre-generated AVAILABLE index in the met directory
# ([meteorology].available = "auto", 2026-09-24): accepted when it matches the directory,
# rejected when a file is missing/empty/added, namelist line rewritten on copy.
#   julia --project=easlee test/test_available_cache.jl
using Dates, Printf
include(joinpath(@__DIR__, "..", "easlee", "src", "io", "met_scan.jl"))
npass = 0; nfail = 0
function check(label, ok, detail="")
    global npass, nfail; ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end
d = mktempdir()
mk(rel) = (p = joinpath(d, rel); mkpath(dirname(p)); write(p, "x"); p)
files = ["20260603/00/gfs.t00z.pgrb2.0p25.f003", "20260603/00/gfs.t00z.pgrb2.0p25.f006",
         "20260603/12/gfs.t12z.pgrb2.0p25.f003"]
foreach(mk, files); write(joinpath(d, "README.txt"), "not met")
idx = joinpath(d, AVAILABLE_INDEX_NAME)
open(idx, "w") do io
    println(io, "/some/other/tree/listvariables/GFS.vars"); println(io, "YYYYMMDD HHMISS"); println(io, "________ ______")
    println(io, "20260603 030000      20260603/00/gfs.t00z.pgrb2.0p25.f003  ON DISC")
    println(io, "20260603 060000      20260603/00/gfs.t00z.pgrb2.0p25.f006  ON DISC")
    println(io, "20260603 150000      20260603/12/gfs.t12z.pgrb2.0p25.f003  ON DISC")
end
check("list_main_met_files skips README and lists 3 files", list_main_met_files(d) == sort(files))
e = read_available_entries(idx)
check("read_available_entries: 3 entries, first date", length(e) == 3 && e[1][1] == DateTime(2026, 6, 3, 3))
ok, why = cached_available_ok(idx, d); check("index without sidecar, nothing newer: accepted", ok, why)
sleep(1.1); mk("20260604/00/gfs.t00z.pgrb2.0p25.f003")
ok, why = cached_available_ok(idx, d); check("no sidecar: a file newer than the index → rescan", !ok, why)
rm(joinpath(d, "20260604/00/gfs.t00z.pgrb2.0p25.f003"))
write_available_sidecar!(idx, d)
check("sidecar written with 3 files + sizes", isfile(sidecar_path(idx)) && length(readlines(sidecar_path(idx))) == 3)
ok, why = cached_available_ok(idx, d); check("with sidecar: matching directory accepted", ok, why)
mk("20260604/00/gfs.t00z.pgrb2.0p25.f003")
ok, why = cached_available_ok(idx, d); check("with sidecar: an added file is detected", !ok, why)
rm(joinpath(d, "20260604/00/gfs.t00z.pgrb2.0p25.f003"))
write(joinpath(d, files[1]), "xy")
ok, why = cached_available_ok(idx, d); check("with sidecar: a changed size is detected", !ok, why)
write(joinpath(d, files[1]), "x")
write(joinpath(d, files[2]), "")
ok, why = cached_available_ok(idx, d); check("an empty listed file is detected", !ok, why)
write(joinpath(d, files[2]), "x"); rm(joinpath(d, files[3]))
ok, why = cached_available_ok(idx, d); check("a missing listed file is detected", !ok, why)
mk(files[3])
# companions: only MAIN files count (ERA5-like tokens)
foreach(mk, ["IFS_EA_MET_20220118_600_1.grb", "IFS_EA_CLD_20220118_600_1.grb", "IFS_PT_MET_20220118_600_1.grb"])
main = list_main_met_files(d; include_only=[r"_MET_"], exclude=[r"IFS_PT_MET"])
check("companion tokens: _MET_ kept, _CLD_ and PT excluded", "IFS_EA_MET_20220118_600_1.grb" in main && !any(occursin("CLD", f) || occursin("PT_MET", f) for f in main))
dst = joinpath(d, "AVAILABLE_1"); copy_available_with_vars(idx, dst, "/my/tree/listvariables/GFS.vars")
check("copy rewrites line 1 with the run's namelist", readlines(dst)[1] == "/my/tree/listvariables/GFS.vars" && length(readlines(dst)) == 6)
println("\n$(npass) passed, $(nfail) failed"); nfail == 0 || exit(1)
