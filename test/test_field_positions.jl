#!/usr/bin/env julia
# test_field_positions.jl — where each mapped 3-D NetCDF field sits (FieldPos,
# _netcdf_field_positions in initialize.jl): WRF `stagger` / *_stag dimensions,
# MesoNH `grid` codes (face and mass dimensions of equal length), the fields'
# accepted positions, and the per-model restrictions.  Synthetic files only.
using Printf, Logging, NetCDF

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-64s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end
errs(f) = try; f(); false; catch; true; end

dir = mktempdir()
Ni, Nj, Km = 4, 3, 5; Nk = Km + 1                      # WRF: Nk = the w-levels

# WRF-like file holding both forms of the winds
wrf = joinpath(dir, "wrfout_test.nc")
for (n, dims, st) in (("U",  ("west_east_stag", Ni + 1, "south_north", Nj, "bottom_top", Km), "X"),
                      ("UC", ("west_east", Ni, "south_north", Nj, "bottom_top", Km), ""),
                      ("V",  ("west_east", Ni, "south_north_stag", Nj + 1, "bottom_top", Km), "Y"),
                      ("VC", ("west_east", Ni, "south_north", Nj, "bottom_top", Km), ""),
                      ("W",  ("west_east", Ni, "south_north", Nj, "bottom_top_stag", Nk), "Z"),
                      ("WC", ("west_east", Ni, "south_north", Nj, "bottom_top", Km), ""),
                      ("T",  ("west_east", Ni, "south_north", Nj, "bottom_top", Km), ""),
                      ("PH", ("west_east", Ni, "south_north", Nj, "bottom_top_stag", Nk), "Z"))
    NetCDF.nccreate(wrf, n, dims..., "Time", 1; atts = Dict("stagger" => st), t = Float32)
end
NetCDF.ncclose(wrf)
pos(f, model, map, ni, nj, nk) = Logging.with_logger(Logging.NullLogger()) do
    _netcdf_field_positions(f, model, collect(keys(map)), map, ni, nj, nk)
end

p = pos(wrf, "WRF", Dict("U" => "Ux", "V" => "Vy", "W" => "W", "T" => "T", "PH" => "ph"), Ni, Nj, Nk)
check("WRF U (stagger X): x-faces, full levels", p["Ux"].x == :face && p["Ux"].y == :mass && p["Ux"].z == :full)
check("WRF V (stagger Y): y-faces", p["Vy"].y == :face && p["Vy"].x == :mass)
check("WRF W (stagger Z): half levels (w-levels)", p["W"].z == :half)
check("WRF T: mass points, full levels", p["T"].x == :mass && p["T"].y == :mass && p["T"].z == :full)
check("WRF PH: half levels", p["ph"].z == :half)
p = pos(wrf, "WRF", Dict("UC" => "Ux", "VC" => "Vy", "WC" => "W"), Ni, Nj, Nk)
check("centred UC, VC: mass points", p["Ux"].x == :mass && p["Vy"].y == :mass)
check("WC on the mass levels: full levels", p["W"].z == :full)
p = pos(wrf, "WRF", Dict("U" => "Ux", "VC" => "Vy", "WC" => "W"), Ni, Nj, Nk)
check("mixed: staggered U with centred V and W", p["Ux"].x == :face && p["Vy"].y == :mass && p["W"].z == :full)
check("T mapped to a staggered U → error", errs(() -> pos(wrf, "WRF", Dict("U" => "T"), Ni, Nj, Nk)))
check("Vy mapped to an x-staggered U → error", errs(() -> pos(wrf, "WRF", Dict("U" => "Vy"), Ni, Nj, Nk)))
check("PH mapped as ph from a full-level field → error", errs(() -> pos(wrf, "WRF", Dict("T" => "ph"), Ni, Nj, Nk)))
check("AROME reader: face winds refused", errs(() -> pos(wrf, "AROME", Dict("U" => "Ux"), Ni, Nj, Nk)))
check("MERRA-2 reader: OMEGA on half levels refused", errs(() -> pos(wrf, "MERRA2", Dict("W" => "levdot"), Ni, Nj, Nk)))

# MesoNH-like file: face and mass dimensions of equal length, `grid` decides
mnh = joinpath(dir, "mesonh_test.nc"); L = 6
for (n, dims, g) in (("UT",  ("ni_u", Ni, "nj_u", Nj, "level", L), 2),
                     ("UM",  ("ni", Ni, "nj", Nj, "level", L), 1),
                     ("VT",  ("ni_v", Ni, "nj_v", Nj, "level", L), 3),
                     ("WT",  ("ni", Ni, "nj", Nj, "level_w", L), 4),
                     ("WM",  ("ni", Ni, "nj", Nj, "level", L), 1),
                     ("THT", ("ni", Ni, "nj", Nj, "level", L), 1),
                     ("TKW", ("ni", Ni, "nj", Nj, "level_w", L), 4))
    NetCDF.nccreate(mnh, n, dims..., "time", 1; atts = Dict("grid" => Int32(g)), t = Float32)
end
NetCDF.ncclose(mnh)
p = pos(mnh, "MESONH", Dict("UT" => "Ux", "VT" => "Vy", "WT" => "levdot", "THT" => "theta", "TKW" => "tke"), Ni, Nj, L)
check("MesoNH UT (grid 2, same length as mass): x-faces", p["Ux"].x == :face)
check("MesoNH VT (grid 3): y-faces", p["Vy"].y == :face)
check("MesoNH WT (grid 4, level_w): half levels", p["levdot"].z == :half)
check("MesoNH THT (grid 1): mass points, full levels", p["theta"].x == :mass && p["theta"].z == :full)
check("MesoNH TKE on level_w: half levels", p["tke"].z == :half)
p = pos(mnh, "MESONH", Dict("UM" => "Ux", "WM" => "levdot"), Ni, Nj, L)
check("MesoNH winds at the mass points (grid 1): no destaggering", p["Ux"].x == :mass)
check("MesoNH w on the mass levels: full levels", p["levdot"].z == :full)

println("\n$npass passed, $nfail failed")
nfail == 0 || exit(1)
