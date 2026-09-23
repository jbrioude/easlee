#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# test_met_crop.jl — domain crop of the met grid (_crop_metconfig! + gather).
#
# Synthetic only (no met archive): builds fake metconfig dicts and full-grid
# fields, crops them, and checks the window, the seam-crossing gather, the
# monotone shifted-origin axis, and the parcel→index mapping continuity across
# the dateline — the properties transport relies on.
# ─────────────────────────────────────────────────────────────────────────────
using Printf

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# a fake global ERA5-like frame: −180 … 179.5 at 0.5°, −90 … 90
mkmc(; Ni=720, Nj=361, lon1=-180.0, dlon=0.5, lat1=-90.0, dlat=0.5, globalx=true) = Dict{Any,Any}(
    "metfileformat" => "GRIB",
    "Ni" => Int32(Ni), "Nj" => Int32(Nj),
    "longitudes" => [lon1 + (i-1)*dlon for i in 1:Ni],
    "latitudes"  => [lat1 + (j-1)*dlat for j in 1:Nj],
    "globalx"    => globalx)

println("="^78)
println("[1] Pacific box 150–230°E (crosses the dateline) on a global −180..180 grid")
println("="^78)
mc = mkmc()
_crop_metconfig!(mc, 150.0, 0.0, 230.0, 50.0; halo=2)
iidx, jr = mc["crop"]
lons = mc["longitudes"]; lats = mc["latitudes"]
check("lon axis monotone ascending across the seam",
      all(diff(lons) .≈ 0.5); info=@sprintf("%.2f…%.2f", lons[1], lons[end]))
check("lon window covers box + halo", lons[1] <= 150.0 - 0.99 && lons[end] >= 230.0 + 0.99;
      info="[$(lons[1]), $(lons[end])]")
check("lat window covers box + halo", lats[1] <= 0.0 - 0.99 && lats[end] >= 50.0 + 0.99;
      info="[$(lats[1]), $(lats[end])]")
check("cropped grid is no longer zonally periodic", mc["globalx"] == false)
check("Ni/Nj match the gather window",
      Int(mc["Ni"]) == length(iidx) && Int(mc["Nj"]) == length(jr))
frac = Int(mc["Ni"]) * Int(mc["Nj"]) / (720 * 361)
check("crop is a small fraction of the full grid", frac < 0.15;
      info=@sprintf("%.1f%%", 100frac))

# the gather must pull the RIGHT columns: fill a full-grid field with
# f(i,j) = lon(i) (mod 360) and check the cropped columns equal the new axis
full = [mod(-180.0 + (i-1)*0.5, 360.0) for i in 1:720, j in 1:361]
sub  = full[iidx, jr]
check("seam-crossing gather pulls the correct source columns",
      all(abs.(sub .- [mod(l, 360.0) for l in lons, _ in 1:length(jr)]) .< 1e-9))

# parcel → fractional index continuity across the dateline, via the SAME
# normalisation transport uses (latlontoindyindx: mod(lon − lonmin, 360))
tolindx(lon) = mod(lon - lons[1], 360.0) / 0.5 + 1.0
i_w = tolindx(179.75)      # just west of the seam
i_e = tolindx(-179.75)     # just east (= 180.25°E)
check("index mapping continuous across the seam", abs(i_e - i_w - 1.0) < 1e-9;
      info=@sprintf("indx(179.75)=%.3f indx(180.25)=%.3f", i_w, i_e))
check("box corners map inside the window",
      1.0 <= tolindx(150.0) <= length(iidx) && 1.0 <= tolindx(230.0) <= length(iidx))
check("a lon far outside the box maps outside the window",
      tolindx(0.0) > length(iidx))

println()
println("="^78)
println("[2] contiguous box (CAPTEX-like, 260–300°E on 0–360 grid) — no seam")
println("="^78)
mc2 = mkmc(lon1=0.0)
_crop_metconfig!(mc2, 260.0, 30.0, 300.0, 50.0; halo=2)
iidx2, jr2 = mc2["crop"]
check("contiguous gather is a plain ascending range",
      all(diff(Int.(iidx2)) .== 1); info="i=$(iidx2[1])…$(iidx2[end])")
check("axis matches the gathered columns",
      mc2["longitudes"][1] ≈ (iidx2[1]-1)*0.5 && mc2["longitudes"][end] ≈ (iidx2[end]-1)*0.5)

println()
println("="^78)
println("[3] lat-only fallback: lon window ≥ full circle keeps the global axis")
println("="^78)
mc3 = mkmc()
_crop_metconfig!(mc3, 0.0, 10.0, 359.9, 60.0; halo=2)   # lon covers the circle
iidx3, jr3 = mc3["crop"]
check("full-circle lon request → lon uncropped, still periodic",
      Int(mc3["Ni"]) == 720 && mc3["globalx"] == true && length(iidx3) == 720)
check("latitude still cropped", Int(mc3["Nj"]) == length(jr3) && Int(mc3["Nj"]) < 361)

println()
println("="^78)
println("[4] non-global regional grid: window clamped to the grid")
println("="^78)
mc4 = mkmc(Ni=161, lon1=260.0, dlon=0.25, Nj=81, lat1=30.0, dlat=0.25, globalx=false)
_crop_metconfig!(mc4, 259.0, 35.0, 275.0, 45.0; halo=2)   # west edge outside grid
iidx4, jr4 = mc4["crop"]
check("window clamped at the grid edge", iidx4[1] == 1 && jr4[1] >= 1)
check("clamped axis starts at the grid origin", mc4["longitudes"][1] ≈ 260.0)

println()
println("="^78)
println("[5] guards")
println("="^78)
mc5 = mkmc(); mc5["metfileformat"] = "NETCDF"
_crop_metconfig!(mc5, 150.0, 0.0, 230.0, 50.0)
check("non-GRIB met: crop skipped with warning (grid untouched)",
      !haskey(mc5, "crop") && Int(mc5["Ni"]) == 720)
err = try _crop_metconfig!(mkmc(), 150.0, 95.0, 230.0, 99.0); false catch; true end
check("empty latitude window errors", err)

println()
println("="^78)
println("[7] top-level interpolation: no out-of-bounds stencil at k = Nk")
println("="^78)
# regression for the build_domainfill segfault: an exactly-integral top-level
# query must not build a k+1 stencil point (zero-weight OOB touch under
# @inbounds).  Convention-agnostic value checks: whatever time slot the
# engine's t-weighting picks, the k = Nk query must equal the k = Nk-1
# integral query + exactly one level step (levels differ by 100 here), for
# BOTH t arguments (the t = 2 case is the one whose OOB touch fell past the
# array end and segfaulted).
zf = zeros(Float32, 4, 3, 5, 2)
for t in 1:2, k in 1:5; zf[:, :, k, t] .= Float32(100k + t); end
f4l = EASLEE.Operators.interpolation4Dlevsingle
for tq in (1.0f0, 2.0f0)
    vtop = f4l(zf, 2.0f0, 2.0f0, 5.0f0, tq)
    vk4  = f4l(zf, 2.0f0, 2.0f0, 4.0f0, tq)
    vmid = f4l(zf, 2.0f0, 2.0f0, 3.5f0, tq)
    check("k = Nk equals k = Nk-1 + one level step (t = $tq)", vtop - vk4 == 100.0f0;
          info="top=$vtop k4=$vk4")
    check("interior half-level query is the level midpoint (t = $tq)",
          vmid - (vk4 - 100.0f0) == 50.0f0; info="mid=$vmid")
end

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || (println("FAIL"); exit(1))
