#!/usr/bin/env julia
# test_track_nudge.jl — observed-track ensemble nudging + wind-bias drift:
#   [1] track loading: centroids, dateline-straddling boxes, lon unwrap,
#   [2] nudge kinematics: exponential relaxation of the ensemble centroid,
#       RIGID shift (internal structure exactly preserved), height-band gating,
#       free-run outside the track span,
#   [3] dateline: ensemble at +179.8 nudged toward a track at −179.8 moves the
#       SHORT way; global-x wrap keeps indx in range,
#   [4] wind-bias drift: Δlon = u·dt/(R cosφ) at each parcel's own latitude.
#
#   julia --project=easlee test/test_track_nudge.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: load_track_nudge, apply_track_nudge!, WindBias,
                         apply_wind_bias!
const ll = EASLEE.Operators.indyindxtolatlon   # (io/tools.jl shadows the name in Main)

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# global 1° met grid (lonmin 0, latmin -90 → indx = lon+1 for lon ≥ 0), 10 levels
G = EASLEE.Types.Rectilinear(Int32(360), Int32(180), Int32(10), 1.0f0, 1.0f0,
        0.0f0, -90.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, true)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
zlev = Array{Float32,4}(undef, 360, 180, 10, 2)
for k in 1:10; zlev[:, :, k, :] .= Float32((k - 1) * 4000.0); end   # 0…36 km
t0 = DateTime(2022, 1, 15, 4, 0)

ll2i(lon, lat) = (Float32(mod(lon, 360.0) + 1.0), Float32(lat + 91.0))
function mkstate(lons, lats, lev)
    st = EASLEE.State.ParcelState(length(lons))
    st.status .= 0x01; st.nest .= Int32(1); st.Dtime .= t0
    for (i, (lo, la)) in enumerate(zip(lons, lats))
        st.indx[i], st.indy[i] = ll2i(lo, la)
        st.indl[i] = Float32(lev)
    end
    return st
end

println("="^78)
println("[1] track loading")
println("="^78)
trk = joinpath(mktempdir(), "track.csv")
write(trk, """
time,lon_west,lon_east,lat_south,lat_north
2022-01-15T04:00:00,182.0,186.0,-22.0,-18.0
2022-01-15T10:00:00,178.0,182.0,-24.0,-20.0
2022-01-15T16:00:00,174.0,178.0,-24.0,-20.0
""")
tn = load_track_nudge(trk, 3600.0, [18000.0, 40000.0])
check("centroids parsed ((182+186)/2 = 184, (-22-18)/2 = -20)",
      tn.clon[1] ≈ 184.0 && tn.clat[1] ≈ -20.0)
check("track lon stays continuous (no wrap jumps)",
      all(abs.(diff(tn.clon)) .< 90.0))
trk2 = joinpath(dirname(trk), "wrap.csv")
write(trk2, "time,lon_west,lon_east,lat_south,lat_north\n" *
            "2022-01-15T04:00:00,178.0,-178.0,-22.0,-18.0\n" *
            "2022-01-15T10:00:00,176.0,-180.0,-22.0,-18.0\n")
tn2 = load_track_nudge(trk2, 3600.0, [18000.0, 40000.0])
check("dateline-straddling box centre = 180 (short arc)",
      abs(mod(tn2.clon[1], 360.0) - 180.0) < 1e-9;
      info=@sprintf("(clon %.2f)", tn2.clon[1]))

println()
println("="^78)
println("[2] nudge kinematics")
println("="^78)
# 100 parcels in a 1°×1° blob at (184, -22.5): obs track at lat -20 → nudge N
lons = 183.5 .+ rand(100); lats = -23.0 .+ rand(100)
st = mkstate(lons, lats, 7.0)          # level 7 → z = 24 km (in band)
rel0 = [(st.indx[i] - st.indx[1], st.indy[i] - st.indy[1]) for i in 1:100]
dt = 900.0
la1, lo1 = ll(G, Float64(st.indy[1]), Float64(st.indx[1]))
nmv, dlo, dla = apply_track_nudge!(st, tn, G, LC, zlev, 1.0, t0, dt)
la2, lo2 = ll(G, Float64(st.indy[1]), Float64(st.indx[1]))
f = 1.0 - exp(-dt / 3600.0)
check("all in-band parcels moved; increment = (obs−ens)·(1−e^(−dt/τ))",
      nmv == 100 && isapprox(dla, (-20.0 - (-22.5)) * f; atol=0.02);
      info=@sprintf("(dlat %.3f vs %.3f)", dla, 2.5f))
check("shift applied to the parcels (lat moved north by dlat)",
      isapprox(la2 - la1, dla; atol=1e-4))
rel1 = [(st.indx[i] - st.indx[1], st.indy[i] - st.indy[1]) for i in 1:100]
check("RIGID shift: internal structure bit-preserved",
      all(rel0[i] == rel1[i] for i in 1:100))
# exponential relaxation over many steps
for _ in 1:40
    apply_track_nudge!(st, tn, G, LC, zlev, 1.0, t0, dt)
end
cla = sum(ll(G, Float64(st.indy[i]), Float64(st.indx[i]))[1] for i in 1:100) / 100
check("centroid relaxes onto the track (40 × 15 min ≫ τ)",
      isapprox(cla, -20.0; atol=0.05); info=@sprintf("(clat %.3f)", cla))
# band gating: parcels below the band don't move
stb = mkstate([184.0], [-22.5], 2.0)   # level 2 → z = 4 km (below band)
x0 = stb.indx[1]
nmvb, _, _ = apply_track_nudge!(stb, tn, G, LC, zlev, 1.0, t0, dt)
check("below-band parcel untouched", nmvb == 0 && stb.indx[1] == x0)
# free-run outside the track span
stc = mkstate([184.0], [-22.5], 7.0)
nmvc, _, _ = apply_track_nudge!(stc, tn, G, LC, zlev, 1.0, t0 + Hour(13), dt)
check("outside the track span → free run (no nudge)", nmvc == 0)

println()
println("="^78)
println("[3] dateline handling")
println("="^78)
trk3 = joinpath(dirname(trk), "dl.csv")
write(trk3, "time,lon_west,lon_east,lat_south,lat_north\n" *
            "2022-01-15T04:00:00,-179.0,-176.0,-22.0,-18.0\n" *
            "2022-01-15T16:00:00,-179.0,-176.0,-22.0,-18.0\n")
tn3 = load_track_nudge(trk3, 600.0, [18000.0, 40000.0])   # short τ
std = mkstate(fill(179.2, 50), -20.0 .+ 0.5 .* rand(50), 7.0)
for _ in 1:30
    apply_track_nudge!(std, tn3, G, LC, zlev, 1.0, t0 + Hour(1), 900.0)
end
clo = 0.0
for i in 1:50
    _, lo = ll(G, Float64(std.indy[i]), Float64(std.indx[i]))
    global clo += 182.5 + (mod(lo - 182.5 + 180.0, 360.0) - 180.0)
end
clo /= 50
check("ensemble crosses the dateline the SHORT way onto the track (−177.5 ≡ 182.5)",
      isapprox(clo, 182.5; atol=0.1); info=@sprintf("(clon %.2f)", clo))
check("global-x wrap keeps indx in [1, Ni+1)",
      all(1.0f0 .<= std.indx[1:50] .< 361.0f0))

println()
println("="^78)
println("[4] wind-bias drift")
println("="^78)
wb = WindBias(2.0, -1.0, 18000.0, 40000.0)    # 2 m/s E, 1 m/s S
ste = mkstate([184.0], [-20.0], 7.0)
la0, lo0 = ll(G, Float64(ste.indy[1]), Float64(ste.indx[1]))
n = apply_wind_bias!(ste, wb, G, LC, zlev, 1.0, 21600.0)   # 6 h in one step
la9, lo9 = ll(G, Float64(ste.indy[1]), Float64(ste.indx[1]))
dlon_exp = 2.0 * 21600.0 / (111320.0 * cosd(-20.0))
dlat_exp = -1.0 * 21600.0 / 111320.0
check("Δlon = u·dt/(R cosφ) at the parcel's own latitude",
      n == 1 && isapprox(lo9 - lo0, dlon_exp; rtol=1e-3);
      info=@sprintf("(%.4f° vs %.4f°)", lo9 - lo0, dlon_exp))
check("Δlat = v·dt/R", isapprox(la9 - la0, dlat_exp; rtol=1e-3))
stf = mkstate([184.0], [-20.0], 2.0)          # below band
check("below-band parcel gets no drift",
      apply_wind_bias!(stf, wb, G, LC, zlev, 1.0, 21600.0) == 0)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
