#!/usr/bin/env julia
# test_polar_stereo.jl — polar-stereographic transport caps
# ([physics].polar_stereo_lat, FLEXPART switchnorth/switchsouth analogue).
#
#   1. polar_stereo_update unit tests on a global ERA5-like 0.5° grid:
#      identity, mid-latitude equivalence with the lat/lon metric, north/south
#      pole crossings and near-pole east motion checked against exact
#      great-circle geometry (haversine distance + spherical triangle).
#   2. advance! integration on a coarse global grid: uniform northward wind
#      moves a cap parcel the right ground distance; a parcel crosses the pole
#      and comes out at lon+180 with the right latitude; south cap symmetric;
#      the lat/lon-path pole-overshoot mirror (switch row close to the pole).
#
# Usage:  julia --project=easlee scripts/test_polar_stereo.jl

using Dates, Printf

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Types: Rectilinear, SigmaPressureLevel, MetFields3D, MetFields2D,
                     NoTurb, NestRunConfig
const OP = EASLEE.Operators

npass = 0; nfail = 0
function check(name, ok)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-64s %s\n", name, ok ? "PASS" : "FAIL")
    return ok
end

const R_E = 6.371e6
havdist(lat1, lon1, lat2, lon2) = 2R_E * asind(min(1.0, sqrt(
    sind((lat2 - lat1) / 2)^2 +
    cosd(lat1) * cosd(lat2) * sind((lon2 - lon1) / 2)^2))) * π / 180
dlon360(a, b) = abs(mod(a - b + 180, 360) - 180)

mkrect(Ni, Nj, dl) = Rectilinear(
    Int32(Ni), Int32(Nj), Int32(4), Float32(dl), Float32(dl),
    -180.0f0, -90.0f0,
    Float32(180 / (dl * R_E * π)), Float32(180 / (dl * R_E * π)),
    1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, true,
    Float32(dl), Float32(dl), true)

ij(g, lat, lon) = OP.latlontoindyindx(g, Float32(lat), Float32(lon))

println("── 1. polar_stereo_update unit tests (global 0.5° grid) ──")
g = mkrect(720, 361, 0.5)

# identity: zero displacement
let (iy0, ix0) = ij(g, 82.3, 47.2)
    iy, ix = OP.polar_stereo_update(g, 1.0f0, iy0, ix0, 0.0f0, 0.0f0)
    check("zero displacement is the identity",
          abs(iy - iy0) < 1e-3 && abs(ix - ix0) < 1e-3)
end

# mid-latitude: matches the lat/lon metric update for a small step
let (iy0, ix0) = ij(g, 76.0, 10.0)
    dxs, dys = 1000.0f0, 500.0f0
    iy, ix   = OP.polar_stereo_update(g, 1.0f0, iy0, ix0, dxs, dys)
    ix_ll    = ix0 + OP.dxtodindx(g, dxs, ix0, iy0)
    iy_ll    = iy0 + OP.dytodindy(g, dys, ix0, iy0)
    check("76°N 1-km step ≈ lat/lon metric update",
          abs(iy - iy_ll) < 2e-4 && abs(ix - ix_ll) < 2e-4)
end

# north pole crossing: 50 km due north from 89.8°N/0°E → other side, lon 180
let (iy0, ix0) = ij(g, 89.8, 0.0)
    iy, ix = OP.polar_stereo_update(g, 1.0f0, iy0, ix0, 0.0f0, 50.0f3)
    lat, lon = OP.indyindxtolatlon(g, iy, ix)
    lat_exp  = 90 - (50e3 / (R_E * π / 180) - 0.2)     # exact meridian arc
    check(@sprintf("N pole crossing: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-3)
    check(@sprintf("N pole crossing: lon %.2f ≈ ±180", lon),
          dlon360(lon, 180.0) < 0.05)
    check("N pole crossing: ground distance preserved",
          abs(havdist(89.8, 0.0, lat, lon) - 50e3) / 50e3 < 5e-3)
end

# near-pole east motion: 5 km due east at 89.9°N — exact spherical triangle
let (iy0, ix0) = ij(g, 89.9, 0.0)
    iy, ix = OP.polar_stereo_update(g, 1.0f0, iy0, ix0, 5.0f3, 0.0f0)
    lat, lon = OP.indyindxtolatlon(g, iy, ix)
    b = 5e3 / R_E                       # chord arc [rad]
    a = 0.1 * π / 180                   # start colatitude [rad]
    colat_exp = acosd(cos(a) * cos(b))
    lon_exp   = atand(tan(b), sin(a))
    check(@sprintf("89.9°N east step: lat %.4f ≈ %.4f", lat, 90 - colat_exp),
          abs(lat - (90 - colat_exp)) < 5e-3)
    check(@sprintf("89.9°N east step: lon %.2f ≈ %.2f", lon, lon_exp),
          dlon360(lon, lon_exp) < 0.2)
end

# south pole crossing: 50 km due south from −89.8°/45°E → lon −135
let (iy0, ix0) = ij(g, -89.8, 45.0)
    iy, ix = OP.polar_stereo_update(g, -1.0f0, iy0, ix0, 0.0f0, -50.0f3)
    lat, lon = OP.indyindxtolatlon(g, iy, ix)
    lat_exp  = -(90 - (50e3 / (R_E * π / 180) - 0.2))
    check(@sprintf("S pole crossing: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-3)
    check(@sprintf("S pole crossing: lon %.2f ≈ −135", lon),
          dlon360(lon, -135.0) < 0.05)
    check("S pole crossing: ground distance preserved",
          abs(havdist(-89.8, 45.0, lat, lon) - 50e3) / 50e3 < 5e-3)
end

println("── 2. advance! integration (global 5° grid, uniform wind) ──")

# One-parcel advance! on a fresh coarse global grid.  Returns final (lat, lon,
# Dtime).  switch_lat sets both cap rows; NoTurb, no mesoscale, ctl = 2.
function run_one(lat0, lon0; u = 0.0f0, v = 0.0f0, ts_s = 3600, switch_lat = 75.0)
    Ni, Nj, Nk, dl = 72, 37, 4, 5.0
    gg  = mkrect(Ni, Nj, dl)
    A() = zeros(Float32, Ni + 1, Nj, Nk, 2)
    zl  = A(); for k in 1:Nk; zl[:, :, k, :] .= 1000.0f0 * (k - 1); end
    m3  = MetFields3D(fill(u, Ni + 1, Nj, Nk, 2), fill(v, Ni + 1, Nj, Nk, 2),
                      A(), A(), A(), A(), zl, A(), A(), A(), A(), A())
    m2  = MetFields2D((zeros(Float32, Ni + 1, Nj, 2) for _ in 1:4)...)
    lev = SigmaPressureLevel(Int32(Nk), zeros(Float32, Nk), zeros(Float32, Nk))
    rc  = NestRunConfig((1, 2), NoTurb(), 2.0f0, Int8(1), (true, true, true),
                        false, 0.16f0, 3600.0f0,
                        0.0f0,                            # kh_ft
                        false, false, NaN32, NaN32,       # v1.1: horiz_tke, vmix_native/rw/ck
                        true, true,
                        Float32(( switch_lat + 90) / dl + 1),
                        Float32((-switch_lat + 90) / dl + 1),
                        false)                            # mo_surface
    iy0, ix0 = ij(gg, lat0, lon0)
    traj = (indx = Float32[ix0], indy = Float32[iy0], indl = Float32[2.0],
            dt = Float32[10.0], uturb = Float32[0], vturb = Float32[0],
            wturb = Float32[0], mesou = Float32[0], mesov = Float32[0],
            mesow = Float32[0], miss2 = Float32[0], miss3 = Float32[0],
            miss4 = Float32[0], Dtime = [DateTime(2025, 1, 1)])
    OP.advance!(traj, Millisecond(ts_s * 1000), rc, m3, m2, gg, lev, 1.0f0, [1])
    lat, lon = OP.indyindxtolatlon(gg, traj.indy[1], traj.indx[1])
    return lat, lon, traj.Dtime[1]
end

# cap parcel, uniform 20 m/s northward for 1 h → 72 km up the meridian
let (lat, lon, dtm) = run_one(85.0, 10.0; v = 20.0f0)
    lat_exp = 85 + 72e3 / (R_E * π / 180)
    check(@sprintf("cap advance: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-3)
    check("cap advance: lon unchanged", dlon360(lon, 10.0) < 0.05)
    check("cap advance: parcel alive, clock advanced",
          dtm == DateTime(2025, 1, 1, 1))
end

# pole crossing through advance!: 80 m/s north from 87.5°N crosses the pole
let (lat, lon, dtm) = run_one(87.5, 0.0; v = 80.0f0)
    lat_exp = 90 - (288e3 / (R_E * π / 180) - 2.5)
    check(@sprintf("advance! N crossing: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-3)
    check(@sprintf("advance! N crossing: lon %.2f ≈ ±180", lon),
          dlon360(lon, 180.0) < 0.05)
    check("advance! N crossing: parcel alive", dtm == DateTime(2025, 1, 1, 1))
end

# south cap symmetric crossing
let (lat, lon, dtm) = run_one(-87.5, 45.0; v = -80.0f0)
    lat_exp = -(90 - (288e3 / (R_E * π / 180) - 2.5))
    check(@sprintf("advance! S crossing: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-3)
    check(@sprintf("advance! S crossing: lon %.2f ≈ −135", lon),
          dlon360(lon, -135.0) < 0.05)
    check("advance! S crossing: parcel alive", dtm == DateTime(2025, 1, 1, 1))
end

# lat/lon-path pole overshoot (switch row very close to the pole): the parcel
# stays on the metric path, overshoots indy > Nj, and must be MIRRORED across
# the pole (FLEXPART HSO/AL fix), not killed.
let (lat, lon, dtm) = run_one(88.0, 0.0; v = 80.0f0, switch_lat = 89.0)
    lat_exp = 180 - (88 + 288e3 / (R_E * π / 180))    # mirror of the overshoot
    check(@sprintf("pole mirror: lat %.4f ≈ %.4f", lat, lat_exp),
          abs(lat - lat_exp) < 5e-2)
    check(@sprintf("pole mirror: lon %.2f ≈ ±180", lon),
          dlon360(lon, 180.0) < 0.05)
    check("pole mirror: parcel alive (not killed at indy ≥ Nj)",
          dtm == DateTime(2025, 1, 1, 1))
end

@printf("\n%d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
