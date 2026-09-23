#!/usr/bin/env julia
# test_umbrella.jl — umbrella-intrusion spreading (VOLCANIC_PLUME.md):
# the gravity-current outflow R(t) = (3λN·Q_V·t²/2π)^(1/3) applied to detrained
# volcanic parcels during transport (the archived met has no eruption outflow).
# Checks:
#   [1] ZERO-TUNING HT closure: with Q_V from the SAME 1-D column validated on
#       the heights (MER·amp/ρ_NBH) and N from the sounding, R(1 h) falls in the
#       observed 150–350 km range; R ∝ t^(2/3); front speed tens of m/s,
#   [2] kinematics: parcels pushed RADIALLY outward, interior ∝ r/R, beyond the
#       front at front speed; mass/vertical untouched; other boxes untouched,
#   [3] dateline wrap on a global grid (shortest-way-around offsets),
#   [4] handoff: no push before t0 or beyond tmax,
#   [5] [release.plume_rise] umbrella_spread TOML round-trip.
#
#   julia --project=easlee test/test_umbrella.jl
using Printf, Random, Statistics, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: UmbrellaIntrusion, umbrella_radius, umbrella_edge_speed,
                         apply_umbrella_spread!, build_plume_profiles,
                         volcanic_column_provenance

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── [1] zero-tuning Hunga Tonga closure ───────────────────────────────────────
println("="^78)
println("[1] R(t) law with the 1-D column's own volume flux (no tuning)")
println("="^78)
z = collect(0.0:200.0:60_000.0); T = similar(z)
for (i, zz) in enumerate(z)
    T[i] = zz <= 16_000 ? 300.0 - 6.5e-3*zz : zz <= 20_000 ? 196.0 : 196.0 + 2.2e-3*(zz-20_000)
end
p = similar(z); p[1] = 101_325.0
for i in 2:length(z); p[i] = p[i-1]*exp(-9.81*(z[i]-z[i-1])/(287.05*0.5*(T[i-1]+T[i]))); end
qv = fill(2e-6, length(z))
MER = 5.0e9; T0 = 700.0; NW = 0.45
prof = build_plume_profiles(MER, T0, NW, z, T, p, qv)
_, _, _, amp, _, _, _, _ = volcanic_column_provenance(MER, T0, NW, z, T, p, qv)
k = clamp(searchsortedlast(z, prof.z_nbh), 2, length(z) - 1)
ρnbh = p[k] / (287.05 * T[k])
θk(kk) = T[kk] * (100000.0 / p[kk])^0.2857
N = sqrt(max(9.81 / θk(k) * (θk(k+1) - θk(k-1)) / (z[k+1] - z[k-1]), 1e-6))
Qv = MER * amp / ρnbh
umb = UmbrellaIntrusion(qv_flux=Qv, nbv=N, t0=DateTime(2022,1,15,4,0,0),
                        vlon=-175.38, vlat=-20.54)
R1 = umbrella_radius(umb, 3600.0); R4 = umbrella_radius(umb, 4*3600.0)
V1 = umbrella_edge_speed(umb, 3600.0)
@printf("  amp %.2f  ρ_NBH %.4f kg/m³  Q_V %.2e m³/s  N %.4f 1/s\n", amp, ρnbh, Qv, N)
@printf("  R(1h) = %.0f km   R(4h) = %.0f km   Ṙ(1h) = %.1f m/s\n", R1/1e3, R4/1e3, V1)
check("R(1h) in the observed 150–350 km range (obs ≈ 200–250)",
      150e3 <= R1 <= 350e3; info=@sprintf("(%.0f km)", R1/1e3))
check("R ∝ t^(2/3)", isapprox(R4 / R1, 4^(2/3); rtol=1e-6))
check("front speed tens of m/s (obs ≈ 35–50)", 25.0 <= V1 <= 65.0;
      info=@sprintf("(%.1f m/s)", V1))
check("R(0) = 0 and Ṙ(0) = 0",
      umbrella_radius(umb, 0.0) == 0.0 && umbrella_edge_speed(umb, 0.0) == 0.0)

# ── [2] kinematics on a synthetic global 1° grid ──────────────────────────────
println()
println("="^78)
println("[2] parcel kinematics (radial, ∝ r/R inside, front speed outside)")
println("="^78)
g = EASLEE.Types.Rectilinear(Int32(360), Int32(181), Int32(10), 1.0f0, 1.0f0,
                -180.0f0, -90.0f0, 1.0f0/111320.0f0, 1.0f0/111320.0f0,
                1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false,
                1.0f0, 1.0f0, true)
# vent at (0°N, 0°E) → indy=91, indx=181; use age 1h (R1e, Rd1e for THIS umb2)
umb2 = UmbrellaIntrusion(qv_flux=Qv, nbv=N, t0=DateTime(2022,1,15,4,0,0),
                        vlon=0.0, vlat=0.0)
R1e = umbrella_radius(umb2, 3600.0); Rd1e = umbrella_edge_speed(umb2, 3600.0)
st = EASLEE.State.ParcelState(5)
st.status .= 0x01; st.relbox .= Int32(7); st.indl .= 3.0f0; st.mass .= 2.0f0
deg_half = (0.5 * R1e) / 111_320.0        # r = R/2 east
deg_out  = (1.5 * R1e) / 111_320.0        # r = 1.5R north
st.indx[1] = 181.0f0 + Float32(deg_half); st.indy[1] = 91.0f0    # inside, east
st.indx[2] = 181.0f0; st.indy[2] = 91.0f0 + Float32(deg_out)     # outside, north
st.indx[3] = 181.0f0; st.indy[3] = 91.0f0                        # at the axis
st.indx[4] = 181.0f0 - Float32(deg_half); st.indy[4] = 91.0f0    # inside, west
st.indx[5] = 200.0f0; st.indy[5] = 120.0f0; st.relbox[5] = Int32(9)  # other box
x0 = copy(st.indx); y0 = copy(st.indy)
dt = 60.0
n = apply_umbrella_spread!(st, g, true, 360, umb2, Set([7]),
                           DateTime(2022,1,15,5,0,0), dt)
check("pushed exactly the in-box off-axis parcels", n == 3; info="($n)")
dx1 = (st.indx[1] - x0[1]) * 111_320.0
# source-fed interior (2026-07-07): u(R/2) = R²/(2·(R/2)·t) = R/t = 1.5·Ṙ —
# the flow that carries vent material to the head (the old u ∝ r profile
# conserved each parcel's front-fraction, freezing the footprint at ~½R)
check("interior parcel (r=R/2) rides the source-fed flow at 1.5·Ṙ·dt",
      isapprox(dx1, 1.5 * Rd1e * dt; rtol=0.05);
      info=@sprintf("(%.0f m vs %.0f m)", dx1, 1.5 * Rd1e * dt))
dy2 = (st.indy[2] - y0[2]) * 111_320.0
check("front-runner (r=1.5R) moves at front speed",
      isapprox(dy2, Rd1e * dt; rtol=0.05);
      info=@sprintf("(%.0f m vs %.0f m)", dy2, Rd1e * dt))
check("westward parcel pushed WEST (radial, not advective)",
      st.indx[4] < x0[4])
check("purely radial: cross-components unchanged",
      st.indy[1] == y0[1] && st.indx[2] == x0[2])
check("axis parcel, other-box parcel untouched; mass/vertical untouched",
      st.indx[3] == x0[3] && st.indx[5] == x0[5] && st.indy[5] == y0[5] &&
      all(st.indl .== 3.0f0) && all(st.mass .== 2.0f0))
# dome saturation: inside r < rsrc the source-fed flow plateaus at u(rsrc)
stD = EASLEE.State.ParcelState(1)
stD.status .= 0x01; stD.relbox .= Int32(7); stD.indl .= 3.0f0; stD.mass .= 1.0f0
deg_dome = (0.5 * umb2.rsrc) / 111_320.0          # r = rsrc/2, inside the dome
stD.indx[1] = 181.0f0 + Float32(deg_dome); stD.indy[1] = 91.0f0
xD0 = stD.indx[1]
apply_umbrella_spread!(stD, g, true, 360, umb2, Set([7]),
                       DateTime(2022,1,15,5,0,0), 1.0)   # 1 s: displacement ≪ r
dxD = (stD.indx[1] - xD0) * 111_320.0
uD  = R1e^2 / (2.0 * 3600.0 * umb2.rsrc)
check("dome parcel (r=rsrc/2) rides the SATURATED outflow u(rsrc)",
      isapprox(dxD, uD * 1.0; rtol=0.05);
      info=@sprintf("(%.0f m/s vs %.0f m/s)", dxD, uD))

# wind shielding: uniform ambient Ux = U; in-cloud parcel with u_r > U gets the
# full resolved-advection cancellation (net eastward motion = (u_r − U)·dt for a
# parcel due-east of the vent); an out-of-cloud parcel (r > 1.2R) is untouched
umbW = UmbrellaIntrusion(qv_flux=Qv, nbv=N, t0=DateTime(2022,1,15,4,0,0),
                         vlon=0.0, vlat=0.0, wind_shield=true)
RW = umbrella_radius(umbW, 3600.0); RdW = umbrella_edge_speed(umbW, 3600.0)
met = (Ux = fill(Float32(10.0), 360, 180, 10, 2),
       Vy = fill(Float32(0.0),  360, 180, 10, 2))
stw = EASLEE.State.ParcelState(2)
stw.status .= 0x01; stw.relbox .= Int32(7); stw.indl .= 3.0f0; stw.mass .= 1.0f0
stw.indx[1] = 181.0f0 + Float32((0.5 * RW) / 111_320.0); stw.indy[1] = 91.0f0
stw.indx[2] = 181.0f0 + Float32((2.0 * RW) / 111_320.0); stw.indy[2] = 91.0f0
xw0 = copy(stw.indx)
apply_umbrella_spread!(stw, g, true, 360, umbW, Set([7]),
                       DateTime(2022,1,15,5,0,0), 60.0; met=met, indt=1.0)
dx_in  = (stw.indx[1] - xw0[1]) * 111_320.0
ur_in  = max(RW^2 / (2*3600.0*max(0.5*RW, umbW.rsrc)), 0.5*RdW)
check("wind shield: in-cloud parcel net = (u_r − U)·dt (full shield, u_r ≫ U)",
      isapprox(dx_in, (ur_in - 10.0) * 60.0; rtol=0.05);
      info=@sprintf("(%.0f m vs %.0f m)", dx_in, (ur_in - 10.0)*60.0))
dx_out = (stw.indx[2] - xw0[2]) * 111_320.0
check("wind shield: beyond 1.2R the resolved wind is NOT cancelled",
      isapprox(dx_out, RdW * 60.0; rtol=0.05);
      info=@sprintf("(%.0f m vs push-only %.0f m)", dx_out, RdW*60.0))

# ── [3] dateline wrap ─────────────────────────────────────────────────────────
println()
println("="^78)
println("[3] dateline wrap (global grid, vent at 179.5°E, parcel at −179.5°E)")
println("="^78)
umb3 = UmbrellaIntrusion(qv_flux=Qv, nbv=N, t0=DateTime(2022,1,15,4,0,0),
                         vlon=179.5, vlat=0.0)
st3 = EASLEE.State.ParcelState(1)
st3.status .= 0x01; st3.relbox .= Int32(1); st3.indl .= 3.0f0
st3.indy[1] = 91.0f0
st3.indx[1] = (( -179.5) - (-180.0)) / 1.0 + 1.0    # = 1.5 (just east of the seam)
x30 = st3.indx[1]
apply_umbrella_spread!(st3, g, true, 360, umb3, Set([1]),
                       DateTime(2022,1,15,5,0,0), 60.0)
moved_m = (st3.indx[1] - x30) * 111_320.0
check("parcel across the seam is pushed EAST by ~front speed (not −359°)",
      0.0 < moved_m < 2.0 * Rd1e * 60.0;
      info=@sprintf("(%.0f m)", moved_m))

# ── [4] handoff / gating ──────────────────────────────────────────────────────
println()
println("="^78)
println("[4] time gating")
println("="^78)
st4 = EASLEE.State.ParcelState(1); st4.status .= 0x01; st4.relbox .= Int32(1)
st4.indx[1] = 190.0f0; st4.indy[1] = 91.0f0
check("no push before t0",
      apply_umbrella_spread!(st4, g, true, 360, umb2, Set([1]),
                             DateTime(2022,1,15,3,0,0), 60.0) == 0)
check("no push beyond tmax (default 12 h)",
      apply_umbrella_spread!(st4, g, true, 360, umb2, Set([1]),
                             DateTime(2022,1,16,4,0,1), 60.0) == 0)

# ── [5] config round-trip ─────────────────────────────────────────────────────
println()
println("="^78)
println("[5] [release.plume_rise] umbrella_spread round-trip")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "volcanic_plume.toml"))
spec.release[1].plume_umbrella = true
spec.release[1].plume_spread_lambda = 0.25
EasleeConfig.validate_runspec(spec)
d = runspec_to_dict(spec)
pr = d["release"][1]["plume_rise"]
check("writes umbrella_spread + spread_lambda",
      pr["umbrella_spread"] === true && pr["spread_lambda"] == 0.25)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse-back preserves both",
      spec2.release[1].plume_umbrella && spec2.release[1].plume_spread_lambda == 0.25)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
