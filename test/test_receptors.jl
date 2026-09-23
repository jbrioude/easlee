#!/usr/bin/env julia
# test_receptors.jl — gridless receptor output (KDE point evaluation):
#   [1] C(x_s) matches the analytic Σ m·G₂(r;h)/Δz for known parcels,
#   [2] layer filter + 3h truncation + R(x_s) mass-weighted mean h,
#   [3] window averaging semantics (mean over nsamp samples),
#   [4] DATEM station file loader, [5] config round-trip + validation.
#
#   julia --project=easlee -t 2 test/test_receptors.jl
using Printf, Random, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: ReceptorSet, receptorset_from_file, receptor_sample!,
                         receptor_flush!, receptor_write

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# synthetic grid: 1° Rectilinear, lon0=0 lat0=0 at index 1 → lat = indy-1, lon = indx-1
G = EASLEE.Types.Rectilinear(Int32(100), Int32(100), Int32(10), 1.0f0, 1.0f0,
        0.0f0, 0.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
zlev = Array{Float32,4}(undef, 100, 100, 10, 2)
for k in 1:10; zlev[:, :, k, :] .= Float32((k - 1) * 100.0); end   # levels every 100 m

println("="^78)
println("[1–2] point evaluation vs analytic, layer filter, truncation, R")
println("="^78)
# station at lat=30, lon=30 (indx=31, indy=31); layer 0–100 m
rs = ReceptorSet(["S1"], [30.0], [30.0], 0.0, 100.0)
st = EASLEE.State.ParcelState(4)
st.status .= 0x01; st.nest .= Int32(1); st.mass .= 2.0f0
# parcel 1: 10 km east of the station, h=20 km, z=50 m (in layer)
doff = 10_000.0 / (111_320.0 * cosd(30.0))
st.indx[1] = Float32(31.0 + doff); st.indy[1] = 31.0f0; st.indl[1] = 1.5f0
st.hband[1] = 20_000.0f0
# parcel 2: at the station, h=5 km, but z=250 m (OUT of layer)
st.indx[2] = 31.0f0; st.indy[2] = 31.0f0; st.indl[2] = 3.5f0
st.hband[2] = 5_000.0f0
# parcel 3: 100 km away, h=20 km → beyond 3h, must not contribute
st.indx[3] = Float32(31.0 + 10 * doff); st.indy[3] = 31.0f0; st.indl[3] = 1.5f0
st.hband[3] = 20_000.0f0
# parcel 4: 5 km north, h=10 km, z=50 m (in layer)
st.indx[4] = 31.0f0; st.indy[4] = Float32(31.0 + 5_000.0 / 111_320.0); st.indl[4] = 1.5f0
st.hband[4] = 10_000.0f0
receptor_sample!(rs, st, G, LC, zlev, 1.0, 1)
G2(r, h) = exp(-r^2 / (2h^2)) / (2π * h^2)
c_ref = (2.0 * G2(10_000.0, 20_000.0) + 2.0 * G2(5_000.0, 10_000.0)) / 100.0
check("C(x_s) = analytic Σ m·G₂/Δz (in-layer, in-support parcels only)",
      isapprox(rs.csum[1], c_ref; rtol=2e-3);
      info=@sprintf("(%.3e vs %.3e)", rs.csum[1], c_ref))
w1 = 2.0 * G2(10_000.0, 20_000.0); w4 = 2.0 * G2(5_000.0, 10_000.0)
R_ref = (w1 * 20_000.0 + w4 * 10_000.0) / (w1 + w4)
check("R(x_s) = mass-weighted mean h of contributors",
      isapprox(rs.rnum[1] / rs.rden[1], R_ref; rtol=2e-3);
      info=@sprintf("(%.0f vs %.0f m)", rs.rnum[1] / rs.rden[1], R_ref))

println()
println("="^78)
println("[3] window averaging")
println("="^78)
csamp1 = rs.csum[1]
receptor_sample!(rs, st, G, LC, zlev, 1.0, 1)      # second identical sample
receptor_flush!(rs, 3600.0)
check("mean over 2 samples = single-sample value",
      isapprox(rs.rows[1][3], csamp1; rtol=1e-12) && rs.rows[1][5] == 2.0)
check("accumulators reset after flush", rs.csum[1] == 0.0 && rs.nsamp == 0)
out = tempname() * ".csv"
receptor_write(rs, out, DateTime(1994, 10, 23, 12))
lines = readlines(out)
check("CSV written (header + 1 row)",
      any(l -> startswith(l, "time_s"), lines) &&
      any(l -> startswith(l, "3600,S1,"), lines))

println()
println("="^78)
println("[4] DATEM station file loader")
println("="^78)
stf = tempname()
open(stf, "w") do io
    println(io, "# stations"); println(io, "101  48.05  -2.01"); println(io, "102  50.10   8.60")
end
rs2 = receptorset_from_file(stf, 0.0, 100.0)
check("id lat lon parsed", rs2.ids == ["101", "102"] &&
      rs2.lat == [48.05, 50.10] && rs2.lon == [-2.01, 8.60])

println()
println("="^78)
println("[5] config round-trip + validation")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "output_etex.toml"))
spec.output.kernel = "adaptive"
spec.output.receptors_file = stf
spec.output.receptors_layer = [0.0, 100.0]
EasleeConfig.validate_runspec(spec)
d = runspec_to_dict(spec)
check("writes [output.receptors] file + layer",
      d["output"]["receptors"]["file"] == stf &&
      d["output"]["receptors"]["layer_top"] == 100.0)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse-back preserves", spec2.output.receptors_file == stf)
# receptors are decoupled from the grid kernel (2026-07-06): with a uniform
# grid the engine refreshes bandwidths for the gridless KDE independently
spec.output.kernel = "uniform"
ok_uni = try EasleeConfig.validate_runspec(spec); true catch; false end
check("receptors accepted with any grid kernel", ok_uni)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
