#!/usr/bin/env julia
# test_cov_kernel.jl — flow-evolved covariance output kernel (kernel="covariance"):
#   [1] P evolution: pure turbulent growth (∇u=0) → cxx = σ0² + 2D·t exact;
#       uniform shear ∂u/∂y=s → analytic FⁿP₀Fⁿᵀ (cxy = σ0²·st, nilpotent —
#       discrete product is EXACT); solid-body rotation → P stays isotropic;
#       init-on-first-touch; growth cap keeps P SPD; nest/status filters
#   [2] deposit (:covariance in _conc_nest!): EXACT conservation with σx≠σy
#       incl. edges, anisotropy visible in the footprint, single-cell limit
#       for a fresh (P=0→1 m) parcel, R(x) = (det P)^(1/4)
#   [3] receptors (usecov): full bivariate value = analytic Gaussian incl.
#       the cxy tilt (the term the separable grid deposit drops)
#   [4] config: [output.covariance] parse (units), round-trip, validation
#   [5] SHAPE mode (mode="shape", no double counting): deposit area = kNN h²
#       with P-anisotropy only (R(x) = h); receptor P rescaled to det = h⁴
#
#   julia --project=easlee -t 2 test/test_cov_kernel.jl
using Printf, Random, Statistics, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: update_parcel_cov!, cov_sigma_eff, _conc_nest!,
                         _out_cell_m, receptor_sample!, ReceptorSet
const ll = EASLEE.Operators.indyindxtolatlon

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# shared synthetic setup (same as test_adaptive_kernel.jl): 1° Rectilinear met
# grid, output = met grid, 10 levels, uniform zlev/rho
G = EASLEE.Types.Rectilinear(Int32(100), Int32(100), Int32(10), 1.0f0, 1.0f0,
        0.0f0, -50.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
NZ = 10
zlev = Array{Float32,4}(undef, 100, 100, NZ, 2)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * 500.0); end
rho = fill(1.0f0, 100, 100, NZ, 2)
zi  = Float32[500.0 * k for k in 1:7]
t0  = DateTime(2025, 1, 1)
my0 = 1.0 / Float64(G.dyconst)                 # metres per met index (y)
mx0 = 1.0 / Float64(G.dxconst)                 # metres per met index (x, equator)

function fresh_state(n)
    st = EASLEE.State.ParcelState(n)
    st.status .= 0x01; st.relbox .= Int32(1); st.nest .= Int32(1)
    st.indl .= 2.5f0; st.mass .= 1.0f0
    st.Dtime .= t0
    return st
end
zmet() = (Ux = zeros(Float32, 100, 100, NZ, 2), Vy = zeros(Float32, 100, 100, NZ, 2))

println("="^78)
println("[1] covariance evolution: P ← F·P·Fᵀ + 2D·dt")
println("="^78)

# ── pure turbulent growth (∇u = 0): cxx = σ0² + n·2σv²·TL·dt, cxy = 0 ────────
met0 = zmet()
st = fresh_state(4)
st.indx .= 51.0f0; st.indy .= 51.0f0
st.nest[3] = Int32(2)          # wrong nest → untouched
st.status[4] = 0x00            # unused slot → untouched
sv = 0.5; tl = 300.0; dt = 90.0; s0 = 100.0
nstep = 40
for _ in 1:nstep
    update_parcel_cov!(st, met0, G, 1.0, 1, dt; sigma0=s0, sigma_v=sv, tl=tl)
end
cxx_ref = s0^2 + nstep * 2 * sv^2 * tl * dt
check("pure growth: cxx = σ0² + 2σv²·TL·t exact",
      isapprox(Float64(st.cxx[1]), cxx_ref; rtol=1e-4) &&
      isapprox(Float64(st.cyy[1]), cxx_ref; rtol=1e-4) &&
      abs(st.cxy[1]) < 1e-3;
      info=@sprintf("(%.4g vs %.4g m²)", st.cxx[1], cxx_ref))
check("nest filter: other-nest parcel untouched", st.cxx[3] == 0.0f0)
check("status filter: unused slot untouched",     st.cxx[4] == 0.0f0)
check("σ_eff = (det P)^(1/4) = √cxx for isotropic P",
      isapprox(cov_sigma_eff(st, 1), sqrt(cxx_ref); rtol=1e-3))

# ── uniform shear ∂u/∂y = s (nilpotent → discrete product exact) ─────────────
metsh = zmet()
s_shear = 1.0e-4
for j in 1:100; metsh.Ux[:, j, :, :] .= Float32(s_shear * (j - 1) * my0); end
st2 = fresh_state(1)
st2.indx .= 51.0f0; st2.indy .= 51.0f0
s0b = 1000.0
nstep2 = 100
for _ in 1:nstep2
    update_parcel_cov!(st2, metsh, G, 1.0, 1, dt; sigma0=s0b, sigma_v=0.0, tl=tl)
end
stt = s_shear * nstep2 * dt                    # total shear s·t = 0.9
check("shear: cxy = σ0²·s·t (tilt grows linearly)",
      isapprox(Float64(st2.cxy[1]), s0b^2 * stt; rtol=2e-3);
      info=@sprintf("(%.4g vs %.4g)", st2.cxy[1], s0b^2 * stt))
check("shear: cxx = σ0²·(1+(s·t)²) (streamwise stretch)",
      isapprox(Float64(st2.cxx[1]), s0b^2 * (1 + stt^2); rtol=2e-3);
      info=@sprintf("(%.4g vs %.4g)", st2.cxx[1], s0b^2 * (1 + stt^2)))
check("shear: cyy = σ0² (cross-stream untouched)",
      isapprox(Float64(st2.cyy[1]), s0b^2; rtol=2e-3))
check("shear: P stays SPD (det > 0)",
      Float64(st2.cxx[1]) * Float64(st2.cyy[1]) - Float64(st2.cxy[1])^2 > 0.0)

# ── solid-body rotation: P stays isotropic (area ~ preserved) ────────────────
metro = zmet()
om = 1.0e-4
for j in 1:100; metro.Ux[:, j, :, :] .= Float32(-om * (j - 1) * my0); end
for i in 1:100; metro.Vy[i, :, :, :] .= Float32(+om * (i - 1) * mx0); end
st3 = fresh_state(1)
st3.indx .= 51.0f0; st3.indy .= 51.0f0     # lat 0 → mx = mx0 (cos = 1)
for _ in 1:nstep2
    update_parcel_cov!(st3, metro, G, 1.0, 1, dt; sigma0=s0b, sigma_v=0.0, tl=tl)
end
grw = (1 + (om * dt)^2)^nstep2             # O(dt²) spurious growth: +0.8 %
check("rotation: P stays isotropic (cxx ≈ cyy, cxy ≈ 0)",
      isapprox(Float64(st3.cxx[1]), Float64(st3.cyy[1]); rtol=1e-3) &&
      abs(st3.cxy[1]) < 1e-3 * st3.cxx[1];
      info=@sprintf("(cxx/cyy %.5f)", st3.cxx[1] / st3.cyy[1]))
check("rotation: only the O(dt²) area growth (+0.8 %)",
      isapprox(Float64(st3.cxx[1]), s0b^2 * grw; rtol=2e-3))

# ── growth cap keeps P SPD ────────────────────────────────────────────────────
st4 = fresh_state(1)
st4.indx .= 51.0f0; st4.indy .= 51.0f0
for _ in 1:400
    update_parcel_cov!(st4, metsh, G, 1.0, 1, 900.0;
                       sigma0=s0b, sigma_v=0.0, tl=tl, sigma_max=5.0e3)
end
det4 = Float64(st4.cxx[1]) * Float64(st4.cyy[1]) - Float64(st4.cxy[1])^2
check("growth cap: cxx ≤ σ_max², P stays SPD",
      st4.cxx[1] <= 5.0f3^2 * 1.001f0 && det4 > 0.0;
      info=@sprintf("(σx %.0f m, det %.3g)", sqrt(st4.cxx[1]), det4))

println()
println("="^78)
println("[2] covariance deposit: conservation, anisotropy, R(x)")
println("="^78)
Niout = 100; Njout = 100; Nk = 7
grid = zeros(Float32, Niout, Njout, Nk, 1, 1)
massarr = zeros(Float64, Niout, Njout, Nk, 1, 1)
resnum = zeros(Float64, Niout, Njout); resden = zeros(Float64, Niout, Njout)
dxm, dym = _out_cell_m(G)

# random cloud incl. edge parcels, log-spread anisotropic covariances
Random.seed!(11)
N3 = 5000
st5 = fresh_state(N3)
st5.indx .= Float32.(1.0 .+ 99.0 .* rand(N3))
st5.indy .= Float32.(1.0 .+ 99.0 .* rand(N3))
st5.cxx  .= Float32.((10 .^ (3 .+ 2 .* rand(N3))).^2)
st5.cyy  .= Float32.((10 .^ (3 .+ 2 .* rand(N3))).^2)
st5.cxy  .= Float32.(0.5 .* sqrt.(st5.cxx .* st5.cyy) .* (2 .* rand(N3) .- 1))
d3 = _conc_nest!(grid, massarr, st5, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                 zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
                 :covariance, resnum, resden)
check("EXACT conservation incl. edges (Σgrid = Σdeposited, 1e-12)",
      isapprox(sum(massarr), d3.deposited_mass; rtol=1e-12);
      info=@sprintf("(Σ %.6f, N_dep %d)", sum(massarr), d3.n_deposited))
check("R(x) accumulators populated and mass-consistent",
      isapprox(sum(resden), d3.deposited_mass; rtol=1e-12) &&
      count(>(0.0), resden) > 1000)

# anisotropic footprint: σx = 3 cells, σy = 1 cell → x-extent > y-extent
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
st6 = fresh_state(1)
st6.indx[1] = 50.5f0; st6.indy[1] = 60.5f0
st6.cxx[1] = Float32((3 * dxm)^2); st6.cyy[1] = Float32((1 * dym)^2)
st6.cxy[1] = 0.0f0
_conc_nest!(grid, massarr, st6, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
            :covariance, resnum, resden)
mcol = dropdims(sum(massarr; dims=(3, 4, 5)); dims=(3, 4, 5))
nx_ext = count(>(1e-6), mcol[:, 61])           # cells along x through the peak
ny_ext = count(>(1e-6), mcol[50, :])           # cells along y
check("anisotropy: σx=3 cells spreads wider in x than σy=1 in y",
      isapprox(sum(massarr), 1.0; rtol=1e-12) && nx_ext > 2 * ny_ext;
      info="($(nx_ext) x-cells vs $(ny_ext) y-cells)")
Rv = resnum[50, 61] / resden[50, 61]
Rref = (Float64(st6.cxx[1]) * Float64(st6.cyy[1]))^0.25
check("R(x) = (det P)^(1/4) (area-equivalent width)",
      isapprox(Rv, Rref; rtol=1e-6);
      info=@sprintf("(%.0f m vs %.0f m)", Rv, Rref))

# fresh parcel (P = 0 → clamps to 1 m) → single-cell deposit
grid .= 0; massarr .= 0
st7 = fresh_state(1)
st7.indx[1] = 50.3f0; st7.indy[1] = 60.7f0     # cxx = 0 (never evolved)
_conc_nest!(grid, massarr, st7, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0], :covariance)
check("fresh parcel (P=0) → single-cell deposit, mass exact",
      count(>(0.0), massarr) == 1 && isapprox(sum(massarr), 1.0; rtol=1e-12))

println()
println("="^78)
println("[3] receptors: full bivariate Gaussian (with the cxy tilt)")
println("="^78)
st8 = fresh_state(1)
st8.indx[1] = 51.0f0; st8.indy[1] = 51.0f0     # lat 0, lon 50; z = 750 m
pa = (15.0e3)^2; pd = (8.0e3)^2; pb = 0.6 * 15.0e3 * 8.0e3
st8.cxx[1] = Float32(pa); st8.cyy[1] = Float32(pd); st8.cxy[1] = Float32(pb)
la_p, lo_p = ll(G, Float64(st8.indy[1]), Float64(st8.indx[1]))
# two stations mirrored in lon: the tilt makes their values differ
st_lat = la_p + 0.05; st_lon1 = lo_p + 0.10; st_lon2 = lo_p - 0.10
rs = ReceptorSet(["S1", "S2"], [st_lon1, st_lon2], [st_lat, st_lat],
                 0.0, 2000.0; hmin=0.0, usecov=true)
receptor_sample!(rs, st8, G, LC, zlev, 1.0, 1)
function bivar(dx, dy)
    det = pa * pd - pb^2
    mah = (pd * dx^2 - 2 * pb * dx * dy + pa * dy^2) / det
    exp(-0.5 * mah) / (2π * sqrt(det)) / 2000.0
end
dy_s = (la_p - st_lat) * 111_320.0
dx1 = mod(lo_p - st_lon1 + 180.0, 360.0) - 180.0; dx1 *= 111_320.0 * cosd(st_lat)
dx2 = mod(lo_p - st_lon2 + 180.0, 360.0) - 180.0; dx2 *= 111_320.0 * cosd(st_lat)
c1_ref = bivar(dx1, dy_s); c2_ref = bivar(dx2, dy_s)
check("bivariate value matches analytic (station 1)",
      isapprox(rs.csum[1], c1_ref; rtol=1e-5);
      info=@sprintf("(%.4g vs %.4g)", rs.csum[1], c1_ref))
check("bivariate value matches analytic (station 2)",
      isapprox(rs.csum[2], c2_ref; rtol=1e-5))
check("the cxy tilt separates the mirrored stations",
      rs.csum[1] / rs.csum[2] > 1.5 || rs.csum[2] / rs.csum[1] > 1.5;
      info=@sprintf("(ratio %.2f, analytic %.2f)",
                    rs.csum[1] / rs.csum[2], c1_ref / c2_ref))
check("R(x_s) = (det P)^(1/4)",
      isapprox(rs.rnum[1] / rs.rden[1], (pa * pd - pb^2)^0.25; rtol=1e-5))
# 3σ ellipse cut: far station reports zero
rs2 = ReceptorSet(["F"], [lo_p + 2.0], [la_p], 0.0, 2000.0; hmin=0.0, usecov=true)
receptor_sample!(rs2, st8, G, LC, zlev, 1.0, 1)
check("beyond the 3σ ellipse → zero", rs2.csum[1] == 0.0)

println()
println("="^78)
println("[4] config: [output.covariance] parse, round-trip, validation")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "era5_conv.toml"); check_paths=false)
spec.output.kernel = "covariance"; spec.output.resolution = true
spec.output.cov_sigma0 = 200.0; spec.output.cov_sigma_v = 0.8
spec.output.cov_tl_s = 600.0; spec.output.cov_stride = 4
EasleeConfig.validate_runspec(spec; check_paths=false)
check("kernel=covariance + resolution=true validates", true)
dct = runspec_to_dict(spec)
cv = get(dct["output"], "covariance", Dict())
check("round-trip writes the [output.covariance] table",
      dct["output"]["kernel"] == "covariance" &&
      cv["sigma0"] == 200.0 && cv["sigma_v"] == 0.8 && cv["stride"] == 4)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(dct)
check("parse-back preserves the covariance parameters",
      spec2.output.kernel == "covariance" && spec2.output.cov_sigma0 == 200.0 &&
      spec2.output.cov_sigma_v == 0.8 && spec2.output.cov_tl_s == 600.0 &&
      spec2.output.cov_stride == 4)
# unit strings parse ("200m", "10min")
dct2 = deepcopy(dct)
dct2["output"]["covariance"] = Dict{String,Any}(
    "sigma0" => "200m", "sigma_v" => 0.8, "tl" => "10min", "stride" => 2)
spec3 = EasleeConfig.ConfigTOML.dict_to_runspec(dct2)
check("unit strings: sigma0=\"200m\", tl=\"10min\"",
      spec3.output.cov_sigma0 == 200.0 && spec3.output.cov_tl_s == 600.0)
spec.output.cov_stride = 0
bad1 = try EasleeConfig.validate_runspec(spec; check_paths=false); false catch; true end
check("stride = 0 rejected", bad1)
spec.output.cov_stride = 1; spec.output.cov_sigma0 = -1.0
bad2 = try EasleeConfig.validate_runspec(spec; check_paths=false); false catch; true end
check("sigma0 ≤ 0 rejected", bad2)

println()
println("="^78)
println("[5] shape mode: statistical area × physical anisotropy")
println("="^78)
# deposit: single parcel, h = 25 km, strong anisotropy Pxx = 16·Pyy (aspect 4)
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
st9 = fresh_state(1)
st9.indx[1] = 50.5f0; st9.indy[1] = 60.5f0                       # cell [50,60]
st9.hband[1] = Float32(dxm)                                      # h = 1 cell
st9.cxx[1] = Float32((40e3)^2); st9.cyy[1] = Float32((10e3)^2)   # aspect 4
st9.cxy[1] = 0.0f0
_conc_nest!(grid, massarr, st9, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
            :covshape, resnum, resden)
mcol9 = dropdims(sum(massarr; dims=(3, 4, 5)); dims=(3, 4, 5))
nx9 = count(>(1e-6), mcol9[:, 60]); ny9 = count(>(1e-6), mcol9[50, :])
check("shape deposit: mass exact, x-extent > y-extent (aspect 4)",
      isapprox(sum(massarr), 1.0; rtol=1e-12) && nx9 > 2 * ny9;
      info="($(nx9) x-cells vs $(ny9) y-cells)")
R9 = resnum[50, 60] / resden[50, 60]
check("shape deposit: R(x) = kNN h (area pinned to the spacing)",
      isapprox(R9, Float64(st9.hband[1]); rtol=1e-6);
      info=@sprintf("(%.0f m)", R9))
# isotropy limit: P isotropic → covshape ≡ adaptive footprint
grid .= 0; massarr .= 0
st9.cxx[1] = Float32((15e3)^2); st9.cyy[1] = Float32((15e3)^2)
d9a = _conc_nest!(grid, massarr, st9, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                  zlev, rho, LC, G, false, 60.0, Float32[1.0f0], :covshape)
m9a = copy(massarr); grid .= 0; massarr .= 0
d9b = _conc_nest!(grid, massarr, st9, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                  zlev, rho, LC, G, false, 60.0, Float32[1.0f0], :adaptive)
check("shape deposit ≡ adaptive when P is isotropic",
      isapprox(m9a, massarr; rtol=1e-10))

# receptor: shape rescaling P → P·h²/√det, R = h
st10 = fresh_state(1)
st10.indx[1] = 51.0f0; st10.indy[1] = 51.0f0
h10 = 12_000.0
pa0 = (30.0e3)^2; pd0 = (8.0e3)^2; pb0 = 0.5 * 30.0e3 * 8.0e3
st10.hband[1] = Float32(h10)
st10.cxx[1] = Float32(pa0); st10.cyy[1] = Float32(pd0); st10.cxy[1] = Float32(pb0)
la10, lo10 = ll(G, Float64(st10.indy[1]), Float64(st10.indx[1]))
rs10 = ReceptorSet(["S"], [lo10 + 0.08], [la10 + 0.03], 0.0, 2000.0;
                   hmin=0.0, usecov=true, shape=true)
receptor_sample!(rs10, st10, G, LC, zlev, 1.0, 1)
det0 = pa0 * pd0 - pb0^2
scl = h10^2 / sqrt(det0)
pas = pa0 * scl; pds = pd0 * scl; pbs = pb0 * scl
dy10 = (la10 - rs10.lat[1]) * 111_320.0
dx10 = (mod(lo10 - rs10.lon[1] + 180.0, 360.0) - 180.0) * 111_320.0 * cosd(rs10.lat[1])
mah10 = (pds * dx10^2 - 2pbs * dx10 * dy10 + pas * dy10^2) / (h10^4)
c10 = exp(-0.5 * mah10) / (2π * h10^2) / 2000.0
check("shape receptor: value = analytic with P·h²/√det",
      isapprox(rs10.csum[1], c10; rtol=1e-4);
      info=@sprintf("(%.4g vs %.4g)", rs10.csum[1], c10))
check("shape receptor: R(x_s) = h",
      isapprox(rs10.rnum[1] / rs10.rden[1], h10; rtol=1e-4))

# config: mode round-trip + validation
spec.output.cov_sigma0 = 50.0; spec.output.cov_stride = 1
spec.output.cov_mode = "shape"
EasleeConfig.validate_runspec(spec; check_paths=false)
dct3 = runspec_to_dict(spec)
spec4 = EasleeConfig.ConfigTOML.dict_to_runspec(dct3)
check("config: mode=\"shape\" round-trips",
      dct3["output"]["covariance"]["mode"] == "shape" &&
      spec4.output.cov_mode == "shape")
spec.output.cov_mode = "bogus"
bad3 = try EasleeConfig.validate_runspec(spec; check_paths=false); false catch; true end
check("config: unknown mode rejected", bad3)
spec.output.cov_mode = "physical"

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
