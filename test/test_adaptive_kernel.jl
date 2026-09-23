#!/usr/bin/env julia
# test_adaptive_kernel.jl — adaptive output kernel + R(x) (ADAPTIVE_KERNEL A+B):
#   [1] erf approximation accuracy (Abramowitz–Stegun 7.1.26),
#   [2] kNN bandwidth: lattice spacing recovered, N^(−1/2) scaling, point-
#       cluster floor (the 3-h rule replacement),
#   [3] adaptive deposit: EXACT mass conservation incl. domain edges
#       (truncated-support renormalisation), single-cell limit at h→floor,
#       R(x) accumulators = h for a single-bandwidth population,
#   [4] uniform kernel unchanged (comparison reference), edge loss documented,
#   [5] output.kernel/resolution config round-trip + validation.
#
#   julia --project=easlee -t 2 test/test_adaptive_kernel.jl
using Printf, Random, Statistics, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: _erf_as, _Φ, update_bandwidths!, _conc_nest!, _out_cell_m

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# shared synthetic setup: global-ish 1° Rectilinear met grid, output = met grid
# (coeffx=coeffy=1, offset 0), 10 levels, uniform zlev/rho
G = EASLEE.Types.Rectilinear(Int32(100), Int32(100), Int32(10), 1.0f0, 1.0f0,
        0.0f0, -50.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
NZ = 10
zlev = Array{Float32,4}(undef, 100, 100, NZ, 2)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * 500.0); end
rho = fill(1.0f0, 100, 100, NZ, 2)
zi  = Float32[500.0 * k for k in 1:7]                      # 7 output layers
t0  = DateTime(2025, 1, 1)

function fresh_state(n)
    st = EASLEE.State.ParcelState(n)
    st.status .= 0x01; st.relbox .= Int32(1); st.nest .= Int32(1)
    st.indl .= 2.5f0; st.mass .= 1.0f0
    st.Dtime .= t0
    return st
end

println("="^78)
println("[1] erf approximation")
println("="^78)
check("erf accuracy (A-S 7.1.26: |err| < 2e-7)",
      abs(_erf_as(0.0)) < 2e-9 &&
      abs(_erf_as(1.0) - 0.8427007929497149) < 2e-7 &&
      abs(_erf_as(-1.0) + 0.8427007929497149) < 2e-7 &&
      abs(_erf_as(2.0) - 0.9953222650189527) < 2e-7)
check("Φ is a CDF (0 at −∞, 1 at +∞, ½ at 0)",
      _Φ(-30.0) < 1e-9 && abs(_Φ(0.0) - 0.5) < 1e-9 && _Φ(30.0) > 1 - 1e-9)

println()
println("="^78)
println("[2] kNN bandwidth estimator")
println("="^78)
# 30×30 lattice, spacing d index units → d·dx_m metres
dxm, dym = _out_cell_m(G)
d_idx = 0.5
st = fresh_state(900)
for (n, (a, b)) in enumerate(Iterators.product(1:30, 1:30))
    st.indx[n] = 20.0f0 + Float32(a * d_idx)
    st.indy[n] = 20.0f0 + Float32(b * d_idx)
end
update_bandwidths!(st, G, 1)
d_m = d_idx * dxm                       # lattice spacing in metres (x); y differs via dym
interior = [n for (n, (a, b)) in enumerate(Iterators.product(1:30, 1:30))
            if 5 <= a <= 25 && 5 <= b <= 25]
h_int = Float64.(st.hband[interior])
# analytic RMS of the 12 nearest square-lattice sites: 4 at d, 4 at d√2, 4 at 2d
d_eff = sqrt(0.5 * (d_idx * dxm)^2 + 0.5 * (d_idx * dym)^2)  # x/y anisotropy mix
h_ref = sqrt((4 * 1 + 4 * 2 + 4 * 4) / 12) * d_eff
check("lattice: interior h uniform (std < 5 %)",
      std(h_int) / mean(h_int) < 0.05;
      info=@sprintf("(cv %.3f)", std(h_int) / mean(h_int)))
check("lattice: h ≈ RMS of 12 nearest sites (±15 %)",
      0.85 * h_ref <= mean(h_int) <= 1.15 * h_ref;
      info=@sprintf("(%.0f m vs %.0f m)", mean(h_int), h_ref))
# N^(−1/2): halve the spacing (×4 parcels in the same area) → h halves
st2 = fresh_state(3600)
for (n, (a, b)) in enumerate(Iterators.product(1:60, 1:60))
    st2.indx[n] = 20.0f0 + Float32(a * d_idx / 2)
    st2.indy[n] = 20.0f0 + Float32(b * d_idx / 2)
end
update_bandwidths!(st2, G, 1)
interior2 = [n for (n, (a, b)) in enumerate(Iterators.product(1:60, 1:60))
             if 10 <= a <= 50 && 10 <= b <= 50]
r = mean(Float64.(st2.hband[interior2])) / mean(h_int)
check("N ×4 → h halves (N^(−1/2) law)", 0.45 <= r <= 0.55;
      info=@sprintf("(ratio %.3f)", r))
# point cluster: all parcels at one spot → h = 1 m floor (3-h rule replacement)
st3 = fresh_state(50)
st3.indx .= 40.0f0; st3.indy .= 40.0f0
update_bandwidths!(st3, G, 1)
check("point release: h collapses to the 1 m floor",
      all(st3.hband[1:50] .== 1.0f0);
      info=@sprintf("(max %.1f m)", maximum(st3.hband[1:50])))

println()
println("="^78)
println("[3] adaptive deposit: conservation, single-cell limit, R(x)")
println("="^78)
Niout = 100; Njout = 100; Nk = 7
grid = zeros(Float32, Niout, Njout, Nk, 1, 1)
massarr = zeros(Float64, Niout, Njout, Nk, 1, 1)
resnum = zeros(Float64, Niout, Njout); resden = zeros(Float64, Niout, Njout)
# random cloud incl. parcels hanging over the domain edge, mixed bandwidths
Random.seed!(11)
N3 = 5000
st4 = fresh_state(N3)
st4.indx .= Float32.(1.0 .+ 99.0 .* rand(N3))
st4.indy .= Float32.(1.0 .+ 99.0 .* rand(N3))
st4.hband .= Float32.(10 .^ (3 .+ 2 .* rand(N3)))          # 1–100 km bandwidths
d3 = _conc_nest!(grid, massarr, st4, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                 zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
                 :adaptive, resnum, resden)
tot_dep = sum(massarr)
check("EXACT conservation incl. edges (Σgrid = Σdeposited, 1e-12)",
      isapprox(tot_dep, d3.deposited_mass; rtol=1e-12);
      info=@sprintf("(Σ %.6f vs %.6f, N_dep %d)", tot_dep, d3.deposited_mass, d3.n_deposited))
check("R(x) accumulators populated and mass-consistent",
      isapprox(sum(resden), d3.deposited_mass; rtol=1e-12) &&
      count(>(0.0), resden) > 1000)
# single-cell limit: one parcel, h at floor, at a cell centre
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
st5 = fresh_state(1)
st5.indx[1] = 50.3f0; st5.indy[1] = 60.7f0; st5.hband[1] = 1.0f0
_conc_nest!(grid, massarr, st5, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
            :adaptive, resnum, resden)
nz_cells = count(>(0.0), massarr)
check("h = floor → single-cell deposition (no 3-h rule needed)",
      nz_cells == 1 && isapprox(sum(massarr), 1.0; rtol=1e-12);
      info="($nz_cells cells)")
# R(x) equals h for a uniform-h population
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
st6 = fresh_state(200)
st6.indx .= Float32.(30.0 .+ 20.0 .* rand(200)); st6.indy .= Float32.(30.0 .+ 20.0 .* rand(200))
st6.hband .= 25_000.0f0
_conc_nest!(grid, massarr, st6, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
            :adaptive, resnum, resden)
Rvals = [resnum[i] / resden[i] for i in eachindex(resnum) if resden[i] > 0]
check("R(x) = h for a single-bandwidth population",
      all(v -> isapprox(v, 25_000.0; rtol=1e-6), Rvals))

println()
println("="^78)
println("[4] uniform kernel path intact (comparison reference)")
println("="^78)
grid .= 0; massarr .= 0
d4 = _conc_nest!(grid, massarr, st4, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                 zlev, rho, LC, G, false, 60.0, Float32[1.0f0], :uniform)
check("uniform kernel deposits (bilinear 4-cell)",
      d4.n_deposited == d3.n_deposited && sum(massarr) > 0.99 * d4.deposited_mass * 0.9)
grid .= 0; massarr .= 0
d5 = _conc_nest!(grid, massarr, st4, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                 zlev, rho, LC, G, false, 60.0, Float32[1.0f0], :none)
check("none = nearest cell, mass exactly counted",
      isapprox(sum(massarr), d5.deposited_mass; rtol=1e-12))

println()
println("="^78)
println("[5] config round-trip + validation")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "era5_conv.toml"); check_paths=false)
spec.output.kernel = "adaptive"; spec.output.resolution = true
EasleeConfig.validate_runspec(spec; check_paths=false)
dct = runspec_to_dict(spec)
check("writes kernel=adaptive + resolution",
      dct["output"]["kernel"] == "adaptive" && dct["output"]["resolution"] === true)
spec2 = EasleeConfig.ConfigTOML.dict_to_runspec(dct)
check("parse-back preserves both",
      spec2.output.kernel == "adaptive" && spec2.output.resolution)
spec.output.kernel = "uniform"
bad1 = try EasleeConfig.validate_runspec(spec; check_paths=false); false catch; true end
check("resolution without adaptive kernel rejected", bad1)
spec.output.resolution = false; spec.output.kernel = "bogus"
bad2 = try EasleeConfig.validate_runspec(spec; check_paths=false); false catch; true end
check("unknown kernel rejected", bad2)

println()
println("="^78)
println("[6] per-population bandwidths (task #35)")
println("="^78)
# two CO-LOCATED populations on the same patch: a dense plume (release 1,
# 3600 parcels at spacing d/4) and a sparse background (release 2, 100
# parcels at spacing 3d).  Pooled kNN would give the background parcels the
# PLUME spacing; per-population must keep them apart.
using .EASLEE.Operators: set_output_popmap!
stp = fresh_state(3700)
for (n, (a, b)) in enumerate(Iterators.product(1:60, 1:60))
    stp.indx[n] = 25.0f0 + Float32(a * d_idx / 4)
    stp.indy[n] = 25.0f0 + Float32(b * d_idx / 4)
    stp.relbox[n] = Int32(1)
end
for (m, (a, b)) in enumerate(Iterators.product(1:10, 1:10))
    n = 3600 + m
    stp.indx[n] = 25.0f0 + Float32(a * 3d_idx)
    stp.indy[n] = 25.0f0 + Float32(b * 3d_idx)
    stp.relbox[n] = Int32(2)
end
# pooled reference
set_output_popmap!(Int[])
update_bandwidths!(stp, G, 1)
h_bg_pooled = mean(Float64.(stp.hband[3601:3700]))
# per-population: boxes 1,2 → populations 1,2
set_output_popmap!([1, 2])
update_bandwidths!(stp, G, 1)
h_bg = mean(Float64.(stp.hband[3601:3700]))
h_pl = mean(Float64.(stp.hband[1:3600]))
# (background parcels inside the dense patch collapse to the plume spacing;
# those outside keep their own — the pooled mean lands well below per-pop)
check("pooled kNN contaminates the background h (the bug)",
      h_bg_pooled < 0.6 * h_bg;
      info=@sprintf("pooled %.0f m vs per-pop %.0f m", h_bg_pooled, h_bg))
check("per-population: background keeps its own spacing",
      h_bg > 5 * h_pl;
      info=@sprintf("bg %.0f m, plume %.0f m", h_bg, h_pl))
# a population never mapped (relbox=0) still gets a bandwidth
stp.relbox[3700] = Int32(0)
update_bandwidths!(stp, G, 1)
check("unmapped relbox → population 0, still gets h", stp.hband[3700] > 0.0f0)
set_output_popmap!(Int[])       # restore for any later sections

println()
println("="^78)
println("[7] physical bandwidth floor (bandwidth_floor = \"dispersion\")")
println("="^78)
using .EASLEE.Operators: _hfloor_taylor, receptor_sample!
# Taylor limits: ballistic σ ≈ σv·t (t ≪ T_L), diffusive σ ≈ σv·√(2 T_L t)
sv = 0.5; tl = 300.0
check("Taylor floor: ballistic limit σ→σv·t",
      isapprox(_hfloor_taylor(3.0, sv, tl), sv * 3.0; rtol=0.01) &&
      _hfloor_taylor(0.0, sv, tl) == 0.0)
check("Taylor floor: diffusive limit σ→σv·√(2·T_L·t)",
      isapprox(_hfloor_taylor(3.0e5, sv, tl), sv * sqrt(2 * tl * 3.0e5); rtol=0.01);
      info=@sprintf("(1 h → %.0f m, 24 h → %.0f m)",
                    _hfloor_taylor(3600.0, sv, tl), _hfloor_taylor(86400.0, sv, tl)))
# kNN pass: a dense point cluster aged 1 h gets the floor, not the 1 m spacing
stf = fresh_state(50)
stf.indx .= 40.0f0; stf.indy .= 40.0f0
stf.Dtimeinit .= t0 - Hour(1)
hf1h = _hfloor_taylor(3600.0, sv, tl)
update_bandwidths!(stf, G, 1; floor_sv=sv, floor_tl=tl)
check("kNN pass: point cluster floored at σ(age = 1 h)",
      all(isapprox.(Float64.(stf.hband[1:50]), hf1h; rtol=1e-4));
      info=@sprintf("(h %.0f m)", hf1h))
update_bandwidths!(stf, G, 1)                       # floor off → old behaviour
check("floor off: point cluster back to the 1 m spacing floor",
      all(stf.hband[1:50] .== 1.0f0))
# deposit: a parcel born AFTER the last kNN refresh (hband = 0) still gets a
# physical footprint at deposit time; conservation stays exact
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
stg = fresh_state(1)
stg.indx[1] = 50.5f0; stg.indy[1] = 60.5f0          # cell centre
stg.Dtimeinit[1] = t0 - Hour(12)                    # σ ≈ 2.5 km ≫ 1° cell? no: ≪ cell
# make the floor span cells: age 12 h → σ = 2546 m vs 111 km cells — use a
# finer-grid config instead: reuse G but scale the floor via a large T_L
svL = 60.0                                          # exaggerated σv to span cells
dL  = _conc_nest!(grid, massarr, stg, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
                  zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
                  :adaptive, resnum, resden, svL, tl)
nzf = count(>(0.0), massarr)
check("deposit: fresh parcel (hband=0) spreads by the age floor",
      nzf > 1 && isapprox(sum(massarr), 1.0; rtol=1e-12);
      info=@sprintf("(%d cells, σ %.0f km)", nzf, _hfloor_taylor(43200.0, svL, tl)/1e3))
grid .= 0; massarr .= 0; fill!(resnum, 0.0); fill!(resden, 0.0)
_conc_nest!(grid, massarr, stg, t0, 1.0f0, Niout, Njout, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 60.0, Float32[1.0f0],
            :adaptive, resnum, resden)
check("deposit: floor off → fresh parcel stays single-cell",
      count(>(0.0), massarr) == 1)
# receptor: floor bounds the single-parcel point-evaluation spike
la1, lo1 = EASLEE.Operators.indyindxtolatlon(G, 40.0, 40.0)
rsF = EASLEE.Operators.ReceptorSet(["X"], [Float64(lo1)], [Float64(la1)], 0.0, 500.0;
                                   fsv=sv, ftl=tl)
rs0 = EASLEE.Operators.ReceptorSet(["X"], [Float64(lo1)], [Float64(la1)], 0.0, 500.0)
sth = fresh_state(1)
sth.indx[1] = 40.0f0; sth.indy[1] = 40.0f0; sth.indl[1] = 1.0f0
sth.Dtimeinit[1] = t0 - Hour(1); sth.hband[1] = 60.0f0     # the ETEX-8021 case
receptor_sample!(rsF, sth, G, LC, zlev, 1.0, 1)
receptor_sample!(rs0, sth, G, LC, zlev, 1.0, 1)
c_ratio = rs0.csum[1] / rsF.csum[1]
check("receptor: 1-h-old 60 m spike bounded by the floor",
      isapprox(c_ratio, (hf1h / 60.0)^2; rtol=1e-3);
      info=@sprintf("(spike ×%.0f tamed)", c_ratio))
# config round-trip
spec7 = read_runspec(joinpath(REPO, "configs", "era5_conv.toml"); check_paths=false)
spec7.output.kernel = "adaptive"; spec7.output.bandwidth_floor = "dispersion"
spec7.output.floor_sigmav = 0.7; spec7.output.floor_tl_s = 600.0
d7 = runspec_to_dict(spec7)
spec7b = EasleeConfig.ConfigTOML.dict_to_runspec(d7)
check("config round-trip: bandwidth_floor + σv + T_L preserved",
      spec7b.output.bandwidth_floor == "dispersion" &&
      spec7b.output.floor_sigmav == 0.7 && spec7b.output.floor_tl_s == 600.0)
spec7.output.bandwidth_floor = "banana"
bad3 = try EasleeConfig.validate_runspec(spec7; check_paths=false); false catch; true end
check("unknown bandwidth_floor rejected", bad3)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
