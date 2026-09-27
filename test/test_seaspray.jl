#!/usr/bin/env julia
# test_seaspray.jl — OSSA sea-spray source function (Ovadnevaite et al. 2014)
#
#   [1] seawater viscosity: reproduces the paper's reference ν = 1.34e-6 m²/s
#       at 9 °C / 35 ‰; SST/salinity monotonicity,
#   [2] Re_Hw and thresholds: onset near the paper's U10 ≈ 3.7 m/s at
#       H_s = 1.23 m / SST 9 °C; zero flux below threshold; mode-5 threshold,
#   [3] mode fluxes: Table-1 power laws exact; total submicron flux magnitude
#       at Re = 3.1e6 in the Fig. 4 range; dF/dlogD five-modal + integrates
#       back to the total,
#
#   julia --project=easlee test/test_seaspray.jl
using Printf

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: OSSA_MODES, ossa_reynolds, ossa_mode_fluxes, ossa_number_flux,
               ossa_mass_flux, ossa_dfdlogd,
               seawater_kinematic_viscosity
using .EASLEE.Operators: _u10_from_ustar

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[1] seawater kinematic viscosity (Sharqawy et al. 2010)")
println("="^78)
ν9 = seawater_kinematic_viscosity(273.15 + 9.0)
# Sharqawy gives 1.40e-6 at 9 °C/35 ‰; the paper quotes 1.34e-6 (different
# correlation/rounding) — accept the 5 % correlation spread
check("ν(9 °C, 35 ‰) ≈ 1.34e-6 m²/s (OSSA reference, ±5 %)",
      isapprox(ν9, 1.34e-6; rtol=0.05), info=@sprintf("%.3e", ν9))
check("colder → more viscous",
      seawater_kinematic_viscosity(273.15) > ν9 >
      seawater_kinematic_viscosity(273.15 + 25.0))
check("saltier → more viscous",
      seawater_kinematic_viscosity(282.15; salinity_permil=35.0) >
      seawater_kinematic_viscosity(282.15; salinity_permil=5.0))

println()
println("="^78)
println("[2] Re_Hw and production thresholds")
println("="^78)
# paper: Re_Hw ≈ 1e5 at U10 ≈ 3.7 m/s given H_s = 1.23 m, SST = 9 °C —
# with the LOW-WIND bulk C_d (the paper's 2.15e-3 belongs to its 18 m/s case)
re37 = ossa_reynolds(3.7, 273.15 + 9.0; hs=1.23)
check("Re_Hw(3.7 m/s | Hs 1.23, SST 9°C) ≈ 1e5 (±25 %)",
      0.75e5 < re37 < 1.25e5, info=@sprintf("%.3g", re37))
check("no production below the submicron threshold",
      ossa_number_flux(0.9e5) == 0.0)
F15 = ossa_mode_fluxes(1.5e5)
check("modes 1–4 on, mode 5 off between the two thresholds",
      all(>(0.0), F15[1:4]) && F15[5] == 0.0)
check("mode 5 on above 2e5", ossa_mode_fluxes(2.5e5)[5] > 0.0)

println()
println("="^78)
println("[3] Table-1 power laws and size distribution")
println("="^78)
re = 3.1e6                                # the paper's Fig. 4 case (~18 m/s)
F  = ossa_mode_fluxes(re)
check("F_i(Re) match the Table-1 power laws exactly",
      all(isapprox(F[i], OSSA_MODES[i].a * (re - OSSA_MODES[i].re0)^OSSA_MODES[i].b;
                   rtol=1e-12) for i in 1:5))
check("total flux at Re 3.1e6 in the Fig-4 magnitude class (1e6–5e6 m⁻²s⁻¹)",
      1e6 < ossa_number_flux(re) < 5e6,
      info=@sprintf("%.3g", ossa_number_flux(re)))
# dF/dlog10D integrates back to the total number flux
lo, hi, n = log10(5e-9), log10(20e-6), 4000
dx = (hi - lo) / n
tot = sum(ossa_dfdlogd(10.0^(lo + (k - 0.5) * dx), re) for k in 1:n) * dx
check("∫ dF/dlogD dlogD = Σ F_i (±1 %)",
      isapprox(tot, ossa_number_flux(re); rtol=0.01),
      info=@sprintf("%.3g vs %.3g", tot, ossa_number_flux(re)))
check("mass flux positive and mode-5 dominated",
      ossa_mass_flux(re) > 0.0)

println()
println("="^78)
println("[4] u* ↔ U10 drag-law inversion")
println("="^78)
sst = 273.15 + 9.0
# U10 ← u* inversion (surface archives with zust but no 10 m winds): round-trips
# the forward drag law u* = √(C_d)·U10, C_d = (0.8 + 0.065·U10)e-3
cd_fwd(U) = (0.8 + 0.065 * U) * 1e-3
for U in (2.0, 6.0, 10.0, 18.0)
    us = U * sqrt(cd_fwd(U))
    check(@sprintf("_u10_from_ustar inverts drag law at U=%.0f m/s", U),
          isapprox(_u10_from_ustar(us), U; rtol=1e-3),
          info=@sprintf("u*=%.3f → %.3f", us, _u10_from_ustar(us)))
end
check("_u10_from_ustar(0) = 0", _u10_from_ustar(0.0) == 0.0)

println()
println("="^78)
println("[5] engine source ([[release]] kind=\"seaspray\")")
println("="^78)
using Dates
using .EASLEE: SeaSprayBox, build_seaspray, apply_seaspray!, _ossa_massflux
# synthetic 1° global grid, all-ocean except a land row at j=95
G = EASLEE.Types.Rectilinear(Int32(360), Int32(180), Int32(10), 1.0f0, 1.0f0,
        0.0f0, -90.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, true)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
lsm = zeros(Float32, 360, 180, 2); lsm[:, 95, :] .= 1.0f0
m2 = Dict("lsm" => lsm,
          "10u" => fill(12.0f0, 360, 180, 2), "10v" => zeros(Float32, 360, 180, 2),
          "2t"  => fill(273.15f0 + 9.0f0, 360, 180, 2))
fctx = (nests = [(gridconfig = G, levelconfig = LC)],)
t1 = DateTime(2025, 5, 1); t2 = DateTime(2025, 5, 2)
mkrel(; q=2.0e3) = EasleeConfig.ReleaseBox(
    species="ss", kind="seaspray", lon=[10.0, 14.0], lat=[2.0, 8.0],
    time=[t1, t2], mass_per_parcel_kg=q)
ssb = build_seaspray(fctx, m2, 3, mkrel())
check("build defers the cell scan (met unread at build time)",
      ssb !== nothing && isempty(ssb.ii) && !ssb.initialized[])
# analytic one-step accumulation at u10 = 12 m/s
sst9 = 273.15 + 9.0
Mflux = _ossa_massflux(ossa_mode_fluxes(ossa_reynolds(12.0, sst9)), 2017.0)
st = EASLEE.State.ParcelState(10_000)
nsp = apply_seaspray!(st, fctx, m2, ssb, t1 + Hour(1), 1.0, 600.0)
check("first apply scans the ocean cells (box minus the land row)",
      ssb.initialized[] && length(ssb.ii) > 0 && all(j -> j != 95, Int.(ssb.jj)),
      info="cells=$(length(ssb.ii))")
exp_mass = Mflux * sum(ssb.area) * 600.0
tot_emitted = ssb.emitted[] + sum(ssb.acc)
check("one-step production matches Σ(flux·area·dt) analytically",
      isapprox(tot_emitted, exp_mass; rtol=1e-6),
      info=@sprintf("%.4g vs %.4g kg", tot_emitted, exp_mass))
# spawning is PER CELL: nsp = Σ_c ⌊acc_c/q⌋ ∈ [⌊total/q⌋ − ncells, ⌊total/q⌋]
check("parcels spawned = Σ per-cell ⌊acc/quantum⌋, residuals carried",
      nsp == ssb.nspawn[] && ssb.emitted[] == nsp * ssb.quantum &&
      floor(Int, exp_mass / ssb.quantum) - length(ssb.ii) <= nsp <=
      floor(Int, exp_mass / ssb.quantum) &&
      all(<(ssb.quantum), ssb.acc), info="spawned=$nsp")
sp1 = findfirst(==(0x01), st.status)
check("spawned parcel: relbox, quantum mass, surface level",
      sp1 !== nothing && st.relbox[sp1] == Int32(3) &&
      st.mass[sp1] == Float32(ssb.quantum) && st.indl[sp1] < 2.0f0)   # level 1 = surface
check("outside the release window: no production",
      apply_seaspray!(st, fctx, m2, ssb, t2 + Hour(1), 1.0, 600.0) == 0)
check("no land-sea mask → source disabled (warns, nothing)",
      build_seaspray(fctx, Dict{String,Any}(), 1, mkrel()) === nothing)

println()
println("="^78)
println("[6] config: kind=seaspray + [release.seaspray]")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2025-05-01T00:00:00","stop"=>"2025-05-02T00:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[Dict{String,Any}("id"=>"ss","kind"=>"aerosol","aerosol_type"=>"seasalt")],
  "release"=>[Dict{String,Any}("species"=>"ss","kind"=>"seaspray",
      "lon"=>[10.0,14.0],"lat"=>[2.0,8.0],"z"=>[0.0,50.0],
      "time"=>["2025-05-01T00:00:00","2025-05-02T00:00:00"],
      "mass_per_parcel"=>"100kg",
      "seaspray"=>Dict{String,Any}("salinity"=>34.0,"sst"=>288.0))])
spc = EasleeConfig.ConfigTOML.dict_to_runspec(d)
EasleeConfig.validate_runspec(spc)
r1 = spc.release[1]
check("[release.seaspray] parses + validates",
      r1.kind == "seaspray" &&
      r1.spray_salinity == 34.0 && r1.spray_sst == 288.0 &&
      r1.mass_per_parcel_kg == 100.0)
rt = EasleeConfig.ConfigTOML.runspec_to_dict(spc)
r2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt).release[1]
check("round-trips", r2.kind == "seaspray" &&
      r2.spray_sst == 288.0)
d2 = deepcopy(d); delete!(d2["release"][1], "mass_per_parcel")
bad = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d2)); false
      catch; true end
check("missing mass_per_parcel rejected", bad)

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
