#!/usr/bin/env julia
# test_ice_physics.jl — transport-phase ice freezing + sedimentation
# ([species.ice], physics/ice_physics.jl):
#   [1] esat_ice: Murphy & Koop (2005) reference values,
#   [2] supersaturated stratospheric parcel loses EXACTLY the analytic
#       (q − q_si)·M_air·(1 − exp(−dt·v_t/h)) per step; frozen budget matches,
#   [3] gates: below rhi_crit no loss; warm boundary-layer parcel no loss;
#       repeated steps stop at RHi = rhi_crit (bounded hysteresis retention),
#   [4] q0-lazy air-mass assignment (no plume ensemble),
#   [5] config: [species.ice] parse / `ice = true` sugar / enabled=false off /
#       round-trip / validation rejects.
#
#   julia --project=easlee test/test_ice_physics.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: IcePhysics, apply_ice!, any_ice, esat_ice, ice_vterm
using .EASLEE.Operators: settling_velocity

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[1] esat_ice (Murphy & Koop 2005)")
println("="^78)
check("e_si(273.15 K) ≈ 611.15 Pa", isapprox(esat_ice(273.15), 611.15; rtol=2e-3),
      info=@sprintf("%.2f Pa", esat_ice(273.15)))
check("e_si(200 K) ≈ 0.1626 Pa", isapprox(esat_ice(200.0), 0.1626; rtol=2e-3),
      info=@sprintf("%.4f Pa", esat_ice(200.0)))
check("e_si(190 K) < e_si(200 K) < e_si(220 K) (monotone)",
      esat_ice(190.0) < esat_ice(200.0) < esat_ice(220.0))

# ── synthetic single-cell met: uniform (T, p) via θ/ρ ────────────────────────
const _Rd = 287.05; const _cp = 1004.0; const _P0 = 1.0e5
function mkctx(T::Float64, p::Float64)
    ρ = p / (_Rd * T)
    θ = T * (_P0 / p)^(_Rd / _cp)
    dims = (6, 6, 6, 2)
    theta = fill(Float32(θ), dims)
    rho   = fill(Float32(ρ), dims)
    (nests = [(metfields3D = (theta = theta, rho = rho),)],)
end
function mkstate(n; mass=1.0e6, airm=0.0)
    st = EASLEE.State.ParcelState(n)
    st.status .= 0x01; st.nest .= Int32(1); st.relbox .= Int32(1)
    st.indx .= 3.0f0; st.indy .= 3.0f0; st.indl .= 3.0f0
    st.mass .= Float32(mass); st.airm .= Float32(airm)
    return st
end
mkice(; rhi=1.5, r_eff=20e-6, h=500.0, q0=0.0) =
    IcePhysics([true], [Float32(rhi)], [Float32(r_eff)],
               [Float32(h)], [Float32(q0)], zeros(1))

println()
println("="^78)
println("[2] analytic single-step loss (195 K, 50 hPa)")
println("="^78)
T0 = 195.0; p0 = 5000.0
esi = esat_ice(T0); qsi = 0.622 * esi / (p0 - esi)
q   = 100.0 * qsi                       # deep supersaturation (fresh plume)
m0  = 1.0e6                             # 1 kt tagged water
Mair = m0 / q
ice = mkice()
st  = mkstate(1; mass=m0, airm=Mair)
dt  = 300.0
# v_t at the PARCEL's T, p — Stokes–Cunningham (matches apply_ice! exactly)
vt  = Float64(settling_velocity(2f0 * Float32(20e-6), 917f0, Float32(T0), Float32(p0)))
check("v_t(20 µm, 195 K, 50 hPa) is cm/s-scale (not the 1.7 m/s plume closure)",
      0.005 <= vt <= 0.30, info=@sprintf("%.3f m/s = %.1f km/day", vt, vt * 86.4))
check("slip correction: v_t grows toward lower pressure",
      settling_velocity(4f-5, 917f0, Float32(T0), 2000f0) >
      settling_velocity(4f-5, 917f0, Float32(T0), 20000f0))
removed = apply_ice!(st, mkctx(T0, p0), ice, dt, [1.0])
dm_ana = (q - qsi) * Mair * (-expm1(-dt * vt / 500.0))
check("removed mass = analytic (q−q_si)·M_air·(1−e^{−dt·v_t/h})",
      isapprox(removed, dm_ana; rtol=2e-4),
      info=@sprintf("%.4g vs %.4g kg", removed, dm_ana))
check("parcel mass decremented consistently",
      isapprox(Float64(st.mass[1]), m0 - dm_ana; rtol=2e-4))
check("frozen accumulator = removed", isapprox(ice.frozen[1], removed; rtol=1e-12))
check("any_ice true / off-box skipped", any_ice(ice))

println()
println("="^78)
println("[3] gates + asymptote")
println("="^78)
# below the nucleation threshold: RHi = 1.2 < 1.5 → untouched
st2 = mkstate(1; mass=1.2 * qsi * Mair, airm=Mair)
r2  = apply_ice!(st2, mkctx(T0, p0), mkice(), dt, [1.0])
check("RHi 1.2 < rhi_crit 1.5: no loss", r2 == 0.0 && st2.mass[1] == Float32(1.2 * qsi * Mair))
# warm BL parcel: q_si huge → subsaturated → untouched
st3 = mkstate(1; mass=m0, airm=Mair)
r3  = apply_ice!(st3, mkctx(288.0, 90000.0), mkice(), dt, [1.0])
check("288 K / 900 hPa parcel: no loss", r3 == 0.0)
# repeated steps: the loss GATES at rhi_crit — the last firing step overshoots
# below it by at most one (q−q_si)·fsh increment, so the retained mixing ratio
# lands in the [q_si, rhi_crit·q_si] band and never goes subsaturated
st4 = mkstate(1; mass=m0, airm=Mair); ic4 = mkice()
for _ in 1:20000
    apply_ice!(st4, mkctx(T0, p0), ic4, dt, [1.0])
end
qf = Float64(st4.mass[1]) / Mair
check("asymptote: q_si ≤ retained q ≤ rhi_crit·q_si (never subsaturated)",
      qsi * 0.999 <= qf <= 1.5 * qsi * 1.001,
      info=@sprintf("q_f/q_si = %.3f", qf / qsi))
check("budget closes: frozen = m0 − m_f",
      isapprox(ic4.frozen[1], m0 - Float64(st4.mass[1]); rtol=1e-6))

println()
println("="^78)
println("[4] q0-lazy air mass (no plume ensemble)")
println("="^78)
q0 = 50.0 * qsi
st5 = mkstate(1; mass=m0, airm=0.0)
ic5 = mkice(q0=q0)
r5  = apply_ice!(st5, mkctx(T0, p0), ic5, dt, [1.0])
check("airm assigned = m/q0 on first touch",
      isapprox(Float64(st5.airm[1]), m0 / q0; rtol=1e-5))
check("loss follows q0 mixing ratio", r5 > 0.0)
st6 = mkstate(1; mass=m0, airm=0.0)          # no q0, no plume → inert
r6  = apply_ice!(st6, mkctx(T0, p0), mkice(q0=0.0), dt, [1.0])
check("no airm + no q0: parcel untouched", r6 == 0.0)

println()
println("="^78)
println("[6] ice co-scavenging ([species.ice_scav])")
println("="^78)
# boxes: 1 = icing water, 2 = partner (uptake 0.6), 3 = partner in a dry cell
ic6 = IcePhysics([true, false, false], Float32[1.5, 1.5, 1.5],
                 Float32[20e-6, 0, 0], Float32[500, 500, 500],
                 Float32[0, 0, 0], zeros(3),
                 Float32[0, 0.6, 0.6], zeros(3))
st6b = mkstate(4; mass=m0, airm=Mair)
st6b.relbox .= Int32[1, 2, 3, 1]
st6b.indx[3] = 5.0f0                        # partner 3: different cell → untouched
st6b.airm[2] = 0f0; st6b.airm[3] = 0f0      # partners need no airm
mp0 = Float64(st6b.mass[2])
r6b = apply_ice!(st6b, mkctx(T0, p0), ic6, dt, [1.0])
# cell shed fraction: two water parcels (1 and 4) in the cell, both shed dm_ana
fsh6 = 2 * dm_ana / (2 * m0)
dm_p = mp0 * 0.6 * fsh6
check("co-located partner loses m·uptake·f_shed",
      isapprox(Float64(mp0 - st6b.mass[2]), dm_p; rtol=5e-4),
      info=@sprintf("%.4g vs %.4g kg", mp0 - Float64(st6b.mass[2]), dm_p))
check("partner in an ice-free cell untouched", st6b.mass[3] == Float32(m0))
check("scavenged accumulator = partner loss",
      isapprox(ic6.scavenged[2], mp0 - Float64(st6b.mass[2]); rtol=1e-6) &&
      ic6.scavenged[3] == 0.0)
check("water budget separate from scavenged budget",
      isapprox(ic6.frozen[1], r6b; rtol=1e-12) && ic6.frozen[2] == 0.0)
check("any_ice true for partner-only config",
      any_ice(IcePhysics([false], Float32[1.5], Float32[0], Float32[500],
                         Float32[0], zeros(1), Float32[0.5], zeros(1))))

println()
println("="^78)
println("[5] config: [species.ice] parse / round-trip / validation")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-15T06:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[Dict{String,Any}("id"=>"H2O_volc","kind"=>"gas",
      "ice"=>Dict{String,Any}("rhi_crit"=>1.4,"r_eff"=>"30um","h_layer"=>"1km","q0"=>0.002))],
  "release"=>[Dict{String,Any}("species"=>"H2O_volc","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>100,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d).species[1]
check("[species.ice] table parses",
      sp.ice && sp.ice_rhi_crit == 1.4 && isapprox(sp.ice_r_eff_m, 30e-6) &&
      sp.ice_h_layer_m == 1000.0 && sp.ice_q0 == 0.002)
d2 = deepcopy(d); d2["species"][1]["ice"] = true
check("`ice = true` sugar enables with defaults",
      EasleeConfig.ConfigTOML.dict_to_runspec(d2).species[1].ice)
d3 = deepcopy(d); d3["species"][1]["ice"] = Dict{String,Any}("enabled"=>false)
check("enabled = false disables (the GUI default)",
      !EasleeConfig.ConfigTOML.dict_to_runspec(d3).species[1].ice)
d3b = deepcopy(d); delete!(d3b["species"][1], "ice")
check("NO ice key → ice stays OFF (the empty-Dict-default regression)",
      !EasleeConfig.ConfigTOML.dict_to_runspec(d3b).species[1].ice)
rt = EasleeConfig.ConfigTOML.runspec_to_dict(EasleeConfig.ConfigTOML.dict_to_runspec(d))
ic = rt["species"][1]["ice"]
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt).species[1]
check("round-trip preserves the ice block",
      ic["enabled"] == true && sp2.ice && sp2.ice_rhi_crit == 1.4 &&
      isapprox(sp2.ice_r_eff_m, 30e-6) && sp2.ice_q0 == 0.002)
d4 = deepcopy(d); d4["species"][1]["ice"]["rhi_crit"] = 0.8
bad = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d4)); false
      catch; true end
check("rhi_crit < 1 rejected", bad)
# ice_scav partner: parse + validation
d5 = deepcopy(d)
push!(d5["species"], Dict{String,Any}("id"=>"SO2","kind"=>"gas",
    "ice_scav"=>Dict{String,Any}("of"=>"H2O_volc","uptake"=>0.3)))
sp5 = EasleeConfig.ConfigTOML.dict_to_runspec(d5).species[2]
EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d5))
check("[species.ice_scav] parses and validates",
      sp5.ice_scav_of == "H2O_volc" && sp5.ice_scav_uptake == 0.3)
d6 = deepcopy(d5); d6["species"][2]["ice_scav"]["of"] = "NOPE"
bad6 = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d6)); false
       catch; true end
check("ice_scav.of must name an ice-enabled species", bad6)
d7 = deepcopy(d5); d7["species"][2]["ice_scav"]["uptake"] = 1.5
bad7 = try EasleeConfig.validate_runspec(EasleeConfig.ConfigTOML.dict_to_runspec(d7)); false
       catch; true end
check("uptake > 1 rejected", bad7)
rt5 = EasleeConfig.ConfigTOML.runspec_to_dict(EasleeConfig.ConfigTOML.dict_to_runspec(d5))
sp5b = EasleeConfig.ConfigTOML.dict_to_runspec(rt5).species[2]
check("ice_scav round-trips",
      sp5b.ice_scav_of == "H2O_volc" && sp5b.ice_scav_uptake == 0.3)

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
