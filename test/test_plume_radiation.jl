#!/usr/bin/env julia
# test_plume_radiation.jl — in-plume LW cooling ([physics.plume_radiation]):
#   [1] analytic one-step cooling: ΔT = −c_rad·(q − q_bg)·dt from the parcel's
#       ice-closure mixing ratio; below-background and airm-less parcels
#       untouched; freed slots reset,
#   [2] descent feedback: dz = dθ/(∂θ/∂z) exact on a linear θ(z) profile;
#       feedback=false leaves indl untouched,
#   [3] config: [physics.plume_radiation] parse / round-trip.
#
#   julia --project=easlee test/test_plume_radiation.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: PlumeRadiation, apply_plume_radiation!, any_prad

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# synthetic met: 10 levels, z = (k−1)·1000 m, θ = 500 + 0.02·z (∂θ/∂z = 0.02),
# uniform p via rho — the module reads theta/rho/zlev/puv
const NZ = 10
dims = (8, 8, NZ, 2)
theta = zeros(Float32, dims); zlev = zeros(Float32, dims)
for k in 1:NZ
    zlev[:, :, k, :] .= Float32((k - 1) * 1000.0)
    theta[:, :, k, :] .= Float32(500.0 + 0.02 * (k - 1) * 1000.0)
end
p0 = 5000.0                                    # 50 hPa
rho  = fill(Float32(p0 / (287.05 * 210.0)), dims)
puv  = fill(Float32(p0), dims)
m3 = (theta = theta, rho = rho, zlev = zlev, puv = puv)
LC = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
fctx = (nests = [(metfields3D = m3, levelconfig = LC,)],
        runconfig = (direction = 1,))          # forward; see [9] for backward
# same context integrating BACKWARD — the radiative tendency must reverse sign
bctx = (nests = [(metfields3D = m3, levelconfig = LC,)],
        runconfig = (direction = -1,))

MW = 28.9647 / 18.0153
c_rad_day = 0.002                              # K/day per ppmv
crad = c_rad_day / 86400.0
q_bg = 5.0
mkpr(fb) = PlumeRadiation([true, false], crad, q_bg, fb, zeros(Float32, 100))
function mkstate()
    st = EASLEE.State.ParcelState(100)
    st.status[1:4] .= 0x01
    st.nest[1:4] .= Int32(1)
    st.relbox[1:4] .= Int32[1, 1, 2, 1]        # 3 in the plume box, 1 elsewhere
    st.indx[1:4] .= 4.0f0; st.indy[1:4] .= 4.0f0; st.indl[1:4] .= 5.0f0
    # parcel 1: 200 ppmv; parcel 2: below background; parcel 4: no airm
    q1 = 200.0 / (MW * 1e6); st.mass[1] = 1.0f6; st.airm[1] = Float32(1.0e6 / q1)
    q2 = 2.0   / (MW * 1e6); st.mass[2] = 1.0f6; st.airm[2] = Float32(1.0e6 / q2)
    st.mass[3] = 1.0f6; st.airm[3] = Float32(1.0e6 / q1)
    st.mass[4] = 1.0f6; st.airm[4] = 0.0f0
    return st
end

println("="^78)
println("[1] analytic one-step cooling")
println("="^78)
dt = 600.0
st = mkstate(); pr = mkpr(false)
pr.trad[50] = 9.0f0                            # stale value in a freed slot
n, mdt = apply_plume_radiation!(st, fctx, pr, dt, [1.0])
dT_exp = -crad * (200.0 - q_bg) * dt
check("ΔT = −c_rad·(q−q_bg)·dt (200 ppmv parcel)",
      isapprox(Float64(pr.trad[1]), dT_exp; rtol=1e-5),
      info=@sprintf("%.4g vs %.4g K", pr.trad[1], dT_exp))
check("below-background parcel untouched", pr.trad[2] == 0.0f0)
check("off-box parcel untouched", pr.trad[3] == 0.0f0)
check("airm-less parcel skipped", pr.trad[4] == 0.0f0 && n == 1)
check("mass-weighted budget", isapprox(mdt, 1.0e6 * dT_exp; rtol=1e-5))
check("freed slot reset", pr.trad[50] == 0.0f0)
check("no feedback: indl untouched", st.indl[1] == 5.0f0)

println()
println("="^78)
println("[2] descent feedback (linear θ profile)")
println("="^78)
st2 = mkstate(); pr2 = mkpr(true)
apply_plume_radiation!(st2, fctx, pr2, dt, [1.0])
# dθ = ΔT·(p0/p)^κ; dz = dθ/0.02; Δindl = dz/1000 m per level
dθ = dT_exp * (1.0e5 / p0)^0.2857
dz_exp = dθ / 0.02
check("parcel descends by dθ/(∂θ/∂z)",
      isapprox(Float64(st2.indl[1] - 5.0f0) * 1000.0, dz_exp; rtol=0.02),
      info=@sprintf("%.2f vs %.2f m", (st2.indl[1] - 5.0) * 1000.0, dz_exp))
check("descent is downward and small per step", dz_exp < 0 && abs(dz_exp) < 50.0)

println()
println("="^78)
println("[4] partner co-cooling (aerosol riding the moist plume)")
println("="^78)
# box 1 = water (produces cooling), box 2 = PARTNER (sulfate: feels it).
# Parcel 3 (box 2) sits in the SAME met cell as the 200-ppmv water parcel 1 →
# cell-mean excess = 195 ppmv (parcel 2 is below background and contributes
# nothing) → the partner cools exactly like the water parcel.
mkprp(fb) = PlumeRadiation([true, false], [false, true], crad, q_bg, fb,
                           zeros(Float32, 100))
st4 = mkstate(); pr4 = mkprp(false)
n4, _ = apply_plume_radiation!(st4, fctx, pr4, dt, [1.0])
check("partner in the water cell cools by the cell-mean excess",
      isapprox(Float64(pr4.trad[3]), dT_exp; rtol=1e-5),
      info=@sprintf("%.4g vs %.4g K", pr4.trad[3], dT_exp))
check("water parcel unchanged by the partner pass",
      isapprox(Float64(pr4.trad[1]), dT_exp; rtol=1e-5))
check("touched count includes the partner", n4 == 2)
st5 = mkstate(); st5.indx[3] = 6.0f0            # partner moved to a DRY cell
pr5 = mkprp(false)
apply_plume_radiation!(st5, fctx, pr5, dt, [1.0])
check("partner in a dry cell untouched (separation switches its cooling off)",
      pr5.trad[3] == 0.0f0)
st6 = mkstate(); pr6 = mkprp(true)              # feedback: partner descends too
apply_plume_radiation!(st6, fctx, pr6, dt, [1.0])
check("partner descends with the plume while co-located",
      st6.indl[3] < 5.0f0 && isapprox(Float64(st6.indl[3]), Float64(st6.indl[1]); rtol=1e-4))
check("any_prad true for partner-only config",
      any_prad(PlumeRadiation([false], [true], crad, q_bg, false, zeros(Float32, 1))))

println()
println("="^78)
println("[5] prognostic cold anomaly T′ ([physics.plume_radiation].tprime_tau)")
println("="^78)
# full positional form: ..., h_ref, tprime_tau, tanom
τ_tp = 6000.0
mkprt(tau) = PlumeRadiation([true, false], falses(2), crad, 0.0, q_bg, true,
                            zeros(Float32, 100), Float32[], 1000.0,
                            tau, zeros(Float32, 100))
relfrac = -expm1(-dt / τ_tp)
st7 = mkstate(); pr7 = mkprt(τ_tp)
pr7.tanom[50] = -9.0f0                          # stale anomaly in a freed slot
apply_plume_radiation!(st7, fctx, pr7, dt, [1.0])
check("first step charges T′ by dT (reservoir, not descent)",
      isapprox(Float64(pr7.tanom[1]), dT_exp; rtol=1e-5),
      info=@sprintf("%.4g vs %.4g K", pr7.tanom[1], dT_exp))
check("no immediate descent (cooling goes to the reservoir)",
      st7.indl[1] == 5.0f0)
check("freed slot T′ reset", pr7.tanom[50] == 0.0f0)
# steady state: T′ → dT/relfrac ≈ −Q_rad·τ
NST = 200
st8 = mkstate(); pr8 = mkprt(τ_tp)
for _ in 1:NST; apply_plume_radiation!(st8, fctx, pr8, dt, [1.0]); end
tinf = dT_exp / relfrac
check("steady state T′ ≈ −Q_rad·τ = dT/(1−e^{−dt/τ})",
      isapprox(Float64(pr8.tanom[1]), tinf; rtol=0.05),
      info=@sprintf("%.4g vs %.4g K", pr8.tanom[1], tinf))
check("relaxation drives descent (indl fell)", st8.indl[1] < 5.0f0)
# energy accounting: total descent ≈ legacy minus the standing anomaly
st9 = mkstate(); pr9 = mkpr(true)               # legacy (τ = 0)
for _ in 1:NST; apply_plume_radiation!(st9, fctx, pr9, dt, [1.0]); end
d_tp  = Float64(5.0f0 - st8.indl[1])
d_leg = Float64(5.0f0 - st9.indl[1])
check("total descent = legacy − standing anomaly (no double counting)",
      isapprox(d_tp, d_leg * (1.0 - tinf / (NST * dT_exp)); rtol=0.05),
      info=@sprintf("T′ %.4g vs legacy %.4g Δindl", d_tp, d_leg))
# a parcel that leaves the plume keeps relaxing: T′ decays, descent continues
stA = mkstate(); prA = mkprt(τ_tp)
apply_plume_radiation!(stA, fctx, prA, dt, [1.0])
t_before = Float64(prA.tanom[1])
q_dry = 2.0 / (MW * 1e6)                        # below background: no new cooling
stA.airm[1] = Float32(1.0e6 / q_dry)
apply_plume_radiation!(stA, fctx, prA, dt, [1.0])
check("out-of-plume parcel: T′ decays by (1−relfrac), still descends",
      isapprox(Float64(prA.tanom[1]), t_before * (1.0 - relfrac); rtol=1e-4) &&
      stA.indl[1] < 5.0f0)
# partner parcels charge their own T′
prB = PlumeRadiation([true, false], [false, true], crad, 0.0, q_bg, false,
                     zeros(Float32, 100), Float32[], 1000.0,
                     τ_tp, zeros(Float32, 100))
stB = mkstate()
apply_plume_radiation!(stB, fctx, prB, dt, [1.0])
check("partner in the water cell charges its own T′",
      isapprox(Float64(prB.tanom[3]), dT_exp; rtol=1e-5))

println()
println("="^78)
println("[6] T′ reaches the saturation gates (the RHi wiring fix)")
println("="^78)
# T, p exactly as the gates compute them: SAME interpolation, SAME inversion
# (indl = 5.0 lands mid-cell → θ = 590, not the level-5 value 580)
θ5 = Float64(EASLEE.Operators.interpolation4Dsingle(m3.theta, 4.0f0, 4.0f0, 5.0f0, 1.0f0))
ρ5 = Float64(EASLEE.Operators.interpolation4Dsingle(m3.rho,   4.0f0, 4.0f0, 5.0f0, 1.0f0))
Tg = max(θ5^(1004.0 / 716.95) * (ρ5 * 287.05 / 1.0e5)^(287.05 / 716.95), 100.0)
pg = ρ5 * 287.05 * Tg
esat = T -> EASLEE.Operators.esat_ice(T)
qsi_g  = 0.622 * esat(Tg) / (pg - esat(Tg))
qsi_g5 = 0.622 * esat(Tg - 5.0) / (pg - esat(Tg - 5.0))
check("test premise: 5 K anomaly separates the two gates (qsi drops)",
      qsi_g5 < 0.9 * qsi_g,
      info=@sprintf("qsi(T)=%.3g qsi(T−5)=%.3g", qsi_g, qsi_g5))
q_test = 0.8 * qsi_g                            # subsaturated ambient, saturated at T−5
@assert q_test > qsi_g5
mkcv() = EASLEE.SpeciesConvert(
    Float32[1.0, 0.0], Int32[2, 0], Float32[1.0, 1.0], zeros(Float64, 2),
    Int8[2, 0], Float32[1.0, 1.5], Float32[0.0, 0.0], Float32[0.0, 0.0])  # SUBLIMATE, σ=0
function mkstc()
    st = EASLEE.State.ParcelState(100)
    st.status[1] = 0x01; st.nest[1] = Int32(1); st.relbox[1] = Int32(1)
    st.indx[1] = 4.0f0; st.indy[1] = 4.0f0; st.indl[1] = 5.0f0
    st.mass[1] = 1.0f6; st.airm[1] = Float32(1.0e6 / q_test)
    return st
end
# SUBLIMATE gate: subsaturated at ambient T → converts; saturated at T+T′ → held
stc = mkstc()
EASLEE.Operators.apply_convert_thermo!(stc, fctx, mkcv(), 600.0, [1.0])
check("no T′: subsaturated parcel sublimates (ICE → daughter box)",
      stc.relbox[1] == Int32(2))
stc2 = mkstc()
ta = zeros(Float32, 100); ta[1] = -5.0f0
EASLEE.Operators.apply_convert_thermo!(stc2, fctx, mkcv(), 600.0, [1.0]; tanom = ta)
check("T′ = −5 K: same parcel is saturated — ice HELD (the §5 stakes)",
      stc2.relbox[1] == Int32(1))
# FREEZE gate ([species.ice]): no shed at ambient, sheds with T′
mkice() = EASLEE.IcePhysics([true], Float32[1.0], Float32[5.0e-6],
                            Float32[500.0], Float32[0.0], zeros(Float64, 1))
sti = mkstc(); ice1 = mkice()
EASLEE.Operators.apply_ice!(sti, fctx, ice1, 600.0, [1.0])
check("apply_ice!, no T′: subsaturated → nothing freezes out",
      sum(ice1.frozen) == 0.0 && sti.mass[1] == 1.0f6)
sti2 = mkstc(); ice2 = mkice()
EASLEE.Operators.apply_ice!(sti2, fctx, ice2, 600.0, [1.0]; tanom = ta)
check("apply_ice!, T′ = −5 K: supersaturated → excess vapour sheds",
      sum(ice2.frozen) > 0.0 && sti2.mass[1] < 1.0f6)

println()
println("="^78)
println("[3] config")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-16T00:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "physics"=>Dict{String,Any}("plume_radiation"=>Dict{String,Any}(
      "species"=>"H2O_volc","partners"=>["sulfate"],
      "c_rad"=>0.003,"q_bg"=>4.0,"feedback"=>false,"tprime_tau"=>"1d")),
  "species"=>[Dict{String,Any}("id"=>"H2O_volc","kind"=>"gas")],
  "release"=>[Dict{String,Any}("species"=>"H2O_volc","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>10,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse (string species sugar + partners)",
      sp.physics.prad_species == ["H2O_volc"] && sp.physics.prad_c_rad == 0.003 &&
      sp.physics.prad_q_bg == 4.0 && !sp.physics.prad_feedback &&
      sp.physics.prad_partners == ["sulfate"])
check("parse tprime_tau duration", sp.physics.prad_tprime_tau == 86400.0)
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
check("round-trip", sp2.physics.prad_species == ["H2O_volc"] &&
      sp2.physics.prad_c_rad == 0.003 && !sp2.physics.prad_feedback &&
      sp2.physics.prad_partners == ["sulfate"] &&
      sp2.physics.prad_tprime_tau == 86400.0)

println("="^78)
println("[9] BACKWARD runs REVERSE the radiative tendency (fixed 2026-08-07)")
println("="^78)
# A layer that cools and SINKS forward in time was WARMER and HIGHER in the
# past.  The timeloop passes |dt|, so without the direction factor a backward
# run cooled — and therefore sank — too: the wrong sign for trajectory work.
stf = mkstate(); prf = mkpr(false)
apply_plume_radiation!(stf, fctx, prf, dt, [1.0])
stb = mkstate(); prb = mkpr(false)
apply_plume_radiation!(stb, bctx, prb, dt, [1.0])
check("forward cools (dT < 0)", prf.trad[1] < 0.0f0,
      info=@sprintf("%+.4g K", prf.trad[1]))
check("backward WARMS (dT > 0)", prb.trad[1] > 0.0f0,
      info=@sprintf("%+.4g K", prb.trad[1]))
check("backward is the exact mirror of forward",
      isapprox(Float64(prb.trad[1]), -Float64(prf.trad[1]); rtol=1e-6),
      info=@sprintf("%+.4g vs %+.4g K", prb.trad[1], -prf.trad[1]))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
