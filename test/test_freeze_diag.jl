#!/usr/bin/env julia
# test_freeze_diag.jl — convert-gate diagnostics (2026-08-05 pre-registration):
#   [1] flux accumulators: `demoted` on a bin demotion (pure relabel, mass
#       unchanged), `gone` on crystal annihilation and off-bottom removal —
#       and `gone` ⊆ `converted` (subset, not a second budget),
#   [2] freeze-event moments: every mode-1 (vapour→ice) conversion records
#       mass-weighted ambient T, applied T′, ln q, ln q_si, ln p,
#   [3] the STAGGERED−PLAIN gate-read offset (fz_mdT): zero on a uniform
#       θ/ρ field (null), and equal to the definitional half-level difference
#       on a vertically graded field — the §0c.2-item-7 convention split,
#       measured where it bites,
#   [4] per-parcel gate state (gate_T/gate_qsi/gate_p) filled when armed via
#       diag_slots, absent (and cost-free) when not,
#   [5] backward compatibility: the 11-arg constructor without diag_slots
#       leaves the gate vectors empty and the physics bit-identical,
#   [6] config: [run].parcel_dump parse + write round-trip.
#
#   julia --project=easlee test/test_freeze_diag.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: apply_convert_thermo!, esat_ice, interpolation4Dsingle

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic met: same harness as test_sectional_ice.jl ─────────────────────
const NZ = 10
dims = (8, 8, NZ, 2)
const T_TARGET = 220.0
const P_TARGET = 2000.0
zlev = zeros(Float32, dims)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * 1000.0); end
rho   = fill(Float32(P_TARGET / (287.05 * T_TARGET)), dims)
theta = fill(Float32(T_TARGET * (1.0e5 / P_TARGET)^0.2857), dims)
m3 = (theta = theta, rho = rho, zlev = zlev, puv = fill(Float32(P_TARGET), dims))
LC = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
fctx = (nests = [(metfields3D = m3, levelconfig = LC,)],)

θ5 = Float64(interpolation4Dsingle(m3.theta, 4.0f0, 4.0f0, 5.0f0, 1.0f0))
ρ5 = Float64(interpolation4Dsingle(m3.rho,   4.0f0, 4.0f0, 5.0f0, 1.0f0))
Tg = max(θ5^(1004.0/716.95) * (ρ5*287.05/1.0e5)^(287.05/716.95), 100.0)
pg = ρ5 * 287.05 * Tg
esi = esat_ice(Tg); qsi = 0.622 * esi / (pg - esi)
Dv = 2.11e-5 * (1.01325e5 / pg) * (Tg / 273.15)^1.94
rs = esi / (461.5 * Tg)
K = Dv * (1.0 - 0.999) * rs / 917.0

# boxes: 1 = vapour (FREEZE → 2), 2 = large ICE (5 µm, demotes → 3),
#        3 = small ICE (1 µm, bottom), 4 = vapour target of sublimation
R_SMALL, R_LARGE = 1.0e-6, 5.0e-6
REDGE2 = sqrt(R_SMALL * R_LARGE)
NB = 4
function mkcv(; diag_slots=0, nslots=20)
    lambda = Float32[1.0, 1.0, 1.0, 0.0]
    target = Int32[2, 4, 4, 0]
    ratio  = Float32[1.0, 1.0, 1.0, 1.0]
    mode   = Int8[1, 2, 2, 0]                    # 1 freeze; 2 sublimate
    rhi    = Float32[1.5, 1.5, 1.5, 1.5]
    lsg    = zeros(Float32, NB)
    sbr    = Float32[0.0, R_LARGE, R_SMALL, 0.0]
    rice   = zeros(Float32, nslots)
    bdn    = Int32[0, 3, 0, 0]
    bre    = Float32[0.0, REDGE2, R_SMALL*sqrt(R_SMALL/R_LARGE), 0.0]
    EASLEE.SpeciesConvert(lambda, target, ratio, zeros(Float64, NB),
                          mode, rhi, lsg, sbr, rice, bdn, bre;
                          diag_slots = diag_slots)
end
function mkst(box::Int; mass=1.0f6, nslots=20, q=0.999*qsi)
    st = EASLEE.State.ParcelState(nslots)
    st.status[1] = 0x01; st.nest[1] = Int32(1); st.relbox[1] = Int32(box)
    st.indx[1] = 4.0f0; st.indy[1] = 4.0f0; st.indl[1] = 5.0f0
    st.mass[1] = mass; st.airm[1] = Float32(Float64(mass) / q)
    return st
end

println("="^78)
println("[1] flux accumulators: demoted / gone")
println("="^78)
dt_big = (R_LARGE^2 - (0.9*REDGE2)^2) / (2K)
cv = mkcv(); st = mkst(2)
apply_convert_thermo!(st, fctx, cv, dt_big, [1.0])
check("demotion increments `demoted[src]` by the parcel's (post-loss) mass",
      st.relbox[1] == Int32(3) &&
      isapprox(cv.demoted[2], Float64(st.mass[1]); rtol=1e-6),
      info=@sprintf("demoted=%.4g, mass=%.4g", cv.demoted[2], st.mass[1]))
check("demotion adds NOTHING to `converted`/`gone` (pure relabel)",
      cv.gone[2] == 0.0 && cv.converted[2] < 1.0e6 * 0.999)
cv2 = mkcv(); st2 = mkst(3)
apply_convert_thermo!(st2, fctx, cv2, R_SMALL^2 / (2K) * 2.0, [1.0])
check("annihilation within a step: gone == converted for that box",
      st2.relbox[1] == Int32(4) &&
      isapprox(cv2.gone[3], cv2.converted[3]; rtol=1e-12) &&
      isapprox(cv2.gone[3], 1.0e6; rtol=1e-9))
# one-bin-per-step: first call demotes 2→3; a second call shrinks past b3's
# bottom edge → off the ladder → vapour + gone.  K depends on the parcel's OWN
# S = q/q_si, and the first call reduced its mass — restore S = 0.999 by
# resetting airm, so the second call's dt lands between the edges (a partial
# loss then removal) instead of annihilating within the step.
cv3 = mkcv(); st3 = mkst(2)
apply_convert_thermo!(st3, fctx, cv3, dt_big, [1.0])
m_after_demote = Float64(st3.mass[1])
st3.airm[1] = Float32(m_after_demote / (0.999 * qsi))
apply_convert_thermo!(st3, fctx, cv3, (Float64(cv3.rice[1])^2 -
                      (0.5*R_SMALL*sqrt(R_SMALL/R_LARGE))^2) / (2K), [1.0])
check("off-bottom removal: gone[bottom] gets the residual (post-loss) mass",
      st3.relbox[1] == Int32(4) && cv3.gone[3] > 0.0 &&
      cv3.gone[3] < 0.05 * m_after_demote,
      info=@sprintf("gone=%.4g of %.4g", cv3.gone[3], m_after_demote))

println()
println("="^78)
println("[2] freeze-event moments (mode 1)")
println("="^78)
q_super = 2.0 * qsi                             # RHi = 2.0 > rhi_crit
cvF = mkcv(); stF = mkst(1; q=q_super)
totF = apply_convert_thermo!(stF, fctx, cvF, 300.0, [1.0])
mF = cvF.converted[1]
check("freeze fired and moments carry its mass",
      stF.relbox[1] == Int32(2) && mF > 0.0 && cvF.fz_mT[1] > 0.0)
check("ambient T moment: fz_mT/m == T_amb",
      isapprox(cvF.fz_mT[1] / mF, Tg; rtol=1e-6),
      info=@sprintf("%.2f K vs %.2f K", cvF.fz_mT[1]/mF, Tg))
check("no T′ passed → fz_mTp == 0", cvF.fz_mTp[1] == 0.0)
check("ln q moment matches the parcel's q",
      isapprox(exp(cvF.fz_mlq[1] / mF), q_super; rtol=1e-5))
check("ln q_si moment matches the gate's q_si",
      isapprox(exp(cvF.fz_mlqsi[1] / mF), qsi; rtol=1e-5))
check("ln p moment matches the gate's p",
      isapprox(exp(cvF.fz_mlp[1] / mF), pg; rtol=1e-5))
# with a cold anomaly: gate T = T_amb + T′, moments split the two
tanom = zeros(Float32, 20); tanom[1] = -4.0f0
cvT = mkcv(); stT = mkst(1; q=q_super)
apply_convert_thermo!(stT, fctx, cvT, 300.0, [1.0]; tanom = tanom)
mT = cvT.converted[1]
check("T′ moment records the anomaly actually applied",
      mT > 0.0 && isapprox(cvT.fz_mTp[1] / mT, -4.0; rtol=1e-5) &&
      isapprox(cvT.fz_mT[1] / mT, Tg; rtol=1e-6),
      info=@sprintf("T'=%.2f K, T_amb=%.2f K", cvT.fz_mTp[1]/mT, cvT.fz_mT[1]/mT))
check("gate q_si is at T_amb + T′ (colder → smaller)",
      exp(cvT.fz_mlqsi[1] / mT) < qsi)

println()
println("="^78)
println("[3] staggered−plain gate offset (fz_mdT)")
println("="^78)
check("UNIFORM field: offset is zero (null test)",
      abs(cvF.fz_mdT[1] / mF) < 1e-6,
      info=@sprintf("%.2e K", cvF.fz_mdT[1]/mF))
# vertically graded θ: staggered read (mass level k+0.5) is warmer than the
# plain read at the same indl — the definitional half-level difference.
# Keep θ near the harness anchor (θ₅ = 672.7 K → T ≈ 220 K at this ρ): a
# large offset pushes T above the `esi < 0.9p` guard and the freeze never
# fires (that mistake produced a silent mG = 0 → NaN in this very test).
thetag = copy(theta)
for k in 1:NZ; thetag[:, :, k, :] .= Float32(645.2 + 5.5 * k); end
m3g = (theta = thetag, rho = rho, zlev = zlev, puv = m3.puv)
gctx = (nests = [(metfields3D = m3g, levelconfig = LC,)],)
θ_stag = Float64(interpolation4Dsingle(m3g.theta, 4.0f0, 4.0f0, 5.0f0, 1.0f0))
θ_plain = Float64(thetag[4, 4, 5, 1])
Tof(θ) = max(θ^(1004.0/716.95) * (ρ5*287.05/1.0e5)^(287.05/716.95), 100.0)
dT_expected = Tof(θ_stag) - Tof(θ_plain)
qsi_g = 0.622 * esat_ice(Tof(θ_stag)) / (ρ5*287.05*Tof(θ_stag) - esat_ice(Tof(θ_stag)))
cvG = mkcv(); stG = mkst(1; q=2.0*qsi_g)
apply_convert_thermo!(stG, gctx, cvG, 300.0, [1.0])
mG = cvG.converted[1]
check("GRADED field: fz_mdT equals the definitional stagger offset",
      mG > 0.0 && isapprox(cvG.fz_mdT[1] / mG, dT_expected; rtol=1e-4),
      info=@sprintf("%.3f K vs %.3f K expected", cvG.fz_mdT[1]/mG, dT_expected))
check("offset has the warm sign when θ increases upward",
      cvG.fz_mdT[1] / mG > 0.0)

println()
println("="^78)
println("[4]/[5] per-parcel gate state + constructor compatibility")
println("="^78)
cvD = mkcv(diag_slots=20); stD = mkst(1; q=q_super)
apply_convert_thermo!(stD, fctx, cvD, 300.0, [1.0])
check("diag_slots > 0: gate_T/gate_qsi/gate_p filled at the gate's values",
      isapprox(Float64(cvD.gate_T[1]), Tg; rtol=1e-5) &&
      isapprox(Float64(cvD.gate_qsi[1]), qsi; rtol=1e-4) &&
      isapprox(Float64(cvD.gate_p[1]), pg; rtol=1e-5))
check("diag_slots = 0 (default): gate vectors empty",
      isempty(cvF.gate_T) && isempty(cvF.gate_qsi) && isempty(cvF.gate_p))
cvE = mkcv(); stE = mkst(1; q=q_super)
totE = apply_convert_thermo!(stE, fctx, cvE, 300.0, [1.0])
check("physics identical with and without per-parcel diagnostics",
      stE.relbox[1] == stD.relbox[1] && isapprox(totE, mF; rtol=1e-12))

println()
println("="^78)
println("[6] config: [run].parcel_dump round-trip")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00",
                          "stop"=>"2022-01-16T00:00:00","parcel_dump"=>"3h"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[Dict{String,Any}("id"=>"X","kind"=>"gas")],
  "release"=>[Dict{String,Any}("species"=>"X","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>10,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
check("parse: parcel_dump read", sp.parcel_dump == "3h")
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
check("round-trip keeps parcel_dump", sp2.parcel_dump == "3h")
check("default stays off (not written when empty)",
      EasleeConfig.ConfigTOML.dict_to_runspec(
          EasleeConfig.ConfigTOML.runspec_to_dict(
              EasleeConfig.ConfigTOML.dict_to_runspec(
                  delete!(deepcopy(d), "run") === nothing ? d : begin
                      dd = deepcopy(d); delete!(dd["run"], "parcel_dump"); dd
                  end))).parcel_dump == "")

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
