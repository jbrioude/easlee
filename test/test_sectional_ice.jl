#!/usr/bin/env julia
# test_sectional_ice.jl — SECTIONAL ice ([species.convert].sectional):
#   [1] the r² law: r'² = r² − 2K·dt integrated exactly, mass fraction kept
#       = (r'/r)³, and the vapour credited to the daughter budget,
#   [2] SIZE SELECTIVITY — the property neither the lamina nor T′ can provide:
#       for one shared (T,p,S) and one dt, a small bin loses a far larger mass
#       fraction than a large one (the saturating exponential loses 100 % of
#       BOTH — that is the ht_full3 null this replaces),
#   [3] bin-to-bin demotion: a parcel whose radius crosses the bin edge moves
#       to the next smaller ICE box; off the bottom of the ladder it becomes
#       vapour and its radius is cleared,
#   [4] size-resolved settling: a shrunken parcel falls slower than its box's
#       nominal r_eff would imply,
#   [5] legacy: empty rice = the exponential path, bit-identical,
#   [6] config: sectional parse + convert-block round-trip (mode / rhi_crit /
#       lamina_sigma / sectional were previously DROPPED by the writer).
#
#   julia --project=easlee test/test_sectional_ice.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: apply_convert_thermo!, apply_species_settling!,
                         esat_ice, interpolation4Dsingle, any_sectional

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-62s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic met (same inversion the gates use) ─────────────────────────────
# UNIFORM θ/ρ chosen so the gates' own inversion returns the HT plume state
# T ≈ 220 K at 20 hPa — the conditions every number in ICE_UPGRADE_DESIGN §0
# is quoted at (5 µm crystal, RHi 0.04 → 0.65 s).
const NZ = 10
dims = (8, 8, NZ, 2)
const T_TARGET = 220.0
const P_TARGET = 2000.0                          # 20 hPa
zlev = zeros(Float32, dims)
for k in 1:NZ; zlev[:, :, k, :] .= Float32((k - 1) * 1000.0); end
rho   = fill(Float32(P_TARGET / (287.05 * T_TARGET)), dims)
theta = fill(Float32(T_TARGET * (1.0e5 / P_TARGET)^0.2857), dims)
puv   = fill(Float32(P_TARGET), dims)
m3 = (theta = theta, rho = rho, zlev = zlev, puv = puv)
LC = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
fctx = (nests = [(metfields3D = m3, levelconfig = LC,)],)

θ5 = Float64(interpolation4Dsingle(m3.theta, 4.0f0, 4.0f0, 5.0f0, 1.0f0))
ρ5 = Float64(interpolation4Dsingle(m3.rho,   4.0f0, 4.0f0, 5.0f0, 1.0f0))
Tg = max(θ5^(1004.0/716.95) * (ρ5*287.05/1.0e5)^(287.05/716.95), 100.0)
pg = ρ5 * 287.05 * Tg
esi = esat_ice(Tg); qsi = 0.622 * esi / (pg - esi)
Dv = 2.11e-5 * (1.01325e5 / pg) * (Tg / 273.15)^1.94
rs = esi / (461.5 * Tg)
Kof(S) = Dv * (1.0 - min(S, 0.999999)) * rs / 917.0       # m²/s
# TWO regimes, because they are physically different:
#  · near-saturation (inside a cooled/lamina parcel): t_sub ≫ dt, crystals
#    shrink gently and the r² law RESOLVES the sizes,
#  · strongly subsaturated ambient: t_sub = 0.03–0.7 s, where the legacy
#    exponential saturates to 1 for EVERY bin and size cannot matter.
const S_NEAR = 0.999
const S_DRY  = 0.04
K = Kof(S_NEAR)
q_test = S_NEAR * qsi
q_dry  = S_DRY  * qsi
@printf("met: T=%.2f K  p=%.1f Pa  qsi=%.4g\n", Tg, pg, qsi)
@printf("  S=%.3f: K=%.3e m2/s  t_sub(1um)=%.1f s  t_sub(5um)=%.1f s\n",
        S_NEAR, K, 1e-12/(2K), 25e-12/(2K))
@printf("  S=%.2f : t_sub(1um)=%.3f s  t_sub(5um)=%.3f s  (exponential saturates)\n",
        S_DRY, 1e-12/(2Kof(S_DRY)), 25e-12/(2Kof(S_DRY)))

# two ICE bins: box1 = small (1 µm), box2 = large (5 µm); box3 = vapour target.
# ladder: box1 is the BOTTOM (bin_down 0), box2 demotes to box1.
R_SMALL, R_LARGE = 1.0e-6, 5.0e-6
REDGE2 = sqrt(R_SMALL * R_LARGE)                 # geometric edge ≈ 2.24 µm
NB = 3
function mkcv(; sectional=true, nslots=20)
    lambda = Float32[1.0, 1.0, 0.0]
    target = Int32[3, 3, 0]                      # both ICE boxes → vapour box 3
    ratio  = Float32[1.0, 1.0, 1.0]
    mode   = Int8[2, 2, 0]                       # SUBLIMATE
    rhi    = Float32[1.5, 1.5, 1.5]
    lsg    = zeros(Float32, NB)                  # sigma = 0 (standing decision)
    sbr    = Float32[R_SMALL, R_LARGE, 0.0]      # seeds / legacy radii
    rice   = sectional ? zeros(Float32, nslots) : Float32[]
    bdn    = Int32[0, 1, 0]                      # box2 → box1; box1 = bottom
    bre    = sectional ? Float32[R_SMALL*sqrt(R_SMALL/R_LARGE), REDGE2, 0.0] :
                         zeros(Float32, NB)
    EASLEE.SpeciesConvert(lambda, target, ratio, zeros(Float64, NB),
                          mode, rhi, lsg, sbr, rice, bdn, bre)
end
function mkst(box::Int; mass=1.0f6, nslots=20, q=q_test)
    st = EASLEE.State.ParcelState(nslots)
    st.status[1] = 0x01; st.nest[1] = Int32(1); st.relbox[1] = Int32(box)
    st.indx[1] = 4.0f0; st.indy[1] = 4.0f0; st.indl[1] = 5.0f0
    st.mass[1] = mass; st.airm[1] = Float32(Float64(mass) / q)
    return st
end

println("="^78)
println("[1] exact r² law: r'² = r² − 2K·dt, mass kept = (r'/r)³")
println("="^78)
dt = 10.0
cv1 = mkcv(); st1 = mkst(2)                      # large bin, stays in its bin
tot1 = apply_convert_thermo!(st1, fctx, cv1, dt, [1.0])
r_exp = sqrt(R_LARGE^2 - 2K*dt)
check("radius follows r' = sqrt(r² − 2K·dt)",
      isapprox(Float64(cv1.rice[1]), r_exp; rtol=1e-5),
      info=@sprintf("%.4f vs %.4f um", cv1.rice[1]*1e6, r_exp*1e6))
f_keep = (r_exp / R_LARGE)^3
check("mass kept = (r'/r)³",
      isapprox(Float64(st1.mass[1]), 1.0e6 * f_keep; rtol=1e-5),
      info=@sprintf("%.5g vs %.5g", st1.mass[1], 1.0e6*f_keep))
# intrinsic invariant (independent of the Float32 radius seeding): every kg
# leaving the parcel is credited to the daughter budget, exactly.
# state.mass is Float32 (ulp ≈ 0.06 kg at 1e6), so the closure can only be
# checked to single precision — the budget itself is accumulated in Float64.
check("mass conserved: m_before = m_after + returned budget",
      isapprox(1.0e6, Float64(st1.mass[1]) + tot1; rtol=1e-6) &&
      isapprox(cv1.converted[2], tot1; rtol=1e-12) &&
      isapprox(tot1, 1.0e6 * (1 - f_keep); rtol=1e-4),
      info=@sprintf("lost %.5g kg, residual %.3g", tot1,
                    1.0e6 - Float64(st1.mass[1]) - tot1))
check("parcel KEEPS its ice identity (no all-or-nothing switch)",
      st1.relbox[1] == Int32(2))

println()
println("="^78)
println("[2] SIZE SELECTIVITY — small bin loses far more than large")
println("="^78)
cvS = mkcv(); stS = mkst(1)                      # 1 µm
cvL = mkcv(); stL = mkst(2)                      # 5 µm
apply_convert_thermo!(stS, fctx, cvS, dt, [1.0])
apply_convert_thermo!(stL, fctx, cvL, dt, [1.0])
lossS = 1.0 - Float64(stS.mass[1]) / 1.0e6
lossL = 1.0 - Float64(stL.mass[1]) / 1.0e6
check("small bin loses a strictly larger mass fraction",
      lossS > lossL,
      info=@sprintf("1um %.1f%% vs 5um %.1f%%", 100lossS, 100lossL))
check("the ratio is large (r² law, not a saturated exponential)",
      lossS / max(lossL, 1e-12) > 3.0,
      info=@sprintf("x%.1f", lossS / max(lossL, 1e-12)))
# In SUBSATURATED ambient air (S = 0.04, dt = 300 s = the production step) the
# legacy exponential has p = 1 − exp(−dt/t_sub) → 1 for BOTH bins: it
# annihilates the 5 µm crystal exactly as fast as the 1 µm one.  That is the
# structural null ht_full3 reproduced over 8 days (r_eff flat).
dt_prod = 300.0
cvSx = mkcv(sectional=false); stSx = mkst(1; q=q_dry)
cvLx = mkcv(sectional=false); stLx = mkst(2; q=q_dry)
apply_convert_thermo!(stSx, fctx, cvSx, dt_prod, [1.0])
apply_convert_thermo!(stLx, fctx, cvLx, dt_prod, [1.0])
check("legacy exponential converts BOTH wholesale (the null it replaces)",
      stSx.relbox[1] == Int32(3) && stLx.relbox[1] == Int32(3))
# the sectional law in the SAME dry air also removes both (correctly — a
# crystal cannot survive RHi 0.04), so the gain is not "ice survives dry air":
# it is that near saturation the sizes SEPARATE instead of all dying at once
cvSd = mkcv(); stSd = mkst(1; q=q_dry)
cvLd = mkcv(); stLd = mkst(2; q=q_dry)
apply_convert_thermo!(stSd, fctx, cvSd, dt_prod, [1.0])
apply_convert_thermo!(stLd, fctx, cvLd, dt_prod, [1.0])
check("sectional in dry air also removes both (no false ice survival)",
      stSd.relbox[1] == Int32(3) && stLd.relbox[1] == Int32(3))

println()
println("="^78)
println("[3] bin-to-bin demotion and the bottom of the ladder")
println("="^78)
# shrink the large bin past the geometric edge in one step
dt_big = (R_LARGE^2 - (0.9*REDGE2)^2) / (2K)
cv3 = mkcv(); st3 = mkst(2)
apply_convert_thermo!(st3, fctx, cv3, dt_big, [1.0])
check("radius crossing the bin edge demotes the parcel to the smaller box",
      st3.relbox[1] == Int32(1) && Float64(cv3.rice[1]) < REDGE2,
      info=@sprintf("box %d, r=%.3f um (edge %.3f)", st3.relbox[1],
                    cv3.rice[1]*1e6, REDGE2*1e6))
check("demotion does NOT change the parcel's mass (loss already applied)",
      st3.mass[1] < 1.0f6 && st3.mass[1] > 0.0f0)
# off the bottom: shrink box1 below its lower edge
dt_gone = (R_SMALL^2 - (0.5*R_SMALL*sqrt(R_SMALL/R_LARGE))^2) / (2K)
cv4 = mkcv(); st4 = mkst(1)
apply_convert_thermo!(st4, fctx, cv4, dt_gone, [1.0])
check("off the bottom of the ladder → vapour box, radius cleared",
      st4.relbox[1] == Int32(3) && cv4.rice[1] == 0.0f0)
# annihilation within a step (r² ≤ 0) → all mass to vapour
cv5 = mkcv(); st5 = mkst(1)
apply_convert_thermo!(st5, fctx, cv5, R_SMALL^2 / (2K) * 2.0, [1.0])
check("r² ≤ 0 within the step → whole parcel to vapour",
      st5.relbox[1] == Int32(3) && cv5.rice[1] == 0.0f0 &&
      isapprox(cv5.converted[1], 1.0e6; rtol=1e-9))

println()
println("="^78)
println("[4] settling reads the parcel's own radius (size-resolved)")
println("="^78)
spp = EASLEE.SpeciesPhysics(zeros(Float32,NB), fill(-1f0,NB), zeros(Float32,NB),
        falses(NB), fill(-1f0,NB), zeros(Float32,NB), zeros(Float32,NB),
        zeros(Float32,NB), 30f0, [zeros(Float64,8,8)],
        Float32[0.0, R_LARGE, 0.0], fill(917f0, NB))
gc_ = EASLEE.Types.Rectilinear(Int32(8), Int32(8), Int32(NZ), 1f0, 1f0,
        0f0, -50f0, Float32(1/111320.0), Float32(1/111320.0),
        1f0, 0f0, 1f0, 0f0, false, false, false, 1f0, 1f0, false)
sctx = (nests = [(metfields3D = m3, levelconfig = LC, gridconfig = gc_)],)
stA = mkst(2); stA.indl[1] = 5.0f0
apply_species_settling!(stA, sctx, spp, 3600.0, [1.0])            # nominal 5 µm
stB = mkst(2); stB.indl[1] = 5.0f0
rsmall = zeros(Float32, 20); rsmall[1] = Float32(0.2 * R_LARGE)   # shrunken
apply_species_settling!(stB, sctx, spp, 3600.0, [1.0]; rice = rsmall)
check("shrunken crystal falls LESS far than the box's nominal r_eff",
      stB.indl[1] > stA.indl[1] && stA.indl[1] < 5.0f0,
      info=@sprintf("indl %.4f (shrunk) vs %.4f (nominal)", stB.indl[1], stA.indl[1]))

println()
println("="^78)
println("[5] legacy path untouched when sectional is off")
println("="^78)
check("any_sectional false with empty rice", !any_sectional(mkcv(sectional=false)))
check("any_sectional true when armed", any_sectional(mkcv()))
cv6 = mkcv(sectional=false); st6 = mkst(2; q=q_dry)
t6 = apply_convert_thermo!(st6, fctx, cv6, dt_prod, [1.0])
check("legacy: whole-parcel switch, mass rescaled by ratio",
      st6.relbox[1] == Int32(3) && isapprox(Float64(st6.mass[1]), 1.0e6; rtol=1e-9) &&
      isapprox(t6, 1.0e6; rtol=1e-9))

println()
println("="^78)
println("[6] config: sectional parse + convert-block round-trip")
println("="^78)
d = Dict{String,Any}("name"=>"x", "task"=>"transport",
  "run"=>Dict{String,Any}("start"=>"2022-01-15T00:00:00","stop"=>"2022-01-16T00:00:00"),
  "meteorology"=>[Dict{String,Any}("model"=>"ERA5","path"=>"/tmp","variables"=>"E")],
  "species"=>[Dict{String,Any}("id"=>"ICE_b1","kind"=>"gas",
      "convert"=>Dict{String,Any}("to"=>"H2O_b1","tau"=>"1h","mode"=>"sublimate",
                                  "lamina_sigma"=>0.0,"sectional"=>true),
      "fall"=>Dict{String,Any}("r_eff"=>"0.87um")),
    Dict{String,Any}("id"=>"H2O_b1","kind"=>"gas",
      "convert"=>Dict{String,Any}("to"=>"ICE_b1","tau"=>"1h","mode"=>"freeze",
                                  "rhi_crit"=>1.2))],
  "release"=>[Dict{String,Any}("species"=>"H2O_b1","lon"=>[0.0,1.0],"lat"=>[0.0,1.0],
      "z"=>[0.0,100.0],"time"=>["2022-01-15T00:00:00","2022-01-15T01:00:00"],
      "particles"=>10,"mass"=>"1kg")])
sp = EasleeConfig.ConfigTOML.dict_to_runspec(d)
ice = sp.species[findfirst(s -> s.id == "ICE_b1", sp.species)]
h2o = sp.species[findfirst(s -> s.id == "H2O_b1", sp.species)]
check("parse sectional + mode + r_eff",
      ice.convert_sectional && ice.convert_mode == "sublimate" &&
      isapprox(ice.settle_r_eff, 0.87e-6; rtol=1e-6))
rt = EasleeConfig.ConfigTOML.runspec_to_dict(sp)
sp2 = EasleeConfig.ConfigTOML.dict_to_runspec(rt)
ice2 = sp2.species[findfirst(s -> s.id == "ICE_b1", sp2.species)]
h2o2 = sp2.species[findfirst(s -> s.id == "H2O_b1", sp2.species)]
check("round-trip keeps sectional AND the thermodynamic gate (was dropped)",
      ice2.convert_sectional && ice2.convert_mode == "sublimate" &&
      h2o2.convert_mode == "freeze" && isapprox(h2o2.convert_rhi_crit, 1.2))
check("round-trip default: sectional off stays off", !h2o2.convert_sectional)

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
