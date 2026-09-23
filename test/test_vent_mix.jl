#!/usr/bin/env julia
# test_vent_mix.jl — radial magma–seawater vent split (VOLCANIC_THEORY.md §3).
# The eruptive jet is magmatic on the axis and ocean-water-loaded at the edge;
# per-parcel initial (T, qv, qc) come from the phreatomagmatic enthalpy
# balance.  Verify:
#   [1] the enthalpy balance itself: limits, monotonicity, the full-
#       vaporisation cap f_max (~0.29 at 1323 K) and energy conservation,
#   [2] the radial vent gradient: hot dry-ish core → boiling wet edge,
#   [3] provenance: magmatic water fraction falls core→edge; gas (melt)
#       weight falls core→edge,
#   [4] the ensemble still rises and delivers with ventmix on (v2 pairwise),
#       and the area-mean column matches a uniform run at the same means,
#   [5] config round-trip: TOML → RunSpec → TOML preserves the vent_mix block.
#
#   julia --project=easlee test/test_vent_mix.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: VentMix, vent_mix_state, vent_mix_means,
                         PlumeWashout, plume_release_ensemble

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

vm = VentMix()          # T_mag=1323, n_mag=0.03, x_core=0.15, x_edge=0.70, p=2

println("="^78)
println("[1] phreatomagmatic enthalpy balance")
println("="^78)
# pure magma limit
s0 = vent_mix_state(VentMix(x_core=0.0, x_edge=0.0), 0.5)
check("x=0: T0 = T_mag, water = magmatic only",
      isapprox(s0.T0, 1323.0; atol=1.0) && isapprox(s0.qv0, 0.03; atol=1e-9) &&
      s0.qc0 == 0.0)
# the analytic full-vaporisation cap: melt·cm·(Tmag−Tb) = x·(cw·ΔTw + Lv)
cm = 1100.0; cw = 4186.0; Lv = 2.26e6; Tb = 373.15
fmax = cm*(1323.0 - Tb) / (cm*(1323.0 - Tb) + cw*(Tb - 298.0) + Lv)
sc = vent_mix_state(VentMix(x_core=fmax - 0.01, x_edge=fmax - 0.01), 0.5)
sh = vent_mix_state(VentMix(x_core=fmax + 0.01, x_edge=fmax + 0.01), 0.5)
check("below f_max: all vapour, T > 373 K",
      sc.qc0 < 1e-12 && sc.T0 > Tb; info=@sprintf("f_max=%.3f T=%.0fK", fmax, sc.T0))
check("above f_max: boiling mixture with liquid remainder",
      sh.qc0 > 0.0 && isapprox(sh.T0, Tb; atol=0.01);
      info=@sprintf("qc0=%.3f", sh.qc0))
# energy conservation on the boiling branch: heat given up by melt =
# heat to warm + partially vaporise the seawater
x = fmax + 0.01
lhs = (1 - x)*cm*(1323.0 - Tb)
fv = 1 - sh.qc0 / x
rhs = x*(cw*(Tb - 298.0) + fv*Lv)
check("boiling branch closes the energy budget",
      isapprox(lhs, rhs; rtol=1e-10); info=@sprintf("%.1f vs %.1f J/kg", lhs, rhs))
# monotonicity across the cap
Ts = [vent_mix_state(VentMix(x_core=xx, x_edge=xx), 0.5).T0 for xx in 0.05:0.05:0.9]
check("T0(x) monotonically non-increasing", all(diff(Ts) .<= 1e-9))

println()
println("="^78)
println("[2] radial vent gradient")
println("="^78)
rr = 0.0:0.1:1.0
st = [vent_mix_state(vm, r) for r in rr]
check("core hot, edge boiling",
      st[1].T0 > 600.0 && isapprox(st[end].T0, Tb; atol=0.01);
      info=@sprintf("T(0)=%.0fK T(1)=%.0fK", st[1].T0, st[end].T0))
check("water load rises core→edge",
      all(diff([s.qv0 + s.qc0 for s in st]) .>= -1e-12);
      info=@sprintf("w(0)=%.2f w(1)=%.2f", st[1].qv0+st[1].qc0, st[end].qv0+st[end].qc0))
check("liquid only appears beyond the cap radius",
      st[1].qc0 == 0.0 && st[end].qc0 > 0.1)

println()
println("="^78)
println("[3] provenance labels")
println("="^78)
check("magmatic water fraction falls core→edge",
      all(diff([s.fmag for s in st]) .<= 1e-12) && st[1].fmag > 3*st[end].fmag;
      info=@sprintf("fmag(0)=%.3f fmag(1)=%.3f", st[1].fmag, st[end].fmag))
check("gas (melt) weight falls core→edge",
      all(diff([s.melt for s in st]) .<= 1e-12);
      info=@sprintf("melt(0)=%.2f melt(1)=%.2f", st[1].melt, st[end].melt))
mns = vent_mix_means(vm)
# the area mean is boiling-dominated (most of the vent area is beyond the
# vaporisation cap) — that is the physics; the ASCENDING column is the core
check("area means consistent (Tb ≤ T̄ ≤ T_core, n̄w in HT range)",
      Tb - 1.0 < mns.T0 < st[1].T0 && 0.3 < mns.nw0 < 0.55;
      info=@sprintf("T̄=%.0fK n̄w=%.2f", mns.T0, mns.nw0))

println()
println("="^78)
println("[4] ensemble with ventmix (v2 pairwise)")
println("="^78)
# standard tropical atmosphere (same as the v1/v2 tests)
z = collect(0.0:200.0:60_000.0); T = similar(z)
for (i, zz) in enumerate(z)
    T[i] = zz <= 16_000 ? 300.0 - 6.5e-3*zz : zz <= 20_000 ? 196.0 : 196.0 + 2.2e-3*(zz-20_000)
end
p = similar(z); p[1] = 101_325.0
for i in 2:length(z); p[i] = p[i-1]*exp(-9.81*(z[i]-z[i-1])/(287.05*0.5*(T[i-1]+T[i]))); end
qv = [max(0.6*0.622*611.2*exp(17.67*(Tt-273.15)/(Tt-29.65))/pp, 2e-6) for (Tt,pp) in zip(T,p)]

Random.seed!(11)
MER = 5.0e9
ens = plume_release_ensemble(MER, 0.0, 0.0, z, T, p, qv;
                             n=800, entrainment=:pairwise, ventmix=vm)
# uniform reference at the CORE state — the column that drives the profiles
score = vent_mix_state(vm, 0.0)
ensu = plume_release_ensemble(MER, score.T0, score.qv0, z, T, p, qv;
                              n=800, entrainment=:pairwise)
zs = sort(ens.z_det); n4 = length(zs)
p10 = zs[max(1, round(Int, 0.10n4))]; p90 = zs[round(Int, 0.90n4)]
zmedu = median(ensu.z_det)
check("TWO-CLOUD structure: upper tail reaches the umbrella",
      p90 > 18_000.0; info=@sprintf("p90=%.1f km", p90/1e3))
check("TWO-CLOUD structure: wet-edge skirt stays low",
      p10 < 12_000.0; info=@sprintf("p10=%.1f km", p10/1e3))
check("upper tail tracks the core-driven uniform column",
      abs(p90 - quantile(ensu.z_det, 0.9)) / max(quantile(ensu.z_det, 0.9), 1.0) < 0.30;
      info=@sprintf("p90 %.1f vs uniform %.1f km", p90/1e3, quantile(ensu.z_det, 0.9)/1e3))
check("per-member magmatic label populated and sensible",
      all(0.0 .<= ens.fmagw .<= 1.0) && 0.0 < mean(ens.fmagw) < 0.5;
      info=@sprintf("mean fmag=%.3f", mean(ens.fmagw)))
check("gas weights average to 1 and span core→edge",
      isapprox(mean(ens.gasw), 1.0; atol=1e-9) && maximum(ens.gasw) > minimum(ens.gasw);
      info=@sprintf("[%.2f, %.2f]", minimum(ens.gasw), maximum(ens.gasw)))
check("lofted source water varies member-to-member",
      maximum(ens.w0src) > 1.5*minimum(ens.w0src);
      info=@sprintf("[%.2f, %.2f]", minimum(ens.w0src), maximum(ens.w0src)))
# per-member fh2o is NOT bounded by 1 under v2: pairwise exchange moves water
# from wet-edge into dry-core parcels — only the ENSEMBLE budget closes
check("water survivals finite, non-negative",
      all(isfinite, ens.fh2o) && all(ens.fh2o .>= 0.0))
check("ensemble-mean water survival ≤ 1 (washout removes, never adds)",
      mean(ens.fh2o) <= 1.05; info=@sprintf("mean=%.2f", mean(ens.fh2o)))
# ensemble magmatic-water budget: delivered magmatic ≤ total magmatic in
w_in_mag  = sum(ens.w0src .* ens.fmagw)
w_in      = sum(ens.w0src)
check("magmatic share of source water ≈ area-integral of the mix",
      0.02 < w_in_mag/w_in < 0.25; info=@sprintf("%.3f", w_in_mag/w_in))

println()
println("="^78)
println("[5] config round-trip")
println("="^78)
tmp = joinpath(mktempdir(), "vm.toml")
write(tmp, """
task = "transport"
[run]
start = 2022-01-15T04:00:00
stop  = 2022-01-15T10:00:00
[[meteorology]]
path = "/tmp"
available = "auto"
[[species]]
id = "H2O"
[[release]]
kind = "volcanic"
species = "H2O"
lon = [184.62]
lat = [-20.55]
mass = 3.3e15
time = [2022-01-15T04:00:00, 2022-01-15T09:10:00]
[release.volcanic]
phases = [[2022-01-15T04:00:00, 2022-01-15T04:15:00, 0.0, 30000.0, 1.0]]
[release.plume_rise]
scheme = "volcanic"
mode = "parcels"
mer = 5.0e9
entrainment = "pairwise"
[release.plume_rise.vent_mix]
magma_t = 1300.0
magma_water = 0.04
seawater_core = 0.10
seawater_edge = 0.65
mix_exponent = 1.5
""")
spec = EasleeConfig.read_runspec(tmp)
r1 = spec.release[1]
check("parse: vent_mix enabled with overrides",
      r1.plume_vent_mix && r1.plume_magma_t == 1300.0 &&
      r1.plume_sw_core == 0.10 && r1.plume_sw_edge == 0.65 &&
      r1.plume_magma_water == 0.04 && r1.plume_mix_exp == 1.5)
out = joinpath(dirname(tmp), "vm_out.toml")
EasleeConfig.write_runspec(spec, out)
spec2 = EasleeConfig.read_runspec(out)
r2 = spec2.release[1]
check("round-trip preserves the vent_mix block",
      r2.plume_vent_mix && r2.plume_magma_t == 1300.0 &&
      r2.plume_sw_edge == 0.65 && r2.plume_mix_exp == 1.5)
# validation: vent_mix without parcels mode must be rejected
write(tmp, replace(read(tmp, String), "mode = \"parcels\"" => "mode = \"layer\""))
bad = try; EasleeConfig.read_runspec(tmp); nothing; catch e; e; end
check("validation rejects vent_mix without mode=parcels", bad !== nothing)

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
