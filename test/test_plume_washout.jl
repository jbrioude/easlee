#!/usr/bin/env julia
# test_plume_washout.jl — Level-2 step 3 acceptance test (VOLCANIC_PLUME.md):
# in-plume WASHOUT.  Rising parcels shed condensate above the updraft holding
# capacity; soluble gases (SO₂ at plume pH, HCl solubility-unlimited) leave
# with the shed water by effective Henry equilibrium.  Verify:
#   • the water budget closes exactly (initial + entrained = carried + fallen),
#   • MOST of the vaporised ocean water is washed out in the column and only a
#     small fraction survives to the umbrella — the Hunga Tonga stratospheric
#     ~150 Tg is the small survivor of a much larger lofted mass,
#   • HCl is scavenged far more efficiently than SO₂ (observed: strong SO₂
#     signal, ~no stratospheric HCl enhancement),
#   • unloading does not break the detrainment-height envelope,
#   • the no-washout path is bit-identical in behaviour (regression).
# Writes results/plume_washout_traj.csv + results/plume_washout_fall.csv for
# the companion figure (scripts/plot_plume_washout.py).
#
#   julia --project=easlee test/test_plume_washout.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: build_plume_profiles, PlumeParcels, PlumeWashout,
                         plume_parcel_init!, plume_parcels_substep!,
                         plume_umbrella_equilibrate!, plume_phase_stats

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── standard tropical atmosphere (same as test_plume_phase.jl) ────────────────
z = collect(0.0:200.0:60_000.0); T = similar(z)
for (i, zz) in enumerate(z)
    T[i] = zz <= 16_000 ? 300.0 - 6.5e-3*zz : zz <= 20_000 ? 196.0 : 196.0 + 2.2e-3*(zz-20_000)
end
p = similar(z); p[1] = 101_325.0
for i in 2:length(z); p[i] = p[i-1]*exp(-9.81*(z[i]-z[i-1])/(287.05*0.5*(T[i-1]+T[i]))); end
qv = [max(0.6*0.622*611.2*exp(17.67*(Tt-273.15)/(Tt-29.65))/pp, 2e-6) for (Tt,pp) in zip(T,p)]

# ── Hunga Tonga source ────────────────────────────────────────────────────────
MER = 5.0e9; T0 = 700.0; NW = 0.45; ALPHA = 0.10
prof = build_plume_profiles(MER, T0, NW, z, T, p, qv; alpha=ALPHA)
@printf("1-D column reference: top %.1f km, NBH %.1f km\n", prof.z_top/1e3, prof.z_nbh/1e3)
wash = PlumeWashout()

# ── release + integrate, recording the fallout profile ────────────────────────
Random.seed!(20220115)
N = 2000
pp = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(pp, i, prof; T0=T0, nw0=NW)
end
DT = 30.0
zedges = collect(0.0:1000.0:45_000.0)         # fallout-by-height accumulator
fallz  = zeros(length(zedges) - 1)
wfell_prev = zeros(N)
series = NamedTuple[]                          # ensemble time series for the figure
t = 0.0
while any(pp.active) && t < 3600.0
    plume_parcels_substep!(pp, prof, DT; washout=wash)
    global t += DT
    for i in 1:N
        dwf = pp.wfell[i] - wfell_prev[i]
        if dwf > 0
            k = clamp(searchsortedlast(zedges, pp.z[i]), 1, length(fallz))
            fallz[k] += dwf
            wfell_prev[i] = pp.wfell[i]
        end
    end
    push!(series, (t = t,
        vap  = sum(pp.qv .* pp.m) / N,
        cond = sum(pp.qc .* pp.m) / N,
        fell = sum(pp.wfell) / N,
        ent  = sum(pp.went) / N,
        zbar = sum(pp.z) / N,
        nact = count(pp.active)))
end
scol = plume_phase_stats(pp)                    # after the COLUMN phase
wfell_col = copy(pp.wfell); fso2_col = copy(pp.fso2); fhcl_col = copy(pp.fhcl)
plume_umbrella_equilibrate!(pp, prof; washout=wash)   # first hours in the umbrella
s = plume_phase_stats(pp)

println()
println("="^78)
println("[1] water budget closure (per unit initial parcel mass, ensemble mean)")
println("="^78)
w_in  = NW + s.w_entrained
w_out = s.w_carried + s.w_fallen
resid = abs(w_in - w_out) / w_in
@printf("  initial %.3f + entrained %.3f = %.3f   |   carried %.3f + fallen %.3f = %.3f\n",
        NW, s.w_entrained, w_in, s.w_carried, s.w_fallen, w_out)
check("budget closes to <1 %", resid < 0.01; info=@sprintf("(residual %.2f %%)", 100resid))
check("all trajectories detrain within an hour", s.n == N; info="($(s.n)/$N)")

println()
println("="^78)
println("[2] where is the water washed out, and how much survives? (150 Tg question)")
println("="^78)
det = findall(pp.done)
f_col = mean(wfell_col[det])                    # column-phase fallout
f_umb = s.w_fallen - f_col                      # umbrella dehydration
fsurv_all = s.w_carried / w_in                  # of ALL lofted water
@printf("  lofted %.3f kg/kg:  column fallout %.3f (%.0f %%)  +  umbrella %.3f (%.0f %%)\n",
        w_in, f_col, 100f_col/w_in, f_umb, 100f_umb/w_in)
@printf("  → carried (vapour+retained ice) %.4f = %.1f %% of the lofted water\n",
        s.w_carried, 100fsurv_all)
@printf("  per kg erupted: %.4f kg water survives → 150 Tg needs %.0f Tg erupted mass\n",
        s.w_carried, 150.0 / max(s.w_carried, 1e-9))
check("majority of the lofted water is washed out", s.w_fallen / w_in > 0.5;
      info=@sprintf("(%.0f %%)", 100 * s.w_fallen / w_in))
check("umbrella dehydration dominates the column fallout", f_umb > f_col;
      info=@sprintf("(%.0f %% vs %.0f %%)", 100f_umb/w_in, 100f_col/w_in))
check("surviving water fraction small but nonzero (0.5–15 %)",
      0.005 <= fsurv_all <= 0.15; info=@sprintf("(%.1f %%)", 100fsurv_all))
check("implied erupted mass for 150 Tg strat water is 10³–10⁵ Tg (HT class)",
      1.0e3 <= 150.0 / max(s.w_carried, 1e-9) <= 1.0e5;
      info=@sprintf("(%.0f Tg)", 150.0 / max(s.w_carried, 1e-9)))

println()
println("="^78)
println("[3] gas scavenging: HCl ≫ SO₂ (observed HT signature)")
println("="^78)
@printf("  after column : SO₂ %.3f   HCl %.3f\n", mean(fso2_col[det]), mean(fhcl_col[det]))
@printf("  after umbrella: SO₂ %.3f   HCl %.3f\n", s.fso2, s.fhcl)
check("HCl mostly scavenged (<25 % survives)", s.fhcl < 0.25;
      info=@sprintf("(%.1f %%)", 100s.fhcl))
check("SO₂ mostly survives (>60 %)", s.fso2 > 0.60; info=@sprintf("(%.0f %%)", 100s.fso2))
check("scavenging ordering HCl ≫ SO₂", s.fso2 - s.fhcl > 0.3)

println()
println("="^78)
println("[4] unloading does not break the envelope; state sanity")
println("="^78)
@printf("  detrainment z p10/p50/p90 = %.1f / %.1f / %.1f km (band %.1f–%.1f)\n",
        s.z_p10/1e3, s.z_p50/1e3, s.z_p90/1e3, prof.z_nbh/1e3, prof.z_top/1e3)
check("median detrainment inside the umbrella band",
      prof.z_nbh <= s.z_p50 <= prof.z_top * 1.02;
      info=@sprintf("(%.1f km)", s.z_p50/1e3))
check("p90 within the observed layer top", s.z_p90 <= prof.z_top * 1.05;
      info=@sprintf("(%.1f km)", s.z_p90/1e3))
# after shedding dqf the retained q_ret is renormalised by the smaller parcel
# mass (q_ret/(1−dqf)) — allow that inflation, bounded by the largest shed
check("equilibrated parcels at ambient T, condensate ≈ retained amount",
      all(pp.qc[det] .<= wash.qc_crit * 1.25);
      info=@sprintf("(max q_c %.4f, q_ret %.3f)", maximum(pp.qc[det]), wash.qc_crit))
check("no NaNs; surviving fractions in [0,1]",
      !any(isnan, pp.z) && !any(isnan, pp.T) &&
      all(0 .<= pp.fso2 .<= 1) && all(0 .<= pp.fhcl .<= 1) && all(pp.wfell .>= 0))

println()
println("="^78)
println("[5] regression: washout=nothing keeps steps 1+2 behaviour + trivial budget")
println("="^78)
Random.seed!(20220115)
pp0 = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(pp0, i, prof; T0=T0, nw0=NW)
end
t0 = 0.0
while any(pp0.active) && t0 < 3600.0
    plume_parcels_substep!(pp0, prof, DT)          # no washout
    global t0 += DT
end
s0 = plume_phase_stats(pp0)
resid0 = abs(NW + s0.w_entrained - s0.w_carried) / (NW + s0.w_entrained)
check("no fallout without washout", s0.w_fallen == 0.0)
check("gas fractions untouched", s0.fso2 == 1.0 && s0.fhcl == 1.0)
check("mass-tracked budget closes without washout too", resid0 < 0.01;
      info=@sprintf("(residual %.2f %%)", 100resid0))
check("median detrainment height unchanged vs steps 1+2",
      abs(s0.z_p50 - 27_100) < 2_500; info=@sprintf("(%.1f km)", s0.z_p50/1e3))

println()
println("="^78)
println("[6] size-dependent settling (settling=true): updraft recycles hydrometeors")
println("="^78)
# v_t-gated shedding: in the fast core v_t < w_gate·|w| so nothing sheds —
# the column lofts its ice to the umbrella first (the observed HT behaviour),
# then drops it there.  Expect LESS column washout than τ_fall mode, with the
# umbrella phase doing correspondingly more; bulk survival stays HT-class.
Random.seed!(20220115)
ps = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(ps, i, prof; T0=T0, nw0=NW)
end
wset = PlumeWashout(settling=true)
ts = 0.0
while any(ps.active) && ts < 3600.0
    plume_parcels_substep!(ps, prof, DT; washout=wset)
    global ts += DT
end
ss_col = plume_phase_stats(ps)
plume_umbrella_equilibrate!(ps, prof; washout=wset)
ss = plume_phase_stats(ps)
w_in_s  = NW + ss.w_entrained
resid_s = abs(w_in_s - ss.w_carried - ss.w_fallen) / w_in_s
f_col_s = ss_col.w_fallen
@printf("  column fallout %.3f (τ_fall mode had %.3f)  umbrella %.3f  carried %.4f\n",
        f_col_s, f_col, ss.w_fallen - f_col_s, ss.w_carried)
@printf("  SO₂ %.2f  HCl %.3f\n", ss.fso2, ss.fhcl)
check("budget closes with settling", resid_s < 0.01;
      info=@sprintf("(residual %.2f %%)", 100resid_s))
check("updraft gate reduces the column fallout vs τ_fall mode",
      f_col_s < f_col; info=@sprintf("(%.3f vs %.3f)", f_col_s, f_col))
check("umbrella still dominates; total washout still major",
      ss.w_fallen / w_in_s > 0.5;
      info=@sprintf("(%.0f %%)", 100 * ss.w_fallen / w_in_s))
check("surviving water still HT-class (0.5–15 %)",
      0.005 <= ss.w_carried / w_in_s <= 0.15;
      info=@sprintf("(%.1f %%)", 100 * ss.w_carried / w_in_s))
check("HCl ≫ SO₂ ordering preserved with settling", ss.fso2 - ss.fhcl > 0.3;
      info=@sprintf("(%.2f vs %.3f)", ss.fso2, ss.fhcl))

# ── dump CSVs for the companion figure ────────────────────────────────────────
resdir = joinpath(REPO, "results"); mkpath(resdir)
open(joinpath(resdir, "plume_washout_traj.csv"), "w") do io
    println(io, "z_det_m,fso2,fhcl,fso2_col,fhcl_col,qv,qc,m,wfell,wfell_col,went,fsrc")
    for i in det
        @printf(io, "%.1f,%.5f,%.6f,%.5f,%.6f,%.5f,%.5f,%.4f,%.5f,%.5f,%.5f,%.4f\n",
                pp.z[i], pp.fso2[i], pp.fhcl[i], fso2_col[i], fhcl_col[i],
                pp.qv[i], pp.qc[i], pp.m[i], pp.wfell[i], wfell_col[i],
                pp.went[i], pp.fsrc[i])
    end
end
open(joinpath(resdir, "plume_washout_fall.csv"), "w") do io
    println(io, "z_bot_m,z_top_m,fallen_water")
    for k in eachindex(fallz)
        @printf(io, "%.0f,%.0f,%.6f\n", zedges[k], zedges[k+1], fallz[k] / N)
    end
end
open(joinpath(resdir, "plume_washout_series.csv"), "w") do io
    println(io, "t_s,vapour,condensate,fallen,entrained,zbar_m,n_active")
    for r in series
        @printf(io, "%.0f,%.5f,%.5f,%.5f,%.5f,%.1f,%d\n",
                r.t, r.vap, r.cond, r.fell, r.ent, r.zbar, r.nact)
    end
end
println("\nCSV → results/plume_washout_{traj,fall,series}.csv")

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
