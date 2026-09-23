#!/usr/bin/env julia
# test_plume_phase_v2.jl — Level-2 v2 acceptance test (VOLCANIC_PLUME.md):
# PAIRWISE-EXCHANGE entrainment (the "Dirac in a crowd" configuration).
# Parcels carry a radial position; only the outer shell engulfs ambient air
# (rate λ/f_shell) and interior parcels mix pairwise with r-adjacent
# neighbours.  Verify:
#   [1] the pairwise pass conserves mass-weighted totals EXACTLY (unit check),
#   [2] a Dirac tracer seeded on the axis spreads through the crowd while its
#       total is conserved (the reference configuration the design is built on),
#   [3] radial structure EMERGES: hot, source-rich core vs cool, diluted edge
#       at the same heights mid-rise — the geometric content v1 cannot have,
#   [4] the column-mean behaviour still matches the 1-D bookkeeping
#       (provenance, detrainment band) — pairwise ≠ different bulk physics,
#   [5] washout still works and the water budget still closes.
# Writes results/plume_v2_snapshots.csv for the cross-section figure
# (scripts/plot_plume_v2.py).
#
#   julia --project=easlee test/test_plume_phase_v2.jl
using Printf, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: build_plume_profiles, PlumeParcels, PlumeWashout,
                         plume_parcel_init!, plume_parcels_substep_v2!,
                         plume_umbrella_equilibrate!, plume_phase_stats,
                         volcanic_column_provenance
using .EASLEE.Operators: _plume_pairwise_pass!, _pp_interp

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── standard tropical atmosphere (same as the v1 tests) ───────────────────────
z = collect(0.0:200.0:60_000.0); T = similar(z)
for (i, zz) in enumerate(z)
    T[i] = zz <= 16_000 ? 300.0 - 6.5e-3*zz : zz <= 20_000 ? 196.0 : 196.0 + 2.2e-3*(zz-20_000)
end
p = similar(z); p[1] = 101_325.0
for i in 2:length(z); p[i] = p[i-1]*exp(-9.81*(z[i]-z[i-1])/(287.05*0.5*(T[i-1]+T[i]))); end
qv = [max(0.6*0.622*611.2*exp(17.67*(Tt-273.15)/(Tt-29.65))/pp, 2e-6) for (Tt,pp) in zip(T,p)]

MER = 5.0e9; T0 = 700.0; NW = 0.45; ALPHA = 0.10
prof = build_plume_profiles(MER, T0, NW, z, T, p, qv; alpha=ALPHA)
_, _, _, _, fsrc1d, flt1d, fft1d, fst1d =
    volcanic_column_provenance(MER, T0, NW, z, T, p, qv; alpha=ALPHA)
@printf("1-D column reference: top %.1f km, NBH %.1f km\n", prof.z_top/1e3, prof.z_nbh/1e3)

println()
println("="^78)
println("[1] pairwise pass conserves mass-weighted totals exactly")
println("="^78)
Random.seed!(7)
pu = PlumeParcels(64)
for i in 1:64
    plume_parcel_init!(pu, i, prof; T0=T0, nw0=NW)
    pu.z[i] = 5_000.0 + 400.0 * rand()          # same z-bin neighbourhood
    pu.T[i] = 300.0 + 200.0 * rand()
    pu.qv[i] = 0.5 * rand(); pu.m[i] = 0.5 + 2.0 * rand()
    pu.fso2[i] = rand()
end
tot(a) = sum(a[i] * pu.m[i] for i in 1:64)
t_T = tot(pu.T); t_q = tot(pu.qv); t_w = tot(pu.w); t_g = sum(pu.fso2)
std0 = std(pu.T)
for _ in 1:200
    _plume_pairwise_pass!(pu, 1.0)
end
check("Σ m·T conserved",  abs(tot(pu.T) - t_T) < 1e-8 * abs(t_T))
check("Σ m·qv conserved", abs(tot(pu.qv) - t_q) < 1e-8 * abs(t_q))
check("Σ m·w conserved",  abs(tot(pu.w) - t_w) < 1e-6 * max(abs(t_w), 1.0))
check("Σ gas load conserved (concentration exchange)",
      abs(sum(pu.fso2) - t_g) < 1e-8 * abs(t_g))
# exchange is r-LOCAL (adjacent partners), so relaxation is diffusive along
# the r-chain, not instantaneous — require substantial, not total, decay
check("pair relaxation homogenises the bin (std halves in 200 passes)",
      std(pu.T) < 0.55 * std0;
      info=@sprintf("(std %.1f → %.1f K)", std0, std(pu.T)))

println()
println("="^78)
println("[2] Dirac tracer: axis-seeded tracer spreads, total conserved")
println("="^78)
Random.seed!(20220115)
N = 2000
pd = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(pd, i, prof; T0=T0, nw0=NW)
    pd.fso2[i] = pd.r[i] < 0.15 ? 1.0 : 0.0     # Dirac: tracer only on the axis
end
n_seed = count(>(0.0), pd.fso2)
g0 = sum(pd.fso2)
t = 0.0
while any(pd.active) && t < 3600.0              # NO washout → conservation exact
    plume_parcels_substep_v2!(pd, prof, 30.0)
    global t += 30.0
end
n_reached = count(>(1e-3), pd.fso2)
check("tracer total conserved through the full rise",
      abs(sum(pd.fso2) - g0) < 1e-6 * g0;
      info=@sprintf("(Σ %.4f vs %.4f)", sum(pd.fso2), g0))
check("tracer spreads far beyond the seeded parcels",
      n_reached > 3 * n_seed; info="($n_seed seeded → $n_reached carry >1e-3)")
check("no parcel exceeds the seed concentration (maximum principle)",
      maximum(pd.fso2 ./ pd.m) <= 1.0 + 1e-9)

println()
println("="^78)
println("[3] radial structure mid-rise: hot source-rich core, cool diluted edge")
println("="^78)
Random.seed!(20220115)
ps = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(ps, i, prof; T0=T0, nw0=NW)
end
plume_parcels_substep_v2!(ps, prof, 120.0)      # 2 min into the rise
act = findall(ps.active)
# compare core vs edge AT THE SAME HEIGHTS (z-bin-matched means)
dT = Float64[]; dS = Float64[]
for zb in 0.0:2_000.0:40_000.0
    sel = [i for i in act if zb <= ps.z[i] < zb + 2_000.0]
    core = [i for i in sel if ps.r[i] < 0.4]
    edge = [i for i in sel if ps.r[i] > 0.8]
    (length(core) >= 10 && length(edge) >= 10) || continue
    push!(dT, mean(ps.T[core]) - mean(ps.T[edge]))
    push!(dS, mean(ps.fsrc[core]) - mean(ps.fsrc[edge]))
end
check("z-matched core–edge contrast sampled in ≥3 height bins",
      length(dT) >= 3; info="($(length(dT)) bins)")
check("core is HOTTER than the edge (mean over bins > 5 K)",
      mean(dT) > 5.0; info=@sprintf("(mean ΔT %.1f K)", mean(dT)))
check("core is more SOURCE-RICH than the edge",
      mean(dS) > 0.02; info=@sprintf("(mean Δf_src %.3f)", mean(dS)))

println()
println("="^78)
println("[4] column-mean behaviour still matches the 1-D bookkeeping")
println("="^78)
Random.seed!(20220115)
pp = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(pp, i, prof; T0=T0, nw0=NW)
end
wash = PlumeWashout()
t2 = 0.0
while any(pp.active) && t2 < 3600.0
    plume_parcels_substep_v2!(pp, prof, 30.0; washout=wash)
    global t2 += 30.0
end
scol = plume_phase_stats(pp)
plume_umbrella_equilibrate!(pp, prof; washout=wash)
s = plume_phase_stats(pp)
@printf("  detrained %d/%d in %.0f min:  z p10/p50/p90 = %.1f / %.1f / %.1f km\n",
        s.n, N, t2/60, s.z_p10/1e3, s.z_p50/1e3, s.z_p90/1e3)
@printf("  ensemble: src %.3f  LT %.3f  FT %.3f  ST %.3f\n", s.fsrc, s.flt, s.fft, s.fst)
@printf("  1-D     : src %.3f  LT %.3f  FT %.3f  ST %.3f\n", fsrc1d, flt1d, fft1d, fst1d)
check("all trajectories detrain within an hour", s.n == N; info="($(s.n)/$N)")
check("median detrainment inside the umbrella band",
      prof.z_nbh <= s.z_p50 <= prof.z_top * 1.02;
      info=@sprintf("(%.1f km)", s.z_p50/1e3))
for (name, a, b) in (("source", s.fsrc, fsrc1d), ("lower-trop", s.flt, flt1d),
                     ("free-trop", s.fft, fft1d), ("strat", s.fst, fst1d))
    check("$name fraction within ±0.10 of the 1-D value", abs(a - b) <= 0.10;
          info=@sprintf("(%.2f vs %.2f)", a, b))
end
check("fractions sum to 1", abs(s.fsrc + s.flt + s.fft + s.fst - 1) < 1e-6)

println()
println("="^78)
println("[5] washout still works under pairwise exchange")
println("="^78)
w_in  = NW + s.w_entrained
w_out = s.w_carried + s.w_fallen
resid = abs(w_in - w_out) / w_in
@printf("  water: lofted %.3f = carried %.4f + fallen %.3f  |  SO₂ %.2f  HCl %.3f\n",
        w_in, s.w_carried, s.w_fallen, s.fso2, s.fhcl)
check("water budget closes to <1 %", resid < 0.01;
      info=@sprintf("(residual %.2f %%)", 100resid))
check("majority of the lofted water washed out", s.w_fallen / w_in > 0.5;
      info=@sprintf("(%.0f %%)", 100 * s.w_fallen / w_in))
check("HCl scavenged ≫ SO₂ (ordering preserved)",
      s.fso2 - s.fhcl > 0.3; info=@sprintf("(%.2f vs %.3f)", s.fso2, s.fhcl))
check("no NaNs; r in [0,1]",
      !any(isnan, pp.z) && !any(isnan, pp.T) && all(0 .<= pp.r .<= 1))

# ── snapshots for the cross-section figure ────────────────────────────────────
Random.seed!(20220115)
pv = PlumeParcels(N)
for i in 1:N
    plume_parcel_init!(pv, i, prof; T0=T0, nw0=NW)
end
resdir = joinpath(REPO, "results"); mkpath(resdir)
open(joinpath(resdir, "plume_v2_snapshots.csv"), "w") do io
    println(io, "t_s,r,z_m,T,Ta,fsrc,active")
    tprev = 0.0
    for tsnap in (60.0, 150.0, 300.0, 600.0)
        plume_parcels_substep_v2!(pv, prof, tsnap - tprev)
        tprev = tsnap
        for i in 1:N
            (pv.active[i] || pv.done[i]) || continue
            Ta = _pp_interp(prof, prof.Ta, pv.z[i])
            @printf(io, "%.0f,%.4f,%.1f,%.1f,%.1f,%.4f,%d\n",
                    tsnap, pv.r[i], pv.z[i], pv.T[i], Ta, pv.fsrc[i], pv.active[i])
        end
    end
end
println("\nCSV → results/plume_v2_snapshots.csv")

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
