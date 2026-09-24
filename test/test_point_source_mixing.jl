"""
test_point_source_mixing.jl
═══════════════════════════
Point-source dispersion test.

Setup
─────
  N_TRAJ = 10 000 Lagrangian parcels drawn from a background Gaussian
  N(0, σ_puff²).  One parcel placed exactly at the origin carries a
  large tracer mixing ratio Q₀ = 1.  All other parcels have q = 0.
  Background: q_bg = 0.

Analytic solution (2-D continuous diffusion)
─────────────────────────────────────────────
  The mixing ratio field from a point source Q₀ at origin spreads as

      c(r, t) = Q₀ / (2π σ²(t)) · exp(−r² / (2 σ²(t)))

  with σ²(t) = 2 K_H t.  The KDE estimate broadens this by h², so we
  compare against the KDE-corrected analytic:

      c_kde(x, 0, t) = Q₀ / (2π (σ²(t) + h²)) · exp(−x² / (2(σ²(t)+h²)))

Tracked quantities
──────────────────
  • Cross-section c(x, 0, t) at several snapshots (KDE vs analytic)
  • Σ qᵢ         — sum of mixing ratios   (conserved in old scheme)
  • Σ qᵢ · mᵢ    — total tracer mass      (should be approx conserved if
                   the exchange is physical; tells us where tracer goes)

Usage
─────
  julia --project=easlee -t 4 test/test_point_source_mixing.jl
"""

using Random, Printf, Statistics
# Plotting is optional: set EASLEE_PLOT=1 to render the summary figure.
# GLMakie is NOT a dependency of the model — install it once in your default
# environment (julia -e 'using Pkg; Pkg.add("GLMakie")') and it is picked up
# through the environment stack.
const SAVE_PLOT = get(ENV, "EASLEE_PLOT", "0") == "1"
SAVE_PLOT && @eval using GLMakie

const REPO = joinpath(@__DIR__, "..")
include(joinpath(REPO, "easlee", "EASLEE.jl"))
using .EASLEE: MixingConfig, apply_mixing_xyz!,
               Sym22, PuffCovarianceState, PuffMixingConfig, apply_puff_mixing!

Random.seed!(42)
mkpath(joinpath(REPO, "results"))

# ─────────────────────────────────────────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────────────────────────────────────────

const K_H       = 1_000.0      # diffusivity [m²/s]
const DT        = 60.0         # time step   [s]
const N_STEPS   = 30           # 30 min total
const Q0        = 1.0          # initial tracer at origin
const N_TRAJ    = 10_000       # total parcels

const SIGMA_PUFF = 10_000.0    # background cloud radius [m]  (>> spreading σ)

# Neighbour kernel width (≈ mean inter-parcel spacing in the puff centre)
# ρ_centre ≈ N/(2π σ_puff²) → spacing ≈ 1/√ρ ≈ 2π^(1/2) σ_puff / √N ≈ 888 m
const DX_NBR     = 1_200.0     # [m] — kernel width for old scheme
const H_KDE      = 1_200.0     # [m] — KDE bandwidth for visualization

# Old-scheme configs — two nsig values to compare
const MCFG3 = MixingConfig(sigma0_h=DX_NBR, K_h=K_H, nsig=3.0, max_neighbors=50)
const MCFG2 = MixingConfig(sigma0_h=DX_NBR, K_h=K_H, nsig=2.0, max_neighbors=50)

# New-scheme config
const RHO_BG      = 1.2
const K_EFF_FAC   = sqrt(2.4)
const DSH_1       = Sym22(2K_H*DT/K_EFF_FAC, 0.0, 2K_H*DT/K_EFF_FAC)
const SIGMA_H0    = DX_NBR
const M0          = RHO_BG * SIGMA_H0^2 * 1.0
const PCFG        = PuffMixingConfig(max_lam_cap=SIGMA_H0^2)

# Snapshot times (steps)
const SNAP_STEPS  = [0, 5, 10, 20, 30]

# ─────────────────────────────────────────────────────────────────────────────
# Initial parcel positions
# ─────────────────────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(7)

# parcel 1: exactly at origin (the point source)
# parcels 2:N_TRAJ: random draw from background Gaussian
x_bg = SIGMA_PUFF .* randn(rng, N_TRAJ - 1)
y_bg = SIGMA_PUFF .* randn(rng, N_TRAJ - 1)
x_all = vcat([0.0], x_bg)
y_all = vcat([0.0], y_bg)
z_all = zeros(N_TRAJ)
act   = trues(N_TRAJ)

# ─────────────────────────────────────────────────────────────────────────────
# 1-D KDE along x-axis  (evaluate at y = 0, summing 2-D Gaussian kernel)
# ─────────────────────────────────────────────────────────────────────────────

const X_EVAL = range(-15_000.0, 15_000.0, length=200)
const X_EVAL_KM = collect(X_EVAL) ./ 1000

"""
    kde_xsection(xs, ys, weights, h) → Vector{Float64}

2-D Gaussian KDE evaluated along the x-axis (y=0):
  ĉ(x, 0) = Σ_k  w_k · G_h(x − x_k, y_k)
(no extra normalization — the kernel 1/(2πh²) is inside G_h)
"""
function kde_xsection(xs, ys, weights, h)
    inv2h2 = 0.5 / h^2
    normk  = 1.0 / (2π * h^2)
    c = zeros(length(X_EVAL))
    for k in eachindex(xs)
        w = weights[k]
        w == 0.0 && continue
        dy2 = ys[k]^2
        dy2 > 9h^2 && continue
        for ix in eachindex(X_EVAL)
            dx = xs[k] - X_EVAL[ix]
            r2 = dx*dx + dy2
            r2 > 9h^2 && continue
            c[ix] += w * exp(-r2 * inv2h2)
        end
    end
    c .* normk
end

# KDE-corrected analytic  (point source convolved with G_h)
function analytic_kde(t, h)
    sig2_eff = max(2K_H*t + h^2, h^2)
    [Q0 / (2π * sig2_eff) * exp(-x^2 / (2*sig2_eff)) for x in X_EVAL]
end

# ─────────────────────────────────────────────────────────────────────────────
# Run both schemes, collecting snapshots
# ─────────────────────────────────────────────────────────────────────────────

q0_vec = zeros(N_TRAJ, 1)
q0_vec[1, 1] = Q0       # single point source at origin

# Storage: Dict(step → (x-section vector, Σq, Σqm, m_vec))
snaps_old3 = Dict{Int, NamedTuple}()
snaps_old2 = Dict{Int, NamedTuple}()
snaps_new  = Dict{Int, NamedTuple}()

# ── Old scheme nsig=3.0 ──────────────────────────────────────────────────────
println("Running old scheme nsig=3.0 (apply_mixing_xyz!) …")
let
    q  = copy(q0_vec)
    yl = fill(1.0, N_TRAJ)
    for step in 0:N_STEPS
        if step in SNAP_STEPS
            xs = kde_xsection(x_all, y_all, q[:, 1], H_KDE)
            snaps_old3[step] = (xsec=xs, sumq=sum(q))
        end
        step == N_STEPS && break
        apply_mixing_xyz!(q, copy(x_all), copy(y_all), copy(z_all),
                          yl, act, DT, MCFG3; merge_packed=false)
    end
end

# ── Old scheme nsig=2.0 ──────────────────────────────────────────────────────
println("Running old scheme nsig=2.0 (apply_mixing_xyz!) …")
let
    q  = copy(q0_vec)
    yl = fill(1.0, N_TRAJ)
    for step in 0:N_STEPS
        if step in SNAP_STEPS
            xs = kde_xsection(x_all, y_all, q[:, 1], H_KDE)
            snaps_old2[step] = (xsec=xs, sumq=sum(q))
        end
        step == N_STEPS && break
        apply_mixing_xyz!(q, copy(x_all), copy(y_all), copy(z_all),
                          yl, act, DT, MCFG2; merge_packed=false)
    end
end

# ── New scheme ──────────────────────────────────────────────────────────────
println("Running new scheme (apply_puff_mixing!) …")
let
    q    = copy(q0_vec)
    cov  = PuffCovarianceState(N_TRAJ; sigma_H0=SIGMA_H0, sigma_z0=1.0, m0=M0)
    dSH  = fill(DSH_1, N_TRAJ)
    dsz2 = zeros(N_TRAJ)
    q_bg = [0.0]

    for step in 0:N_STEPS
        if step in SNAP_STEPS
            xs   = kde_xsection(x_all, y_all, q[:, 1], H_KDE)
            sumq = sum(q)
            sumqm = sum(q[i,1] * cov.m[i] for i in 1:N_TRAJ)
            snaps_new[step] = (xsec=xs, sumq=sumq, sumqm=sumqm,
                               m_center=cov.m[1], q_center=q[1,1])
        end
        step == N_STEPS && break
        apply_puff_mixing!(cov, q, dSH, dsz2, x_all, y_all, z_all,
                           act, RHO_BG, q_bg, PCFG)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary table
# ─────────────────────────────────────────────────────────────────────────────
println()
println("=" ^ 80)
println("  POINT-SOURCE DISPERSION  (Q₀=$Q0 at origin, N_TRAJ=$N_TRAJ parcels)")
println("  K_H=$(Int(K_H)) m²/s  dt=$(Int(DT)) s  σ_puff=$(Int(SIGMA_PUFF/1000)) km")
println("=" ^ 80)
@printf("  %-6s  %-8s  %-12s  %-12s  %-10s  %-10s  %-10s\n",
        "step", "t[min]",
        "OLD3 Σq", "OLD2 Σq", "NEW Σq", "NEW Σqm", "analytic")
for step in SNAP_STEPS
    t      = step * DT
    o3     = snaps_old3[step]
    o2     = snaps_old2[step]
    nw     = snaps_new[step]
    sig2t  = 2K_H*t
    c_peak = t > 0 ? Q0/(2π*sig2t) : Inf
    @printf("  %-6d  %-8.1f  %-12.4f  %-12.4f  %-10.4f  %-10.4f  %-10s\n",
            step, t/60,
            o3.sumq, o2.sumq, nw.sumq, nw.sumqm,
            t > 0 ? @sprintf("%.4e", c_peak) : "∞")
end
println("=" ^ 80)
println()
println("  Σq = sum of mixing ratios (= Q₀ if conserved)")
println("  Σqm = total tracer mass Σᵢ qᵢ·mᵢ")

# ─────────────────────────────────────────────────────────────────────────────
# Figure (optional — numeric checks above are complete)
# ─────────────────────────────────────────────────────────────────────────────
# ── assertions (review 2026-09-23: the test printed a table and always passed) ──
# The peer-to-peer schemes conserve the sum of mixing ratios exactly; the puff scheme
# dilutes into the background by design, so its Σq drifts — bounded here as a regression
# guard on the values recorded on 2026-09-23 (1.070 at 30 min).
nfail = 0
for step in SNAP_STEPS
    for (lab, sn) in (("OLD3", snaps_old3[step]), ("OLD2", snaps_old2[step]))
        abs(sn.sumq - Q0) < 1e-6 || (global nfail += 1; println("FAIL: $lab Σq = $(sn.sumq) at step $step (expected $Q0)"))
    end
    nw = snaps_new[step]
    (isfinite(nw.sumq) && Q0 <= nw.sumq < 1.10 * Q0) || (global nfail += 1; println("FAIL: NEW Σq = $(nw.sumq) at step $step (expected within [1, 1.10] × Q0)"))
    (isfinite(nw.sumqm) && nw.sumqm > 0) || (global nfail += 1; println("FAIL: NEW Σqm not finite/positive at step $step"))
end
println(nfail == 0 ? "ALL PASS (point-source mixing: $(length(SNAP_STEPS)) snapshots checked)" : "FAIL: $nfail checks")
nfail == 0 || exit(1)

SAVE_PLOT || exit(0)

colors_steps = [:grey70, :dodgerblue, :steelblue, :royalblue, :navy]

fig = Figure(size=(1400, 1000), fontsize=13)

Label(fig[0, :],
      "Point-source dispersion  |  Q₀=$Q0 at origin  |  " *
      "N_TRAJ=$N_TRAJ in background Gaussian (σ=$(Int(SIGMA_PUFF/1000)) km)  |  " *
      "K_H=$(Int(K_H)) m²/s  dt=$(Int(DT)) s",
      fontsize=14, font=:bold)

# ── Row 1: Cross-sections — old nsig=3 | old nsig=2 | new scheme ─────────────
ax1 = Axis(fig[1, 1:2],
           title="Old  nsig=3.0  — cross-section  c(x,0,t)",
           xlabel="x [km]", ylabel="mixing ratio  q")
ax2 = Axis(fig[1, 3:4],
           title="Old  nsig=2.0  — cross-section  c(x,0,t)",
           xlabel="x [km]", ylabel="mixing ratio  q")
ax_new = Axis(fig[1, 5:6],
              title="New scheme  (puff_mixing!)  — cross-section  c(x,0,t)",
              xlabel="x [km]", ylabel="mixing ratio  q")

for (k, step) in enumerate(SNAP_STEPS)
    t   = step * DT
    col = colors_steps[k]
    ana = analytic_kde(t, H_KDE)
    lbl = "t=$(Int(t÷60)) min"
    lines!(ax1, X_EVAL_KM, snaps_old3[step].xsec; color=col, linewidth=2.5, label=lbl)
    lines!(ax1, X_EVAL_KM, ana;                    color=col, linewidth=1, linestyle=:dash)
    lines!(ax2, X_EVAL_KM, snaps_old2[step].xsec; color=col, linewidth=2.5, label=lbl)
    lines!(ax2, X_EVAL_KM, ana;                    color=col, linewidth=1, linestyle=:dash)
    lines!(ax_new, X_EVAL_KM, snaps_new[step].xsec; color=col, linewidth=2.5, label=lbl)
    lines!(ax_new, X_EVAL_KM, ana;                   color=col, linewidth=1, linestyle=:dash)
end
for ax in (ax1, ax2, ax_new)
    lines!(ax, [NaN],[NaN]; color=:black, linewidth=2.5, label="Scheme (KDE)")
    lines!(ax, [NaN],[NaN]; color=:black, linewidth=1, linestyle=:dash, label="Analytic")
    axislegend(ax; position=:rt, labelsize=9, nbanks=2)
end

# ── Row 2: Normalised cross-sections ─────────────────────────────────────────
ax3 = Axis(fig[2, 1:2], title="Old nsig=3 — normalised",
           xlabel="x [km]", ylabel="q / peak(q)")
ax4 = Axis(fig[2, 3:4], title="Old nsig=2 — normalised",
           xlabel="x [km]", ylabel="q / peak(q)")
ax5 = Axis(fig[2, 5:6], title="New — normalised",
           xlabel="x [km]", ylabel="q / peak(q)")

for (k, step) in enumerate(SNAP_STEPS)
    step == 0 && continue
    t   = step * DT
    col = colors_steps[k]
    ana = analytic_kde(t, H_KDE);  pk_a = maximum(ana)
    lbl = "t=$(Int(t÷60)) min"
    for (ax, snaps) in ((ax3, snaps_old3), (ax4, snaps_old2), (ax5, snaps_new))
        xs  = snaps[step].xsec;  pk = maximum(xs)
        pk > 0 && lines!(ax, X_EVAL_KM, xs./pk;    color=col, linewidth=2.5, label=lbl)
        pk_a > 0 && lines!(ax, X_EVAL_KM, ana./pk_a; color=col, linewidth=1, linestyle=:dash)
    end
end
for ax in (ax3, ax4, ax5);  axislegend(ax; position=:rt, labelsize=9);  end

# ── Row 3: Conservation + central parcel ─────────────────────────────────────
ax6 = Axis(fig[3, 1:3],
           title="Tracer conservation (Σq / Q₀)",
           xlabel="time [min]", ylabel="Σq / Q₀")

t_steps   = [s*DT/60 for s in SNAP_STEPS]
sumq_old3 = [snaps_old3[s].sumq for s in SNAP_STEPS]
sumq_old2 = [snaps_old2[s].sumq for s in SNAP_STEPS]
sumq_new  = [snaps_new[s].sumq  for s in SNAP_STEPS]
sumqm_new = [snaps_new[s].sumqm for s in SNAP_STEPS]

hlines!(ax6, [1.0]; color=(:black,0.3), linewidth=1, linestyle=:dot, label="Q₀")
lines!(ax6, t_steps, sumq_old3./Q0; color=:steelblue, linewidth=2.5, label="Old nsig=3  Σq")
lines!(ax6, t_steps, sumq_old2./Q0; color=:dodgerblue, linewidth=2.5,
       linestyle=:dash, label="Old nsig=2  Σq")
lines!(ax6, t_steps, sumq_new./Q0;  color=:crimson, linewidth=2.5, label="New  Σq")
lines!(ax6, t_steps, sumqm_new./(Q0*M0); color=:orangered, linewidth=2.5,
       linestyle=:dashdot, label="New  Σqm/(Q₀m₀)")
for (y, col) in ((sumq_old3./Q0,:steelblue),(sumq_old2./Q0,:dodgerblue),
                 (sumq_new./Q0,:crimson),(sumqm_new./(Q0*M0),:orangered))
    scatter!(ax6, t_steps, y; color=col, markersize=7)
end
axislegend(ax6; position=:lb, labelsize=10)

ax7 = Axis(fig[3, 4:6],
           title="KDE at x=0 vs analytic peak  (old nsig=2 vs nsig=3)",
           xlabel="time [min]", ylabel="c(0,0,t) / Q₀")

t_cont   = range(1, N_STEPS*DT, length=200)
ana_peak = [Q0/(2π*(2K_H*t + H_KDE^2)) for t in t_cont]
mid      = div(length(X_EVAL), 2)

q_ctr3 = [snaps_old3[s].xsec[mid] for s in SNAP_STEPS]
q_ctr2 = [snaps_old2[s].xsec[mid] for s in SNAP_STEPS]

lines!(ax7, t_cont./60, ana_peak./Q0; color=:black, linewidth=2,
       linestyle=:dash, label="Analytic (KDE-corr.)")
lines!(ax7, t_steps, q_ctr3./Q0; color=:steelblue, linewidth=2.5, label="Old nsig=3 KDE")
lines!(ax7, t_steps, q_ctr2./Q0; color=:dodgerblue, linewidth=2.5,
       linestyle=:dash, label="Old nsig=2 KDE")
scatter!(ax7, t_steps, q_ctr3./Q0; color=:steelblue, markersize=7)
scatter!(ax7, t_steps, q_ctr2./Q0; color=:dodgerblue, markersize=7)
axislegend(ax7; position=:rt, labelsize=10)

# ─────────────────────────────────────────────────────────────────────────────
outpath = joinpath(REPO, "results", "point_source_mixing.png"); mkpath(dirname(outpath))
mkpath(dirname(outpath))
save(outpath, fig; px_per_unit=2)
println("Saved → $outpath")
