"""
test_two_puff_mixing.jl
═══════════════════════
Two-puff mixing comparison.  Each puff is represented by N_TRAJ = 10 000
randomly-sampled Lagrangian parcels drawn from the initial Gaussian.

Four methods compared
─────────────────────
  1. Analytic         exact PDE solution (two independent Gaussians)
  2. Stochastic       10 000 Brownian-motion walkers per puff (no mixing scheme)
  3. Old scheme       apply_mixing_xyz!  on 10 000 random parcels per puff
  4. New scheme       apply_puff_mixing! on the same parcels

Key questions
─────────────
  • Does each mixing scheme reproduce the analytic mixing fraction
    f(x,t) = 1/(1 + exp(2·d·x/σ²(t)))  between the two puffs?
  • Does each scheme mix with the background (q_bg = 0) or only with
    the other parcels in the ensemble?

Physical setup
──────────────
  Puff centres: (x₁,0) = (−7500 m, 0),  (x₂,0) = (+7500 m, 0)
  σ₀ = 5 000 m,  K_H = 1 000 m²/s,  dt = 60 s,  N_steps = 20

Concentration estimation from parcels
──────────────────────────────────────
  After mixing the parcel positions do not change (pure mixing, no advection).
  The concentration field is recovered by a Gaussian KDE:
      ĉ₁(x,y) = Σ_k  q[k,1] · G_h(x−x_k, y−y_k)
  where G_h is a normalised 2-D Gaussian kernel with bandwidth h.
  This is normalised so that the initial peak matches the analytic value.

Usage
─────
  julia --project=easlee -t 4 test/test_two_puff_mixing.jl
"""

using Random, LinearAlgebra, Printf, Statistics
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
# Physical parameters
# ─────────────────────────────────────────────────────────────────────────────

const K_H      = 1_000.0
const DT       = 60.0
const N_STEPS  = 20
const SIGMA_0  = 5_000.0
const X1       = -7_500.0
const X2       =  7_500.0
const D        = 0.5*(X2-X1)
const T_FINAL  = N_STEPS * DT
const SIG2_T   = SIGMA_0^2 + 2.0*K_H*T_FINAL

const N_TRAJ   = 10_000     # parcels per puff

# KDE bandwidth: Silverman's rule for 2-D
const H_KDE    = sqrt(SIG2_T) * N_TRAJ^(-1/6)   # ≈ 1130 m

# Evaluation grid for concentration estimates
const NX_EVAL  = 60
const LH       = 30_000.0
const DX_EVAL  = 2*LH / NX_EVAL
xs_eval = [(-LH + (i-0.5)*DX_EVAL) for i in 1:NX_EVAL]
ys_eval = [(-LH + (j-0.5)*DX_EVAL) for j in 1:NX_EVAL]
xs_km   = xs_eval ./ 1000
ys_km   = ys_eval ./ 1000

# Old-scheme config (kernel width = H_KDE ≈ neighbouring parcel scale)
const SIGMA_NBR = H_KDE
const MCFG = MixingConfig(sigma0_h=SIGMA_NBR, K_h=K_H, nsig=3.0, max_neighbors=50)

# New-scheme config
const SIGMA_H0_PUFF = SIGMA_NBR
const SIGMA_Z0_PUFF = 1.0
const RHO_BG        = 1.2
const M0_PUFF       = RHO_BG * SIGMA_H0_PUFF^2 * SIGMA_Z0_PUFF
const K_EFF_FACTOR  = sqrt(2.4)
const DSH_STEP      = Sym22(2.0*K_H*DT/K_EFF_FACTOR, 0.0, 2.0*K_H*DT/K_EFF_FACTOR)
const PCFG          = PuffMixingConfig(max_lam_cap=SIGMA_H0_PUFF^2)

# ─────────────────────────────────────────────────────────────────────────────
# Initial parcel positions: draw from each puff's Gaussian
# ─────────────────────────────────────────────────────────────────────────────

rng = Random.MersenneTwister(42)

function sample_puff(cx, cy, n)
    x = cx .+ SIGMA_0 * randn(rng, n)
    y = cy .+ SIGMA_0 * randn(rng, n)
    x, y
end

x1p, y1p = sample_puff(X1, 0.0, N_TRAJ)   # puff-1 parcel positions
x2p, y2p = sample_puff(X2, 0.0, N_TRAJ)   # puff-2 parcel positions

# Combined parcel arrays for the mixing schemes (20k total)
x_all  = vcat(x1p, x2p)
y_all  = vcat(y1p, y2p)
z_all  = zeros(2*N_TRAJ)
act_all = trues(2*N_TRAJ)

# Initial mixing ratios: puff-1 parcels carry species 1; puff-2 carry species 2
q0 = zeros(2*N_TRAJ, 2)
q0[1:N_TRAJ,       1] .= 1.0   # puff-1: tracer 1
q0[N_TRAJ+1:end,   2] .= 1.0   # puff-2: tracer 2

# ─────────────────────────────────────────────────────────────────────────────
# Gaussian KDE: estimate concentration field from weighted parcel positions
# ─────────────────────────────────────────────────────────────────────────────

"""
    kde_field(xs_p, ys_p, weights, h)

Gaussian KDE estimate of the concentration field.

Normalization: ∫∫ ĉ dxdy = 2π·σ₀² · (Σ weights / N_TRAJ)
so that a unit-weight run of N_TRAJ parcels from the initial Gaussian
gives the same spatial integral as the analytic c field.

The un-normalised kernel accumulates Σ_k w_k · exp(−r²/(2h²)).
Final factor: σ₀²/(N_TRAJ·h²) converts to concentration units.
"""
function kde_field(xs_p, ys_p, weights, h)
    inv2h2 = 0.5 / h^2
    norm   = SIGMA_0^2 / (N_TRAJ * h^2)
    c = zeros(NX_EVAL, NX_EVAL)
    for k in eachindex(xs_p)
        w = weights[k]
        w == 0.0 && continue
        for iy in 1:NX_EVAL
            dy2 = (ys_p[k] - ys_eval[iy])^2
            dy2 > 9h^2 && continue
            for ix in 1:NX_EVAL
                dx = xs_p[k] - xs_eval[ix]
                r2 = dx*dx + dy2
                r2 > 9h^2 && continue
                c[ix, iy] += w * exp(-r2 * inv2h2)
            end
        end
    end
    c .* norm
end

# ─────────────────────────────────────────────────────────────────────────────
# Analytic solution
# ─────────────────────────────────────────────────────────────────────────────

function c_analytic_grid(cx, cy)
    [let r2=(xs_eval[ix]-cx)^2+(ys_eval[iy]-cy)^2
         (SIGMA_0^2/SIG2_T)*exp(-r2/(2*SIG2_T))
     end
     for ix in 1:NX_EVAL, iy in 1:NX_EVAL]
end

f_analytic_x(x) = 1.0 / (1.0 + exp(2.0*D*x/SIG2_T))

ac1 = c_analytic_grid(X1, 0.0)
ac2 = c_analytic_grid(X2, 0.0)

# ─────────────────────────────────────────────────────────────────────────────
# Stochastic reference: pure Brownian motion (no mixing scheme)
# ─────────────────────────────────────────────────────────────────────────────
println("Running stochastic reference (pure Brownian motion, $(N_TRAJ) walkers/puff) …")

function brownian_final(xs0, ys0)
    σ_step = sqrt(2.0*K_H*DT)
    x = copy(xs0);  y = copy(ys0)
    for _ in 1:N_STEPS
        x .+= σ_step * randn(rng, length(x))
        y .+= σ_step * randn(rng, length(y))
    end
    x, y
end

xb1, yb1 = brownian_final(x1p, y1p)
xb2, yb2 = brownian_final(x2p, y2p)

# KDE with unit weights (all walkers carry equal tracer)
sc1 = kde_field(xb1, yb1, ones(N_TRAJ), H_KDE)
sc2 = kde_field(xb2, yb2, ones(N_TRAJ), H_KDE)

# ─────────────────────────────────────────────────────────────────────────────
# Old scheme: apply_mixing_xyz! on 20 000 random parcels
# ─────────────────────────────────────────────────────────────────────────────
println("Running old scheme (apply_mixing_xyz! on $(2*N_TRAJ) random parcels) …")

q_old = copy(q0)
let yl = fill(1.0, 2*N_TRAJ)
    for _ in 1:N_STEPS
        apply_mixing_xyz!(q_old, copy(x_all), copy(y_all), copy(z_all),
                          yl, act_all, DT, MCFG; merge_packed=false)
    end
end

# Concentration estimate: KDE weighted by updated mixing ratios
old_c1 = kde_field(x_all, y_all, q_old[:, 1], H_KDE)
old_c2 = kde_field(x_all, y_all, q_old[:, 2], H_KDE)

# ─────────────────────────────────────────────────────────────────────────────
# New scheme: apply_puff_mixing! on the same 20 000 parcels
# ─────────────────────────────────────────────────────────────────────────────
println("Running new scheme (apply_puff_mixing! on $(2*N_TRAJ) random parcels) …")

q_new = copy(q0)
let
    cov  = PuffCovarianceState(2*N_TRAJ; sigma_H0=SIGMA_H0_PUFF,
                                          sigma_z0=SIGMA_Z0_PUFF, m0=M0_PUFF)
    dSH  = fill(DSH_STEP, 2*N_TRAJ)
    dsz2 = fill(0.0, 2*N_TRAJ)
    q_bg = [0.0, 0.0]
    for _ in 1:N_STEPS
        apply_puff_mixing!(cov, q_new, dSH, dsz2, x_all, y_all, z_all,
                           act_all, RHO_BG, q_bg, PCFG)
    end
end

new_c1 = kde_field(x_all, y_all, q_new[:, 1], H_KDE)
new_c2 = kde_field(x_all, y_all, q_new[:, 2], H_KDE)

# ─────────────────────────────────────────────────────────────────────────────
# Mixing fraction  f = c₁/(c₁+c₂),  set to NaN where signal too weak
# ─────────────────────────────────────────────────────────────────────────────

function mixing_frac(c1, c2; thresh=5e-3)
    [let s=c1[ix,iy]+c2[ix,iy]; s>thresh ? c1[ix,iy]/s : NaN end
     for ix in 1:NX_EVAL, iy in 1:NX_EVAL]
end

thresh = 0.005 * maximum(ac1)
af  = mixing_frac(ac1,   ac2;  thresh=thresh)
sf  = mixing_frac(sc1,   sc2;  thresh=thresh)
of  = mixing_frac(old_c1, old_c2; thresh=thresh)
nf  = mixing_frac(new_c1, new_c2; thresh=thresh)

# ─────────────────────────────────────────────────────────────────────────────
# Statistics: total tracer & RMSE of mixing fraction along x-axis (y≈0)
# ─────────────────────────────────────────────────────────────────────────────

iy_mid = NX_EVAL ÷ 2

a_xf   = [f_analytic_x(xs_eval[ix]) for ix in 1:NX_EVAL]
s_xf   = sf[:,  iy_mid]
o_xf   = of[:,  iy_mid]
n_xf   = nf[:,  iy_mid]

mask   = .!isnan.(a_xf) .& .!isnan.(o_xf) .& .!isnan.(n_xf) .& .!isnan.(s_xf)

rmse_s = sqrt(mean((s_xf[mask] .- a_xf[mask]).^2))
rmse_o = sqrt(mean((o_xf[mask] .- a_xf[mask]).^2))
rmse_n = sqrt(mean((n_xf[mask] .- a_xf[mask]).^2))

total_ref = sum(ac1)
@printf("\n")
println("=" ^ 70)
println("  TWO-PUFF MIXING  ($(2*N_TRAJ) random parcels, $(N_STEPS) steps, t=$(Int(T_FINAL)) s)")
println("  σ₀=$(Int(SIGMA_0)) m  x₁=$(Int(X1)) m  x₂=$(Int(X2)) m  K_H=$(Int(K_H)) m²/s")
println("=" ^ 70)
@printf("  %-28s  %9s  %9s  %9s\n", "Metric", "Stoch", "Old", "New")
@printf("  %-28s  %9.4f  %9.4f  %9.4f\n", "Total c₁ / analytic",
        sum(sc1)/total_ref, sum(old_c1)/total_ref, sum(new_c1)/total_ref)
@printf("  %-28s  %9.4f  %9.4f  %9.4f\n", "Total c₂ / analytic",
        sum(sc2)/total_ref, sum(old_c2)/total_ref, sum(new_c2)/total_ref)
@printf("  %-28s  %9.4e  %9.4e  %9.4e\n", "RMSE(f) along x-axis",
        rmse_s, rmse_o, rmse_n)
println("=" ^ 70)
println()
println("  Background mixing:")
println("  Old (mixing_xyz): NO — only peer-to-peer parcel exchange;")
println("                    parcels with no neighbours don't mix at all.")
println("  New (puff_mixing): YES — shell fraction not covered by neighbours")
println("                    is diluted into q_bg = [0,0].")
println()

# ─────────────────────────────────────────────────────────────────────────────
# Figure (optional — numeric checks above are complete)
# ─────────────────────────────────────────────────────────────────────────────
SAVE_PLOT || exit(0)

fig = Figure(size=(1400, 1100), fontsize=13)

Label(fig[0, :],
      "Two-puff mixing  |  $(N_TRAJ) random parcels/puff  |  " *
      "t = $(N_STEPS)×$(Int(DT)) s = $(Int(T_FINAL/60)) min  |  " *
      "K_H=$(Int(K_H)) m²/s  |  σ₀=$(Int(SIGMA_0/1000)) km  |  separation=$(Int((X2-X1)/1000)) km",
      fontsize=14, font=:bold)

titles = ["Analytic", "Stochastic (Brownian)", "Old: mixing_xyz", "New: puff_mixing"]

# ── Row 1: c₁ field ──────────────────────────────────────────────────────────
clim_c = (0.0, maximum(ac1) * 1.05)
for (j, (title, field)) in enumerate(zip(titles,
        [ac1, sc1, old_c1, new_c1]))
    ax = Axis(fig[1, j], title=title, xlabel="x [km]",
              ylabel=(j==1 ? "y [km]" : ""), yticklabelsvisible=(j==1))
    heatmap!(ax, xs_km, ys_km, field; colormap=:thermal, colorrange=clim_c)
end
Colorbar(fig[1, 5]; colormap=:thermal, limits=clim_c, label="c₁  (puff-1 tracer)")
for j in 1:4; colsize!(fig.layout, j, Relative(0.22)); end

# ── Row 2: mixing fraction ────────────────────────────────────────────────────
for (j, (title, field)) in enumerate(zip(titles,
        [af, sf, of, nf]))
    ax = Axis(fig[2, j], title=(j==1 ? "Mixing fraction  f = c₁/(c₁+c₂)" : ""),
              xlabel="x [km]", ylabel=(j==1 ? "y [km]" : ""),
              yticklabelsvisible=(j==1))
    heatmap!(ax, xs_km, ys_km, field; colormap=:RdBu, colorrange=(0.0, 1.0),
             nan_color=(:grey70, 0.4))
    scatter!(ax, [X1/1000, X2/1000], [0.0, 0.0]; color=:yellow, markersize=8,
             strokewidth=1, strokecolor=:black)
end
Colorbar(fig[2, 5]; colormap=:RdBu, limits=(0.0, 1.0), label="f = c₁/(c₁+c₂)")

# ── Row 3: x-axis cross-sections ─────────────────────────────────────────────
ax_c = Axis(fig[3, 1:2],
            title="x cross-section  c₁ and c₂  (y ≈ 0 km)  — KDE estimate",
            xlabel="x [km]", ylabel="concentration")
lines!(ax_c, xs_km, ac1[:, iy_mid]; color=:navy,     linewidth=2.5, linestyle=:dash, label="c₁ analytic")
lines!(ax_c, xs_km, ac2[:, iy_mid]; color=:firebrick,linewidth=2.5, linestyle=:dash, label="c₂ analytic")
lines!(ax_c, xs_km, sc1[:, iy_mid]; color=(:dodgerblue,0.8), linewidth=1.5, linestyle=:dot, label="c₁ stochastic")
lines!(ax_c, xs_km, sc2[:, iy_mid]; color=(:salmon,0.8),     linewidth=1.5, linestyle=:dot, label="c₂ stochastic")
lines!(ax_c, xs_km, old_c1[:, iy_mid]; color=:steelblue, linewidth=2, label="c₁ old")
lines!(ax_c, xs_km, old_c2[:, iy_mid]; color=:tomato,    linewidth=2, label="c₂ old")
lines!(ax_c, xs_km, new_c1[:, iy_mid]; color=:royalblue, linewidth=2, linestyle=:dashdot, label="c₁ new")
lines!(ax_c, xs_km, new_c2[:, iy_mid]; color=:orangered, linewidth=2, linestyle=:dashdot, label="c₂ new")
axislegend(ax_c; position=:ct, nbanks=2, labelsize=10)

ax_f = Axis(fig[3, 3:5],
            title="x cross-section  f = c₁/(c₁+c₂)  (y ≈ 0 km)  — logistic analytic",
            xlabel="x [km]", ylabel="f  (mixing fraction)")
lines!(ax_f, xs_km, a_xf;  color=:black,     linewidth=3, linestyle=:dash, label="Analytic (logistic)")
lines!(ax_f, xs_km, replace(s_xf, NaN=>NaN); color=(:grey50,0.8), linewidth=1.5, linestyle=:dot, label="Stochastic")
lines!(ax_f, xs_km, replace(o_xf, NaN=>NaN); color=:steelblue, linewidth=2.5, label="Old: mixing_xyz")
lines!(ax_f, xs_km, replace(n_xf, NaN=>NaN); color=:crimson,   linewidth=2.5, label="New: puff_mixing")
hlines!(ax_f, [0.5]; color=(:black, 0.25), linewidth=1)
axislegend(ax_f; position=:rt)
text!(ax_f, -27.0, 0.90;
      text="RMSE(f) on x-axis:\n" *
           "  Stochastic: $(round(rmse_s, sigdigits=2))\n" *
           "  Old:         $(round(rmse_o, sigdigits=2))\n" *
           "  New:         $(round(rmse_n, sigdigits=2))",
      fontsize=11, align=(:left, :top))

outpath = joinpath(REPO, "results", "two_puff_mixing.png")
mkpath(dirname(outpath))
save(outpath, fig; px_per_unit=2)
println("Saved → $outpath")
