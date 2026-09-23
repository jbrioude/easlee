"""
plot_plume_comparison.jl
═════════════════════════
Side-by-side visual comparison of plume shape produced by the two mixing
algorithms after N_STEPS time steps.

Panels
──────
  (row 1) 2-D concentration heatmaps — old | new | new−old
  (row 2) Radial profiles c(r) — old, new, analytic theory
  (row 3) Normalised residual |c - c_theory| / c_theory along x-axis

Output
──────
  results/plume_shape_comparison.png

Usage
─────
  julia -t 4 scripts/plot_plume_comparison.jl
"""

using Random, LinearAlgebra, Printf, Statistics
using GLMakie   # not an EASLEE dep — install once in your default env: Pkg.add("GLMakie")

const REPO = joinpath(@__DIR__, "..", "..")
include(joinpath(REPO, "easlee", "EASLEE.jl"))
using .EASLEE: MixingConfig, apply_mixing_xyz!,
               Sym22, PuffCovarianceState, PuffMixingConfig, apply_puff_mixing!

Random.seed!(1234)
mkpath(joinpath(REPO, "results"))

# ─────────────────────────────────────────────────────────────────────────────
# Physical parameters  (identical to compare_mixing.jl)
# ─────────────────────────────────────────────────────────────────────────────

const K_H        = 1_000.0
const DT         = 60.0
const N_STEPS    = 20
const NX         = 100
const LH         = 50_000.0
const DX         = 2.0 * LH / NX
const N          = NX * NX
const SIGMA_SRC  = 5_000.0
const Z0         = 1_000.0

const SIGMA0_H   = DX
const MCFG = MixingConfig(sigma0_h=SIGMA0_H, K_h=K_H, nsig=3.0, max_neighbors=50)

const SIGMA_H0_PUFF = DX
const SIGMA_Z0_PUFF = 1.0
const K_EFF_RATIO   = 2.4
const K_EFF_FACTOR  = sqrt(K_EFF_RATIO)
const DSH_STEP   = Sym22(2.0*K_H*DT/K_EFF_FACTOR, 0.0, 2.0*K_H*DT/K_EFF_FACTOR)
const DSIGMA2_Z  = 0.0
const RHO_BG     = 1.2
const M0_PUFF    = RHO_BG * SIGMA_H0_PUFF^2 * SIGMA_Z0_PUFF
const PCFG = PuffMixingConfig(max_lam_cap=SIGMA_H0_PUFF^2)

# Grid positions: parcel k = (jj-1)*NX + ii  →  x=(ii-0.5)*DX-LH, y=(jj-0.5)*DX-LH
x_grid = [(((k-1) % NX) + 0.5) * DX - LH for k in 1:N]
y_grid = [(((k-1) ÷ NX) + 0.5) * DX - LH for k in 1:N]
z_grid = fill(Z0, N)

# Axis tick values in km for heatmaps
const xs_km = range(-LH/1000, LH/1000, length=NX)
const ys_km = range(-LH/1000, LH/1000, length=NX)

# ─────────────────────────────────────────────────────────────────────────────
# Initial tracer: 2-D Gaussian c(x,y) = exp(−r²/(2σ_src²))
# ─────────────────────────────────────────────────────────────────────────────

init_conc() = [exp(-(x_grid[i]^2 + y_grid[i]^2) / (2.0 * SIGMA_SRC^2)) for i in 1:N]

# Analytic solution after time t
function c_theory_2d(x, y, t)
    sig2 = SIGMA_SRC^2 + 2.0 * K_H * t
    (SIGMA_SRC^2 / sig2) * exp(-(x^2 + y^2) / (2.0 * sig2))
end

# ─────────────────────────────────────────────────────────────────────────────
# Warm-up (compile everything before timing)
# ─────────────────────────────────────────────────────────────────────────────
println("Compiling …")

let
    c   = init_conc()
    q   = reshape(copy(c), N, 1)
    yl  = fill(1.0, N)
    act = trues(N)
    apply_mixing_xyz!(q, copy(x_grid), copy(y_grid), copy(z_grid), yl, act, DT, MCFG;
                      merge_packed=false)
end
let
    q   = reshape(init_conc(), N, 1)
    cov = PuffCovarianceState(N; sigma_H0=SIGMA_H0_PUFF, sigma_z0=SIGMA_Z0_PUFF, m0=M0_PUFF)
    act = trues(N)
    dSH = fill(DSH_STEP, N);  dsz2 = fill(DSIGMA2_Z, N);  q_bg = [0.0]
    apply_puff_mixing!(cov, q, dSH, dsz2, x_grid, y_grid, z_grid, act, RHO_BG, q_bg, PCFG)
end
GC.gc()
println("Ready.\n")

# ─────────────────────────────────────────────────────────────────────────────
# Run old algorithm
# ─────────────────────────────────────────────────────────────────────────────
println("Running old algorithm (apply_mixing_xyz!) …")
q_old = let
    q   = reshape(init_conc(), N, 1)
    yl  = fill(1.0, N)
    act = trues(N)
    for _ in 1:N_STEPS
        apply_mixing_xyz!(q, copy(x_grid), copy(y_grid), copy(z_grid), yl, act, DT, MCFG;
                          merge_packed=false)
    end
    q[:, 1]
end

# ─────────────────────────────────────────────────────────────────────────────
# Run new algorithm
# ─────────────────────────────────────────────────────────────────────────────
println("Running new algorithm (apply_puff_mixing!) …")
q_new = let
    q   = reshape(init_conc(), N, 1)
    cov = PuffCovarianceState(N; sigma_H0=SIGMA_H0_PUFF, sigma_z0=SIGMA_Z0_PUFF, m0=M0_PUFF)
    act = trues(N)
    dSH  = fill(DSH_STEP, N)
    dsz2 = fill(DSIGMA2_Z, N)
    q_bg = [0.0]
    for _ in 1:N_STEPS
        apply_puff_mixing!(cov, q, dSH, dsz2, x_grid, y_grid, z_grid,
                           act, RHO_BG, q_bg, PCFG)
    end
    q[:, 1]
end

# Analytic field at final time
T_FINAL = N_STEPS * DT
q_theory = [c_theory_2d(x_grid[i], y_grid[i], T_FINAL) for i in 1:N]

# ─────────────────────────────────────────────────────────────────────────────
# Reshape to 2-D grids for heatmaps
# Grid layout: parcel k = (jj-1)*NX + ii  →  grid[ii, jj]
# ─────────────────────────────────────────────────────────────────────────────
field_old    = reshape(q_old,    NX, NX)
field_new    = reshape(q_new,    NX, NX)
field_theory = reshape(q_theory, NX, NX)

# Normalize to peak=1 for shape-only comparison
norm_old    = field_old    ./ maximum(field_old)
norm_new    = field_new    ./ maximum(field_new)
norm_theory = field_theory ./ maximum(field_theory)

diff_norm = norm_new .- norm_theory  # shape residual (amplitude-free)

# ─────────────────────────────────────────────────────────────────────────────
# Radial profiles  — bin starting from r=DX/2 to avoid grid aliasing at r≈0
# ─────────────────────────────────────────────────────────────────────────────
r = sqrt.(x_grid.^2 .+ y_grid.^2)

N_BINS   = 50
r_min    = DX * 0.5
r_max    = 3.0 * sqrt(SIGMA_SRC^2 + 2.0 * K_H * T_FINAL)
r_edges  = range(r_min, r_max, length=N_BINS+1)
r_mid    = [(r_edges[k] + r_edges[k+1]) / 2 for k in 1:N_BINS]
r_mid_km = r_mid ./ 1000

function radial_profile(q_vec)
    counts = zeros(N_BINS)
    sums   = zeros(N_BINS)
    for i in 1:N
        ri = r[i]
        (ri < r_min || ri >= r_max) && continue
        b = searchsortedlast(r_edges, ri)
        b = clamp(b, 1, N_BINS)
        sums[b]   += q_vec[i]
        counts[b] += 1
    end
    [counts[b] > 0 ? sums[b] / counts[b] : NaN for b in 1:N_BINS]
end

rp_old    = radial_profile(q_old)
rp_new    = radial_profile(q_new)
rp_theory = [c_theory_2d(r_mid[b], 0.0, T_FINAL) for b in 1:N_BINS]

# Normalised radial profiles (peak=1) — isolates shape from amplitude
rp_old_n    = rp_old    ./ maximum(filter(!isnan, rp_old))
rp_new_n    = rp_new    ./ maximum(filter(!isnan, rp_new))
rp_theory_n = rp_theory ./ maximum(rp_theory)

# x-axis cut (iy = NX÷2 row, all ix)
iy_mid  = NX ÷ 2
xcut_km = collect(xs_km)
xcut_old    = field_old[:,    iy_mid]
xcut_new    = field_new[:,    iy_mid]
xcut_theory = field_theory[:, iy_mid]

# ─────────────────────────────────────────────────────────────────────────────
# Total-mass comparison
# ─────────────────────────────────────────────────────────────────────────────
total_old    = sum(q_old)
total_new    = sum(q_new)
total_theory = sum(q_theory)

# ─────────────────────────────────────────────────────────────────────────────
# Figure  (4 rows)
# ─────────────────────────────────────────────────────────────────────────────
fig = Figure(size=(1350, 1200), fontsize=13)

T_MIN = N_STEPS * DT / 60
Label(fig[0, :],
      "Plume shape after $(N_STEPS) steps × $(Int(DT)) s = $(Int(T_MIN)) min  " *
      "(K_h = $(Int(K_H)) m²/s,  N = $(N) parcels,  $(Threads.nthreads()) threads)",
      fontsize=15, font=:bold)

# ── Row 1: absolute 2-D heatmaps (shared colour scale = theory range) ────────
clims = (0.0, maximum(q_theory) * 1.02)
cmap  = :thermal

ax1 = Axis(fig[1,1], title="Old: apply_mixing_xyz!",
           xlabel="x [km]", ylabel="y [km]", aspect=1)
ax2 = Axis(fig[1,2], title="New: apply_puff_mixing!",
           xlabel="x [km]", ylabel="y [km]", aspect=1)
ax3 = Axis(fig[1,3], title="Theory  (analytic)",
           xlabel="x [km]", ylabel="y [km]", aspect=1)

hm1 = heatmap!(ax1, xs_km, ys_km, field_old;    colormap=cmap, colorrange=clims)
       heatmap!(ax2, xs_km, ys_km, field_new;    colormap=cmap, colorrange=clims)
       heatmap!(ax3, xs_km, ys_km, field_theory; colormap=cmap, colorrange=clims)
Colorbar(fig[1, 4], hm1, label="c  (mixing ratio)")

# ── Row 2: peak-normalised 2-D heatmaps + shape residual ────────────────────
nclims = (0.0, 1.02)

ax4 = Axis(fig[2,1], title="Old  (normalised to peak)",
           xlabel="x [km]", ylabel="y [km]", aspect=1)
ax5 = Axis(fig[2,2], title="New  (normalised to peak)",
           xlabel="x [km]", ylabel="y [km]", aspect=1)
ax6 = Axis(fig[2,3], title="Shape residual: norm(New) − norm(Theory)",
           xlabel="x [km]", ylabel="y [km]", aspect=1)

hm4 = heatmap!(ax4, xs_km, ys_km, norm_old;  colormap=cmap, colorrange=nclims)
       heatmap!(ax5, xs_km, ys_km, norm_new;  colormap=cmap, colorrange=nclims)
dmax = max(maximum(abs.(diff_norm)), 1e-6) * 1.05
hm6  = heatmap!(ax6, xs_km, ys_km, diff_norm; colormap=:RdBu, colorrange=(-dmax, dmax))
Colorbar(fig[2, 4], hm4, label="c / peak(c)")
Colorbar(fig[2, 5], hm6, label="Δ(normalised c)")

# ── Row 3: radial profiles (absolute + normalised) ───────────────────────────
ax7 = Axis(fig[3, 1:2],
           title="Radial profile  c̄(r)  — absolute",
           xlabel="r [km]", ylabel="c  (mixing ratio)")
lines!(ax7, r_mid_km, rp_theory; color=:black,    linewidth=2.5, linestyle=:dash, label="Theory")
lines!(ax7, r_mid_km, rp_old;    color=:steelblue, linewidth=2,                   label="Old")
lines!(ax7, r_mid_km, rp_new;    color=:crimson,   linewidth=2,                   label="New")
axislegend(ax7, position=:rt)

ax8 = Axis(fig[3, 3],
           title="Radial profile  — normalised to peak",
           xlabel="r [km]", ylabel="c / peak(c)")
lines!(ax8, r_mid_km, rp_theory_n; color=:black,    linewidth=2.5, linestyle=:dash, label="Theory")
lines!(ax8, r_mid_km, rp_old_n;    color=:steelblue, linewidth=2,                   label="Old")
lines!(ax8, r_mid_km, rp_new_n;    color=:crimson,   linewidth=2,                   label="New")
axislegend(ax8, position=:rt)

# ── Row 4: x-axis cross-section (absolute + normalised) ──────────────────────
ax9 = Axis(fig[4, 1:2],
           title="x-axis cross-section  (y ≈ 0 km)  — absolute",
           xlabel="x [km]", ylabel="c  (mixing ratio)")
lines!(ax9, xcut_km, xcut_theory; color=:black,    linewidth=2.5, linestyle=:dash, label="Theory")
lines!(ax9, xcut_km, xcut_old;    color=:steelblue, linewidth=2,                   label="Old")
lines!(ax9, xcut_km, xcut_new;    color=:crimson,   linewidth=2,                   label="New")
axislegend(ax9, position=:rt)

ax10 = Axis(fig[4, 3],
            title="x-axis cross-section  — normalised to peak",
            xlabel="x [km]", ylabel="c / peak(c)")
lines!(ax10, xcut_km, xcut_theory ./ maximum(xcut_theory);
       color=:black,    linewidth=2.5, linestyle=:dash, label="Theory")
lines!(ax10, xcut_km, xcut_old    ./ maximum(xcut_old);
       color=:steelblue, linewidth=2,                   label="Old")
lines!(ax10, xcut_km, xcut_new    ./ maximum(xcut_new);
       color=:crimson,   linewidth=2,                   label="New")
axislegend(ax10, position=:rt)

# ─────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────────
rmse_old = sqrt(mean((q_old .- q_theory).^2))
rmse_new = sqrt(mean((q_new .- q_theory).^2))
bias_old = mean(q_old .- q_theory)
bias_new = mean(q_new .- q_theory)
peak_old = maximum(q_old)
peak_new = maximum(q_new)
peak_ref = maximum(q_theory)

rmse_shape_old = sqrt(mean((norm_old .- norm_theory).^2))
rmse_shape_new = sqrt(mean((norm_new .- norm_theory).^2))

println("=" ^ 65)
println("  PLUME SHAPE STATISTICS  (t = $(Int(T_FINAL)) s = $(Int(T_MIN)) min)")
println("=" ^ 65)
@printf("  %-28s  %9s  %9s  %9s\n", "Metric",        "Old",    "New",    "Theory")
@printf("  %-28s  %9.4f  %9.4f  %9.4f\n", "Peak c",
        peak_old, peak_new, peak_ref)
@printf("  %-28s  %9.4f  %9.4f  %9.4f\n", "Total Σc  (mass proxy)",
        total_old, total_new, total_theory)
@printf("  %-28s  %9.4f  %9.4f\n",  "Total / theory",
        total_old/total_theory, total_new/total_theory)
@printf("  %-28s  %9.2e  %9.2e\n",  "RMSE vs theory (absolute)",
        rmse_old, rmse_new)
@printf("  %-28s  %9.2e  %9.2e\n",  "Bias vs theory (absolute)",
        bias_old, bias_new)
@printf("  %-28s  %9.2e  %9.2e\n",  "RMSE vs theory (shape only)",
        rmse_shape_old, rmse_shape_new)
println("=" ^ 65)
println("  Note: 'shape only' metrics use fields normalised to peak=1.")
println("  New algo exchanges mass with background (q_bg=0) → total drops.")

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
outpath = joinpath(REPO, "results", "plume_shape_comparison.png")
save(outpath, fig; px_per_unit=2)
println("\nSaved → $outpath")
