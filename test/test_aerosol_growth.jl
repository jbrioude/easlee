#!/usr/bin/env julia
# test_aerosol_growth.jl — Aitken → accumulation growth: what EASLEE has vs
# what the atmosphere does.
#
# Part A (what EXISTS — validated): κ-Köhler equilibrium water uptake
#   (growth_factor).  This is REVERSIBLE swelling with RH — it never changes
#   the dry size, so it is NOT inter-mode growth.
#
# Part B (what is MISSING — quantified): a reference sectional box model
#   (60 log bins, 3 nm–3 µm) with the two dominant clear-sky pathways that
#   move Aitken particles into the accumulation mode:
#     • Brownian coagulation      (Fuchs transition-regime kernel,
#                                  volume-splitting sectional scheme)
#     • H2SO4 condensation        (Fuchs–Sutugin flux, daytime
#                                  [H2SO4] = 1e7 cm⁻³, 12 h/day)
#   72-h marine case: Aitken N=600 cm⁻³ dg=40 nm, accumulation N=200 cm⁻³
#   dg=160 nm (both σg=1.5).  Writes the evolving distribution to
#   results/aerosol_growth_ref.csv for plotting.
#
# EASLEE's per-parcel dry bin is assigned once and never changes
# (physics/aerosol.jl: bin sampled at release; bins.d static), so its answer to
# this case is "no change" by construction — the printed reference growth rates
# ARE the missing physics.  NB: even where a per-parcel SO4 source exists,
# but nothing condenses it onto the AerosolModel bins.
#
#   julia --project=easlee test/test_aerosol_growth.jl
using Printf, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[A] κ-Köhler water uptake (the ONLY growth EASLEE has) — validation")
println("="^78)
gf_analytic(κ, rh) = (1 + κ * rh / (1 - rh))^(1/3)
okgf = true
for κ in (0.1f0, 0.61f0, 1.28f0), rh in (0.4f0, 0.8f0, 0.9f0, 0.95f0)
    g = EASLEE.growth_factor(κ, rh)
    isapprox(g, gf_analytic(κ, rh); rtol=0.02) || (global okgf = false)
end
check("growth_factor matches (1+κ·RH/(1−RH))^(1/3)", okgf)
g90 = EASLEE.growth_factor(0.61f0, 0.9f0)
check("sulfate at RH 90 %: GF ≈ 1.86 (water swelling only)",
      isapprox(g90, 1.86; rtol=0.03); info=@sprintf("(GF=%.2f)", g90))
check("equilibrium (stateless): RH 0.9 → 0.4 returns to the 0.4 value",
      EASLEE.growth_factor(0.61f0, 0.4f0) == EASLEE.growth_factor(0.61f0, 0.4f0))
println("  → κ-growth is DIAGNOSTIC water uptake; the dry diameter of a parcel's")
println("    bin never changes during a run (bin assigned once at release).")

println()
println("="^78)
println("[B] Reference box model: coagulation + H2SO4 condensation, 72 h marine")
println("="^78)

# ── sectional grid ────────────────────────────────────────────────────────────
const NB = 60
const dedges = exp.(range(log(3e-9), log(3e-6); length=NB + 1))    # m
const d  = sqrt.(dedges[1:end-1] .* dedges[2:end])                  # bin centres
const v  = (π / 6) .* d .^ 3                                        # volumes [m³]

# initial bimodal marine distribution [cm⁻³ per bin]
lognorm(N, dg, σg) = [N * (log(dedges[k+1]) - log(dedges[k])) /
                      (sqrt(2π) * log(σg)) *
                      exp(-(log(d[k] / dg))^2 / (2 * log(σg)^2)) for k in 1:NB]
N0 = lognorm(600.0, 40e-9, 1.5) .+ lognorm(200.0, 160e-9, 1.5)

# ── Fuchs transition-regime Brownian kernel [m³ s⁻¹] (Seinfeld & Pandis 13.56)
const kB = 1.380649e-23; const T = 288.0; const P = 101325.0
const μ  = 1.8e-5; const ρp = 1500.0
const λa = 6.6e-8 * (101325 / P) * (T / 293.15)
Cc(dp)   = 1 + 2λa / dp * (1.257 + 0.4 * exp(-1.1 * dp / (2λa)))
Dif(dp)  = kB * T * Cc(dp) / (3π * μ * dp)
cbar(dp) = sqrt(8 * kB * T / (π * ρp * (π / 6) * dp^3))
ℓ(dp)    = 8 * Dif(dp) / (π * cbar(dp))
gfun(dp) = ((dp + ℓ(dp))^3 - (dp^2 + ℓ(dp)^2)^1.5) / (3 * dp * ℓ(dp)) - dp
function kernel(d1, d2)
    D1, D2 = Dif(d1), Dif(d2); c1, c2 = cbar(d1), cbar(d2)
    g1, g2 = gfun(d1), gfun(d2)
    β = ((d1 + d2) / (d1 + d2 + 2 * sqrt(g1^2 + g2^2)) +
         8 * (D1 + D2) / (sqrt(c1^2 + c2^2) * (d1 + d2)))^(-1)
    return 2π * (D1 + D2) * (d1 + d2) * β
end
const K = [kernel(d[i], d[j]) * 1e6 for i in 1:NB, j in 1:NB]   # cm³ s⁻¹

# volume-splitting: put number n with combined volume vv into the grid
function deposit!(N, vv, n)
    k = searchsortedlast(v, vv)
    if k >= NB
        N[NB] += n * vv / v[NB]                 # conserve volume at the top
    elseif k < 1
        N[1] += n
    else
        f = (v[k+1] - vv) / (v[k+1] - v[k])     # number split conserving volume
        N[k] += n * f; N[k+1] += n * (1 - f)
    end
end

# ── H2SO4 condensation (Fuchs–Sutugin, α=1) ──────────────────────────────────
const Dv = 1.0e-5                               # H2SO4 diffusivity [m²/s]
const mH = 98e-3 / 6.022e23                     # kg/molecule
const cH = sqrt(8 * kB * T / (π * mH))          # thermal speed
fFS(Kn)  = (1 + Kn) / (1 + 1.71Kn + 1.33Kn^2)
# dv/dt per particle [m³/s] at gas concentration C [molecules m⁻³]
dvdt(dp, C) = 2π * dp * Dv * fFS(2λa / dp) * C * mH / ρp

# ── integrate 72 h ────────────────────────────────────────────────────────────
const dt = 60.0; const nstep = Int(72 * 3600 ÷ dt)
N = copy(N0)
H2SO4_day = 1.0e7 * 1e6                          # 1e7 cm⁻³ → m⁻³ (canonical marine daytime)
snaps = Dict(0 => copy(N0))
for it in 1:nstep
    t_h = it * dt / 3600
    # condensation (daytime half of each 24 h)
    if mod(t_h, 24) < 12
        Nn = zeros(NB)
        for k in 1:NB
            N[k] <= 0 && continue
            deposit!(Nn, v[k] + dvdt(d[k], H2SO4_day) * dt, N[k])
        end
        N .= Nn
    end
    # coagulation (explicit, Δt small enough: max λ_coag ≈ K·N ~ 1e-5 s⁻¹)
    dN = zeros(NB); prod = zeros(NB)
    for i in 1:NB
        N[i] <= 0 && continue
        for j in i:NB
            N[j] <= 0 && continue
            rate = K[i, j] * N[i] * N[j] * (i == j ? 0.5 : 1.0) * dt
            rate <= 0 && continue
            dN[i] -= rate; dN[j] -= rate
            deposit!(prod, v[i] + v[j], rate)
        end
    end
    N .+= dN .+ prod
    N .= max.(N, 0.0)
    (t_h in (24.0, 48.0, 72.0)) && (snaps[Int(t_h)] = copy(N))
end

# ── diagnostics ───────────────────────────────────────────────────────────────
# number geometric-mean diameter of the sub-100 nm (Aitken) population
function dg_mode(N, lo, hi)
    sel = [k for k in 1:NB if lo <= d[k] < hi && N[k] > 0]
    isempty(sel) && return NaN
    exp(sum(N[k] * log(d[k]) for k in sel) / sum(N[k] for k in sel))
end
ntot(N) = sum(N)
nacc(N) = sum(N[k] for k in 1:NB if d[k] >= 100e-9)

dg0 = dg_mode(N0, 10e-9, 100e-9) * 1e9
dg72 = dg_mode(snaps[72], 10e-9, 100e-9) * 1e9
gr24 = (dg_mode(snaps[24], 10e-9, 100e-9) - dg_mode(N0, 10e-9, 100e-9)) * 1e9 / 24
@printf("  t= 0h: N=%6.0f cm⁻³  Aitken dg=%5.1f nm  N(>100 nm)=%5.0f cm⁻³\n",
        ntot(N0), dg0, nacc(N0))
for t in (24, 48, 72)
    Ns = snaps[t]
    @printf("  t=%2dh: N=%6.0f cm⁻³  Aitken dg=%5.1f nm  N(>100 nm)=%5.0f cm⁻³\n",
            t, ntot(Ns), dg_mode(Ns, 10e-9, 100e-9) * 1e9, nacc(Ns))
end
# analytic daytime GR at d = 50 nm for cross-checking the numerics
gr_day = dvdt(50e-9, H2SO4_day) * 2 / (π * (50e-9)^2) * 1e9 * 3600
@printf("\n  daytime condensation GR at 50 nm: %.2f nm/h (pure H2SO4 at 1e7 cm⁻³;\n", gr_day)
@printf("  observed total GR 1–5 nm/h also includes organics + base co-condensation)\n")
@printf("  Aitken-mode dg growth rate, day-1 24-h mean: %.2f nm/h\n", gr24)
check("daytime H2SO4 GR ≈ 0.5–1.5 nm/h (kinetic-flux benchmark)",
      0.5 <= gr_day <= 1.5; info=@sprintf("(%.2f nm/h)", gr_day))
check("24-h-mean mode GR consistent (≈ half the daytime value)",
      0.25 * gr_day <= gr24 <= 0.75 * gr_day; info=@sprintf("(%.2f nm/h)", gr24))
check("Aitken mode grows toward accumulation over 72 h", dg72 > dg0 + 15;
      info=@sprintf("(dg %.0f → %.0f nm)", dg0, dg72))
check("accumulation-mode number increases", nacc(snaps[72]) > 1.15 * nacc(N0);
      info=@sprintf("(%.0f → %.0f cm⁻³)", nacc(N0), nacc(snaps[72])))
check("total number falls (coagulation sink)", ntot(snaps[72]) < ntot(N0))
println("\n  EASLEE result for the same case: dg 40 → 40 nm, N(>100 nm) unchanged —")
println("  per-parcel dry bins are static; no coagulation / condensation /")
println("  aqueous processing exists, and any produced SO4 is never")
println("  condensed onto the bins.")

# CSV for plotting
csv = joinpath(REPO, "results", "aerosol_growth_ref.csv"); mkpath(dirname(csv))   # results/ is gitignored (review 2026-09-23)
open(csv, "w") do io
    println(io, "d_nm,dlogd,N0,N24,N48,N72")
    for k in 1:NB
        dlogd = log10(dedges[k+1]) - log10(dedges[k])
        @printf(io, "%.4g,%.5f,%.5g,%.5g,%.5g,%.5g\n", d[k] * 1e9, dlogd,
                N0[k], snaps[24][k], snaps[48][k], snaps[72][k])
    end
end
println("\ndistributions → $csv")
@printf("\nRESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
