# Model justifications & references (paper support)

Purpose: every non-obvious modelling choice made in EASLEE's volcanic-plume /
estimator work, with (a) the literature basis, (b) the internal
validation evidence (test or run), and (c) an honest flag where a value is a
placeholder or an internal construction rather than a published result.
Companion docs: `CONFIG.md` (volcanic release options), the estimator
sections in `ETEX_RESULTS.md`.

Status of citations: entries marked **[verified 2026-07-07]** were checked
against the source abstract online; the rest are standard results cited from
memory — **verify page/volume details before submission**.

---

## 1. Volcanic eruption column (1-D moist plume, `plume_rise.scheme="volcanic"`)

| choice | basis |
|---|---|
| Entrainment closure u_e = α·w, α ≈ 0.09–0.12 | Morton, Taylor & Turner (1956), Proc. R. Soc. A 234 — the foundational entrainment hypothesis. |
| Column thermodynamics: mixture of pyroclasts + vapour + entrained moist air, latent-heat release by saturation adjustment | Woods (1988), Bull. Volcanol. 50 (dry dynamics); Woods (1993), JGR 98 (moist convection — latent heat can dominate loft for water-rich plumes). |
| Phreatomagmatic case: high vent water fraction (HT: `water_frac = 0.45`, T₀ = 700 K after magma–water interaction) drives the loft; near-dry column collapses | Koyaguchi & Woods (1996), JGR 101 (external-water effects on column dynamics); the HT-specific submarine-eruption water source is discussed in Van Eaton et al. (2023), GRL 50 **[verified]**. Internal evidence: dry column collapses at 4 km vs 40 %-vapour lofting to 29 km (`validation/test_volcanic_column.jl`, 7/7). |
| MER ↔ height cross-check | Mastin et al. (2009), JVGR 186 — empirical H–MER. Internal: untuned column reproduces Mastin ×0.83–1.04. |
| HT source parameters MER 5×10⁹ kg/s, T₀ 700 K, 45 % vapour | Chosen so the computed injection layer matches the observed cloud: model [19.9, 34.0] km vs observed [20, 33] km (stereo); plume top 57–58 km reported by Proud et al. (2022), Science 378 **[verified via Van Eaton search]** is the transient overshoot, not the persistent layer. **Internal calibration choice** — document as such. |
| Suzuki vertical profile for satellite-table releases | Suzuki (1983), in *Arc Volcanism*; standard tephra-model profile. Coefficients vs the paper verified for the Sofiev BB scheme separately (Sofiev et al. 2012, ACP 12). |

## 2. Umbrella intrusion (gravity-current spreading, `umbrella_spread=true`)

| choice | basis |
|---|---|
| Front law R(t) = (3λNQ_V t²/2π)^(1/3), Ṙ = (2/3)R/t | Woods & Kienle (1994), JVGR 62; Sparks et al. (1997), *Volcanic Plumes*, ch. 11; identical closure used in ash models (Costa et al. 2013, GRL 40 — density-driven umbrella transport). Derivation: volume conservation V = Q_V t = πR²h̄ + intrusion front condition Ṙ = λNh. |
| Intrusion constant λ = 0.2 | Costa et al. (2013), GRL 40 (and refs therein); standard range 0.1–0.3. |
| Q_V = MER·amp/ρ(NBH), N from the met sounding | Internal: both from the SAME validated 1-D column (no free parameter beyond λ). Zero-tuning result: R(1 h) = 173 km vs observed ≥200 km radius at 1 h; Van Eaton et al. (2023), GRL 50 **[verified]**: umbrella >400 km across within 1 h, ~600 km diameter within 2 h — our digitised stereo envelope (equiv. radius 242 km at +1 h, 276 km at +2 h) is consistent. |
| **Interior radial flow u(r) = Q_V/(2πrh̄) = R²/(2rt)** (source-fed), saturated inside the collapse dome r < r_src, with the front-speed floor Ṙ·min(r/R,1) | Steady radial continuity — the same advective umbrella wind used by tephra models (Costa et al. 2013). Internal finding that motivated it (2026-07-07 0.25° verdict run): the previous self-similar u ∝ r conserves each parcel's front fraction exactly (df/dt = 0), so vent-released material never reaches the head — parcel footprint tracked only ~½R. With the source-fed term the parcel mass median tracks ~0.9R (measured). |
| Collapse-dome radius r_src = 40 km (HT) | Satellite-imaged central dome scale of the HT plume (Van Eaton et al. 2023 imagery; Proud et al. 2022). Config knob `umbrella_r_source`; generic default 10 km. **Scale choice, not a fitted constant** — u(r_src, +10 min) ≈ 40–60 m/s follows, matching the observed early front speeds (35–50 m/s class). |
| **Feeding cutoff + coasting**: constant-flux law only while the column feeds (t ≤ t_feed, mass-weighted release-end age); afterwards constant-volume coasting R ∝ t^(1/3) with the same front condition (V frozen: dR/dt = λNV/(πR²) ⇒ R³ ∝ t) | Gravity-current regime change, e.g. Ungarish (2009), *An Introduction to Gravity Currents and Intrusions* (constant-flux vs constant-volume regimes). Internal evidence: the observed stereo radius follows the transition quantitatively — 272 km at +1.5 h × (4.5/1.5)^(1/3) = 392 km vs the observed 389 km plateau at +4.5 h.  Verdict run (v3, 2026-07-07, HT t_feed = 1.93 h): analytic R(6 h) = 391 km vs observed 380–389; parcel annulus at +2 h spans 150–400 km (obs equivalent radius 276–309), mass tracking ~0.9R during feeding.  **Metric note for the paper**: percent-of-max footprint radii under-read an annulus, and quantiles about the mass centroid are invalid once the ring drifts (the centroid sits on the dense arc) — report maps + the analytic law + quantiles about the ring centre. |
| Momentum bookkeeping: the radial push is imposed (pressure-driven intrusion), not taken from parcel w | Parcel models carry no pressure forces; the imposed-flow approach is the standard remedy (Costa et al. 2013; FALL3D/Ash3d umbrella options — Mastin et al. 2014, JVGR). |
| **Initial plume width**: emergent parcels are placed AREA-UNIFORMLY across the collapse dome (radius = `umbrella_r_source`, HT 40 km) rather than at the vent point | The umbrella does not start from a point: the overshoot dome collapse delivers material dome-wide (Van Eaton et al. 2023 imagery — central dome ~40 km scale).  Before this the model built the annulus from r ≈ 0 (column width, km-scale) and under-read the first half-hour by construction. |
| **Wind shielding inside the intrusion** (`umbrella_wind_shield`): the resolved-wind displacement of in-cloud parcels (r ≤ 1.2R) is cancelled by the factor s = min(1, u_r/\|u_amb\|) — full shielding while the outflow dominates, fading automatically as Ṙ decays; no time constant | A gravity current of this intensity EXPELS the ambient flow — superimposing archived-met advection on the imposed radial flow double-counts the ambient wind inside the current.  Observational basis: the HT umbrella spread as a vent-centred circle for the first hours while the outflow (40–90 m/s) dwarfed the ~10–25 m/s stratospheric winds, and only began drifting as the current weakened (Van Eaton et al. 2023).  The s-factor makes the handoff to ambient advection self-regulating.  Turbulent displacements are NOT shielded. |

## 3. In-plume washout (`PlumeWashout`, two-phase)

| choice | basis |
|---|---|
| Henry-equilibrium gas scavenging into condensate; effective Henry constants with dissociation | Sander (2015), ACP 15 — Henry's-law compilation (H_eff values); Seinfeld & Pandis (2016), ch. 7. |
| SO₂ survives (~99 %): S(IV) dissociation suppressed at plume pH ≈ 2 (co-erupted HCl/H₂SO₄ acidify the condensate) | Standard aqueous S(IV) equilibria (Seinfeld & Pandis 2016); volcanic-plume acidity argument as in Tabazadeh & Turco (1993), Science 260. Consistent with the observed HT stratospheric SO₂ 0.41±0.02 Tg (Millán et al. 2022, GRL 49 **[verified]**) for an at-vent load of the same order. |
| HCl mostly scavenged (~2 % survives): solubility-unlimited uptake (H_eff ~ 2×10⁶ M/atm at plume pH) | Tabazadeh & Turco (1993), Science 260 — HCl removal in supercooled plume water; **caveat**: Textor et al. (2003), JGR 108 **[verified]** showed ice-dominated (dry-magmatic) columns can carry 25 %+ of HCl to the stratosphere on ice. HT was the opposite regime — a water-flooded phreatomagmatic column — and the OBSERVED HT stratospheric HCl was "unexceptional" (Millán et al. 2022), i.e. small next to the water/SO₂ anomaly: the liquid-scavenging-dominant treatment is the appropriate branch here. Document the regime dependence. |
| Umbrella dehydration to the ~210 K ambient dominates the water budget (69 % of removal); condensate shed above holding capacity with τ_fall; size-dependent settling with updraft recycling | Internal construction (validated: water budget closes exactly; the loft-then-drop ice curtain reproduced; `test_plume_washout.jl` 22/22). Physical anchors: ambient-T saturation at the NBH; Rogers & Yau (1989) for fall-speed parameterisations. |
| Delivered water ≈ 5–6 % of at-vent 3300 Tg → 183–213 Tg vs observed 146±5 Tg | Observation: Millán et al. (2022), GRL 49 **[verified]**; Vömel et al. (2022), Science 377 (radiosonde detection of the moist layers). The at-vent 3300 Tg is an **internal inversion** (chosen so delivered ≈ observed given the model's survival fraction) — state it as such; independent at-vent estimates vary. |

## 4. Trajectory model / estimators (paper §4 material)

| choice | basis |
|---|---|
| Langevin turbulence: Thomson well-mixed criterion; Hanna scheme parameterisation | Thomson (1987), JFM 180; Hanna (1982), in *Atmospheric Turbulence and Air Pollution Modelling*; audited faithful to FLEXPART 11 (Pisso et al. 2019, GMD 12; Bakels et al. 2024, GMD — FLEXPART 11). |
| FLEXPART-class validation metrics (ATMES-II / RANK) | Mosca et al. (1998); Stohl et al. (1998); RANK composite as in FLEXPART evaluations (Hegarty et al. 2013, JAMC — CAPTEX multi-model RANKs 2.3–2.8). Internal: EASLEE CAPTEX RANK ≈ 2.5, ETEX framed-engine RANK 2.75. |
| Adaptive output kernel: per-parcel truncated Gaussian, bandwidth = kNN(12) RMS spacing | kNN density estimation: Loftsgaarden & Quesenberry (1965), Ann. Math. Stat.; variable-bandwidth KDE: Breiman, Meisel & Purcell (1977), Technometrics 19; Abramson (1982), Ann. Stat. 10; Silverman (1986), ch. 5.3. Internal: N^(−1/2) law measured (exponent −0.50…−0.56 on ETEX; unit lattice 0.500 exact); mass conservation exact incl. domain edges. |
| Dispersion bandwidth floor σ(age): Taylor growth σ² = 2σ_v²T_L²(t/T_L − 1 + e^(−t/T_L)) | Taylor (1921), Proc. Lond. Math. Soc. 20. Internal A/B verdict: neutral on the (renormalised) grid deposit, harmful on the (non-renormalised) receptor point-KDE — floor is grid-only; documented in `ETEX_RESULTS.md`. |
| Observed-track ensemble nudging (rigid-shift centroid relaxation) | Newtonian relaxation / observation nudging: Stauffer & Seaman (1990), MWR 118 — applied here to the ensemble CENTROID only (structure-preserving), free-running outside the observed window. Internal: `test_track_nudge.jl` 14/14; HT demo lat RMS 1.12°→0.95°. |
| Mesoscale velocity fluctuations (OU process, fturbmeso·σ) | FLEXPART's mesoscale scheme: Stohl et al. (2005), ACP 5 (§ mesoscale velocity fluctuations, after Maryon 1998). |
| Backward SRR normalisation | Seibert & Frank (2004), ACP 4. |
| Covariance output kernel (`kernel="covariance"`): per-parcel 2×2 horizontal covariance P evolved by the resolved deformation P ← (I+∇u·Δt)·P·(I+∇u·Δt)ᵀ plus turbulent growth 2σ_v²T_LΔt per axis; P read only at evaluation (grid: axis widths, exact separable conservation; receptors: full bivariate form incl. the cxy tilt); R(x) = (det P)^(1/4) | Puff moment-tensor evolution under resolved velocity gradients + turbulence: SCIPUFF second-order closure — Sykes & Henn (1995), JAM 34; Sykes et al. (2004), SCIPUFF Tech. Doc. Hybrid particle–puff carriers: de Haan & Rotach (1998), QJRMS 124 (puff–particle model); Hurley (1994), Atmos. Environ. 28 (PARTPUFF). Tangent-map deformation in Lagrangian transport (F = I+∇u·Δt): Konopka et al. (2004), JGR 109 (CLaMS deformation-driven mixing). Turbulent-growth diffusivity (Taylor limit D = σ_v²T_L): Taylor (1921). Internal: `test_cov_kernel.jl` 27/27 — nilpotent-shear analytic P exact, rotation-invariance, conservation 1e-12 incl. edges, receptor bivariate value = analytic (tilt ratio 2.63 = analytic 2.63). |
