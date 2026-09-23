# EASLEE — Convective transport schemes

Status: 2026-06-12.  Two parcel-redistribution schemes live in
`easlee/src/physics/convection.jl`, selected by what the met dataset
provides; wet deposition (`wetdep.jl`) shares the same driving fields.

| Scheme | Driving fields | Datasets | Entry point |
|---|---|---|---|
| **Archived mass-flux** (preferred) | updraught/downdraught mass flux + detrainment rates | ERA5 + CLD companion (full set), MERRA2 (`CMFMC`, updraught only) | `convective_redistribution!` |
| **Diagnostic CAPE-gated** (fallback) | CAPE (+ CIN), T/q/p profiles | GFS, any source archiving CAPE | `diagnostic_convection!` |

Convection-permitting drivers (MesoNH 2 km, AROME 1.3 km) resolve convective
motion explicitly — neither scheme should be enabled there.

---

## 1. Archived mass-flux scheme

Monte-Carlo redistribution driven by the IFS/GEOS archived convection
diagnostics, per parcel and per time step (operator-split after advection):

**Fields** (ERA5 `IFS_EA_CLD_*` companion files; MERRA2 `mst_Ne` collection):

| EASLEE | ERA5 | MERRA2 | Meaning |
|---|---|---|---|
| `umf` | `avg_umf` | `CMFMC` | updraught mass flux M_u ≥ 0 [kg m⁻² s⁻¹] |
| `dmf` | `avg_dmf` | — | downdraught mass flux M_d ≤ 0 |
| `udr` | `avg_udr` | — | updraught detrainment D_u [kg m⁻³ s⁻¹] |
| `ddr` | `avg_ddr` | — | downdraught detrainment D_d |

**Algorithm** (per parcel; well-mixed *by construction* — a uniform-in-air-mass
tracer is an exact fixed point):

1. *Entrainment / detrainment* come straight from the **flux divergence** of the
   interface mass fluxes `M_u^if[k] = ½(M_u[k]+M_u[k+1])`:
   `entrain_up[k] = max(M_u^if[k]−M_u^if[k−1],0)·Δt` (air drawn out of layer k),
   `detrain_up[k] = max(M_u^if[k−1]−M_u^if[k],0)·Δt` (air deposited into layer k);
   mirror for the downdraught.  The capture probability is `entrain/dm`.  Because
   the conveyor is built from the flux itself, its flux through every interface
   equals the archived `M_u` exactly — the prerequisite for up- and down-transport
   to cancel.  (The archived `udr`/`ddr` are no longer needed; the divergence is
   entrainment − detrainment.)
2. *Captured parcels* ride the draught and detrain in layer `l` with probability
   `1 − M_u^if[l]/M_u^if[l−1]` (the flux fraction lost there), detraining
   unconditionally where the flux dies (cloud top).  This keeps the physical
   **non-local** transport — sub-cloud air reaches cloud top in one step — so the
   scheme is advective, not diffusive.
3. *Compensating subsidence* is the **exact mass-conserving residual**, not a
   velocity.  With `m1[k]` the air-mass image of the conveyor and
   `S[k] = Σ_{l≤k}(dm[l]−m1[l])` the cumulative deficit below interface k, every
   parcel then descends across interface `c−½` with probability
   `S[c−1]/(dm[c]+S[c−1])`.  The composite (conveyor ∘ subsidence) maps the
   air-mass vector `dm` onto itself in expectation, so a uniform tracer stays
   uniform while the conveyor carries the asymmetric transport.

MERRA2 runs with `umf` only (no downdraught) — a usable degradation, not a
failure mode.

**Validation** (well-mixed convection test, ERA5 2025-05-31,
strongest-`umf` column, 200 k parcels initialised uniform in air mass, Δt = 15
min, 4 h of convective churning):

- BL mixing ratio stays uniform: **BL std(χ) 0.0003 → 0.021**, a flat fixed point
  (every BL level 0.91–1.08 around 1.00; the surface holds 0.98).  The previous
  velocity-subsidence scheme plateaued at **0.155** with the sub-cloud drained to
  χ ≈ 0.5;
- the violation was shown to be a *structural* bias (a Δt = 15/5/1 min sweep all
  plateaued at the same ≈0.15–0.16), not a time-step artefact — hence the
  flux-consistent rebuild rather than sub-cycling.

**Entrainment/detrainment validation** (`validation/test_conv_entrain_detrain.jl`,
one-step Monte-Carlo, 2 M parcels pinned to the strongest-umf ERA5
column): parcel-derived interface
crossings and origin/destination histograms are compared with the ARCHIVED
fields — including `avg_udr`, which the scheme never reads:

- parcel up-flux `F_up(k)` and down-flux `F_dn(k)` (downdraught + compensating
  subsidence) both converge to the archived `Mu_if(k)` as Δt→0:
  R² = 0.85 / 0.97 / **0.997** and Σ-ratio 81 % / 94 % / **101 %** at
  Δt = 900/300/60 s.  The deficit at large Δt is intra-step up-then-subside
  cancellation (nets to zero), not a bias: `F_up ≈ F_dn` to < 1 % at every Δt,
  so up- and down-transport cancel exactly;
- riser origins (entrainment) match the flux-divergence entrainment,
  R² = **0.993** at Δt = 60 s (identical vs the archived-implied `E = D + ΔMu`);
  riser destinations (detrainment) match the flux-divergence profile,
  R² = **0.980**.  Correlation with the raw archived `avg_udr·Δz` is limited by
  construction (R² ≈ 0.42): the IFS archives *gross* detrainment, while the
  flux divergence — the scheme's contract — is the *net* E−D per layer, and the
  two must differ wherever entrainment and detrainment coexist in a layer;
- χ stays uniform after the step (std 0.002–0.005 across Δt).

**Earlier bugs found during convection validation** (all fixed, all pre-existing):

1. `src/gpu/rng.jl` used `^=` for XOR — in Julia `^` is exponentiation, so
   the xorshift state collapsed to 1 and every "uniform" draw returned 1e-7
   (and the GPU `randn` always +3.0).  Affected all previous GPU-turbulence
   runs.
2. `zlev = +Inf` at the model top for full-137-level ERA5 (top half-level has
   p = 0 → `log(p/0)` in the hydrostatic integration); floored at 0.1 Pa.

---

## 2. Diagnostic CAPE-gated scheme (no archived fluxes)

MPTRAC-style "extreme convection parameterization" (Hoffmann et al.) — the
cheap alternative to FLEXPART's Emanuel scheme, which computes its own fluxes
from T/q profiles:

1. **Event gate**: archived `CAPE ≥ cape_min` (default 100 J/kg) and, when
   archived, `CIN ≤ cin_max` (default 50 J/kg).
2. **Equilibrium level**: θe of a lowest-model-level parcel (Bolton 1980) vs
   the environment's *saturated* θe profile; EL = topmost level (≤ 17 km)
   where the lifted parcel is still buoyant.  Specific humidity is
   reconstructed from the typed fields via `Tv = p/(R_d·ρ)` — no additional
   3-D inputs needed.  Events shallower than `zel_min` (2 km) are ignored.
3. **Redistribution**: parcels below the EL are re-placed **uniformly in air
   mass** over [0, z_EL] with probability `p = 1 − exp(−Δt/τ_mix)`;
   `τ_mix = 0` (default) is the MPTRAC extreme limit — the convective column
   mixes completely every time step.  Uniform-in-mass redistribution
   preserves the well-mixed criterion trivially.

This is an **upper bound** on convective vertical transport: every gated
column mixes to its EL.  Use `τ_mix` (e.g. 3–6 h) to temper it, and prefer
the archived-flux scheme whenever the dataset provides one.

**Validation** (well-mixed convection test, GFS 0.25° 2026-06-07
12z, strongest column CAPE = 4223 J/kg, 200 k parcels, Δt = 15 min, 2 h):

- the EL trigger **fires reliably on GFS** — 94 % of parcels moved at step 1
  (mean |Δlevel| = 5.6), the case the diagnostic scheme is meant for;
- the boundary layer stays well mixed: **BL std(χ) 0.0001 → 0.038**, a flat
  plateau (not a drift), comparable to the mass-flux scheme's 0.021;
- one caveat: the **surface layer (z = 0) sits ~10 % low** (χ ≈ 0.90) while
  every level above is spot-on — a bottom-edge convention in the `cm[1]` layer
  mass, stable, not a leak.

**Validation on ERA5** (after the qv fix below; `configs/era5_cape.toml`,
2025-05-31, max-CAPE cell 6568 J/kg, 200 k parcels, Δt = 15 min, 4 h):

- the scheme now fires on ERA5: 93.5 % of column parcels moved in one step
  (one-step test, `validation/test_conv_entrain_detrain.jl`), EL diagnosed at
  16.8 km, **zero** parcels above the EL touched, destinations uniform in air
  mass (χ std 0.011);
- well-mixed: **BL std(χ) 0.0003 → 0.027 flat plateau** over 4 h (every BL
  level 0.91–1.07), comparable to the mass-flux scheme's 0.021.  Before the
  fix this test passed *vacuously* (0 % moved).

**ERA5 EL fix (2026-07-02).**  On hybrid-level ERA5 the scheme used to gate but
find no EL (CAPE 6568 J/kg at the strongest cell, **0 % moved**).  Root cause:
the ERA5 reader stores **dry** density (`rho = p/(R_d·T)`, `initialize.jl`), so
reconstructing the source-parcel humidity from `Tv = p/(R_d·ρ)` returns exactly
q ≡ 0 → θe_parcel ≈ θ → never exceeds the saturated environment θe.  (GFS
worked because its reader stores moist ρ from Tv.)  Fix: `DiagConvFields` now
carries the **archived `qv`** (3-D dict passed in `build_run`), used for the
source-parcel θe; the Tv reconstruction remains only as a fallback for readers
that map no humidity.  `ERA5_CAPE.vars` maps `qv = q` + `cape = cape`.

---

## 3. Usage

```julia
# Archived fluxes (ERA5+CLD, MERRA2) — wet deposition shares the fields:
convfields = [EASLEE.ConvFields(rt.metfields3D[n]) for n in 1:rt.nest]
wd = EASLEE.WetDepState(rt.ctx)
step_model!(state, ctx, date; convfields=convfields,
            wetdepparams=EASLEE.WetDepParams(), wetdepstate=wd)

# CAPE-gated diagnostic (GFS, ERA5 with cape+qv mapped) — pass the 3-D dict too
# so the archived qv feeds the source-parcel θe (required on ERA5, dry-ρ reader):
diagconv = [EASLEE.DiagConvFields(rt.metfields2D[n], rt.metfields3D[n]) for n in 1:rt.nest]
step_model!(state, ctx, date; diagconv=diagconv,
            diagconvparams=EASLEE.DiagConvParams(tau_mix=3600.0f0))
```

Both constructors return `nothing` when the namelist doesn't map the needed
fields, and `step_model!` skips the scheme for that nest.  When both kwargs
are passed, the archived-flux scheme takes precedence.

Time-mean archived fluxes are grid-cell means, consistent with parcels
representing grid-mean air masses; probabilities scale linearly with the
model time step.  The mass-flux scheme is applied in backward runs unchanged
(this Monte-Carlo form is its own adjoint); wet deposition is skipped in
backward mode.

---

## 4. Known limitations / next steps

- **Instantaneous redistribution (FLEXPART-style).** The conveyor is a
  redistribution *within* the step, not a velocity-resolved ascent: where the
  updraught flux grows with height the detrainment probability is zero, so an
  entrained parcel can ride from cloud base to its detrainment level — up to
  cloud top — in a single step.  This is correct for the grid-mean *tracer*
  redistribution (and matches FLEXPART's Forster scheme), but a given parcel's
  intermediate positions are skipped, which matters for trajectory-time
  diagnostics (in-cloud residence time, age spectra).
  A velocity-limited variant (cap the ride at `w_up·Δt`, `w_up = Mu/(ρ·σ_up)`,
  carrying an "in-transit" flag across steps, with the `m1`/`S` bookkeeping
  rebuilt from the capped conveyor to keep well-mixedness) is a possible
  `convection.resolved_ascent` option — a different physics contract, not yet
  built.  NB the full `Mu(z)` profile is still fully used: its *rising* part
  sets entrainment (how many parcels enter per level), its *falling* part sets
  detrainment (where they land).
- Mass-flux/detrainment levels are treated as collocated with the other 3-D
  fields (half/full-level distinction neglected — second-order at 137 levels).
- The compensating-subsidence descent is downward-only (`S` clamped ≥ 0): it is
  *exact* for the usual updraught-dominated column, and simply applies no
  compensation (no error) in the rare net-downdraught column — adding an upward
  branch for that case is a possible refinement.
- The diagnostic scheme uses a single undiluted source parcel from the lowest
  model level; no entrainment dilution of the ascent (consistent with its
  role as an upper bound).
- GFS wet deposition uses the pseudo-3-D `prate` column (below cloud top),
  not a true precipitation-flux profile.
- Possible refinement: gate the diagnostic scheme on `prate > 0` to restrict
  events to actively precipitating columns.
