# Biomass-burning plume rise (injection height)

> **Status: implemented standalone (config-driven, single fire from the TOML).**
> `easlee/src/physics/plume_rise.jl` (Sofiev default + Freitas option + θ-profile PBL
> diagnostic); wired into the release seeding via `[[release]].[plume_rise]`. Validated:
> 500 MW fills the PBL ([31,1437] m), 5 GW lofts above it ([1739,2840] sofiev /
> [2097,3022] freitas). Sofiev coefficients verified vs the paper (2026-07-02, Eq. 13).
> Inventory (many fires) wiring = refactor slice 4.
>
> **Tested end to end (2026-07-02, `validation/test_plume_rise.jl`, ERA5 sounding at the
> era5_conv release cell, H_abl = 1438 m):** both schemes monotone over FRP 10 MW–20 GW;
> Sofiev 500 MW stays in the PBL ([0, 1088] m), 5 GW lofts ([1438, 2008] m); doubling FT
> stability lowers the Sofiev top (1088→714 m).  Wired path: 100 % of 20 k released
> parcels land inside — and fill — the computed injection layer, for both schemes at
> 0.5/5 GW.  **Caveat (simplified Freitas):** with `fire_area` fixed the low-FRP response
> is weak (top 1.9 km at 10 MW vs 3.3 km at 20 GW) because the top-hat column's initial
> radius R₀=√(A/π) and fixed w₀ = 1 m/s dominate over buoyancy for weak fires — the full
> Freitas PRM ties w₀ to the fire's convective energy flux.  Prefer Sofiev (default) for
> FRP-inventory releases; treat Freitas as the intense-fire/pyroCb option.
>
> ```toml
> [[release]]
>   [release.plume_rise]
>   scheme = "sofiev"   # none | sofiev | freitas
>   frp    = 5.0e8      # W ;  fire_area = 2.0e6  # m² (freitas)
>   bottom = "pbl"      # pbl | ground | fraction:<f>
> ```

## Original proposal

Goal: given a fire's **FRP** (fire radiative power), **fire size**, and the ambient
state, decide the **injection height(s)** over which the emitted parcels are released,
instead of dumping them at the surface. Wires into the slice-4 `inventory` release kind.

## Landscape (what exists, incl. newer than Freitas-2007)

| scheme | inputs | cost | notes |
|---|---|---|---|
| **Freitas et al. 2007/2010** (1-D PRM) | fire area + convective heat flux (biome LUT or from FRP) + ambient profile | **high** — a 1-D entraining-plume PDE per fire/timestep | physical; MISR-validated; the workhorse in WRF-Chem/BRAMS/FLEXPART. Gives a full plume *profile* (detrainment layer → top). |
| **Sofiev et al. 2012 / IS4FIRES 2013** | **FRP** + PBL height + free-trop Brunt–Väisälä N² | **very low** — one algebraic formula | semi-empirical, CAPE-analogy, fitted to MISR; used operationally in SILAM; the natural match for FRP fire products (GFAS/VIIRS/MODIS). |
| **Briggs** (classic) | buoyancy flux (from heat flux) + stability/wind | low | point-source stack formula; sometimes applied to fires; simplest. |
| **Paugam et al. 2015/2016** (PRMv2) | remote-sensing-optimised Freitas (FRP-derived fire properties) | high | review + improved 1-D PRM driven by satellite inputs. |
| **2024 refinements** (Ma 2024 CAM-chem; Lu 2024 E3SM; GEOS 2024) | improved 1-D PRM, interactive ESM coupling, pyroCb | high | mostly target extreme **pyro-convection**; refinements of the Freitas-type PRM. |

**Takeaway:** the main "newer than Freitas" option is **Sofiev (2012/2013)** — a fast
FRP-based algebraic injection-height formula. Post-2020 work largely *refines the 1-D
PRM for pyroCb* rather than replacing the Freitas/Sofiev pair. Data-driven/ML height
estimators exist in the literature but aren't standard in LPDMs.

## Proposed for EASLEE: Sofiev default + Freitas option

Two schemes behind one config switch (consistent with the engine refactor — a
`[release.plume_rise]` block on an `inventory` release; no script):

```toml
[[release]]
kind      = "inventory"
from_file = "gfas_20250520.nc"      # gridded FRP / emissions
  [release.plume_rise]
  scheme    = "sofiev"      # none | sofiev | freitas | briggs
  frp       = "frp"         # inventory field: per-fire/cell FRP [W]
  fire_area = "area"        # [m²]  (freitas/briggs only)
  profile   = "uniform"     # uniform | gaussian — how parcels fill the plume depth
  bottom    = "pbl"         # pbl | ground | fraction:<f>  — lower edge of the layer
```

### Default — Sofiev et al. (2012) semi-empirical formula
```
H_p = α·H_abl + β·(FRP / P_f0)^γ · exp( −δ · N²_FT / N²_0 )
```
- `H_abl` = boundary-layer height (ERA5 `blh`); `N²_FT` = Brunt–Väisälä² in the free
  troposphere (from the θ profile above the PBL); `FRP` per fire from the inventory.
- Reference scales `P_f0 = 1 MW`, `N²_0 = 2.5e-4 s⁻²`.
- Coefficients (fitted to MISR): **α=0.24, β=170 m, γ=0.35, δ=0.6**
  — ✓ verified against Sofiev et al. (2012), Eqs. (10)+(13) (2026-07-02).  The paper's
  two-step variant (separate FT-plume parameters, Eq. 19) scores slightly worse
  overall (~64 % vs ~66 % within 500 m) — the one-step form implemented here is
  the paper's better-scoring default.
- Strong stable layers (large `N²_FT`) suppress rise (exp term → small) → plume stays
  near/in the PBL — matches the MISR finding that most plumes don't escape the PBL.

Cheap (one evaluation per fire) → scales to a whole FRP inventory.

### High-fidelity option — Freitas et al. (2007/2010) 1-D PRM
Integrate the 1-D entraining-plume equations (vertical velocity, buoyancy, continuity)
forced by the fire's convective heat flux (from `fire_area` × biome heat flux, or from
FRP via a radiative→convective fraction) in the ambient sounding; the plume top is where
`w → 0`. Returns a **detrainment layer** (neutral-buoyancy level → top), not just a top.
More physical for intense fires / pyroCb; costs a small ODE integration per fire.

**What "simplified" means in EASLEE's implementation:** the full Freitas PRM carries
water vapour + a bulk cloud-microphysics closure (cloud water / rain / ice) whose
condensation releases latent heat inside the column — the term that powers pyroCb
deepening.  The EASLEE version is the **dry** MTT (Morton–Taylor–Turner) core only: no
qv equation, no microphysics, no latent-heat boost (flagged "v2" in `plume_rise.jl`).
NB the PRM's microphysics is a *thermodynamic* closure — it does not use smoke-aerosol
number to activate CCN/IN (aerosol–microphysics coupling of that kind lives in
interactive ESM versions, e.g. Lu 2024 E3SM); dropping it removes latent heating, not an
aerosol budget.  Smoke CCN/IN properties belong to the released species (species DB
`ccn_aero`/`in_aero`), not to the injection-height calculation.

### Releasing the parcels
Both schemes yield an injection layer `[z_bottom, H_p]` (Freitas: NBL→top; Sofiev:
`bottom`→`H_p`). The engine distributes the cell's emitted parcels over that layer per
`profile` (uniform, or Gaussian peaked near the top). This is a release-placement option
— it reuses the existing parcel-seeding machinery; nothing else in transport changes.

### Inputs & wiring
- From the **inventory**: per-cell `FRP` (W), `fire_area` (m²), emission mass/flux.
- From the **met**: PBL height (`blh`), θ profile (→ `N²_FT`), and full sounding (Freitas).
- New code: `src/physics/plume_rise.jl` with `plume_top_sofiev(...)` and
  `plume_layer_freitas(...)`; called from the `inventory` release seeding in `build_run`/
  `run_spec!`. Default `scheme="none"` (surface release) preserves current behaviour.

## Sources
- Freitas et al. (2007), ACP — 1-D plume-rise model.
- Sofiev et al. (2012), ACP, "Evaluation of the smoke-injection height from wild-land
  fires using remote-sensing data"; Sofiev et al. (2013), GMD — IS4FIRES.
- Paugam et al. (2016), ACP — review + PRMv2 from remote sensing.
- Ma et al. (2024), JAMES — improved 1-D PRM for pyro-convection in CAM-chem.
- Lu et al. (2024), JAMES — interactive fire plume-rise in E3SMv1.
