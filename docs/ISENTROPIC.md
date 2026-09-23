# EASLEE — Isentropic & hybrid vertical transport

Status: 2026-06-13.  Implementation: `easlee/src/physics/isentropic.jl`;
validation: `validation/test_isentropic.jl`.

## Why

Kinematic vertical velocities (etadot/ω) in reanalyses carry assimilation
noise that produces spurious vertical dispersion in the stratosphere — the
validation run measured a kinematic p95 |Δθ| drift of ~40 K in 3 h at
450–550 K.  Transport there is far more accurate on **isentropes** with an
explicit **diabatic** cross-isentropic drift (CLaMS; MPTRAC), which is
~1 K/day and smooth.

## Scheme

Hybrid after Mahowald et al. (2002) / MPTRAC, blended in PARCEL space:

1. **Prognostic θ_p per parcel.**  In the isentropic regime θ_p is never
   re-interpolated from the grid — no cumulative vertical numerical
   diffusion (the CLaMS property).
2. **Along-isentrope motion** — horizontal advection with **u, v archived ON
   the isentropes** (ERA5 `IFS_PT_MET_*` files: 16 θ-levels 265–850 K,
   3-hourly, global 0.25°), read by a self-contained bracket reader
   (`PTWinds`).  Sampling model-level winds at the isentrope's level index is
   NOT equivalent: the horizontal part of that interpolation runs at constant
   model level and aliases ∇θ|_level·∂u/∂θ into the wind (largest at the
   vortex edge and jets).  The PT winds replace the fraction w of the
   horizontal displacement; without PT files the scheme falls back to
   model-level winds with this caveat documented.
3. **Across-isentrope (diabatic) motion** — the potential-temperature
   tendency, evaluated per parcel:
       θ̇ = (θ_p/T) · (∂T/∂t)_diab
   from `avg_ttlwr + avg_ttswr` (ERA5 archived radiative tendencies, the
   leading diabatic term above the tropopause).  When the TOTAL parameterised
   tendency is downloaded, set `ttpm = mttpm` in the namelist (replaces the
   radiative pair — one-line change, no code).
4. **Vertical placement** — after each step the parcel's level index is
   re-inverted from the local, time-interpolated θ ladder at
   θ_p (+ θ̇·Δt), over the full column (θ is monotonic above the boundary
   layer; bracketing the search to the post-kinematic index is forbidden —
   a single etadot-noise excursion would clamp the inversion and ratchet the
   parcel off its isentrope; found during validation).
5. **Blend** — w(σ): 0 for σ = p/p_s > 0.3 (pure kinematic), 1 for σ < 0.1
   (pure isentropic-diabatic), cosine ramp between:
       indl ← indl_kin + w·(indl_θ − indl_kin)
   Both limits are mass-consistent and w is smooth in space/time → the
   well-mixed criterion holds through the tropopause (verified: σ-histogram
   of blend-zone parcels unchanged to ×1.05 over 3 h).

Backward runs: θ̇ enters with the signed Δt (diabatic drift reverses);
the placement/inversion machinery is direction-agnostic.

## Inputs

| Field | Source | Role |
|---|---|---|
| `avg_ttlwr`, `avg_ttswr` | ERA5 CLD companion (model levels) | diabatic θ̇ (radiation) |
| `mttpm` (optional) | ERA5 download | total parameterised θ̇ |
| u, v, pres on θ | `IFS_PT_MET_*` (PT) files | along-isentrope advection |

The 16-level PT files are NOT used as the transport grid (their vertical
spacing, 15–70 K, is coarser than the 137-level grid at these heights, and
they carry no tendency); they supply the winds only.

## Usage

```julia
hyb = HybridVert(length(state.status), rt;
                 sigma_r=0.3, sigma_t=0.1,         # blend bounds
                 diabatic=true,                     # θ̇ on/off (validation)
                 ptdir=..., ptpattern="IFS_PT_MET") # isentropic winds
for date in timesteps
    read_met_step!(rt, date)
    step_model!(state, ctx, date; hybridvert=hyb)
end
```

## Validation (`test_isentropic.jl`)

- archived heating rates: LW ≈ −0.8 K/day, SW up to +3 K/day at θ≈500 K;
- **adiabatic limit** (θ̇ forced 0): parcels must conserve θ against a
  kinematic baseline that drifts several K in 3 h;
- **diabatic**: drift at the integrated-heating scale, not kinematic noise;
- **blend zone**: no pile-up/void at the kinematic/isentropic seam;
- |u_θ − u_ml| diagnostic quantifies the constant-level interpolation error
  the PT winds remove.
