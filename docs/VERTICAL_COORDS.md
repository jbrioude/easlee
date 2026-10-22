# EASLEE — Vertical coordinates and transport

Status: 2026-10-22 (stacked layout contract; readers since commit b778cec).

## The core idea: transport in level-INDEX space

The transport kernel never works in a named physical coordinate.  A parcel's
vertical position is the **continuous level index** `indl` ∈ [1, Nk)
(1 = the ground), height is linear in `indl` between integer slots, and every
met source is reduced at read time to fields on two interleaved sets of
levels — the **stacked layout**:

| class | fields | where slot k sits | correct read |
|---|---|---|---|
| H (half levels) | `zlev`, `levdot`, `tke` (and `tls`) | AT half level k, height `zlev[k]`; `zlev[1] = 0` is the ground | at `indl = k`: `interpolation4Dlevsingle`, `get_vertical_profile*` paired with `zlev[k]`, `ztoindl` / `indltoz` |
| F (full levels) | `Ux`, `Vy`, `T`, `puv`, `rho`, `theta`, `qv`, `cc`, `tpf`, `umf`, `dmf`, `udr`, `ddr`, `ttlwr`, `ttswr`, `ttpm` | BETWEEN half levels k−1 and k, read at `indl = k − ½`; height `zf(k) = ½(zlev[k−1] + zlev[k])` for k ≥ 2 | the stacked read `interpolation4Dsingle` (`_split_z_stagger`: index z + ½ for z ≥ 1.5, 2z − 1 below), `interpolateU/V` |

| Ingredient | Meaning | Used by |
|---|---|---|
| `zlev(x,y,k,t)` | height above ground of HALF level k [m] | `ztoindl` (release, turbulence, reflection), output binning, profiles |
| `levdot(x,y,k,t)` | **index** velocity dindl/dt [1/s] on half level k; 0 at the ground (k = 1) and the top (k = Nk) — coordinate surfaces | advection |
| `puv`, `rho`, `T`, `theta`, winds, … | FULL level k, between `zlev[k−1]` and `zlev[k]` | everything else |

**Slot 1 of a class-F field holds SURFACE values**, read at `indl = 1` (the
ground) and blended linearly to slot 2 at `indl = 1.5`: U, V = slot 2 (constant
extension; the Monin–Obukhov option then shapes the wind below the lowest full
level), T = 2 m temperature, `puv` = surface pressure, `rho`/`theta` derived from
those, `qv` = slot 2, convective fluxes and `cc` = 0, `tpf` = slot 2.

Rules for any consumer:

- The layer between half levels k and k+1 (`indl` ∈ [k, k+1), heights
  [`zlev[k]`, `zlev[k+1]`)) holds full level **k+1**: its air mass is
  `rho[k+1]·(zlev[k+1] − zlev[k])`, its wind slot k+1 — the stacked read at the
  layer's mid index k + ½ (boundary inflow, domain fill, well-mixed checks).
- A profile that pairs heights with class-F values uses the full-level heights
  `zf(k)` (with `zf(1) = zlev[1]`, the ground, for the surface slot) — e.g. the
  plume-rise and volcanic-column soundings (`_full_level_heights` in
  `core/operators.jl`).
- A profile that pairs `zlev[k]` (half levels) with ρ, θ, dρ/dz must evaluate
  them AT the half levels: ρ_h(k) = ½(ρ[k] + ρ[k+1]) for 1 < k < Nk (ρ[1], the
  ground value, for k = 1), dρ/dz_h(k) = (ρ[k+1] − ρ[k]) / (zf(k+1) − zf(k)).
- Equivalently: a class-F field at an arbitrary level index z is read with the
  stacked interpolation (`interpolation4Dsingle`), a class-H field with the plain
  one (`interpolation4Dlevsingle`).  **Any code that pairs `zlev[k]` (or a level
  index k) with a class-F element `[k]` is half a level wrong.**

`ztoindl` inverts the local `zlev` profile (monotonic by construction), so the
kernel is agnostic to whether the levels are hybrid σ-p, eta, isobaric,
fixed heights, or terrain-following — **the coordinate lives entirely in the
reader**, in how `zlev` and `levdot` are built.

## Per-dataset reduction to index space

Every reader delivers the stacked layout above.  Datasets that archive heights
and vertical velocity on the full levels have them moved to the half levels
between them (`_stack_half_levels!`: `zlev[s]` = midpoint of full levels s and
s+1, 0 at the ground, the top half a layer above the top full level;
`_hybrid_p_levdot_stacked!`: ω → the level-index rate of the coordinate velocity
on the half levels).  Terrain-following readers use
`_terrain_following_levdot_stacked!` (w minus the slope term, on the half levels).

| Dataset | Native coordinate | `zlev` (half levels) | `levdot` (half levels) |
|---|---|---|---|
| ERA5 / IFS | 137 hybrid σ-p | ak/bk interface heights, hydrostatic integration of T_v (top half-level p=0 floored at 0.1 Pa) | η̇ (archived on full levels) averaged onto the half levels ÷ the local η spacing |
| MERRA2 | 72 eta (GEOS) | archived `H` − PHIS/g on the full levels → `_stack_half_levels!` | OMEGA → `_hybrid_p_levdot_stacked!` on archived `PL` |
| **NCEP GFS** | **33 isobaric levels** | `gh` − orog, monotonicity enforced below terrain → `_stack_half_levels!` | ω → `_hybrid_p_levdot_stacked!` (ω/Δp on an isobaric half level; the V·∇p of the below-ground ramp removed) |
| AROME ATMO-HDF | 22 fixed heights AGL (20–3000 m), optional isobaric stitch | the heights (and stitched geopotential heights) → `_stack_half_levels!` | 0 (not disseminated) — turbulence only |
| AROME NetCDF | hybrid Ak/Bk | Ak/Bk interface heights (top interface extrapolated) | w on the interfaces → `_terrain_following_levdot_stacked!` |
| MesoNH | Gal-Chen ZHAT/ZS | w-level heights ZHAT·(1 − ZS/ztop) (`zlev[s]` = w-level s+1; w-level 1 is the −10 m ghost) | WT on the w-levels, C-grid slope term |
| WRF | eta (mass-flux archived) | ph+phb geopotential on the w-levels | coordinate velocity on the w-levels (W, or the time-averaged mass-coupled flux) |

The contract test is `validation/check_met_readers.jl` (hydrostatic balance with
the full level of slot s between `zlev[s−1]` and `zlev[s]`).

### Note on GFS isobaric levels

There is no dedicated `PressureLevel` machinery in the transport kernel — and
none is needed: the isobaric ladder is just another irregular level set.  What
matters and is handled explicitly:

1. **level VALUE → index mapping** at read time (`leveltype = isobaric` in the
   namelist builds the descending-pressure ladder; `isobaricInPa`
   stratospheric duplicates are skipped);
2. **zlev from geopotential height minus orography**, with monotonicity
   enforced where pressure surfaces intersect terrain (GFS extrapolates
   u/v/T/q below ground; those slots become a ramp between the surface and the
   first level above ground), then moved to the half levels;
3. **ω → index velocity on the half levels** (`_hybrid_p_levdot_stacked!`) —
   0 at the ground and the top; until 2026-10-22 it was 2ω/Δp on the full
   levels with a fixed ladder spacing;
4. vertical interpolation between isobaric levels is linear in index — the
   same approximation as between hybrid levels (both ladders are irregular
   in z); 33 levels are coarse above 100 hPa, which is a property of the
   dataset, not of the scheme.

End-to-end check: `validation/test_gfs_transport.jl` (parcel displacement vs
wind, vertical sanity).

## The hybrid kinematic/isentropic vertical option

For stratospheric work the kinematic `levdot` is replaced (smoothly) by
isentropic-diabatic motion — see `ISENTROPIC.md`.  The blend happens in
PARCEL space, not in the fields: below σ_r = 0.3 pure kinematic, above
σ_t = 0.1 pure isentropic, cosine ramp between, preserving the well-mixed
criterion through the tropopause.

### ⚠ Any other process that moves a parcel vertically MUST update `theta_p`

Under the hybrid/isentropic option the parcel's level index is **re-derived
every step from its prognostic potential temperature**:

```julia
indl_th = _theta_to_indl(th, th_tgt, 2, nz)      # th_tgt = theta_p + thdot·dt
indl    = indl + w * (indl_th - indl)            # isentropic.jl
```

so a vertical displacement written only to `state.indl` is **undone on the next
step** — and with w → 1 in the stratosphere the cancellation is *exact*, i.e.
the process becomes a silent no-op rather than a degraded one.

Two processes have hit this:

| process | status |
|---|---|
| `plume_radiation` (radiative descent) | handled from the start — takes `thetap` |
| `apply_species_settling!` (gravitational fall) | **was broken until 2026-07-31** |

The correct treatment depends on the physics:

- **Cross-isentropic motion** (sedimentation: the particle falls *through* the
  air, it does not carry its isentrope) → rebase `theta_p` on the ambient θ of
  the level reached.
- **Diabatic motion of the air itself** (radiative cooling) → change `theta_p`
  by the heating rate; the parcel then finds its own new level.

If you add any vertical displacement, ask which of the two it is, pass
`hybridvert.theta_p`, and verify with an early check that the displacement
survives a step.

## Adding a new vertical coordinate

Write the dataset's `met*fieldcalculation!` so that it fills `zlev` (monotonic
heights above ground on the HALF levels, `zlev[1] = 0`), `levdot` (index
velocity on the same half levels, 0 at the ground and the top) and the
class-F fields on the full levels stacked between them, with the surface values
in slot 1 (see the table at the top) — nothing else in the model needs to
know.  Run `validation/check_met_readers.jl` on it.  The `SigmaPressureLevel` level type is the generic container
(`ztoindl` only reads `zlev`); the special types (`PressureLevel`,
`SigmaGeopotentialLevel`, …) in `types/level.jl` exist for future dispatch but
are not required by the current kernel.
