# EASLEE — Vertical coordinates and transport

Status: 2026-06-13.

## The core idea: transport in level-INDEX space

The transport kernel never works in a named physical coordinate.  A parcel's
vertical position is the **continuous level index** `indl` (1 = surface slot,
Nk = top), and every met source is reduced at read time to the same three
ingredients:

| Ingredient | Meaning | Used by |
|---|---|---|
| `zlev(x,y,k,t)` | height AGL of level k [m] | `ztoindl` (release, turbulence, reflection), output binning, profiles |
| `levdot(x,y,k,t)` | **index** velocity dindl/dt [1/s] | advection |
| `puv`, `rho`, `T`, `theta` | thermodynamics on the same k | everything else |

`ztoindl` inverts the local `zlev` profile (monotonic by construction), so the
kernel is agnostic to whether the levels are hybrid σ-p, eta, isobaric,
fixed heights, or terrain-following — **the coordinate lives entirely in the
reader**, in how `zlev` and `levdot` are built.

## Per-dataset reduction to index space

| Dataset | Native coordinate | `zlev` construction | `levdot` construction |
|---|---|---|---|
| ERA5 / IFS | 137 hybrid σ-p | ak/bk hydrostatic integration of T_v (top half-level p=0 floored at 0.1 Pa) | etadot ÷ ∂p/∂k (local layer spacing) |
| MERRA2 | 72 eta (GEOS) | archived `H` − PHIS/g (no integration) | OMEGA · 2/(p[k+1]−p[k−1]) on archived `PL` |
| **NCEP GFS** | **33 isobaric levels** | `gh` − orog, monotonicity enforced below terrain | ω · 2/Δp on the FIXED ladder (1000…1 hPa); never uses sp at k=1 (degenerates over terrain) |
| AROME ATMO-HDF | 22 fixed heights AGL (20–3000 m) | the heights themselves | 0 (not disseminated) — turbulence only |
| MesoNH | Gal-Chen ZHAT/ZS | ZHAT·(1 − ZS/ztop) | WT ÷ local Δz |
| WRF | eta (mass-flux archived) | ph+phb geopotential | AVGFLX_WWM ÷ column mass · (Nk−1) |

### Note on GFS isobaric levels

There is no dedicated `PressureLevel` machinery in the transport kernel — and
none is needed: the isobaric ladder is just another irregular level set.  What
matters and is handled explicitly:

1. **level VALUE → index mapping** at read time (`leveltype = isobaric` in the
   namelist builds the descending-pressure ladder; `isobaricInPa`
   stratospheric duplicates are skipped);
2. **zlev from geopotential height minus orography**, with monotonicity
   enforced where pressure surfaces intersect terrain (GFS extrapolates
   u/v/T/q below ground; those pseudo-levels are never reached by parcels
   because `zlev ≥ 0` with ground reflection);
3. **ω → index velocity on the FIXED ladder spacing** — using the
   surface-pressure slot in the centred difference would blow up over
   high terrain (found and fixed during validation);
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
| `apply_species_settling!` (gravitational fall) | **was broken until 2026-07-31** (`fb16b1c`) |

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
heights AGL) and `levdot` (index velocity) — nothing else in the model needs
to know.  The `SigmaPressureLevel` level type is the generic container
(`ztoindl` only reads `zlev`); the special types (`PressureLevel`,
`SigmaGeopotentialLevel`, …) in `types/level.jl` exist for future dispatch but
are not required by the current kernel.
