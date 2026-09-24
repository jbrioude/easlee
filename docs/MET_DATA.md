# EASLEE — Eulerian met data reference

Status: 2026-06-12.  This file documents every meteorological dataset EASLEE
can ingest, what each provides, what is derived, and what each enables
(turbulence, convection, wet deposition).

---

## How a dataset is wired in

Each dataset is described by a **variables namelist** in `listvariables/`,
selected per run by the config's `meteorology.variables` field (see
CONFIG.md).  The namelist maps input variable names to EASLEE
field names and declares dataset properties:

```
#info
model = GFS                      # selects the derived-field calculation
format = grib                    # grib | netcdf
leveltype = isobaric             # hybrid (default) | isobaric | height
multistep = true                 # several valid times per GRIB file
companions = __HP1__>__SP1__     # filename patterns for split datasets
#3D
Ux = u                           # EASLEE name = input name
...
#2D
sp = sp
...
```

- Fields listed with a **dummy mapping** (`rho = rho`) are allocated but never
  read — they are filled by the model-specific calculation in
  `easlee/src/io/metcalc_models.jl` (or `metgribfieldcalculation!` for
  hybrid-level ECMWF data).
- **companions** handles datasets whose fields are split over several files at
  one valid time (substitution applied to the file name).
- Vertical convention everywhere after loading: `k = 1` surface slot
  (`zlev = 0`), `k = 2 … Nk` ascending.  All vertical motion is expressed as
  `levdot` [level-index/s].

A run selects a dataset entirely from its TOML `[[meteorology]]` block
(`model`, `path`, `variables`, `available = "auto"`); example configs for
each dataset are in `configs/`.

---

## Supported datasets

### ERA5 model levels (ECMWF), `listvariables/ERA5.vars` / `ERA5_CONV.vars`
- **Files**: GRIB, one valid time per file. `ERA5_CONV` additionally reads the
  `IFS_EA_CLD_*` companion of each `IFS_EA_MET_*` file
  (`companions = _MET_>_CLD_`).
- **Grid**: regular lat-lon (tested: 0.25°, 60°N–60°S, global wrap).
- **Vertical**: 137 hybrid levels; `zlev` from ak/bk hydrostatic integration
  (top half-level p=0 floored at 0.1 Pa); `sp` from `lnsp` when needed.
- **Direct fields**: u, v, etadot, T, q, lnsp/sp, 2t, 2d, z, lsm, zust→ustar,
  blh→pblh.  CLD companion: `avg_umf/dmf/udr/ddr` (draught mass fluxes &
  detrainment), `avg_tpf` (3-D precipitation flux), `cc`.
- **Enables**: full turbulence (Hanna from archived ustar/blh), **convection**
  (complete mass-flux set), **wet deposition** (3-D tpf + cc).

### MERRA2 (NASA GEOS), `listvariables/MERRA2.vars`
- **Files**: NetCDF4, 8 valid times per daily file (tavg3, 01:30+3h…).
  Companion `companions = _asm_Nv.>_mst_Ne.` (and `trb_Ne` if mapped).
- **Grid**: regular lat-lon 0.625°×0.5° (576×361), global wrap.
- **Vertical**: 72 eta levels, archived **top-down** (flipped by the loader);
  `PL` (mid-level pressure) and `H` (mid-level height) archived directly —
  no hydrostatic integration; `zlev = H − PHIS/g`.
- **Direct fields**: U, V, T, QV, OMEGA→levdot, PL, H, CLOUD→cc, PS, PHIS;
  mst_Ne: `CMFMC`→umf (updraught mass flux on the 73 edges),
  `PFLCU+PFICU+PFLLSAN+PFILSAN`→tpf (convective + large-scale rain + ice).
- **Derived**: rho, theta, drhodz, levdot from OMEGA/(∂p/∂k);
  ustar/ol/pblh from profile fallbacks (no flux collection on disk —
  add `tavg1_2d_flx_Nx` and map `ustar`/`pblh` for archived values).
- **Enables**: **convection** (umf only: entrainment from ∂M/∂z, detrainment
  flux-consistent — no downdraughts), **wet deposition** (3-D tpf + cc),
  turbulence from derived surface parameters.

### NCEP GFS 0.25° pgrb2, `listvariables/GFS.vars`
- **Files**: GRIB2, one valid time per file (`date/cycle/gfs.tCCz.pgrb2.0p25.fFFF`).
- **Grid**: regular lat-lon 0.25° (1440×721), global wrap.
- **Vertical**: 33 isobaric levels (1000→1 hPa); the 8 `isobaricInPa`
  stratospheric duplicates are skipped.  `zlev = gh − orog`, monotonicity
  enforced where pressure levels run below terrain.
- **Direct fields**: u, v, t, q, w(ω)→levdot, gh, tcdc→cc (GRIB.jl alias of
  tcc), sp, orog, 2t, 2d, `avg_utaua/avg_vtaua` (stress), `avg_ishf`
  (sensible heat flux), lsm, fsr→z0, prate.
- **Derived**: rho/theta/drhodz; levdot = ω·2/Δp on the fixed ladder;
  ustar from stress, ol from heat flux, pblh from bulk-Richardson profile;
  pseudo-3-D tpf from `prate` distributed below cloud top.
- **Enables**: turbulence (full derived chain), **wet deposition**
  (pseudo-3-D tpf + cc), **diagnostic convection**: no archived mass fluxes,
  but `cape`/`cin` are mapped — `diagnostic_convection!` (MPTRAC-style
  extreme convection: CAPE-gated, parcels below the θe-derived equilibrium
  level redistributed uniformly in air mass; `DiagConvParams.tau_mix`
  softens the extreme limit).  Pass `diagconv=[DiagConvFields(rt.metfields2D[n]) …]`
  to `step_model!`.

### AROME ATMO-HDF (Météo-France 0.025°), `listvariables/AROME_HDF.vars`
- **Files**: GRIB2, **7 valid times per file** (`multistep = true`; messages
  filtered by validity time).  Packages at the same valid time are companions:
  `__HP1__>__SP1__`, `__HP1__>__SP2__`, `__HP1__>__SP3__`.
- **Grid**: regular lat-lon 0.025° (1121×717), France domain.  Corners outside
  the native Lambert domain are bitmap-missing → stored as NaN (≈17 % of the
  box) and ignored by the physics kernels.
- **Vertical**: 22 fixed height-above-ground levels, **20–3000 m** —
  boundary-layer dataset.  `zlev` is exactly those heights.
- **Direct fields**: HP1: u, v, t, r (RH %), pres; SP1: 2t, 2r, prmsl;
  SP2: blh→pblh; SP3: iews/inss→stress, sshf (sign flipped to +up).
- **Derived**: qv from RH (Magnus), rho/theta/drhodz; ustar from stress,
  ol from sshf.  **levdot = 0** — the dissemination has no vertical velocity;
  vertical transport is turbulence-only.
- **Enables**: turbulence with archived blh + real fluxes.  No convection.
  Wet deposition **not yet enabled**: `tp` is accumulated-from-start;
  per-step differencing is not implemented.
- **Caveats**: domain top 3 km (parcels are capped there); BL studies only.
  Accumulated/averaged SP3 fields (stress, heat fluxes) have **no step-0
  message** — at the cycle's analysis time the surface parameters come from
  the derived fallbacks (log-law ustar, Rib-based L) and switch to the
  archived values from forecast hour 1 onward.

### MesoNH (`_reduced` files), `listvariables/MESONH.vars`
- **Files**: NetCDF4, one valid time per file.
- **Grid**: **Mercator** conformal projection (RPK = 0) — new grid type in
  `types/grid.jl` + `met/coords.jl` (exact lon/lat ↔ index, scale-factor-aware
  displacements).  Lambert (RPK ≠ 0) and rotated (BETA ≠ 0) MesoNH domains are
  not wired yet and error out explicitly.
- **Vertical**: Gal-Chen terrain-following: `zlev = ZHAT·(1 − ZS/ztop)`;
  MesoNH ghost level 1 becomes the surface slot.
- **Direct fields**: UT, VT, WT→levdot, THT→theta, PABST→puv, RVT→rv,
  **TKET→tke (native TKE — no Hanna needed)**, ZS.
- **Derived**: T from θ/Exner, qv from rv, rho, drhodz, levdot from WT/Δz;
  pblh from bulk-Richardson profile, ustar from log-law (for the Hanna
  option-b path and the entrainment-zone scheme).
- **Caveats**: UT/VT are on C-grid u/v points but treated as mass points
  (½-cell offset, ~1 km at dx = 2 km); no precipitation fields in the
  `_reduced` files → no wet deposition; `BOR21.1.00-01.002_reduced.nc` on disk
  is corrupt (bad compressed chunks) and must be excluded from the met-file
  listing.

### WRF and AROME NetCDF
Unchanged in v1.1 (`listvariables/WRF.vars`, `AROME.vars`): Lambert
grids, `inputdimensions`-based allocation, AVGFLX mass-flux winds for WRF.

---

## Reader guards and unit detection (review 2026-09-04 A3)

### Surface-flux units: sshf, ewss/nsss (`initialize.jl`, `_sfcflux_convention!`)
The Obukhov length and (when `ustar` is not archived) the friction velocity
are DERIVED in `metgribfieldcalculation!` from the archived sensible heat flux
`sshf` and surface stress `ewss/nsss` for every name map that carries them
without `ustar`/`ol`.  In this tree that derivation runs on the hybrid-level
GRIB path only (`_load_met_file!`: `leveltype = hybrid`, model ≠ AROME), i.e.
**`listvariables/ERA5.vars`** (the CAPTEX/ETEX configs and every other config
on that map).  `WRF.vars` / `WRFv2.vars` (`sshf = HFX`) and `AROME.vars`
(`sshf`) map the same names, but their calculations (`checkmetfields!`,
`metaromefieldcalculation!`) do not derive `ol` from `sshf` today; the detector
below is model-agnostic so it also covers them (instantaneous W/m²) if they
are ever routed through the derivation.  `ERA5_CAPTEX`, `ERA5_CONV`,
`ERA5_FULL` and `ERA5_CAPE` read `zust`/`blh` and carry no `sshf`, so they
never enter it.

The **GRIB header cannot be trusted** for the unit: ERA5 files carry
`stepType=accum`, `units = J m**-2` / `N m**-2 s` whether the values are raw
MARS accumulations (J/m² over the step) or flex_extract de-accumulated
per-second values (W/m², N/m²).  The archive on the reference workstation
(`/data/ECMWF/ERA5/EAg/2013/EAg13091700`) has `sshf ∈ [−405, 223]` and
`|τ| ≤ 4` under an "accumulated" header — physical only as W/m² and N/m².
Dividing those by 3600 (the old unconditional rule) made |L| 3600× too large
(neutral surface layer everywhere); taking a raw accumulated archive as
instantaneous would make u* ~60× too large.  WRF `HFX` and AROME fluxes are
instantaneous W/m² and would be broken the same way by a ÷3600.

The reader therefore decides **from the field magnitude, once per run**
(first file carrying non-zero fluxes), with no model-name branch:
`max |sshf| > 5000 → J/m² accumulated`, else W/m²; `max |τ| > 50 → N s/m²
accumulated`, else N/m².  The accumulation period is `endStep − startStep`
of the sshf message (header); when it is 0, 3600 s is assumed with a warning.
The decision is logged once:
`surface fluxes: sshf treated as <W/m² | J/m² accumulated over N s> (max |sshf| = …)`.
Override with the namelist key `#info sfcflux_units = auto | per-s | accum`
(same pattern as `etadot_units`).

### Other reader guards added in the same review
- **Latitude scan direction**: GRIB fields are flipped to an ascending latitude
  axis only when `jScansPositively == 0`; the axis is asserted strictly
  ascending after the flip (`gribinfo`).  The decode writes Float32 into the
  slot in one pass (missing-value test inline, flip as index reversal) instead
  of two Float64 temporaries per message; bit-identical output.
- **Incomplete GRIB files**: per-file level counts per 3-D variable
  (`GribScan`).  The first file sets the reference; a later file with fewer
  levels (truncated file, absent companion, level subset) warns naming file,
  variable and levels found/expected — the slot would otherwise silently keep
  the previous file's field.  Zero levels for T/U/V is an error.
- **Run start outside the archive**: the first met bracket is clamped to
  `[1, nfile]` and a start before/after `[filedate[1], filedate[end]]` errors
  naming the start date and the span (was a BoundsError / MethodError).
- **End of archive with prefetch**: `readmetfiles_prefetch!` now holds the
  last bracket ("nothing to read (end of met archive — holding the last
  frame)") like the earlier reader instead of indexing `filename[nfiles+1]`.
- **WRF `Times`**: the time-slice search loops over the time dimension and
  errors when the requested date is not in the file (it previously walked the
  19 characters of the first date string).

---

## Capability matrix

| | ERA5(+CLD) | MERRA2 | GFS | AROME-HDF | MesoNH | WRF |
|---|---|---|---|---|---|---|
| zlev source | ak/bk integr. | archived H | gh − orog | fixed heights | Gal-Chen | ph/phb |
| levdot | etadot | OMEGA | ω | — (0) | WT | AVGFLX_WWM |
| tke | Hanna | Hanna | Hanna | Hanna | **native** | Hanna |
| ustar | zust / stress | log-law* | stress | stress | log-law* | derived |
| Obukhov L | sshf | Rib* | heat flux | sshf | Rib* | derived |
| pblh | blh | Rib profile* | Rib profile* | blh | Rib profile* | derived |
| convection (umf/dmf/D) | ✓ full | umf only | CAPE-gated diag.* | — | resolved | — |
| wet dep (3-D tpf) | ✓ | ✓ | pseudo-3-D* | — (tp accum) | — | — |
| land-sea mask | lsm | — | lsm | — | — | — |
| roughness z0 | — | — | fsr | — | — | — |

`*` = derived fallback implemented in `metcalc_models.jl`
(`_surface_turbulence!` chain: stress → log-law ustar; heat flux → Rib-based
L; archived → Rib-profile pblh).

---

## Adding a new dataset

1. Create `listvariables/<NAME>.vars` mapping its variables (`#info` model
   name, format, leveltype, companions as needed).
2. If the grid/levels need new handling, add a `netcdfinfo_<NAME>` (NetCDF) or
   extend the GRIB ladder logic in `gridcheck_onefile!`, and add a
   `met<name>fieldcalculation!` in `metcalc_models.jl` deriving
   rho/theta/zlev/levdot/drhodz and calling `_surface_turbulence!`.
3. Hook the model name into `_load_met_file!` (initialize.jl).
4. Point a TOML config's `[[meteorology]]` block at the data
   (`available = "auto"`) and smoke-test a short transport run.
