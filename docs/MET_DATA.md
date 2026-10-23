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
- **Edge fields**: CMFMC and the PF* fluxes are on the layer edges; each
  mid-layer slot gets the mean of its two edges (slot 1 = the surface edge).
  Used raw until 2026-10-13, they sat half a layer too high.
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
  stratospheric duplicates are skipped.  `zlev = gh − orog`.  Levels below the
  terrain (p ≥ surface pressure) carry GFS's extrapolation: they are rebuilt as a
  ramp from the surface to the first level above ground (heights evenly spaced,
  T linear, u/v/q/cc those of the first level above ground, p log-linear, no
  vertical motion) — until 2026-10-13 they were
  kept 1 m apart above the ground and the noon PBL over high terrain came out
  at 50 m.
- **Duplicate 2-D messages**: `prate` (average and instantaneous), `cape`/`cin`
  (surface and three layers) — the surface, interval-mean message is kept
  (`_grib2d_score`; the choice is logged), deterministically.  eccodes names
  differ by version (GRIB.jl's bundled eccodes: `tcdc`; eccodes ≥ 2.4x: `tcc`);
  a renamed field that matches no message now warns.
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
- **Vertical**: Gal-Chen terrain-following.  ZHAT are the W-level heights;
  mass level k sits at (ZHAT[k] + ZHAT[k+1])/2 · (1 − ZS/ztop) (2026-10-13; it
  was put at ZHAT[k], half a layer low).  MesoNH ghost level 1 becomes the
  surface slot.
- **Direct fields**: UT, VT, WT→levdot, THT→theta, PABST→puv, RVT→rv,
  **TKET→tke (native TKE — no Hanna needed)**, ZS.
- **Derived**: T from θ/Exner, qv from rv, rho, drhodz, levdot from WT
  averaged onto the mass levels / Δz;
  pblh from bulk-Richardson profile, ustar from log-law (for the Hanna
  option-b path and the entrainment-zone scheme).
- **Caveats**: UT/VT are averaged from their C-grid faces to the mass points
  (2026-10-13; used as mass-point values before); no precipitation fields in the
  `_reduced` files → no wet deposition; `BOR21.1.00-01.002_reduced.nc` on disk
  is corrupt (bad compressed chunks) and must be excluded from the met-file
  listing.

### WRF (`wrfout`), `listvariables/WRFv2.vars`
Reader `metwrffieldcalculation!` (`io/metcalc_models.jl`, rewritten 2026-10-13
after Yue Jia's report — the previous path filled only rho).
- **Grid**: Lambert conformal (MAP_PROJ = 1).  Displacements include the map
  factor (`lambert_mapfac`); winds stay grid-relative (U, V along the grid axes).
- **Vertical**: a surface slot (PSFC, T2, U10/V10) under the mass levels, each
  at its own height: the midpoint of its two PH+PHB interfaces, above the terrain.
- **Fields**: T = (θ′ + 300)·(P/10⁵)^(2/7) with P = P + PB (WRF's R/cp); puv = P;
  q from the QVAPOR mixing ratio; U/V averaged from their faces to the mass
  points; levdot from W averaged onto the mass levels; rho, theta, drhodz.
- **Turbulence**: UST, PBLH archived; Obukhov L from HFX (positive upward in
  WRF); Hanna tke.
- **Not available**: wet deposition (no precipitation flux is mapped; RAINC/RAINNC
  are accumulated from the start and would need per-step differencing),
  mass-flux convection (WRF does not archive it).  `WRF.vars` (AVGFLX variant)
  stops with a message.
- **Tested** on synthetic wrfout files built from GFS with WRF's conventions
  (`validation/make_synthetic_wrfout.py`, `check_wrf_synthetic.jl`) — there was
  no WRF output on the reference machine.

### AROME (Météo-France NetCDF, `AR_T<HHMM>P_<yyyymmdd>_<hh>.nc`), `listvariables/AROME.vars`
- **Times**: from the file attributes `date` + `forecast` (the valid hour); the
  run hour (T<HHMM>P) makes overlapping runs resolve to the freshest one.
- **Vertical**: `Ak`/`Bk` are the layer interfaces (top-down, the ground not
  listed): mass-level pressure = mean of its interfaces, height integrated
  hypsometrically (2026-10-13; the interfaces were used as mass levels before,
  hydrostatic balance off by −2.3 %).
- Lambert grid (map factor applied); native TKE (`tke_kind = e`).

---

## Reader guards and unit detection (review 2026-09-04 A3)

### Surface-flux units: sshf, ewss/nsss (`initialize.jl`, `_sfcflux_convention!`)
The Obukhov length and (when `ustar` is not archived) the friction velocity
are DERIVED in `metgribfieldcalculation!` from the archived sensible heat flux
`sshf` and surface stress `ewss/nsss` for every name map that carries them
without `ustar`/`ol`.  In this tree that derivation runs on the hybrid-level
GRIB path only (`_load_met_file!`: `leveltype = hybrid`, model ≠ AROME), i.e.
**`listvariables/ERA5.vars`** (no shipped configuration uses it; the CAPTEX/ETEX
configurations use `ERA5_CAPTEX`).  WRF (`WRFv2.vars`: `shf = HFX`) derives L
in `_derive_surface_params!` from a positive-upward flux; `AROME.vars` (`sshf`)
does not derive `ol` from `sshf` today; the detector
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

**Accumulation windows, every GRIB dataset (2026-10-22).** A 2-D field with
`stepType = accum` or `avg` whose header gives a window (`startStep < endStep`:
raw MARS/CDS hourly files, IFS or AROME-HDF fields accumulated from the
forecast start, GFS 0-3/0-6/6-9/6-12 h buckets) is converted to the mean rate
over the last archived step: rate = (I(e) − I(e′)) / (e − e′), with I = value
(accum) or value·(e − s0) (avg), and e′ the previous step of the same run found
in the same file or in a file of the same name pattern in the same directory
(log line `window → rate: …`). Without such a step the mean over the whole
window is used (exact for a one-step window). flex_extract files (step 0) are
untouched. Whether an `accum` header can be believed is decided **once per met
source**, under the assumption that the files of one directory and name pattern
follow one convention: `sfcflux_units = accum` → the header is believed;
`per-s` → the values are taken as rates and the windows ignored; `auto` → from
the heat-flux / stress fields of the first file that carries them over a window
(above the thresholds above → accumulated; all below → the values were
de-accumulated under the original header, warned). A declaration that
contradicts the magnitudes is applied, with a warning. `avg` fields are rates in
every case; only their window is converted. Fields converted from their window
are not divided again by the period above (before this, `sfcflux_units = accum`
on raw hourly files divided twice: u* 60× too small). Verified on one ERA5 hour
re-encoded four ways (flex_extract, raw 1-h accumulations, accumulations from
the forecast start, an accumulated header over per-second values) under the
three declarations: u* and L agree with the flex_extract reading to the GRIB
packing precision wherever the declaration is right
(`/data/jbrioude/met_samples/fluxconv/check.jl`).

### Heat-flux SIGN (2026-10-13)
The Obukhov length L = −ρ·cp·T·u*³/(κ·g·H) needs H positive UPWARD (a heated
surface gives L < 0).  ECMWF archives `sshf` positive DOWNWARD (−236 W/m² at
noon over the Sahara in a flex_extract ERA5 file), and the hybrid-GRIB path used
it raw: every daytime surface layer was stable and the derived pblh had an
inverted diurnal cycle (reported by Yue Jia).  The derivation now flips it
(`#info sshf_positive = down`, the default on that path; `up` for a dataset
already positive upward) and logs the median L over heating and cooling cells
once per run.  GFS `avg_ishf`, WRF `HFX` and AROME-HDF (flipped on read) are
positive upward and go through `_derive_surface_params!` unchanged.

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

| | ERA5(+CLD) | MERRA2 | GFS | AROME-HDF | AROME NetCDF | MesoNH | WRF |
|---|---|---|---|---|---|---|---|
| zlev source | ak/bk integr. | archived H | gh − orog (ramp below ground) | fixed heights | Ak/Bk interfaces integr. | Gal-Chen, w-level midpoints | PH+PHB interface midpoints |
| levdot | etadot | OMEGA | ω | — (0) | NH_dW | WT → mass levels | W → mass levels |
| tke | Hanna | Hanna | Hanna | Hanna | **native** | **native** | Hanna |
| ustar | zust / stress | log-law* | stress | stress | log-law* | log-law* | UST |
| Obukhov L | sshf (sign flipped) | Rib* | heat flux | sshf | Rib* | Rib* | HFX |
| pblh | blh | Rib profile* | Rib profile* | blh | TKE threshold* | TKE threshold* | PBLH |
| convection (umf/dmf/D) | ✓ full | umf only (edges → mid) | CAPE-gated diag.* | — | — | resolved | — |
| wet dep (3-D tpf) | ✓ | ✓ (edges → mid) | pseudo-3-D* | — (tp accum) | — | — | — |
| land-sea mask | lsm | — | lsm | — | — | — | LANDMASK |
| roughness z0 | — | — | fsr | — | — | — | — |

`*` = derived fallback implemented in `metcalc_models.jl`
(`_surface_turbulence!` chain: stress → log-law ustar; heat flux → Rib-based
L; archived → Rib-profile pblh).

---

Wet deposition needs only `tpf` (and `cc`), convection only `umf`: since
2026-10-13 they are built separately, so a dataset with precipitation but no
mass flux (GFS) gets wet deposition, and a requested-but-impossible process warns.

## Checking a reader

Five scripts in `validation/` (2026-10-13) — run them on a new dataset, and
after any reader change:
- `check_met_readers.jl <config>`: layout, ranges, θ and ρ consistency,
  hydrostatic balance between levels (catches heights and pressures that belong
  to different levels), tke filled, u*/pblh ranges, Obukhov L negative at local
  noon and positive at local midnight over land, pblh deeper at noon, tpf/cc.
- `check_output_mapping.jl <config>`: one parcel at a time through both deposit
  paths, the mass must land in the output cell containing the parcel.
- `make_synthetic_wrfout.py` + `check_wrf_synthetic.jl`: the WRF reader against
  known truth, and one advection step against the wind on the sphere.  The
  synthetic W is diagnosed as WRF does (η̇ from mass continuity, then the
  geopotential equation), so it is consistent with the synthetic winds.
- `check_vertical_kinematics.jl <config> mesonh|arome|merra2` (60-s timestep,
  start at a met time): one step without turbulence; a parcel's height
  (MERRA-2: pressure) must change by w·dt (ω·dt) read from the file.
- `check_wellmixed_met.jl <config> HOURS NPARCELS [lon0 lon1 lat0 lat1]`: the
  well-mixed criterion on the real fields.  Parcels seeded uniformly in air
  mass must stay uniform (χ = parcel mass / air mass = 1) under the model's own
  transport, convection and deposition off; `WM_TURB=none` for advection alone.
  Reports χ by height band in an interior that outside air cannot reach, its
  block-to-block spread against Poisson noise, and the two lowest bands split by
  the relief of the block.  The only test here that does not reuse the model's
  own understanding of the fields — it measures whether the winds, the vertical
  velocity, the density and the levels belong together.

### Well-mixed results, 2026-10-13

χ in the lowest 100 m and the worst other band below 6 km, after the time shown
(start ≈ 1 at seeding; Poisson noise of a band mean 1–3 %).  `old` = the v1.1
archive (7790a03), `new` = after the 2026-10-13 commits.  Advection only unless
marked "turb." (Hanna turbulence on).

| data (grid) | test | old | new |
|---|---|---|---|
| ERA5 η̇ (1°) | 4 h, 15 min | 0.989 / 0.981 | 0.989 / 0.981 |
| ERA5, turb. | 4 h, 15 min | 0.995 / 0.985 | 1.008 / 0.985 |
| MERRA-2 (0.625°×0.5°) | 4 h, 15 min | **1.778** / 0.817 (spread 8×) | 1.012 / 0.983 |
| MERRA-2, turb. | 4 h, 15 min | **1.331** / 1.049 | 1.016 / 0.987 |
| GFS (0.5°) | 4 h, 15 min | 1.019 / 1.036 | **1.137** / 1.030 |
| GFS, turb. | 4 h, 15 min | **1.339** / 0.967 | 1.038 / 1.019 |
| WRF synthetic (30 km) | 2 h, 15 min | no density (NaN) | 1.032 / 0.994 |
| WRF synthetic, turb. | 2 h, 15 min | — | 1.012 / 1.010 |
| MesoNH (2 km) | 3 h, 1 min | 1.049 / 1.025 (spread 2× aloft) | 1.052 / 1.009 |
| AROME (1.3 km) | 0.5 h, 1 min | 1.115 / **0.845** | 1.045 / 1.030 |
| AROME | 0.5 h, 15 min | 1.053 / **0.840** | **1.155** / 1.047 |
| AROME, turb. | 0.5 h, 15 min | crash (segfault) | **1.135** / 1.045 |

Over STEEP terrain (block relief > 500 m) the lowest 100 m is not well mixed,
old or new:

| | flat blocks 0–100 / 100–300 m | steep blocks 0–100 / 100–300 m |
|---|---|---|
| MesoNH 3 h, 1 min, old | 1.030 / 1.003 | 1.337 / **0.326** |
| MesoNH 3 h, 1 min, new | 1.014 / 1.001 | **2.933** / 0.839 |
| AROME 0.5 h, 1 min, new | 1.010 / 1.018 | **1.258** / 1.104 |
| AROME 0.5 h, 15 min, new | 1.071 / 1.046 | **1.574** / 1.059 |

Open (not fixed):
1. **Steep km-scale terrain.**  The coordinate velocity is computed from the
   archived w with centred slopes of the level surfaces, while the parcels
   follow the bilinear surfaces between grid points; at 1–2 km over relief of
   several hundred metres per cell the two disagree and air accumulates in the
   lowest 100 m.  Not caused by the time step, the surface-layer wind
   (`surface_layer_wind = false` changes AROME by < 0.002) or the ground value
   of levdot (0 instead of ½·levdot[2]: identical on AROME).  Candidate fix: a
   coordinate velocity diagnosed from mass continuity on EASLEE's own grid, as
   η̇ is for ERA5.
2. **Time step at km scale.**  The horizontal step is Euler-forward with the
   wind at the start position (the CTL sub-steps refine only the vertical).
   AROME, 0–100 m after 0.5 h: 1.045 (1 min), 1.064 (3 min), 1.078 (5 min),
   1.155 (15 min).  The shipped AROME configs use 15 min.
3. **GFS near the ground** (+14 % in 4 h without turbulence, +4 % with): the
   isobaric levels plus the surface ramp give a very uneven index spacing near
   the ground (the first level above ground can be 1 hPa or 25 hPa up).
   Unchanged at 1-min steps and without the surface-layer wind; a
   mass-consistent ramp velocity made it worse (0.80 at 100–300 m).  Candidate
   fix: interpolate GFS onto hybrid σ–p levels and use the MERRA-2 path.

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
