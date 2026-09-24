# EASLEE — NetCDF output grid

Status: 2026-06-13.  Code: `easlee/src/output/conc.jl` (accumulation),
`easlee/src/io/output.jl` (NetCDF init + write).  Validation:
`validation/test_backward_srr.jl`.

## TOML representation (`[output.grid]`)

FLEXPART's `OUTGRID` file is ambiguous about units: `OUTLON0`/`DXOUT` mean **degrees** for
lat/lon mets (ERA5 `0.25`) but **metres** for projected mets (WRF `24000`), with
the choice inferred from the first nest's projection, and the vertical is always
height-AGL. The TOML makes **both axes explicit**.

**Horizontal** — `crs` selects the coordinate system (and thus the units):
```toml
[output.grid]
crs = "lonlat"          # "lonlat" (deg) | "native" (metres of a nest's projection)
# nest = 1              # native only: which nest's projection defines x/y
# two interchangeable field styles (use whichever is natural):
origin  = [60.0, -55.0]
spacing = [0.25, 0.25]
size    = [120, 100]
# ── or the scalar style ──
# lon0 = 60.0 ; lat0 = -55.0 ; dx = 0.25 ; dy = 0.25 ; nx = 120 ; ny = 100
```
For `crs="native"`, `origin`/`lon0,lat0` are `x0,y0` and `spacing`/`dx,dy` are
metres in nest `nest`'s projection.

**Vertical** — `kind` + `levels` (units follow `kind`):
```toml
  [output.grid.vertical]
  kind   = "height_agl"   # height_agl (m) | pressure (hPa) | model_level (idx) | theta (K)
  levels = [0, 500, 1000, 2000, 3000, 6000, 10000]
```
A bare `heights = [...]` is still accepted and means `kind = "height_agl"`.

**Status:** the config is fully parsed/validated (both field styles, all four
`vkind`s, `crs`, `nest`). Runtime binning currently still uses **height-AGL** and
infers the horizontal CRS from the met projection, so the default cases
(`lonlat→lonlat`, `native→native`) work today; a build-time `@warn` fires for a
non-`height_agl` `vkind`. Wiring the four vertical binners and the
`projected→lonlat` horizontal case in `conc.jl`/`initoutgrid` is the
next step (forward-compat keys `CRS`/`OUTNEST`/`OUTVKIND` are already emitted).

## What the grid is

A 5-D accumulator `[Niout, Njout, Nkout, Nc, Nage]` defined by the
`[output.grid]` config (regular lon/lat box + a list of layer-edge heights).  Two quantities are
accumulated every step by `conccalc!` (deposit happens BEFORE advance, so the
position recorded is the cell the parcel occupies at the step's start):

| Accumulator | Per deposit | Forward meaning | Backward meaning |
|---|---|---|---|
| `massarray` | `Σ m_p·w` | mass → **concentration** | residence → **sensitivity** |
| `concarray` | `Σ (m_p/ρ)·w` | volume → **mixing ratio** | → **mixing-ratio sensitivity** |

`w = 1` forward; `w = Δt/M_released` backward (see below).

## Output kinds — FLEXPART `IOUT` semantics

| IOUT | Written |
|---|---|
| 1, 9 | concentration only |
| 2, 10 | mass mixing ratio only |
| 3, 11 | both |

- **Concentration** `spec{NNN}_conc = Σ m_p /(V·n_samples)`  [kg m⁻³].
- **Mixing ratio** `spec{NNN}_mr = Σ (m_p/ρ_p)/(V·n_samples)`  [kg kg⁻¹]
  (mass of tracer per mass of air: `m/ρ` is the tracer's air-equivalent
  volume, `/V` the air volume fraction).

`n_samples` is the number of output sub-samples in the averaging interval.

## Time-framed output (engine, 2026-07-03)

Forward `run_spec!` runs that fit at least one full `output.interval` write a
TIME-RESOLVED `grid_<name>.nc` with FLEXPART `LOUTSTEP/LOUTAVER/LOUTSAMPLE`
semantics: parcels are deposited only within the trailing `output.average`
window of each `output.interval`, at the `output.sample` cadence, and each
frame is frozen at the interval end:

- `conc(lon, lat, height, time)` = Σ m_p /(V·n_samples)  [kg m⁻³]
  (summed over release boxes), when `quantity` includes concentration;
- `mr(lon, lat, height, time)`   = Σ (m_p/ρ_p)/(V·n_samples)  [kg kg⁻¹],
  when `quantity` includes mixing ratio;
- `time` = seconds since run start at each averaging-window END; the
  `interval_s/average_s/sample_s` used are global attributes.

Constraints: `sample`, `average`, `interval` must be integer multiples of
`run.timestep`, `sample ≤ average ≤ interval`.  Backward runs (time-integrated
sensitivity) and runs shorter than one interval keep the single-frame
writer below.  This is what makes transport runs directly comparable with
time-resolved observations — see `configs/output_etex.toml` +
`validation/check_etex.jl` (ETEX self-check) and ETEX_RESULTS.md.

## Per species and per release box — `[output].group`

- `SPECNUM_REL` (from the release configuration) sets the FLEXPART species
  number of each release box (single value → all boxes; a list → one per box).  Parsed into
  `ReleaseSchedule.species` and recorded on `OutputConfig.boxspecies`.
- **What the accumulator's 4th dimension is** is chosen by `[output] group`
  (COMMAND `OUTPUTGROUP`), because the grids cost 12 B/cell **per column**
  (Float32 conc + Float64 mass ≈ 18 MB/column on the HT grid — per-box columns
  OOM-killed a 3400-box file release at 70 GB):
  - `group = "species"` (**default for forward runs**): one column per
    species; `OutputConfig.outclass` maps each release box to its species
    column at deposit time (`conc.jl` / `conc_ka.jl`), so grid memory is
    independent of release-box count.  `conc_box`/`mr_box` are **not**
    written (they would duplicate the `spec{NNN}_*` variables).
  - `group = "release"` (alias `"box"`; **default for backward runs**, where
    sensitivities must be attributable to individual boxes): identity
    mapping, bit-identical to the behaviour before `group` existed, same
    variable set.
- With `IOUTPUTFOREACHRELEASE = 1` the writer produces:
  - **per species**: `spec{NNN}_conc` / `spec{NNN}_mr`, summed over the
    columns that share species NNN (`OutputConfig.classspecies`);
  - **per box** (release grouping only): `conc_box` / `mr_box` over the
    `numpoint` dimension, plus `release_species(numpoint)` mapping each
    column to its species.
- With `IOUTPUTFOREACHRELEASE = 0` (`Nc = 1`) everything aggregates into a
  single class, written as `spec{NNN}_*` if one species, else `specall_*`.
- Backward normalisation (`relmass`) is **per column**: under species grouping
  the released mass is summed over the boxes of each species.
- The original `conc` variable (units label `ng/m3`,
  `concarray·1e12/n_samples/V`) is always written for backward compatibility
  with the `validate_*` scripts (per-column when `Nc > 1`).
- Tests: `test/test_output_grouping.jl`.

## Backward runs — source-receptor sensitivity

Set `LDIRECT = -1`.  The loop runs from `enddate` to `startdate`; parcels are
released at the *receptor* and advected upwind.  Each deposit carries the
residence weight `w = Δt/M` (Seibert & Frank 2004), so the grids become
**time-integrated** (no `n_samples` averaging):

```
sens_ijk    = (Δt/M) Σ_t Σ_{p∈cell ijk} m_p      / V_ijk        [s m⁻³]
sens_mr_ijk = (Δt/M) Σ_t Σ_{p∈cell ijk} m_p/ρ_p  / V_ijk        [s kg⁻¹]
```

Variables: `spec{NNN}_sens`, `spec{NNN}_sens_mr` (and `sens_box`,
`sens_mr_box` per box).  Units attributes are written on each variable.

**Interpretation.** `sens_ijk` is the sensitivity of the receptor
concentration to an emission in cell `ijk`: multiply by an emission flux
density and sum to predict the receptor signal.  Because every released
parcel carries `m_p = M/N`, the grid satisfies the residence-time identity

```
Σ_ijk sens_ijk · V_ijk = (1/N) Σ_p τ_p
```

— the **mean time the released particles spend inside the output domain**.
This identity is the basis of `test_backward_srr.jl`: the grid sum is compared
against the residence time tracked independently from the parcel positions
(non-circular), and the two agree to round-off.  Parcels that leave the
(finite) output box stop contributing — correctly, since the sensitivity is
only defined where the grid is.

## Validation summary

Forward (2 boxes / 2 species / IOUT=3):
- mass closure `Σ conc·V·n = Σ massarray` to machine precision per species;
- mixing-ratio closure `Σ mr·V·n = Σ(m/ρ)` likewise;
- per-box variable ≡ matching species aggregate; `release_species` correct;
- the original `conc` variable present.

Backward (`test_backward_srr.jl`):
- grid `Σ sens·V` equals the independently tracked mean in-domain residence
  time;
- mixing-ratio sensitivity equals `Σ sens·V` divided through by air density
  (cell-by-cell), i.e. `sens_mr = sens/ρ`;
- parcels move UPWIND (footprint extends opposite the mean flow).

## Axis registration — cell CENTRES (2026-09-24)

Both NetCDF writers store the horizontal axes at the **cell centres**:
`longitude[i] = lon0 + dlon·(i − ½)` in the single-frame writer (`io/output.jl`,
`initnetcdf!`, variables `longitude`/`latitude`, also `lon2d`/`lat2d` on a native CRS) and
`lon[i] = lon0 + (i − ½)·dlon` in the time-framed writer (`output/chemconc.jl`, `lon`/`lat`).
`lon0`/`lat0` in the config remain the grid's south-west corner; cell i spans
`[lon0 + dlon·(i−1), lon0 + dlon·i]`. Until 2026-09-24 the single-frame writer stored the
west/south EDGES, and every bilinear reader that took the axis values as nodes sampled half
a cell to the south-west (the CAPTEX validation, paper Table 3). Files written before that
date by `initnetcdf!` need `+ dlon/2` before sampling; the time-framed files never did.

## Backward runs and convection (2026-09-24)

The mass-flux scheme applies the ADJOINT of its forward kernel in backward runs
(`physics/convection.jl`: reversed subsidence, then reversed conveyor — the mass-weighted
transpose, `test/test_convection_adjoint.jl`). The `cape` scheme mixes toward uniform and is
its own adjoint.
