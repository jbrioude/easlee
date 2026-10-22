# EASLEE transport distribution — changes

## v1.1.1 (2026-10-22) — verified meteorological interfaces

The interfaces to the meteorological models were verified on real data with tests
that do not depend on observations (hydrostatic and kinematic consistency, the
well-mixed criterion on real fields, the position of each field in the files),
following points reported by Yue Huang on the WRF interface.  The tests are in
`validation/`.

Changes to the ERA5 results of the model description paper:
- **Gridded output on regular met grids**: deposits were offset by
  (met spacing / output spacing − 1) cells.  ETEX (0.1° output on 0.25° ERA5, 72 h,
  n = 2629): RANK 2.76 → 2.26, entirely from station 8021, 20 km from the release
  (2.64 → 2.67 without it).  CAPTEX (aligned grids) and the station output
  (`[output.receptors]`) are unchanged.
- **Surface-layer stability without an archived heat flux** (ERA5_CAPTEX, ERA5_CONV,
  ERA5_FULL, ERA5_CAPE, used for ETEX and CAPTEX): the Obukhov length now comes from
  a bulk Richardson number; the surface layer was treated as neutral.  Together with
  the vertical-placement changes below: ETEX RANK 2.26 → 2.21, CAPTEX 2.43 → 2.47.
  The conclusions of the paper are unchanged.

Meteorological interfaces:
- One vertical layout for every reader (heights and vertical velocity on the half
  levels, the other fields on the full levels between them); the position of each
  field (staggered or centred u, v, w) is read from the file, checked and logged.
- Vertical velocity as the coordinate velocity on terrain-following and hybrid levels.
- WRF: thermodynamic fields, staggered winds, Lambert map factor, time-averaged
  mass-coupled fluxes (`listvariables/WRFv2_avgflx.vars`).
- MesoNH, AROME NetCDF: level and wind positions, valid times, grid size from the file.
- AROME ATMO-HDF: accumulated fluxes, isobaric stitch, vertical velocity from mass
  continuity (`#info vertical_velocity`).
- GFS: isobaric levels remapped to a hybrid σ–p coordinate; averaged fields converted
  to rates.  MERRA-2: edge and mid-layer fields.
- ERA5 and other GRIB data: accumulation windows and units decided once per met source
  (`#info sfcflux_units`); η̇ level placement decided from the data
  (`#info etadot_levels`); `sshf` sign (`ERA5.vars`).
- Output: `[output.grid] crs` (lonlat or native Lambert); edge cells of the GPU deposit.

Physics:
- Vertical placement consistent with the readers' layout in the turbulence profiles,
  the surface-layer wind, convection, dry and wet deposition, boundary inflow and
  domain fill.
- Wet deposition independent of the convective mass flux.
- Turbulence step kept inside the model column.

Other: download time limits for large files (MERRA-2); default kernel width of the
pairwise mixing on regular grids.

Configuration: on km-scale grids use synchronisation steps of a few minutes
(`docs/CONFIG.md`); the AROME configurations use 5 min.

## 2026-09-23 — pre-release review fixes

- **Web GUI security.** The server binds to `127.0.0.1` (override with `EASLEE_GUI_HOST`);
  every POST route refuses requests without `Content-Type: application/json` or with a
  non-local `Host`; the Julia console (`/api/eval`) is disabled unless the GUI is started
  with `EASLEE_GUI_REPL=1`. The Chemistry section is no longer offered by the GUI schema.
- **`[run].seed`** (`Random.seed!` at run start; 0 = unseeded). Reproducible bit for bit for
  the same thread count and advance backend (`docs/CONFIG.md`).
- **Chemistry leftovers refused at validation**: `task = "chemistry"`, `[chemistry] enabled
  = true` and `[chemistry.emissions]` now fail in `validate_runspec` with the sentence that
  names the missing module, instead of an `UndefVarError` during the run; the task enum
  offers `transport` and `validation` only.
- **Tests**: the three tests that write into the gitignored `results/` create it;
  `test_point_source_mixing.jl` and `test_two_puff_mixing.jl` now assert (conservation and
  totals within tolerance) instead of only printing; the three slow well-mixed tests
  (10–45 min each) run only with `EASLEE_SLOW_TESTS=1`, so the default suite is 35 tests
  in minutes; `.github/workflows/ci.yml` runs it on Julia 1.10 and 1.12.
- **Suite green again after the GPU port**: `test_v11_transport_ports.jl` expected the
  retired "horizontal_tke is CPU-only" refusal and the old unsupported-switch list (now it
  checks the explicit-backend refusal instead); `io/vmix_namelist.jl` stays this
  distribution's 3-tuple reader (the developer tree's 4-tuple Kz variant had been copied in
  and broke `test_vertical_mixing_c.jl`).
- **Metadata**: personal path removed from `configs/arome_hdf_france.toml`; README test count;
  `CITATION.cff` date and a `repository-code` placeholder to fill with the GitLab URL (and the
  Zenodo DOI after the first release); note in `docs/CONFIG.md` that checkpoints use Julia's
  `Serialization` and must be resumed with the Julia version that wrote them.

## 2026-09-23 — GPU/KA advance: the last threads-only switches ported, two kernel bugs fixed, switches passed as types

Ported from the developer tree (commits a0d6026, 7fd48d0, 4345524, 2bdf8f8, 77e0d3f,
8a12d39, 0e46fe6; chemistry left out).  Files: `physics/advance_ka.jl`,
`physics/advance_ka_gpu.jl`, `physics/turbulence_4d.jl`, `physics/turbulence_hanna.jl`,
`physics/advance.jl`, `met/interp.jl`, `core/operators.jl`,
`core/config.jl`, `types/config.jl`, `config/{runspec,schema,toml_io}.jl`, `build_run.jl`.

- **`[physics.turbulence] horizontal_hanna`** (new): Hanna (1982) horizontal Langevin
  (σ_u, σ_v, T_L from u*, h, L, w*) on Hanna-fill drivers, the in-PBL horizontal term
  FLEXPART uses; `hanna_h`/`hanna_wstar` in `turbulence_hanna.jl`; surface parameters
  sampled once per sync step at the entry position on both advances.
- **Nothing is threads-only any more.** On-the-fly Hanna (`turbulence_hanna_4d`, the
  device twin of `turbulence_hanna_v1`), `mesoscale`, `kh_ft`, `horizontal_tke`,
  `horizontal_hanna` and `[physics] surface_layer_wind` all run on the KA/GPU path; the
  mesoscale OU state joins the resident device set (12 → 15 arrays), the surface
  arrays (u*, h, L) are uploaded when any of them is armed.
- **No silent downgrade.** When `EASLEE_ADVANCE=ka_gpu|ka_cpu` is set explicitly and
  the backend cannot be honoured (CUDA not functional, 2-D/1-D wind mask), the run now
  ERRORS instead of running the threads advance under a GPU label (`_ADVANCE_REQUESTED`,
  `advance_dispatch!`).  The auto default (no `EASLEE_ADVANCE`) still falls back.  The
  two `build_run` refusals for horizontal_tke/horizontal_hanna on KA are retired.
- **Two KA kernel bugs fixed** (found by the well-mixed criterion, `test/test_wellmixed_*_4d.jl`):
  the top model layer was unreachable in `interpolation4D_irregzsingle` (loop to Nkmax,
  fallback Nkmax − 1e-3), and the KA Thomson drift term was not aligned with the CPU's
  (wshare·∂tke/∂z/σ_w with the SIGW_FLOOR guard).
- **Physics switches reach the kernel as `Val{Bool}`** (types, not values): the GPU
  compiler specialises on argument types, so runtime Bools compiled every branch into
  one kernel and cost +44 % on the advance with all switches OFF (ETEX 1 M: 116 → 167 s).
  Back to 118 s.  Any new kernel switch must go in as a type.
- Validation record: `docs/BENCHMARKS.md` §3b (all five switches ON, threads vs ka_gpu:
  2.3e-3 in mass, correlation 0.99993; per-switch pairs; well-mixed PBL and Thomson);
  benchmark matrix for the paper in §0d
  (`scripts/bench_matrix.sh`, `scripts/bench_matrix_table.py`, `scripts/ka_compare_grids.py`,
  `configs/bench_*`).  Tests added: `test_ka_switches_unit.jl`, `test_hanna_gpu_4d.jl`,
  `test_wellmixed_pbl_4d.jl`, `test_wellmixed_thomson_4d.jl`.
- `test/test_ka_parity.jl` updated to the new dispatch semantics (mesoscale and on-the-fly
  Hanna are supported on KA; an explicit `EASLEE_ADVANCE` with a 2-D wind mask errors; the
  auto default warns once): 33/33.  Note for this distribution: `NestRunConfig` has no
  `vmix_kz` field here (the developer tree's Kz port added one), so positional
  constructions of it differ between the trees by that argument.
- Integration check of the port on zion (2026-09-23): `EASLEE_ADVANCE=ka_gpu PREFETCH=1`
  on `configs/bench_gpu_etex_1M.toml` — `mode = :ka_gpu`, no override, advection 120.9 s,
  met read 19.6 s, TOTAL(run) 171.0 s (developer tree the same morning: 168 s).  Tests in
  this tree: `test_ka_switches_unit.jl` 5/5, `test_hanna_gpu_4d.jl` 12/12,
  `test_ka_parity.jl` 33/33, `test_wellmixed_pbl.jl` all pass.
- Not ported (chemistry, outside this distribution): sulfate dilution, HNO3-on-ice
  scavenging and isotherm, the coarse-ice sweep accumulator (`state.csweep`), NAT ceiling.

## 2026-09-13 — v1.1 (transport and met-reading fixes from the developer tree)

Version 1.0.0 → **1.1.0** (`CITATION.cff`).  Everything below is transport,
met reading or output; the developer tree's chemistry changes of the same
period are not part of this distribution.  Each hunk in the source carries a
`v1.1` comment naming the change.

Documentation and met-reader checks (2026-09-19)

- **`docs/MIXING_SCHEMES.md` removed.**  It framed `physics/mixing.jl` as the
  "old scheme", a port of a Fortran EASLEE, set against a "new scheme" on a TODO
  list — a distinction that does not exist in this distribution and was
  misleading.  The part worth keeping, the derivation of the pairwise exchange
  (cell list, Gaussian pair weights, second-moment renormalisation to `K_h`,
  conservative exchange, θ-limiter, coverage-based background relaxation), is
  now in the `[physics.mixing]` section of `docs/CONFIG.md`, where the
  configuration it describes already lived.  The routines themselves are
  unchanged.
- **Dangling references fixed.**  Four docs pointed at `test/<file>.jl` for
  tests that ship in `validation/`; three pointed at documents that are not
  part of this distribution.
- **`docs/ETEX_RESULTS.md` and `docs/CAPTEX_RESULTS.md` now carry a dated
  warning.**  Their tables were measured before the 2026-09-04 surface-flux
  correction (`sshf`, `ewss`, `nsss`), which feeds u*, the Obukhov length and
  therefore the PBL, so they do not describe the physics this version produces.
  The current-code measurements available are quoted alongside, with their
  configurations stated; a like-for-like re-run of the published ETEX config is
  still to be done.
- **FLEXPART, WRF and AROME are no longer called "legacy".**  `TIMESTEPS.md`
  carried the `LSYNCTIME`/`LOUTSTEP`/`LOUTAVER`/`LOUTSAMPLE`/`CTL`/`IFINE` table
  under a "Mapping (legacy → EASLEE TOML)" heading with a "Legacy" column; it is
  now a section of its own, "Coming from FLEXPART", framed as a cross-reference
  from FLEXPART's `COMMAND` to the TOML key that plays the same role.
  `OUTPUT_GRID.md` called FLEXPART's `OUTGRID` legacy, and `MET_DATA.md` and
  `ARCHITECTURE.md` called the WRF and AROME reader paths legacy when the point
  was only that v1.1 did not change them.  `ARCHITECTURE.md`'s "Legacy I/O"
  layer is EASLEE's own untyped Dict-based code, but it sat next to "reads
  FLEXPART-style input files" where it read as a judgement on FLEXPART; it is
  now "Untyped I/O", which is more descriptive anyway.  Where the word genuinely
  meant EASLEE's own earlier behaviour it now says so.
- **Sanity check on the converted vertical velocity** (`_levdot_sanity!` in
  `easlee/src/io/initialize.jl`).  The existing detector judges the magnitude of
  the ARCHIVED eta-dot; this one judges the RESULT, so an archive in neither
  convention — or a convention forced wrongly through `#info etadot_units` —
  is caught even when the input-side test cannot see it.  Runs once per met
  source, warns when max|w| over the interior levels leaves 1e-3…50 m/s.  On
  the ETEX files the correct reading gives max|w| = 0.457 m/s and the forced
  wrong one 5.06e-6 m/s, a factor 9e4.  See `docs/ETADOT_UNITS.md`;
  `test/test_etadot_units.jl` is now 20 checks.

Installation

- **No `Manifest.toml` is shipped** (engine or GUI), as since 2026-09-05: a
  manifest resolved on Julia 1.12 fails to precompile on 1.10, which is the
  install error reported from the HPC.  `Pkg.instantiate()` resolves against the
  `[compat]` bounds in the two `Project.toml` files instead.  Julia RECREATES a
  manifest as soon as anyone runs in the tree, so `make_tarball.sh` excludes
  them at build time and fails loudly if one slips in — deleting them once is
  not enough.  Build the distribution with `sh easlee_v1_1_shared/make_tarball.sh`.

Met reading and nested domains

- **AROME: `pblh` was identically zero for every run** (`io/initialize.jl`,
  `metaromefieldcalculation!`).  AROME archives neither `ustar` nor `pblh`
  (`listvariables/AROME.vars` leaves both blank) and the native-TKE path never
  ran the derivation the Hanna-fill paths do, so the field stayed at its
  allocation zeros: the `[tracers]` PBL-contact tag never fired (it floors at
  50 m), the trajectory CSV's `pblh` column was 0, and the `kh_ft`
  free-troposphere gate treated everything above 50 m as free troposphere.
  Now derived exactly as the other native-e driver (MesoNH) does it — log-law
  `ustar`, bulk-Richardson `ol`, TKE-threshold `pblh`.  The 3-D `tke` is
  untouched.  A Thomson run is unchanged by this and only gains the diagnostic.
- **Per-nest vertical-mixing closure** (`types/config.jl`, `core/config.jl`,
  `advance.jl`, `runtime.jl`).  Nests may legitimately declare different
  closures — AROME native `e` inside ERA5 Hanna fill (`tke = 0.5·σ_w²`) — and
  applying one nest's CBR constants to the other's field misread the archived
  `tke` by a factor `2·r_w`.  Each nest now carries its own
  (`native_e, r_w, c_k`) from its own variables namelist; the former "nests
  declare conflicting closures" error is an informational message.  A parcel
  transferred between nests has its stored turbulent velocities zeroed (they
  were drawn under the other nest's σ; the OU process re-establishes the right
  variance within one τ_L).

Transport

- **`[physics.turbulence] horizontal_tke`** (`turbulence.jl` `tke_partition_h`,
  `advance.jl`): TKE-driven horizontal Langevin velocity on native-e drivers.
  The Langevin route dispersed parcels vertically only (`uout = vout = 0`), so
  with `kh_ft = 0` and `[physics.mixing]` off there was **no horizontal
  dispersion at all**: on an AROME marine back-trajectory the cross-wind σ of a
  5000-parcel cloud went 0.19 → 0.15 → 0.15 → 0.18 km at ages 0/1/3/6 h.
  σ_u² = σ_v² = 2·e·(1−r_w)/2 from the same archived `e` as the vertical step,
  τ_Lh = K_h/σ_u²; refused on Hanna-fill met, with `[physics.mixing]` armed,
  without a turbulence scheme, or on the KA/GPU advance.  Off by default;
  off reproduces the previous output.  See CONFIG.md.
- `turbulence_hanna_v1`: guard against a rare non-finite / runaway normalised
  velocity (|w'| > 20 σ), reset to a fresh draw, counted in
  `_HANNA_WP_RESETS`, warned once.
- Release-slot allocation: the free-slot scan restarted from 1 for every
  release box every step (O(boxes × used slots)); a monotone cursor now
  resumes where the last scan stopped (share review B2c, previously applied
  in the developer tree only).  The per-step release count no longer
  materialises a BitVector over the whole run.
- Plume rise / plume phase: the Bolton liquid saturation fit (valid −35…+35 °C,
  1.7–1.9× e_si at 200–215 K) is replaced by Murphy & Koop over ice below 235 K
  and over liquid above (`ice_physics.jl` `esat_mk`, `esat_liq_mk`); the
  plume-phase ice saturation now calls the single `esat_ice` definition.

Output

- **NetCDF georeference of a projected output grid** (`types/output.jl`,
  `io/output.jl`, `io/initialize.jl`).  `initnetcdf` wrote the grid axes as
  `lon0 + dlon·(i−1)` under hard-coded `degrees_east` / `degrees_north`.  On
  `[output.grid] crs = "native"` over a Lambert nest those numbers are
  **projection metres**, so every external reader placed the grid at 0°N 0°E
  with 1300° spacing.  Such a grid is now written as CF `x`/`y` in metres plus
  2-D auxiliary `longitude`/`latitude` sampled from the nest's own LON/LAT;
  data variables carry `coordinates = "longitude latitude"`.  A lon/lat output
  grid is unchanged bit for bit.
- **`[output.trajectory] format = "demovor"`**: the DEMOVOR/AROME column
  set (no header, `%.5g`) plus an eighth `pblh_m` column.  Rows are sampled on
  the parcel's own nest: the writer used to read `ctx.nests[1]` for every
  parcel, so a parcel that had left an inner domain was placed with the inner
  grid's projection — silently, since the bogus indices usually still fell
  inside its bounds.

Tests

- `test/test_v11_transport_ports.jl` (new, synthetic data only): the
  horizontal-TKE partition and OU stationary variance, the Murphy–Koop laws,
  the free-slot cursor, config keys through the RunSpec and a TOML round
  trip, the NestRunConfig fields and CPU-only switch, the trajectory header by
  format.
- Two test-side positional `NestRunConfig(...)` constructors brought to the
  full field order (`test/test_ka_parity.jl`; `validation/test_polar_stereo.jl`,
  which was a 12-argument form predating `kh_ft`/`mo_surface` and was already
  broken before this release).

Validation (2026-09-13, zion, Julia 1.12.6)

- Test suite: **39/39** (`test/runtests.jl`, 8 threads).
- ETEX regression: `configs/output_etex.toml` (400 k parcels, 72 h, 16 threads)
  scored with `validation/check_etex.jl` reproduces the documented baseline
  **RANK 2.75** to the last digit (R 0.660 vs 0.651, FMS 55.8 vs 57.2 %, within
  run-to-run scatter); no Hanna-guard warning fired.  BENCHMARKS.md §0c.
- End-to-end trajectory smoke on real ERA5 (ETEX met, 6 h, 20 k parcels,
  `format = "demovor"`): 535 rows, every row 8 columns, `pblh_m` non-zero
  (mean 745 m, max 895 m), TKE sampled at the parcel.
- Not exercised here (needs an AROME/MesoNH archive): `horizontal_tke` on live
  met, the AROME `pblh` derivation on live met, and a Lambert `crs = "native"`
  NetCDF.  Their logic is the same code as the developer tree, where they were
  run on the TROPIVOC AROME 1.3 km batch; only their unit-level contracts are
  pinned in `test_v11_transport_ports.jl`.

## 2026-09-05

Installation

- No `Manifest.toml` is shipped any more (engine and GUI). `Pkg.instantiate()`
  resolves against the `[compat]` bounds now present in both `Project.toml`
  files, on Julia 1.10 or newer. Tested from a clean depot on Julia 1.12.6 and
  1.10.12. (A manifest resolved on 1.12 failed to precompile on 1.10, and the
  GUI manifest pinned a networking library whose precompile self-test used a
  certificate that expired on 2026-08-06.) See HOWTO §1.
- Every run prints one startup line saying whether the GPU is actually used
  (`GPU backend enabled: <device>; CUDA initialized`) or why not. Being inside
  a GPU allocation is not using it; the line settles the question. HOWTO §1.3.
- Unknown or misspelled configuration keys are now an error listing the key and
  the nearest known one (`strict_keys = false` downgrades it to a warning).
- `test/runtests.jl` passes the thread count to each test; 38 tests.

Physics and numerics (transport results change only where noted)

- `EASLEE_ADVANCE=ka_cpu|ka_gpu`: the vertical mean wind was applied
  (ctl+1)/2 times per step, and the default `hanna` turbulence scheme ran as no
  turbulence. Both fixed; the KA path now reproduces the threads path (8-h
  CAPTEX release-1 check: mean plume height 511 m vs 510 m, layer mass
  fractions within 0.003; before the fix 42 % of the mass sat below 100 m and
  a tail reached 4 km). CPU-only switches (`surface_layer_wind`,
  `turbulence.mesoscale`, `kh_ft`, on-the-fly Hanna) now force the threads
  advance with a warning instead of being dropped silently. Met fields are
  uploaded to the GPU only when they change (was: every step).
- Threads advance: the sub-step loop is exact (it ran ctl+1 sub-steps for
  ~20 % of (Δt, ctl) pairs; the shipped configs were not among them); the
  staggered vertical read at the model top no longer touches level nz+1;
  `physics.ctl` must be an integer; the Monin-Obukhov correction is skipped
  above the lowest level (same result, less work).
- Aerosols: `dry_deposition` / `wet_deposition` and `[physics.deposition]` are
  honoured for aerosol species (deposition used to run regardless); gas wet
  deposition now also runs in mixed gas + aerosol runs (was silently off); two
  aerosol species with different size parameters are refused (one bin set
  was shared silently). `species.molar_mass` and `species.refractive_index`
  now reach the physics. `[species.settling]` is the documented table name
  (`[species.fall]` still accepted).
- Met reader: surface-flux units (sshf, ewss/nsss) are detected from the field
  magnitude instead of assuming hourly accumulations (affects only the `ERA5`
  name map, no shipped config); end-of-archive holds the last frame on the
  prefetch path (was a crash); a run start outside the archive names the
  dates; the WRF time-slice lookup no longer loops over a date-string length;
  latitude flip follows `jScansPositively`; incomplete GRIB files warn per
  variable; GRIB decode no longer makes two Float64 copies per message.
- Output: culled domain-fill parcels are not deposited one extra step.
- Example config `captex_rel1.toml` now points at the 1983 dates (CAPTEX-83).

Every item above was re-read in the source and, where numeric, reproduced
with a small script before being fixed.
in the developer tree.

## 2026-09-24 — validation drivers: forward/backward consistency

- `validation/validate_captex_run.jl`: the backward convolution now restricts each
  sensitivity frame to its overlap with the release window (Seibert & Frank 2004); the
  previous whole-run integral × emission rate overestimated the backward concentration
  ~3× at CAPTEX release 2. The time loop starts at the run-start step (the first release
  step was skipped, 1/13 of the mass), and stations are sampled at the cell centres
  (the single-frame NetCDF axis holds the cell edges). No model change.
- `docs/CAPTEX_RESULTS.md`: backward section rewritten accordingly.

## 2026-09-24 — output axis at cell centres; adjoint convection in backward runs

- `easlee/src/io/output.jl`, `easlee/src/io/initialize.jl`: the single-frame NetCDF axes
  (`longitude`/`latitude`, `lon2d`/`lat2d`) are now the CELL CENTRES (`lon0 + dlon·(i−½)`),
  as in the time-framed writer; before they were the west/south edges and bilinear readers
  sampled half a cell off. Readers of older files: add half a cell.
- `easlee/src/physics/convection.jl`: backward runs apply the adjoint (mass-weighted
  transpose) of the mass-flux kernel instead of the forward kernel; the per-parcel kernel is
  now `_convect_column(rng, k0, …, backward)`; `test/test_convection_adjoint.jl` (5 checks).
- `docs/OUTPUT_GRID.md` updated.

## 2026-09-25 — the download tools are Julia, so Julia is the only prerequisite

- `download/` is now Julia: `gfs.jl`, `merra2.jl`, `era5_archive.jl`,
  `easlee_available.jl`, `repack_grib.jl`, sharing `cli.jl` (flag parsing and a
  throttled progress printer).  The Python versions are gone.  Nothing outside
  Julia is needed any more, which matters most on Windows, where Python is
  usually absent.  Verified against the Python output: `easlee_available.jl`
  writes byte-identical `AVAILABLE` and `AVAILABLE.files`, `repack_grib.jl`
  produces byte-identical GRIB, `gfs.jl` builds identical URLs, and
  `era5_archive.jl` round-trips a 2.1 GB ETEX tar to an identical index.
- `repack_grib.jl` takes `grib_set`/`grib_count` from `eccodes_jll`, the
  artifact GRIB.jl already uses, so a system ecCodes install is no longer
  needed.  `era5_archive.jl` unpacks with the `Tar` standard library and
  `p7zip_jll`, and lifts a wrapping directory out of the archive so both tar
  layouts land in the same place.
- The web GUI's Data tab spawns these scripts with the easlee project instead of
  `python3`; `EASLEE_PYTHON` is retired.
- `easlee/Project.toml` declares `Downloads`, `Tar`, `eccodes_jll` and
  `p7zip_jll`, all of which were already in the tree indirectly.
- Analysis helpers now exist in both languages, so post-processing can be done
  in either: `scripts/ka_compare_grids.jl`, `scripts/bench_matrix_table.jl` and
  `tools/verify_nc.jl` join their `.py` twins.

## 2026-09-24 — data download tools, pre-generated AVAILABLE index, GUI Data tab, AROME configs

- `download/`: `gfs.jl` (NOMADS + NOAA/AWS archive, analysis+forecast or date range),
  `merra2.jl` (GES DISC, Earthdata token/netrc), `era5_archive.jl` (prepared ETEX/CAPTEX
  archives, links in `links.toml`), `easlee_available.jl` (index writer).  Julia, run
  with the easlee project.  Every script writes `<dir>/AVAILABLE` + `AVAILABLE.files`.
- Model: `[meteorology].available = "auto"` first uses `<path>/AVAILABLE` when it still
  describes the directory (sidecar file list + sizes, or no newer file), else scans and
  writes the index there for the next run; `"scan"` forces a rescan
  (`io/met_scan.jl`, `build_run.jl _resolve_available`, `test/test_available_cache.jl`).
  Measured: 0.3 s instead of a 21 s scan on a 26-file ERA5 day.
- GUI: **Data** tab (`POST /api/met/download`, `/status`, `/cancel`) runs the scripts and
  can set `[[meteorology]]` of the Config tab to the downloaded directory.
- Configs: `arome_meteofrance.toml` (native Météo-France NetCDF, model `AROME`) added
  beside `arome_hdf_france.toml` (ATMO Hauts-de-France GRIB2, model `AROME_HDF`).

## 2026-09-25 — mineral-dust source (Kok et al. 2014), sea-spray level fix

- `[[release]] kind = "dust"` (`easlee/src/physics/dust_kok.jl`): Kok et al. (2014) vertical
  dust flux over the land cells of the box from the met friction velocity, Shao & Lu (2000)
  threshold with the Fécan et al. (1999) soil-moisture increase, constant or mapped clay and
  bare fractions (`[release.dust]`); species `size.distribution = "kok"` = the Kok (2011)
  emitted size distribution. Checkpointed; run-end "dust (Kok 2014) budget" line.
  `docs/DUST_KOK.md`, `configs/example_dust_kok.toml`, `test/test_dust_kok.jl` (19 checks).
- Sea-spray source: spawned parcels were placed at level index Nk − 0.2, which is the model
  TOP (level 1 is the surface after the reader's level shift); now the lowest layer (1.2),
  with the jittered position kept inside the met grid. Any sea-spray run made before this
  date released its spray at the top of the atmosphere.

## 2026-09-26 — AOD optics validation without transport (validation/check_aod_optics_cams.jl)

- The CAMS aerosol-class optics (table, single-mode and lognormal mass-extinction
  coefficients, κ-Köhler/Mie wet enhancement) moved to `easlee/src/physics/cams_optics.jl`
  (transport side, included by Operators).
- `validation/check_aod_optics_cams.jl`: reads one CAMS analysis (aermr01–20 + aerext532) and
  the ERA5 file of the same hour (t, q, lnsp), computes the 550-nm column AOD with EASLEE's
  optics at every grid point and compares it with the CAMS 532-nm extinction column — no
  parcels, no chemistry, ~100 s. 20 Sept 2025 15 UTC, Africa / SW Indian Ocean window:
  lognormal σ 2, 3 radii → ratio 1.11, r 0.75 (PASS: 0.7 ≤ ratio ≤ 1.4, r > 0.7);
  `--dry --size single` → 0.55, r 0.96 (the pre-growth optics).

## 2026-09-27 — sea-spray source: single OSSA flux evaluation per cell

- `easlee/src/physics/seaspray.jl`: the flux is the published Ovadnevaite et al. (2014) function
  evaluated once per ocean cell and step; a development-only option and its config key were
  removed. Spawned parcels now sit in the surface layer (fix of 2026-09-25); the test checks it.
- `test/test_source_mass_quantum.jl`: the spawned parcel mass plus the unspawned remainder
  equals the flux integrated independently over the cells and steps, for both sources.

## 2026-10-02 — vertical velocity to level rate was four times too small (WRF, AROME)

Reported by **Yue Huang**, who spotted it in the WRF reader.  It was in the
AROME reader as well.

A vertical velocity in m/s becomes a level-index rate by dividing by the LOCAL
level spacing, which for an interior level is the centred half span,
`(z[k+1] - z[k-1]) / 2`.  Both readers divided by the full two-level span and
then by two:

```julia
W ./ (zlev[:,:,3:Nk] - zlev[:,:,1:Nk-2]) / 2.f0     # = W / (2·span)
lev[:,:,k+1] / dz2 / 2.0f0                          # = w / (2·dz2)
```

`/` is left-associative in Julia, so the 2 divided the result rather than the
spacing, and the rate came out **four times too small** at every interior level.
The one-sided formulas these same readers use at k = 1 and k = Nk were correct,
so the interior disagreed with its own boundaries by a factor of four.  The
third path that does this conversion (`metcalc_models.jl`, the `WT` reader) had
it right all along, which is what made the other two diagnosable.

Affected: `WRFv2.vars` (maps `W`) and `AROME.vars` (maps `NH_dW`).  ERA5, GFS,
MERRA-2 and the WRF archived-etadot branch never took these paths.  WRF and
AROME runs through them were under-ventilated vertically in the interior.

The conversion now lives in one place, `levrate_from_w` in `initialize.jl`, and
both readers call it.  The duplication is what let the two copies drift from the
correct one.  Note that a centred DERIVATIVE on the same stencil (drho/dz,
dsigw/dz) divides by the FULL span and must not carry the factor; those sit on
the adjacent lines and are unchanged.

New test `test/test_levdot_spacing.jl` (11 checks): uniform and stretched
columns, agreement between an interior level and the one-sided boundaries, a
regression guard asserting the old form was wrong by exactly four, and a check
that the centred derivative still divides by the full span.  Against the old
code 10 of its 11 checks fail.
