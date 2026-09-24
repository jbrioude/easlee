# EASLEE — install and run

EASLEE is a Lagrangian particle dispersion model (LPDM) in the FLEXPART
lineage, written in Julia. This distribution is the **transport engine**: met
ingestion, trajectories, turbulence, convection, deposition, aerosols,
inter-trajectory mixing and the gridded/receptor output estimators. It does not
include the reactive-chemistry module.

Everything is driven by **one TOML file per run**. There is one entry point —
you never edit a script.

Meteorology supported here: **WRF**, ERA5, NCEP GFS, MERRA-2, AROME (Météo-France NetCDF: `configs/arome_meteofrance.toml`; ATMO Hauts-de-France GRIB2: `configs/arome_hdf_france.toml`), MesoNH.
Configs to start from live in [`configs/`](configs/):

| file | what it is |
|---|---|
| `configs/wrf_single.toml` | WRF only, one domain |
| `configs/wrf_era5_nest.toml` | WRF as the inner nest, ERA5 outside it |
| `configs/wrf_gfs_nest.toml` | WRF as the inner nest, NCEP GFS 0.25° outside it |
| `configs/volcanic_mixing.toml` + `volcanic_nomix.toml` | a single-factor twin isolating what `[physics.mixing]` does |
| `configs/era5_conv.toml`, `gfs.toml`, `merra2.toml`, `arome_hdf_france.toml`, `mesonh.toml` | one per other supported driver |

---

## 1. Install

### 1.1 Linux / macOS

**Julia** (1.10 or newer; 1.12 is what this tree is developed on):

```bash
curl -fsSL https://install.julialang.org | sh     # installs `juliaup` + `julia`
exec $SHELL -l                                    # pick up the new PATH
julia --version
```

Then, from the root of this directory, install the model's dependencies once:

```bash
cd /path/to/easlee_share
julia --project=easlee -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Optional — only if you want the web interface:

```bash
julia --project=webgui -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
```

First install downloads and compiles a few hundred MB (mostly `NetCDF`, `GRIB`
and `CairoMakie` for the GUI); expect 5–20 minutes. Nothing is installed
system-wide: `--project=easlee` keeps everything in this directory's
environment plus `~/.julia`.

**No `Manifest.toml` is shipped.** `Pkg.instantiate()` resolves the newest
package versions compatible with the `[compat]` bounds in
`easlee/Project.toml` / `webgui/Project.toml` on *your* Julia, so a stale pin
cannot break the install and Julia 1.10, 1.11 and 1.12 all resolve their own
set (a manifest resolved on 1.12 fails to precompile on 1.10). The
combination this release was tested with (Julia 1.12.6 and 1.10.12,
2026-09-05) was: Atomix 1.1.3, GRIB 0.4.0, Interpolations 0.16.3,
KernelAbstractions 0.9.42, NetCDF 0.12.2, StructArrays 0.7.3; GUI: CairoMakie
0.15.11, HTTP 2.1.1 (Reseau 1.4.1), JSON3 1.14.3, PackageCompiler 2.4.1. If a
later resolve ever misbehaves, pin these with `Pkg.add(name=..., version=...)`.

> **`easlee/` is a Julia environment, not a package.** You never write
> `using EASLEE`. The run scripts `include` the engine source directly;
> `--project=easlee` only supplies the dependencies.

### 1.2 Windows

Double-click **`install_easlee.bat`**. It installs Julia (via `winget`) if
needed, then the model and GUI libraries. If it just installed Julia, close the
window and double-click the file once more so Windows picks up the new `PATH`.

Optionally double-click **`build_sysimage.bat`** afterwards: it compiles the
model + GUI into `easlee_sysimage.dll` so the interface starts in seconds
rather than a minute. It takes 10–25 minutes and runs once.

Everything below works the same on Windows; use `\` in paths and
`scripts\easlee.jl`. See [`docs/WINDOWS.md`](docs/WINDOWS.md) for the
click-only route.

**If the GUI step fails with `tls handshake failed: tls: certificate has
expired`** (in `Reseau`, while precompiling `HTTP`): versions of the
networking library Reseau before 1.3.5 (2026-08-07) ran a precompile self-test
against a built-in certificate that expired on 2026-08-06. A fresh resolve
picks a fixed version; if you still see it your local package registry is
stale, so run once

```bash
julia -e 'import Pkg; Pkg.Registry.update()'
julia --project=webgui -e 'import Pkg; Pkg.update("Reseau"); Pkg.precompile()'
```

and re-run the installer. Nothing in the model itself is affected: the
transport engine (`--project=easlee`) does not use this library.

### 1.3 GPU (optional)

CPU-only is the default and needs nothing extra. To use an NVIDIA GPU for the
advection kernel, install CUDA **once in your default Julia environment** (not
in `easlee`), then set an environment variable:

```bash
julia -e 'using Pkg; Pkg.add("CUDA")'
EASLEE_ADVANCE=ka_gpu julia --project=easlee scripts/easlee.jl configs/wrf_single.toml
```

If CUDA is missing or non-functional the run falls back to the CPU path with a
warning — it never fails because of the GPU.

**Check that the GPU is really used.** Being inside a GPU allocation (on an
HPC, `nvidia-smi` lists the card) and using it are two different things: the
GPU is used only if `EASLEE_ADVANCE=ka_gpu` reaches the Julia process. Every
run prints one status line right after the run header, so look for it in the
log:

```
[ Info: GPU backend enabled: NVIDIA A30, 24 GiB (CUDA runtime 12.x, driver 12.y); CUDA initialized — particle advance on GPU, 16 CPU threads
```

or, when nothing was requested,

```
[ Info: GPU: not requested — CPU run on 16 threads (set EASLEE_ADVANCE=ka_gpu for the GPU particle advance)
```

A `Warning:` variant says why a requested GPU could not be used (CUDA not
installed, or no usable device/driver). In a batch script, `export` the
variable before the `julia` line (or write it on the same line as above) — a
variable set only in your login shell does not reach a Slurm/PBS job. Once the
line says "CUDA initialized", `nvidia-smi` shows a Julia process holding a few
hundred MiB on the card.

---

## 2. Check the install

```bash
julia --project=easlee test/runtests.jl          # 38 self-contained tests,
                                                 # synthetic data only — no met
                                                 # archive needed, so this runs
                                                 # anywhere
julia --project=easlee test/runtests.jl mixing   # or a subset by name
```

To validate a config without running it (this parses **and** validates every
field, and reports *all* problems at once):

```bash
julia --project=easlee -e 'include("easlee/src/config/EasleeConfig.jl");
                           using .EasleeConfig; read_runspec("configs/wrf_single.toml")'
```

---

## 3. Run from the command line

One entry point, one config:

```bash
julia --project=easlee scripts/easlee.jl configs/wrf_single.toml
```

The config's top-level `task` selects what happens (`"transport"` — the default
— or `"validation"` for the CAPTEX/ETEX DATEM tracer experiments), and
`[run].threads` sets the worker count (the launcher passes `-t` for you).

Output lands in `[output].path`:

| file | contents |
|---|---|
| `grid_<name>.nc` | the gridded concentration / mixing-ratio field, one time frame per `[output].interval` |
| `receptors_<name>.csv` | gridless time series at station coordinates, if `[output.receptors].file` is set |
| `traj_<name>.csv` | per-parcel trajectories, if `[output.trajectory].enabled` |

**Lower-level runner.** `scripts/easlee.jl` is a thin dispatcher; if it ever
gets in the way, run the engine directly — same config, same result:

```bash
julia -t 8 --project=easlee scripts/run_config.jl configs/wrf_single.toml
```

**Long runs.** Julia buffers stdout heavily, so an empty log does *not* mean a
dead run — check progress from the number of time frames in the output NetCDF,
or from the process CPU time. On Linux, prefer a resource-capped unit over a
bare `nohup`:

```bash
systemd-run --user --unit=easlee-wrf --collect \
  -p MemoryHigh=20G -p MemoryMax=24G \
  --working-directory=/path/to/easlee_share \
  sh -c 'julia -t 8 --project=easlee scripts/run_config.jl configs/wrf_single.toml \
         > /tmp/wrf.log 2>&1'
```

---

## 3b. Get meteorological data

No archive at hand?  `download/` fetches public data and writes the `AVAILABLE` index the
model needs, so the first run starts immediately:

```bash
python3 download/gfs.py --analysis 2026-06-03T00 --hours 48 --out ./data/GFS        # NCEP GFS, no account
python3 download/merra2.py --start 2022-01-18 --end 2022-01-20 --out ./data/MERRA2  # MERRA-2, Earthdata account
python3 download/era5_archive.py --dataset etex --out ./data/ERA5_ETEX              # ERA5 of the ETEX-1 case
```

then point `[[meteorology]]` at the directory (`model`, `path`, `variables`, `available =
"auto"`); see `download/README.md`.  The GUI's **Data** tab does the same with a form and a
**Use in config** button.  `available = "auto"` reuses `<path>/AVAILABLE` whenever it still
matches the directory (files added/removed/resized trigger a rescan) and writes it after a
scan; `available = "scan"` forces the scan.

## 4. Run from the web interface

```bash
./gui.sh                    # Linux/macOS  (= julia --threads=auto --project=webgui webgui/easlee_app.jl)
```

On Windows: double-click **`EASLEE.vbs`** (no console window) or
**`run_gui.bat`** (console visible; close it to stop).

Then open <http://localhost:8080>. Tabs:

- **Config** — an auto-generated form built from the model's own field schema,
  so it always matches the code: every option, its unit, its default and its
  allowed values (hover a row for the description). Edit → **Save** → **RUN**.
  **Load TOML** opens any `configs/*.toml` for editing; **Download TOML** saves
  your edits back to a file you can re-run from the command line.
- **Figures** — column tracer, vertical cross-section, parcel positions and raw
  met fields, refreshed live as the run steps.
- **Console** — a Julia REPL attached to the running model.
- **Data** — download NCEP GFS (one analysis + forecast, or a date range of cycles and
  steps), MERRA-2 (date range, Earthdata token) or a prepared ERA5 archive into a
  directory, with the `AVAILABLE` index written; **Use in config** points the Config tab's
  `[[meteorology]]` at it.  Needs `python3` (`EASLEE_PYTHON` to override).

The GUI and the CLI drive the *same* engine from the *same* TOML — anything you
build in the form runs identically with `scripts/easlee.jl`.

---

## 5. Pointing EASLEE at WRF data

```toml
[[meteorology]]
model     = "WRF"
path      = "/data/WRF/my_run/"   # directory of wrfout_* files
variables = "WRFv2"               # → listvariables/WRFv2.vars
available = "auto"                # scan `path`, index every valid time
crop      = "off"
```

- `path` is a **directory**; `available = "auto"` reads each `wrfout_*` file and
  builds the time index. Files with several valid times each are handled.
- `variables` names a mapping file in [`listvariables/`](listvariables/), which
  translates the archive's variable names to EASLEE's field names. Two exist for
  WRF: **`WRFv2`** (recommended: `U`/`V`/`W` staggered winds, `PH`+`PHB` →
  geopotential height, `PSFC`, `PBLH`, `UST`, `HFX`, map factors) and the older
  `WRF` (the `AVGFLX_*` mass-flux variant). Use `WRFv2` unless you know you need
  the other.
- The run window must be covered by the files in `path`.
- Your wrfout must carry the fields listed in `WRFv2.vars`. If a field is
  missing, either add it to the WRF registry/output and rerun WRF, or comment
  the line out in a **copy** of the `.vars` file and point `variables` at that
  copy — EASLEE derives what it can (`ρ`, `θ`, `∂ρ/∂z`, and the surface
  turbulence chain) and warns about what it cannot.

What WRF gives you and what it does not (full matrix in
[`docs/MET_DATA.md`](docs/MET_DATA.md)):

| | WRF | ERA5 (+CLD) | GFS 0.25° |
|---|---|---|---|
| height of levels | `PH`+`PHB` | ak/bk integration | `gh − orog` |
| vertical motion | `W` | `etadot` | `ω` |
| turbulence | Hanna, from derived `u*`/`L`/PBLH | Hanna, archived `zust`/`blh` | Hanna, derived |
| convection | resolved (none parameterised) | **archived mass fluxes** | **CAPE-gated diagnostic** |
| wet deposition | — | ✓ 3-D precipitation flux | pseudo-3-D from `prate` |

That table is the reason to nest: a convection-permitting WRF domain resolves
what ERA5 has to parameterise, but it ends at its own boundary.

---

## 6. Nesting WRF inside ERA5 or GFS

List **one `[[meteorology]]` table per nest, finest first**:

```toml
[[meteorology]]          # nest 1 — WRF (finest)
model = "WRF"  ; path = "/data/WRF/my_run/"  ; variables = "WRFv2" ; available = "auto"

[[meteorology]]          # nest 2 — the global driver
model = "ERA5" ; path = "/data/ERA5/2021/"   ; variables = "ERA5_CONV" ; available = "auto"
```

How it behaves:

- **At release**, each parcel is assigned to the *finest* nest whose horizontal
  domain contains it.
- **On exit**, a parcel leaving its nest is handed to the next coarser nest that
  contains it — same lat/lon, vertical level index re-clamped to that nest's
  ladder — and keeps travelling. Only a parcel outside *every* nest is
  terminated.
- **Hand-off is one-way**: a parcel that later re-enters the WRF footprint stays
  on the coarse driver. So size the fine domain around the source, and don't
  expect re-capture.
- The nests need **not** share a grid, a projection or a level ladder — Lambert
  terrain-following WRF and lat/lon hybrid-level ERA5 coexist. They **do** need
  to cover the same run window.
- The **output grid is always geographic (lon/lat degrees)**, whatever
  projection each nest uses; parcels are mapped nest grid → lat/lon → output
  cell at deposit time, so both nests write into one field. Size the output box
  around the *far-field* plume, not around the WRF domain.
- Physics options apply run-wide but degrade gracefully per nest: in
  `wrf_era5_nest.toml`, `convection = "massflux"` acts on parcels sitting in the
  ERA5 nest (which has the archived mass fluxes) and is a no-op inside WRF,
  where convection is resolved.

`configs/wrf_era5_nest.toml` and `configs/wrf_gfs_nest.toml` are complete
worked examples — start from one of them rather than assembling this by hand.

---

## 7. The two features to know about

### 7.1 Inter-trajectory mixing (`[physics.mixing]`)

By default each parcel is a closed box: it keeps its tracer mass, and only its
*position* disperses. `scheme = "pairwise"` adds Fickian tracer exchange between
co-located parcels of the same species, plus coverage-based entrainment of
ambient air at the plume edge (an interior parcel, fully surrounded by
neighbours, sees none).

```toml
[physics.mixing]
scheme        = "pairwise"
K_h           = 1000.0     # horizontal diffusivity of the exchange [m²/s]
sigma0_h      = ""         # kernel half-width; empty = met grid spacing
sigma_v       = ""         # vertical pair half-width; empty = the turbulent yl
background    = true       # edge entrainment toward ambient air
nsig          = 2.0        # pair-search radius in σ₀ units
max_neighbors = 30         # cap on partners per parcel per pass
every         = 1          # apply every N sync steps
```

> **The one rule:** a turbulent K either **moves** parcels or **exchanges mass**
> between them — never both, or you double-count the same diffusion. So keep
> `[physics.turbulence].kh_ft = 0` whenever `[physics.mixing]` is on, and vice
> versa. Details and the scheme's derivation: the `[physics.mixing]` section of
> [`docs/CONFIG.md`](docs/CONFIG.md).

### 7.2 The kNN adaptive output kernel (`[output].kernel`)

How parcel mass becomes a gridded concentration:

| `kernel` | estimator |
|---|---|
| `"none"` | nearest cell |
| `"uniform"` | FLEXPART's 4-cell bilinear — the comparison reference |
| `"adaptive"` | **per-parcel Gaussian whose bandwidth is the local inter-trajectory spacing**, measured from the k nearest neighbours each output sample |
| `"covariance"` | per-parcel 2×2 covariance evolved by wind shear + turbulent growth, read at evaluation |

`"adaptive"` is the kNN estimator: effective resolution follows the local
trajectory density instead of a fixed cell size, and there is no FLEXPART-style
3-hour age rule. Add `resolution = true` to also write `R(x)`, the mass-weighted
kernel width per cell — read that field to see where the run is genuinely
resolved and where the grid is interpolating between a handful of parcels.

```toml
[output]
kernel          = "adaptive"
resolution      = true
bandwidth_floor = "dispersion"   # honest footprint for just-released parcels,
floor_sigma_v   = 0.5            # whose kNN bandwidth is still 0
floor_tl        = "5min"
```

More in [`docs/OUTPUT_GRID.md`](docs/OUTPUT_GRID.md).

---

## 8. Where to look next

| topic | file |
|---|---|
| every config field, block by block | [`docs/CONFIG.md`](docs/CONFIG.md) |
| what each met dataset provides / enables | [`docs/MET_DATA.md`](docs/MET_DATA.md) |
| the inter-trajectory mixing scheme | [`docs/CONFIG.md`](docs/CONFIG.md) (`[physics.mixing]`) |
| output grid, kernels, memory | [`docs/OUTPUT_GRID.md`](docs/OUTPUT_GRID.md) |
| code layout | [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) |
| vertical coordinates (kinematic / hybrid / isentropic) | [`docs/VERTICAL_COORDS.md`](docs/VERTICAL_COORDS.md), [`docs/ISENTROPIC.md`](docs/ISENTROPIC.md) |
| convection schemes | [`docs/CONVECTION.md`](docs/CONVECTION.md) |
| aerosols (size bins, settling, deposition, AOD) | [`docs/AEROSOLS.md`](docs/AEROSOLS.md) |
| species database | [`docs/SPECIES_DB.md`](docs/SPECIES_DB.md) |
| timestep hierarchy | [`docs/TIMESTEPS.md`](docs/TIMESTEPS.md) |
| ETEX / CAPTEX validation results | [`docs/ETEX_RESULTS.md`](docs/ETEX_RESULTS.md), [`docs/CAPTEX_RESULTS.md`](docs/CAPTEX_RESULTS.md) |
| Windows click-to-run | [`docs/WINDOWS.md`](docs/WINDOWS.md) |

---

## 9. Troubleshooting

**"config not found" / paths.** Run from the root of this directory; config
paths in the TOML (`[[meteorology]].path`, `[output].path`) are used as given,
so absolute paths are safest.

**`read_runspec` throws a list of errors.** That is the point — it validates
every field and reports all problems at once, before anything expensive starts.
Fix them top to bottom.

**"unknown configuration key … did you mean …".** Every TOML key is checked
against the schema; a misspelled or retired key is refused rather than
silently ignored (it used to be). Rename it as suggested, or put
`strict_keys = false` at the top of the file to load anyway with a warning.

**A config key that seems to do nothing.** Check the run log: the engine logs
each armed mechanism at startup (`plume rise enabled`, `mixing`, `prognostic air
mass armed`, …). If your setting produced no such line, it was not reached.

**The run dies during met reading.** The window in `[run]` must be fully covered
by files in `[[meteorology]].path`. With `available = "auto"` the startup log
prints how many valid times were indexed and their range — compare that to
`[run].start`/`stop`.

**Memory.** Grid memory is dominated by the output columns. Keep
`[output].group = "species"` (the forward default: one column per species,
independent of the number of release boxes). `group = "release"` gives one
column per box — needed for per-box backward sensitivities, but it scales with
box count and can be very large.

**Two big runs at once.** Two concurrent multi-GB NetCDF writers have been seen
to fail with `NetCDF error code -101` in the per-frame reopen. Chain long runs
rather than running them side by side.

**Checks that need real data** — reader tests, convection, PBL height, the
CAPTEX/ETEX tracer-experiment scoring — live in `validation/` rather than
`test/`, because each one loads a met archive from a path you have to set.
See [`validation/README.md`](validation/README.md).

**Slow GUI start.** Expected on the first launch (Julia compiles as it goes).
On Windows, `build_sysimage.bat` removes it. On Linux the equivalent is
`julia --project=webgui webgui/build_sysimage.jl`.
