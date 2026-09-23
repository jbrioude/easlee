# EASLEE — Architecture & File Reference

## Directory layout

```
easlee_v1.0/
  easlee/              ← model package (everything Julia needs to run)
    EASLEE.jl          ← module root
    runtime.jl         ← entry-point helpers (build, read met, step)
    Project.toml
    Manifest.toml
    src/
      core/            ← EASLEE module internals (orchestration)
      types/           ← data-type definitions
      physics/         ← advection & turbulence kernels
      output/          ← concentration accumulation
      met/             ← coordinate transforms & interpolation
      gpu/             ← GPU-safe utilities
      io/              ← file I/O (untyped layer)
  scripts/             ← entry point: easlee.jl (task dispatcher) + run_config.jl
  configs/             ← example TOML run configs
  listvariables/       ← met variables namelists (<NAME>.vars)
  docs/                ← documentation (this file and the other *.md)
  test/                ← self-contained test scripts
  validation/          ← met-archive validation scripts (CAPTEX/ETEX, DATEM, readers)
  tools/               ← generation + plotting utilities
  webgui/              ← web GUI (schema-driven config form + run/plot)
  assets/              ← static assets
  juliaout/            ← model output (netCDF, created at run time)
```

---

## Two-layer design

The code has two layers that talk to each other:

| Layer | Files | Role |
|---|---|---|
| **Untyped I/O** | `src/io/initialize.jl`, `src/io/tools.jl`, `src/io/output.jl` | Reads FLEXPART-style input files, allocates met arrays, writes netCDF. Uses plain `Dict` and `Any`-typed structs internally. |
| **Typed EASLEE** | everything in `src/core/`, `src/types/`, `src/physics/`, `src/output/`, `src/met/`, `src/gpu/` | Physics, advection, turbulence, concentration accumulation. Works only with concrete typed structs → zero dynamic dispatch in hot loops. |

`runtime.jl` is the bridge: it calls the untyped I/O layer to initialise and read met, then wires the results into the typed EASLEE context.

---

## One model step — what happens

```
read_met_step!(rt, currentdate)          # untyped layer: load met files into arrays
    └─ readmetfiles!()                   # src/io/initialize.jl

step_runtime!(rt, currentdate)           # typed: one full physics step
    └─ step_model!()                     # src/core/timeloop.jl
        ├─ release_parcels!()            # activate new parcels this step
        ├─ initialize_new_parcels!()     # set vertical level, turbulence state
        ├─ deposit_to_grid!()            # accumulate mass → concentration array
        ├─ advance_active_parcels!()     # move parcels (advection + turbulence)
        └─ reassign_escaped_parcels!()   # nest-switching for multi-nest runs
```

Output is written to netCDF by the untyped writer (`src/io/output.jl`) at the end of the run, or at each `outstep` if triggered from the script.

---

## File-by-file reference

### `easlee/EASLEE.jl` — module root

Defines the `EASLEE` module. Loads submodules in dependency order:
`Types → State → Config → MetAdapter → Operators → TimeLoop`

No logic here — just `include` and `using` statements.

---

### `easlee/runtime.jl` — entry-point helpers

Three functions used by every script:

| Function | What it does |
|---|---|
| `build_runtime(input)` | Reads grid/level/output config via the untyped layer, then builds the typed `ModelContext` and allocates `ParcelState`. Returns a named tuple `rt` with everything. |
| `read_met_step!(rt, date)` | Calls the untyped `readmetfiles!` for each nest. Met arrays are updated in-place; the typed structs see the new data automatically (shared memory). |
| `step_runtime!(rt, date)` | Thin wrapper around `step_model!`. |

---

### `src/core/state.jl` — parcel & release data structures

| Type | Description |
|---|---|
| `ParcelState{F,I}` | Struct-of-arrays for all parcels. One vector per field (position `indx/indy/indl`, turbulent velocity `uturb/vturb/wturb`, `mass`, `Dtime`, `status`, …). Default `F=Float32`. |
| `ReleaseSchedule{F,I}` | Release box geometry and timing (lon/lat bounds, z bounds, release times, number of parcels, mass). |
| `sync_dtime_ns!(state)` | Copies `Dtime` → `Dtime_ns` (Int64 ms) so GPU kernels can compare dates without `DateTime` objects. |

**Status codes for `state.status[i]`:** `0x00` = unused · `0x01` = active · `0x02` = deposited · `0x03` = outside domain.

---

### `src/core/config.jl` — context construction

| Type | Description |
|---|---|
| `NestContext{G,L,M3,M2}` | Everything about one nest: grid, level, met fields, run config. Fully parametric → compiler specialises on exact types. |
| `ModelContext` | Vector of `NestContext` + release schedule + output config + global run config. This is `rt.ctx`. |

| Function | What it does |
|---|---|
| `build_context(...)` | Converts untyped Dicts into typed structs and assembles `ModelContext`. Called once by `build_runtime`. |
| `build_release_schedule(...)` | Parses the `RELEASES` block from the input Dict into a typed `ReleaseSchedule`. |
| `to_typed_grid(g, globalx)` | Converts a untyped `Lambert` or `Rectilinear` object to the typed version. |
| `to_typed_level(l)` | Converts a untyped `sigmapressurelevel` object to `SigmaPressureLevel`. |

---

### `src/core/metadapter.jl` — time interpolation index

| Function | What it does |
|---|---|
| `met_time_index_array(date, nests)` | Returns a `Vector{Float32}` of length `nnest`. Each element is the fractional index between the two bracketing met files for the current date (used for 4D linear interpolation). |

---

### `src/core/operators.jl` — per-step physics operations

This module includes all the physics kernels (`src/physics/`, `src/met/`, `src/output/`, `src/gpu/`) and exposes five high-level functions:

| Function | CPU/GPU | What it does |
|---|---|---|
| `release_parcels!(state, ctx, date)` | **CPU** | Activates new parcels from the `ReleaseSchedule`. Distributes them randomly inside each release box. Assigns to finest nest that contains the point. |
| `initialize_new_parcels!(state, ctx, date, t_array)` | **CPU** | For each newly released parcel: finds the vertical level index, sets initial turbulence state (`wturb`, `dt`, `miss2`, `miss3`). |
| `advance_active_parcels!(state, ctx, date, t_array)` | **CPU multi-thread** | Calls `advance!` → `updatepos!` for each nest. The parcel loop inside `updatepos!` uses `Threads.@threads :static`. |
| `deposit_to_grid!(ctx, state)` | **CPU** | Calls `conccalc!` to accumulate parcel mass into the 5-D concentration array. Matches only parcels where `Dtime == currentdate`. |
| `reassign_escaped_parcels!(state, ctx, date)` | **CPU** | After advance, parcels that left their nest are converted back to lat/lon and assigned to the next coarser nest (multi-nest runs only). |

There are also GPU-path variants (used by the KA/GPU advance,
`EASLEE_ADVANCE=ka_gpu`):

| Function | CPU/GPU | What it does |
|---|---|---|
| `initialize_new_parcels_gpu!(backend, gs, gpu_met, ...)` | **GPU kernel** | KA version of parcel initialisation. |
| `deposit_conc_gpu!(backend, ...)` | **GPU kernel** | KA version of concentration accumulation. |
| `advance_parcel_kernel!(backend, blocksize)` | **GPU kernel** | KA version of advection + turbulence. |

---

### `src/core/timeloop.jl` — step orchestration

| Function | What it does |
|---|---|
| `step_model!(state, ctx, date; do_output, rng, …)` | Calls the operator functions in the correct order for one time step. |

Order matters: **deposit before advance** — `conccalc!` matches parcels where `Dtime == currentdate`; `advance!` then increments `Dtime` by `timestep`.

---

### `src/types/` — data type definitions

| File | Types defined |
|---|---|
| `grid.jl` | `Rectilinear`, `Lambert`, `AbstractGrid`; turbulence tags `Thomson`, `Hanna` |
| `level.jl` | `SigmaPressureLevel`, `AbstractLevel` |
| `met.jl` | `MetFields3D{A}`, `MetFields2D{B}`, `MetConfig`, `MetSnapshot` |
| `config.jl` | `RunConfig`, `NestRunConfig` |
| `output.jl` | `OutputConfig` |
| `parcel.jl` | Lightweight `ParcelState` used by the untyped layer (not the active one — `src/core/state.jl` is used) |

`MetFields3D` is the key type: it wraps 12 four-dimensional arrays `(Ni, Nj, Nk, 2)` for the two bracketing time slices. On CPU these are `Array{Float32,4}`; on GPU they become `CuArray{Float32,4}` — the struct is parametric so all physics code works unchanged on both.

`ConvFields` (also in `met.jl`) is the optional companion for convection/wet-deposition runs: 6 arrays (`umf`, `dmf`, `udr`, `ddr`, `tpf`, `cc`) with the same `(Ni, Nj, Nk, 2)` layout. Built from the untyped met Dict with `ConvFields(rt.metfields3D[n])` (shares array references; returns `nothing` when the variables file doesn't list the convective fields). The ERA5 GRIB reader fills them automatically from the `IFS_EA_CLD_*` companion of each `IFS_EA_MET_*` file, and `step_model!` applies convection + wet deposition when the `convfields` kwarg is passed (see CONVECTION.md).

---

### `src/physics/` — advection and turbulence

| File | CPU/GPU | Key functions |
|---|---|---|
| `turbulence.jl` | Both | `turbulence_v3` — main turbulence scheme (Thomson/BL89 or Hanna/no-turb). `turbulent_length_scale_v2`, `tke_partition` — building blocks. Uses pre-extracted 1D vertical profiles. |
| `turbulence_4d.jl` | **GPU** | `turbulence_v3_4d` — same physics as `turbulence_v3` but reads directly from 4D arrays without profile allocation (GPU-safe, no heap alloc in kernel). |
| `advance.jl` | **CPU** | `advance!(...)` → `updatepos!(...)` — advects all active parcels with `Threads.@threads :static`. Inner loop: sub-steps (`ctl` per timestep), mean wind + turbulence, boundary conditions, `Dtime += timestep`. |
| `advance_ka.jl` | **GPU** | `advance_parcel_kernel!` — KernelAbstractions kernel, one parcel per GPU thread. Uses `turbulence_v3_4d`. Also `_advance_one!` (inlined inner body). |
| `convection.jl` | **CPU multi-thread** | `convective_redistribution!` — Monte-Carlo parcel redistribution from archived IFS convective mass fluxes (ERA5 CLD companion files: `avg_umf`, `avg_dmf`, `avg_udr`, `avg_ddr`). Entrainment reconstructed from mass-flux continuity; captured parcels walk up/down the draught and detrain with probability D·Δz/M; uncaptured parcels subside with the compensating environment velocity −(M_u+M_d)/ρ, preserving the well-mixed criterion. |
| `wetdep.jl` | **CPU multi-thread** | `wet_deposition!` — in-cloud/below-cloud scavenging from the 3-D precipitation flux `avg_tpf` and cloud fraction `cc`. Λ = A·I^B (Hertel et al. 1995) with sub-grid precipitating-fraction correction fp^(1−B); mass decays as exp(−ΛΔt); removed mass accumulated on a per-nest 2-D grid (`WetDepState`). Parameters in `WetDepParams`. Skipped in backward runs. |
| `aerosol.jl` | **CPU multi-thread** | Size-resolved aerosol transport. `AerosolBins` (lognormal or explicit bins; per-bin Mie Q_ext from a Bohren-Huffman series validated against the Wiscombe benchmarks; refractive index by species — dust/sulfate/BC/OC/sea salt/ash). `apply_settling!` (Stokes-Cunningham + Clift-Gauvin transition regime), `dry_deposition_aerosol!` (resistance model, MO-similarity ra + Slinn rb), `wet_deposition_aerosol!` (size-dependent: Greenfield-gap collection efficiency below cloud, activation fraction in cloud), `compute_aod` (column optical depth on the output grid). Per-parcel size bin in `AerosolModel`; pass `aerosol=AerosolModel(...)` to `step_model!`. |

**Turbulence schemes** (set in `NestRunConfig.turbulence`):

| Tag | Scheme | Notes |
|---|---|---|
| `NoTurb()` | No turbulence | Pure advection — mean wind only. Reference/validation runs. |
| `Hanna(false)` | Hanna (option a) | Pre-computed TKE from Hanna; runs Thomson Langevin path. |
| `Hanna(true)` | Hanna (option b) | On-the-fly Hanna: σ_w and T_L analytic at every sub-step. |
| `Thomson(BL89, read_tls)` | Thomson/BL89 | Full TKE-based stochastic turbulence. AROME default. |

---

### `src/output/` — concentration accumulation

| File | CPU/GPU | Key functions |
|---|---|---|
| `conc.jl` | **CPU** | `conccalc!` — outer loop over nests, calls `_conc_nest!`. Deposits mass: maps parcel `(indx, indy, indl)` to output grid cell `(i, j, k)`, adds `mass / (rho * dz)` to `concarray`. |
| `conc_ka.jl` | **GPU** | `conc_parcel_kernel!` — KA kernel. Uses `Atomix.@atomic` for thread-safe writes. `deposit_conc_gpu!` is the launcher wrapper. |

The output grid is 5-D: `[Niout, Njout, Nkout, Nc, nage]` where `Nc` = number of release boxes and `nage` = number of age classes.

---

### `src/met/` — meteorological interpolation

| File | CPU/GPU | Key functions |
|---|---|---|
| `interp.jl` | Both | `interpolate4D` — 4D trilinear interpolation (i,j,k,t). `interpolation4Dsingle` — same for a single scalar. `get_vertical_profile!` / `get_vertical_profile` — extract a 1D vertical column. `interpolation1Dsingle` — 1D linear interpolation on a profile. `interpolateU/V/W` — wind-component interpolation with staggered-grid offsets. |
| `coords.jl` | Both | `latlontoindyindx` — geographic → grid index for `Rectilinear` and `Lambert`. `dxtodindx` / `dytodindy` — physical displacement (m) → grid-index increment. `indyindxtolatlon` — inverse. `ztoindl` / `indltolev` / `levtoindl` — vertical coordinate transforms. |

All functions are `@inline` and fully type-specialised — they are hot-path code called thousands of times per parcel per step.

---

### `src/gpu/rng.jl` — GPU-safe random numbers

Implements a per-parcel xorshift64 random number generator. Each parcel carries its own `rng_state::UInt64`. Functions:

| Function | What it does |
|---|---|
| `_xorshift64!(state)` | One xorshift64 step, returns a uniform `[0,1)` Float32. |
| `_randn_xorshift64!(state)` | Box-Muller transform on top of xorshift64 → standard normal Float32. |
| `_ensure_rng_states!(state, idx)` | Seeds `rng_state` from Julia's global RNG for newly released parcels. |
| `_sync_dtime_from_ns!(state, idx)` | Converts `Dtime_ns` (Int64) back to `DateTime` after a GPU kernel step. |

---

### Multi-source met input (namelists in `listvariables/`)

Each met dataset is described by a **variables namelist** selected by the
config's `meteorology.variables` field. The namelist maps input variable names to EASLEE
field names (`#3D` / `#2D` sections) and carries dataset properties in `#info`:

| `#info` key | Meaning |
|---|---|
| `model` | Selects the derived-field calculation (`ERA5`, `WRF`, `AROME`, `MERRA2`, `GFS`, `AROME_HDF`, `MESONH`) |
| `format` | `grib` or `netcdf` |
| `companions` | Filename-pattern substitutions for datasets that split fields across files at one valid time, e.g. `_MET_>_CLD_` (ERA5 conv) or `_asm_Nv.>_mst_Ne.` (MERRA2) |
| `leveltype` | Vertical coordinate of 3-D GRIB messages: `hybrid` (default), `isobaric` (GFS), `height` (AROME-HDF) |
| `multistep` | `true` when one GRIB file holds several valid times (AROME-HDF) — messages filtered by validity time |

Supported datasets and their derived-field paths (`src/io/metcalc_models.jl`):

| Dataset | Grid | Vertical | Notes |
|---|---|---|---|
| ERA5 (model levels) | regular lat-lon | 137 hybrid | ak/bk hydrostatic integration; optional CLD companion → convection/wetdep fields |
| MERRA2 | regular lat-lon 0.625°×0.5° | 72 eta | `PL`/`H` archived → no integration; `CMFMC`→umf, `PFL*/PFI*`→tpf (convection + wetdep work directly) |
| NCEP GFS | regular lat-lon 0.25° | 33 isobaric | zlev = gh − orog; ω→levdot on the fixed ladder; ustar/L/PBLH derived from archived stress + heat flux |
| AROME ATMO-HDF | regular lat-lon 0.025° | 22 height-AGL (20–3000 m) | BL-only dataset; no vertical velocity (levdot = 0, turbulence-driven); HP1 + SP1 companion |
| MesoNH | **Mercator** (RPK = 0) | Gal-Chen `ZHAT`/`ZS` | native TKE (`TKET`); θ (`THT`) direct; new `Mercator` grid type in `types/grid.jl` + `met/coords.jl` |
| WRF / AROME (NetCDF) | Lambert | eta / hybrid | reader unchanged in v1.1 |

### `src/io/initialize.jl` — untyped initialisation (Main scope)

Called once at startup. The FLEXPART-style inputs it reads (`pathnames`,
`COMMAND`, `RELEASES`, `OUTGRID`, `AVAILABLE`) are **temporary files emitted
internally by `build_run` from the TOML spec** — the user never writes them.
Key functions:

| Function | What it does |
|---|---|
| `readflexpartinput(pathnames)` | Reads the `pathnames` file, then reads `COMMAND`, `RELEASES`, `OUTGRID`. Returns the top-level `input` Dict. |
| `gridcheck_onefile!(path, input)` | Reads the first met file for a nest, determines grid type (AROME/ERA5), fills `metconfig[n]` Dict. |
| `initializerun(input, metconfig)` | Builds the global `runconfig` Dict from `COMMAND`. |
| `initoutgrid(runconfig, input, metconfig)` | Builds the `outputconfig` Dict from `OUTGRID`. Allocates the 5-D `concarray`. |
| `initgridconfig(metconfig, outputconfig, n)` | Builds the untyped `Lambert` or `Rectilinear` struct for nest `n`. |
| `initlevelconfig(metconfig, outputconfig, n)` | Builds the untyped `sigmapressurelevel` struct. |
| `initarrays(metconfig, gridconfig, levelconfig)` | Allocates the met arrays as `Dict{String, Array{Float32,4}}`. |
| `readmetfiles!(nestctx, date, dir, ...)` | Reads one or two met files bracketing `date` into the allocated arrays. Updates `t_index` in `NestRunConfig`. |

---

### `src/io/tools.jl` — utility functions (Main scope)

Miscellaneous helpers: `interpolate4D` (the untyped version used inside `initialize.jl`), coordinate helpers used by the reader, and `latlontoindyindx` (the untyped version — the typed version in `src/met/coords.jl` is used by the physics code).

---

### `src/io/output.jl` — netCDF output (Main scope)

| Function | What it does |
|---|---|
| `initnetcdf!(oc, filename)` | Creates the output netCDF file and writes coordinate variables. |
| `writenetcdf!(oc, filename, date)` | Writes the accumulated `concarray` to the file for one output time step. |
| `reset_accumulators!(oc)` | Zeros `concarray` and `massarray` between output intervals. |

---

## CPU vs GPU summary

| Component | CPU path | GPU path |
|---|---|---|
| **File I/O & initialisation** | ✅ always CPU | — |
| **Met reading** | ✅ always CPU | — |
| **Parcel release** | ✅ CPU | — (CPU even in GPU runs) |
| **Parcel initialisation** | ✅ CPU (`initialize_new_parcels!`) | ✅ GPU (`initialize_new_parcels_gpu!`) |
| **Advection** | ✅ CPU multi-thread (`updatepos!`, `Threads.@threads`) | ✅ GPU kernel (`advance_parcel_kernel!`) |
| **Turbulence** | ✅ CPU (`turbulence_v3`, 1D profiles) | ✅ GPU (`turbulence_v3_4d`, direct 4D access) |
| **Concentration** | ✅ CPU (`conccalc!`) | ✅ GPU kernel (`conc_parcel_kernel!`) |
| **Nest switching** | ✅ CPU | — (CPU even in GPU runs) |
| **netCDF output** | ✅ always CPU | — |

### Which path a run uses

Every run goes through `scripts/easlee.jl <config.toml>` →
`scripts/run_config.jl` (`read_runspec` → `build_run` / `run_spec!`). The CPU
multi-thread advance is the default; `EASLEE_ADVANCE=ka_gpu` selects the
KernelAbstractions GPU advance (falls back to CPU without a CUDA device).

### How the GPU path works

The GPU path is a manual orchestration loop — it does **not** call `step_runtime!`. Instead, each step it:
1. Calls `release_parcels!` on CPU
2. Uploads new parcel state to GPU
3. Launches `initialize_new_parcels_gpu!` (KA kernel)
4. Launches `deposit_conc_gpu!` (KA kernel)
5. Launches `advance_parcel_kernel!` (KA kernel)
6. Copies `Dtime_ns` back to CPU for the next step's parcel selection

The met arrays and output grids persist on the GPU across steps; only new-parcel fields travel CPU→GPU each step.

---

## Key data flow

```
TOML config  (read_runspec → RunSpec; build_run emits the
    │         internal FLEXPART-style temp inputs)
    ▼
readflexpartinput()          → input Dict
    │
    ▼
build_runtime(input)
    ├─ gridcheck_onefile!()  → metconfig[n] Dict (grid geometry)
    ├─ initializerun()       → runconfig Dict    (COMMAND)
    ├─ initoutgrid()         → outputconfig Dict (OUTGRID + concarray)
    ├─ initgridconfig()      → untyped Lambert/Rectilinear
    ├─ initlevelconfig()     → untyped sigmapressurelevel
    ├─ initarrays()          → metfields3D[n] Dict (empty Float32 arrays)
    └─ build_context()       → ModelContext (all typed structs)
           │
           └─ NestContext × nnest
                ├─ MetConfig      (grid metadata)
                ├─ gridconfig     (typed Lambert or Rectilinear)
                ├─ levelconfig    (typed SigmaPressureLevel)
                ├─ metfields3D    (typed MetFields3D — shares arrays with Dict)
                ├─ metfields2D    (typed MetFields2D)
                └─ nestrunconfig  (NestRunConfig: t_index, turbulence, ctl, direction)

─── time loop ────────────────────────────────────────────────────────────────

read_met_step!(rt, date)     fills metfields3D arrays in-place; updates t_index
step_runtime!(rt, date)  [or step_model! directly]
    ├─ release_parcels!        → state.indx/indy/nest/mass/Dtime set for new parcels
    ├─ initialize_new_parcels! → state.indl/wturb/dt/miss set
    ├─ deposit_to_grid!        → ctx.output.concarray += parcel contributions
    ├─ advance_active_parcels! → state.indx/indy/indl/uturb/wturb/Dtime updated
    └─ reassign_escaped_parcels! → state.nest/indx/indy updated for escapees
```
