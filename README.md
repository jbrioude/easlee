<p align="center">
  <img src="assets/easlee_logo.png" alt="EASLEE logo" width="180">
</p>

# EASLEE — transport distribution

An atmospheric **Lagrangian particle dispersion model** (LPDM) written in
**Julia**, in the lineage of FLEXPART. Everything is driven by a **single TOML
config** (or the web GUI) through one engine and one entry point — there are no
per-task scripts to edit.

```bash
julia --project=easlee scripts/easlee.jl configs/wrf_single.toml
```

**→ Start with [HOWTO.md](HOWTO.md)**: install (Linux and Windows), run from the
command line or the web interface, point EASLEE at WRF data, and nest WRF inside
ERA5 or GFS.

> **Scope of this distribution.** This is the **transport engine**. The
> reactive-chemistry module (gas-phase mechanisms, solvers, photolysis,
> heterogeneous chemistry, emissions inventories, adjoint) is not part of it.
> Inter-trajectory **mixing** and the **kNN adaptive output kernel** *are*
> included.

**→ What changed in v1.1** (transport, met reading, output — no chemistry):
[CHANGES.md](CHANGES.md).

## Features

- **Transport**: forward & backward (source–receptor sensitivity) on **WRF**,
  ERA5, NCEP GFS, MERRA-2, AROME(-HDF) and Meso-NH meteorology, with arbitrary
  **nesting** (finest domain first; parcels are handed to the coarser driver on
  exit); polar-stereographic handling at the poles on global grids.
- **Turbulence** (Hanna, Thomson, driver-native TKE where archived) and
  **convection** (archived mass-flux, diagnostic CAPE-gated), **dry/wet
  deposition**, radioactive **decay**, gravitational **settling**.
- **Inter-trajectory mixing** (`[physics.mixing]`): Fickian tracer exchange
  between co-located parcels plus coverage-based edge entrainment — the
  mass-exchange alternative to a horizontal random walk.
  See the `[physics.mixing]` section of [docs/CONFIG.md](docs/CONFIG.md).
- **Output estimators**: FLEXPART-uniform, **kNN adaptive kernel** (bandwidth =
  local inter-trajectory spacing), and a flow-evolved covariance kernel; plus a
  per-cell resolution field `R(x)`.
  See [docs/OUTPUT_GRID.md](docs/OUTPUT_GRID.md).
- **Aerosols**: lognormal / binned size distributions, κ-Köhler growth, Mie AOD;
  OSSA sea-spray source (Ovadnevaite wave-Reynolds flux, mean/gust split);
  in-flight species conversion.
- **Vertical coordinates**: kinematic, hybrid and isentropic.
- **Plume rise**: biomass-burning (Sofiev FRP default, Freitas 1-D option) and a
  volcanic module (1-D moist eruption column + Lagrangian plume phase with
  two-phase washout).
- **Species database**: FLEXPART species with molar mass, dry/wet deposition and
  decay ([docs/SPECIES_DB.md](docs/SPECIES_DB.md)).
- **Wind dimensionality** (3D / 2D / 1D) to isolate vertical mixing and convection.
- **Warm start**: full-state checkpoints (`[run] checkpoint`/`restart`) — extend
  a finished run without re-simulating.
- **Output**: time-framed gridded NetCDF (FLEXPART LOUTSTEP/LOUTAVER semantics),
  gridless receptor time series, per-trajectory CSV; flexible grids and vertical
  coordinates (height / pressure / model level / theta).
- **Web GUI** for building configs, running, and plotting.
- Optional **GPU** particle advance (`EASLEE_ADVANCE=ka_gpu`; CUDA is installed
  only if you use it). Every run prints one startup line saying whether the
  GPU is actually in use (`GPU backend enabled: <device>; CUDA initialized`)
  or not, and why — see `HOWTO.md` §1.3.

Validated against tracer experiments: **CAPTEX-83** (RANK ≈ 2.5) and **ETEX-1**
(competitive with FLEXPART 11) — see [docs/CAPTEX_RESULTS.md](docs/CAPTEX_RESULTS.md)
and [docs/ETEX_RESULTS.md](docs/ETEX_RESULTS.md).

## Layout

```
easlee/          the engine (Julia source) + its run environment
scripts/         easlee.jl (entry point) and run_config.jl (lower-level runner)
configs/         example TOML run configs, including the WRF ones
listvariables/   met variable namelists (<NAME>.vars), one per dataset
docs/            reference documentation
test/            self-contained tests (synthetic data, no met archive needed)
validation/      checks that need a real met archive + CAPTEX/ETEX DATEM scoring
webgui/          the web interface
tools/           plotting and profiling utilities
```

## Quick check

```bash
julia --project=easlee -e 'import Pkg; Pkg.instantiate()'   # once
julia --project=easlee test/runtests.jl                     # 30 tests, no data needed
```

## Licence and citation

GPL-3.0-or-later ([LICENSE](LICENSE)). If you use EASLEE in your research,
please cite it as described in [CITATION.cff](CITATION.cff).
