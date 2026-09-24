<p align="center">
  <img src="assets/easlee_logo.png" alt="EASLEE logo" width="180">
</p>

# EASLEE — the transport release

**Where does the air go?** Give EASLEE some weather data and a release — a
chimney, a wildfire, a volcano, an accidental leak — and it follows millions of
imaginary air parcels through the atmosphere to show where the material
travels, how it spreads out, and where it comes back down. Run it the other way
round and it answers the reverse question: for a measurement taken at one
place and time, which upwind regions did that air come from?

Technically, it is a **Lagrangian particle dispersion model** in the lineage of
FLEXPART, written in **Julia**. One run is described by a **single TOML file**
(or built in the web interface), and one command starts it — there are no
per-task scripts to edit.

```bash
julia --project=easlee scripts/easlee.jl configs/wrf_era5_nest.toml
```

**→ Start with [HOWTO.md](HOWTO.md)**: installing it on Linux and Windows,
running it from the command line or the web interface, pointing it at WRF data,
and nesting WRF inside ERA5 or GFS.

> **What is in this release.** This is the **transport engine** only. The
> reactive-chemistry module — gas-phase mechanisms, chemical solvers,
> photolysis, heterogeneous chemistry, emission inventories, the adjoint — is
> not included here. Mixing between trajectories and the adaptive output kernel
> *are* included.

**→ What changed in v1.1** (transport, weather-data reading and output — no
chemistry): [CHANGES.md](CHANGES.md).

## What it can do

- **Weather data.** Runs on **WRF**, ERA5, NCEP GFS, MERRA-2, AROME and
  Meso-NH. Grids can be **nested**: put a fine WRF domain inside a global ERA5
  or GFS field, list the finest first, and parcels leaving the fine domain are
  handed to the coarser one automatically. Regular latitude–longitude, Lambert
  conformal and polar-stereographic grids are all handled, including runs that
  cross the poles.

- **Both directions.** *Forward*: release material at a source and follow it
  downwind. *Backward*: start from an instrument and compute which upwind areas
  could have contributed to what it measured — the source–receptor sensitivity
  that inverse modelling needs.

- **Small-scale motion.** Random turbulent jostling (Hanna and Thomson
  schemes), and deep convection — the fast updraughts inside storm clouds —
  either read from the weather model or diagnosed from the instability it
  carries.

- **Loss along the way.** Material is removed by dry deposition at the surface,
  washout by rain, radioactive decay, and gravitational settling of heavy
  particles.

- **Turning parcels into a map.** The classic method averages parcels onto a
  fixed grid. EASLEE also offers an adaptive kernel that widens where the
  trajectories are sparse and tightens where they are dense (**kNN**), and
  reports a per-cell resolution field *R(x)* so you can see which parts of the
  map are well sampled and which are mostly noise. See
  [docs/OUTPUT_GRID.md](docs/OUTPUT_GRID.md).

- **Aerosols.** Particle size distributions (lognormal or binned), the swelling
  of particles as humidity rises (κ-Köhler), and optical depth from Mie theory.
  Sea spray is emitted from the wave state (Ovadnevaite), and species can
  convert from one form into another in flight.

- **The vertical.** EASLEE follows whatever vertical coordinate the weather
  model itself uses. Motion can be kinematic, hybrid, or isentropic (following
  surfaces of constant potential temperature).

- **Plume rise.** Smoke rises above a fire before it starts drifting (Sofiev's
  fire-radiative-power scheme by default, Freitas's 1-D column as an option).
  Volcanoes get their own treatment: a 1-D moist eruption column feeding a
  Lagrangian plume phase with two-phase washout.

- **Species database.** Each chemical species carries its molar mass, its dry
  and wet deposition behaviour, and its decay rate
  ([docs/SPECIES_DB.md](docs/SPECIES_DB.md)).

- **Wind dimensionality.** Run with full 3-D winds, or deliberately restrict
  the wind to 2-D or 1-D to isolate what vertical mixing and convection are
  contributing.

- **Warm start.** Checkpoints save the complete state of a run, so a finished
  run can be extended rather than simulated again from the beginning.

- **Output.** Gridded NetCDF in time frames (the same averaging semantics as
  FLEXPART's LOUTSTEP/LOUTAVER), time series at receptor points with no grid at
  all, and per-trajectory CSV. Output grids can be in height, pressure, model
  level, or potential temperature.

- **Web interface** for building a config, launching the run, and plotting the
  result.

- **GPU, optionally.** Set `EASLEE_ADVANCE=ka_gpu` to move the particle advance
  onto the GPU; CUDA is installed only if you ask for it. Every run prints one
  line at startup saying whether the GPU is genuinely being used
  (`GPU backend enabled: <device>; CUDA initialized`) or not, and why — see
  `HOWTO.md` §1.3.

EASLEE has been checked against two real tracer-release experiments, where a
measured gas was deliberately released and tracked by a ground network:
**CAPTEX-83** (RANK ≈ 2.5) and **ETEX-1**, where it is competitive with
FLEXPART 11. See [docs/CAPTEX_RESULTS.md](docs/CAPTEX_RESULTS.md) and
[docs/ETEX_RESULTS.md](docs/ETEX_RESULTS.md).

## What is in the box

```
easlee/          the engine (Julia source) and its run environment
scripts/         easlee.jl (the entry point) and run_config.jl (lower-level runner)
configs/         example TOML run configs, including the WRF ones
listvariables/   which variables to read from each weather dataset (<NAME>.vars)
docs/            reference documentation
test/            self-contained tests (synthetic data, no weather archive needed)
validation/      checks that need a real weather archive, plus CAPTEX/ETEX scoring
webgui/          the web interface
tools/           plotting and profiling utilities
```

## Check that it works

```bash
julia --project=easlee -e 'import Pkg; Pkg.instantiate()'   # once
julia --project=easlee test/runtests.jl                     # 30 tests, no data needed
```

## Contact

EASLEE is developed by **Jérôme Brioude** at the Laboratoire de l'Atmosphère et
des Cyclones (LACy), Université de La Réunion. Questions about using the model,
bug reports, and requests for the chemistry module are all welcome:

**jerome.brioude@univ-reunion.fr**

If you hit a problem that others might hit too, opening an issue on this GitLab
project is more useful than an email, since the answer stays visible.

## Licence and citation

GPL-3.0-or-later ([LICENSE](LICENSE)). If you use EASLEE in your research,
please cite it as described in [CITATION.cff](CITATION.cff).
