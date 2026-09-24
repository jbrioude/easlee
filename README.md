![](assets/easlee_logo.png "EASLEE logo")

# EASLEE : an Efficient AtmoSpheric Lagrangian Easy-to-usE particle dispersion model

Lagrangian particle dispersion models (LPDMs) are widely used tools in the atmospheric sciences, from the interpretation of in situ measurements to source–receptor analysis and the inverse modelling of surface fluxes. Their practical use, however, still requires a level of technical expertise (compiling Fortran code, per-dataset preprocessors, configuration fragmented over several files) that is unrelated to the physics of the problem being studied.

**EASLEE** (Efficient AtmoSpheric Lagrangian Easy-to-usE model) is a new LPDM in the lineage of FLEXPART, written entirely in the **Julia** language, and intended to remove that barrier. Because Julia is compiled just in time, no compilation step is exposed to the user: the same source code runs unmodified under **Linux, macOS and Windows**, and a run is defined by a single TOML configuration file or through a web-based graphical interface. Forward and backward (source–receptor) transport are both supported.

```bash
julia --project=easlee scripts/easlee.jl configs/wrf_era5_nest.toml
```

**→ Start with [HOWTO.md](HOWTO.md)**: installing it on Linux, macOS and Windows, running it from the command line or the web interface, pointing it at WRF data, and nesting WRF inside ERA5 or GFS.

> **What is in this release.** This is the **transport engine** only. The reactive-chemistry module (gas-phase mechanisms, chemical solvers, photolysis, heterogeneous chemistry, emission inventories, the adjoint) is not included here. Mixing between trajectories and the adaptive output kernel *are* included.

**→ What changed in v1.1** (transport, meteorological input and output, no chemistry): [CHANGES.md](CHANGES.md).

## What it can do

- **Meteorological input.** EASLEE ingests WRF, ERA5, NCEP GFS, MERRA-2, AROME and Meso-NH fields through a common level-index vertical transport formulation, so a new dataset is added by describing its variables in a small namelist rather than by writing new transport code. Domains can be **nested**: a fine WRF domain inside a global ERA5 or GFS field, finest listed first, with parcels handed to the coarser driver on exit. Regular latitude–longitude, Lambert conformal and polar-stereographic grids are handled, including transport across the poles of global grids.

- **Both directions.** *Forward*: release material at a source and follow it downwind. *Backward*: start from a receptor and compute which upwind areas could have contributed to what it measured. This is the source–receptor sensitivity that inverse modelling requires.

- **Turbulence and convection.** Boundary-layer turbulence through the Hanna and Thomson well-mixed Langevin schemes, and a mass-flux convection scheme that is well mixed by construction, with a CAPE-gated diagnostic scheme as a fallback for datasets that archive no fluxes.

- **Removal processes.** Dry and wet deposition, gravitational settling, and radioactive decay.

- **Gridded output.** The output uses an adaptive kernel in which each parcel deposits its mass over a support set by the local inter-parcel spacing, so that the spatial resolution of a simulation is determined by the number of trajectories rather than by the output grid. A companion resolution field *R(x)* reports this true local resolution at every grid cell. See [docs/OUTPUT_GRID.md](docs/OUTPUT_GRID.md).

- **Aerosols.** Size-resolved distributions (lognormal or binned) with κ-Köhler hygroscopic growth, size-dependent settling and scavenging, and an aerosol optical depth diagnostic computed from Mie extinction efficiencies. Sea spray is emitted from the wave state (Ovadnevaite), and species can convert from one form into another in flight.

- **Vertical transport.** EASLEE follows the native vertical coordinate of the driving model. Motion can be kinematic, hybrid, or isentropic–kinematic, the last of these intended for stratospheric transport.

- **Plume rise.** Biomass burning driven by fire radiative power (Sofiev by default, Freitas's 1-D column as an option), and a volcanic module: a 1-D moist eruption column feeding a Lagrangian plume phase with two-phase washout.

- **Species database.** Each species carries its molar mass, its dry and wet deposition behaviour and its decay rate ([docs/SPECIES_DB.md](docs/SPECIES_DB.md)).

- **Wind dimensionality.** Full 3-D winds, or a deliberate restriction to 2-D or 1-D to isolate the contribution of vertical mixing and convection.

- **Warm start.** Full-state checkpoints let a finished run be extended rather than simulated again from the beginning.

- **Output formats.** Gridded NetCDF in time frames (the averaging semantics of FLEXPART's LOUTSTEP/LOUTAVER), gridless receptor time series, and per-trajectory records. Output grids may be in height, pressure, model level or potential temperature.

- **Web interface** for building a configuration, launching the run and plotting the result.

- **GPU, optionally.** Set `EASLEE_ADVANCE=ka_gpu` to move the particle advance onto the GPU. CUDA is installed only if you ask for it. Every run prints one line at startup stating whether the GPU is genuinely in use (`GPU backend enabled: <device>; CUDA initialized`) or not, and why. See `HOWTO.md` §1.3.

The model is evaluated against the **ETEX-1** and **CAPTEX-83** tracer experiments, for which EASLEE obtains composite ranks of **2.76** and **2.51** respectively, comparable to other established LPDMs. See [docs/ETEX_RESULTS.md](docs/ETEX_RESULTS.md) and [docs/CAPTEX_RESULTS.md](docs/CAPTEX_RESULTS.md). It is also fast: the one-million-parcel, two-day ETEX-1 reference run completes in 12 min on 16 threads of a workstation, and in 5 min on its GPU ([docs/BENCHMARKS.md](docs/BENCHMARKS.md)).

## What is in the box

```text
easlee/          the engine (Julia source) and its run environment
scripts/         easlee.jl (the entry point) and run_config.jl (lower-level runner)
configs/         example TOML run configs, including the WRF ones
listvariables/   which variables to read from each meteorological dataset (<NAME>.vars)
download/        scripts that fetch public meteorology (GFS, MERRA-2, ERA5 archives)
data/            meteorology and validation observations (gitignored, see data/README.md)
docs/            reference documentation
test/            self-contained tests (synthetic data, no meteorological archive needed)
validation/      checks that need a real archive, plus CAPTEX/ETEX DATEM scoring
webgui/          the web interface
tools/           plotting and profiling utilities
```

## Check that it works

```bash
julia --project=easlee -e 'import Pkg; Pkg.instantiate()'   # once
julia --project=easlee test/runtests.jl                     # 42 tests, no data needed
```

## Authors and acknowledgements

EASLEE was written by **Jérôme Brioude**, with contributions from **Rémi Augerai**. The co-authors of the model description paper are warmly thanked for the discussions that shaped the model and for testing it on real cases.

The reader is referred to Stohl et al. (2005) for the physics and numerics inherited from FLEXPART.

## Contact

EASLEE is developed at the Laboratoire de l'Atmosphère et des Cyclones (LACy, UMR 8105 CNRS, Université de La Réunion, Météo-France), Saint-Denis, La Réunion, France. Questions about using the model, bug reports, and requests for the chemistry module are all welcome:

**jerome.brioude@univ-reunion.fr**

If you hit a problem that others might hit too, opening an issue on this GitLab project is more useful than an email, since the answer stays visible.

## Licence and citation

GPL-3.0-or-later ([LICENSE](LICENSE)). If you use EASLEE in your research, please cite it as described in [CITATION.cff](CITATION.cff).
