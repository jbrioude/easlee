# validation/ — checks that need real data

Everything here reads a **meteorological archive** (or observations) from disk,
so these do *not* run out of the box: edit the `path` in the config each one
loads, or point it at your own archive. The self-contained tests that need no
data at all live in [`../test/`](../test) and are the ones to run first
(`julia --project=easlee test/runtests.jl`).

Run any of them directly:

```bash
julia --project=easlee validation/test_era5_io.jl
```

| file | what it checks | data it needs |
|---|---|---|
| `test_era5_io.jl` | ERA5 GRIB reader: fields, levels, hydrostatic `zlev` | ERA5 model-level GRIB |
| `test_gfs_transport.jl` | GFS GRIB2 reader + a short transport run | NCEP GFS 0.25° pgrb2 |
| `test_isentropic.jl` | hybrid / isentropic vertical transport | ERA5 (+ PT winds) |
| `test_mesoscale.jl` | mesoscale sub-grid wind fluctuations | ERA5 |
| `test_polar_stereo.jl` | polar-stereographic handling on a global grid | global ERA5 |
| `test_aerosol.jl` | size-resolved settling / deposition / AOD in a run | ERA5 (`configs/example_dust.toml`) |
| `test_antimeridian.jl` | dateline-crossing output grid | ERA5 (`configs/era5_conv.toml`) |
| `test_backward_srr.jl` | backward source–receptor sensitivity | ERA5 (`configs/era5_backward.toml`) |
| `test_conv_entrain_detrain.jl` | mass-flux convection entrainment/detrainment | ERA5 + CLD companions |
| `test_diag_convection.jl` | CAPE-gated diagnostic convection | NCEP GFS (`configs/gfs.toml`) |
| `test_kh_ft.jl` | free-atmosphere horizontal random walk (`kh_ft`) | ERA5 over the ETEX domain (`configs/test_wind.toml`) |
| `test_pblh_diagnostic.jl` | bulk-Richardson PBL height | ERA5 (`configs/era5_pacific.toml`) |
| `test_plume_rise.jl` | fire / volcanic plume-rise injection heights | ERA5 (`configs/era5_conv.toml`) |
| `test_wind_dim.jl` | wind dimensionality (3D / 2D / 1D) | ERA5 over the ETEX domain (`configs/test_wind.toml`) |

## Tracer-experiment validation (DATEM)

`validate_captex.jl` + `check_captex.jl` / `check_etex.jl` score a run against
the CAPTEX-83 and ETEX-1 measurement archives (ATMES-II statistics via
`datem.jl`). These are driven from a config with `task = "validation"`:

```bash
julia --project=easlee scripts/easlee.jl configs/validate_etex.toml
```

Published scores for reference: [`../docs/CAPTEX_RESULTS.md`](../docs/CAPTEX_RESULTS.md),
[`../docs/ETEX_RESULTS.md`](../docs/ETEX_RESULTS.md).

## check_aod_optics_cams.jl — aerosol optics vs a CAMS analysis (2026-09-26)

Optics-only check: EASLEE's per-class AOD550 from the CAMS mixing ratios of one analysis
against the CAMS 532-nm extinction column of the same file, RH from the ERA5 file of the same
hour. No transport, no chemistry. PASS = ratio in [0.7, 1.4] and r > 0.7 (reference: 1.11 / 0.75).
