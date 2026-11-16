# ETEX-1 validation of EASLEE

**Experiment.** First European Tracer Experiment release (Nodop et al., 1998).
340 kg of PMCH released near Monterfil, Brittany (48.058 °N, 2.008 °W), over
11 h 50 min from 16:00 UTC on 23 October 1994, sampled at 168 surface stations
across Europe over the following days.

**Meteorology.** ERA5 hourly, hybrid model levels (`IFS_EA`) with the surface
companion (`IFS_SFC`), 23–27 October 1994. The domain crosses the prime meridian.
The GRIB header reports 350°→30° but ecCodes decodes the axis to −10°…30°,
which is the convention EASLEE uses internally, and the validation driver detects
it from the decoded `distinctLongitudes` axis.

**Configuration.** The reference run follows the FLEXPART ETEX-1 setup. Forward
transport, release at 8–20 m AGL, 4 × 10⁵ parcels, 10 min synchronisation step,
CTL = 10, IFINE = 10, Thomson turbulence, mass-flux convection, 0.1° output grid
(399 × 299 cells, −10…30 °E × 40…70 °N), 0–100 m surface layer, 3 h time frames
matching the sampling windows (`configs/etex_n800k.toml` with
`particles = 400000`). Model concentrations are paired with the observations at
every station and 3 h window over the first 72 h of sampling (n = 2629, zeros
included) and scored with the ATMES-II statistics, summarised by
`RANK = R² + (1 − |FB|/2) + FMS/100 + (1 − KSP/100)`, range 0–4.

## Reference run (uniform kernel, 72 h)

Scores of v1.1.2. They replace those of the model description preprint
(RANK 2.76), which were computed before two corrections listed in `CHANGES.md`:
the position of the gridded deposits on a regular meteorological grid, and the
surface-layer stability, which now follows the ERA5 surface sensible heat flux.

| statistic | all stations (n = 2629) | without station 8021 (n = 2607) | ideal |
|---|---:|---:|---:|
| observed mean | 109.2 pg m⁻³ | 99.2 pg m⁻³ | |
| model mean | 425.6 pg m⁻³ | 96.4 pg m⁻³ | obs. |
| FB | +1.18 | −0.03 | 0 |
| FA2 | 0.34 | 0.35 | 1 |
| FA5 | 0.69 | 0.69 | 1 |
| R (per window) | 0.56 | 0.49 | 1 |
| FMS | 56.2 % | 56.3 % | 100 % |
| KSP | 9.1 % | 9.1 % | 0 % |
| **RANK** | **2.19** | **2.70** | 4 |

The parcels are drawn at random, so two runs are not identical: the RANK of
repeated runs of this configuration varies by about 0.01.

**Station 8021** (20 km from the release, 22 of the 2629 pairs) carries
the whole difference between the two columns. Its observed mean is about
1300 pg m⁻³ and the model gives about 39 000 pg m⁻³ there, thirty times more.
With it the pooled bias, correlation and NMSE describe that one station; without
it the model mean is within 3 % of the observed mean. The plume footprint (FMS)
and the concentration distribution (KSP) do not depend on it. In the preprint the
deposits were shifted by 1.5 cells of the 0.1° grid, which moved the near-source
filament off this station and hid the over-prediction (RANK 2.76 with it,
2.64 without).

Extending the run to the four sampling days with 10⁶ parcels (n = 3083) gives
RANK 2.06 (R 0.43, FB +1.16, FMS 55.6 %, KSP 9.9 %), and 2.58 without station
8021. The preprint value of 2.41 is superseded.

**Turbulence scheme.** The preprint reported the Hanna scheme within 0.01 of the
Thomson scheme. That comparison has not been repeated since the corrections.

## Comparison with FLEXPART v11 (per-station dosage)

Bakels et al. (2024) score ETEX-1 on the per-station time-integrated dosage
(n = number of stations). Scored the same way, the four-day run gives

| metric | EASLEE (n = 159) | FLEXPART v11 (n = 168) |
|---|---:|---:|
| FMS | 0.80 | 0.68 |
| R | 0.70 | 0.84 |
| FA2 | 0.47 | 0.74 |
| FA5 | 0.80 | 0.87 |
| NMSE | 333 | 0.43 |
| FB | +1.16 | |

EASLEE places the tracer in the right regions, but the per-station amplitudes
are noisier. FB and NMSE are set by station 8021.

## Output estimators: uniform kernel, adaptive kernel, gridless receptors

The model offers two independent choices. The first is where the concentration
is evaluated: on the output grid (the station value is interpolated from the
four surrounding cells) or directly at the station coordinates, without a grid
(`[output.receptors]`). The second, for the grid only, is how a parcel is
deposited: over the four nearest cells (`kernel = "uniform"`, as in FLEXPART) or
with a Gaussian whose width follows the local parcel density
(`kernel = "adaptive"`). The gridless evaluation always uses the adaptive width.

RANK of the 72 h run, 4 × 10⁵ parcels:

| estimator | all stations | model mean | without 8021 | model mean |
|---|---:|---:|---:|---:|
| grid, uniform kernel | 2.19 | 426 | 2.70 | 96 |
| grid, adaptive kernel (1 run) | 2.23 | 436 | 2.70 | 93 |
| gridless receptors (5 runs) | 1.83–2.10 | 554–4066 | 2.64–2.70 | 87–103 |

Model means in pg m⁻³; observed 109.2 (all) and 99.2 (without 8021).

Without station 8021 the three estimators agree within the spread between
identical runs: on ETEX-1 the choice of estimator does not change the score. At
station 8021 all three over-predict. The gridless value there changes by a factor
of seven from one run to the next, because the point evaluation resolves
individual parcels 20 km from the source (bandwidth of the order of 1 km); the
all-station receptor score is therefore not reproducible at this parcel number.

The preprint attributed the contrast between the uniform and adaptive kernels
(RANK 2.76 against 2.59) to a filament passing 0.2° north of station 8021. That
displacement was the offset of the gridded deposits; it is corrected and the
attribution no longer holds.

**Parcel-number sweep.** The native resolution reported at the stations improves
from a median of about 12 km at 5 × 10⁴ parcels to 2.5 km at 8 × 10⁵, following
the N⁻¹ᐟ² law of the model description paper (Eq. 15). Resolution is set by the
number of trajectories, not by the output grid. The gridded scores of the sweep
and of the 0.2° twin (`configs/etex_n*.toml`, `configs/etex_dx02.toml`) were
computed before the correction of the deposit position and have not been
repeated.

**Gridless receptors.** The point evaluation at the station coordinates reports
a median native resolution of 4 km (p10 1 km, p90 21 km) on the 400 k run.
The receptor bandwidth is kept decoupled from the dispersion floor of the grid
deposit, since the point evaluation does not renormalise and a widened
bandwidth in a dense field over-counts. `min_bandwidth` in `[output.receptors]`
bounds the point-evaluation variance as an explicit experiment (default 0).

The estimator calibration tests (`test/test_estimator_calibration.jl`) cover
the slab normalisation, a prescribed-puff reproduction, the station-8021
filament case, deposit mass conservation on the real ETEX grid, the window
semantics and the longitude conventions. Run them before trusting any score.

## Reproduction

Engine run with time-framed output, then scoring of both estimators:

```bash
julia -t 8 --project=easlee run.jl configs/etex_n800k.toml   # particles = 400000
julia --project=easlee validation/check_etex.jl
```

Per-record pairs are written to `results/etex/etex_station_check.csv`.
The validation task of the driver gives the same statistics in one step:

```bash
CAPTEX_OBS=<ETEX dir> CAPTEX_MET=<ERA5 1994 dir> \
  TURB=thomson CONV=massflux MODE=run DIRECTION=forward RELEASE=all \
  NP=400000 RUNDAYS=3 julia -t 8 --project=easlee validation/validate_captex.jl
```

## Files

- `results/etex/etex_paired.txt` — observed and modelled concentration
  (pg m⁻³) for every pair.
- `results/etex/etex_scatter.png` — observed against modelled, log–log, with
  1:1 and factor-of-2 lines.
- `results/etex/etex_station_check.csv` — per-record pairs of `check_etex.jl`.
- `results/estimator_2x2/` — figures of the estimator comparison
  (`tools/plotting/plot_estimator_2x2.py`).
- `configs/etex_n*.toml`, `configs/etex_dx02.toml` — parcel-number sweep and
  0.2° twin.
- `validation/check_etex.jl`, `validation/validate_captex.jl`,
  `validation/datem.jl` — scoring and drivers.
