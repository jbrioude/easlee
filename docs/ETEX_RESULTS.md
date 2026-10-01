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

## Reference run (uniform kernel, 72 h, n = 2629)

| statistic | EASLEE | ideal |
|---|---:|---:|
| model mean (obs. 109.2 pg m⁻³) | 151.1 pg m⁻³ | obs. |
| FB | +0.32 | 0 |
| FA2 | 0.35 | 1 |
| FA5 | 0.66 | 1 |
| R (per window) | 0.66 | 1 |
| FMS | 56.6 % | 100 % |
| KSP | 8.4 % | 0 % |
| **RANK** | **2.76** | 4 |

The composite rank is at the upper end of the range reported for FLEXPART,
HYSPLIT and similar models on ETEX-1. The plume footprint overlap (FMS) and the
concentration distribution (KSP) are both strong. Extending the run to the full
four sampling days with 10⁶ parcels (n = 3083) gives RANK 2.41, with a moderate
positive bias carried by the late, dilute windows.

**Turbulence scheme.** The Hanna scheme gives pooled scores within 0.01 of
composite rank of the Thomson scheme.

## Comparison with FLEXPART v11 (per-station dosage)

Bakels et al. (2024) score ETEX-1 on the per-station time-integrated dosage
(n = number of stations). Scored the same way, the four-day run gives

| metric | EASLEE | FLEXPART v11 |
|---|---:|---:|
| FMS | 0.80 | 0.68 |
| R | 0.74 | 0.84 |
| FA2 | 0.41 | 0.74 |
| FB | comparable | low |

EASLEE places the tracer in the right regions and matches the spatial pattern,
but the per-station amplitudes are noisier. The residual points to the
calibration of the near-source vertical dilution (boundary layer depth,
turbulence intensity, depth of the surface sampling layer) rather than to a
transport-position error.

## Output estimators: uniform kernel, adaptive kernel, gridless receptors

The same 72 h run scored with the adaptive kernel instead of the uniform
deposit is essentially unbiased (FB = −0.04) but its correlation drops to 0.33
(RANK 2.59). A pair-level attribution by distance from the source shows that the
entire contrast between the two estimators, in correlation and in bias, comes
from one station 20 km from the source (station 8021). The modelled filament
passes about 0.2° north of the observed track. The uniform deposit spreads the
filament over the neighbouring cells, fills the gap and counts as a detection,
whereas the adaptive kernel, whose local bandwidth there (0.1–1 km) resolves the
filament, reports the miss. Beyond 150 km from the source the two estimators
are statistically identical (band means within 2 %, equal per-band
correlations).

Two distinct points follow. The near-source azimuth error is a transport or
meteorology issue, not a kernel one. And the adaptive kernel separates model
error from representation error where the uniform deposit conflates them.

**Parcel-number sweep.** With the adaptive kernel the gridded scores are
independent of the parcel count from 5 × 10⁴ to 8 × 10⁵ parcels, while the
native resolution reported at the stations improves from a median of about
12 km to 2.5 km, following the N⁻¹ᐟ² law of the model description paper
(Eq. 15). Resolution is set by the number of trajectories, not by the output
grid. A 0.2° twin of the 400 k run reproduces the station-sampled fields beyond
150 km (correlation 0.95–0.98) and differs only at the near-source station,
where the 22 km cells bridge the gap to the misplaced filament.

**Gridless receptors.** The adaptive-kernel point evaluation at the station
coordinates reports a median native resolution of 4 km (p10 1 km, p90 21 km)
on the 400 k run. Its bias and NMSE improve monotonically with the parcel count
as the single-parcel spikes near the source shrink with the kNN bandwidth.
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
