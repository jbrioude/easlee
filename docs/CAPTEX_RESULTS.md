# CAPTEX-83 validation of EASLEE

**Experiment.** Cross-Appalachian Tracer Experiment (Ferber et al., 1986). Seven
PMCH releases from Dayton, Ohio, and Sudbury, Ontario, in September and October
1983, sampled at about 84 surface stations.

**Meteorology.** ERA5 hourly, 1983, hybrid model levels (`IFS_EA`) with the
surface companion (`IFS_SFC`), regional North-America crop (260–300 °E).

**Configuration.** Forward runs, one per release, 0–20 m AGL release, Thomson
turbulence, mass-flux convection, 0.25° output grid, 0–100 m surface layer. Each
release runs for the three days covered by its meteorological archive (four for
release 5). Scoring follows the ATMES-II / DATEM statistics over every station ×
sampling window with zeros included, summarised by the composite rank of Draxler
et al. (2001), `RANK = R² + (1 − |FB|/2) + FMS/100 + (1 − KSP/100)`, range 0–4.

## Results per release and pooled

| release | N | R | FB | FA2 | FMS % | KSP % | RANK |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 395 | 0.380 | +0.66 | 0.18 | 50.0 | 11.1 | 2.20 |
| 2 | 400 | 0.821 | −0.00 | 0.26 | 69.1 | 18.0 | 3.18 |
| 3 | 404 | 0.665 | +1.35 | 0.13 | 33.7 | 12.6 | 1.98 |
| 4 | 367 | 0.558 | +0.52 | 0.15 | 54.1 | 7.6 | 2.52 |
| 5 (4 days) | 653 | 0.653 | −0.36 | 0.33 | 20.7 | 13.5 | 2.32 |
| 6 | 390 | −0.038 | −1.71 | 0.18 | 8.6 | 6.4 | 1.17 |
| 7 | 288 | 0.147 | +0.13 | 0.12 | 34.2 | 6.6 | 2.23 |
| **pooled** | **2897** | **0.462** | **+0.19** | **0.21** | **41.5** | **10.2** | **2.43** |

The pooled RANK lies within the multi-model CAPTEX range of 2.3–2.8 reported by
Hegarty et al. (2013) for FLEXPART, HYSPLIT and STILT. On the per-station
time-integrated dosage (n = 80) the pooled figure of merit in space is 0.96, the
correlation 0.64 and FB +0.19.

The spread between releases is large, as expected for individual synoptic cases.
Release 2 reaches RANK 3.18 (R = 0.82) while release 6 is strongly
under-predicted (R = −0.04). Removing release 6 from the pool changes the RANK
by less than 0.02. Per-release scores are plotted in Fig. 6 of the model
description paper.

**Turbulence scheme.** Repeating the pooled evaluation with the Hanna scheme
changes the RANK by less than 0.01 and every per-release score by at most 0.01.
At the synoptic scale of CAPTEX the boundary-layer turbulence parameterisation
is not the limiting factor, as found on ETEX-1.

## Output estimators

CAPTEX is the control case for the ETEX-1 estimator attribution
(`ETEX_RESULTS.md`). No station lies in the dense near-source zone (the nearest
are beyond 100 km), and the uniform and adaptive output kernels give
indistinguishable scores. In a release-1 experiment the gridless receptor
evaluation is the least biased estimator (FB +0.37, against +0.65 for the
gridded product of the same run) at a median native resolution of 7–8 km on the
0.25° grid.

## Forward–backward consistency

Backward runs from the six sampling windows of the strongest release-2 receptor
(station 318) and the two strongest windows of station 316, with source and
receptor defined as the same output cell and the sensitivity convolved with the
release-2 emission over the release hours (Seibert and Frank, 2004), reproduce
the forward concentrations of the same runs with r = 0.998 and a pooled
backward-to-forward ratio of 0.88, without any scaling. See
`scripts/captex_fwdbwd_check.jl` and the per-receptor outputs under
`results/captex_bwd_*`.

## Reproduction

```bash
CAPTEX_OBS=<CAPTEX dir> CAPTEX_MET=<ERA5 1983 dir> \
  TURB=thomson CONV=massflux MODE=run DIRECTION=forward RELEASE=all \
  julia -t 8 --project=easlee validation/validate_captex.jl
```

Set `TURB=hanna` for the Hanna comparison. Backward runs use
`DIRECTION=backward` with the receptor list of `validate_captex_run.jl`.

## Files

- `results/captex_paired_thomson.txt`, `results/captex_paired_hanna.txt` —
  release, observed and modelled concentration (pg m⁻³) for every pair.
- `results/captex_fwd_rel{1..7}/` — per-release run output.
- `results/captex_bwd_*` — backward runs of the consistency check.
- `validation/validate_captex.jl`, `validation/validate_captex_run.jl` — drivers.
- `validation/check_captex.jl` — scoring of both estimators from an engine run.
- `validation/datem.jl` — shared DATEM / ATMES-II statistics.
