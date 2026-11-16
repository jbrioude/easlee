# CAPTEX-83 validation of EASLEE

**Experiment.** Cross-Appalachian Tracer Experiment (Ferber et al., 1986). Seven
PMCH releases from Dayton, Ohio, and Sudbury, Ontario, in September and October
1983, sampled at about 84 surface stations.

**Meteorology.** ERA5 hourly, 1983, hybrid model levels (`IFS_EA`) with the
surface companion (`IFS_SFC`), regional North-America crop (260–300 °E).

**Configuration.** Forward runs, one per release, 0–20 m AGL release, 2 × 10⁵
parcels, Thomson turbulence, mass-flux convection, 0.25° output grid, 0–100 m
surface layer, three days per release (`RUNDAYS=3`). Scoring follows the
ATMES-II / DATEM statistics over every station × sampling window with zeros
included, summarised by the composite rank of Draxler et al. (2001),
`RANK = R² + (1 − |FB|/2) + FMS/100 + (1 − KSP/100)`, range 0–4.

## Results per release and pooled

Scores of v1.1.3. The model description preprint gave a pooled RANK of
2.43; the change comes from the surface-layer stability, which now follows the
ERA5 surface sensible heat flux.

| release | N | R | FB | FA2 | FMS % | KSP % | RANK |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 395 | 0.405 | +0.48 | 0.30 | 51.4 | 10.1 | 2.33 |
| 2 | 400 | 0.827 | −0.45 | 0.25 | 67.5 | 19.2 | 2.94 |
| 3 | 404 | 0.597 | +1.25 | 0.17 | 32.6 | 16.1 | 1.90 |
| 4 | 367 | 0.645 | +0.20 | 0.23 | 54.1 | 9.0 | 2.77 |
| 5 | 368 | 0.916 | +0.33 | 0.33 | 32.4 | 13.6 | 2.86 |
| 6 | 390 | −0.034 | −1.73 | 0.09 | 8.8 | 5.6 | 1.17 |
| 7 | 288 | 0.109 | −0.22 | 0.12 | 36.1 | 6.2 | 2.20 |
| **pooled** | **2612** | **0.462** | **−0.00** | **0.24** | **43.1** | **9.0** | **2.55** |

The pooled model mean is 210.3 pg m⁻³ against 211.3 observed. The pooled RANK
lies within the multi-model CAPTEX range of 2.3–2.8 reported by Hegarty et al.
(2013) for FLEXPART, HYSPLIT and STILT. On the per-station time-integrated dosage
(n = 80) the figure of merit in space is 0.95, the correlation 0.66, FA2 0.51,
FA5 0.79 and FB −0.00.

**Protocol of the preprint.** Table 3 of the preprint pooled the three-day runs
with a four-day run of release 5 (n = 2897, RANK 2.43). The same protocol now
gives RANK 2.48 (R 0.45, FB −0.09, FMS 41.6 %, KSP 9.4 %); the four-day release 5
alone scores 2.34 (n = 653, R 0.65, FB −0.34).

**Releases 5, 6 and 7 are not independent.** The seven releases use the same
tracer, and releases 5 (26 October, Sudbury), 6 (28 October, Dayton, 32 kg) and
7 (29 October, Sudbury, 183 kg) follow each other within three days. A run of
one release is scored against samples that contain the tracer of the others:
the fourth day of release 5 and most of the release 6 samples measure the plume
of release 7. The rows 5 to 7 of the table, and the apparent under-prediction of
release 6, should not be read as the skill of the model on each release.
Simulating the three releases together gives RANK 2.30 (n = 715, R 0.45,
FB +0.18, observed mean 180, model mean 215 pg m⁻³), and pooling this run with
releases 1–4 gives RANK 2.46 (n = 2281, R 0.49, FB +0.33, FMS 46.3 %, KSP 8.3 %).

The spread between releases 1 to 4 remains large, as expected for individual
synoptic cases (RANK 1.90 to 2.94).

**Turbulence scheme.** The preprint reported the Hanna scheme within 0.01 of the
Thomson scheme. That comparison has not been repeated since the corrections.

## Output estimators

CAPTEX is the control case for the ETEX-1 estimator comparison
(`ETEX_RESULTS.md`). No station lies in the dense near-source zone (the nearest
are beyond 100 km), and the uniform and adaptive output kernels give
indistinguishable scores. The gridded deposits were not affected by the position
correction listed in `CHANGES.md`, the output and meteorological grids having the same
spacing. In a release-1 experiment made before the corrections the gridless receptor
evaluation was the least biased estimator (FB +0.37, against +0.65 for the
gridded product of the same run) at a median native resolution of 7–8 km.

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
  TURB=thomson CONV=massflux MODE=run DIRECTION=forward RELEASE=all RUNDAYS=3 \
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
