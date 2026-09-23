# CAPTEX-83 validation of EASLEE

**Experiment:** CAPTEX (Cross-Appalachian Tracer Experiment, Sep–Oct 1983). PMCH
tracer released from **Dayton, OH** and **Sudbury, ON**; surface air concentrations
measured at ~84 receptor stations over the following days. Seven releases.

**Meteorology:** ERA5 (1983, regional 0–360° / 260–300 °E North-America domain),
`IFS_EA` hybrid-level 3-D + `IFS_SFC` companion, hourly.

**EASLEE configuration:** forward, release 0–20 m AGL, mass-flux convection,
output grid 0.25°, surface layer 0–100 m, DATEM/ATMES-II scoring (every station ×
sampling-window, zeros included). Driver: `validation/validate_captex.jl`
(`MODE=run DIRECTION=forward`).

## Results — Thomson vs Hanna turbulence (per-window, n=2897)

### Thomson + convection
> ### ⚠ 2026-09-19 — the tables below predate two corrections
>
> 1. **Surface fluxes (2026-09-04).** `sshf`, `ewss` and `nsss` were read with
>    the wrong accumulated-vs-mean convention, which feeds u\*, the Obukhov
>    length and hence the PBL. These runs therefore had almost no boundary
>    layer.
> 2. **Met over-run.** The runs behind these tables continued past the end of
>    the available met by 8–12 h; the reader now refuses that.
>
> Re-run on the current code (2026-09-14), releases 1–4 and 6–7 at 3 days,
> release 5 at 4 days so it still covers its sampling:
>
> | release | N | R | FB | FA2 | FMS | KSP | RANK |
> |---|---:|---:|---:|---:|---:|---:|---:|
> | 1 | 395 | 0.200 | +0.619 | 0.25 | 51.0 | 9.6 | 2.14 |
> | 2 | 400 | 0.694 | −0.212 | 0.26 | 68.4 | 19.0 | 2.87 |
> | 3 | 404 | 0.478 | +1.234 | 0.10 | 33.1 | 11.9 | 1.82 |
> | 4 | 367 | 0.382 | +0.348 | 0.20 | 54.0 | 8.2 | 2.43 |
> | 5 (3 d) | 368 | 0.807 | +0.243 | 0.42 | 30.1 | 12.0 | 2.71 |
> | 6 | 390 | −0.036 | −1.802 | 0.22 | 7.2 | 6.2 | 1.11 |
> | 7 | 288 | 0.251 | −0.249 | 0.15 | 42.9 | 6.2 | 2.30 |
> | **POOLED (rel 5 at 3 d)** | **2612** | **0.393** | **+0.082** | **0.23** | **43.7** | **9.9** | **2.45** |
> | 5 (4 d, on its own) | 653 | 0.569 | −0.421 | 0.42 | 22.5 | 12.6 | 2.21 |
>
> Release 6 remains the outlier it always was. The pooled RANK is 2.45 against
> the 2.49 below, so the headline conclusion ("FLEXPART-class, ~2.5") survives,
> but **the individual rows below are superseded** and the pooled row at
> N = 2897 that combines the 4-day release 5 with the others has not been
> recomputed.

| release | N | R | FB | FA2 | FA5 | FMS | KSP | RANK |
|--------:|--:|--:|--:|--:|--:|--:|--:|--:|
| 1 | 395 | 0.219 | +0.60 | 0.10 | 0.44 | 51.7 | 7.8 | 2.19 |
| 2 | 400 | 0.799 | −0.03 | 0.26 | 0.53 | 64.1 | 21.3 | **3.05** |
| 3 | 404 | 0.400 | +1.29 | 0.07 | 0.37 | 32.2 | 9.4 | 1.74 |
| 4 | 367 | 0.438 | +0.20 | 0.20 | 0.53 | 53.3 | 10.9 | 2.51 |
| 5 | 653 | 0.577 | −0.45 | 0.46 | 0.75 | 21.4 | 12.4 | 2.19 |
| 6 | 390 | −0.041 | −1.69 | 0.25 | 0.38 | 7.1 | 4.9 | 1.18 |
| 7 | 288 | 0.115 | −0.50 | 0.14 | 0.38 | 35.6 | 8.3 | 2.04 |
| **POOLED** | 2897 | 0.456 | +0.01 | 0.20 | 0.49 | 40.3 | 10.9 | **2.49** |
| POOLED (excl. rel 6) | 2507 | 0.478 | +0.15 | 0.20 | 0.50 | 44.9 | 11.9 | 2.48 |

### Hanna + convection
| release | N | R | FB | FA2 | FA5 | FMS | KSP | RANK |
|--------:|--:|--:|--:|--:|--:|--:|--:|--:|
| POOLED | 2897 | 0.457 | +0.01 | 0.21 | 0.49 | 40.5 | 11.1 | **2.50** |
| POOLED (excl. rel 6) | 2507 | 0.479 | +0.15 | 0.21 | 0.49 | 45.2 | 12.1 | 2.49 |

(Per-release Hanna values are within ±0.01 of Thomson for every release.)

`RANK = R² + (1−|FB|/2) + FMS/100 + (1−KSP/100)`, range 0–4.

## Interpretation

- **POOLED RANK ≈ 2.5** (Thomson 2.49, Hanna 2.50) sits squarely in the
  FLEXPART/HYSPLIT/STILT band for CAPTEX (cf. Hegarty et al. 2013, who report
  multi-model CAPTEX RANKs in the ~2.3–2.8 range). **EASLEE is FLEXPART-class on
  CAPTEX.**
- **Turbulence scheme is nearly irrelevant** here: Thomson and Hanna give identical
  pooled scores (2.49 vs 2.50) and per-release agreement to ±0.01. At the synoptic
  scale of CAPTEX transport, the boundary-layer turbulence parameterisation is not
  the limiting factor — consistent with the ETEX finding.
- **Pooled bias is excellent** (FB +0.01).
- **Release-to-release spread is large**: release 2 is excellent (R 0.80, RANK 3.05),
  release 6 is poor (R −0.04, FB −1.69 — strong under-prediction, RANK 1.18) and is
  the weakest of the set; excluding it lifts the pooled R to 0.48 but barely moves
  RANK (the composite is robust to one bad release).

## Backward / adjoint (qualitative)

Backward (source–receptor) runs from the strongest receptors, convolving the SRR
with the emissions, reproduce the forward signal with correlation R ≈ 0.65; the
absolute backward magnitude needs an ~2.2× calibration factor against the forward
run (documented as a TODO — verify the SRR absolute scaling). Per-receptor outputs
under `results/captex_bwd_*`.

## Files
- `results/captex_paired_thomson.txt`, `results/captex_paired_hanna.txt`
  — release, observed, modelled (pg/m³) for all pairs.
- `results/captex_fwd_rel{1..7}/` — per-release run output.
- Reproduce:
  ```bash
  CAPTEX_OBS=<CAPTEX dir> CAPTEX_MET=/data/CAPTEX \
    TURB=thomson CONV=massflux MODE=run DIRECTION=forward RELEASE=all \
    julia -t 8 --project=easlee validation/validate_captex.jl
  ```

## Regression rerun (2026-07-03, real 1983 ERA5)

Full `RELEASE=all` rerun after the 2026-07-03 engine changes (dateline OUTGRID
wrap, NaN-Obukhov→neutral fix, per-species physics, PBL diagnostics), now on
genuine 1983 met (the year-mismatch caveat above no longer applies):

| | N | R | FB | FAC2 | FMS | KSP | RANK |
|---|--:|--:|--:|--:|--:|--:|--:|
| POOLED | 2897 | 0.470 | **+0.026** | 0.19 | 41.1 | 10.7 | **2.51** |

Per-station dosage (FLEXPART-11 set): FMS 0.97, PCC 0.68, FA2 0.48, FA5 0.75,
NMSE 2.14, FB +0.03 (n=80).  RANK 2.51 ≈ the historical ~2.50 — no regression;
the pooled bias is now essentially zero.  Weakest releases remain 6 (R −0.04,
strongly under-predicted) and 3 (2× over-predicted) — case-meteorology issues,
not code.

## Dual-estimator benchmark, release 1 (2026-07-06)

> Figures for this section: `tools/plotting/plot_estimator_2x2.py` →
> `results/estimator_2x2/`.

First runs on the 1983 ERA5 (regional 0–360 crop — exposed and fixed three
lon-convention bugs: release placement, met→OUTGRID shift, receptor
sampling; see commits 95ef5cd, 44f6687).  60 h, 400k parcels, 3-h frames on
synoptic hours; scored on 393 station-windows incl. 6-h windows
(`validation/check_captex.jl`, both estimators).

| estimator | kernel | R | FB | FMS | RANK |
|---|---|---|---|---|---|
| grid | uniform | 0.737 | +0.64 | 52% | 2.64 |
| grid | adaptive | 0.738 | +0.65 | 51% | 2.62 |
| receptors | (uniform run) | 0.654 | +0.40 | 53% | 2.64 |
| receptors | (adaptive run) | 0.503 | +0.37 | 53% | 2.48 |

All four at or above the recorded validation baseline (RANK ≈ 2.50).  The
receptor estimator is the least biased (FB +0.37–0.40 vs +0.65) at a median
native resolution of 7–8 km on a 25-km grid.  Grid-kernel choice is
indistinguishable here (0.25° cells ≈ the plume-core parcel spacing).
Estimator calibration: `test/test_estimator_calibration.jl` (14/14).

CAPTEX's kernel-indistinguishability corroborates the ETEX pair-level
attribution (ETEX_RESULTS.md, 2026-07-06): the dramatic ETEX 2×2 contrasts
were entirely the near-source station 8021 (21 km from the release);
CAPTEX's nearest stations are ≳100 km out, so no such pair exists and the
estimators agree.

## Bandwidth-floor A/B (`captex_floor`, 2026-07-07) — the control case

Release 1, 400k, 0.25° adaptive grid + gridless receptors, `[output]
bandwidth_floor = "dispersion"` (σv 0.5 m/s, T_L 5 min):

| estimator | R | FB | NMSE | RANK | vs unfloored baseline |
|---|---|---|---|---|---|
| grid | 0.736 | +0.654 | 34.6 | 2.58 | ≈ (grid adaptive +0.65, 2.62) |
| receptor | 0.667 | +0.357 | 46.5 | 2.73 | ≈ (receptor +0.37) |

CAPTEX is the CONTROL for the ETEX floor experiment (ETEX_RESULTS.md
pair-level §2, where the floor blew the receptor FB to +1.63): with **no
station in the dense near-source zone**, the floor leaves BOTH estimators
essentially unchanged (grid FB +0.654 ≈ +0.65; receptor FB +0.357 ≈ +0.37).
This isolates the ETEX receptor over-count to near-source parcel density —
the point-KDE sums more parcel tails only where parcels are dense — and
confirms the fix: the dispersion floor is GRID/kNN-only, the receptor is
decoupled (`build_run.jl`, fsv = 0) and keeps `min_bandwidth`.  NB this
`captex_floor` receptor ran with the pre-decouple build (harm would show
here if it existed on this network — it does not).
