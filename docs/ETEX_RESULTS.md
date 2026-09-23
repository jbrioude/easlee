# ETEX-1 validation of EASLEE

**Experiment:** ETEX-1 (European Tracer Experiment, release 1). Continuous PMCH
release from **Monterfil, Brittany** (48.058 °N, 2.0083 °W), starting **16:00 UTC
23 Oct 1994**, **340 kg over 11 h 50 m**. Surface air concentrations measured at
**168 stations** across Europe over the following ~90 h.

**Meteorology:** ERA5 (`IFS_EA` hybrid-level 3-D + `IFS_SFC` companion), hourly,
`/data/ETEX/19941023`–`19941027`. The ETEX domain crosses the prime meridian:
the GRIB *header* reports 350°→30° (0–360) but ecCodes *decodes* to **−10°…30°**
(−180…180), which is the convention EASLEE uses internally. The validation driver
now auto-detects this from the decoded `distinctLongitudes` axis (commit `ab09618`).

**EASLEE configuration:**
- Turbulence: **Thomson** (well-mixed Langevin).
- Convection: **mass-flux** (well-mixed-by-construction scheme).
- 100 000 particles; 4-day forward transport; output grid −10…30 °E × 40…70 °N
  at 0.25°; surface layer 0–100 m; 3-h average sampling at every station/window.

## Two configurations

| run | particles | release | grid | step / CTL / IFINE |
|-----|----------:|---------|------|--------------------|
| **A** (first) | 100 000 | 0–20 m | 0.25° | 15 min / 5 / 4 |
| **B** (FLEXPART-match) | 1 000 000 | 8–20 m | 0.1° | 15 min / **10 / 10** |

Config **B** mirrors the FLEXPART ETEX-1 setup (1 M trajectories, 8–20 m AGL,
CTL=10, IFINE=10, 0.1° output grid).

> ### ⚠ 2026-09-19 — the tables below predate the surface-flux fix
>
> Every number in this file was measured before the 2026-09-04 correction to
> the archived surface fluxes (`sshf`, `ewss`, `nsss`: accumulated vs mean
> convention), which feeds u\*, the Obukhov length and therefore the PBL.
> Before that fix the boundary layer was effectively absent, so these runs are
> not the physics the current code produces, and **the RANK values here should
> not be quoted for the present version**.
>
> The most recent measurement on the current code, ETEX-1 at 0.1°, 400 k
> parcels, 72 h from 1994-10-23T12:00, Thomson turbulence with mass-flux
> convection, CTL 10, IFINE 10, adaptive output kernel, 3-h frames on a
> 399 × 299 grid at 0.1° (`configs/etex_n800k.toml` with `particles = 400000`),
> scored with `validation/check_etex.jl` over 2629 station-windows:
> **R 0.660, FB +0.327, NMSE 74.4, FA2 0.35,
> FA5 0.66, FMS 56.2 %, KSP 8.6, RANK 2.75.** That is a different
> configuration from config B below (400 k rather than 1 M parcels, a 72-h
> window rather than the full run), so it replaces none of the rows here; it
> is recorded so the gap is visible. A like-for-like re-run of config B is
> still to be done.

## Results (ATMES-II / FLEXPART-11 suite)

Paired over every station × sampling-window (zeros included): **n = 3083**
(obs > 0 in 938).

| statistic | A: 0.25°/100k | **B: 0.1°/1M (FLEXPART-match)** | target |
|-----------|------:|------:|--------|
| model mean (obs 96.2) | 188.8 | **83.7** pg/m³ | obs |
| **FB**   | +0.650 | **−0.139** | 0 |
| **FOEX** | −37.3 % | −38.2 % | 0 |
| **FA2**  | 0.349 | 0.341 | ↑ |
| **FA5**  | 0.669 | **0.702** | ↑ |
| **NMSE** | 468 | **42.1** | ↓ |
| **R**    | 0.425 | 0.335 | ↑ |
| **FMS**  | 54.6 % | 54.8 % | ↑ |
| **KSP**  | 10.3 % | 11.5 % | ↓ |
| **RANK** | 2.299 | **2.476** | ↑ (0–4) |

The FLEXPART-grade settings **removed the ~2× over-prediction** (model mean
189→84, FB +0.65→−0.14) and cut **NMSE 468→42** — the higher 8–20 m release plus
the finer grid and 1 M particles spread the near-source peak realistically (max
modelled conc 102 000→12 800 pg/m³). R and FA2 dip slightly: the classic
**finer-grid "double penalty"** — small position errors now fall in the wrong
0.1° cell. Net **RANK 2.48 ≈ EASLEE's CAPTEX (~2.50)** and FLEXPART-class.

### (config A) detail

| statistic | value | meaning (target) |
|-----------|------:|------------------|
| obs mean  | 96.2 pg/m³ | — |
| model mean| 188.8 pg/m³ | — |
| **FB**    | **+0.650** | factor-of-bias (0 = unbiased) |
| **FOEX**  | **−37.3 %** | factor of exceedance (0 best) |
| **FA2**   | **0.349** | fraction within factor 2 (↑) |
| **FA5**   | **0.669** | fraction within factor 5 (↑) |
| **NMSE**  | **468** | normalised mean square error (↓) |
| **PCC (R)** | **0.425** | Pearson correlation (↑) |
| **FMS**   | **54.6 %** | figure of merit in space (↑) |
| **KSP**   | **10.3 %** | Kolmogorov–Smirnov parameter (↓) |
| **RANK**  | **2.299** | composite 0–4 (↑) |

`RANK = R² + (1−|FB|/2) + FMS/100 + (1−KSP/100)`.

## Comparison to FLEXPART 11 (per-station dosage)

FLEXPART v11 (Bakels et al. 2024, GMD, Table 1) reports ETEX-1 statistics on the
**per-station time-integrated dosage** (n = #stations), *not* per-3h-window. Scored
the same way (config B, threads, n=159 stations):

| metric | EASLEE | FLEXPART 11 | note |
|--------|------:|------:|------|
| **FMS** (plume footprint) | **0.81** | 0.68 | EASLEE better |
| **rank correlation** | **0.79** | 0.84 (PCC) | comparable |
| FB (bias) | −0.14 | low | comparable |
| FA5 | 0.77 | 0.87 | close |
| FA2 | 0.41 | 0.74 | EASLEE trails |
| linear Pearson | 0.37 | — | outlier-dominated; use rank |

EASLEE is **competitive with FLEXPART 11 on ETEX-1**: it *beats* FLEXPART on the
plume footprint (which stations are hit, FMS 0.81 vs 0.68), matches it on the
spatial rank-pattern of dosage (0.79 vs 0.84) and on bias, and trails only on the
strict within-factor-2 amplitude (FA2 0.41 vs 0.74) — i.e. EASLEE puts the right
total tracer in the right places but the per-station *amounts* are noisier. That
residual is a dilution / vertical-mixing calibration target (PBL depth, Thomson
intensity), not a transport-position error. (Raw linear Pearson 0.37 is misleading
— a few high-dosage stations dominate it; the rank/log correlations, 0.79 / 0.53,
are the meaningful ones and what FLEXPART-class metrics use.)

## Interpretation

- **RANK 2.30** sits in the "good model" band of the ATMES-II ETEX-1
  intercomparison and is comparable to FLEXPART/HYSPLIT-class results — and to
  EASLEE's own CAPTEX RANK (≈2.50).
- **Spatial overlap (FMS 55 %) and distribution match (KSP 10 %) are strong** —
  the plume reaches the right stations and the concentration *distribution*
  matches the observations well.
- **The model over-predicts by ~2×** (FB +0.65; model mean 189 vs obs 96), while
  FOEX is −37 % — i.e. the modelled plume exceeds the obs in only ~13 % of pairs
  but is much higher where it does hit. This is the classic LPDM signature of a
  plume that is slightly **too peaky / under-dispersed in the core** while missing
  some low-concentration edges (also why FA2 is 0.35 despite FA5 0.67). The over-
  prediction is a calibration target: surface-sampling-layer depth (0–100 m vs the
  true ~ground-level intake), near-source vertical mixing, and release mass/timing.

## Self-check via the engine's time-framed grid output (2026-07-03)

The plain-transport engine (`run_spec!`) now writes TIME-FRAMED grids (FLEXPART
LOUTSTEP/LOUTAVER semantics — one frame per `output.interval`, mean over the
trailing `output.average` window), so ETEX can be checked WITHOUT the validation
task: run `configs/output_etex.toml`, then `validation/check_etex.jl` samples the
24 3-h frames bilinearly at every station/window (0–100 m layer) and scores with
the same `datem_stats` (shared `validation/datem.jl`).  Per-record pairs land in
`results/etex/etex_station_check.csv` for independent inspection/plotting.

Setup: 400 k parcels, 8–20 m release, 0.1° grid, Δt 10 min, CTL/IFINE 10/10,
Thomson + mass-flux, 72 h (→ 2629 station-windows; the 4-day runs above have
~3083).  **Result (pooled, zeros included):**

| statistic | framed engine (this check) | B: validation task (1M, 96 h) |
|-----------|------:|------:|
| model mean (obs 109.2 / 96.2) | 149.1 | 83.7 pg/m³ |
| **R**    | **0.651** | 0.335 |
| **FB**   | +0.309 | −0.139 |
| **FA2**  | 0.32  | 0.341 |
| **FA5**  | 0.65  | 0.702 |
| **FMS**  | **57.2 %** | 54.8 % |
| **KSP**  | **9.3 %**  | 11.5 % |
| **RANK** | **2.75** | 2.476 |

The per-window correlation more than doubles (0.34 → 0.65) and RANK rises to
2.75 — the best EASLEE ETEX score so far — confirming the trajectory →
concentration pipeline is sound once the output is properly time-averaged
per sampling window.  (The runs are not a strict A/B: 400 k vs 1 M parcels,
72 vs 96 h, Δt 10 vs 15 min; FB +0.31 partly reflects the shorter window
dropping late low-concentration pairs.)  The earlier "noisy/unusable" gridded
output was the missing LOUTSTEP/LOUTAVER semantics — the engine used to write a
single whole-run average, structurally incomparable with 3-h observations.

## Files
- `results/etex/etex_paired.txt` — release, observed, modelled (pg/m³) for all pairs.
- `results/etex/etex_scatter.png` — obs-vs-model log-log scatter with 1:1 and ±factor-2.
- Driver: `validation/validate_captex.jl` (`MODE=run`) + `validation/validate_captex_run.jl`.
  Reproduce:
  ```bash
  CAPTEX_OBS=<ETEX dir> CAPTEX_MET=/data/ETEX \
    TURB=thomson CONV=massflux MODE=run DIRECTION=forward RELEASE=all \
    NP=100000 RUNDAYS=4 julia -t 4 --project=easlee validation/validate_captex.jl
  ```

## Regression rerun of config B (2026-07-03, after the engine changes)

Same setup as B (1 M parcels, 8–20 m release, 0.1°, Thomson + mass-flux, 4 days)
rerun after the 2026-07-03 changes (dateline OUTGRID wrap, NaN-Obukhov→neutral,
per-species physics, uniform kernel active):

| statistic | rerun | B (baseline) |
|-----------|------:|------:|
| **R** (per-window, n=3083) | **0.491** | 0.335 |
| FB | +0.581 | −0.139 |
| FA2 | 0.33 | 0.341 |
| **FMS** | **56.7 %** | 54.8 % |
| KSP | 11.1 % | 11.5 % |
| RANK | 2.41 | 2.476 |
| dosage FMS / PCC / FA2 / FA5 (n=159) | 0.80 / **0.74** / 0.41 / 0.75 | 0.81 / 0.37 (linear) / 0.41 / 0.77 |

Correlations improved markedly (per-window R 0.34→0.49; per-station linear PCC
0.37→0.74, now approaching FLEXPART-11's 0.84) with FMS/KSP slightly better —
but the amplitude bias flipped to +0.58 (model ≈1.8× obs), pulling RANK to
2.41 ≈ baseline.  Candidate causes for the bias shift: the NaN-Obukhov fix
(NaN-L columns now mix as NEUTRAL instead of stable — more surface
concentration retained in the sampling layer) and the uniform output kernel;
follow-up: re-check the 0–100 m sampling-layer sensitivity.  Overall: no
regression in skill; the bias calibration item stays open.

## Dual-estimator benchmark: gridded vs gridless receptors (2026-07-06)

> Figures for this section: `tools/plotting/plot_estimator_2x2.py` →
> `results/estimator_2x2/`.

Same engine run (`configs/output_etex_uni.toml` / `output_etex_ab.toml`,
400k parcels, 72 h), two output estimators scored against the SAME 2629
station-windows (`validation/check_etex.jl`, both sections):

- **gridded** — 0.1° time-framed grid, bilinear at the stations;
- **receptors** — `[output.receptors]`: the adaptive-KDE point evaluation AT
  the station coordinates (no output grid; bandwidth = kNN spacing of the
  surrounding trajectories, validated against brute force in
  `test/test_estimator_calibration.jl`, 14/14).

| estimator | kernel | R | FB | FMS | RANK | notes |
|---|---|---|---|---|---|---|
| grid | uniform | 0.651 | +0.30 | 56% | 2.74 | = recorded baseline, reproduced twice |
| receptors | (gridless) | 0.642 | +1.07 | 59% | 2.36 | native resolution median 4 km (p10 1, p90 21) |
| grid | adaptive | 0.322 | −0.07 | 58% | 2.55 | mass now conserved (task #40); unbiased but smoother |
| receptors | (gridless, adaptive run) | 0.177 | +0.29 | 59% | 2.37 | |

Reading the 2×2 (SUPERSEDED — see the pair-level attribution below): the
apparent bias–variance trade-off (adaptive unbiased/smoother, uniform
sharp/biased) is NOT a field-level property.  A pair-level diff of the two
runs shows every contrast in this table — R AND FB — comes from the
near-source station 8021 alone.

### Pair-level attribution (2026-07-06): the 2×2 contrasts are ONE station

Joining the 2629 station-window pairs of the two 400k runs and splitting by
distance from the source:

| band | n | obs mean | uniform (R) | adaptive (R) |
|---|---|---|---|---|
| 0–150 km | 74 | 387 | 1626 (0.83) | **0.6** (0.80) |
| 150–400 km | 249 | 143 | 475 (0.73) | 467 (0.73) |
| 400–800 km | 549 | 87 | 87 (0.41) | 88 (0.42) |
| 800–2000 km | 1702 | 103 | 60 (0.51) | 61 (0.52) |

Beyond 150 km the two estimators are statistically IDENTICAL (means within
2 %, per-band R equal, Spearman on obs>0 pairs 0.584 vs 0.570; grid totals
conserved to 0.995–1.018 per frame).  The entire pooled contrast — R
0.651→0.322 AND FB +0.30→−0.07 (Δmean 46.2 of 46.2 pg/m³) — is the 74
near-source pairs, dominated by the four giant-obs windows of station 8021.

Mechanism (2-D field dump around 8021, frames 2–4): in BOTH runs the
modelled filament runs along ~48.25 °N; the station is at 48.07 °N, ~0.18°
south.  The uniform path smears the filament twice — bilinear deposit
(±1 cell) then bilinear station sampling (±1 cell) — bridging the ~20 km
gap and scoring 27k–38k pg/m³ against obs 0.9k–12.6k (hit, with the right
timing).  The adaptive kernel (h = kNN spacing ≈ 0.1–1 km ≪ 11 km cells)
deposits the filament into its own cells and reports the miss honestly.
Pearson-with-zeros is leveraged by exactly those pairs; per-band R is
unchanged.  The N-sweep confirms: adaptive grid R is 0.320 / 0.315 / 0.322
at 50k / 100k / 400k — N-invariant, i.e. a deterministic near-source
transport mismatch, not estimator noise.

Two real issues, correctly separated by the estimator pair:
1. **Model near-source azimuth error** (transport/met, NOT kernel): the
   filament passes ~0.2° north of the observed track at 21 km range.
   Candidates: night-stable-BL directional shear vs the 8–20 m release
   sampling, ERA5 low-level veer.  Open.
2. **kNN bandwidth ≠ physical plume width** (kernel — a `dispersion`
   bandwidth floor was TESTED and PARTLY REJECTED, 2026-07-07): the
   receptor's single-parcel spikes (one 0.85 g parcel at h ≈ 60 m ⇒
   ~750 000 pg/m³ point value; receptor R swings 0.64↔0.18 between
   realisations on the luck of one spike) come from h measuring SAMPLING
   spacing.  `[output] bandwidth_floor = "dispersion"` floors per-parcel h
   at the Taylor growth σ(age) = σv·age → σv·√(2·T_L·age) (σv = 0.5 m/s,
   T_L = 5 min → 0.7 km @ 1 h, 3.6 km @ 24 h).  The `etex_floor` A/B
   settled it:
   - on the **grid** the floor is NEUTRAL (RANK 2.57 vs 2.55 unfloored,
     R/FB unchanged) — the renormalised deposit conserves each parcel's
     mass, so a wider h just spreads that fixed mass over more of the 0.1°
     cells, which already exceed most floored h;
   - on the **receptor** it is HARMFUL (FB +0.29 → +1.63, NMSE 605 →
     8625; 8021 conc ×350).  The point-KDE does NOT renormalise
     (C = Σ mᵢ·G₂ integrates to mᵢ over ALL space), so widening h near a
     source — where parcel density is high — makes each station sum over
     many more parcels' tails and over-count.  Flooring bounds an ISOLATED
     spike (test_adaptive_kernel.jl [7], 28/28, single-parcel ×138) but
     inflates a DENSE field.
   Resolution: the floor is now GRID/kNN-only; the receptor is DECOUPLED
   (`build_run.jl`, fsv = 0) and keeps its honest native-resolution report
   plus the fixed `min_bandwidth` knob.  The real receptor fix — the
   age-growing PER-PARCEL footprint σ(age)/√N_local (absolute dispersion
   ÷ the local parcel count, which the point-KDE needs and the plain σ(age)
   over-estimates by exactly that √N) — stays open.

### The near-source station 8021 (the receptor's story in one station)

Station 8021 sits ~20 km from the Rennes source.  Splitting the receptor
statistics:

| set | n | obs mean | model mean | R | FB |
|---|---|---|---|---|---|
| all stations | 2629 | 109 | 358 | 0.642 | +1.07 |
| excl. 8021 | 2607 | 99 | 82 | 0.454 | −0.19 |

One station carries the ENTIRE overprediction — and a large share of the
correlation.  The model plume's timing over 8021 is exactly right (obs peak
12.6 ng/m³ in our +9 h window), but the point evaluation reports it 20×
too strong: at 20 km range the parcel cloud is a ~200 m filament whose kNN
bandwidths (50–110 m) honestly resolve a plume the MODEL has not dispersed
enough (near-field relative dispersion is not yet represented — the
per-parcel covariance growth the puff machinery tracks would supply it).
The 0.1° grid scores "better" there by smearing the same filament to zero.
Interpretation: the receptor estimator separates MODEL near-field error
from REPRESENTATION error; the grid conflates them.  `min_bandwidth` in
`[output.receptors]` bounds the point-evaluation variance as an explicit
experiment (default 0 — width always from the trajectories).

Estimator calibration (closed-form, seconds — run before trusting any
score): `test/test_estimator_calibration.jl` — slab normalization,
prescribed-puff reproduction (±10% documented adaptive peak/flank tilt,
aggregate mass exact), the deterministic 8021 filament case, the real-ETEX
mapping deposit conservation (task #40), window semantics, 0–360/−180 lon
conventions.

## N-sweep + grid invariance (#31 wrap-up, 2026-07-07 overnight)

Same 72-h adaptive-kernel setup at five particle counts (0.1° grid,
`configs/etex_n*.toml`) plus a 0.2° twin at 400k (`etex_dx02.toml`); each run
scored on the same 2629 station-windows, gridded AND receptor
(`/tmp` sweep logs; receptor R(x_s) medians recomputed unrounded from the
receptors_*.csv files):

| N | grid R | grid FB | grid RANK | recept R | recept FB | recept NMSE | R(x_s) p50 / p90 |
|---|---|---|---|---|---|---|---|
| 50k | 0.320 | −0.067 | 2.56 | 0.232 | +0.82 | 1836 | 11.7 / 56 km |
| 100k | 0.315 | −0.051 | 2.54 | 0.158 | +0.64 | 1744 | 8.7 / 44 km |
| 200k | 0.322 | −0.073 | 2.55 | 0.181 | +0.24 | 588 | 5.6 / 30 km |
| 400k | 0.322 | −0.069 | 2.55 | 0.177 | +0.29 | 605 | 3.7 / 22 km |
| 800k | 0.321 | −0.066 | 2.56 | 0.190 | +0.15 | 374 | 2.5 / 14 km |
| 400k @ 0.2° | 0.526 | +1.09 | 2.22 | 0.164 | +0.33 | 731 | 3.7 / 21 km |

1. **The N^(−1/2) law holds on the real plume**: receptor-native resolution
   medians 11.7 → 2.5 km over ×16 trajectories (measured exponent −0.56 on
   the p50, −0.50 exactly on the p90).  At 800k the 0.1° grid (7–11 km
   cells) is fully resolved at the stations — resolution is bought with
   trajectories, not with the output grid.
2. **Gridded statistics are N-invariant from 50k up** (R 0.315–0.322,
   FB −0.05…−0.07, RANK 2.54–2.56): the kernel adapts its width to the
   sampling density, so the gridded product is already converged at low N —
   and the R gap vs the uniform kernel is N-independent, sealing the
   pair-level attribution above (a deterministic near-source transport miss,
   not estimator noise).
3. **Receptor FB/NMSE improve monotonically with N** (+0.82 → +0.15;
   1836 → 374): the near-source single-parcel spikes shrink as kNN h
   collapses.  (The dispersion bandwidth-floor was tried as an N-independent
   bound and REJECTED for the receptor — it over-counts a dense field; see
   the pair-level section.  The open fix is the √N-corrected per-parcel
   footprint.)
4. **Grid invariance holds exactly where representativeness allows**: at
   0.2° vs 0.1° the station-sampled fields correlate 0.95–0.98 beyond
   150 km with near-equal means and per-band R (0.77/0.48/0.52 vs
   0.73/0.42/0.52).  The pooled numbers diverge (R 0.526, FB +1.09) ENTIRELY
   inside 150 km: the 22-km cells bridge the ~20 km gap between the
   misplaced filament and station 8021 (band mean 9921 vs 0.6 pg/m³, obs
   387) — cell-size representativeness at a near-source station, not a
   deposit artefact.  The receptor estimator is unchanged across grids
   (R 0.164/0.177, FB +0.33/+0.29, R(x_s) 3.7 km both) — the
   model-error vs representation-error separation working as designed.

## Estimator wall-time benchmark (2026-07-07, 400k / 72 h / 6 threads)

Same machine, same met window, sequential runs; `[run] benchmark = true`
per-group tables at the end of each log (persisted copies in
`results/etex/bench_etex_{uniform,adaptive,meso}400k.log`).

| run (kernel) | TOTAL(run) | output group | output alloc | advection | notes |
|---|---|---|---|---|---|
| output_etex_uni (uniform) | 3053 s | **3.2 s** (0.1 %) | 61 MiB | 1857 s | bilinear deposit only |
| etex_meso400k (uniform + meso) | 2937 s | 2.8 s (0.1 %) | 62 MiB | 1813 s | OU adds ~nothing |
| output_etex_ab (adaptive kNN) | 2827 s | **522.9 s** (18.5 %) | **1.35 TiB** | 1746 s | kNN refresh + Gaussian deposit |
| etex_cov (covariance, physical) | 2382 s | **4.0 s** (0.2 %) | 71 MiB | 1808 s | P-deposit, NO kNN pass; `covkernel` 16.7 s (P evolution, 433 calls) |
| etex_cov_shape (covariance, shape) | 2338 s | 21.5 s (0.9 %) | 350 MiB | 1763 s | kNN size × P orientation (post-fix kNN); `covkernel` 15.9 s |

Reading it honestly:
- **TOTAL(run) is NOT the comparison metric** — run-to-run scatter in
  advection/metread (±5–8 %) is larger than the estimator cost (the adaptive
  TOTAL came out *lower* than uniform despite +520 s of output).  Compare the
  `output` / `covkernel` groups.
- **The adaptive-kNN output machinery costs ~520 s ≈ 9 min** against 3 s for
  uniform: 24 kNN bandwidth refreshes (once per 3-h frame, `bw_stride` = 6)
  plus 144 Gaussian-stencil deposits (30-min sample cadence).  The 1.35 TiB
  of allocations sit in the kNN pass (neighbour search at 400k).  This is the
  POST-fix number — before the ±32-cell stencil cap and the per-frame
  `bw_stride` the same runs took HOURS in the output loop.
- The two covariance runs complete the split: their `output` group is the
  Gaussian deposit WITHOUT any kNN pass, and `covkernel` isolates the
  P-evolution cost (per sync step, 433 calls) — together they attribute the
  520 s between neighbour search and stencil deposit.

**kNN allocation bug FOUND + FIXED (2026-07-07, user-spotted via this
table).**  Synthetic profiling (400k filament):
the adaptive DEPOSIT is allocation-free (0.09 s/call, 0 B) — the whole
1.35 TiB / ~520 s was the kNN pass, at 167 GiB / 69–94 s PER REFRESH.
`Profile.Allocs`: 97 % of bytes at conc.jl:243–257 (inner distance loop) as
boxed Float64/Int64 — `xmin`/`ymin` are reassigned in the bbox scan, the
`@threads` closure captures them as Core.Box, the cell index goes dynamic
and poisons j → x[j] (every op heap-allocates).  Fix: single-assignment
rebinding before the closure (conc.jl, 3 lines).  AFTER: **0.85 s /
0.012 GiB per refresh (~90× faster, ~14 000× fewer bytes, GC 9 %→0)**;
test_adaptive_kernel.jl 28/28 unchanged.  ETEX projection: adaptive output
~520 s → **≈ 35 s** (24×0.9 s kNN + 144×0.09 s deposit); the uniform run's
untimed receptor-kNN (~its 609 s group-sum gap) collapses identically.
The table rows above are PRE-fix numbers.

**Covariance-kernel accuracy (2026-07-08, both chains complete, 400k/72h,
2629 station-windows incl. zeros).**  Baselines from the 2×2 above: grid
uniform R 0.651 / FB +0.30 / FMS 56 % / RANK 2.74, adaptive R 0.322 / −0.07 /
58 % / 2.55; receptors uniform 0.642 / +1.07 / 2.36, adaptive 0.177 / +0.29 /
2.37.

| kernel | estimator | R | FB | NMSE | FMS | RANK |
|---|---|---|---|---|---|---|
| covariance physical | grid | 0.325 | −0.064 | 51.4 | 55.4 % | 2.54 |
| covariance physical | receptors | 0.601 | **+1.635** | **8608** | 54.1 % | **1.99** |
| covariance shape | grid | 0.323 | −0.067 | 50.9 | 59.7 % | **2.57** |
| covariance shape | receptors | 0.391 | +0.401 | 502 | 59.0 % | 2.43 |

Reading it:
- **On the grid both modes reproduce the adaptive-kNN field** (R ≈ 0.32,
  RANK 2.54–2.57 vs adaptive 2.55) — no gain, no harm.
- **The double-counting prediction is CONFIRMED on the receptor.**  Physical
  mode carries P as the *full* width, adding relative dispersion on top of the
  ensemble's absolute spread → receptor model-mean 1088 vs 109 pg/m³ observed
  (10× over), FB +1.64, NMSE 8608, RANK 1.99 — the worst of any estimator.
  Physical mode is falsified for point evaluation.
- **Shape mode is the best receptor estimator to date** (RANK 2.43 > uniform
  2.36 and adaptive 2.37; R 0.39 vs adaptive's 0.18, NMSE 502 vs 605), because
  it pins kernel *area* to the kNN spacing and takes only *anisotropy/tilt*
  from P.  It supersedes the σ(age)/√N_local receptor-footprint idea — but it
  STILL needs the kNN pass, so it does not deliver the "no-kNN" speedup that
  motivated physical mode.
- **Cost is negligible either way**: `covkernel` (P evolution) is 0.7 % of
  wall (~16 s of ~2380 s); both runs ≈ 40 min, on par with adaptive.
