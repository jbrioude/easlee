# EASLEE benchmarks

Two reference machines are used throughout this document.  They are named here so
that any table below can be read without knowing which one produced it.

| name | CPU | GPU | software |
|---|---|---|---|
| **workstation** | AMD Threadripper 3970X, 32 cores / 64 threads, 220 GB RAM | NVIDIA RTX 4080 16 GB | Debian 13, Julia 1.12 |
| **laptop** | 20-core mobile CPU | NVIDIA RTX A1000 6 GB | Julia with a functional CUDA |

Unless stated otherwise every run uses the same protocol: ETEX-1 ERA5 meteorology
(−10…30 °E × 40…70 °N), the PMCH release from Monterfil, forward transport, Thomson
turbulence, **convection off**, IFINE = 5, output grid 0.25°.  Before each timing the
sustained clocks were checked to be fully boosted (4.1–4.3 GHz all-core CPU,
2.6–2.8 GHz GPU SM), so none of the numbers below are deflated by a power-saving
governor.  Machines were otherwise idle and runs were sequential.

## 0. Workstation benchmark (2026-07-17): thread-scaling curve + RTX 4080

Protocol as in §1 (validate_captex.jl, NP=1000000 RUNDAYS=2 TURB=thomson CONV=off
IFINE=5 DXOUT=0.25 RELEASE=1).  Meteorology: a fresh MARS retrieval 19941023–27
with the surface companions.  Observations: NOAA DATEM, with emit-t1.txt
reconstructed from the documented release (340 kg, 48.058 °N, −2.0083 °E,
16:00 UTC + 11 h 50 m).  Wall clock is for the whole run including Julia startup.

| backend | threads | wall clock | speedup vs 1 t | notes |
|---------|--------:|-----------:|---------------:|-------|
| `threads` | 1  | 1:51:24 (6684 s) | 1.00× | single core @ 4.3 GHz |
| `threads` | 2  | 57:58 (3478 s) | 1.92× | |
| `threads` | 4  | 32:15 (1935 s) | 3.45× | |
| `threads` | 8  | 20:36 (1236 s) | 5.41× | |
| `threads` | 16 | 12:10 (730 s)  | 9.16× | laptop 16 t = 27:30, workstation 2.26× faster |
| `threads` | 32 | 7:49 (469 s)   | 14.3× | all physical cores; parallel eff. 45 % |
| `threads` | 64 | 6:22 (382 s)   | 17.5× | SMT step +1.23× |
| `ka_gpu` (RTX 4080) | 4 | 7:44 (464 s) | — | first run, incl. CUDA kernel compile |
| `ka_gpu` (RTX 4080) | 4 | **4:54 (294 s)** | 22.7× | repeat (warm cache); **651 MiB / 16 GB**, SM 2.79 GHz, util 100 % peak |

Reading:
- **RTX 4080 = 2.48× the 16-thread CPU** (294 vs 730 s) and **1.30× the best CPU
  row** (64 t, 382 s).  vs the laptop A1000 (15:24) the 4080 whole-run is **3.1×**
  — far below the ~7× FP32 FLOP ratio, as expected: the shared CPU work (GRIB
  met reads, deposit, NetCDF, release) and host↔device copies bound the whole-run
  speedup by Amdahl (§1 headroom items unchanged).
- **Scaling curve** flattens exactly as the cost-centre split predicts: near-ideal
  to 4 t (1.92×/1.80× per doubling), ~1.6× per doubling in the 8–32 t range as the
  serial met/IO fraction grows, +23 % from SMT.
- **Scores identical across all CPU rows** (POOLED N=1720, R 0.486–0.487, RANK
  2.22) — same physics, only wall clock varies.  `ka_gpu` scores RANK 2.10/2.11,
  the same ka-vs-threads dispersal offset documented in §2 (laptop: 2.04 vs 2.15).
  Note absolute scores differ slightly from the laptop table (RANK 2.22 vs 2.15).
  Controls rule out every input-side cause:
  - **Met**: laptop subsample (6 files, 19941027 1800 steps 7–9) **bit-identical**
    to the workstation's MARS retrieval (max|Δ|=0 on 8325 fields; laptop-only extra `chnk`
    unused).  **Emit**: laptop emit-t1.txt value-identical to the reconstruction.
  - **RNG/sampling**: 3 identical 16 t runs → RANK 2.21–2.22, R ±0.000, FMS ±0.3
    (and all 7 thread counts scored 2.22) → sampling spread ±0.005, 10× too small.
  - **12-bit jpeg repack** of the meteorological archive: RANK 2.22,
    **no score impact**; met-read cost +26 s of 756 s at 16 t (~3.5 %); 2.8×
    smaller on disk.  The archive recipe is safe.
  - **Code drift**: the tree that recorded the laptop table, re-run on the
    workstation, also gives RANK 2.21.  Not the code.
  The shift is systematic and backend-independent (+0.06–0.07 on threads AND
  ka_gpu), and is most probably a Julia-version or hardware systematic between
  the two machines.  The timing columns are unaffected by it.

## 0c. v1.1 transport ports — ETEX regression (2026-09-13)

Same config (`configs/output_etex.toml`, 400 k parcels, 72 h, 16 threads,
Thomson, massflux), scored with `validation/check_etex.jl` (the pooled
2629-station-window scorer of ETEX_RESULTS.md — NOT the validation-task RANK of
the §0/§0b tables, which is a different scorer and not comparable):

| tree | R | FB | FA2 | FA5 | FMS | KSP | RANK | wall |
|------|--:|---:|----:|----:|----:|----:|-----:|-----:|
| documented baseline (ETEX_RESULTS.md, framed engine) | 0.651 | +0.309 | 0.32 | 0.65 | 57.2 % | 9.3 | **2.75** | — |
| v1.1 shared | 0.660 | +0.326 | 0.32 | 0.66 | 55.8 % | 8.3 | **2.75** | 11.6 min |

RANK reproduced to the last digit; the components move within run-to-run RNG
scatter.  Expected: none of the v1.1 ports touches an ERA5 (Hanna-fill,
lon/lat-grid, single-nest) run — `horizontal_tke` is off and would be refused
on ERA5 anyway, the per-nest closure collapses to the run-global one on a
single nest, the projected-grid axes only engage on `crs = "native"`, and the
Hanna runaway guard did not fire (no warning in the run log).  What DID change
for this run is the release-slot cursor and the release step count, which are
allocation order only.

## 0b. Workstation re-run after the September 2026 kernel fixes (2026-09-05)

Same machine, same protocol and knobs as §0 (validate_captex.jl, NP=1000000
RUNDAYS=2 TURB=thomson CONV=off IFINE=5 DXOUT=0.25 RELEASE=1), mid-run CPU clocks
4.0–4.1 GHz on the loaded cores, GPU SM 2.82 GHz at 100 % utilisation, machine
otherwise idle, runs sequential.  New knob `MO=false` switches the
Monin-Obukhov surface-layer wind off so that threads and `ka_gpu` run the
SAME physics (the KA kernel does not implement it; since 2026-09-04
`advance_dispatch!` would otherwise force the threads path).

| backend | threads | MO | wall clock | RANK | July (§0) |
|---------|--------:|:--:|-----------:|-----:|-----------|
| `ka_gpu` (RTX 4080), cold | 4 | off | 290 s | 2.19 | 464 s / 2.11 |
| `ka_gpu` (RTX 4080), warm | 4 | off | **282 s** | 2.20 | 294 s / 2.11 |
| `threads` | 16 | off | 607 s | 2.21 | — |
| `threads` | 32 | off | 408 s | 2.21 | — |
| `threads` | 16 | on (production) | **606 s** | 2.20 | 730 s / 2.22 |

Reading:
- **Speed-up held**: RTX 4080 = **2.15× the 16-thread CPU** and **1.45× the
  32-thread CPU** with identical physics (July: 2.48× vs 16 t, but against a
  CPU row that carried the MO work the GPU skipped).  The kernel now runs the
  Hanna/Thomson turbulence it used to drop (A2) and the entry/NaN guards (A4),
  and the met stencil is uploaded only when it changes (B1) instead of every
  step: the two effects roughly cancel on the GPU row (294 → 282 s), so the
  whole-run time is still bounded by the shared CPU work (met read, deposit,
  output, release) and the per-step state round trip.
- **CPU rows 17 % faster than July** (730 → 606 s at 16 t) with unchanged
  scores: the Monin-Obukhov correction is now evaluated only below the lowest
  model level instead of two extra 4-D interpolations per parcel per sub-step
  everywhere (B2a); MO on vs off is now free (606 vs 607 s).
- **Parity offset closed**: `ka_gpu` scores 2.19–2.20 vs threads 2.20–2.21,
  within the ±0.005 sampling spread measured in §0 (July: 2.11 vs 2.22, the
  "ka disperses more" offset of §2, which was the (ctl+1)/2 vertical-wind bug
  A1 plus the missing turbulence A2).  `ka_*` can now serve as a regression
  oracle for the threads path when the CPU-only switches are off.
- Cold vs warm GPU: 8 s (kernel compile cache).  GPU memory 619 MiB / 16 GB.

## 0d. Benchmark matrix for the paper (workstation, 2026-09-23; supersedes §0 where they disagree)

Same protocol as §0 (ETEX, 1 M particles, 2 days, 15-min sync steps, Thomson turbulence,
convection off, `surface_layer_wind = false`, uniform output kernel), run by
`scripts/bench_matrix.sh` on an otherwise idle machine with the CPU fully boosted, one
row after the other, each in its own batch unit.  The **GPU rows were re-measured on the
Val-switch fix** (see "the GPU regression" below; the first GPU pass read
295 / 293 / 210 / 211 s).  Numbers are **TOTAL(run)** from
`[run].benchmark = true` — run only, no Julia start-up (start-up is the difference to the
"unit wall" column of the raw summary, 60–90 s). The tables are regenerated by
`scripts/bench_matrix_table.py`. `PREFETCH=1` sets the met prefetch loader; `-t 16,1`
gives it an interactive thread (rows `gpu_t16i*`).

### Table A — every row (run-only seconds; ETEX 1 M parcels, 2 days; GPU rows and the 32-thread pair re-measured on the Val-switch fix)

| row | backend | threads | prefetch | TOTAL | metread | advection | output kernel | release | untimed | vs 1 t | vs 16 t |
|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| cpu_t1 | threads | 1 | off | 5910 | 94.3 | 5795.4 | 9.8 | 6.6 | 4.2 | 1.00× | 0.10× |
| cpu_t2 | threads | 2 | off | 3055 | 84.2 | 2949.4 | 9.7 | 7.1 | 4.2 | 1.93× | 0.19× |
| cpu_t4 | threads | 4 | off | 1683 | 109.4 | 1551.7 | 9.9 | 7.5 | 4.3 | 3.51× | 0.34× |
| cpu_t8 | threads | 8 | off | 975 | 105.5 | 847.6 | 10.0 | 7.9 | 4.2 | 6.06× | 0.59× |
| cpu_t16 | threads | 16 | off | 576 | 108.5 | 444.5 | 10.1 | 8.6 | 4.4 | 10.26× | 1.00× |
| cpu_t32 | threads | 32 | off | 374 | 107.7 | 243.2 | 10.1 | 8.5 | 4.4 | 15.81× | 1.54× |
| cpu_t64 | threads | 64 | off | 316 | 116.9 | 175.2 | 10.2 | 9.6 | 4.4 | 18.68× | 1.82× |
| cpu_t16_pf | threads | 16 | on | 558 | 14.0 | 512.4 | 10.8 | 15.8 | 4.5 | 10.60× | 1.03× |
| cpu_t32_pf | threads | 32 | on | 360 | 14.1 | 317.2 | 10.6 | 13.4 | 4.5 | 16.42× | 1.60× |
| cpu_t64_pf | threads | 64 | on | 295 | 14.3 | 253.6 | 10.4 | 12.0 | 4.6 | 20.04× | 1.95× |
| gpu_t4 | ka_gpu | 4 | off | 249 | 109.3 | 117.7 | 10.2 | 7.8 | 4.5 | 23.70× | 2.31× |
| gpu_t16i | ka_gpu | 16 | off | 245 | 104.8 | 117.3 | 10.3 | 8.1 | 4.5 | 24.12× | 2.35× |
| gpu_t4_pf | ka_gpu | 4 | on | 168 | 16.8 | 120.0 | 10.4 | 16.1 | 4.5 | 35.22× | 3.43× |
| gpu_t16i_pf | ka_gpu | 16 | on | 165 | 14.4 | 118.9 | 10.5 | 16.5 | 4.8 | 35.80× | 3.49× |

### Table B — per process: time [s] and speedup (vs 1 thread; in brackets vs 16 threads)

| process | 1 t | 16 t | 64 t | GPU (4 t) | GPU + prefetch | scaling exponent (CPU, no prefetch) |
|---|---:|---:|---:|---:|---:|---:|
| met fields read + process | 94.3 | 108.5 (0.87×) | 116.9 (0.81×; [0.93×]) | 109.3 (0.86×; [0.99×]) | 16.8 (5.61×; [6.45×]) | t ∝ N^+0.06 (7 pts, 1–64 t) |
| advection + turbulence | 5795.4 | 444.5 (13.04×) | 175.2 (33.07×; [2.54×]) | 117.7 (49.24×; [3.78×]) | 120.0 (48.31×; [3.70×]) | t ∝ N^-0.86 (7 pts, 1–64 t) |
| output kernel (deposit_to_grid!) | 9.8 | 10.1 (0.97×) | 10.2 (0.96×; [0.99×]) | 10.2 (0.96×; [1.00×]) | 10.4 (0.95×; [0.98×]) | t ∝ N^+0.01 (7 pts, 1–64 t) |
| release | 6.6 | 8.6 (0.76×) | 9.6 (0.68×; [0.89×]) | 7.8 (0.85×; [1.11×]) | 16.1 (0.41×; [0.53×]) | t ∝ N^+0.09 (7 pts, 1–64 t) |
| untimed (frame dump, reassign, sync) | 4.2 | 4.4 (0.97×) | 4.4 (0.96×; [0.99×]) | 4.5 (0.94×; [0.97×]) | 4.5 (0.93×; [0.96×]) | t ∝ N^+0.01 (7 pts, 1–64 t) |
| TOTAL(run) | 5910.3 | 576.0 (10.26×) | 316.4 (18.68×; [1.82×]) | 249.4 (23.70×; [2.31×]) | 167.8 (35.22×; [3.43×]) | t ∝ N^-0.72 (7 pts, 1–64 t) |
| convection | off in the ETEX protocol — see the full-physics table (2026-07-09): 4.7× at 6 threads | | | | | |
| dry / wet deposition | off in the ETEX protocol — full-physics table: wet dep 1.5× at 6 threads | | | | | |


### Reading

- **Thread scaling.** The advance (advection + turbulence, the only parallel cost centre)
  scales as N^−0.86 over 1–64 threads: 13.0× at 16 t, 33× at 64 t (32 physical cores,
  parallel efficiency 49 % at 32 t, +1.39× from SMT 32 → 64). The whole run scales as
  N^−0.72 because everything else is serial: met read 94–117 s, output kernel 10 s, release
  7–10 s, frame dump and bookkeeping 4 s — a **~125 s floor** that is 2 % of the run at
  1 thread and 40 % at 64. The met read even slows slightly with threads (94 → 117 s,
  N^+0.06): it allocates 14 GiB over 193 calls and the GC has more to stop.
- **GPU.** The RTX 4080 advance is 118 s against 444 s on 16 threads (**3.8×**) and
  175 s on 64 threads (**1.5×**) — one consumer GPU beats the full 32-core CPU on the
  parallel part. Whole run: 2.3× vs 16 threads, 1.27× vs 64. Giving the prefetch loader an
  interactive thread (`-t 16,1`) changes nothing on a quiet machine (the 16 % of 21 Sept was
  measured under load).
- **Prefetch.** On the CPU it hides the met read (108 → 14 s) but moves the decode into
  the advance's wall time (+68 s at 16 t, +74 s at 32 t, +78 s at 64 t: the loader competes
  with the workers for cores and memory bandwidth), net −3 % / −4 % / −7 %. The 32-thread
  pair was re-measured as a clean pair on the idle machine: 374 s
  (reproducing the matrix row to 0.1 %) against 360 s with prefetch. On the GPU the CPU is idle while the
  kernel runs, so the read is hidden for free: **168 s, −33 %, the best row** — 3.4× the
  16-thread CPU and 1.76× the best CPU row (64 t + prefetch, 295 s). The remaining GPU run
  is 71 % advance (kernel + host↔device copies) and 29 % serial floor (met read 17 s,
  release 16 s, output kernel 10 s, untimed 5 s).
- **The GPU regression, found by this matrix and fixed the same morning.** The first GPU
  pass read 167 s in the advection group where the same config had given 128–137 s on
  21 September (under load) and 294 s whole-run in July. A three-point bisect on the quiet
  machine: before the five-switch port **116.3 s**, the port itself **167.5 s**, head
  **167.0 s**, so the whole +44 % is the port and the later KA fixes cost nothing. Cause: the switches (`use_hanna`, `mesoscale`, `kh_ft`,
  `horizontal_tke`, `horizontal_hanna`, `surface_layer_wind`) reached the kernel as runtime
  `Bool`s; the GPU compiler specialises on argument TYPES, so every branch was compiled
  into the one kernel and paid in registers and occupancy even with all five off. Fix
  (`advance_ka.jl`, 23 Sept): the switches are passed as `Val{Bool}` and unwrapped with
  `_kb`, so each configuration compiles its own kernel with the dead branches eliminated.
  Re-measured: **117.7 s** (`gpu_t4`), i.e. the pre-port number; unit tests
  `test_ka_switches_unit.jl` 5/5 and `test_hanna_gpu_4d.jl` 12/12; the all-five
  threads-vs-GPU validation and the well-mixed 4-D tests were rerun on the fix (results
  in §3). The lesson for the paper's method section: **a KA kernel must receive physics
  switches as types, not values.**
- **Against the July table (§0)** the tree is not slower on the CPU (1 t 5979 s wall vs
  6684; 16 t 661 vs 730; 64 t 390 vs 382) and the fixed GPU row is faster (321 s wall vs
  294 + start-up; run-only 249 s).
- **What this protocol does not measure**: convection and deposition are off in ETEX (the
  2026-07-09 full-physics table below has them at 6 threads: 4.7× and 1.5×), the output
  kernel here is the uniform one (the adaptive kernel is the super-linear cost centre of
  `BENCHMARKS_COSTCENTRES.md` and has no GPU path), and the NetCDF frame dump is inside the
  4 s "untimed" remainder, not a group of its own. For the chemistry runs the picture is
  different again: the CB05/strat solver is ~70 % of the wall clock and GC-bound
  (`BENCHMARKS_COSTCENTRES.md`).

## 1. Advance backend — GPU vs CPU (performance)

The parcel advance can run on three backends, selected by `EASLEE_ADVANCE`:
`threads` (default, multithreaded CPU), `ka_cpu` (KernelAbstractions on CPU),
`ka_gpu` (same kernel on `CUDABackend()`). Wiring: `Operators.advance_dispatch!`.

**1 000 000 particles, 2-day transport** (wall clock incl. ~40 s Julia startup):

| backend | threads | wall clock | speedup vs CPU | GPU util | GPU mem |
|---------|--------:|-----------:|---------------:|---------:|--------:|
| `threads` (CPU) | 16 | **27:30** | 1.0× | — | — |
| `ka_gpu` (RTX A1000) | 4 | **15:24** | **1.78×** (1.82× net of startup) | peak 100 % | 432 MiB / 6 GB |

Notes:
- The GPU win is **wall-clock for the whole run**, not advance-only. Per-step the
  GPU also pays met + state host↔device copies, and the shared CPU work (GRIB met
  reads, `deposit_to_grid!`, NetCDF output, release) is unchanged — so by Amdahl
  the total speedup (~1.8×) is well below the advance-kernel speedup alone.
- **Memory is a non-issue**: the entire ERA5 met stencil + 1 M parcels used only
  **432 MiB of 6 GB**. Much larger particle counts / domains fit easily.
- **Headroom** (not yet done): upload met only on met-refresh (not every step),
  gather only the *active* parcels instead of the full arrays, and deposit on the
  GPU (`output/conc_ka.jl` already exists). These target the copy + CPU-deposit
  bottlenecks and should push the speedup well past 2×.

## 2. Backend correctness / parity

**ETEX-1, 50 000 particles, 1-day** (stochastic — independent RNG per backend):

| backend | N | R | FB | FAC2 | FMS | KSP | RANK |
|---------|--:|--:|--:|--:|--:|--:|--:|
| `threads` (CPU, production) | 617 | 0.539 | +1.44 | 0.19 | 46.1 | 7.6 | 1.96 |
| `ka_cpu` | 617 | 0.521 | +1.29 | 0.17 | 35.0 | 9.2 | 1.88 |
| `ka_gpu` | 617 | 0.516 | +1.29 | 0.16 | 36.3 | 9.6 | 1.89 |

**1 000 000 particles, 2-day** (the perf-benchmark runs):

| backend | N | R | FB | FAC2 | FMS | KSP | RANK |
|---------|--:|--:|--:|--:|--:|--:|--:|
| `threads` | 1720 | 0.476 | +1.08 | 0.25 | 57.5 | 11.2 | 2.15 |
| `ka_gpu`  | 1720 | 0.432 | +0.76 | 0.16 | 40.0 | 16.4 | 2.04 |

Reading the parity:
- **`ka_gpu` ≡ `ka_cpu`** (same kernel) — confirmed: FB identical, RANK 1.89 vs 1.88.
  The GPU implementation is numerically correct relative to the KA kernel.
- **`ka_*` ≈ `threads` but not identical** (tables above, pre-2026-09-04 code).
  The ~0.1 RANK offset was traced to two KA-kernel bugs (vertical mean wind
  applied (ctl+1)/2 times; default Hanna turbulence run as no turbulence) and
  is closed: see §0b (2.19–2.20 vs 2.21 at 1 M parcels) and
  `test/test_ka_parity.jl`.  The validated threads path remains the reference;
  `ka_*` now reproduces it when the CPU-only switches (§3) are off.

## 3. Known limitations of the KA / GPU advance

- **Turbulence**: Thomson, `Hanna(false)` (pre-computed σ_w/T_L fields — the
  ERA5/GFS/MERRA-2/WRF default, run through the same Thomson body as the CPU)
  or none.  Hanna on-the-fly (`turbswitch = true`), the mesoscale OU term,
  `kh_ft > 0` and the Monin-Obukhov surface layer are `threads`-only:
  `advance_dispatch!` forces the threads path and warns once naming the switch
  (before 2026-09-04 the KA kernel silently dropped them).
- **IFINE is hard-coded to 5** in the KA wrapper (ignores the config / `IFINE` env).
- **Convection** is unaffected: it is a separate operator-split step applied after
  the advance regardless of backend, so `ka_*` runs still include mass-flux
  convection.
- The GPU path keeps one resident buffer set **per nest** and uploads the 8 met
  arrays only when they changed — detected from the met bracket
  (`nestrunconfig.t_index`, written by `readmetfiles!`/prefetch when a slot is
  refilled) plus a sampled fingerprint (~1024 strided elements per array, bit
  hash) that catches in-place edits leaving `t_index` alone.  The 12 state
  arrays are still uploaded/downloaded every step (gather of the active set is
  the remaining headroom item of §1).
- **Per-parcel RNG**: seeded from (slot, first-advance time, run seed) through
  splitmix64; `EASLEE_SEED=<integer>` changes the ensemble (default: fixed
  constant → reproducible).  Before 2026-09-04 the seed was the slot index
  alone, so successive occupants of a slot replayed identical noise.

### KA-kernel review fixes (2026-09-04) — `test/test_ka_parity.jl`

Synthetic-met parity test (threads vs `ka_cpu`, and `ka_gpu` with
`EASLEE_GPU=1`), all asserting against the CPU path:

| fix | symptom before | check |
|-----|----------------|-------|
| A1 vertical wind applied (ctl+1)/2 times (accumulated `dlev` re-added each sub-step) | 3x the descent at ctl=5 | per-parcel `indl` identical to threads (max\|Δ\| = 0) |
| A5 Float32 `while inc_time < Δt` sub-step loop ran ctl+1 sub-steps for ~20 % of (Δt, ctl) | 10/9 at Δt=3 s, ctl=9 | Δindl = W·Δt within 0.5 %, identical to threads |
| A2 `Hanna(false)` demoted to NoTurb | zero vertical spread in the PBL | std(Δindl) ka/threads = 1.02–1.04 (1 km PBL, 2 h) |
| A4 no entry / NaN guards | OOB read (`unsafe_trunc` on index 0), NaN positions | OOB parcels untouched with Dtime held; NaN met cell → zero displacement |
| A7 slot-only RNG seed | replayed noise on slot reuse | reused slot gets a different stream |
| A5 met re-uploaded every step, one buffer set for all nests | copy bound, realloc per step on nests | in-place met change detected on the GPU (Δ doubles) |

### 3b. Validation record for the ported switches (2026-09-23)

- **2026-09-21: NOTHING is threads-only any more.** `mesoscale`, `kh_ft`,
  `horizontal_tke`, `horizontal_hanna` and `[physics] surface_layer_wind` all run
  on the KA path (`advance_ka.jl`; the mesoscale OU state joins the resident set,
  so the device state went 12 → 15 arrays). `_ka_unsupported_switch` returns
  nothing; the two `build_run.jl` refusals for horizontal_tke/horizontal_hanna are
  retired. Validated with ALL FIVE ON AT ONCE (ETEX, 200 k particles, 12 h, each
  leg's backend confirmed from its own log):

  | comparison | total mass | correlation |
  |---|---|---|
  | ka_cpu vs ka_gpu | 1.09e-3 apart | 0.999963 |
  | threads vs ka_cpu | 3.18e-4 apart | 0.999951 |
  | threads vs ka_gpu | 1.40e-3 apart | 0.999962 |

  All three agree to about one part in a thousand in mass and 0.99995 in pattern —
  the spread between backends is the ensemble's own Monte-Carlo noise, not a
  physics difference.
  **Re-validated 2026-09-23 on the Val-switch fix** (the switches now reach the kernel
  as types, §0b): threads vs ka_gpu, all five on, ETEX 200 k / 12 h, each leg's backend
  confirmed from its log, **2.34e-3 apart in mass, correlation 0.999935**
  (`scripts/ka_compare_grids.py`, outputs `bench_all5_{threads,ka_gpu}`).
- **Well-mixed criterion on the device-side Hanna (2026-09-22,
  `test/test_wellmixed_pbl_4d.jl`)** — the physics test for a Langevin scheme,
  independent of any backend comparison: 100 k parcels started uniform in air mass
  over 0–4 km, 24 h of turbulence only. Max deviation from uniform 4.3 % (convective,
  u* 0.35, h 1200 m, L −50), 3.8 % (stable, u* 0.15, h 300, L 80), 4.3 % (neutral,
  u* 0.30, h 800), against a 10.9 % tolerance; the PBL air-mass fraction is
  conserved to 0.1–0.3 % in every case. No pile-up or void at the PBL top.
- **Per-switch validation (2026-09-22, each switch ALONE, threads vs ka_gpu, ETEX
  200 k / 12 h):** `surface_layer_wind`, `mesoscale`, `kh_ft` and `horizontal_hanna`
  each ran on both backends with correlation 0.99987–0.99990. **`horizontal_tke`
  could NOT be validated**: `build_run` refuses it on ETEX because it requires a
  native-e driver (AROME/MesoNH `tke_kind = e`) and ERA5 archives a Hanna-filled
  0.5σ_w² — so that one switch is ported but untested end-to-end until an AROME
  configuration is run. Also observed in these pairs: a systematic ~1.4–1.9 % mass
  excess on the KA side relative to `threads` (same sign in all four); see the
  control note below for whether it predates the ports.
- **Two pre-existing KA bugs found by the 2026-09-22 double-check, both fixed and
  validated by the well-mixed criterion (`test/test_wellmixed_thomson_4d.jl`,
  `test/test_wellmixed_pbl_4d.jl`):**
  1. *The top model layer was unreachable on the KA path.* `interpolation4D_irregzsingle`
     looped `Nkmin:(Nkmax-1)`, never bracketing the top layer, so every height inside
     it fell through to the fallback and was placed on the layer's lower boundary
     (the CPU's `interpolation1Dsingleirreg` loops `2:Nz`). Well-mixed test on a
     4 km column: KA Thomson **14.4 %** max deviation with a −20 % void in the top
     bin and zero parcels above z(Nk) − 50 m, in all three stability regimes; CPU
     3.0–3.7 %. Fixed (loop to `Nkmax`): **KA 4.9 %**, PASS everywhere; the
     device-side Hanna went from 3.8–4.3 % to 2.3–4.3 %.
  2. *Inconsistent Langevin drift term in the KA Thomson body.* It differenced the
     CLAMPED σ_w at the bracketing levels; the CPU had been corrected to
     (wshare/σ_w)·∂tke/∂z (the form consistent with σ_w² = 2·wshare·tke) and the KA
     twin had not. Aligned.
  **The remaining KA-vs-threads mass offset is pre-existing, bounded and undiagnosed.**
  All-off control (pure advection + Thomson fill, ETEX 200 k / 12 h), KA relative to
  threads: old code **+1.72 %**; + drift fix +1.61 %; + top-layer fix **+1.50 %**. It is
  EARLY-TIME, not cumulative (per 3-h frame: +2.0–2.3 % in frames 2–3, +0.9 % in
  frame 4), and `ka_cpu ≈ ka_gpu` (0.18 % apart), so it is algorithmic, not a device
  effect. Ruled out: the KA random generator (mean +0.0004, std 0.998, kurtosis 2.92 —
  identical to `clamp(randn(), ±3)`), active-set selection, escape/kill logic, and
  new-parcel initialisation (the GPU initialiser is never called; both paths use the
  CPU one). Next step is instrumentation — dump the first steps of identical parcels on
  both paths — not more reading. FB/RANK scores never saw it: they are insensitive to a
  2 % mass offset.
- **A requested backend is now honoured or refused, never swapped.** If
  `EASLEE_ADVANCE` is set explicitly and anything would force a different path
  (an unported switch, 2-D/1-D wind, CUDA not functional), the run ERRORS. Only the
  unset default still falls back. Three separate silent-downgrade sites existed
  before this (`_ka_unsupported_switch`, two guards in `build_run.jl`, and the
  CUDA check) — each found by tripping over it.
- **Turbulence**: Thomson, none, and — since 2026-09-21 — **Hanna on the fly**
  (`turbulence_hanna_4d`, the device twin of `turbulence_hanna_v1`; the surface
  parameters ustar/pblh/ol join the resident met set). Validated on ETEX,
  200 k particles, 12 h, both legs confirmed on the KA path by their log lines:
  total mass 1.05e-3 apart, correlation 0.999964, max difference 6.3e-3 of the
  field maximum. Not bitwise — FMA contraction and different `exp`/`sqrt`/`cbrt`
  implementations move the last bit, and a 12-hour Langevin trajectory amplifies
  that into individual parcels landing elsewhere, as it already does for Thomson.
- **IFINE is hard-coded to 5** in the KA wrapper (ignores the config / `IFINE` env).
- **Convection** is unaffected: it is a separate operator-split step applied after
  the advance regardless of backend, so `ka_*` runs still include mass-flux
  convection.
- The current GPU path uploads met + full state each step (correctness-first); see
  §1 "headroom" for the optimisations that remain.


## How to reproduce

Two TOMLs now carry the protocol, so the benchmark needs no observation files
and no environment variables beyond the backend switch (`[run].threads` reaches
the child julia as `-t` through scripts/easlee.jl):

```bash
EASLEE_ADVANCE=ka_gpu  julia --project=easlee scripts/easlee.jl configs/bench_gpu_etex_1M.toml
EASLEE_ADVANCE=threads julia --project=easlee scripts/easlee.jl configs/bench_cpu_etex_1M.toml
```

Run each arm twice and quote the second (the first GPU run pays the CUDA kernel
compile). `[run].benchmark = true` prints the per-group split; `advection` is
the group the backend changes, `TOTAL(run)` is the number quoted above.

**Three settings in those TOMLs exist only to keep the arms on the same
physics** — the KA kernel does not implement them and `advance_dispatch!` falls
back to threads if any is on, which would silently benchmark threads against
threads: `physics.surface_layer_wind = false` (schema default is **true**),
`physics.turbulence.scheme = "thomson"` (Hanna on-the-fly is threads-only), and
`output.kernel = "uniform"` (adaptive/covariance refuse ka_gpu). Also
threads-only: `mesoscale`, `kh_ft > 0`, `horizontal_tke`, 2-D/1-D wind masks.
Verify in the log that `advance backend selected mode = :ka_gpu` appears AND
that no "overriding EASLEE_ADVANCE=... with threads" warning follows it — the
mode line is printed before the override, so it alone proves nothing (checked
2026-09-18: with `surface_layer_wind = true` the override fires).

The §0/§1 rows above predate `surface_layer_wind` (added 2026-08-13), so both
arms then ran without it — which is what these TOMLs reproduce.

Do not shrink the run to make it quick: at small N the GPU loses. A 3-h,
20 000-parcel smoke test on the RTX 4080 spent 12.9 s in advection against 2.5 s
on 8 CPU threads (kernel compile + per-step host<->device copies dominate).

The original env-var route still works:

```bash
CAPTEX_OBS=<ETEX obs dir> CAPTEX_MET=<ETEX met dir> EASLEE_ADVANCE=ka_gpu MO=false \
  TURB=thomson CONV=off IFINE=5 MODE=run DIRECTION=forward RELEASE=1 \
  NP=1000000 RUNDAYS=2 DXOUT=0.25 julia -t 4 --project=easlee validation/validate_captex.jl
```

(`MO=false` was added to validate_captex_run.jl for exactly the reason above; it
is not in the July command line because the switch did not exist yet.)

A standalone copy of these instructions, paired with the meteorology, is in
`docs/README_GPU_BENCHMARK.md`.
