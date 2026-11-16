# EASLEE benchmarks

Machine: laptop, 20-core CPU, **NVIDIA RTX A1000 6 GB** (Laptop), Julia, CUDA functional.
All runs below: ETEX-1 ERA5 met (−10…30 °E × 40…70 °N), PMCH release from Monterfil,
forward, Thomson turbulence, **convection off**, IFINE=5, output grid 0.25°.

## 0. Server benchmark — zion (2026-07-17): thread-scaling curve + RTX 4080

Machine: **zion**, AMD Threadripper 3970X (32 c / 64 t), **RTX 4080 16 GB**, 220 GB RAM,
Debian 13, Julia 1.12.6.  Governor **schedutil** (no sudo → performance mode not
settable); measured sustained clocks 4.1–4.3 GHz all-core CPU, 2.6–2.8 GHz GPU SM —
i.e. fully boosted, timings NOT deflated ([[benchmark-governor-check]] satisfied).
Same protocol as §1 (validate_captex.jl, NP=1000000 RUNDAYS=2 TURB=thomson CONV=off
IFINE=5 DXOUT=0.25 RELEASE=1); met = fresh MARS retrieval 19941023–27 (same params
as the 2025 archive pipeline + SFC companions; shortNames avg_* verified), obs =
NOAA DATEM (emit-t1.txt reconstructed from documented release: 340 kg, 48.058 °N
−2.0083 °E, 16:00 UTC + 11 h 50 m).  Whole-run wall clock incl. Julia startup;
machine otherwise idle; runs sequential.  Logs: ~/zion_bench/ on zion.

| backend | threads | wall clock | speedup vs 1 t | notes |
|---------|--------:|-----------:|---------------:|-------|
| `threads` | 1  | 1:51:24 (6684 s) | 1.00× | single core @ 4.3 GHz |
| `threads` | 2  | 57:58 (3478 s) | 1.92× | |
| `threads` | 4  | 32:15 (1935 s) | 3.45× | |
| `threads` | 8  | 20:36 (1236 s) | 5.41× | |
| `threads` | 16 | 12:10 (730 s)  | 9.16× | laptop 16 t = 27:30 → zion 2.26× faster |
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
  Controls run on zion 2026-07-17 rule out every input-side cause:
  - **Met**: laptop subsample (6 files, 19941027 1800 steps 7–9) **bit-identical**
    to the zion MARS retrieval (max|Δ|=0 on 8325 fields; laptop-only extra `chnk`
    unused).  **Emit**: laptop emit-t1.txt value-identical to the reconstruction.
  - **RNG/sampling**: 3 identical 16 t runs → RANK 2.21–2.22, R ±0.000, FMS ±0.3
    (and all 7 thread counts scored 2.22) → sampling spread ±0.005, 10× too small.
  - **12-bit jpeg repack** (script_era5.sh recipe, /data/ETEX_JPEG): RANK 2.22 —
    **no score impact**; met-read cost +26 s of 756 s at 16 t (~3.5 %); 2.8×
    smaller on disk.  The archive recipe is safe.
  - **Code drift**: worktree at 4d84253 (the commit that recorded the laptop
    table, pre-719c3ec k=Nk stencil fix) → RANK 2.21 on zion.  Not the code.
  The shift is systematic and backend-independent (+0.06–0.07 on threads AND
  ka_gpu).  OPEN: rerun the §How-to-reproduce command verbatim on the laptop —
  if it now yields ≈2.22, the 2.15/2.04 rows came from different knobs (or
  uncommitted state) at record time; if ≈2.15, it is a real Julia
  1.12.5-vs-1.12.6 / hardware systematic.  Timing columns are unaffected.

## 0b. Benchmark matrix for the paper — zion, 2026-09-23 (supersedes §0 where they disagree)

Same protocol as §0 (ETEX, 1 M particles, 2 days, 15-min sync steps, Thomson turbulence,
convection off, `surface_layer_wind = false`, uniform output kernel), run by
`scripts/bench_matrix.sh` on a QUIET machine (no other julia; the chemistry runs had
finished 3 h earlier), CPU rows at commit `1c75d62` (KA fixes `8a12d39` included), **GPU
rows re-measured the same morning on the Val-switch fix** (see "the GPU regression"
below — the first GPU pass at `1c75d62` read 295 / 293 / 210 / 211 s, all four rows
kept in `bench_matrix_20260922/`), governor
`schedutil` (fully boosted, trap #7 satisfied), one row after the other, each in its own
systemd unit with `HDF5_USE_FILE_LOCKING=FALSE`. Numbers are **TOTAL(run)** from
`[run].benchmark = true` — run only, no Julia start-up (start-up is the difference to the
"unit wall" column of the raw summary, 60–90 s). Raw logs and summary:
`/data/jbrioude/easlee_results/bench_matrix_20260922/`; tables regenerated by
`scripts/bench_matrix_table.py`. `PREFETCH=1` sets the met prefetch loader; `-t 16,1`
gives it an interactive thread (rows `gpu_t16i*`).

### Table A — every row (run-only seconds; ETEX 1 M parcels, 2 days; `bench_matrix_20260923_merged` = CPU rows of 03:02–07:24 + GPU rows re-measured 08:42–09:01 on the Val-switch fix + the 32-thread pair re-measured 10:04–10:18 (Jerome))

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
  pair was re-measured as a clean pair on the quiet machine at Jerome's request: 374 s
  (reproducing the matrix row to 0.1 %) against 360 s with prefetch. On the GPU the CPU is idle while the
  kernel runs, so the read is hidden for free: **168 s, −33 %, the best row** — 3.4× the
  16-thread CPU and 1.76× the best CPU row (64 t + prefetch, 295 s). The remaining GPU run
  is 71 % advance (kernel + host↔device copies) and 29 % serial floor (met read 17 s,
  release 16 s, output kernel 10 s, untimed 5 s).
- **The GPU regression, found by this matrix and fixed the same morning.** The first GPU
  pass read 167 s in the advection group where the same config had given 128–137 s on
  21 September (under load) and 294 s whole-run in July. A three-point bisect on the quiet
  machine (`/tmp/bisect_gpu.log`): `4345524` (before the five-switch port) **116.3 s**,
  `2bdf8f8` (the port) **167.5 s**, head **167.0 s** — the whole +44 % is the port, the KA
  fixes of `8a12d39` cost nothing. Cause: the switches (`use_hanna`, `mesoscale`, `kh_ft`,
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

## 0c. Cost-centre speed-ups with convection ON — zion, 2026-10-01 (the paper's Table 5)

Jerome (1 Oct 2026): *"it would be more informative if the speed-up factors for 16 or 32 CPU
threads and GPU are shown for each category instead of the time spent in seconds."*  §0b gives
every cost centre except convection, which the published ETEX protocol switches off.  The ETEX
ERA5 carries `avg_umf`/`avg_dmf`/`cc`, so the same protocol was rerun with
`[physics.convection] scheme = "massflux"` (`configs/bench_{cpu,gpu}_etex_1M_conv.toml`,
switch verified at RunSpec level before launch) at 1 / 16 / 32 threads and on the GPU —
sequential, one systemd unit per row, quiet machine, governor `schedutil`, commit `70a426f`.
Raw logs and summary: `/data/jbrioude/easlee_results/bench_conv_20261001/`; the table is
regenerated by `scripts/bench_costcentre_table.py --conv-only`.

| cost centre | 1 t [s] | share @ 16 t | 16 t | 32 t | GPU (4 t) |
|---|---:|---:|---:|---:|---:|
| meteorological read + preparation | 86.9 | 17 % | 0.77× | 0.74× | 0.80× |
| advection + turbulence | 5937.2 | 69 % | **12.94×** | **22.65×** | **41.60×** |
| convection (mass-flux redistribution) | 482.1 | 11 % | 6.73× | 8.88× | 3.06× |
| output gridding (deposit + uniform kernel) | 9.9 | 1 % | 1.01× | 1.00× | 1.03× |
| release, bookkeeping, other | 7.4 | 1 % | 0.84× | 0.76× | 0.88× |
| **whole run** | **6527.3** | 100 % | **9.81×** | **14.30×** | **15.14×** |

Absolute totals: 6527 s (1 t), 665 s (16 t), 456 s (32 t), 431 s (GPU).  Against the
convection-off matrix of §0b the same rows are 5910 / 576 / 374 / 249 s, i.e. convection adds
8–73 % of wall clock depending on where the advance has gone.

### Reading

- **Convection parallelises, but only half as well as the advance**: 6.7× at 16 threads and
  8.9× at 32 against 12.9× and 22.7× for the advance.  It is `Threads.@threads :static` over
  the parcels of each nest like the advance, but each parcel walks a *column* (profile
  extraction, flux build, the conveyor and the subsidence walk), so it is memory-bound on the
  3-D met arrays rather than compute-bound, and it keeps a per-thread buffer set.
- **Convection is the bottleneck of a GPU run.**  The mass-flux redistribution has no device
  path — it runs on the CPU with whatever threads the GPU row is given (4 here) — so its GPU
  column is a 4-thread CPU figure, 3.06×, and it becomes **37 % of the GPU run** (157.5 s of
  431 s) against 11 % at 16 threads.  The whole-run GPU speed-up therefore falls from 23.7×
  (convection off) to 15.1×.  Porting convection, or simply giving the GPU row more CPU
  threads (`-t 16` instead of `-t 4`), is the next lever for full-physics GPU runs.
- **Everything else is serial and flat**, as in §0b: met read 0.74–0.80× (it *slows* with
  threads, 87 → 117 s, GC on 14 GiB of allocation), output kernel 1.00×, release 0.8×.  With
  the asynchronous prefetch the met read is hidden (§0b: 5.6×, and 35.2× whole-run on the GPU
  with convection off); the convection matrix was run without prefetch so that the four rows
  differ by thread count alone.
- **The practical reading for a user**: the whole run gains ~10× at 16 threads and ~14× at 32;
  beyond that the serial floor (met read + output + bookkeeping ≈ 130 s) and the poorly
  scaling convection dominate, which is why the 64-thread and GPU rows of §0b add little once
  full physics is on.

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
- **Headroom — REVISED 2026-09-21 after measuring, and the first item is already done.**
  *Met upload only on refresh* was implemented on 2026-09-04 (review A5,
  `advance_ka_gpu.jl`): device met buffers are resident, a bracket change forces a
  full re-upload, and a ~1024-element strided fingerprint per array catches in-place
  edits so only changed arrays are copied. The other two items were sized with
  `[run] benchmark = true` on this very config:

  | group | `-t 4` | `-t 16,1` | share |
  |---|---:|---:|---:|
  | **metread** | **225.0 s** | **176.1 s** | **53–57 %** |
  | advection (GPU kernel + copies) | 137.0 | 128.5 | 34–38 % |
  | release | 14.5 | 12.1 | 4 % |
  | **output (deposit)** | **13.2** | **12.2** | **3.3 %** |
  | TOTAL | 397.9 | 335.8 | |

  So **the GPU deposit is worth at most 3.3 %** — and only for `kernel = "uniform"`;
  the adaptive kernel (the super-linear cost centre, 272.9 s at 20 M parcels in
  `BENCHMARKS_COSTCENTRES.md`) has no KA implementation at all. **Reading the met is
  the run.** Two findings on it: `prefetch.jl` spawns its loader on the `:interactive`
  pool, which `-t 4` does not provide, and giving it one (`-t 16,1`) cut metread by
  49 s and the whole run by 16 %; and the met read allocates **14 GiB** over 193 calls
  (~74 MB per call), i.e. it decodes into fresh arrays rather than preallocated
  buffers. Those are the two levers worth pulling, not the deposit.

  ⚠ Both rows above were measured while three 8-day chemistry runs held 48 of the
  machine's 64 threads, so the ABSOLUTE times are inflated against the 294 s in the
  table above (which was measured on a quiet machine). The comparison between the two
  columns is indicative, not clean — repeat it on a quiet machine before quoting.

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
- **`ka_*` ≈ `threads` but not identical.** The KA kernel is a *reimplementation*
  of the Thomson turbulence advance; it disperses slightly more than the
  production threads path (lower FB / RANK by ~0.1). Differences come from the
  turbulence sub-stepping detail, not RNG alone. For science, the validated
  threads path remains the reference; `ka_*` is for GPU acceleration studies.

## 3. Known limitations of the KA / GPU advance

- ⚠ **`[physics] surface_layer_wind` DEFAULTS TO `true`, and it forces the threads
  advance.** Any config that does not set it to `false` silently runs on the CPU
  whatever `EASLEE_ADVANCE` says — the only evidence is a `Warning: ... overriding
  EASLEE_ADVANCE=... with threads` line in the log. **Check the log for
  `mode = :ka_gpu` before believing a GPU timing.** (This invalidated a first
  attempt at the Hanna validation on 2026-09-21: both "backends" were the CPU.)
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
  (`scripts/ka_compare_grids.py`; outputs `bench_all5_{threads,ka_gpu}`, the 22 Sept
  ones archived as `*_PREFIX_20260922`).
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
CAPTEX_OBS=<ETEX dir> CAPTEX_MET=/data/ETEX EASLEE_ADVANCE=ka_gpu MO=false \
  TURB=thomson CONV=off IFINE=5 MODE=run DIRECTION=forward RELEASE=1 \
  NP=1000000 RUNDAYS=2 DXOUT=0.25 julia -t 4 --project=easlee validation/validate_captex.jl
```

(`MO=false` was added to validate_captex_run.jl for exactly the reason above; it
is not in the July command line because the switch did not exist yet.)

## Full-physics per-process benchmark (2026-07-09, domain-fill chemistry)

Domain-fill CB05×CAMS chemistry on a cropped 20°×15° Pacific box
(`[meteorology].crop="auto"` → 85×65×138 of the global 0.25° ERA5 grid),
250k parcels, 4×5-min sync steps spanning one ERA5 file rotation, chemistry
every step, ALL production physics on: hanna turbulence, mass-flux convection
+ wet dep, pairwise chem mixing, LUT photolysis + cloud τ, CAMS aerosols/het/
AOD.  20-core box, RTX A1000.  Per-step = group total / 5 calls; chemistry has
~±7 % run-to-run solver variance (adaptive steps).

| process (bench group)      | CPU 1 thr | CPU 6 thr | scaling | share (6 thr) |
|----------------------------|-----------|-----------|---------|---------------|
| chemistry (CB05 solver + J)| 59.4 s/st | 27.3 s/st | 2.2×    | 78 %          |
| met read + field calc      | 13.8      | 5.8       | 2.4×    | 17 %          |
| advection (incl. turb)     | 3.9       | 0.85      | 4.6×    | 2.4 %         |
| convection (mass-flux)     | 1.77      | 0.38      | 4.7×    | 1.1 %         |
| CAMS reseed                | 0.29      | 0.29      | 1×      | 0.8 %         |
| chem extract (T/p/cld-τ)   | 0.26      | 0.11      | 2.3×    | 0.3 %         |
| pairwise mixing            | 0.068     | 0.068     | 1×      | 0.2 %         |
| wet deposition             | 0.095     | 0.064     | 1.5×    | 0.2 %         |
| output (streamed frame)    | 0.16/frame| 0.17/frame| —       | 0.5 %         |

Derived (A/B toggles — these live inside hot kernels, not separable by timer):
- **turbulence** = Δadvection(hanna − none) = **0.77 s/step** — 90 % of the
  advection cost, 2.2 % of the run.
- **LUT photolysis** = Δchemistry(lut − analytic) ≈ **0 within the ±7 % solver
  variance** (the analytic run came out *slower*); the TUV table lookup is
  effectively free.  Cloud-τ extraction: 0.08 s/step.

GPU (`EASLEE_ADVANCE=ka_gpu`, advance only — CB05 chemistry is CPU by design):
- at 250k parcels GPU advection is **2.5× slower** than 6 CPU threads
  (10.7 vs 4.25 s per 5 calls: per-step met+state upload dominates at this
  size); the §1 result (~1.8× GPU win) needs ≥ ~1M parcels.  Everything else
  unchanged → hybrid GPU-advance/CPU-chemistry only pays at 1M+.
- GPU *chemistry*: three findings.  (i) Install the optional stack as
  `Pkg.add([... "DiffEqGPU@3.15" ...])` — the LATEST DiffEqGPU does not
  precompile against the SciMLBase 3.x this model pins (the run then warns and
  falls back to CPU).  (ii) The ENGINE deliberately forces CPU chemistry
  (`gpu_thresh = typemax(Int)`, build_run.jl): the engine feeds per-parcel
  H2O/lat/lon into the solver — scientifically required — and the GPU kernels
  only support scalar H2O, so GPU routing would silently degrade the
  chemistry.  (iii) Even routed directly (micro-benchmark, identical 250k
  input, bit-identical output), bb12-GPU is 1.28× SLOWER than 6 CPU threads on
  this card at 100 % GPU utilisation — the RTX A1000 executes Float64 at 1/32
  of FP32, capping any Float64 stiff solve at ~CPU parity.  GPU chemistry
  needs an FP64-strong card (A100-class) or mixed-precision kernels, plus
  per-parcel-H2O support, to be worth switching on.  (bb12 CPU: 2.7 s/step vs
  CB05's 27.3 — mechanism size is worth 10× before any GPU.)
- Overlapping GPU advance with CPU chemistry is not worth building: advance is
  3 % of the step — the overlap ceiling.  The lever with actual headroom is the
  chemistry solver's allocation churn (30 GiB/step at 250k parcels — the GC
  serialisation behind the 2.2× thread scaling); cutting it would speed EVERY
  configuration.

Peak RSS: 6.5–7.6 GB (cropped; +1 GB with the CUDA context loaded).

## CB05 Jacobian sparsity scan (2026-07-09) — sparse-GPU-kernel feasibility

FD Jacobian of `rhs_cb05!` at six states (BL/FT/UT × day/night + polluted BL),
composition from `cb05_background()` ±50 %:

| pattern | nnz | of 49² | UMFPACK LU fill | per-thread factors |
|---|---|---|---|---|
| dense (reference) | 2401 | 100 % | — | 19.2 KB |
| structural union | 404 | 16.8 % | 555 | **4.3 KB** |
| numerically-active union | 277 | 11.5 % | 371 | 2.9 KB |

- Pattern is REGIME-STABLE (380–403 structural nnz per state) → one compiled
  KPP-style unrolled sparse-LU pattern serves all regimes, no divergence.
- Block structure: one 41-species strongly-connected core + 8 feed-forward
  singletons (substitution-eliminable).  No small-block decomposition —
  everything couples through OH/NOx; the win is sparsity (4.4×), not blocks.
- Feasibility: 4.3 KB factors + 1.2 KB state per parcel; batched SoA layout →
  ~5.5 GB device memory at 1M parcels (fits a 16 GB card).  In Float64 a
  sparse CB05 kernel on a 1/64-FP64 consumer card is a modest win; mixed
  precision is the transformative step.

Per-parcel H2O/lat/lon in GPU chemistry (same date): the bb12 GPU path now
consumes the per-parcel vectors (they enter the parameter build; kernel
unchanged) — verified vs CPU to solver tolerance (max rel 4e-3, mean 5e-4)
with a 41× max-rel effect vs the old scalar H2O.  The engine gate is now
mechanism-dependent: bb12 routes to GPU under EASLEE_GPU=1; all other
mechanisms stay CPU-forced until their kernels support the per-parcel inputs.

## RB23 thread-scaling re-measure (2026-07-10)

*(Mode caveat: measured under the POWERSAVE governor — post-reboot default,
noticed later that day.  The 2026-07-09 tables above are performance-mode.
Thread-scaling and allocation numbers are mode-insensitive; the absolute
s/call are deflated and must not be compared across the sections.)*

Same domain-fill CB05×CAMS Pacific setup (crop auto, 250k parcels, 4×5-min
sync steps, all production physics), but the config's own 15-min chemistry
step (2 chemistry calls; per-CALL numbers below solve dt = 900 s — do not
compare per-step against the 2026-07-09 FBDF table, which solved dt = 300 s
every step) and the now-default `solver="rb23"`:

| chemistry (rb23) | 1 thr | 6 thr | scaling |
|---|---|---|---|
| time / call | 52.5 s | 19.7 s | **2.7×** |
| alloc / call | 2.13 GiB | 2.13 GiB | — |

- Allocation churn is down ~14× vs the FBDF path (2.1 vs ~30 GiB per call) on
  top of the 3–5× absolute speedup, but 6-thread scaling only moved 2.2× →
  2.7×: **solver allocations no longer explain the gap**.  The remaining
  serialisation is the still-nonzero GC pressure (2 GiB/call is not free) +
  memory bandwidth of the per-parcel dense 49×49 LU workspaces.
- metread scales 3.4× (17.8 → 5.2 s/step); everything else is noise at this
  scale.
- Consequence for the roadmap: with chemistry now 19.7 s / 15-min call at
  250k on 6 threads, further CPU alloc work is LOW priority; the bigger
  remaining lever for chemistry-dominated runs is the CB05 sparse-FP32 GPU
  stage 2 (above), or simply the rb23 CPU path as-is.

## 24-h Pacific chemistry: FBDF-CPU vs rb23-fp32-GPU — sanity check (2026-07-11)

Same config (1M parcels, 24 h, domainfill CB05 × CAMS, chemistry every
15 min) as the production run of record (results/chem_pacific/run_24h_v1.log,
FBDF fp64 CPU, performance governor, uncapped): chem_pacific_gpu32.toml =
solver rb23 + precision fp32 + EASLEE_GPU=1.  NOT a clean benchmark — the
rerun carried two handicaps: 16 of 24 sim-hours under the powersave governor
(CPU 400 MHz floor, GPU 210/2100 MHz; switched mid-run) and a systemd
MemoryHigh=36G cap with a 47.4 GB peak → reclaim throttling.  Clean rerun
deferred (with the RTX 4080 columns).

- Wall: 30 382 s → **13 982 s (2.2×)** despite the handicaps.  Segment paces
  9.6 (powersave) vs 7.9 (performance) min per sim-hour → clean-performance
  extrapolation ≈ 3.2–3.5 h (~2.5×), likely better uncapped.
- Chemistry proper: **1737 s (12.4 % of the run, 97 calls — ~18 s per 1M-parcel
  call, mostly under the throttled GPU)**.  The old run predates the chemistry
  bench group; subtracting the shared non-chemistry cost estimates its
  chemistry at ~18 ks → **chemistry speedup ≈ 10×** (rb23 3–5× × GPU, as
  predicted by the component benchmarks).
- Science: species-by-species vs the old output — O3 corr 0.990 / ratio
  0.9998, CO 0.992 / 1.0002, CH4 0.995 / 1.0000; the inert CH4 correlation IS
  the parcel-position noise floor (unseeded turbulence RNG), so the reacting
  species are indistinguishable at the resolvable level.  Radical/NOx totals
  move a few % (OH −0.9 %, HO2 −4.2 %, NO +3.9 %) — the FBDF-legacy-defaults
  outlier + fp32 floor, far inside the mechanism's ±15–35 % CAMS envelope.
  O3 final-frame mean profile ratio 1.000 ± 0.001 at every level.
- FINDING for the next optimisation pass: BOTH runs spend ~10.5 ks
  (~39 s/step, now 75 % of the accelerated run) in per-step work outside every
  bench group (suspects: output/conc deposit at 1M, domainfill inject/cull,
  GC/reclaim; the advection group's 21 s total says transport is largely
  untracked too).  Add bench groups before the clean benchmark.

## The untracked 39 s/step — named (2026-07-11, instrumented rerun)

2-h slice of the Pacific GPU-fp32 run with new bench groups (:chemgrid =
chemconc_sample!, :dfill = domainfill inject+cull, :sync), performance
governor, MemoryHigh=46G (no reclaim), 1M parcels.  TOTAL 1398 s = init
~860 s (CAMS seeding + met scan + build) + 24 × ~21.5 s/step, FULLY
attributed:

| per step | s | share of step |
|---|---|---|
| **chemgrid (chemconc_sample!)** | **14.2** | **66 %** |
| chemistry (fp32 GPU, ÷ steps) | 3.8 | 17 % |
| metread | 1.6 | 7 % |
| chem_extract | 1.2 | 6 % |
| cams_reseed | 0.5 | 2 % |
| dfill / advection / rest | 0.4 | 2 % |

The 24-h run's "39 s/step untracked" was this ×(powersave + memory-reclaim).
chemconc_sample! loops 49 species × 1M parcels EVERY transport step and
allocates 1.3 GiB/step.  Fixes, in order of value: (1) sample the chemistry
grid at the output.sample cadence like the transport grid (immediate ~3×);
(2) kill the per-call allocations; (3) thread it.  With chemgrid fixed the
projected 24-h wall is ~init + 24×(7 s/step)×12 ≈ 1.0–1.2 h — the original
1.5 h estimate was right once the REAL bottleneck is removed.
Note advection at 1M really is 0.08 s/step (simple memory-bound loop) — the
earlier suspicion was unfounded.
