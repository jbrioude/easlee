# Paper tables from the 2026-09-23 benchmark matrix (for Sect. 5 of EASLEE v1.0 GMD)

Source: `docs/BENCHMARKS.md` §0b (workstation, Threadripper 3970X 32c/64t, RTX 4080; ETEX-1, 1 M parcels, 2 days, ERA5 hourly, Thomson, convection off, otherwise idle machine, CPU fully boosted). Times are RUN-ONLY (the model's own `TOTAL(run)`); Julia start-up adds ~1 min to every row and is excluded so that the speed-ups are the model's. The July Table 4 was whole-run including start-up: 111.4 / 32.3 / 12.2 / 7.8 / 6.4 / 4.9 min; the new run-only rows are 98.5 / 28.0 / 9.6 / 6.2 / 5.3 / 4.2 min, i.e. the same machine and code family within the start-up offset except the GPU rows, which are re-measured on the kernel fix of 23 Sept.

## Table 4 (replacement) — run time and speed-up, reference ETEX-1 configuration

| Backend | Threads | Prefetch | Run time | Speed-up vs 1 thread | Speed-up vs 16 threads |
|---|---:|---|---:|---:|---:|
| Workstation CPU | 1 | off | 98.5 min | 1.0× | 0.10× |
| Workstation CPU | 2 | off | 50.9 min | 1.9× | 0.19× |
| Workstation CPU | 4 | off | 28.1 min | 3.5× | 0.34× |
| Workstation CPU | 8 | off | 16.2 min | 6.1× | 0.59× |
| Workstation CPU | 16 | off | 9.6 min | 10.3× | 1.00× |
| Workstation CPU | 32 | off | 6.2 min | 15.8× | 1.54× |
| Workstation CPU (SMT) | 64 | off | 5.3 min | 18.7× | 1.82× |
| Workstation CPU | 16 | on | 9.3 min | 10.6× | 1.03× |
| Workstation CPU | 32 | on | 6.0 min | 16.4× | 1.60× |
| Workstation CPU (SMT) | 64 | on | 4.9 min | 20.0× | 1.95× |
| Workstation GPU (RTX 4080) | 4 | off | 4.2 min | 23.7× | 2.31× |
| Workstation GPU (RTX 4080) | 4 | on | 2.8 min | 35.2× | 3.43× |

Laptop rows (20-core, RTX A1000) unchanged from July: 27.5 min at 16 threads, 15.4 min GPU (1.78×), 432 MiB.

## Table 5 (new) — per-process time and speed-up on the same runs

| Process | 1 thread | 16 threads | 64 threads | GPU (4 CPU threads) | GPU + prefetch | CPU scaling exponent |
|---|---:|---:|---:|---:|---:|---:|
| Met fields read + preparation | 94 s | 108 s (0.9×) | 117 s (0.8×) | 109 s (0.9×) | 17 s (5.6×) | N^+0.06 |
| Advection + turbulence | 5795 s | 444 s (13.0×) | 175 s (33.1×) | 118 s (49.2×) | 120 s (48.3×) | N^−0.86 |
| Output kernel (deposit to grid) | 10 s | 10 s (1.0×) | 10 s (1.0×) | 10 s (1.0×) | 10 s (0.9×) | N^+0.01 |
| Release | 7 s | 9 s (0.8×) | 10 s (0.7×) | 8 s (0.8×) | 16 s (0.4×) | N^+0.09 |
| Frame dump, bookkeeping | 4 s | 4 s (1.0×) | 4 s (1.0×) | 4 s (0.9×) | 4 s (0.9×) | N^+0.01 |
| Whole run | 5910 s | 576 s (10.3×) | 316 s (18.7×) | 249 s (23.7×) | 168 s (35.2×) | N^−0.72 |

Convection and deposition are off in the ETEX protocol; the July full-physics profiling (Table 5 of the draft: 500 k parcels, 6 h, Hanna, mass-flux, 6 threads) remains the source for those cost centres (convection 31 % of the loop, 4.7× at 6 threads).

## Text points for Sect. 5

- The advance scales as N^−0.86 from 1 to 64 threads (13× at 16, 33× at 64; parallel efficiency 49 % at 32 physical cores, +1.4× from SMT). The whole run scales as N^−0.72 because the met read (94–117 s), the output kernel (10 s), the release (7–10 s) and the frame dump (4 s) are serial: a ~125 s floor, 2 % of the run at 1 thread and 40 % at 64.
- The GPU advance (118 s) is 3.8× the 16-thread advance and 1.5× the 64-thread one: one consumer GPU beats the 32-core CPU on the parallel part; the whole run is 2.3× vs 16 threads and 1.27× vs 64.
- Prefetch hides the met read only where compute overlaps it. On the CPU the single-threaded decode competes with the workers (a shared worker thread delays the advance barrier; its 14 GiB of allocation pause every thread at garbage collection), so the read drops 108 → 14 s but the advance rises by 68–78 s: net −3 % (16 t), −4 % (32 t), −7 % (64 t). On the GPU the CPU is idle during the kernel and the full read is hidden: 249 → 168 s, −33 %, the best row at 3.4× the 16-thread CPU. The loader reads one file ahead in one task and the decode is single-threaded, so a read-dominated workload (few parcels, many files) gains from parallel reading, not from prefetch as it stands.
- Method note: physics switches must reach a KernelAbstractions kernel as types (`Val{Bool}`), not as runtime Bools — the GPU compiler specialises on types and otherwise compiles every branch into one kernel (+44 % on the advance with all switches off, found by this matrix and fixed).
