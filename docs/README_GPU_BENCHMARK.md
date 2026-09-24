# EASLEE GPU-vs-CPU speed benchmark — run it yourself

This is the exact protocol behind the GPU rows of `docs/BENCHMARKS.md`, packaged
so it needs nothing but this directory and an EASLEE checkout.

## What is in here

```
bench_gpu_etex_1M.toml      GPU arm   (4 CPU threads feeding the device)
bench_cpu_etex_1M.toml      CPU arm   (16 threads; identical in every other line)
met/19941023 19941024 19941025    ERA5 met, ETEX-1 domain (−10…30 °E, 40…70 °N)
```

4.2 GB of met. It is plain ERA5 GRIB as retrieved from MARS (`IFS_EA_*` model
levels + `IFS_SFC_*` surface companions, 12-hourly files, 13 chunks each).

## Setup

1. Copy the two TOMLs into `configs/` of your EASLEE tree.
2. Edit `[[meteorology]].path` in both to wherever you put `met/`.
3. The GPU arm needs CUDA working in your Julia environment. EASLEE loads CUDA
   lazily, so the CPU arm works without it.

## Run

```bash
EASLEE_ADVANCE=ka_gpu  julia --project=easlee scripts/easlee.jl configs/bench_gpu_etex_1M.toml
EASLEE_ADVANCE=threads julia --project=easlee scripts/easlee.jl configs/bench_cpu_etex_1M.toml
```

`scripts/easlee.jl` reads `[run].threads` from the TOML and passes it to the
child julia as `-t`, so the thread count is in the config, not on the command
line. Change `threads` in `bench_cpu_etex_1M.toml` to whatever core count you
want to quote.

**Run each arm twice and report the SECOND run.** The first GPU run pays the
CUDA kernel compilation (about 3 minutes on the RTX 4080 here).

The number to compare is the **whole-run wall clock**, and `[run].benchmark =
true` also prints a per-group breakdown at the end:

```
── benchmark (per-group totals) ──────────────────────
  group             time[s]        %        alloc     calls
  metread            ...
  advection          ...          <- the only group the backend changes
  output             ...
  TOTAL(run)         ...
```

`advection` is the GPU-vs-CPU comparison proper; `TOTAL(run)` is the number
quoted in the paper, and it is much less impressive than the advection ratio
because GRIB met reads, the deposit, and NetCDF output are the same CPU work in
both arms (Amdahl).

## The trap to check for

Three settings in the TOML exist only to keep the two arms on the **same
physics**. The KernelAbstractions/GPU kernel does not implement them, and EASLEE
falls back to the CPU advance if any is on — you would then benchmark threads
against threads and measure a 1.0× "speedup":

| setting in the TOML | value required for the GPU arm |
|---|---|
| `physics.surface_layer_wind` | `false` (schema **default is `true`**) |
| `physics.turbulence.scheme` | `"thomson"` (Hanna on-the-fly is threads-only) |
| `output.kernel` | `"uniform"` (`adaptive`/`covariance` refuse `ka_gpu`) |

Also threads-only: `physics.turbulence.mesoscale`, `kh_ft > 0`,
`horizontal_tke`, and 2-D/1-D wind masks.

So check the log for both of these:

```
┌ Info: advance backend selected
└   mode = :ka_gpu
```

and the **absence** of

```
┌ Warning: ... is implemented by the threads advance only — overriding EASLEE_ADVANCE=ka_gpu with threads
```

The `mode = :ka_gpu` line is printed before the override happens, so it alone
does not prove the GPU ran. The warning is the decisive line.

## Sanity check on the parcel count

Do not shrink the run to make it quick: at small N the GPU **loses**. A 3-hour,
20 000-parcel smoke test on the RTX 4080 spends 12.9 s in advection against 2.5 s
on 8 CPU threads — kernel compile plus the per-step host↔device copy dominate.
The 1 000 000-parcel, 2-day configuration here is the smallest size at which the
comparison means anything.

## Reference numbers (workstation: Threadripper 3970X 32c/64t, RTX 4080 16 GB)

1 000 000 parcels, 2 days, whole-run wall clock including Julia startup:

| backend | threads | wall clock |
|---|---:|---:|
| `threads` | 16 | 12:10 (730 s) |
| `threads` | 64 | 6:22 (382 s) |
| `ka_gpu` | 4 | **4:54 (294 s)** |

RTX 4080 = 2.48× the 16-thread CPU, 1.30× the best CPU row. On a laptop RTX
A1000 6 GB the same run took 15:24 against 27:30 on 16 CPU threads (1.78×).

Two caveats on these reference numbers. They were measured on 2026-07-17, before
`surface_layer_wind` existed (2026-08-13), so both arms then ran without it —
which is what the TOMLs here reproduce. And the GPU path currently uploads the
met and the full parcel state every step (correctness first), so the whole-run
speedup is bounded well below the advance-kernel speedup.
