# EASLEE time steps — naming convention

EASLEE distinguishes two kinds of "time step", which a single `COMMAND`-style
namelist tends to hold side by side:

- **Cadence / operator-split steps** — *when* to do something (advance, sample, write,
  react, convect). Every one is an integer multiple of a single master clock.
- **Accuracy sub-steps** — *how finely* to integrate *within* one transport step. These
  are subdivisions rather than cadences, which is why `ctl` and `ifine` live under
  `[physics]` here and not with the output cadences.

## Coming from FLEXPART

EASLEE is in the FLEXPART lineage and the two models solve the same problem, so
if you know FLEXPART's `COMMAND` this table is the fastest way in. It maps each
`COMMAND` entry to the TOML key that plays the same role. The names differ
because EASLEE groups the keys by what they control (`[run]`, `[output]`,
`[physics]`) rather than by file; the meanings are the same.

| Concept | FLEXPART `COMMAND` | EASLEE field | Default | Constraint |
|---|---|---|---|---|
| Master sync clock | `LSYNCTIME` | `run.timestep` | `15min` | base unit |
| Grid write interval | `LOUTSTEP` | `output.interval` | `1h` | × timestep |
| Averaging window | `LOUTAVER` | `output.average` | = interval | × timestep, ≤ interval, multiple of sample |
| Sampling cadence | `LOUTSAMPLE` | `output.sample` | = timestep | × timestep, ≤ average |
| Turbulence accuracy | `CTL` | `physics.ctl` | `5` | adaptive vs Lagrangian timescale T_L (`<0` = off) |
| Vertical refine | `IFINE` | `physics.ifine` | `4` | integer ≥ 1 |
| Convection split | — | `physics.convection.timestep` | = run.timestep | × timestep |

The last row has no FLEXPART counterpart: convection there runs on the
transport step, whereas EASLEE lets it have its own cadence.

## Worked example — the output cadence

*"Output every hour, snapshot the trajectories every 15 min, average over the hour"* →
```toml
[run]
timestep = "15min"
[output]
interval = "1h"      # write a frame each hour   (LOUTSTEP)
average  = "1h"      # average over the whole hour (LOUTAVER)
sample   = "15min"   # bin parcels onto the grid every 15 min (LOUTSAMPLE)
```
→ 4 snapshots accumulated and averaged into one written frame. **sample ≤ average ≤ interval**,
each a multiple of `run.timestep`, and `average` a multiple of `sample`.

## Accuracy sub-steps (CTL / IFINE)

These are *not* cadences — they refine the integration inside one `timestep`:
- `physics.ctl`: the Langevin turbulence step is reduced to `T_L / ctl` (adaptive; larger
  ctl ⇒ finer steps in turbulent regions). Governs horizontal **and** vertical turbulence.
- `physics.ifine`: integer extra subdivision of the **vertical** step (`dt_vert = dt_turb / ifine`).

## Process steps

- `physics.convection.timestep` — convective redistribution cadence. Default = `run.timestep`.

## Auto vs user-set
- **User sets:** `run.timestep`, `output.interval`.
- **Auto-derived (overridable):** `output.average = interval`, `output.sample = timestep`,
  `ctl = 5`, `ifine = 4`,
  `physics.convection.timestep = run.timestep`.
- **Validated:** every cadence is an integer multiple of `run.timestep`;
  `sample | average | interval`; `ifine ≥ 1`. (See `validate_runspec`.)

## Status
Config + validation implemented: `run.timestep`, `output.{interval,average,sample}`,
`physics.{ctl,ifine}`, `physics.convection.timestep`. `convection.timestep` is parsed/validated and
defaults to current behaviour; a convection sub-cadence inside `step_model!` is a follow-up.
