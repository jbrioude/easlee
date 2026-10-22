# ERA5 etadot units — the bug, and its independent confirmation

## The bug (fixed 2026-08-12)

`easlee/src/io/initialize.jl` converted the archived vertical velocity to a
level-index rate by dividing by the pressure thickness:

    lev /= (puv[k+1] - puv[k])          # _levdrho_chunk!

That is flex_extract's convention, where `calc_etadot` PRODUCES η̇·∂p/∂η in
Pa/s. The raw ERA5 archive is not that field: **param 77 is η̇ in s⁻¹**. The
division is then dimensionally wrong — s⁻¹ ÷ Pa is s⁻¹Pa⁻¹, not a level rate —
and leaves the resolved vertical advection too small by ∂p/∂η ≈ 10⁵.

Affected: every ERA5 model-level (hybrid GRIB) run — ETEX, CAPTEX, Réunion, HT
transport. HT runs in hybrid/isentropic mode are least affected above the blend
region, where the level index is re-derived from `theta_p`; below σ_r the
kinematic levdot dominates and was ~0. Not affected: the MERRA-2/GFS paths
(ω [Pa/s] → correct conversion) or AROME/WRF (w in m/s).

## Three independent confirmations

1. **GRIB header.** `grib_ls -p units` on param 77 returns `s**-1`.
2. **Dimensional.** Only a Pa/s input makes `field ÷ Δp` a level-index rate.
3. **ERA5's own ω (param 135, Pa/s), retrieved for this check.** In hybrid
   coordinates ω = B·∂p_s/∂t + η̇·∂p/∂η, so regressing the archived param 77
   against ω must give a slope of ∂p/∂η ≈ 10⁵ if it is η̇ [s⁻¹], or 1 if it is
   already Pa/s. Least squares through the origin, 19 481 points per level,
   1994-10-23 12:00 UTC:

   | GRIB level | slope ω/η̇ | r | mean abs η̇ | mean abs ω |
   |---|---|---|---|---|
   | 60 | 1.012e5 | 1.000 | 1.25e-7 s⁻¹ | 1.27e-2 Pa/s |
   | 100 | 4.33e4 | 0.704 | 2.27e-6 | 1.44e-1 |
   | 130 | 6.82e4 | 0.278 | 6.63e-7 | 1.38e-1 |

   Level 60 (upper troposphere / lower stratosphere) is the clean test: slope
   1.01e5, correlation 1.000. Lower levels degrade as expected, where ω also
   carries B·∂p_s/∂t and terrain effects.

## Both conventions are supported (2026-08-12)

The original `÷ Δp` was NOT simply wrong: it is correct for **flex_extract**
output, which is what `calc_etadot` produces and what that code path was written
for. Making the Δη conversion unconditional would have broken flex_extract
input by the same 1e5 in the other direction. The reader now handles both:

- **`#info etadot_units`** in the variables namelist: `auto` (default) | `s-1`
  (raw MARS/ERA5 dissemination) | `pa/s` (flex_extract).
- **`auto`** decides ONCE per met source, on the first file (a file set is either flex_extract output or raw dissemination, never mixed), from the field magnitude. eta-dot is ~1e-6 s^-1
  in the mid-troposphere; eta-dot*dp/deta is ~1e-1 Pa/s — five orders apart, so
  a threshold at 1e-2 sits two orders clear of both populations. The detector
  uses the MAX over interior levels (not the mean, so a quiescent domain cannot
  flip it) and logs its decision with the statistic at the point of use.
- **Known limit, by design:** an all-zero or vanishingly small field is read as
  s^-1. That is why the explicit override exists; declare it in the namelist for
  any dataset where the magnitude test could be ambiguous.

Pinned by `test/test_etadot_units.jl` (10 checks), including that both
conventions describing the SAME physical motion produce the same level-index
rate.

## A second, independent check on the RESULT (2026-09-19)

The detector above judges the magnitude of the ARCHIVED field.  A second check,
`_levdot_sanity!`, judges the RESULT: once the divisor has been applied, the
level-index rate times the local layer thickness is a physical vertical
velocity, and a convention that was picked wrongly — or forced wrongly through
`#info etadot_units`, or an archive that is in neither convention — shows up as
a `w` orders of magnitude away from the synoptic range whatever the input
looked like.  That is precisely the case the magnitude test cannot cover on its
own, since it reads an all-zero or vanishingly small field as s⁻¹ by
construction.

It runs **once per met source**, on the same file as the detector, and costs one
pass over the middle third of the levels.  The plausible band is deliberately
wide — max|w| between **1e-3 and 50 m/s** over the interior levels — because it
exists to catch a factor ~1e5, not to police the meteorology.  Outside the band
it emits a warning naming the convention that was applied and the namelist key
that overrides it.

Measured on the ETEX-1 ERA5 files:

| run | max\|w\| | mean\|w\| | outcome |
|---|---:|---:|---|
| `etadot_units = auto` → s⁻¹ (correct) | 0.457 m/s | 1.52 cm/s | info, passed |
| `etadot_units = pa/s` forced (wrong) | 5.06e-6 m/s | 1.49e-7 m/s | **warning** |

The two differ by 9.0e4, i.e. the dp/deta factor, which is what a flex_extract
archive read as raw dissemination (or the reverse) costs.  Pinned by
`test/test_etadot_units.jl`.

## After the fix

Mean |w| = 1.9 cm/s at 880 hPa, 1.8 cm/s at 246 hPa, peaks 0.36 m/s — the
synoptic range; previously ~10⁵ smaller, i.e. no resolved vertical transport.

Consequences measured so far: spurious parcel clustering falls sharply — the
spread of the log parcel-density, sd(ln c), went from GROWING 0.145 → 0.290 to
CONTRACTING 0.145 → 0.083 over 48 h — and ETEX gains in structure while losing
in bias (docs/ETEX_RESULTS.md: R +0.006, KSP −1.0, FB +0.053, RANK flat).

Data: `/data/ETEX_VOD/omega/IFS_EA_W.grb` (param 135, 2 dates × 2 bases × 2
steps × 137 levels), MARS request `ea_w.req` on the ECMWF side.

## Levels (2026-10-22): the two conventions are not on the same levels either

The same continuity argument fixes the flux through every interface from u, v and
p_s alone: F(k+½) = B(k+½)·ΣD − Σ_{j≤k} D_j, D_j = ∇·(v_j Δp_j).  Correlating the
archived value at GRIB level k with F above (k−½), below (k+½) and with their mean
(`validation/check_etadot_levels.py`):

| archive | best match | slope at levels 136 … 70 | rms(137)/rms(136) |
|---|---|---|---|
| flex_extract η̇·∂p/∂η (EAg 2010-11-04 22 UTC, ERA5 global 2022-01-15 03 UTC) | interface k+½ below the label: **half levels**, level 137 = the ground | 1.08-1.19 | 0.053, 0.058 |
| raw MARS η̇ (ETEX IFS_EA 1994-10-24 06 UTC +10 h) | mean of k±½: **full levels** | 1.06-1.14 | 0.319 |

Loaded bottom-up and unshifted, slot k is read at zlev[k], the half level below
full level k: right for flex_extract, half a layer low for raw MARS, whose values
are averaged to the half levels.  `_etadot_levels!` decides once per met source
from the lowest-level ratio (threshold 0.15); `#info etadot_levels = full | half`
skips the test.  Until 2026-10-22 raw MARS archives (ETEX, CAPTEX) were read half
a layer low; for a few hours on 2026-10-22 (b778cec..02858a1, never released)
flex_extract archives were averaged, i.e. read half a layer high.
