# EASLEE — configuring and running

A run is one TOML file.  No FLEXPART `COMMAND`/`RELEASES`/`OUTGRID`/`pathnames`
and no hand-written `AVAILABLE` file are needed.

## Quick start

```bash
# copy a config that is close to what you want, edit it, then run it:
julia --project=easlee scripts/easlee.jl my_run.toml
```

The config's top-level `task` field (`"transport"` — default — or
`"validation"`) selects what the run does; `[run].threads` sets the worker count
(the launcher passes `-t` for you).  To check a config without running it:

```julia
julia --project=easlee -e 'include("easlee/src/config/EasleeConfig.jl");
                           using .EasleeConfig; read_runspec("my_run.toml")'
# read_runspec validates every field and throws with a list of all problems
```

### Unknown keys (`strict_keys`)

Every key of the TOML is checked against the schema
(`easlee/src/config/schema.jl`): a key the model does not know — a typo
(`[physics.mixing] Kh = 500`), a key from another block, a retired name — is
an **error**, listing each unknown key with the nearest schema key:

```
unknown configuration key in my_run.toml — not in the schema (...), the model would silently ignore it:
  - physics.mixing.Kh  (did you mean "physics.mixing.K_h"?)
```

`strict_keys = false` at the top level downgrades this to a warning (the
config then loads with the unknown keys ignored).  Keys the reader accepts
in TOML-only spellings (`[species.settling]` as a table, `[release.chem_init]`
species names, `[chemistry.init.sources.<name>]`, the `origin/spacing/size`
grid vectors, …) are listed in `schema_extra_paths()` in the same file.
Test: `test/test_config_keys.jl`.

GUI: `julia --threads=auto --project=webgui webgui/easlee_app.jl` → the **Config**
tab is an auto-generated form (built from the model schema); edit, **Save**, **RUN**.

Ready examples are in [configs/](configs/).

## The config file

Grouped by concern, every field has a default; durations/lengths/masses accept
human units (`15min`, `2.0um`, `1kg`).

```toml
name = "my_run"

[run]
direction = "forward"          # forward | backward (source-receptor)
start     = 2025-05-31T07:00:00
stop      = 2025-05-31T10:00:00
timestep  = "15min"
benchmark = false              # per-group @time report at the end

[[meteorology]]                # one table per nest, finest first
model     = "ERA5"             # ERA5 | MERRA2 | GFS | AROME_HDF | MESONH | WRF | AROME
path      = "/data/ERA5/20250531/"
variables = "ERA5_CONV"        # → listvariables/ERA5_CONV.vars
available = "auto"             # scan the directory; or an AVAILABLE file path
crop      = "auto"             # domain crop, see below: auto | off | [lon0, lon1, lat0, lat1]

[[species]]
id   = "dust"
kind = "aerosol"               # gas | aerosol
aerosol_type = "dust"          # dust|sulfate|bc|oc|seasalt|ash|water (κ + refractive index)
density = 2600.0
[species.size]
distribution = "lognormal"     # lognormal | bins
dg = "2.0um" ; sigma_g = 2.0 ; nbins = 8

[[release]]
species   = "dust"             # or a list: ["SO2","NOx"]
lon = [76.5, 78.0] ; lat = [-42.0, -38.5] ; z = ["200m", "1500m"]
time      = [2025-05-31T07:00:00, 2025-05-31T08:00:00]
particles = 20000
mass      = "1kg"

[physics]
turbulence = "hanna"           # none | hanna | thomson
[physics.convection]
scheme = "massflux"            # off | massflux | cape
[physics.vertical]
mode = "kinematic"             # kinematic | hybrid | isentropic
[physics.deposition]
wet = true

[output]
quantity    = "both"           # concentration | mixing_ratio | both
per_release = true ; per_species = true ; aod = true
group       = "species"        # output-column grouping: "species" (one column
                               # per species — grid memory independent of the
                               # release-box count) | "release"/"box" (one per
                               # box; needed for per-box backward sensitivities).
                               # Omit = auto: species forward, release backward.
                               # See docs/OUTPUT_GRID.md.
path        = "juliaout/output_dust/"
interval    = "1h" ; sample = "15min"
[output.grid]
lon0 = 60.0 ; lat0 = -55.0 ; nx = 120 ; ny = 100 ; dx = 0.25 ; dy = 0.25
heights = [0, 500, 1000, 2000, 3000, 4000, 6000, 8000, 10000, 12000]
```

### Warm start (`[run] checkpoint` / `restart`)

```toml
[run]
checkpoint = "12h"     # save the full state every 12 h AND at run end ("end" = end only)
# restart  = "results/myrun/checkpoint_myrun.jls"   # resume from a saved state
```

The checkpoint carries the complete parcel state, the seeded mask and the
physics accumulators (ice, sea spray).  To EXTEND a
finished run: keep the config identical, set a later `[run].stop`, point
`restart` at the checkpoint — the run resumes at the saved step instead of
re-simulating (releases never double-fire: the per-step quota only covers
remaining window steps).  Output files of the resumed segment get an
`_rMMDDHH` suffix so the original grids are kept; frames carry true
timestamps, so segments concatenate trivially.  Same-config continuation
only (pool ≥ saved; forward runs).  Tests: `test/test_checkpoint.jl`.

### Output estimators (gridding kernel, receptors, resolution)

```toml
[output]
kernel     = "adaptive"        # none | uniform (FLEXPART 4-cell) | adaptive
resolution = true              # write R(x) = mass-weighted kernel width per cell
# GRID/kNN-ONLY physical bandwidth floor: per-parcel deposit h ≥ Taylor
# σ(age) = σv·age → σv·√(2·T_L·age); gives just-released parcels (hband=0 until
# the next kNN refresh) an honest footprint.  NOT applied to the receptor — a
# non-renormalised point-KDE over-counts a dense near-source field (ETEX 8021
# A/B, ETEX_RESULTS.md §pair-level).  Neutral on the grid at benchmark N.
bandwidth_floor = "dispersion" # off (default) | dispersion
floor_sigma_v   = 0.5          # m/s → σ ≈ 0.7 km @ 1 h, 3.6 km @ 24 h
floor_tl        = "5min"
[output.receptors]             # gridless KDE time series at station coords
file = "/path/to/stations.txt" # DATEM list: id lat lon; → receptors_<name>.csv
layer_top = "100m"             # sampling layer (bottom 0 m)
# min_bandwidth = "5km"        # fixed representativeness floor (off by default)
```

### Sea-spray source (`[[release]] kind = "seaspray"`)

```toml
[[species]]
id   = "ss"
kind = "aerosol"
aerosol_type = "seasalt"

[[release]]
species   = "ss"
kind      = "seaspray"           # met-driven OSSA source (Ovadnevaite et al. 2014)
lon = [160.0, 200.0] ; lat = [-40.0, -10.0]   # ocean cells inside this box
time = [2025-05-21T03:00:00, 2025-05-22T03:00:00]
mass_per_parcel = "50kg"         # REQUIRED: the spawn quantum
  [release.seaspray]
  gust_sigma = 1.5               # sub-grid/gust wind σ [m/s]; 0 = mean-field only
  # sst      = 288.0             # K; omit = met 2t per cell
  # hs       = 2.0               # m; omit = Pierson–Moskowitz from U10
  # cd       = 1.5e-3            # omit = bulk (0.8+0.065·U10)e-3
```

The five-modal OSSA flux (wave Reynolds number Re_Hw = u*·H_s/ν_w — wind,
wave state, SST/salinity through the seawater viscosity) is evaluated every
transport step on the ocean cells (`lsm < 0.5`) of the box, from the archived
`10u/10v` (+ `2t`); the accumulated dry-mass per cell spawns surface parcels
of one quantum each.  With `gust_sigma > 0` the flux is integrated over
N(U10, σ) (Gauss–Hermite) — the TRANSIENT contribution beyond the mean
fields, dominant below ~10–15 m/s mean wind where the production threshold
truncates the wind PDF.  The run-end budget line reports the mean-field and
gust-integrated emissions separately (`transient_enhancement` = their ratio).
Physics + mode table: `easlee/src/physics/seaspray.jl`; tests:
`test/test_seaspray.jl`.

### In-flight species conversion (`[species.convert]`)

```toml
[[species]]
id = "SO2"
kind = "gas"
molar_mass = 64.06           # g/mol; needed for the auto mass_ratio when the id is not in the species DB
  [species.convert]
  to  = "SO4a"          # daughter species (defined below)
  tau = "25d"           # e-folding; fit against the observed SO2 decay
  # mass_ratio = 1.5    # omit = auto M_to/M_from from the molar masses

[[species]]
id = "SO4a"
kind = "aerosol"
aerosol_type = "sulfate"
molar_mass = 96.06
density = 1770.0
```

Stochastic Lagrangian conversion: each parcel switches species with
probability 1−exp(−dt/τ) per step (per-parcel deterministic RNG stream),
rescaling its mass by M_to/M_from — mass-conserving in expectation, and the
daughter parcels get the full aerosol machinery (bins, settling, deposition,
`output.aod` Mie AOD fields) from the moment of conversion via auto-created
shadow release boxes.  The HT use: SO2 → stratospheric sulfate → model AOD
vs the Maïdo lidar.  Budget logged at run end; carried through warm-start
checkpoints.  Tests: `test/test_species_convert.jl`.

### Ice freezing + sedimentation for tagged water (`[species.ice]`)

```toml
[[species]]
id   = "H2O_volc"
kind = "gas"
  [species.ice]
  rhi_crit = 1.5             # freeze only above this RHi (homogeneous nucleation, Koop)
  r_eff    = "20um"          # effective ice radius → Stokes–Cunningham v_t at parcel T,p
  h_layer  = "500m"          # settling escape scale (loss rate v_t/h_layer)
  # q0     = 0.002           # release mixing ratio [kg/kg]; omit = from the plume ensemble
```

Transport-phase ice microphysics for a tagged water species — the
first-48-hour stratospheric loss the eruption-column washout cannot see
(it ends at detrainment).  Every sync step, each parcel evaluates its water
mixing ratio q = mass/M_air against saturation over ICE (Murphy & Koop 2005)
at the met (T, p); above `rhi_crit` × q_si the excess condenses and sediments
out on the fall time scale `h_layer`/v_t — removed from the tagged budget
(the stratosphere→troposphere ice flux; totals logged at run end).  Retention
is bounded in the [q_si, rhi_crit·q_si] band (nucleation-threshold
hysteresis).  M_air (the air mass the parcel represents) comes from the plume
ensemble's detrainment water fraction (`plume_rise mode="parcels"`) or from
the prescribed `q0`.  Consistent with `[physics.mixing]` (mass exchange
between equal-air-mass parcels is mixing-ratio exchange).  Presence of the
table enables the process (`ice = true` sugar works; the GUI writes an
explicit `enabled` key).  Co-scavenging BY the ice: a partner species
declares

```toml
[[species]]
id = "SO2"
kind = "gas"
  [species.ice_scav]
  of     = "H2O_volc"        # must name an ice-enabled species
  uptake = 0.3               # retained-in-ice fraction (HNO3 ≈ 1; SO2 smaller)
```

and loses `mass · uptake · f_shed(cell)` wherever the named species sheds
ice in the same met cell — the removed totals are logged per box at run end
next to the water number.
Tests: `test/test_ice_physics.jl`.

### Gravitational fall of a condensed gas-kind species (`[species.settling]`)

```toml
[[species]]
id   = "ICE_b1"
kind = "gas"                 # a tagged condensed population, not an aerosol size bin
  [species.settling]
  r_eff = "5um"              # effective particle radius → Stokes–Cunningham v_t at the parcel's T,p (0 = off)
  rho   = 917.0              # particle density [kg/m³] (ice)
```

Gives the ice/ash shadow population of a `[species.convert]` chain its own
falling trajectory, so it sublimates where it ARRIVES rather than where it
formed.  `[species.fall]` is accepted as an older alias of the same table;
the run log says `species settling armed ([species.settling])`.  For an
aerosol-kind species the same table form carries the size-bin settling toggle
as `enabled = true|false` (the scalar `settling = true` still works).

### In-plume radiative cooling (`[physics.plume_radiation]`)

```toml
[physics.plume_radiation]
species  = "H2O_volc"       # excess-water species (needs the [species.ice] air-mass closure)
c_rad    = 0.002            # K/day per ppmv excess — CALIBRATE (see below)
q_bg     = 5.0              # ppmv background subtracted before cooling
feedback = true             # buoyancy-equilibrium descent of the cooled parcel
partners = ["sulfate"]     # species that FEEL the co-located water cooling without
                           # producing it: partner parcels get the met-cell mean
                           # excess-q cooling + the same buoyancy descent — the
                           # aerosol/water SEPARATION emerges as the aerosol settles
                           # out of the moist layer and its co-cooling shuts off
tprime_tau = "1d"          # OPTIONAL: relaxation time of the prognostic cold
                           # anomaly T′.  Cooling charges T′ (the parcel is
                           # colder than ambient AT ITS LEVEL); T′ relaxes into
                           # the buoyancy descent; the [species.convert] and
                           # [species.ice] saturation gates evaluate
                           # e_si(T_ambient + T′).  Omit/0 = previous: instant
                           # equilibrium, and the self-cooling NEVER reaches
                           # the gates (the hybrid vertical keeps the parcel at
                           # ambient T by construction).  Steady state
                           # T′ ≈ −Q_rad·τ; physical scale = the mesoscale
                           # buoyancy adjustment (L/h)/N ≈ 0.5–1 day.
```

Linearised cooling-to-space from the parcel's own excess water:
dT/dt = −c_rad·(q − q_bg), with q = mass/airm.  Each parcel accumulates its
ΔT_rad along the trajectory (mass-weighted stats at run end — compare to the
few-K in-plume cold anomalies in the HT radiosonde analyses); with
`feedback` the cooled parcel sinks through the ambient θ profile by
dz = dθ/(∂θ/∂z) — the slow radiative descent of the water layer.  The
default c_rad gives −0.3 K/day at 150 ppmv excess (≈ −3 K over ten days —
the observed anomaly class), but it is a one-line closure, not a radiation
code: calibrate against published heating rates before citing.
Tests: `test/test_plume_radiation.jl`.

### Inter-trajectory mixing during transport (`[physics.mixing]`)

```toml
[physics.mixing]
scheme   = "pairwise"          # off (default) | pairwise
K_h      = 1000.0              # horizontal eddy diffusivity of the exchange [m²/s]
sigma0_h = ""                  # pair kernel half-width; empty/0 = met grid spacing
sigma_v  = "300m"              # vertical pair half-width; empty/0 = turbulence yl
background = true              # coverage-based edge entrainment of ambient air
every    = 1                   # apply every N transport steps
```

By default every parcel keeps its tracer mass forever — turbulence moves
positions but never exchanges composition between trajectories.
`scheme = "pairwise"` turns on Fickian tracer-mass exchange between
co-located same-species parcels (Gaussian pair weights, second-moment
renormalisation to `K_h`, θ-limiter, exact mass conservation — derived below),
plus — with `background = true` —
a coverage-based relaxation toward zero-tracer ambient air: parcels fully
surrounded by plume neighbours exchange only among themselves, while edge and
isolated parcels entrain ambient air at rate `(1−cov)·K_h/σ₀²` (the physical
leak of tracer into untracked air).  Applies to release parcels only
(`relbox > 0`; domainfill background parcels are skipped).  In a measured
44-h counterfactual at K_h = 1000 m²/s the plume footprint and compactness were
unchanged, peaks ~10 % smoother, and ~9 % of the mass was entrained into ambient
air.  Config tests: `test/test_transport_mixing.jl`.

**How the exchange is computed** (`easlee/src/physics/mixing.jl`).

1. *Cell list.*  A 2-D spatial hash bins the active parcels into cells of side
   `nsig × σ₀_h` (default `nsig = 3`), so the neighbours of a parcel are looked
   up in the 3 × 3 block around it.  Cost O(N), not O(N²).
2. *Pair weights.*  For every pair (i, j) inside `nsig × σ₀_h`,

   ```
   λ⁰_ij = α · K_h / (π σ_h⁴) · exp( −r_h²/(2σ_h²) − r_v²/(2σ_v²) )
   ```

   with `r_h`, `r_v` the horizontal and vertical separations, `σ_h = σ₀_h`, and
   `σ_v = ½(yl_i + yl_j)` from the Lagrangian length scale of the last
   advection step.  `max_neighbors` (default 50) caps the pair list so that a
   point release cannot make the cost quadratic.
3. *Second-moment renormalisation.*  The weights are rescaled so that
   `Σ_j λ_ij r²_h,ij = 4 K_h` for every parcel.  A Taylor expansion of
   `dq_i/dt = Σ_j λ_ij (q_j − q_i)` gives `½ Σ_j λ_ij (Δx² + Δy²) ∇²q`, so the
   rescaling makes the scheme reproduce `∂q/∂t = K_h ∇²q` for any smooth parcel
   distribution, whatever λ⁰ was.
4. *Conservative exchange.*  Each pair exchanges `F_ij = λ_ij (q_j − q_i)`,
   added to i and subtracted from j, so `Σ_i q_i` is conserved exactly.
   Threaded, with per-thread accumulators summed at the end.
5. *θ-limiter.*  After the fluxes are accumulated, the largest `θ ∈ [0,1]` with
   `q_i + θ·dq_i ≥ 0` for every parcel and every species is applied.
   Positivity is guaranteed and conservation preserved.
6. *Background relaxation* (`background = true`).  From the pre-renormalisation
   second moment `mom2_i = Σ_j λ⁰_ij r²_h,ij` come a coverage fraction
   `cov_i = min(mom2_i σ₀_h² / (4 α K_h), 1)` and a rate
   `β_i = (1 − cov_i) K_h / σ₀_h²`.  A parcel in the interior has `cov ≈ 1` and
   exchanges only with its neighbours; an edge or isolated parcel has `cov < 1`
   and relaxes toward the ambient value at `β_i`.

A turbulent K either moves parcels or exchanges mass between them, never both.
Keep `[physics.mixing]` off whenever a horizontal trajectory term such as
`kh_ft` is armed; `build_run` refuses the combination.

### Observed-track nudging + wind-bias drift (`[physics.nudge]`, `[physics.wind_bias]`)

```toml
[physics.nudge]                # put a cloud where the observations say (HT)
file = "results/ht_cloud_envelope.csv"  # time, lon_west, lon_east, lat_south, lat_north
tau  = "3h"                    # centroid relaxation time
z    = [18000.0, 40000.0]      # only parcels in this height band [m]

[physics.wind_bias]            # manual stratospheric-wind tuning knob
u = 0.0                        # constant eastward drift [m/s] in the band
v = 0.0                        # northward
z = [18000.0, 40000.0]
```

`nudge` rigid-shift relaxes the in-band parcel ensemble's horizontal CENTROID
toward the time-interpolated observed track — every in-band parcel gets the
SAME increment, so the cloud's internal structure and relative dispersion are
untouched; only the advection bias is corrected.  Outside the track's time
span the run is free (nudge while observed, free-run after — e.g. nudge to
the Jan 15–16 stereo/MLS track, then let ERA5 carry the cloud to the
Australian radiosondes).  `wind_bias` adds a constant per-parcel drift in the
band (the archived met is never modified).  Dateline-safe (unwrapped-lon
arithmetic); `test/test_track_nudge.jl` (14/14).

## Meteorology — automatic AVAILABLE

`available = "auto"` scans the met directory at run time and builds the
date→file mapping itself:

- reads the authoritative valid time from each file's header (GRIB
  validity / NetCDF time), with a filename fallback;
- filters companion files (CLD, PT, SP1…) — only main files are listed;
- recurses subdirectories (e.g. GFS `date/cycle/`);
- de-duplicates overlapping forecast cycles, keeping the freshest run.

Point `available` at an explicit `AVAILABLE` file to override.

## Meteorology — domain crop (`crop`)

Loading a global archive for a run confined to a box wastes memory: with
`crop = "auto"` (the default) a **domain-fill** run loads only the release box
**+ 2 grid cells** of every met field.  This is *exact*
for domain-fill — the per-step cull already discards any parcel leaving the
box, so met outside box + halo is provably never touched — and cuts resident
memory by the grid-area ratio (a Pacific box on global 0.25° ERA5: ~10 % of
the grid, ~3–4× lower peak RSS).  A box crossing the dateline is handled
(the cropped axis continues monotonically past 180°).

Ordinary release runs are **never auto-cropped** — trajectories can travel
arbitrarily far.  If you know yours stay inside a region, set an explicit
window `crop = [lon0, lon1, lat0, lat1]`: it applies to any run, with
**regional-domain semantics — parcels reaching the crop edge are terminated
there**, exactly as if the met archive ended at that boundary.  `crop = "off"`
disables cropping entirely.  The output grid must lie inside the cropped
window (checked at startup).  GRIB met only; other formats run uncropped with
a warning.

See [MET_DATA.md](MET_DATA.md) for per-dataset details and the variables
namelists in `listvariables/`.

## Releases

A release defines where/when/what particles enter (forward) or the receptor
(backward).  `geometry` selects the shape (the total `particles`/`mass` are
shared across the resulting point sources):

| geometry | lon / lat | meaning |
|---|---|---|
| `box` (default) | `[a,b]` / `[c,d]` | a lon/lat/z volume |
| `point` | scalar / scalar | a single point, z range `[z1,z2]` |
| `profile` | scalar / scalar | one release at each level in `z = [...]` |
| `line` | vertex list | sample `npoints` along the polyline (ship/flight track) |
| `polygon` | vertex list | sample `npoints` uniformly inside (city/fire perimeter) |

```toml
[[release]]
geometry  = "line"
lon = [76.0, 78.0, 80.0] ; lat = [-40.0, -39.0, -41.0] ; z = ["300m", "1500m"]
npoints   = 100
time = [...] ; particles = 20000 ; mass = "1kg" ; species = "tracer"
```

Or read locations from a file instead of a geometry:

- **a text file** of locations: `from_file = "locs.txt"`, `file_format =
  "latlonlev"` (loose `lat lon level [weight]`, comma or space, comments ok);
- **a gridded NetCDF emission inventory** (CAMS-GLOB-* style: lat × lon × time,
  one variable per sector in kg m⁻² s⁻¹): `from_file = "...nc"`, `file_format =
  "inventory"`.  Sector fluxes are summed, area-weighted, cropped to the
  `lon`/`lat` box, and thresholded to the dominant cells; total `mass` is
  renormalised over the kept cells (emission conserved).  Raise `particles`
  to resolve more (weaker) cells.

`species` may be a list to release the same locations for several species.
See [OUTPUT_GRID.md](OUTPUT_GRID.md) for the output variables (per species /
per box, concentration / mixing ratio, backward sensitivity).

### Release `kind`

Besides the plain emission above (`kind = "release"`, the default), a release
can be one of three kinds:

```toml
[[release]]
kind = "release"        # release | domainfill | inventory
```

- `domainfill` — fill the lon/lat/z **box** with `particles` parcels at t=0
  (uniform in lon/lat, model levels drawn ∝ air mass ρ·Δz), each carrying
  `box air mass / particles` kg.  Every step the 4 box edges are replenished
  proportional to the incoming mass flux ρ·v·Δz·Δt, and parcels that leave the
  box (horizontally or above `z[end]`) are culled, keeping the population
  ≈ `particles`.  Forward-only; one domainfill box per run; the box may cross
  the dateline (`lon = [150, 230]`).
  E.g. `configs/era5_domainfill.toml`.
- `inventory` — release from a gridded emission inventory: shorthand for
  `from_file = "...nc"` + `file_format = "inventory"`.
- `volcanic` — time-resolved eruption column at the vent (lon/lat scalars,
  z ASL).  `[release.volcanic]` holds the satellite-constrained injection
  table `phases = [[t0, t1, z_bottom, z_top, frac], …]` (or `from_file` CSV)
  and the per-phase vertical mass profile `profile = "suzuki"` (shape
  `suzuki_a`, peak at (1−1/A)·z_top) | `"umbrella"` (Gaussian at
  z_b + 0.85·Δz) | `"uniform"`.  Mass accepts Tg-scale units ("0.4Tg").
  Mastin H↔MER helpers: `mastin_height_km` / `mastin_mer_kg_s`.

## Physics

| Group | Options | Docs |
|---|---|---|
| turbulence | none / hanna / thomson (+ mesoscale) | — |
| convection | off / massflux / cape (`timestep` is accepted but NOT implemented: convection runs every `run.timestep`; a different value only warns) | [CONVECTION.md](CONVECTION.md) |
| vertical | kinematic / hybrid / isentropic | [ISENTROPIC.md](ISENTROPIC.md), [VERTICAL_COORDS.md](VERTICAL_COORDS.md) |
| deposition (+ aerosol) | dry / wet, settling, AOD, κ-growth | [AEROSOLS.md](AEROSOLS.md) |

### Horizontal turbulence from an archived TKE field (`horizontal_tke`)

The Langevin schemes disperse parcels **vertically only**: `advance.jl` sets
`uout = vout = 0`, so the sole horizontal spreading is resolved-flow
deformation.  In coherent flow that means no spreading at all — measured on an
AROME 1.3 km back-trajectory in laminar trade wind, the cross-wind σ of a
5000-parcel cloud stayed at 0.15–0.19 km over 60 km of travel.

A **native-e** driver (one whose `listvariables/*.vars` declares
`tke_kind = e`: AROME native `e`, MesoNH) archives the full
e = ½(σ_u²+σ_v²+σ_w²), so the horizontal variance is already in the profile the
vertical step reads.  `horizontal_tke` runs a horizontal OU velocity off it, on
the same sub-step:

```toml
[physics.turbulence]
scheme         = "thomson"
horizontal_tke = true      # default false
```

Closure (`turbulence.jl` `tke_partition_h`), with `vmix_rw` / `vmix_ck` from the
same `#info` block as the vertical step:

    σ_u² = σ_v² = 2·e·r_h,  r_h = (1 − r_w)/2      (so r_w + 2·r_h = 1 exactly)
    K_h  = C_K·ℓ·√e,  τ_Lh = K_h/σ_u²  (capped at 600 s)
    u' ← u'·r + σ_u·√(1−r²)·N(0,1),   r = exp(−Δt/τ_Lh)

The horizontal share is whatever the vertical share leaves over — **not an
independent tunable**.  The OU form is exact for constant coefficients over any
Δt, and `u'`/`v'` persist across sync steps in `uturb`/`vturb`.

`build_run` refuses the key when the met is a Hanna-fill driver (its `tke` array
is 0.5·σ_w² and carries no horizontal information), when no turbulence scheme is
armed, when `[physics.mixing]` is also on (one K either moves parcels or
exchanges mass, never both), and on the KA/GPU advance, which does not implement
it.  The alternatives remain `[physics.turbulence] kh_ft` (FLEXPART's constant
free-troposphere `d_trop = 50 m² s⁻¹` memoryless random walk, applied above the
PBL top only) and the pairwise `[physics.mixing]` — pick one route.

In a **nested** run each nest uses its own closure (`NestRunConfig.vmix_*`,
read from that nest's variables namelist): AROME native `e` inside an ERA5
Hanna-fill outer domain is the supported case, and `horizontal_tke` is then
active on the AROME nest only.  A parcel crossing between nests has its stored
turbulent velocities reset, since they were drawn under the other nest's σ.

### Trajectory CSV layout (`[output.trajectory].format`)

```toml
[output.trajectory]
enabled = true
format  = "full"      # full (default) | demovor
```

`full` writes a headered CSV — `time,id,lon,lat,height_m,nest,mass_kg,age_s,
rho_kgm3,pblh_m`.  `demovor` writes the DEMOVOR/AROME column set with no
header and `%.5g` reals: `yyyy-mm-dd,HH:MM:SS,id,lon,lat,height_m,tke,pblh_m`
(the eighth column, `pblh_m`, is appended after the seven unchanged DEMOVOR
fields so the PBL / free-troposphere split can be computed from the text —
`height_m` is AGL, so `height_m ≤ pblh_m` is the boundary-layer test).  Rows
are sampled on the parcel's **own** nest; a parcel that has left an inner
domain is placed with the outer grid's projection, not the inner one's.

### Mesoscale velocity fluctuations

The PBL schemes only represent boundary-layer turbulence; the archived winds
also carry an unresolved SUB-GRID/mesoscale component in all layers
(FLEXPART's usig/vsig/wsig).  Enable the FLEXPART-equivalent per-parcel OU
velocity — memory `exp(−2Δt/T_met)` (T_met = met archive interval), amplitude
`fturbmeso × σ` of the 2×2×2×2 surrounding grid-point winds, vertical in
level-index space with ground reflection:

```toml
[physics.turbulence]
scheme    = "hanna"
mesoscale = true       # default false
fturbmeso = 0.16       # FLEXPART default; empirical (Stohl et al. 2005) —
                       # consider lowering for fine/hourly met (double-counting)
```

Validated in `validation/test_mesoscale.jl` (neighbourhood σ analytic, OU
stationary std = fturbmeso·σ, lag-1 autocorrelation = r).  CPU advance only —
the KA/GPU path ignores it (warned at build).

### Polar-stereographic transport (`polar_stereo_lat`)

On a lat/lon grid the horizontal metric degenerates at the poles: the
east-displacement → Δlon conversion carries 1/cos(lat) (→ ∞), and a parcel
crossing the pole cannot be represented by a monotone latitude step (it was
killed at the top grid row).  Following FLEXPART (switchnorth/switchsouth),
poleward of a switch latitude the position update runs in a pole-tangent
stereographic plane — conformal, no singularity, exact pole crossing (the
parcel comes out at lon+180°):

```toml
[physics]
polar_stereo_lat = 70.0   # default 75 (FLEXPART value); ≥ 90 disables
```

Auto-enabled per pole and per nest only where meaningful: a zonally-periodic
Rectilinear met grid whose latitude rows reach that pole (e.g. the global
ERA5 archive `/data/ERA5_full`, −90…90, variables namelist `ERA5_FULL` with
`companions = _EA_>_SFC_`).  Parcels equatorward of the switch that overshoot
the pole row in one step are mirrored across the pole (FLEXPART HSO/AL fix),
not killed.  Works on both the CPU and the KA/GPU advance.  Demo:
`configs/era5_polar.toml` (release at 85°N); validated in
`validation/test_polar_stereo.jl` (identity, lat/lon-metric equivalence at 76°,
N/S pole crossings against exact great-circle geometry, mirror path).

## Benchmark

`[run] benchmark = true` prints a per-group table at the run end (met read,
advection, convection, wet/dry deposition, settling, output …):
time, %, allocations, calls.  Useful for spotting the bottleneck (often I/O
when prefetch is off).

## How it works

The schema (`easlee/src/config/schema.jl`) is the single source of truth: it
drives TOML parsing, CLI `--help`, the GUI form, and validation.  `build_run`
(`easlee/build_run.jl`) turns a spec into the runtime and physics objects,
reusing the tested ingestion layer.
