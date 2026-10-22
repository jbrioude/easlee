
# ── RunSpec — declarative, typed run specification ────────────────────────────
#
# The single source of truth for a run.  Parsed from TOML (config/toml_io.jl),
# imported from FLEXPART files (config/flexpart_import.jl), or built by the GUI;
# turned into the runtime objects by build_run (config/build_run.jl).
#
# Every field has a default so a minimal run is a few lines.  Durations are
# stored in seconds, lengths in metres, masses in kg (parsed from human units
# at load time — see Units).  Grouped by concern, not by legacy file.

module RunSpecs

using Dates

export InitSource, InitSpecies
export RunSpec, MetSource, Species, ReleaseBox, Physics, Turbulence, Convection,
       VerticalSpec, Deposition, TransportMixing, OutputSpec, OutputGrid, Chemistry,
       validate_runspec, parse_refractive_index

# ── Meteorology (one per nest, finest first) ──────────────────────────────────
Base.@kwdef mutable struct MetSource
    model::String                       # ERA5, MERRA2, GFS, AROME_HDF, MESONH, WRF, AROME
    path::String                        # met directory
    variables::String                   # listvariables/<variables>.vars  (or absolute path)
    available::String = "auto"          # "auto" → scan path, else AVAILABLE file path
    # Domain crop: load only the release box + 2 halo cells of each met field.
    # "auto" = crop when the run has a domainfill release (exact there — the
    # per-step cull keeps parcels in the box); "off" = never; "explicit" with
    # crop_box = [lon0, lon1, lat0, lat1] forces a window on any run (parcels
    # reaching the crop edge are terminated — regional-domain semantics).
    crop::String = "auto"               # auto | off | explicit
    crop_box::Vector{Float64} = Float64[]
    # Refuse to run if the met has a coverage gap (missing files) larger than
    # this WITHIN the run window — a gap that big almost always means a missing
    # met file, not an intended overrun, and would be silently interpolated /
    # held-constant across it.  Small overruns (≤ this) still hold the last
    # frame.  Set to 0 to disable the guard (a deliberate sparse-archive run).
    max_met_gap_s::Float64 = 6 * 3600.0   # 6 h
end

# ── Species / material (fully inline) ─────────────────────────────────────────
Base.@kwdef mutable struct Species
    id::String                          # referenced by release boxes
    kind::String = "gas"                # "gas" | "aerosol"
    specnum::Int = 0                    # optional FLEXPART species number (metadata)
    # aerosol-only (ignored for gases)
    size_distribution::String = "lognormal"  # "lognormal" | "bins"
    dg_m::Float64 = 1.0e-6              # lognormal median diameter [m]
    sigma_g::Float64 = 2.0             # geometric std-dev
    nbins::Int = 8
    bin_edges_m::Vector{Float64} = Float64[]  # explicit bins (size_distribution="bins")
    density::Float64 = 2000.0          # particle density [kg/m³]
    kappa::Float64 = -1.0              # κ-Köhler (<0 → default from aerosol_type)
    aerosol_type::String = "dust"      # dust|sulfate|bc|oc|seasalt|ash|water (κ + refractive index defaults)
    refractive_index::String = "auto"  # "auto" → from aerosol_type, or "1.53-0.003i"
    settling::Bool = true
    # process toggles common to both kinds
    dry_deposition::Bool = false
    wet_deposition::Bool = false
    decay_halflife_s::Float64 = 0.0    # 0 → stable
    # material properties (0/-1 → filled from the species DB by `id`, see species_db.jl)
    molar_mass::Float64 = 0.0          # g/mol  (0 → from DB, else 29)
    reldiff::Float64  = -1.0           # gas dry-dep diffusivity ratio D(H2O)/D(gas)
    henry::Float64    = -1.0           # gas dry-dep Henry's-law constant [M/atm]
    f0::Float64       = -1.0           # gas dry-dep reactivity (0–1)
    dryvel::Float64   = -1.0           # prescribed gas dry-dep velocity [m/s] (<0 → resistance model)
    weta_gas::Float64 = -1.0           # below-cloud scavenging A
    wetb_gas::Float64 = -1.0           # below-cloud scavenging B
    # ── [species.ice]: transport-phase freezing + sedimentation ─────────────
    # Tagged water species: above ice_rhi_crit × saturation over ICE
    # (Murphy–Koop) the excess vapour condenses and sediments out on the fall
    # time scale h_layer/v_t(r_eff) — the first-48-h stratospheric loss the
    # eruption-column washout cannot see.  The parcel mixing ratio comes from
    # the plume ensemble's detrainment water fraction (plume_rise
    # mode="parcels") or from ice_q0 [kg/kg] otherwise.
    ice::Bool = false
    ice_rhi_crit::Float64 = 1.5        # homogeneous-nucleation RHi threshold (Koop)
    ice_r_eff_m::Float64 = 20.0e-6     # effective ice radius → Stokes–Cunningham v_t(T,p) per parcel
    ice_h_layer_m::Float64 = 500.0     # settling escape scale [m]
    ice_q0::Float64 = 0.0              # prescribed release mixing ratio [kg/kg]; 0 = plume
    # ── [species.ice_scav]: co-scavenging BY the ice of another species ─────
    # This species rides the co-located freezing water (`of` names an
    # ice-enabled species): dm = m · uptake · f_shed(cell), uptake ∈ (0,1] the
    # retained-in-ice fraction (HNO3 ≈ 1; SO2/HCl smaller, T-dependent).
    ice_scav_of::String = ""
    ice_scav_uptake::Float64 = 0.0
    # ── [species.convert]: in-flight species conversion (e.g. SO2 → sulfate) ─
    # Each parcel of this species switches to `to` with probability
    # 1 − exp(−dt/τ) per step (stochastic, mass-conserving in expectation);
    # the mass is rescaled by mass_ratio (0 = auto, M_to/M_from from the
    # species molar masses — SO2→sulfate ≈ 96/64).  The daughter gets its own
    # (auto-created) release boxes, so aerosol physics / deposition / AOD
    # apply from the moment of conversion.
    convert_to::String = ""
    convert_tau_s::Float64 = 0.0
    convert_mass_ratio::Float64 = 0.0
    # THERMODYNAMIC gate on the conversion ([species.convert].mode):
    #   ""          → plain τ (legacy)
    #   "freeze"    → only where q > rhi_crit·q_si   (vapour → ice)
    #   "sublimate" → only where q < q_si            (ice → vapour)
    # The pair makes the phase change reversible, so ice that falls into warmer
    # subsaturated air gives its water back instead of being lost for good.
    convert_mode::String = ""
    convert_rhi_crit::Float64 = 1.5
    # sub-grid lamina width: relative spread of q about the parcel mean.
    # 0 = legacy step gate (bit-identical).  See species_physics.jl.
    convert_lam_sigma::Float64 = 0.0
    # SECTIONAL ice ([species.convert].sectional): carry a per-PARCEL crystal
    # radius, shrink it on the exact r² = r₀² − 2K·dt law (mass fraction lost
    # = 1 − (r'/r)³), demote the parcel to the next smaller ICE bin when the
    # radius crosses the bin edge, and let settling read the same radius.
    # This is the only mechanism that produces SIZE SELECTIVITY: the lamina
    # and T′ both act on (1−S), which scales every bin alike (ht_full3, 8 d:
    # r_eff flat 3.10–3.44 µm).  Set on the sublimating (ICE) species.
    convert_sectional::Bool = false
    # ── [species.settling]: gravitational fall of a CONDENSED species ────────
    # Gives the ice/ash shadow population its own falling trajectory: it
    # sublimates where it ARRIVES, not where it formed.
    settle_r_eff::Float64 = 0.0        # effective particle radius [m] (0 = off)
    settle_rho::Float64 = 917.0        # particle density [kg/m³] (ice)
end

# ── Release ───────────────────────────────────────────────────────────────────
# A release defines where/when/what particles enter the run.  In forward runs
# it is an emission; in backward runs it is the receptor.  Geometry is either a
# lon/lat/z box (default), or a list of point locations read from a text file /
# emission inventory (`from_file`).
Base.@kwdef mutable struct ReleaseBox
    species::String                     # Species.id (single; multi handled at parse time)
    kind::String = "release"            # release | domainfill | inventory | volcanic
                                        #   release    — parcels emitted from the geometry over the window
                                        #   domainfill — fill the lon/lat/z box ∝ air mass at t=0, replenish
                                        #                the box edges from the inflow flux, cull leavers
                                        #   inventory  — gridded emission inventory (from_file, NetCDF or text)
                                        #   volcanic   — time-resolved eruption column at (lon,lat):
                                        #                injection phases (t0,t1,z_bottom,z_top,frac) with a
                                        #                Suzuki/umbrella vertical mass profile (see [release.volcanic])
    geometry::String = "box"            # box | point | line | polygon | profile
    lon::Vector{Float64} = Float64[]    # box [lon1,lon2]; point scalar; line/polygon vertices
    lat::Vector{Float64} = Float64[]
    z::Vector{Float64}   = [0.0, 100.0] # box/line/polygon [z1,z2]; profile = level list
    zkind::String = "agl"               # "agl" | "asl" | "pressure"
    time::Vector{DateTime} = DateTime[] # [start, stop]
    particles::Int = 10000              # TOTAL particles for this release
    mass_kg::Float64 = 1.0              # TOTAL released mass
    mass_per_parcel_kg::Float64 = 0.0   # parcel mass quantum (>0 → particles = ⌈mass/q⌉,
                                        # overrides `particles`; consistent resolution
                                        # across releases of very different mass)
    npoints::Int = 100                  # sampling density for line / polygon
    comment::String = ""
    # File / inventory source (overrides the geometry when set)
    from_file::String = ""              # path to a locations / inventory file
    file_format::String = "latlonlev"   # latlonlev | lonlatlev | inventory (lat lon value)
    # Volcanic column (kind="volcanic"): satellite-constrained injection phases.
    # Each phase = (start, stop, z_bottom [m ASL], z_top [m ASL], mass fraction);
    # within a phase the mass is distributed over the column by `volc_profile`.
    volc_phases::Vector{Tuple{DateTime,DateTime,Float64,Float64,Float64}} =
        Tuple{DateTime,DateTime,Float64,Float64,Float64}[]
    volc_from_file::String = ""         # CSV t0,t1,z_bottom,z_top,frac (alt. to phases)
    # Optional per-phase source pulses (multi-blast scenarios — e.g. the Hunga
    # Tonga TWO CLOUDS: climactic MER ~5e9–2e10 → stratospheric umbrella;
    # sustained MER ~1e8 → tropopause-capped cloud).  One (mer, T0, water_frac)
    # per phase, parsed from phase rows of length 8; empty → the plume_* values
    # apply to every phase.
    volc_pulses::Vector{NTuple{3,Float64}} = NTuple{3,Float64}[]
    volc_profile::String   = "suzuki"   # suzuki | umbrella | uniform
    volc_suzuki_a::Float64 = 4.0        # Suzuki (1983) shape A (larger → top-heavier)
    volc_suzuki_l::Float64 = 1.0        # Suzuki λ exponent
    # Plume rise (injection height computed at release time from the met column)
    # — overrides z over the plume layer.  Biomass burning: sofiev|freitas (FRP).
    # Volcanic: scheme="volcanic" — 1-D moist eruption column driven by the
    # source energetics (MER, mixture temperature, vaporised-water fraction).
    plume_scheme::String   = "none"     # none | sofiev | freitas | volcanic
    plume_frp_w::Float64   = 0.0        # fire radiative power [W]
    plume_fire_area::Float64 = 0.0      # active-fire area [m²] (freitas)
    plume_mer::Float64     = 0.0        # volcanic: mass eruption rate [kg/s] (total mixture)
    plume_temp_k::Float64  = 800.0      # volcanic: mixture T after magma–water interaction [K]
    plume_water_frac::Float64 = 0.05    # volcanic: vapour mass fraction at the vent
    plume_profile::String  = "uniform"  # uniform | gaussian — parcel fill over the layer
    plume_bottom::String   = "pbl"      # pbl | ground | fraction:<f>
    # mode="parcels" (volcanic only): each released trajectory runs the
    # Lagrangian plume phase (buoyant rise + entrainment + washout) through the
    # met column and enters transport at its EMERGENT detrainment height, its
    # mass scaled by the washout survival of the box species (SO2/HCl/H2O by
    # id; others unscaled).  The [release] mass is then the AT-VENT mass, not
    # the delivered stratospheric mass.
    plume_mode::String     = "layer"    # layer | parcels
    plume_washout::Bool    = true       # parcels mode: in-plume + umbrella washout
    plume_entrainment::String = "meanfield"  # parcels mode: meanfield | pairwise
                                        # (pairwise = v2: edge engulfment + parcel-
                                        # parcel exchange → hot core / cool edge)
    # Umbrella-intrusion spreading (parcels mode): after detrainment the cloud
    # spreads as a gravity current R(t)=(3λNQ_V t²/2π)^(1/3) fed by the
    # column's volume flux (from the SAME 1-D physics) — the archived met has
    # no eruption outflow, so without this the first hours are far too narrow.
    plume_umbrella::Bool = false
    plume_spread_lambda::Float64 = 0.2  # intrusion constant λ (Costa et al. 2013)
    plume_umbrella_windshield::Bool = false  # shield in-cloud parcels from the
                                        # resolved wind by s = min(1, u_r/|u_amb|)
                                        # (density current expels the ambient flow)
    plume_umbrella_rsource::Float64 = 10000.0  # source/collapse-dome radius [m]:
                                        # the source-fed interior flow saturates
                                        # inside it (HT-class ≈ 40 km)
    # Handoff to the resolved met: the radial push stops once Ṙ < spread_umin or
    # past spread_tmax.  These END the gravity current, so they SET the cloud's
    # final width — for HT the 5 m/s default fires at 8 h (R = 432 km) and the
    # 12 h cap never bites.  Physically umin is not "the wind is faster" (the
    # wind ADVECTS the cloud, it does not stop the current spreading) but "the
    # intrusion is too thin/weak to keep expanding"; keep it below the ambient
    # |v| at the neutral level if the meridional width matters.
    plume_spread_umin::Float64 = 5.0    # [m/s]
    plume_spread_tmax::Float64 = 43200.0  # [s]
    # Umbrella retention (post-detrainment washout).  `qc_crit` is an UPDRAFT
    # holding capacity; the umbrella has w ≈ 0, so setting `umb_r_med > 0`
    # replaces it by sedimentation of a lognormal crystal spectrum through the
    # umbrella (see PlumeWashout).  This is the parameter that SETS the mass
    # delivered to the stratosphere — for HT, qc_crit = 5 g/kg ≡ r ≈ 105 µm.
    plume_qc_crit::Float64  = 5.0e-3    # [kg/kg]
    # Vapour cap at the umbrella equilibration [× ice saturation at the
    # ambient T]: retained vapour above cap·q_si leaves as the COARSE-mode
    # budget (fcoarse), to be released as a large-crystal species (ICE_C*)
    # instead of arriving as mislabelled supersaturation.  Inf = off.
    plume_vapor_cap::Float64 = Inf
    plume_umb_r_med::Float64 = 0.0      # [m] 0 ⇒ legacy qc_crit behaviour
    plume_umb_r_sg::Float64 = 2.0
    plume_umb_h::Float64    = 2000.0    # [m]
    plume_umb_t::Float64    = 3600.0    # [s]
    plume_solid_r::Float64  = 0.0       # [m] coarse-ash fallout radius, 0 ⇒ off
    # z-dependent spread amplification (backward-footprint constrained):
    # u_r × (1 + zgain·ramp(z; z0→z1)).  0 = off.
    plume_spread_zgain::Float64 = 0.0
    plume_spread_zg_z0::Float64 = 26.0e3
    plume_spread_zg_z1::Float64 = 30.0e3
    # azimuth-dependent spread amplification (inversion-wing constrained):
    # u_r × (1 + g·cos²θ_N·ramp(z; az_z0→az_z1)), g = azgain_n/_s per hemisphere
    # of the vent.  spread_azgain sets both; _n/_s override.  az_z1 ≤ az_z0 ⇒
    # no ramp (all z).  0 = off.
    plume_spread_azgain_n::Float64 = 0.0
    plume_spread_azgain_s::Float64 = 0.0
    plume_spread_az_z0::Float64 = 0.0
    plume_spread_az_z1::Float64 = 0.0
    plume_spread_hmin::Float64 = 0.0    # [m] PHYSICAL handoff: stop once the
                                        # intrusion has thinned to h = V/(πR²)
                                        # below this (no longer a distinct
                                        # layer).  0 ⇒ use umin/tmax only.
    # Radial magma–seawater vent split ([release.plume_rise.vent_mix], parcels
    # mode): the jet is magmatic on the axis, ocean-water-loaded at the edge —
    # x(r) = x_core + (x_edge−x_core)·r^p, per-parcel (T,qv,qc) from the
    # phreatomagmatic enthalpy balance (VOLCANIC_THEORY.md §3).  Overrides
    # plume_temp_k / plume_water_frac (the column runs on the area means).
    plume_vent_mix::Bool = false
    plume_magma_t::Float64 = 1323.0     # magma temperature [K]
    plume_magma_water::Float64 = 0.03   # magmatic volatile wt fraction of melt
    plume_sw_core::Float64 = 0.15       # seawater mass fraction on the axis
    plume_sw_edge::Float64 = 0.70       # … at the edge
    plume_mix_exp::Float64 = 2.0        # radial profile exponent
    # ── [release.seaspray] (kind = "seaspray"): OSSA met-driven source ───────
    # The box reserves the species/relbox; parcels spawn from the Ovadnevaite
    # et al. (2014) flux over the ocean cells inside the box (lsm < 0.5).
    # mass_per_parcel is REQUIRED (the spawn quantum).
    spray_sst::Float64 = -1.0           # SST [K]; < 0 → met 2t per cell
    spray_salinity::Float64 = 35.0      # ‰ (viscosity → Re_Hw)
    spray_hs::Float64 = -1.0            # significant wave height [m]; < 0 → Pierson–Moskowitz
    spray_cd::Float64 = -1.0            # drag coefficient; < 0 → bulk (0.8+0.065·U10)e-3
    spray_rho::Float64 = 2017.0         # dry sea-salt density [kg/m³]
    # ── [release.dust] (kind = "dust"): Kok et al. (2014) mineral-dust source ──
    # The box reserves the species/relbox; parcels spawn from the Kok flux over
    # the LAND cells of the box (lsm ≥ 0.5) from u* (met "ustar").  The species
    # should use size.distribution = "kok" (Kok 2011 emitted size distribution).
    # mass_per_parcel is REQUIRED (the spawn quantum).
    dust_clay::Float64 = 0.2            # clay mass fraction of the soil (0–1)
    dust_bare::Float64 = 1.0            # bare (erodible) fraction of the cell (0–1)
    dust_w_pct::Float64 = 0.0           # gravimetric soil moisture [%] when the met has no swvl1
    dust_d_m::Float64 = 127.0e-6        # threshold grain diameter [m] (Shao & Lu 2000)
    dust_cd0::Float64 = 4.4e-5          # Kok 2014 C_d0
    dust_ce::Float64 = 2.0              # Kok 2014 C_e
    dust_calpha::Float64 = 2.7          # Kok 2014 C_α
    dust_clay_file::String = ""         # NetCDF lon/lat map of the clay fraction (optional)
    dust_bare_file::String = ""         # NetCDF lon/lat map of the bare fraction (optional)
    # ── [release.chem_init] — release-scoped chemistry initial condition ──────
    # Parcels born from THIS release start with the plume composition, while
    # every other parcel keeps the [chemistry.init] background — the spatial
    # eruption injection the global so2_init/H2O_ppb scalars cannot express.
    # Keys are intensive (mixing ratio / mass loading), applied at each parcel's
    # own (T,P):  "so2_ppb" [ppb]→cs.so2, "so4_ugm3" [µg/m³]→cs.so4, "h2o_ppm"
    # [ppmv]→cs.h2o (and pins it), any mechanism species name → its column [ppb].
    chem_init::Dict{String,Float64} = Dict{String,Float64}()
    # Volcanic plume: seed the prognostic chemistry water DIRECTLY from the column
    # model's detrained water fraction (qdet = mass/airm) instead of a fixed
    # h2o_ppm — the thermodynamically-determined stratospheric water the plume
    # delivers.  Set via [release.chem_init].h2o_from_plume = true (needs
    # plume_rise mode=parcels, which sets airm, and [chemistry].h2o_prognostic).
    chem_h2o_plume::Bool = false
end

# ── Physics ───────────────────────────────────────────────────────────────────
Base.@kwdef mutable struct Turbulence
    scheme::String = "hanna"            # "none" | "hanna" | "thomson"
    hanna_onthefly::Bool = false        # Hanna option b
    sigw_ft::Float64 = 0.1              # free-troposphere σ_w [m/s]
    # FLEXPART free-atmosphere horizontal diffusivity D_h [m²/s] (Stohl et al.
    # 2005: d_trop = 50).  Random walk √(2KΔt) above the PBL, trajectory route
    # only — [physics.mixing] must stay off when this is on (double counting).
    # 0 = off (legacy: no horizontal turbulence above the PBL at all).
    kh_ft::Float64 = 0.0
    # TKE-driven horizontal turbulence: an OU (Langevin) horizontal velocity
    # from the SAME archived e the vertical step uses, σ_u²=σ_v²=2·e·(1−r_w)/2.
    # Native-e drivers only (AROME native e, MesoNH) — a Hanna-fill driver
    # archives w-variance only and build_run refuses the option there.
    # Trajectory route: [physics.mixing] must stay off (double counting).
    horizontal_tke::Bool = false
    # Hanna (1982) horizontal Langevin in the PBL — σ_u, σ_v, T_Lu, T_Lv from the
    # same (u*, w*, h, L) the vertical Hanna fill uses; FLEXPART's in-PBL
    # horizontal turbulence.  Hanna-fill drivers only (ERA5/GFS/MERRA-2/
    # AROME-HDF); native-e nests use horizontal_tke instead.  Trajectory route:
    # [physics.mixing] must stay off (double counting).  kh_ft may stay armed
    # for the free troposphere (FLEXPART d_trop), as in FLEXPART.
    horizontal_hanna::Bool = false
    # Mesoscale velocity fluctuations (FLEXPART usig/vsig/wsig): per-parcel OU
    # velocity with memory exp(−2Δt/T_met) and amplitude fturbmeso·σ of the
    # surrounding grid-point winds — sub-grid wind uncertainty in ALL layers
    # (the PBL schemes above only cover boundary-layer turbulence).
    mesoscale::Bool = false
    fturbmeso::Float64 = 0.16           # FLEXPART default; empirical (Stohl 2005) —
                                        # consider lowering for fine/hourly met
end

Base.@kwdef mutable struct Convection
    scheme::String = "off"              # "off" | "massflux" | "cape"
    cape_min::Float64 = 100.0
    cin_max::Float64 = 50.0
    tau_mix_s::Float64 = 0.0
    timestep_s::Float64 = 0.0           # convection op-split step; 0 → = run.timestep
end

Base.@kwdef mutable struct VerticalSpec
    mode::String = "kinematic"          # "kinematic" | "hybrid" | "isentropic"
    sigma_r::Float64 = 0.3
    sigma_t::Float64 = 0.1
    diabatic::Bool = true
    pt_winds::String = "auto"           # "auto" | "off" | directory
    pt_pattern::String = "IFS_PT_MET"
end

Base.@kwdef mutable struct Deposition
    dry::Bool = false                   # global toggles (species flags can refine)
    wet::Bool = false
end

# ── Inter-trajectory mixing during transport ([physics.mixing]) ──────────────
# Pairwise Fickian exchange of tracer mass between co-located parcels of the
# SAME species (physics/mixing.jl, derived in docs/CONFIG.md), plus a
# coverage-based relaxation toward zero-tracer ambient air: parcels fully
# surrounded by plume neighbours (cov ≈ 1) exchange only among themselves;
# edge / isolated parcels (cov < 1) entrain ambient air at rate
# (1−cov)·K_h/σ₀².  Independent of the stochastic turbulence scheme (which
# moves positions but never exchanges tracer between trajectories).
Base.@kwdef mutable struct TransportMixing
    scheme::String = "off"              # "off" | "pairwise"
    K_h::Float64 = 1000.0               # horizontal diffusivity [m²/s]
    sigma0_h::Float64 = 0.0             # kernel half-width [m]; 0 = met grid spacing
    sigma_v::Float64 = 0.0              # vertical pair half-width [m]; 0 = turbulence yl
    nsig::Float64 = 2.0                 # pair-search radius in σ₀_h units
    max_neighbors::Int = 30
    background::Bool = true             # coverage-based edge entrainment of ambient air
    every::Int = 1                      # apply every N transport steps
end

Base.@kwdef mutable struct Physics
    turbulence::Turbulence = Turbulence()
    # ── [physics.plume_radiation]: LW cooling from excess plume water ───────
    # Linearised cooling-to-space dT/dt = −c_rad·(q_ppmv − q_bg) per parcel
    # (needs the ice-closure airm); optional buoyancy-equilibrium descent.
    # c_rad in K/day per ppmv — CALIBRATE against published heating rates.
    prad_species::Vector{String} = String[]   # empty = off
    prad_c_rad::Float64 = 0.002               # K/day per ppmv excess
    # LW opacity [1/(ppmv*m)] for the thickness-aware saturation; 0 = legacy
    # linear law + hard cap.  See physics/plume_radiation.jl.
    prad_lw_kappa::Float64 = 0.0
    prad_h_ref::Float64 = 1000.0              # fallback layer thickness [m]
    prad_q_bg::Float64 = 5.0                  # ppmv
    prad_feedback::Bool = true                # radiative descent
    # Relaxation time [s] of the prognostic cold anomaly T′: cooling charges
    # T′ (the parcel is colder than ambient AT ITS LEVEL), T′ relaxes into the
    # buoyancy descent, and the freeze/sublimate gates evaluate e_si at
    # T_ambient + T′.  0 = legacy (instant equilibrium, no anomaly state).
    # Physical estimate: mesoscale adjustment (L/h)/N ≈ 0.5–1 day.
    prad_tprime_tau::Float64 = 0.0
    prad_partners::Vector{String} = String[]  # species that FEEL the co-located
                                              # cooling without producing it (the
                                              # sulfate riding the moist plume) —
                                              # aerosol/water separation emerges
    # [physics].prognostic_airm : the air mass each parcel represents is
    # refreshed from the LOCAL parcel spacing (kNN) instead of being frozen at
    # release, so q = mass/airm (ice closure + plume_radiation) DILUTES as the
    # cloud disperses — read from the trajectories, never prescribed.
    prognostic_airm::Bool = false
    prognostic_airm_every::Int = 12     # transport steps between refreshes
    # [tracers]: passive transport-history air-mass tags (strat residence, …)
    tags_enabled::Bool = false
    tags_strat_theta::Float64 = 380.0         # θ threshold for "stratospheric" [K]
    tags_tropo_origin::Bool = false           # plume-origin (tropospheric) air-fraction
                                              # tag: 1.0 on volcanic-release parcels,
                                              # diluted toward 0 by transport mixing →
                                              # gridded f_tropo (needs mixing enabled)
    convection::Convection = Convection()
    vertical::VerticalSpec = VerticalSpec()
    deposition::Deposition = Deposition()
    mixing::TransportMixing = TransportMixing()
    ctl::Float64 = 5.0
    ifine::Int = 4
    # ── Polar-stereographic transport ([physics].polar_stereo_lat) ──────────
    # Poleward of ±this latitude, horizontal parcel-position updates run in a
    # pole-tangent stereographic plane instead of lat/lon (FLEXPART
    # switchnorth/switchsouth: no 1/cos(lat) blow-up, exact pole crossing).
    # Auto-enabled per pole on global Rectilinear met grids whose latitude
    # rows reach that pole (e.g. /data/ERA5_full, −90…90).  ≥ 90 disables.
    polar_stereo_lat::Float64 = 75.0
    # Monin-Obukhov surface-layer wind ([physics].surface_layer_wind): the
    # reader constant-extends the lowest model level's wind to the ground, so
    # below ~20 m a parcel carries the lowest-level wind instead of the
    # similarity profile (1.8x too fast at 2 m over land).  On = rescale onto
    # u(z) = U1 - (u*/k)[ln(z1/z) - psi_m(z1/L) + psi_m(z/L)], anchored on the
    # model's own u* (physics/surface_layer.jl).
    surface_layer_wind::Bool = true
    # ── Observed-track ensemble nudging ([physics.nudge], HT thread) ─────────
    # Rigid-shift relaxation of the in-band parcel ensemble's horizontal
    # CENTROID toward a time-interpolated observed cloud track (CSV: time,
    # lon_west, lon_east, lat_south, lat_north — the ht_cloud_envelope.csv
    # format).  Preserves the cloud's internal structure (every in-band parcel
    # gets the SAME increment); corrects only the advection bias; free-runs
    # outside the track's time span.  Off when file is empty.
    nudge_file::String = ""
    nudge_tau_s::Float64 = 10800.0      # relaxation time [s]
    nudge_z::Vector{Float64} = [18000.0, 40000.0]   # height band [m]
    # ── Constant wind-bias correction ([physics.wind_bias]) ─────────────────
    # Per-parcel drift u/v [m/s] added inside the height band — the manual
    # stratospheric-wind tuning knob (applied to parcel motion, the archived
    # met is untouched).  Off when both are 0.
    windbias_u::Float64 = 0.0
    windbias_v::Float64 = 0.0
    windbias_z::Vector{Float64} = [18000.0, 40000.0]
end

# ── Output ────────────────────────────────────────────────────────────────────
# Horizontal CRS is EXPLICIT (legacy inferred it from the met projection):
#   crs = "lonlat"  → lon0/lat0 in degrees, dx/dy in degrees
#   crs = "native"  → lon0/lat0 = x0/y0 in metres of nest `nest`'s projection,
#                     dx/dy in metres
# Vertical is one of four coordinates; `heights` holds the level VALUES in the
# units implied by `vkind` (m AGL | hPa | model-level index | K).
Base.@kwdef mutable struct OutputGrid
    crs::String = "lonlat"          # "lonlat" | "native"
    nest::Int = 1                   # (native) which nest's projection defines x/y metres
    lon0::Float64 = -180.0          # origin x  (degrees if lonlat, metres if native)
    lat0::Float64 = -90.0           # origin y
    nx::Int = 360
    ny::Int = 180
    dx::Float64 = 1.0               # spacing x (degrees if lonlat, metres if native)
    dy::Float64 = 1.0               # spacing y
    vkind::String = "height_agl"    # height_agl | pressure | model_level | theta
    heights::Vector{Float64} = Float64[0, 100, 500, 1000, 3000, 10000]  # level values (units per vkind)
end

# Per-trajectory CSV dump (one row per sampled parcel per interval): position +
# time, mass + age, chemistry species (if chem on), and met at the parcel.
Base.@kwdef mutable struct TrajectoryOutput
    enabled::Bool = false
    interval_s::Float64 = 0.0          # row cadence; 0 → = output.sample
    max_parcels::Int = 10000           # subsample cap (0 → all parcels)
    species::Vector{String} = String[] # chem columns (empty → all run species)
    # Row layout:
    #   "full"    — headered CSV: time,id,lon,lat,height_m,nest,mass_kg,age_s,rho,pblh
    #   "demovor" — the legacy DEMOVOR/AROME column set, no header, %.5g reals:
    #               yyyy-mm-dd,HH:MM:SS,id,lon,lat,height_m,tke,pblh_m
    format::String = "full"            # full | demovor
end

Base.@kwdef mutable struct OutputSpec
    grid::OutputGrid = OutputGrid()
    quantity::String = "concentration" # "concentration" | "mixing_ratio" | "both"
    name::String = ""                  # output file base name (empty → run name)
    # Output gridding kernel (ADAPTIVE_KERNEL design):
    #   "none"     — nearest-cell deposition
    #   "uniform"  — FLEXPART uniform kernel (4-cell bilinear; comparison reference)
    #   "adaptive" — per-parcel truncated Gaussian, bandwidth = kNN spacing:
    #                resolution is bought with trajectories, not the grid; the
    #                3-h age rule disappears (h→0 at a fresh point release).
    kernel::String = "uniform"
    resolution::Bool = false           # write R(x) = mass-weighted mean kernel width
                                       # per cell (needs kernel = "adaptive"):
                                       # R < Δx_out → grid resolves the run;
                                       # R > Δx_out → sampling noise, raise N.
    # Gridless receptor output ([output.receptors], needs kernel="adaptive"):
    # pointwise KDE evaluation at station coordinates — the field's value at
    # the coordinate at the model's NATIVE resolution (no cell-average
    # representativeness error), plus R(x_s) from the same sum.  One time
    # series per station (sample cadence, window average — FLEXPART
    # semantics), written to <output.path>/receptors_<name>.csv.
    receptors_file::String = ""        # DATEM station list: "id lat lon" per line
    receptors_layer::Vector{Float64} = [0.0, 100.0]  # sampling layer [m] (all stations)
    receptors_hmin::Float64 = -1.0     # bandwidth floor [m] for the point KDE
                                       # (subgrid-footprint variance control);
                                       # -1 → half the output grid spacing, 0 → pure
                                       # sample-point estimator
    # Physical bandwidth floor (adaptive kernel + receptors): the kNN hband is
    # the SAMPLING spacing; in a young dense plume it collapses to metres —
    # far below the width of atmosphere a parcel represents after turbulent
    # transport (single-parcel receptor spikes, one-cell filaments).
    # "dispersion" floors per-parcel h at the Taylor growth
    # σ(age) = σv·age → σv·√(2·T_L·age) using floor_sigma_v / floor_tl.
    bandwidth_floor::String = "off"    # "off" | "dispersion"
    floor_sigmav::Float64 = 0.5        # [m/s] effective relative-dispersion velocity
    floor_tl_s::Float64 = 300.0        # [s] Lagrangian time scale of the growth
    # Flow-evolved covariance kernel (kernel = "covariance", [output.covariance]):
    # each parcel carries a 2×2 horizontal covariance P [m²] evolved every sync
    # step by the resolved deformation P ← (I+∇u·dt)·P·(I+∇u·dt)ᵀ (wind shear
    # stretches/rotates it) plus turbulent growth 2·σv²·TL·dt per axis, and P
    # is read only when a concentration is evaluated (grid deposit: axis σx/σy;
    # receptors: full bivariate Gaussian).  The width is the PHYSICAL footprint
    # of the parcel (puff-style), not the kNN statistical spacing — no kNN pass.
    # mode="physical": σ from P directly — P is ADDED width on top of the
    # ensemble's absolute dispersion (double-counts relative dispersion; kept
    # as the falsifiable variant, and it needs no kNN pass).
    # mode="shape": kernel AREA pinned to the kNN spacing hband² (statistical,
    # vanishes with N — no double counting), ANISOTROPY from P (shear-aligned
    # smoothing along the filament); receptors rescale the full tilted P to
    # det = h⁴.  Needs the kNN pass (same cadence as kernel="adaptive").
    cov_mode::String = "physical"      # "physical" | "shape"
    cov_sigma0::Float64 = 50.0         # [m] release footprint σ₀ (P init = σ₀²·I)
    cov_sigma_v::Float64 = 0.5         # [m/s] unresolved relative-dispersion velocity
    cov_tl_s::Float64 = 300.0          # [s] Lagrangian time scale (D = σv²·TL)
    cov_sigma_max::Float64 = 5.0e5     # [m] per-axis growth cap
    cov_stride::Int = 1                # evolve P every N sync steps (dt_eff = N·dt)
    receptors::Vector{@NamedTuple{id::String, lon::Float64, lat::Float64}} =
        @NamedTuple{id::String, lon::Float64, lat::Float64}[]   # inline stations
    per_release::Bool = true
    per_species::Bool = true
    # [output].group — what the accumulator's 4th dimension is:
    #   "species" : one output column per species (grid memory independent of
    #               release-box count — a per-box grid is ~18 MB/box on the HT
    #               grid and a 3400-box file release OOM'd at 70 GB)
    #   "release" : one column per release box (legacy; required to attribute
    #               backward sensitivities to individual boxes; "box" = alias)
    #   ""        : auto — "species" forward, "release" backward
    group::String = ""
    aod::Bool = false
    # Passive-reference comparison ([output.reference], CHEM_INIT.md §6): the
    # name of a [chemistry.init] source to compare the reacting chemistry
    # against.  When set, the engine carries a FROZEN passive twin of every
    # parcel's initial composition (advected, never reacted) and writes, per
    # species, X (reacted) / X_pass (transport-only) / X_ref (the source field
    # gridded) into the chemistry NetCDF, plus a per-frame bias/rmse line.
    # Empty = off.  Replaces the bespoke chem_pass/accumulate_all in the
    # retired chem_cams_pacific.jl.
    reference::String = ""
    # Per-parcel chemistry-only tendency ([output.tendencies]): a list of
    # mechanism species whose in-situ photochemical production/loss is
    # accumulated per parcel along its trajectory (transport excluded) and
    # gridded as <species>_tend [ppb] in the chemistry NetCDF.  Empty = off.
    tendencies::Vector{String} = String[]
    path::String = "juliaout/output/"
    interval_s::Float64 = 3600.0
    average_s::Float64 = 3600.0
    sample_s::Float64 = 900.0
    trajectory::TrajectoryOutput = TrajectoryOutput()
end

# ── Chemistry initialisation (CHEM_INIT.md) ───────────────────────────────────
# A named external field source ([chemistry.init.sources.<name>]).
Base.@kwdef mutable struct InitSource
    name::String                        # referenced by init.species / default
    kind::String = "grib"               # grib | netcdf | constant
    dir::String = ""                    # file directory (grib/netcdf)
    units::String = "mmr"               # mmr (kg/kg) | vmr | ppb
    vars::Dict{String,String} = Dict{String,String}()   # species → file variable
                                        # (mechanism default map applies otherwise)
end

# One per-species assignment ([[chemistry.init.species]]).
Base.@kwdef mutable struct InitSpecies
    name::String                        # mechanism species name
    from::String                        # a source name, or a derivation:
                                        #   correlation | steady_state | family |
                                        #   constant | profile
    predictor::String = ""              # correlation: predictor species
    curve::String = ""                  # correlation: x,y CSV (clamped at ends)
    total::String = ""                  # family: the family-total species/source
    regime::String = "any"              # steady_state: any | day (night → skip)
    value::Float64 = 0.0                # constant [ppb]
    profile::String = ""                # profile: z/p CSV
end

# ── Chemistry ─────────────────────────────────────────────────────────────────
# Gas-phase chemistry applied to active parcels after transport, via the
# ChemState/ChemConfig layer (chemistry/coupler.jl).  `mechanism` selects the
# reaction set; `init = "cams"` seeds parcel composition from a CAMS GRIB
# (chemistry/cams_background.jl) so EASLEE can be compared against CAMS with the
# same chemistry.  Disabled by default — a normal passive run ignores this block.
Base.@kwdef mutable struct Chemistry
    enabled::Bool = false
    mechanism::String = "reducedqss"    # reduced | reducedqss | bb12 | cb05 | strat
    timestep_s::Float64 = 3600.0        # chemistry op-split step [s] (must be × run.timestep)
    photolysis_timestep_s::Float64 = 0.0 # J-value refresh step; 0 → = chemistry.timestep
    init::String = "background"         # "background" (clean defaults) | "cams" (GRIB)
    boundary::Bool = true               # re-seed newly released parcels from CAMS (init=cams)
    lat_deg::Float64 = 45.0             # fallback latitude when per-parcel unavailable
    h2o_ppb::Float64 = 1.0e4            # water vapour mixing ratio [ppb]
    doy::Int = 0                        # day of year; 0 → derived from run start
    # Solver tolerances: defaults aligned with ChemConfig (radicals sit near/below
    # a 1e-3 abstol — OH ~1e-4 ppb, night NO3 ~5e-6 ppb — so the old 1e-5/1e-3
    # TOML defaults silently degraded radical accuracy; see PHOTOLYSIS_TUV.md §Audit.
    # abstol 1e-6 per the 2026-07-09 truth study: 24-h HO2 error 31 % at 1e-5 vs
    # 0.1 % at 1e-6, at equal-or-lower rb23 cost).
    reltol::Float64 = 1.0e-6
    abstol::Float64 = 1.0e-6
    # solver arithmetic for the bb12/cb05-rb23 fast paths: "fp64" (default) or
    # "fp32" (tolerances auto-clamped: bb12 ≥1e-4, cb05 ≥1e-5; docs/GPU_CHEMISTRY.md)
    precision::String = "fp64"
    # CB05 integrator: "rb23" (analytic-Jacobian Rosenbrock23, cb05_fast.jl —
    # 3–5× faster than FBDF at equal tolerances and at least as accurate;
    # GPU_CHEMISTRY.md §8.2) or "fbdf" (legacy generic fallback).  Ignored by
    # the other mechanisms (they have their own fast paths).
    solver::String = "rb23"
    # photolysis (CB05 mechanism): TUV lookup table or analytic cos-SZA
    photolysis::String = "analytic"     # "analytic" | "lut"
    lut_path::String = ""               # LUT NetCDF (empty → TUV/cb05_photolysis_lut.nc)
    o3col_DU::Float64 = 300.0           # LUT O3 column [DU]; CAMS-derived per-parcel when init=cams
    albedo::Float64 = 0.05              # LUT surface UV albedo
    cloud::Bool = true                  # cloud correction of J from the archived cloud
                                        # fraction (cc, CLD companions); clear-sky if absent
    # CAMS coupling (init = "cams")
    cams_dir::String = ""               # dir of CAMS_YYYYMMDD_HHMM.grib, or a single .grib path
    cams_species::Vector{String} = String[]   # CAMS shortNames to load (empty → all)
    # Surface emissions ([chemistry.emissions]) — inject gridded CAMS-GLOB-* inventory
    # fluxes into the mixing ratios of the domain-fill parcels in the surface layer.
    # Requires a [[release]] kind="domainfill".
    emissions_dir::String = ""          # inventory directory ("" → no emissions)
    emissions_files::Vector{String} = String[]  # basename glob patterns (empty → all *.nc)
    emissions_res::Float64 = 0.4        # working-grid resolution [deg]
    emissions_inject_m::Float64 = 150.0 # inject into parcels below this height [m]
    emissions_nox_no::Float64 = 0.9     # mole fraction of NOx emitted as NO (rest NO2)
    # diurnal cycle applied to the daily/monthly-mean flux (local-solar-time
    # scaling, 24-h mean preserved): "none" | "gaussian" | "mu2011"
    emissions_diurnal::String = "none"
    emissions_diurnal_peak::Float64 = 13.5   # peak local solar time [h] (gaussian)
    emissions_diurnal_sigma::Float64 = 2.0   # width [h] (gaussian)
    # Lightning NOx ([chemistry.lightning]): Price & Rind (1992) flash rate from
    # the archived updraft-mass-flux cloud top → NO into UT parcels.  `yield` is
    # the tunable mol-NO-per-flash (calibrate against the ~1 ppb FT-O3 signal).
    lightning_enabled::Bool = false
    lightning_yield::Float64 = 250.0         # mol NO per flash
    lightning_umf_thresh::Float64 = 0.01     # updraft mass flux threshold [kg m-2 s-1]
    lightning_h_min::Float64 = 5.0           # min cloud-top height for lightning [km]
    lightning_ut_lo::Float64 = 6000.0        # inject into parcels between these [m]
    lightning_ut_hi::Float64 = 14000.0
    # Parcel–parcel composition exchange during transport ([chemistry.mixing]):
    # without it every parcel is a CLOSED chemical box (no dilution-driven
    # chemistry — a concentrated plume never titrates/recovers).  "pairwise"
    # bins active seeded parcels by met grid cell × model level and relaxes
    # random r-adjacent pairs toward their mass-weighted mean with
    # μ = 1−exp(−dt/τ) — conservative in Σ m·c per species (same form as the
    # plume-phase v2 / puff-mixing engines).
    het_awet::Float64 = 0.0             # [chemistry.het] awet: background wet aerosol
                                        # surface density [m²/m³] for HO2/N2O5/SO2
                                        # uptake (marine ≈ 1e-4); 0 = off
    # SWV-prognostic mode (HISTORICC Task 2): freeze the met re-prescription of
    # per-parcel H2O after each slot's first seeding — the parcel carries its own
    # water budget (HT injection, CH4-oxidation source) along the trajectory
    h2o_prognostic::Bool = false
    # In-situ stratospheric H2O source from CH4 oxidation (strat mechanisms):
    # per-parcel CH4 decays against the solved OH/O(1D)/Cl; each CH4 returns
    # 2 H2O (Jones & Pyle 1986: H2O + 2·CH4 ≈ const).  Requires h2o_prognostic.
    ch4_oxidation::Bool = false
    # Sulfate-SAD-driven het (Analysis C / Evan 2023): parcel carries its own
    # sulfate (+SO2 → sulfate), whose wet surface area sets the strat het rates.
    sulfate_sad::Bool = false
    sulfate_gmr_m::Float64 = 0.35e-6
    sulfate_sigma_g::Float64 = 1.3
    sulfate_kappa::Float64 = 0.5
    sulfate_so2_tau_s::Float64 = 17.1*86400.0
    sulfate_so2_oh::Bool = false
    # PHYSICS-carried ice → chemistry ([species.convert] coupling): crystal
    # r_eff for A = 3·m/(ρ·r_eff), and the self-shielding length of the
    # Beer–Lambert τ = (A/2)·L on the strat photolysis.  ice_shield=0 ⇒ off.
    ice_reff_m::Float64   = 1.0e-5
    ice_shield_m::Float64 = 0.0               # actual SO2+OH+M reaction (JPL, prognostic OH)
    sulfate_so2_oh_feedback::Bool = false      # EXPERIMENTAL (OH is fast-QSS; see coupler.jl)
    so4_init_kgm3::Float64 = 0.0
    so2_init_kgm3::Float64 = 0.0
    # PSC / cold-point DEHYDRATION sink on the prognostic water: freeze-out of
    # the excess above RHi=psc_rhi, sedimenting on psc_tau — the cold-point entry
    # dehydration + the dominant polar-vortex SWV sink.  Requires h2o_prognostic.
    psc_dehydration::Bool = false
    psc_rhi::Float64 = 1.0
    psc_tau_s::Float64 = 6.0 * 3600.0
    # NAT DENITRIFICATION: permanent HNO3 removal below T_NAT (Hanson &
    # Mauersberger 1988 equilibrium) on the large-NAT sedimentation timescale —
    # delays chlorine deactivation, prolonging springtime ozone loss.
    psc_denitrification::Bool = false
    psc_nat_tau_s::Float64 = 5.0 * 86400.0
    mixing::String = "off"              # off | pairwise
    mixing_tau_s::Float64 = 1800.0      # relaxation timescale τ [s]
    mixing_cells::Int = 1               # neighbourhood size [met grid cells]
    mixing_levels::Int = 1              # neighbourhood size [model levels]
    # [chemistry.aerosols]: carry the CAMS aerosol fields (aermr01–20) on the
    # domain-fill parcels (Lagrangian aerosol transport), derive the per-parcel
    # wet surface area for the heterogeneous chain (replaces the scalar
    # het_awet), and output a 550-nm AOD field per frame.
    aerosols::Bool = false              # read + carry CAMS aermr* per parcel
    aerosols_het::Bool = true           # per-parcel A_wet → HO2/N2O5/SO2 uptake
    aerosols_aod::Bool = true           # AOD550 output field per frame
    aerosols_rh_growth::Bool = true     # RH-grown MEC in the overhead AOD550 column
                                        # (per-level κ-Köhler + wet Mie); false = dry v1
    # [chemistry.init] — declarative species initialisation (CHEM_INIT.md).
    # Sources name external fields; per-species assignments pick a source or a
    # derivation.  Legacy init="cams" + [chemistry.cams] expands to this.
    init_default::String = ""           # source for species not listed (""=zeros)
    init_sources::Vector{InitSource} = InitSource[]
    init_species::Vector{InitSpecies} = InitSpecies[]
end

# ── Top-level run ─────────────────────────────────────────────────────────────
Base.@kwdef mutable struct RunSpec
    name::String = "easlee_run"
    task::String = "transport"          # "transport" | "chemistry" | "validation"
                                        # — what scripts/easlee.jl dispatches to.
    direction::String = "forward"       # "forward" | "backward"
    wind::String = "3D"                 # 3D | 2D-lonlat | 2D-latalt | 2D-lonalt | 1D-vertical
    # ── Warm start ([run] checkpoint/restart) ────────────────────────────────
    # checkpoint = "end" (write <output.path>/checkpoint_<name>.jls at run end)
    # or a period ("6h": every 6 h AND at the end).  restart = path to a
    # checkpoint: the run resumes at the saved step with the saved parcels,
    # chemistry and physics accumulators — extend a finished run by setting a
    # later [run].stop instead of re-simulating.  Same-config continuation
    # only (pool ≥ saved, same releases/species); output files get an
    # _rMMDDHH suffix so the original grids are never overwritten.
    checkpoint::String = ""             # "" (off) | "end" | duration
    restart::String = ""                # checkpoint file to resume from
    # Per-parcel diagnostic snapshots (duration, e.g. "3h"; "" = off): periodic
    # dumps of parcel state + convert-gate state + flux accumulators to
    # <output.path>/snapshots/, for the freeze/demotion attribution and the
    # sectional-τ instrument.  Diagnostic only — never read back by the model.
    parcel_dump::String = ""            # "" (off) | duration
                                        # — freeze the suppressed parcel index/indices.
    start::DateTime = DateTime(2000,1,1)
    stop::DateTime = DateTime(2000,1,1,1)
    timestep_s::Float64 = 900.0
    benchmark::Bool = false             # per-group @time-style timing report
    seed::Int = 0                       # [run].seed: Random.seed! at run start (0 = unseeded)
    meteorology::Vector{MetSource} = MetSource[]
    species::Vector{Species} = Species[]
    release::Vector{ReleaseBox} = ReleaseBox[]
    physics::Physics = Physics()
    output::OutputSpec = OutputSpec()
    chemistry::Chemistry = Chemistry()
end

# ── Validation ────────────────────────────────────────────────────────────────
# Is `a` an integer multiple of base `b` (within rounding)?  (b ≤ 0 → false)
_ismult(a::Real, b::Real) = b > 0 && abs(a - round(a / b) * b) < 1e-6 * max(b, 1)

"""
    parse_refractive_index(s) → ComplexF64 | nothing

`species.refractive_index` string ("1.53-0.003i", "1.53-0.003im", "1.53") →
complex m = n − k·i (the AEROSOL_REFRACTIVE sign convention).  `nothing` for
"auto" (→ the aerosol_type default) and for an unparseable string — callers
that need to tell the two apart test `== "auto"` first.
"""
function parse_refractive_index(s::AbstractString)   # review 2026-09-04 A7 item 3
    t = replace(strip(String(s)), r"\s+" => "")
    lowercase(t) == "auto" && return nothing
    t = replace(replace(t, r"im$" => "i"), r"i$" => "im")
    return tryparse(ComplexF64, t)
end

# review 2026-09-04 A7 item 4: build_run constructs ONE AerosolBins from the FIRST
# aerosol species and shares it with every aerosol box — the bin-defining
# properties must therefore agree across aerosol species until per-species bins exist.
_aerosol_bin_key(sp::Species) = (sp.size_distribution,
    sp.size_distribution == "bins" ? sp.bin_edges_m : Float64[],
    sp.size_distribution == "bins" ? (0.0, 0.0, 0) : (sp.dg_m, sp.sigma_g, sp.nbins),
    sp.density, sp.aerosol_type, sp.kappa, sp.refractive_index)

function validate_runspec(s::RunSpec; check_paths::Bool=true)
    errs = String[]
    s.direction in ("forward", "backward") ||
        push!(errs, "direction must be forward|backward (got \"$(s.direction)\")")
    s.stop > s.start || push!(errs, "stop must be after start")
    s.timestep_s > 0 || push!(errs, "timestep must be > 0")
    isempty(s.meteorology) && push!(errs, "at least one [[meteorology]] source is required")
    isempty(s.release)     && push!(errs, "at least one [[release]] box is required")

    specids = Set(sp.id for sp in s.species)
    for sp in s.species
        sp.kind in ("gas", "aerosol") ||
            push!(errs, "species \"$(sp.id)\": kind must be gas|aerosol")
        sp.kind == "aerosol" && !(sp.size_distribution in ("lognormal", "bins", "kok")) &&
            push!(errs, "species \"$(sp.id)\": size.distribution must be lognormal|bins|kok")
        sp.kind == "aerosol" && sp.size_distribution == "bins" && length(sp.bin_edges_m) < 2 &&
            push!(errs, "species \"$(sp.id)\": bins distribution needs ≥2 bin_edges_m")
        sp.kind == "aerosol" && sp.refractive_index != "auto" && parse_refractive_index(sp.refractive_index) === nothing &&
            push!(errs, "species \"$(sp.id)\": refractive_index must be \"auto\" or n-ki like \"1.53-0.003i\" (got \"$(sp.refractive_index)\")")   # review 2026-09-04 A7 item 3
        sp.molar_mass >= 0 || push!(errs, "species \"$(sp.id)\": molar_mass must be ≥ 0 g/mol (0 = from the species DB)")
        if sp.ice
            sp.ice_rhi_crit >= 1.0 ||
                push!(errs, "species \"$(sp.id)\": ice.rhi_crit must be ≥ 1")
            sp.ice_r_eff_m > 0 || push!(errs, "species \"$(sp.id)\": ice.r_eff must be > 0")
            sp.ice_h_layer_m > 0 || push!(errs, "species \"$(sp.id)\": ice.h_layer must be > 0")
            sp.ice_q0 >= 0 || push!(errs, "species \"$(sp.id)\": ice.q0 must be ≥ 0 (kg/kg)")
        end
        if !isempty(sp.convert_to)
            sp.convert_tau_s > 0 ||
                push!(errs, "species \"$(sp.id)\": convert.tau must be > 0")
            sp.convert_to != sp.id ||
                push!(errs, "species \"$(sp.id)\": convert.to must differ from the species itself")
            any(t -> t.id == sp.convert_to, s.species) ||
                push!(errs, "species \"$(sp.id)\": convert.to must name a defined species (got \"$(sp.convert_to)\")")
            sp.convert_mass_ratio >= 0 ||
                push!(errs, "species \"$(sp.id)\": convert.mass_ratio must be ≥ 0 (0 = auto from molar masses)")
        end
        if sp.ice_scav_uptake != 0.0 || !isempty(sp.ice_scav_of)
            0.0 < sp.ice_scav_uptake <= 1.0 ||
                push!(errs, "species \"$(sp.id)\": ice_scav.uptake must be in (0, 1]")
            tgt = findfirst(t -> t.id == sp.ice_scav_of, s.species)
            (tgt !== nothing && s.species[tgt].ice) ||
                push!(errs, "species \"$(sp.id)\": ice_scav.of must name an ice-enabled species (got \"$(sp.ice_scav_of)\")")
            sp.ice && push!(errs, "species \"$(sp.id)\": a species cannot be both ice-enabled and an ice_scav partner")
        end
    end
    # review 2026-09-04 A7 item 4: one shared AerosolBins → aerosol species must agree
    let aer = filter(sp -> sp.kind == "aerosol", s.species)
        if length(aer) > 1 && length(unique(_aerosol_bin_key.(aer))) > 1
            push!(errs, "aerosol species " * join(["\"$(sp.id)\"" for sp in aer], ", ") *
                  " differ in size distribution / bin edges (dg, sigma_g, nbins) / density / aerosol_type / kappa / refractive_index, " *
                  "but build_run shares ONE AerosolBins (defined by the first aerosol species) across all aerosol boxes — " *
                  "per-species bins are not implemented; make them identical or run them in separate configs")
        end
    end
    for (k, r) in enumerate(s.release)
        r.species in specids ||
            push!(errs, "release[$k]: species \"$(r.species)\" not defined in [[species]]")
        length(r.time) == 2 || push!(errs, "release[$k]: time must be [start, stop]")
        r.kind in ("release", "domainfill", "inventory", "volcanic", "seaspray", "dust") ||
            push!(errs, "release[$k]: kind must be release|domainfill|inventory|volcanic|seaspray|dust")
        if r.kind == "dust"
            r.mass_per_parcel_kg > 0 ||
                push!(errs, "release[$k] (dust): mass_per_parcel must be set (the spawn quantum)")
            (0.0 <= r.dust_clay <= 1.0 && 0.0 <= r.dust_bare <= 1.0) ||
                push!(errs, "release[$k] (dust): clay and bare must be fractions in [0, 1]")
            r.dust_w_pct >= 0 || push!(errs, "release[$k] (dust): soil_moisture must be ≥ 0 %")
            (r.dust_d_m > 0 && r.dust_cd0 > 0) || push!(errs, "release[$k] (dust): threshold_d and cd0 must be > 0")
            for (f, nm) in ((r.dust_clay_file, "clay_file"), (r.dust_bare_file, "bare_file"))
                isempty(f) || !check_paths || isfile(f) || push!(errs, "release[$k] (dust): $nm not found: $f")
            end
        end
        if r.kind == "seaspray"
            r.mass_per_parcel_kg > 0 ||
                push!(errs, "release[$k] (seaspray): mass_per_parcel must be set (the spawn quantum)")
            r.spray_salinity > 0 ||
                push!(errs, "release[$k] (seaspray): salinity must be > 0 ‰")
        end
        r.plume_scheme in ("none", "sofiev", "freitas", "volcanic") ||
            push!(errs, "release[$k]: plume_rise.scheme must be none|sofiev|freitas|volcanic")
        r.plume_mode in ("layer", "parcels") ||
            push!(errs, "release[$k]: plume_rise.mode must be layer|parcels")
        r.plume_mode == "parcels" && r.plume_scheme != "volcanic" &&
            push!(errs, "release[$k]: plume_rise.mode=\"parcels\" requires scheme=\"volcanic\"")
        r.plume_entrainment in ("meanfield", "pairwise") ||
            push!(errs, "release[$k]: plume_rise.entrainment must be meanfield|pairwise")
        if r.plume_vent_mix
            r.plume_mode == "parcels" ||
                push!(errs, "release[$k]: plume_rise.vent_mix requires mode=\"parcels\"")
            0.0 <= r.plume_sw_core <= r.plume_sw_edge <= 0.95 ||
                push!(errs, "release[$k]: vent_mix needs 0 ≤ seawater_core ≤ seawater_edge ≤ 0.95")
            r.plume_magma_t > 500.0 ||
                push!(errs, "release[$k]: vent_mix magma_t must be > 500 K")
            0.0 <= r.plume_magma_water <= 0.2 ||
                push!(errs, "release[$k]: vent_mix magma_water must be in [0, 0.2]")
        end
        if r.kind == "volcanic"
            (length(r.lon) >= 1 && length(r.lat) >= 1) ||
                push!(errs, "release[$k]: volcanic needs the vent lon/lat")
            (isempty(r.volc_phases) && isempty(r.volc_from_file)) &&
                push!(errs, "release[$k]: volcanic needs [release.volcanic] phases or from_file")
            !check_paths || isempty(r.volc_from_file) || isfile(r.volc_from_file) ||
                push!(errs, "release[$k]: volcanic from_file not found: $(r.volc_from_file)")
            for (i, p) in enumerate(r.volc_phases)
                p[1] < p[2] || push!(errs, "release[$k]: volcanic phase $i: start must precede stop")
                p[3] < p[4] || push!(errs, "release[$k]: volcanic phase $i: z_bottom must be < z_top")
                p[5] > 0    || push!(errs, "release[$k]: volcanic phase $i: mass fraction must be > 0")
            end
            r.volc_profile in ("suzuki", "umbrella", "uniform") ||
                push!(errs, "release[$k]: volcanic profile must be suzuki|umbrella|uniform")
            continue
        end
        if r.kind == "domainfill"
            isempty(r.from_file) ||
                push!(errs, "release[$k]: domainfill uses the lon/lat/z box, not from_file")
            r.geometry == "box" ||
                push!(errs, "release[$k]: domainfill needs geometry=\"box\"")
            length(r.lon) == 2 || push!(errs, "release[$k]: domainfill lon must be [lon1, lon2]")
            length(r.lat) == 2 || push!(errs, "release[$k]: domainfill lat must be [lat1, lat2]")
            isempty(r.z) && push!(errs, "release[$k]: domainfill needs z (fill-layer top = z[end])")
            s.direction == "forward" ||
                push!(errs, "release[$k]: domainfill is forward-only (boundary inflow)")
            continue
        elseif r.kind == "inventory"
            isempty(r.from_file) &&
                push!(errs, "release[$k]: kind=\"inventory\" requires from_file")
        end
        if !isempty(r.from_file)
            !check_paths || isfile(r.from_file) ||
                push!(errs, "release[$k]: from_file not found: $(r.from_file)")
            r.file_format in ("latlonlev","lonlatlev","inventory") ||
                push!(errs, "release[$k]: file_format must be latlonlev|lonlatlev|inventory")
        else
            r.geometry in ("box","point","line","polygon","profile") ||
                push!(errs, "release[$k]: geometry must be box|point|line|polygon|profile")
            if r.geometry == "box"
                length(r.lon) == 2 || push!(errs, "release[$k]: box lon must be [lon1, lon2]")
                length(r.lat) == 2 || push!(errs, "release[$k]: box lat must be [lat1, lat2]")
                length(r.z)   == 2 || push!(errs, "release[$k]: box z must be [z1, z2]")
            elseif r.geometry in ("point","profile")
                (length(r.lon) >= 1 && length(r.lat) >= 1) ||
                    push!(errs, "release[$k]: $(r.geometry) needs lon and lat")
                r.geometry == "profile" && length(r.z) < 1 &&
                    push!(errs, "release[$k]: profile needs a z level list")
            else  # line | polygon
                (length(r.lon) >= 2 && length(r.lon) == length(r.lat)) ||
                    push!(errs, "release[$k]: $(r.geometry) needs matching lon/lat vertex lists (≥2)")
                r.geometry == "polygon" && length(r.lon) < 3 &&
                    push!(errs, "release[$k]: polygon needs ≥3 vertices")
            end
        end
    end
    s.physics.turbulence.scheme in ("none","hanna","thomson") ||
        push!(errs, "turbulence.scheme must be none|hanna|thomson")
    # review 2026-09-04 A7 item 7: the sub-step count dt/ctl is exact only for an integer ctl
    (isfinite(s.physics.ctl) && isinteger(s.physics.ctl) && s.physics.ctl >= 1) ||
        push!(errs, "physics.ctl must be a positive integer value (got $(s.physics.ctl)): the turbulence sub-step is run.timestep/ctl and a non-integer ctl leaves a partial step")
    # review 2026-09-04 A7 item 7: "one K moves parcels OR exchanges mass, never both".  A WARNING
    # (not an error, unlike kh_ft in build_run) because the shipped volcanic_mixing.toml arms both
    # deliberately and the two act at different scales: the mesoscale OU velocity carries the
    # unresolved variance of the ARCHIVED winds (σ of the surrounding grid points, memory = met
    # interval), the pairwise K_h is an intra-plume exchange at the kernel scale — a partial
    # overlap in sub-grid horizontal dispersion, not the same K armed twice.
    s.physics.turbulence.mesoscale && s.physics.mixing.scheme != "off" &&
        @warn "[physics.turbulence].mesoscale = true AND [physics.mixing].scheme = \"$(s.physics.mixing.scheme)\": both represent unresolved horizontal dispersion (one moves parcels, one exchanges mass) — a partial double count; keep [physics.mixing] off unless the overlap is intended"
    let m = s.physics.mixing
        m.scheme in ("off","pairwise") ||
            push!(errs, "physics.mixing.scheme must be off|pairwise")
        if m.scheme == "pairwise"
            m.K_h > 0        || push!(errs, "physics.mixing.K_h must be > 0")
            m.sigma0_h >= 0  || push!(errs, "physics.mixing.sigma0_h must be ≥ 0 (0 = auto)")
            m.sigma_v >= 0   || push!(errs, "physics.mixing.sigma_v must be ≥ 0 (0 = turbulence yl)")
            m.every >= 1     || push!(errs, "physics.mixing.every must be ≥ 1")
        end
    end
    s.physics.convection.scheme in ("off","massflux","cape") ||
        push!(errs, "convection.scheme must be off|massflux|cape")
    s.physics.vertical.mode in ("kinematic","hybrid","isentropic") ||
        push!(errs, "vertical.mode must be kinematic|hybrid|isentropic")
    s.output.quantity in ("concentration","mixing_ratio","both") ||
        push!(errs, "output.quantity must be concentration|mixing_ratio|both")
    g = s.output.grid
    # degree values that cannot be degrees are projection metres given without
    # crs = "native" (the WRF/AROME confusion): say so instead of producing
    # negative cell areas or an empty grid (2026-10-13).  Only values that are
    # clearly not degrees are refused: a pole-centred global grid legitimately
    # has its last row half a cell past 90° (configs/merra2.toml: lat0 = −90,
    # dy = 0.5, ny = 361).
    if g.crs == "lonlat"
        (g.dx > 0 && g.dy > 0 && g.dx <= 360 && g.dy <= 180 && abs(g.lat0) <= 90 + g.dy) ||
            push!(errs, "output.grid (crs = \"lonlat\") is not a longitude/latitude box: lat0=$(g.lat0), " *
                        "dy=$(g.dy), ny=$(g.ny), dx=$(g.dx), nx=$(g.nx) — if these are projection metres " *
                        "of a Lambert nest (WRF, AROME), set [output.grid] crs = \"native\"")
    end
    g.crs in ("lonlat","native") ||
        push!(errs, "output.grid.crs must be lonlat|native")
    g.vkind in ("height_agl","pressure","model_level","theta") ||
        push!(errs, "output.grid.vkind must be height_agl|pressure|model_level|theta")
    g.nest >= 1 || push!(errs, "output.grid.nest must be ≥ 1")
    (g.nx > 0 && g.ny > 0) || push!(errs, "output.grid nx/ny must be > 0")
    isempty(g.heights) && push!(errs, "output.grid needs at least one vertical level")
    s.output.kernel in ("none", "uniform", "adaptive", "covariance") ||
        push!(errs, "output.kernel must be none|uniform|adaptive|covariance (got \"$(s.output.kernel)\")")
    s.output.bandwidth_floor in ("off", "dispersion") ||
        push!(errs, "output.bandwidth_floor must be off|dispersion (got \"$(s.output.bandwidth_floor)\")")
    s.output.bandwidth_floor == "dispersion" && (s.output.floor_sigmav <= 0.0 || s.output.floor_tl_s <= 0.0) &&
        push!(errs, "output.bandwidth_floor=\"dispersion\" needs floor_sigma_v > 0 and floor_tl > 0")
    s.output.resolution && !(s.output.kernel in ("adaptive", "covariance")) &&
        push!(errs, "output.resolution=true requires output.kernel=\"adaptive\" or \"covariance\" (R(x) is the kernel-width field)")
    if s.output.kernel == "covariance"
        s.output.cov_mode in ("physical", "shape") ||
            push!(errs, "output.covariance.mode must be physical|shape (got \"$(s.output.cov_mode)\")")
        s.output.cov_sigma0 > 0.0   || push!(errs, "output.covariance.sigma0 must be > 0")
        s.output.cov_sigma_v >= 0.0 || push!(errs, "output.covariance.sigma_v must be ≥ 0")
        s.output.cov_tl_s > 0.0     || push!(errs, "output.covariance.tl must be > 0")
        s.output.cov_stride >= 1    || push!(errs, "output.covariance.stride must be ≥ 1")
    end
    has_rec = !isempty(s.output.receptors) || !isempty(s.output.receptors_file)
    # (receptors work with any grid kernel: the engine refreshes the per-parcel
    #  bandwidths for the gridless KDE independently when the grid isn't adaptive)
    !check_paths || isempty(s.output.receptors_file) || isfile(s.output.receptors_file) ||
        push!(errs, "output.receptors_file not found: $(s.output.receptors_file)")
    length(s.output.receptors_layer) == 2 && s.output.receptors_layer[1] < s.output.receptors_layer[2] ||
        push!(errs, "output.receptors layer must be [z_bottom, z_top] with z_bottom < z_top")

    # observed-track nudging / wind-bias drift
    if !isempty(s.physics.nudge_file)
        !check_paths || isfile(s.physics.nudge_file) ||
            push!(errs, "physics.nudge.file not found: $(s.physics.nudge_file)")
        s.physics.nudge_tau_s > 0 || push!(errs, "physics.nudge.tau must be > 0")
        (length(s.physics.nudge_z) == 2 && s.physics.nudge_z[1] < s.physics.nudge_z[2]) ||
            push!(errs, "physics.nudge.z must be [z_bottom, z_top] with z_bottom < z_top")
    end
    (length(s.physics.windbias_z) == 2 && s.physics.windbias_z[1] < s.physics.windbias_z[2]) ||
        push!(errs, "physics.wind_bias z must be [z_bottom, z_top] with z_bottom < z_top")

    # This distribution has no reactive-chemistry module (review 2026-09-23): refuse the
    # chemistry task, chemistry.enabled and chemistry.emissions here, with the same sentence
    # build_run uses, instead of failing later with an UndefVarError.
    _NOCHEM = "the reactive-chemistry module is not included in this EASLEE distribution (transport + aerosol only)"
    s.task == "chemistry" && push!(errs, "task = \"chemistry\": $_NOCHEM")
    s.chemistry.enabled && push!(errs, "[chemistry] enabled = true: $_NOCHEM — remove the [chemistry] block or set enabled = false")
    isempty(s.chemistry.emissions_dir) || push!(errs, "[chemistry.emissions]: $_NOCHEM — the emission inventories are read by that module")
    s.task == "chemistry" && !s.chemistry.enabled &&
        push!(errs, "task=\"chemistry\" needs a [chemistry] block (enabled = true) — it is the domain-fill chemistry comparison run")

    if s.chemistry.enabled
        c = s.chemistry
        c.mechanism in ("reduced","reducedqss","bb12","cb05","strat","strat_brox") ||
            push!(errs, "chemistry.mechanism must be reduced|reducedqss|bb12|cb05|strat")
        c.init in ("background","cams","sources") ||
            push!(errs, "chemistry.init must be background|cams|sources")
        c.ch4_oxidation && !c.h2o_prognostic &&
            push!(errs, "chemistry.ch4_oxidation requires chemistry.h2o_prognostic = true " *
                        "(the produced H2O must accumulate on the parcel, not be " *
                        "re-prescribed from met every step)")
        c.ch4_oxidation && !startswith(c.mechanism, "strat") &&
            push!(errs, "chemistry.ch4_oxidation needs a stratospheric mechanism (strat|strat_brox)")
        c.sulfate_sad && !startswith(c.mechanism, "strat") &&
            push!(errs, "chemistry.sulfate_sad needs a stratospheric mechanism (strat|strat_brox)")
        c.psc_dehydration && !c.h2o_prognostic &&
            push!(errs, "chemistry.psc_dehydration requires chemistry.h2o_prognostic = true")
        c.psc_dehydration && !startswith(c.mechanism, "strat") &&
            push!(errs, "chemistry.psc_dehydration needs a stratospheric mechanism (strat|strat_brox)")
        c.psc_denitrification && !startswith(c.mechanism, "strat") &&
            push!(errs, "chemistry.psc_denitrification needs a stratospheric mechanism (strat|strat_brox)")
        # [chemistry.init] declarative block (CHEM_INIT.md)
        if c.init == "sources" || !isempty(c.init_species) || !isempty(c.init_sources)
            srcnames = Set(src.name for src in c.init_sources)
            length(srcnames) == length(c.init_sources) ||
                push!(errs, "chemistry.init: duplicate source names")
            for src in c.init_sources
                src.kind in ("grib","netcdf","constant") ||
                    push!(errs, "chemistry.init source \"$(src.name)\": kind must be grib|netcdf|constant")
                src.units in ("mmr","vmr","ppb") ||
                    push!(errs, "chemistry.init source \"$(src.name)\": units must be mmr|vmr|ppb")
                src.kind != "constant" && isempty(src.dir) &&
                    push!(errs, "chemistry.init source \"$(src.name)\": dir required")
            end
            isempty(c.init_default) || c.init_default in srcnames ||
                push!(errs, "chemistry.init.default \"$(c.init_default)\" is not a declared source")
            derivs = ("correlation","steady_state","family","constant","profile")
            names = Set{String}()
            for sp in c.init_species
                sp.name in names && push!(errs, "chemistry.init: species \"$(sp.name)\" assigned twice")
                push!(names, sp.name)
                sp.from in srcnames || sp.from in derivs ||
                    push!(errs, "chemistry.init species \"$(sp.name)\": from must be a source or $(join(derivs, "|"))")
                if sp.from == "correlation"
                    isempty(sp.predictor) &&
                        push!(errs, "chemistry.init species \"$(sp.name)\": correlation needs predictor")
                    isempty(sp.curve) ? push!(errs, "chemistry.init species \"$(sp.name)\": correlation needs curve") :
                        (!check_paths || isfile(sp.curve) || push!(errs, "chemistry.init species \"$(sp.name)\": curve not found: $(sp.curve)"))
                end
                sp.from == "family" && isempty(sp.total) &&
                    push!(errs, "chemistry.init species \"$(sp.name)\": family needs total")
                sp.from == "profile" && check_paths && !isfile(sp.profile) &&
                    push!(errs, "chemistry.init species \"$(sp.name)\": profile not found: $(sp.profile)")
                sp.regime in ("any","day") ||
                    push!(errs, "chemistry.init species \"$(sp.name)\": regime must be any|day")
            end
        end
        c.photolysis in ("analytic","lut") ||
            push!(errs, "chemistry.photolysis must be analytic|lut")
        c.aerosols && !(c.init == "cams" || any(src -> src.kind == "grib", c.init_sources)) &&
            push!(errs, "chemistry.aerosols needs a CAMS grib init source (the aermr fields come from the gribs)")
        if c.photolysis == "lut" && !isempty(c.lut_path) && check_paths && !isfile(c.lut_path)
            push!(errs, "chemistry.lut_path not found: $(c.lut_path)")
        end
        c.timestep_s > 0 || push!(errs, "chemistry.timestep must be > 0")
        _ismult(c.timestep_s, s.timestep_s) ||
            push!(errs, "chemistry.timestep must be a multiple of run.timestep")
        (c.photolysis_timestep_s == 0 || _ismult(c.photolysis_timestep_s, s.timestep_s)) ||
            push!(errs, "chemistry.photolysis_timestep must be a multiple of run.timestep")
        if c.init == "cams"
            isempty(c.cams_dir) &&
                push!(errs, "chemistry.init=\"cams\" requires chemistry.cams_dir")
            !check_paths || isdir(c.cams_dir) || isfile(c.cams_dir) ||
                push!(errs, "chemistry.cams_dir not found: $(c.cams_dir)")
        end
        c.mixing in ("off", "pairwise") ||
            push!(errs, "chemistry.mixing.scheme must be off|pairwise")
        c.mixing == "pairwise" && c.mixing_tau_s <= 0 &&
            push!(errs, "chemistry.mixing.tau must be > 0")
        (c.mixing_cells >= 1 && c.mixing_levels >= 1) ||
            push!(errs, "chemistry.mixing.cells/levels must be ≥ 1")
        if !isempty(c.emissions_dir)
            !check_paths || isdir(c.emissions_dir) ||
                push!(errs, "chemistry.emissions.dir not found: $(c.emissions_dir)")
            any(r -> r.kind == "domainfill", s.release) ||
                push!(errs, "chemistry.emissions requires a [[release]] kind=\"domainfill\" (air-mass parcels)")
            c.emissions_res > 0 || push!(errs, "chemistry.emissions.res must be > 0")
            (0.0 <= c.emissions_nox_no <= 1.0) ||
                push!(errs, "chemistry.emissions.nox_no must be in [0, 1]")
            c.emissions_diurnal in ("none", "gaussian", "mu2011") ||
                push!(errs, "chemistry.emissions.diurnal must be none|gaussian|mu2011")
            c.emissions_diurnal_sigma > 0 ||
                push!(errs, "chemistry.emissions.diurnal_sigma must be > 0")
        end
        # [output.reference]: the named init source must be declared
        if !isempty(s.output.reference)
            srcnames = Set(src.name for src in c.init_sources)
            s.output.reference in srcnames ||
                push!(errs, "output.reference \"$(s.output.reference)\" is not a declared [chemistry.init] source")
        end
    elseif !isempty(s.output.reference)
        push!(errs, "output.reference needs chemistry.enabled = true (it compares the reacting run to its init source)")
    end

    # ── time-step hierarchy: every cadence is a multiple of run.timestep, and
    #    output.sample | output.average | output.interval ──────────────────────
    o = s.output
    _ismult(o.sample_s,  s.timestep_s) || push!(errs, "output.sample must be a multiple of run.timestep")
    _ismult(o.average_s, s.timestep_s) || push!(errs, "output.average must be a multiple of run.timestep")
    _ismult(o.interval_s,s.timestep_s) || push!(errs, "output.interval must be a multiple of run.timestep")
    o.sample_s  <= o.average_s  || push!(errs, "output.sample must be ≤ output.average")
    o.average_s <= o.interval_s || push!(errs, "output.average must be ≤ output.interval")
    _ismult(o.average_s, o.sample_s)   || push!(errs, "output.average must be a multiple of output.sample")
    (s.physics.convection.timestep_s == 0 || _ismult(s.physics.convection.timestep_s, s.timestep_s)) ||
        push!(errs, "physics.convection.timestep must be a multiple of run.timestep")
    # review 2026-09-04 A7 item 6: parsed + validated, consumed nowhere in easlee/
    s.physics.convection.timestep_s > 0 && s.physics.convection.timestep_s != s.timestep_s &&
        @warn "physics.convection.timestep = $(s.physics.convection.timestep_s) s is not implemented: convection runs every run.timestep ($(s.timestep_s) s); the key is accepted but has no effect"
    s.physics.ifine >= 1 || push!(errs, "physics.ifine must be ≥ 1")

    isempty(errs) || error("Invalid run configuration:\n  - " * join(errs, "\n  - "))
    return s
end

end # module RunSpecs
