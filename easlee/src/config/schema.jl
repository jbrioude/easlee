
# ── Schema — field metadata driving the GUI form, CLI help, and validation ────
# The RunSpec structs own the types and defaults; this table adds the human
# layer (enum options, help, units, grouping, visibility) so the GUI can render
# a form generically and the CLI can produce --help.  ONE table → form + CLI +
# docs, so they never drift.

module Schema

using ..RunSpecs

export Field, schema_fields, schema_json, schema_extra_paths, schema_known_paths

struct Field
    path::String          # dotted path; "[]" marks an array-of-tables element
    type::String          # "string" | "number" | "int" | "bool" | "datetime" | "enum" | "list"
    default::Any
    group::String
    help::String
    unit::String
    options::Vector{String}     # enum values
    visible_if::String          # e.g. "species[].kind==aerosol"
end
Field(path, type, default, group, help; unit="", options=String[], visible_if="") =
    Field(path, type, default, group, help, unit, options, visible_if)

# The full field list.  Order = display order within each group.
function schema_fields()::Vector{Field}
    F = Field[]
    push!(F, Field("name", "string", "easlee_run", "Run", "Run name (output file prefix)"))
    push!(F, Field("task", "enum", "transport", "Run", "What to run: plain transport or DATEM validation (the chemistry task is not part of this distribution)";
                   options=["transport","validation"]))
    push!(F, Field("run.direction", "enum", "forward", "Run", "Integration direction"; options=["forward","backward"]))
    push!(F, Field("run.start", "datetime", "", "Run", "Simulation start (UTC)"))
    push!(F, Field("run.stop",  "datetime", "", "Run", "Simulation stop (UTC)"))
    push!(F, Field("run.timestep", "string", "15min", "Run", "Synchronisation time step"; unit="duration"))
    push!(F, Field("run.wind", "enum", "3D", "Run", "Wind dimensionality: 2D/1D freeze parcel motion along the dropped axes (the frozen index stays at its release value) — reduced-dimension transport studies";
                   options=["3D","2D-lonlat","2D-latalt","2D-lonalt","1D-vertical"]))
    push!(F, Field("run.benchmark", "bool", false, "Run", "Per-group @time-style timing report at run end"))
    push!(F, Field("run.seed", "int", 0, "Run", "Random seed (Random.seed! at run start; 0 = unseeded). Bit-for-bit reproducibility also needs the same thread count and the same advance backend"))
    push!(F, Field("run.checkpoint", "string", "", "Run", "Warm-start checkpoint: \"end\" = save the full parcel/chemistry state at run end; a duration (e.g. \"6h\") also saves periodically — resume later with run.restart instead of re-simulating"))
    push!(F, Field("run.restart", "string", "", "Run", "Path to a checkpoint file: resume the run from the saved step (same config, later stop; output files get an _r suffix)"))
    push!(F, Field("run.parcel_dump", "string", "", "Run", "Per-parcel diagnostic snapshots: a duration (e.g. \"3h\") writes parcel + convert-gate state and flux accumulators to <output.path>/snapshots/ at that cadence; \"\" = off. Diagnostic only, never read back by the model"))
    # review 2026-09-04 A7 item 1/6: keys the reader/launcher consume but the schema lacked
    push!(F, Field("run.threads", "int", 8, "Run", "Worker threads passed as `julia -t` by the launcher (scripts/easlee.jl); the engine itself does not read it"))
    push!(F, Field("strict_keys", "bool", true, "Run", "Unknown-key policy: true = a TOML key that is not in this schema is an ERROR (it would be silently ignored otherwise); false = warn and load"))

    push!(F, Field("meteorology[].model", "enum", "ERA5", "Meteorology", "Driving model";
                   options=["ERA5","MERRA2","GFS","AROME_HDF","MESONH","WRF","AROME"]))
    push!(F, Field("meteorology[].path", "string", "", "Meteorology", "Directory of met files"))
    push!(F, Field("meteorology[].variables", "string", "", "Meteorology", "Variables namelist (listvariables/<name>.vars)"))
    push!(F, Field("meteorology[].available", "string", "auto", "Meteorology", "AVAILABLE file path; 'auto' = use <path>/AVAILABLE when present and up to date (written by download/*.py or by an earlier run), else scan the path and write it; 'scan' = always rescan"))
    push!(F, Field("meteorology[].crop", "string", "auto", "Meteorology", "Domain crop: 'auto' (crop to the domainfill box + 2 cells; exact — the cull keeps parcels inside), 'off', or [lon0, lon1, lat0, lat1] (any run; parcels are terminated at the crop edge)"))
    push!(F, Field("meteorology[].max_gap", "string", "6h", "Meteorology", "Refuse to run when the met archive has a coverage gap larger than this inside the run window (a missing file, not an intended overrun); 0 = disable the guard"; unit="duration"))   # review 2026-09-04 A7 item 6

    push!(F, Field("species[].id", "string", "spec1", "Species", "Identifier (referenced by releases)"))
    push!(F, Field("species[].kind", "enum", "gas", "Species", "Material kind"; options=["gas","aerosol"]))
    push!(F, Field("species[].specnum", "int", 0, "Species", "FLEXPART species number (optional metadata; 0 = none)"))
    push!(F, Field("species[].molar_mass", "number", 0.0, "Species", "Molar mass [g/mol]; 0 = from the species DB by id (else unknown). Needed for the automatic [species.convert] mass_ratio = M_to/M_from"; unit="g/mol"))   # review 2026-09-04 A7 item 3
    push!(F, Field("species[].aerosol_type", "enum", "dust", "Species", "Aerosol type (κ + refractive index defaults)";
                   options=["dust","sulfate","bc","oc","seasalt","ash","water"], visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].size.distribution", "enum", "lognormal", "Species", "Size distribution";
                   options=["lognormal","bins"], visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].size.dg", "string", "1um", "Species", "Lognormal median diameter"; unit="length", visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].size.sigma_g", "number", 2.0, "Species", "Geometric standard deviation"; visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].size.nbins", "int", 8, "Species", "Number of size bins"; visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].size.edges", "list", Float64[], "Species", "Explicit bin edge diameters (distribution=bins)"; unit="length", visible_if="species[].size.distribution==bins"))
    push!(F, Field("species[].density", "number", 2000.0, "Species", "Particle density"; unit="kg/m3", visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].kappa", "number", -1.0, "Species", "κ-Köhler hygroscopicity (<0 = default from type)"; visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].refractive_index", "string", "auto", "Species", "Refractive index n−ki (\"auto\" = from type, or e.g. 1.53-0.003i)"; visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].settling", "bool", true, "Species", "Gravitational settling (aerosol size bins). TOML also accepts a [species.settling] TABLE {r_eff, rho, enabled} for the condensed-species fall of a gas-kind tracer (ice/ash shadow population) — see schema_extra_paths"; visible_if="species[].kind==aerosol"))
    push!(F, Field("species[].dry_deposition", "bool", false, "Species", "Dry deposition"))
    push!(F, Field("species[].wet_deposition", "bool", false, "Species", "Wet deposition"))
    push!(F, Field("species[].decay_halflife", "string", "", "Species", "Radioactive half-life (blank = stable)"; unit="duration"))
    push!(F, Field("species[].ice.enabled", "bool", false, "Species", "Transport-phase ice microphysics for a tagged water species: above rhi_crit × saturation over ICE (Murphy–Koop) the excess vapour condenses and sediments out of the stratospheric budget — the first-48-h loss the eruption-column washout cannot see"))
    push!(F, Field("species[].ice.rhi_crit", "number", 1.5, "Species", "Freeze only above this RHi. 1.5 = HOMOGENEOUS nucleation (Koop) — correct only for a plume with no ice nuclei. When efficient IN are present (volcanic ash, mineral dust) deposition nucleation starts near RHi 1.05–1.2 and 1.5 is unphysical: use ~1.1"; visible_if="species[].ice.enabled==true"))
    push!(F, Field("species[].ice.r_eff", "string", "20um", "Species", "Effective ice radius → Stokes–Cunningham fall speed at the parcel's T,p (~6 cm/s at 25 km for 20 µm)"; unit="length", visible_if="species[].ice.enabled==true"))
    push!(F, Field("species[].ice.h_layer", "string", "500m", "Species", "Settling escape scale: condensate leaves the layer at rate v_t/h_layer"; unit="length", visible_if="species[].ice.enabled==true"))
    push!(F, Field("species[].ice.q0", "number", 0.0, "Species", "Prescribed release mixing ratio [kg/kg] defining the parcel air mass (q = mass/airm); 0 = from the plume ensemble's detrainment water fraction (plume_rise mode=parcels)"; visible_if="species[].ice.enabled==true"))
    push!(F, Field("species[].ice_scav.of", "string", "", "Species", "Ice co-scavenging: this species rides the co-located freezing water — name an ice-enabled species; the parcel loses mass · uptake · f_shed(cell) wherever that species sheds ice"))
    push!(F, Field("species[].ice_scav.uptake", "number", 0.0, "Species", "Retained-in-ice fraction (0 = off; HNO3 ≈ 1, SO2/HCl smaller and T-dependent)"; visible_if="species[].ice_scav.of!="))
    push!(F, Field("species[].convert.to", "string", "", "Species", "In-flight species conversion: each parcel switches to this species with probability 1−exp(−dt/τ) per step (e.g. SO2 → a sulfate aerosol species; the daughter gets aerosol physics/AOD from the moment of conversion)"))
    push!(F, Field("species[].convert.tau", "string", "", "Species", "Conversion e-folding time (e.g. \"25d\" for stratospheric SO2 → sulfate)"; unit="duration", visible_if="species[].convert.to!="))
    push!(F, Field("species[].convert.mass_ratio", "number", 0.0, "Species", "Mass rescaling on conversion (0 = auto: M_to/M_from from the species molar masses; SO2→sulfate ≈ 1.5)"; visible_if="species[].convert.to!="))
    push!(F, Field("species[].convert.mode", "enum", "", "Species", "Thermodynamic gate: freeze (only where q > rhi_crit·q_si) / sublimate (only where q < q_si). The pair makes the phase change reversible"; options=["", "freeze", "sublimate"], visible_if="species[].convert.to!="))
    push!(F, Field("species[].convert.rhi_crit", "number", 1.5, "Species", "Freezing threshold on RHi (mode = freeze)"; visible_if="species[].convert.to!="))
    push!(F, Field("species[].convert.lamina_sigma", "number", 0.0, "Species", "Sub-grid lamina width: relative lognormal spread of q about the parcel mean; 0 = step gate (bit-identical legacy)"; visible_if="species[].convert.to!="))
    push!(F, Field("species[].convert.sectional", "bool", false, "Species", "SECTIONAL ice: per-parcel crystal radius shrinking on the exact r² law, bin-to-bin demotion when it crosses the bin edge, and size-resolved settling — the only mechanism that gives SIZE SELECTIVITY (the lamina and T′ scale all bins alike). Set on the sublimating ICE species"; visible_if="species[].convert.to!="))

    push!(F, Field("release[].species", "string", "spec1", "Release", "Species id released"))
    push!(F, Field("release[].kind", "enum", "release", "Release", "release = emit from the geometry; domainfill = fill the box ∝ air mass + edge inflow; inventory = gridded emission inventory (from_file); volcanic = time-resolved eruption column at the vent; seaspray = met-driven OSSA sea-spray source over the ocean cells of the box (Ovadnevaite et al. 2014; needs mass_per_parcel)";
                   options=["release","domainfill","inventory","volcanic","seaspray"]))
    push!(F, Field("release[].seaspray.sst", "number", -1.0, "Release", "Sea-surface temperature [K] for the wave Reynolds number; -1 = met 2t per cell"; visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].seaspray.salinity", "number", 35.0, "Release", "Salinity [‰] (seawater viscosity → Re_Hw)"; visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].seaspray.hs", "number", -1.0, "Release", "Significant wave height [m]; -1 = fully-developed Pierson–Moskowitz H_s = 0.0246·U10²"; visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].seaspray.cd", "number", -1.0, "Release", "Wave drag coefficient; -1 = bulk (0.8 + 0.065·U10)×10⁻³"; visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].seaspray.gust_sigma", "number", 0.0, "Release", "Sub-grid/gust wind σ [m/s]: the flux integrates over N(U10, σ) — the transient contribution beyond the mean fields (largest below ~10-15 m/s mean wind); 0 = mean-field only"; visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].seaspray.rho_dry", "number", 2017.0, "Release", "Dry sea-salt density for the mass flux"; unit="kg/m3", visible_if="release[].kind==seaspray"))
    push!(F, Field("release[].chem_init.so2_ppb", "number", 0.0, "Release", "Release-scoped chemistry IC: plume SO2 [ppb] on THIS release's parcels only (→ per-parcel sulfate SAD); every other parcel keeps the [chemistry.init] background. Needs [chemistry].sulfate_sad"; unit="ppb"))
    push!(F, Field("release[].chem_init.so4_ugm3", "number", 0.0, "Release", "Release-scoped chemistry IC: plume sulfate loading [µg/m³] on this release's parcels (→ SAD). Needs [chemistry].sulfate_sad"; unit="ug/m3"))
    push!(F, Field("release[].chem_init.h2o_ppm", "number", 0.0, "Release", "Release-scoped chemistry IC: plume water vapour [ppmv] on this release's parcels; pins the prognostic water so met stops re-prescribing it. Needs [chemistry].h2o_prognostic"; unit="ppmv"))
    push!(F, Field("release[].volcanic.profile", "enum", "suzuki", "Release", "Vertical mass profile of each eruption phase";
                   options=["suzuki","umbrella","uniform"], visible_if="release[].kind==volcanic"))
    push!(F, Field("release[].volcanic.suzuki_a", "number", 4.0, "Release", "Suzuki shape A (larger → mass closer to column top)"; visible_if="release[].volcanic.profile==suzuki"))
    push!(F, Field("release[].volcanic.suzuki_lambda", "number", 1.0, "Release", "Suzuki shape λ (peak sharpness)"; visible_if="release[].volcanic.profile==suzuki"))
    push!(F, Field("release[].volcanic.from_file", "string", "", "Release", "CSV injection table: t0, t1, z_bottom, z_top, frac (alternative to inline phases)"; visible_if="release[].kind==volcanic"))
    push!(F, Field("release[].plume_rise.scheme", "enum", "none", "Release", "Injection height from the met column at release time: fire FRP (sofiev/freitas) or a 1-D moist eruption column (volcanic)";
                   options=["none","sofiev","freitas","volcanic"]))
    push!(F, Field("release[].plume_rise.frp", "number", 0.0, "Release", "Fire radiative power (fire schemes: sofiev/freitas)"; unit="W", visible_if="release[].plume_rise.scheme!=none"))
    push!(F, Field("release[].plume_rise.fire_area", "number", 0.0, "Release", "Burning area (freitas: initial plume radius R₀ = √(A/π))"; unit="m²", visible_if="release[].plume_rise.scheme!=none"))
    push!(F, Field("release[].plume_rise.profile", "enum", "uniform", "Release", "Parcel fill over the computed injection layer";
                   options=["uniform","gaussian"], visible_if="release[].plume_rise.scheme!=none"))
    push!(F, Field("release[].plume_rise.bottom", "string", "pbl", "Release", "Injection layer bottom: pbl | ground | fraction:<f> (fraction of the plume top)"; visible_if="release[].plume_rise.scheme!=none"))
    push!(F, Field("release[].plume_rise.mer", "number", 0.0, "Release", "Volcanic: mass eruption rate (total mixture)"; unit="kg/s", visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.temperature", "number", 800.0, "Release", "Volcanic: mixture temperature after magma–water interaction"; unit="K", visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.water_frac", "number", 0.05, "Release", "Volcanic: vapour mass fraction at the vent (phreatomagmatic seawater)"; visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.mode", "enum", "layer", "Release", "layer: inject over the computed [NBH, top]; parcels: each trajectory rises through the Lagrangian plume phase and enters at its emergent detrainment height (mass scaled by washout survival — release mass is then the AT-VENT mass)";
                   options=["layer","parcels"], visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.washout", "bool", true, "Release", "Parcels mode: in-plume + umbrella washout (SO2/HCl Henry scavenging, water fallout)"; visible_if="release[].plume_rise.mode==parcels"))
    push!(F, Field("release[].plume_rise.entrainment", "enum", "meanfield", "Release", "Parcels mode: meanfield (v1, each parcel dilutes toward the unperturbed ambient) or pairwise (v2, edge engulfment + parcel-parcel exchange → hot core / cool edge)";
                   options=["meanfield","pairwise"], visible_if="release[].plume_rise.mode==parcels"))
    push!(F, Field("release[].plume_rise.umbrella_spread", "bool", false, "Release", "Umbrella gravity-current spreading during transport: R(t)=(3λNQ_V t²/2π)^(1/3) fed by the column volume flux — the archived met has no eruption outflow"; visible_if="release[].plume_rise.mode==parcels"))
    push!(F, Field("release[].plume_rise.spread_lambda", "number", 0.2, "Release", "Intrusion constant λ (Costa et al. 2013)"; visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.umbrella_r_source", "string", "10km", "Release", "Source-fed interior radius: inside it parcels ride the collapse-dome outflow u(r) = Q_V/(2πrh̄); emergent parcels start area-uniform across this dome"; unit="length", visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.qc_crit", "number", 0.005, "Release", "Hydrometeor holding capacity of the UPDRAFT [kg/kg]: condensate above this is shed during ascent. Also the legacy umbrella retention — prefer umbrella_r_ice for the umbrella"; visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.vapor_cap", "number", 1.0e30, "Release", "Umbrella vapour cap [× ice saturation at ambient T]: retained vapour above cap·q_si becomes the COARSE-mode budget (fcoarse, printed in the ensemble log line), to be released as a large-crystal ICE_C* species; unset/Inf = off"; visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.umbrella_r_ice", "string", "0m", "Release", "Mass-median crystal radius of the umbrella spectrum. >0 replaces the qc_crit umbrella retention by SEDIMENTATION (physical: the umbrella has w≈0). SETS the mass delivered to the stratosphere — HT: qc_crit=5 g/kg ≡ r≈105 µm"; unit="length", visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.umbrella_r_ice_sigma", "number", 2.0, "Release", "Geometric standard deviation of the lognormal crystal spectrum (a single radius is knife-edge: 50 µm→52 % retained, 200 µm→0 %)"; visible_if="release[].plume_rise.umbrella_r_ice!=0m"))
    push!(F, Field("release[].plume_rise.umbrella_shed_depth", "string", "2km", "Release", "Umbrella thickness the crystals settle through"; unit="length", visible_if="release[].plume_rise.umbrella_r_ice!=0m"))
    push!(F, Field("release[].plume_rise.umbrella_shed_time", "string", "1h", "Release", "Time over which the umbrella sheds"; unit="duration", visible_if="release[].plume_rise.umbrella_r_ice!=0m"))
    push!(F, Field("release[].plume_rise.tephra_r", "string", "0m", "Release", "Mass-median coarse-ash radius for TEPHRA fallout during ascent (mode=parcels): solids shed where v_t beats the updraft gate, unloading the parcel so it detrains lower and colder. 0 = legacy (ash rides to detrainment, forcing a hot high arrival)"; unit="length", visible_if="release[].plume_rise.scheme==volcanic"))
    push!(F, Field("release[].plume_rise.spread_zgain", "number", 0.0, "Release", "Altitude-dependent umbrella amplification: u_r × (1+zgain·ramp(z; z0→z1)). The single-level law under-spreads the upper cloud (HT backward footprints: 79 % of filaments beyond 4° of the core missing, 62-66 % above 26 km). 0 = off"; visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.spread_zgain_z0", "string", "26km", "Release", "Ramp start of the z-dependent amplification"; unit="length", visible_if="release[].plume_rise.spread_zgain!=0"))
    push!(F, Field("release[].plume_rise.spread_zgain_z1", "string", "30km", "Release", "Ramp end (full amplification at/above)"; unit="length", visible_if="release[].plume_rise.spread_zgain!=0"))
    push!(F, Field("release[].plume_rise.spread_azgain", "number", 0.0, "Release", "Azimuth-dependent umbrella amplification: u_r × (1+azgain·cos²θ), θ = azimuth from north — meridional quadrants spread faster, zonal axis untouched (HT inversion wings). Sets both hemispheres; spread_azgain_n / spread_azgain_s override per side. 0 = off"; visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.spread_azgain_n", "number", 0.0, "Release", "Azimuthal amplification for the NORTHERN hemisphere of the vent (overrides spread_azgain on that side)"; visible_if="release[].plume_rise.umbrella_spread==true"))   # review 2026-09-04 A7 item 6
    push!(F, Field("release[].plume_rise.spread_azgain_s", "number", 0.0, "Release", "Azimuthal amplification for the SOUTHERN hemisphere of the vent (overrides spread_azgain on that side)"; visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.spread_azgain_ramp_z0", "string", "0km", "Release", "Optional z-ramp start for the azimuthal amplification (z1 ≤ z0 = no ramp, full at all z)"; unit="length", visible_if="release[].plume_rise.spread_azgain!=0"))
    push!(F, Field("release[].plume_rise.spread_azgain_ramp_z1", "string", "0km", "Release", "z-ramp end (full azimuthal amplification at/above)"; unit="length", visible_if="release[].plume_rise.spread_azgain!=0"))
    push!(F, Field("release[].plume_rise.spread_hmin", "string", "0m", "Release", "PHYSICAL handoff: stop the umbrella once the intrusion has thinned to h = V/(πR²) below this — it is no longer a distinct layer. Preferred over spread_umin (HT: h≈350 m tracks the MLS width on every day). 0 = disabled"; unit="length", visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.spread_umin", "number", 5.0, "Release", "Handoff: stop the radial push once Ṙ < this [m/s]. ENDS the gravity current, so it SETS the final cloud width (HT: 5 m/s fires at 8 h, R=432 km). Keep below the ambient |v| at the neutral level if meridional width matters"; unit="m/s", visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.spread_tmax", "string", "12h", "Release", "Handoff: stop the radial push past this age"; unit="duration", visible_if="release[].plume_rise.umbrella_spread==true"))
    push!(F, Field("release[].plume_rise.umbrella_wind_shield", "bool", false, "Release", "Cancel the resolved-wind displacement of in-cloud parcels by s = min(1, u_r/|u_amb|): the density current expels the ambient flow, so met advection would double-count; hands off smoothly as the current decays"; visible_if="release[].plume_rise.umbrella_spread==true"))
    # release[].plume_rise.vent_mix (bool OR [.. .vent_mix] subtable: magma_t,
    # magma_water, seawater_core, seawater_edge, mix_exponent) is TOML-only:
    # the flat form would silently enable it by collecting subtable defaults.
    push!(F, Field("release[].geometry", "enum", "box", "Release", "Source geometry";
                   options=["box","point","line","polygon","profile"], visible_if="release[].kind!=inventory"))
    push!(F, Field("release[].lon", "list", [0.0,0.0], "Release", "Longitude: box [a,b]; point scalar; line/polygon vertices"; unit="deg"))
    push!(F, Field("release[].lat", "list", [0.0,0.0], "Release", "Latitude: box [a,b]; point scalar; line/polygon vertices"; unit="deg"))
    push!(F, Field("release[].npoints", "int", 100, "Release", "Sampling density for line/polygon"; visible_if="release[].geometry!=box"))
    push!(F, Field("release[].z", "list", [0.0,100.0], "Release", "Height box [z1, z2]"; unit="length"))
    push!(F, Field("release[].zkind", "enum", "agl", "Release", "Vertical reference"; options=["agl","asl","pressure"]))
    push!(F, Field("release[].time", "list", ["",""], "Release", "Release window [start, stop]"; unit="datetime"))
    push!(F, Field("release[].particles", "int", 10000, "Release", "Total number of particles"))
    push!(F, Field("release[].mass", "string", "1kg", "Release", "Total released mass"; unit="mass"))
    push!(F, Field("release[].mass_per_parcel", "string", "", "Release", "Mass quantum: the parcel count follows the released mass, n = ⌈mass/quantum⌉ (overrides particles; empty = use particles)"; unit="mass"))
    push!(F, Field("release[].from_file", "string", "", "Release", "Locations/inventory file (overrides the lon/lat box)"))
    push!(F, Field("release[].file_format", "enum", "latlonlev", "Release", "Columns in the locations file";
                   options=["latlonlev","lonlatlev","inventory"], visible_if="release[].from_file!="))
    push!(F, Field("release[].comment", "string", "", "Release", "Free-text label for this release"))

    push!(F, Field("physics.turbulence.scheme", "enum", "hanna", "Physics", "Turbulence scheme"; options=["none","hanna","thomson"]))
    push!(F, Field("physics.turbulence.onthefly", "bool", false, "Physics", "Hanna option b: compute σ on the fly (vs profile lookup)"; visible_if="physics.turbulence.scheme==hanna"))
    push!(F, Field("physics.turbulence.sigw_ft", "number", 0.1, "Physics", "Free-troposphere σ_w"; unit="m/s", visible_if="physics.turbulence.scheme!=none"))
    push!(F, Field("physics.turbulence.kh_ft", "number", 0.0, "Physics", "FLEXPART free-atmosphere horizontal diffusivity D_h (d_trop = 50): random walk √(2KΔt) above the PBL; 0 = off. Trajectory route — keep [physics.mixing] OFF when set (one K moves parcels OR exchanges mass, never both)"; unit="m²/s", visible_if="physics.turbulence.scheme!=none"))
    push!(F, Field("physics.turbulence.mesoscale", "bool", false, "Physics", "Mesoscale velocity fluctuations (FLEXPART usig/vsig/wsig): per-parcel OU velocity, amplitude fturbmeso × σ of the surrounding grid-point winds, memory = met archive interval — sub-grid wind uncertainty in all layers"))
    push!(F, Field("physics.turbulence.horizontal_tke", "bool", false, "Physics", "TKE-driven HORIZONTAL turbulence: OU velocity with memory from the same archived e as the vertical step (sigma_u^2 = sigma_v^2 = 2*e*(1-vmix_rw)/2, tau_Lh = K_h/sigma_u^2). Native-e drivers ONLY (AROME native e, MesoNH); a Hanna-fill driver archives w-variance only and build_run refuses it. Trajectory route - keep [physics.mixing] OFF"; visible_if="physics.turbulence.scheme!=none"))
    push!(F, Field("physics.turbulence.horizontal_hanna", "bool", false, "Physics", "Hanna (1982) HORIZONTAL Langevin in the PBL: sigma_u, sigma_v, T_Lu, T_Lv from the same u*, w*, h, L as the vertical Hanna fill — FLEXPART's in-PBL horizontal turbulence (grid components, no wind-frame rotation). Hanna-fill drivers ONLY (ERA5/GFS/MERRA-2/AROME-HDF); a native-e nest uses horizontal_tke. Fades to zero across the entrainment zone; kh_ft remains the free-troposphere term. Trajectory route - keep [physics.mixing] OFF"; visible_if="physics.turbulence.scheme!=none"))
    push!(F, Field("physics.turbulence.fturbmeso", "number", 0.16, "Physics", "Mesoscale amplitude factor (FLEXPART default 0.16; empirical — consider lowering for fine/hourly met)"; visible_if="physics.turbulence.mesoscale==true"))
    push!(F, Field("physics.polar_stereo_lat", "number", 75.0, "Physics", "Polar-stereographic transport (FLEXPART switchnorth/switchsouth): poleward of ±this latitude parcel positions advance in a pole-centred stereographic plane — no 1/cos(lat) blow-up, exact pole crossing. Auto-enabled per pole when the (global rectilinear) met grid reaches it, e.g. ERA5 −90…90. 90 = off"))
    push!(F, Field("physics.convection.scheme", "enum", "off", "Physics", "Convection (massflux needs CLD fields; cape is diagnostic)"; options=["off","massflux","cape"]))
    push!(F, Field("physics.convection.timestep", "string", "", "Physics", "Convection step (empty → = run.timestep)"; unit="duration", visible_if="physics.convection.scheme!=off"))
    push!(F, Field("physics.convection.cape_min", "number", 100.0, "Physics", "CAPE event threshold (cape scheme)"; unit="J/kg", visible_if="physics.convection.scheme==cape"))
    push!(F, Field("physics.convection.cin_max", "number", 50.0, "Physics", "Skip columns more strongly capped than this (cape scheme)"; unit="J/kg", visible_if="physics.convection.scheme==cape"))
    push!(F, Field("physics.convection.tau_mix", "string", "", "Physics", "Convective mixing timescale (cape scheme; empty/0 = extreme limit)"; unit="duration", visible_if="physics.convection.scheme==cape"))
    push!(F, Field("physics.vertical.mode", "enum", "kinematic", "Physics", "Vertical transport"; options=["kinematic","hybrid","isentropic"]))
    push!(F, Field("physics.vertical.sigma_r", "number", 0.3, "Physics", "Blend bottom σ=p/ps (kinematic below)"; visible_if="physics.vertical.mode==hybrid"))
    push!(F, Field("physics.vertical.sigma_t", "number", 0.1, "Physics", "Blend top σ=p/ps (isentropic above)"; visible_if="physics.vertical.mode==hybrid"))
    push!(F, Field("physics.vertical.diabatic", "bool", true, "Physics", "Diabatic θ̇ drift"; visible_if="physics.vertical.mode==hybrid"))
    push!(F, Field("physics.vertical.pt_winds", "string", "auto", "Physics", "Isentropic θ̇ winds: \"auto\" | \"off\" | met directory"; visible_if="physics.vertical.mode!=kinematic"))
    push!(F, Field("physics.vertical.pt_pattern", "string", "IFS_PT_MET", "Physics", "Filename pattern for θ̇ met files"; visible_if="physics.vertical.mode!=kinematic"))
    push!(F, Field("physics.surface_layer_wind", "bool", true, "Physics", "Monin-Obukhov surface-layer wind: below the lowest model level the archived wind is constant-extended to the ground (a parcel at 2 m carries the ~20 m wind); on = rescale onto the similarity profile anchored at that level, using the run own u* and Obukhov length"))
    push!(F, Field("physics.deposition.dry", "bool", false, "Physics", "Dry deposition (global)"))
    push!(F, Field("physics.deposition.wet", "bool", false, "Physics", "Wet deposition (global)"))
    push!(F, Field("physics.mixing.scheme", "enum", "off", "Physics", "Inter-trajectory mixing during transport: off = each parcel keeps its tracer mass (positions still disperse via turbulence); pairwise = Fickian tracer exchange between co-located same-species parcels + coverage-based entrainment of ambient air at the cloud edge (interior parcels, fully surrounded, see none)"; options=["off","pairwise"]))
    push!(F, Field("physics.mixing.K_h", "number", 1000.0, "Physics", "Horizontal diffusivity of the exchange"; unit="m²/s", visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.sigma0_h", "string", "", "Physics", "Horizontal kernel half-width (empty/0 = met grid spacing)"; unit="length", visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.sigma_v", "string", "", "Physics", "Vertical pair half-width (empty/0 = turbulence Lagrangian length yl — ~10 m in the stratosphere)"; unit="length", visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.background", "bool", true, "Physics", "Coverage-based edge entrainment: parcels not fully surrounded relax toward zero-tracer ambient air at rate (1−cov)·K_h/σ₀²"; visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.nsig", "number", 2.0, "Physics", "Pair-search radius in σ₀ units (cell-list cutoff)"; visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.max_neighbors", "int", 30, "Physics", "Cap on exchange partners per parcel per pass"; visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.mixing.every", "int", 1, "Physics", "Apply every N transport steps"; visible_if="physics.mixing.scheme==pairwise"))
    push!(F, Field("physics.ctl", "number", 5.0, "Physics", "Lagrangian time-scale factor (sub-steps): must be a positive INTEGER value (the sub-step count dt/ctl is exact only then)"))
    push!(F, Field("physics.ifine", "int", 4, "Physics", "Turbulence sub-steps per CTL interval"))
    # review 2026-09-04 A7 item 6: read by toml_io.jl, were missing here
    push!(F, Field("physics.prognostic_airm", "bool", false, "Physics", "Refresh the air mass each parcel represents from the LOCAL parcel spacing (kNN) instead of freezing it at release, so q = mass/airm dilutes as the cloud disperses (ice closure + plume_radiation)"))
    push!(F, Field("physics.prognostic_airm_every", "int", 12, "Physics", "Transport steps between air-mass refreshes"; visible_if="physics.prognostic_airm==true"))

    push!(F, Field("physics.plume_radiation.species", "list", String[], "Physics", "In-plume LW radiative cooling from excess water: species ids it applies to (needs the [species.ice]-style per-parcel air mass); empty = off"))
    push!(F, Field("physics.plume_radiation.c_rad", "number", 0.002, "Physics", "Cooling coefficient [K/day per ppmv excess] — linearised cooling-to-space; CALIBRATE against published heating rates / sonde anomalies"; visible_if="physics.plume_radiation.species!="))
    push!(F, Field("physics.plume_radiation.partners", "list", String[], "Physics", "Species that FEEL the co-located water cooling without producing it (aerosol riding the plume)"))
    push!(F, Field("physics.plume_radiation.q_bg", "number", 5.0, "Physics", "Background water [ppmv] subtracted before cooling"; visible_if="physics.plume_radiation.species!="))
    push!(F, Field("physics.plume_radiation.lw_kappa", "number", 0.0, "Physics", "LW opacity [1/(ppmv·m)] for the thickness-aware saturation of the cooling; 0 = legacy linear law + hard cap (physics/plume_radiation.jl)"; visible_if="physics.plume_radiation.species!="))   # review 2026-09-04 A7 item 6
    push!(F, Field("physics.plume_radiation.h_ref", "string", "1000m", "Physics", "Fallback layer thickness for the lw_kappa saturation when no local thickness is available"; unit="length", visible_if="physics.plume_radiation.lw_kappa!=0"))
    push!(F, Field("physics.plume_radiation.feedback", "bool", true, "Physics", "Buoyancy-equilibrium descent: the cooled parcel sinks by dz = dθ/(∂θ/∂z) — the radiative descent of the water layer"; visible_if="physics.plume_radiation.species!="))
    push!(F, Field("physics.plume_radiation.tprime_tau", "string", "", "Physics", "Relaxation time of the prognostic cold anomaly T′ (e.g. \"1d\"): cooling charges T′, T′ relaxes into the descent, and the freeze/sublimate gates see e_si(T_ambient+T′). Empty/0 = legacy instant buoyancy equilibrium (no anomaly)"; visible_if="physics.plume_radiation.species!="))
    push!(F, Field("tracers.enabled", "bool", false, "Tracers", "Passive transport-history air-mass tags: each parcel carries diagnostics of how it got where it is (complements chemistry + the adjoint). Currently: stratospheric residence (STE)"))
    push!(F, Field("tracers.strat_theta", "number", 380.0, "Tracers", "Potential-temperature threshold [K] for 'stratospheric' air (θ≈380 K = tropical overworld/tropopause) — sets the strat_frac / age_strat tags"; visible_if="tracers.enabled==true"))
    push!(F, Field("tracers.tropo_origin", "bool", false, "Tracers", "Plume-origin air-fraction tag: 1.0 on volcanic-release parcels, diluted toward 0 by transport mixing → gridded f_tropo (needs [physics.mixing])"; visible_if="tracers.enabled==true"))   # review 2026-09-04 A7 item 6
    push!(F, Field("physics.nudge.file", "string", "", "Physics", "Observed-track ensemble nudging (HT): CSV (time, lon_west, lon_east, lat_south, lat_north) — the in-band parcel ensemble's centroid is rigid-shift relaxed toward the track (structure preserved); free-runs outside the track span"))
    push!(F, Field("physics.nudge.tau", "string", "3h", "Physics", "Nudging relaxation time"; unit="duration", visible_if="physics.nudge.file!="))
    push!(F, Field("physics.nudge.z_bottom", "string", "18km", "Physics", "Nudge only parcels in this height band — bottom"; unit="length", visible_if="physics.nudge.file!="))
    push!(F, Field("physics.nudge.z_top", "string", "40km", "Physics", "Nudge height band — top"; unit="length", visible_if="physics.nudge.file!="))
    push!(F, Field("physics.wind_bias.u", "number", 0.0, "Physics", "Constant eastward drift [m/s] added to parcel motion inside the wind-bias height band (stratospheric wind tuning; archived met untouched)"))
    push!(F, Field("physics.wind_bias.v", "number", 0.0, "Physics", "Constant northward drift [m/s] inside the band"))
    push!(F, Field("physics.wind_bias.z_bottom", "string", "18km", "Physics", "Wind-bias height band — bottom"; unit="length"))
    push!(F, Field("physics.wind_bias.z_top", "string", "40km", "Physics", "Wind-bias height band — top"; unit="length"))
    push!(F, Field("output.quantity", "enum", "concentration", "Output", "Quantity written"; options=["concentration","mixing_ratio","both"]))
    push!(F, Field("output.name", "string", "", "Output", "Output file base name (grid_<name>.nc); empty = the run name"))   # review 2026-09-04 A7 item 6
    push!(F, Field("output.kernel", "enum", "uniform", "Output", "Gridding kernel: none = nearest cell; uniform = FLEXPART 4-cell bilinear (comparison reference); adaptive = per-parcel Gaussian with kNN-spacing bandwidth (resolution set by trajectory count, no 3-h age rule); covariance = per-parcel 2×2 covariance evolved by wind shear + turbulent growth, read only at evaluation (physical footprint, no kNN pass)";
                   options=["none","uniform","adaptive","covariance"]))
    push!(F, Field("output.resolution", "bool", false, "Output", "Write the R(x) true-resolution field (mass-weighted mean kernel width per cell): R < grid spacing → resolved; R > grid spacing → raise N, not the grid"; visible_if="output.kernel==adaptive"))
    push!(F, Field("output.covariance.mode", "enum", "physical", "Output", "physical = σ from P directly (adds relative dispersion on top of the ensemble — double-counts; no kNN pass); shape = kernel area from the kNN spacing, anisotropy/tilt from P (no double counting; needs the kNN pass)";
                   options=["physical","shape"], visible_if="output.kernel==covariance"))
    push!(F, Field("output.covariance.sigma0", "length", 50.0, "Output", "Covariance kernel: release footprint σ₀ (P init = σ₀²·I)"; unit="length", visible_if="output.kernel==covariance"))
    push!(F, Field("output.covariance.sigma_v", "number", 0.5, "Output", "Covariance kernel: unresolved relative-dispersion velocity σv [m/s] (turbulent growth D = σv²·T_L per axis)"; visible_if="output.kernel==covariance"))
    push!(F, Field("output.covariance.tl", "string", "5min", "Output", "Covariance kernel: Lagrangian time scale T_L of the turbulent growth"; unit="duration", visible_if="output.kernel==covariance"))
    push!(F, Field("output.covariance.sigma_max", "string", "500km", "Output", "Per-axis growth cap on the parcel covariance (keeps P bounded in persistent deformation)"; unit="length", visible_if="output.kernel==covariance"))
    push!(F, Field("output.covariance.stride", "int", 1, "Output", "Covariance kernel: evolve P every N sync steps (deformation is smooth on the met time scale — 4–8 is safe and cheaper)"; visible_if="output.kernel==covariance"))
    push!(F, Field("output.bandwidth_floor", "enum", "off", "Output", "Physical bandwidth floor: per-parcel kernel width ≥ Taylor dispersion σ(age) — the width of atmosphere a parcel represents; bounds near-source single-parcel receptor spikes and gives just-released parcels an honest footprint";
                   options=["off","dispersion"], visible_if="output.kernel==adaptive"))
    push!(F, Field("output.floor_sigma_v", "number", 0.5, "Output", "Dispersion-floor velocity σv [m/s]: σ(age) = σv·age → σv·√(2·T_L·age); 0.5 m/s → ≈0.7 km at 1 h, 3.6 km at 24 h"; visible_if="output.bandwidth_floor==dispersion"))
    push!(F, Field("output.floor_tl", "string", "5min", "Output", "Dispersion-floor Lagrangian time scale T_L (ballistic→diffusive crossover)"; unit="duration", visible_if="output.bandwidth_floor==dispersion"))
    push!(F, Field("output.receptors.file", "string", "", "Output", "Gridless receptor stations (DATEM list: id lat lon) — pointwise KDE time series at the station coordinates, no cell-average representativeness error; writes receptors_<name>.csv with C and R per station (adaptive or covariance kernel)"))
    push!(F, Field("output.receptors.layer_bottom", "string", "0m", "Output", "Receptor sampling layer bottom (AGL)"; unit="length", visible_if="output.receptors.file!="))
    push!(F, Field("output.receptors.layer_top", "string", "100m", "Output", "Receptor sampling layer top (AGL)"; unit="length", visible_if="output.receptors.file!="))
    push!(F, Field("output.receptors.min_bandwidth", "string", "", "Output", "Fixed representativeness floor on the receptor KDE bandwidth (empty = off; e.g. 5km tames near-source single-parcel spikes at the cost of native resolution)"; unit="length", visible_if="output.receptors.file!="))
    push!(F, Field("output.trajectory.enabled", "bool", false, "Output", "Per-parcel trajectory CSV: position, mass/age, chemistry and met sampled at the parcel"))
    push!(F, Field("output.trajectory.interval", "string", "", "Output", "Trajectory sampling cadence (empty = every sync step)"; unit="duration", visible_if="output.trajectory.enabled==true"))
    push!(F, Field("output.trajectory.max_parcels", "int", 10000, "Output", "Random subsample cap — at most this many parcels are written"; visible_if="output.trajectory.enabled==true"))
    push!(F, Field("output.trajectory.species", "list", String[], "Output", "Restrict the chemistry columns to these species ids (empty = all)"; visible_if="output.trajectory.enabled==true"))
    push!(F, Field("output.trajectory.format", "enum", "full", "Output", "Row layout: full = headered CSV (position, mass/age, met); demovor = the legacy DEMOVOR/AROME columns, no header, %.5g reals — yyyy-mm-dd,HH:MM:SS,id,lon,lat,height_m,tke,pblh_m";
                   options=["full","demovor"], visible_if="output.trajectory.enabled==true"))
    push!(F, Field("output.per_release", "bool", true, "Output", "Keep per-release-box dimension"))
    push!(F, Field("output.per_species", "bool", true, "Output", "Per-species variables"))
    push!(F, Field("output.group", "enum", "", "Output", "Output-column grouping: species (one column per species — grid memory independent of release-box count), release (one per box — needed for per-box backward sensitivities). Empty = auto: species forward, release backward"; options=["", "species", "release"]))
    push!(F, Field("output.aod", "bool", false, "Output", "Write AOD(550nm) for aerosol species"))
    push!(F, Field("output.reference", "string", "", "Output", "Passive-reference comparison: name of a [chemistry.init] source (e.g. \"cams\") to compare the reacting chemistry against — writes X (reacted) / X_pass (transport-only twin) / X_ref (source field) per species + a per-frame bias/rmse line"))
    push!(F, Field("output.tendencies", "list", "", "Output", "Per-parcel chemistry-only tendency: list of mechanism species (e.g. [\"O3\",\"NO2\"]) whose in-situ photochemical production/loss is accumulated along each trajectory (transport excluded) and gridded as <species>_tend [ppb] in the chemistry NetCDF — the per-trajectory attribution of chemistry vs transport"))
    push!(F, Field("output.path", "string", "juliaout/output/", "Output", "Output directory"))
    push!(F, Field("output.interval", "string", "1h", "Output", "Output interval"; unit="duration"))
    push!(F, Field("output.average", "string", "1h", "Output", "Output averaging window"; unit="duration"))
    push!(F, Field("output.sample", "string", "15min", "Output", "Output sampling step"; unit="duration"))
    push!(F, Field("output.grid.crs", "enum", "lonlat", "Output grid", "Horizontal CRS: lonlat (deg) or a nest's native projection (m)"; options=["lonlat","native"]))
    push!(F, Field("output.grid.nest", "int", 1, "Output grid", "Native CRS: which nest's projection"; visible_if="output.grid.crs==native"))
    push!(F, Field("output.grid.vertical.kind", "enum", "height_agl", "Output grid", "Vertical coordinate"; options=["height_agl","pressure","model_level","theta"]))
    push!(F, Field("output.grid.vertical.levels", "list", Float64[], "Output grid", "Output level boundaries in the vertical coordinate's units (m AGL / hPa / model level / K); empty = use heights"))
    push!(F, Field("output.grid.lon0", "number", -180.0, "Output grid", "Origin x: lower-left longitude (deg) or x0 (m, native)"; unit="deg"))
    push!(F, Field("output.grid.lat0", "number", -90.0, "Output grid", "Lower-left latitude"; unit="deg"))
    push!(F, Field("output.grid.nx", "int", 360, "Output grid", "Longitude cells"))
    push!(F, Field("output.grid.ny", "int", 180, "Output grid", "Latitude cells"))
    push!(F, Field("output.grid.dx", "number", 1.0, "Output grid", "Longitude spacing"; unit="deg"))
    push!(F, Field("output.grid.dy", "number", 1.0, "Output grid", "Latitude spacing"; unit="deg"))
    push!(F, Field("output.grid.heights", "list", [0.0,100.0,500.0,1000.0,3000.0,10000.0], "Output grid", "Layer edge heights"; unit="length"))

    push!(F, Field("chemistry.enabled", "bool", false, "Chemistry", "Apply gas-phase chemistry after transport"))
    push!(F, Field("chemistry.mechanism", "enum", "reducedqss", "Chemistry", "Reaction set (cb05 = C-IFS/CAMS)";
                   options=["reduced","reducedqss","bb12","cb05","strat","strat_brox"], visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.timestep", "string", "1h", "Chemistry", "Chemistry operator-split step (multiple of run.timestep)"; unit="duration", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.photolysis_timestep", "string", "", "Chemistry", "J-value refresh step (empty → = chemistry.timestep)"; unit="duration", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.init", "enum", "background", "Chemistry", "Initial composition source";
                   options=["background","cams"], visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.boundary", "bool", true, "Chemistry", "Re-seed newly released parcels from CAMS"; visible_if="chemistry.init==cams"))
    push!(F, Field("chemistry.lat_deg", "number", 45.0, "Chemistry", "Fallback latitude for photolysis"; unit="deg", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.h2o_ppb", "number", 1.0e4, "Chemistry", "Water vapour mixing ratio"; unit="ppb", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.precision", "enum", "fp64", "Chemistry", "Solver arithmetic for the bb12/cb05-rb23 fast paths (fp32: single precision, tolerances auto-clamped — bb12 1e-4, cb05 1e-5; see docs/GPU_CHEMISTRY.md)";
                   options=["fp64","fp32"], visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.solver", "enum", "rb23", "Chemistry", "CB05 integrator: rb23 = analytic-Jacobian Rosenbrock23 (3-5x faster, at least as accurate); fbdf = legacy generic path";
                   options=["rb23","fbdf"], visible_if="chemistry.mechanism==cb05"))
    push!(F, Field("chemistry.doy", "int", 0, "Chemistry", "Day of year (0 → from run start)"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.reltol", "number", 1.0e-6, "Chemistry", "ODE solver relative tolerance"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.abstol", "number", 1.0e-6, "Chemistry", "ODE solver absolute tolerance (radicals: HO2 ~1e-4 ppb, night NO3 ~5e-6 ppb — 1e-5 leaves them uncontrolled; don't loosen for CAMS work)"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.cams.dir", "string", "", "Chemistry", "CAMS GRIB directory or file"; visible_if="chemistry.init==cams"))
    push!(F, Field("chemistry.cams.species", "list", String[], "Chemistry", "CAMS shortNames to load (empty → all)"; visible_if="chemistry.init==cams"))
    push!(F, Field("chemistry.photolysis", "enum", "analytic", "Chemistry", "Photolysis: analytic cos-SZA or TUV lookup table (cb05)";
                   options=["analytic","lut"], visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.lut_path", "string", "", "Chemistry", "Photolysis LUT NetCDF (empty → TUV/cb05_photolysis_lut.nc)"; visible_if="chemistry.photolysis==lut"))
    push!(F, Field("chemistry.o3col_DU", "number", 300.0, "Chemistry", "O3 column for the LUT (CAMS-derived per-parcel when init=cams)"; unit="DU", visible_if="chemistry.photolysis==lut"))
    push!(F, Field("chemistry.albedo", "number", 0.05, "Chemistry", "Surface UV albedo for the LUT"; visible_if="chemistry.photolysis==lut"))
    push!(F, Field("chemistry.ice_reff", "string", "10um", "Chemistry", "Crystal r_eff of the PHYSICS-carried ICE species for the heterogeneous surface A = 3·m/(ρ·r_eff) — the chemistry reads the transported ice instead of re-deriving an equilibrium"; unit="length", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.ice_shield", "string", "0m", "Chemistry", "Self-shielding length of the plume for the strat photolysis: τ = (A/2)·L, J·exp(−τ/cosSZA). Core (dense) and edge (thin) parcels then see different J. 0 = clear-sky"; unit="length", visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.h2o_prognostic", "bool", false, "Chemistry", "SWV-prognostic: parcel keeps its own H2O budget after first met seeding (HT water studies)"))
    push!(F, Field("chemistry.ch4_oxidation", "bool", false, "Chemistry", "In-situ H2O from CH4 oxidation (strat mechanisms; 2 H2O per CH4; needs h2o_prognostic)"; visible_if="chemistry.h2o_prognostic==true"))
    push!(F, Field("chemistry.cloud", "bool", true, "Chemistry", "Cloud correction of photolysis from the archived 3-D cloud fraction (CLD companions; clear-sky when absent)"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.aerosols.enabled", "bool", false, "Chemistry", "Carry the CAMS aerosol fields (aermr01–20) on the parcels: Lagrangian aerosol transport, per-parcel het surface area, AOD550 output. Needs init=cams"; visible_if="chemistry.init==cams"))
    push!(F, Field("chemistry.aerosols.het", "bool", true, "Chemistry", "Per-parcel wet aerosol surface drives HO2/N2O5/SO2 uptake"; visible_if="chemistry.aerosols.enabled==true"))
    push!(F, Field("chemistry.aerosols.aod", "bool", true, "Chemistry", "Write the 550 nm AOD field per output frame"; visible_if="chemistry.aerosols.enabled==true"))
    push!(F, Field("chemistry.aerosols.rh_growth", "bool", true, "Chemistry", "RH-grown mass-extinction in the overhead AOD550 column (per-level κ-Köhler growth + wet Mie enhancement); false = dry v1 MEC"; visible_if="chemistry.aerosols.enabled==true"))
    push!(F, Field("chemistry.het.awet", "number", 0.0, "Chemistry", "Prescribed background wet-aerosol surface density for heterogeneous uptake (HO2, N2O5, SO2 …); 0 = off (or per-parcel A_wet from aerosols.het). Marine BL ≈ 1e-4"; unit="m²/m³", visible_if="chemistry.enabled==true"))
    # review 2026-09-04 A7 item 6: stratospheric sulfate / PSC / lightning keys read by toml_io.jl, were missing here
    push!(F, Field("chemistry.sulfate_sad", "bool", false, "Chemistry", "Parcel carries its own sulfate (+SO2 → sulfate) whose wet surface area sets the strat heterogeneous rates (strat mechanisms)"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.sulfate_gmr", "string", "0.35um", "Chemistry", "Sulfate lognormal geometric mean radius"; unit="length", visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.sulfate_sigma_g", "number", 1.3, "Chemistry", "Sulfate lognormal geometric standard deviation"; visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.sulfate_kappa", "number", 0.5, "Chemistry", "Sulfate κ-Köhler hygroscopicity"; visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.sulfate_so2_tau", "string", "17.1d", "Chemistry", "SO2 → sulfate e-folding time (prescribed; ignored when sulfate_so2_oh is on)"; unit="duration", visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.sulfate_so2_oh", "bool", false, "Chemistry", "SO2 oxidation by the solved OH (JPL SO2+OH+M) instead of the prescribed τ"; visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.sulfate_so2_oh_feedback", "bool", true, "Chemistry", "EXPERIMENTAL: let the SO2 sink consume OH (OH is fast-QSS; see coupler.jl)"; visible_if="chemistry.sulfate_so2_oh==true"))
    push!(F, Field("chemistry.so4_init", "number", 0.0, "Chemistry", "Global initial sulfate loading [kg/m³] on every parcel (prefer [release.chem_init].so4_ugm3)"; unit="kg/m3", visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.so2_init", "number", 0.0, "Chemistry", "Global initial SO2 loading [kg/m³] on every parcel (prefer [release.chem_init].so2_ppb)"; unit="kg/m3", visible_if="chemistry.sulfate_sad==true"))
    push!(F, Field("chemistry.psc_dehydration", "bool", false, "Chemistry", "PSC / cold-point dehydration sink on the prognostic water: freeze-out above RHi = psc_rhi, sedimenting on psc_tau (needs h2o_prognostic, strat mechanism)"; visible_if="chemistry.h2o_prognostic==true"))
    push!(F, Field("chemistry.psc_rhi", "number", 1.0, "Chemistry", "RHi threshold of the dehydration sink"; visible_if="chemistry.psc_dehydration==true"))
    push!(F, Field("chemistry.psc_tau", "string", "6h", "Chemistry", "Sedimentation time scale of the dehydration sink"; unit="duration", visible_if="chemistry.psc_dehydration==true"))
    push!(F, Field("chemistry.psc_denitrification", "bool", false, "Chemistry", "NAT denitrification: permanent HNO3 removal below T_NAT on psc_nat_tau (strat mechanism)"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.psc_nat_tau", "string", "5d", "Chemistry", "Large-NAT sedimentation time scale"; unit="duration", visible_if="chemistry.psc_denitrification==true"))
    push!(F, Field("chemistry.lightning.enabled", "bool", false, "Chemistry", "Lightning NOx: Price & Rind (1992) flash rate from the archived updraft-mass-flux cloud top → NO into UT parcels"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.lightning.yield", "number", 250.0, "Chemistry", "mol NO per flash (tunable)"; visible_if="chemistry.lightning.enabled==true"))
    push!(F, Field("chemistry.lightning.umf_thresh", "number", 0.01, "Chemistry", "Updraft mass-flux threshold [kg m⁻² s⁻¹] defining the convective cloud top"; visible_if="chemistry.lightning.enabled==true"))
    push!(F, Field("chemistry.lightning.h_min", "number", 5.0, "Chemistry", "Minimum cloud-top height [km] for lightning"; visible_if="chemistry.lightning.enabled==true"))
    push!(F, Field("chemistry.lightning.ut_lo", "string", "6km", "Chemistry", "Inject NO into parcels above this height"; unit="length", visible_if="chemistry.lightning.enabled==true"))
    push!(F, Field("chemistry.lightning.ut_hi", "string", "14km", "Chemistry", "… and below this height"; unit="length", visible_if="chemistry.lightning.enabled==true"))
    # [chemistry.init] as a TABLE (default/sources.<name>/[[species]] with kind/
    # from/predictor/curve/…) is TOML-only: keyed sub-tables and nested arrays
    # don't fit the flat form — the init enum above covers the common sugar.
    push!(F, Field("chemistry.mixing.scheme", "enum", "off", "Chemistry", "Parcel–parcel composition exchange during transport: off = each parcel is a closed chemical box; pairwise = conservative mass-weighted exchange between co-located parcels (in-plume chemistry needs this)";
                   options=["off","pairwise"], visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.mixing.tau", "duration", "30min", "Chemistry", "Mixing relaxation timescale τ (μ = 1−exp(−dt/τ) per chemistry step)"; visible_if="chemistry.mixing.scheme==pairwise"))
    push!(F, Field("chemistry.mixing.cells", "int", 1, "Chemistry", "Neighbourhood size, met grid cells"; visible_if="chemistry.mixing.scheme==pairwise"))
    push!(F, Field("chemistry.mixing.levels", "int", 1, "Chemistry", "Neighbourhood size, model levels"; visible_if="chemistry.mixing.scheme==pairwise"))
    push!(F, Field("chemistry.emissions.dir", "string", "", "Chemistry", "Surface-emission inventory directory (CAMS-GLOB-* NetCDFs); empty = no emissions. Requires a domainfill release"; visible_if="chemistry.enabled==true"))
    push!(F, Field("chemistry.emissions.files", "list", String[], "Chemistry", "Inventory basename glob patterns (empty → all *.nc in dir)"; visible_if="chemistry.emissions.dir!="))
    push!(F, Field("chemistry.emissions.res", "number", 0.4, "Chemistry", "Emission working-grid resolution"; unit="deg", visible_if="chemistry.emissions.dir!="))
    push!(F, Field("chemistry.emissions.inject_depth", "string", "150m", "Chemistry", "Inject emissions into parcels below this height"; unit="length", visible_if="chemistry.emissions.dir!="))
    push!(F, Field("chemistry.emissions.nox_no", "number", 0.9, "Chemistry", "Mole fraction of NOx emitted as NO (rest NO2)"; visible_if="chemistry.emissions.dir!="))
    push!(F, Field("chemistry.emissions.diurnal", "enum", "none", "Chemistry", "Diurnal cycle applied to the daily/monthly-mean flux by local solar time (24-h mean preserved): none = flat; gaussian = afternoon-peaked Gaussian; mu2011 = tabulated savanna cycle (Mu et al. 2011). Matters for photochemistry — fires peak in phase with peak OH/photolysis"; options=["none","gaussian","mu2011"], visible_if="chemistry.emissions.dir!="))
    push!(F, Field("chemistry.emissions.diurnal_peak", "number", 13.5, "Chemistry", "Peak local solar time [h] for the gaussian diurnal cycle"; visible_if="chemistry.emissions.diurnal==gaussian"))
    push!(F, Field("chemistry.emissions.diurnal_sigma", "number", 2.0, "Chemistry", "Width [h] of the gaussian diurnal cycle"; visible_if="chemistry.emissions.diurnal==gaussian"))

    # review 2026-09-04 A7 item 1: [validation] (task = "validation") is consumed by the
    # launcher (scripts/easlee.jl → env → validation/validate_captex.jl), not by the reader
    push!(F, Field("validation.obs_dir", "string", "", "Validation", "DATEM observation directory (ETEX / CAPTEX)"; visible_if="task==validation"))
    push!(F, Field("validation.mode", "enum", "run", "Validation", "run = simulate + score; selftest = score an existing output"; options=["run","selftest"], visible_if="task==validation"))
    push!(F, Field("validation.direction", "enum", "forward", "Validation", "Forward plume vs backward receptor-footprint scoring"; options=["forward","backward"], visible_if="task==validation"))
    push!(F, Field("validation.release", "string", "all", "Validation", "Which tracer release(s) to run: all, or a release number"; visible_if="task==validation"))
    push!(F, Field("validation.maxback", "int", 8, "Validation", "Backward mode: maximum number of receptor stations run"; visible_if="validation.direction==backward"))
    push!(F, Field("validation.rundays", "int", 4, "Validation", "Simulation length after each release [days]"; visible_if="task==validation"))
    push!(F, Field("validation.dxout", "number", 0.25, "Validation", "Output grid spacing [deg]"; unit="deg", visible_if="task==validation"))
    push!(F, Field("validation.relz", "list", [0.0, 20.0], "Validation", "Release layer [z_bottom, z_top] m AGL (FLEXPART-11 used 8–20 m)"; unit="m", visible_if="task==validation"))
    push!(F, Field("validation.ctl", "number", 5.0, "Validation", "CTL for the validation run (overrides physics.ctl; positive integer value)"; visible_if="task==validation"))
    push!(F, Field("validation.ifine", "int", 4, "Validation", "IFINE for the validation run (overrides physics.ifine)"; visible_if="task==validation"))
    return F
end

# ── TOML-only spellings ───────────────────────────────────────────────────────
# review 2026-09-04 A7 item 1: keys the TOML reader (toml_io.jl) accepts that the
# flat GUI form cannot render (keyed sub-tables, wildcards, vector aliases, a
# table under a key the form treats as a scalar).  They are part of the schema
# for the unknown-key check; "*" matches any single path segment, "[]" marks an
# array-of-tables element as in `Field.path`.  Path => help.
function schema_extra_paths()::Vector{Pair{String,String}}
    Pair{String,String}[
        # [species.settling] TABLE (preferred) / [species.fall] (legacy alias, review item 2):
        # gravitational fall of a CONDENSED gas-kind species (ice/ash shadow population)
        "species[].settling.r_eff"  => "Effective particle radius of the falling condensed species (0 = off); Stokes–Cunningham v_t at the parcel's T,p",
        "species[].settling.rho"    => "Particle density [kg/m³] of the falling condensed species (ice: 917)",
        "species[].settling.enabled"=> "Aerosol size-bin settling toggle when [species.settling] is written as a table (default true)",
        "species[].fall.r_eff"      => "Legacy alias of species[].settling.r_eff",
        "species[].fall.rho"        => "Legacy alias of species[].settling.rho",
        # [release.volcanic] phases + [release.chem_init] free species keys
        "release[].volcanic.phases"        => "Injection phases [[t0, t1, z_bottom, z_top, frac (, mer, temperature, water_frac)], …]",
        "release[].chem_init.h2o_from_plume" => "Seed the prognostic chemistry water from the column model's detrained water fraction (needs plume_rise mode=parcels + [chemistry].h2o_prognostic)",
        "release[].chem_init.*"            => "Any mechanism species name → its initial column [ppb] on this release's parcels",
        # [release.plume_rise.vent_mix]: `vent_mix = true` or a sub-table
        "release[].plume_rise.vent_mix"               => "Radial magma–seawater vent split (true = defaults, or a [..vent_mix] table)",
        "release[].plume_rise.vent_mix.magma_t"       => "Magma temperature [K]",
        "release[].plume_rise.vent_mix.magma_water"   => "Magmatic volatile weight fraction of the melt",
        "release[].plume_rise.vent_mix.seawater_core" => "Seawater mass fraction on the axis",
        "release[].plume_rise.vent_mix.seawater_edge" => "Seawater mass fraction at the edge",
        "release[].plume_rise.vent_mix.mix_exponent"  => "Radial profile exponent",
        # vector aliases of the scalar height-band keys
        "physics.nudge.z"     => "Height band [z_bottom, z_top] in metres (alias of z_bottom/z_top)",
        "physics.wind_bias.z" => "Height band [z_bottom, z_top] in metres (alias of z_bottom/z_top)",
        # output grid vector form + receptor stations + table forms of reference/tendencies
        "output.grid.origin"  => "[lon0, lat0] (alias of the scalar keys)",
        "output.grid.spacing" => "[dx, dy] (alias of the scalar keys)",
        "output.grid.size"    => "[nx, ny] (alias of the scalar keys)",
        "output.receptors.station[].id"  => "Inline receptor station id ([[output.receptors.station]])",
        "output.receptors.station[].lon" => "Inline receptor station longitude [deg]",
        "output.receptors.station[].lat" => "Inline receptor station latitude [deg]",
        "output.receptors[].id"  => "Inline receptor station id (plain [[output.receptors]] form)",
        "output.receptors[].lon" => "Inline receptor station longitude [deg]",
        "output.receptors[].lat" => "Inline receptor station latitude [deg]",
        "output.reference.source"  => "Table form of output.reference: the [chemistry.init] source name",
        "output.reference.compare" => "Table form of output.reference (alias of source)",
        "output.tendencies.species" => "Table form of output.tendencies: the species list",
        # [chemistry.init] declarative block (CHEM_INIT.md): keyed sources + species assignments
        "chemistry.init.default"           => "Source for species not listed in [[chemistry.init.species]] (\"\" = zeros)",
        "chemistry.init.sources.*.kind"    => "Source kind: grib | netcdf | constant",
        "chemistry.init.sources.*.dir"     => "Source file directory",
        "chemistry.init.sources.*.units"   => "Source units: mmr | vmr | ppb",
        "chemistry.init.sources.*.vars.*"  => "Species → file variable name override",
        "chemistry.init.species[].name"      => "Mechanism species name",
        "chemistry.init.species[].from"      => "A source name, or correlation | steady_state | family | constant | profile",
        "chemistry.init.species[].predictor" => "correlation: predictor species",
        "chemistry.init.species[].curve"     => "correlation: x,y CSV",
        "chemistry.init.species[].total"     => "family: the family-total species/source",
        "chemistry.init.species[].regime"    => "steady_state: any | day",
        "chemistry.init.species[].value"     => "constant [ppb]",
        "chemistry.init.species[].profile"   => "profile: z/p CSV",
    ]
end

# Every accepted dotted path (Field paths + TOML-only spellings), "[]" stripped
# — the flattened form the unknown-key check compares against.
schema_known_paths()::Vector{String} =
    unique(vcat([replace(f.path, "[]" => "") for f in schema_fields()],
                [replace(first(p), "[]" => "") for p in schema_extra_paths()]))

# JSON-serialisable schema for the GUI (plain Dicts; the server adds JSON3).
function schema_json()
    [Dict("path"=>f.path, "type"=>f.type, "default"=>f.default, "group"=>f.group,
          "help"=>f.help, "unit"=>f.unit, "options"=>f.options, "visible_if"=>f.visible_if)
     for f in schema_fields()]
end

end # module Schema
