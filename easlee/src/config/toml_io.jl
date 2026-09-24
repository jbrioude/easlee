
# ── TOML ↔ RunSpec ────────────────────────────────────────────────────────────
# Canonical on-disk format.  Reading applies unit parsing (durations→s,
# lengths→m, masses→kg) and fills defaults; writing emits a commented, grouped
# file.  TOML.jl is in the Julia stdlib.

module ConfigTOML

using Dates
using TOML
using ..RunSpecs
using ..Units
using ..SpeciesDB: get_species

export read_runspec, write_runspec, runspec_to_dict, dict_to_runspec
export unknown_config_keys, check_config_keys

# ── helpers ───────────────────────────────────────────────────────────────────
_get(d, k, def) = haskey(d, k) ? d[k] : def
# vent_mix accepts `vent_mix = true` (defaults) or a [ ... .vent_mix] subtable
_ventmix(pr) = _get(pr, "vent_mix", Dict{String,Any}()) isa AbstractDict ?
               _get(pr, "vent_mix", Dict{String,Any}()) : Dict{String,Any}()
_dt(x::DateTime) = x
_dt(x::AbstractString) = DateTime(x)          # ISO 8601
_dt(x::Date) = DateTime(x)

# A TOML datetime may come back as DateTime, Date, or String depending on the
# value; normalise.
_asdt(x) = x isa DateTime ? x : x isa Date ? DateTime(x) : DateTime(String(x))

# Parse a unit value only when present AND non-blank.  An empty string (a GUI
# field left blank, or an explicit "" in TOML) means "unset" → return the
# default, rather than throwing in parse_duration_s/parse_length_m.
_blank(x) = x isa AbstractString && isempty(strip(String(x)))
_durget(d, k, def) = (haskey(d, k) && !_blank(d[k])) ? parse_duration_s(d[k]) : def
_lenget(d, k, def) = (haskey(d, k) && !_blank(d[k])) ? parse_length_m(d[k]) : def

# ── dict → RunSpec ────────────────────────────────────────────────────────────
function dict_to_runspec(d::AbstractDict)::RunSpec
    s = RunSpec()
    s.name      = String(_get(d, "name", s.name))
    s.task      = lowercase(String(_get(d, "task", s.task)))   # transport|chemistry|validation
    run = _get(d, "run", Dict{String,Any}())
    s.direction  = String(_get(run, "direction", s.direction))
    s.wind       = String(_get(run, "wind", s.wind))
    haskey(run, "start") && (s.start = _asdt(run["start"]))
    haskey(run, "stop")  && (s.stop  = _asdt(run["stop"]))
    s.timestep_s = _durget(run, "timestep", s.timestep_s)
    s.benchmark = Bool(_get(run, "benchmark", false))
    s.seed      = Int(_get(run, "seed", 0))
    s.checkpoint = String(_get(run, "checkpoint", ""))
    s.restart    = String(_get(run, "restart", ""))
    s.parcel_dump = String(_get(run, "parcel_dump", ""))

    # meteorology (array of tables)
    for m in _get(d, "meteorology", Any[])
        push!(s.meteorology, MetSource(
            model     = String(_get(m, "model", "ERA5")),   # tolerant of half-filled GUI rows;
            path      = String(_get(m, "path", "")),        # real gaps surface at validate/build
            variables = String(_get(m, "variables", "")),
            available = String(_get(m, "available", "auto")),
            # crop accepts "auto" | "off" | [lon0, lon1, lat0, lat1]
            crop      = _get(m, "crop", "auto") isa AbstractVector ? "explicit" :
                        lowercase(String(_get(m, "crop", "auto"))),
            crop_box  = _get(m, "crop", "auto") isa AbstractVector ?
                        Float64.(collect(m["crop"])) : Float64[],
            max_met_gap_s = _durget(m, "max_gap", 6 * 3600.0),
        ))
    end

    # species
    for sp in _get(d, "species", Any[])
        size = _get(sp, "size", Dict{String,Any}())
        # [species.ice]: table {rhi_crit, r_eff, h_layer, q0} or `ice = true`;
        # ABSENT key must yield nothing (an empty-Dict default silently
        # enabled ice for every species — caught by test_ice_physics [5])
        icg = _get(sp, "ice", nothing)
        icg isa Bool && (icg = icg ? Dict{String,Any}() : nothing)
        # review 2026-09-04 A7 item 2: the condensed-species fall table is
        # [species.settling] (preferred; matches runspec/build_run/the run log) with
        # [species.fall] kept as an alias.  `settling = true|false` stays the
        # aerosol size-bin toggle; a table carries it as `enabled`.
        stl = _get(sp, "settling", nothing)
        fall = stl isa AbstractDict ? stl : _get(sp, "fall", Dict{String,Any}())
        push!(s.species, Species(
            id                = String(_get(sp, "id", "spec")),
            kind              = String(_get(sp, "kind", "gas")),
            specnum           = Int(_get(sp, "specnum", 0)),
            size_distribution = String(_get(size, "distribution", "lognormal")),
            dg_m              = _lenget(size, "dg", 1.0e-6),
            sigma_g           = Float64(_get(size, "sigma_g", 2.0)),
            nbins             = Int(_get(size, "nbins", 8)),
            bin_edges_m       = haskey(size, "edges") ? Float64[parse_length_m(e) for e in size["edges"] if !_blank(e)] : Float64[],
            density           = Float64(_get(sp, "density", 2000.0)),
            kappa             = Float64(_get(sp, "kappa", -1.0)),
            aerosol_type      = String(_get(sp, "aerosol_type", "dust")),
            refractive_index  = String(_get(sp, "refractive_index", "auto")),
            settling          = stl isa AbstractDict ? Bool(_get(stl, "enabled", true)) : Bool(something(stl, true)),
            molar_mass        = Float64(_get(sp, "molar_mass", 0.0)),   # review 2026-09-04 A7 item 3: was documented, never parsed
            dry_deposition    = Bool(_get(sp, "dry_deposition", false)),
            wet_deposition    = Bool(_get(sp, "wet_deposition", false)),
            decay_halflife_s  = _durget(sp, "decay_halflife", 0.0),
            # a hand-written [species.ice] table enables by presence; the GUI
            # writes `enabled` explicitly (its defaults must not enable)
            ice               = icg !== nothing && Bool(_get(icg, "enabled", true)),
            ice_rhi_crit      = icg === nothing ? 1.5 : Float64(_get(icg, "rhi_crit", 1.5)),
            ice_r_eff_m       = icg === nothing ? 20.0e-6 : _lenget(icg, "r_eff", 20.0e-6),
            ice_h_layer_m     = icg === nothing ? 500.0 : _lenget(icg, "h_layer", 500.0),
            ice_q0            = icg === nothing ? 0.0 : Float64(_get(icg, "q0", 0.0)),
            ice_scav_of       = String(_get(_get(sp, "ice_scav", Dict{String,Any}()), "of", "")),
            ice_scav_uptake   = Float64(_get(_get(sp, "ice_scav", Dict{String,Any}()), "uptake", 0.0)),
            convert_to        = String(_get(_get(sp, "convert", Dict{String,Any}()), "to", "")),
            convert_tau_s     = _durget(_get(sp, "convert", Dict{String,Any}()), "tau", 0.0),
            convert_mode      = lowercase(String(_get(_get(sp, "convert", Dict{String,Any}()), "mode", ""))),
            convert_rhi_crit  = Float64(_get(_get(sp, "convert", Dict{String,Any}()), "rhi_crit", 1.5)),
            convert_lam_sigma = Float64(_get(_get(sp, "convert", Dict{String,Any}()), "lamina_sigma", 0.0)),
            convert_sectional = Bool(_get(_get(sp, "convert", Dict{String,Any}()), "sectional", false)),
            settle_r_eff      = _lenget(fall, "r_eff", 0.0),
            settle_rho        = Float64(_get(fall, "rho", 917.0)),
            convert_mass_ratio = Float64(_get(_get(sp, "convert", Dict{String,Any}()), "mass_ratio", 0.0)),
        ))
    end

    # releases — `species` may be a single id or a list (→ one release each)
    for r in _get(d, "release", Any[])
        sp = _get(r, "species", "")
        specs = sp isa AbstractVector ? String.(sp) : [String(sp)]
        # lon/lat may be a scalar (point/profile) or a list (box/line/polygon)
        aslist(x) = x === nothing ? Float64[] : (x isa AbstractVector ? Float64.(x) : [Float64(x)])
        pr = _get(r, "plume_rise", Dict{String,Any}())
        kind = lowercase(String(_get(r, "kind", "release")))
        # a gridded/text inventory release always uses the inventory column layout
        fmt = String(_get(r, "file_format", kind == "inventory" ? "inventory" : "latlonlev"))
        kind == "inventory" && (fmt = "inventory")
        # volcanic column ([release.volcanic]): phases = [[t0, t1, z_bottom, z_top, frac], …]
        vo = _get(r, "volcanic", Dict{String,Any}())
        vphases = Tuple{DateTime,DateTime,Float64,Float64,Float64}[]
        vpulses = NTuple{3,Float64}[]
        for p in _get(vo, "phases", Any[])
            length(p) >= 5 || error("[release.volcanic] phase needs [t0, t1, z_bottom, z_top, frac] (+ optional mer, temperature, water_frac)")
            push!(vphases, (_asdt(p[1]), _asdt(p[2]),
                            parse_length_m(p[3]), parse_length_m(p[4]), Float64(p[5])))
            length(p) >= 8 && push!(vpulses, (Float64(p[6]), Float64(p[7]), Float64(p[8])))
        end
        isempty(vpulses) || length(vpulses) == length(vphases) ||
            error("[release.volcanic] pulse params (cols 6-8) must be given for ALL phases or none")
        # parcel allocation: an explicit count, or a mass QUANTUM — with
        # `mass_per_parcel` the count follows the released mass (n = ⌈M/q⌉),
        # so multi-release configs get a consistent kg-per-parcel resolution
        rmass = (haskey(r, "mass") && !_blank(r["mass"])) ? parse_mass_kg(r["mass"]) : 1.0
        qpp   = (haskey(r, "mass_per_parcel") && !_blank(r["mass_per_parcel"])) ?
                parse_mass_kg(r["mass_per_parcel"]) : 0.0
        # seaspray: mass_per_parcel is the SPAWN quantum and `particles` the
        # parcel-pool cap — never derive the count from the (unused) mass
        rpart = kind == "seaspray" ? Int(_get(r, "particles", 10000)) :
                qpp > 0.0 ? max(1, ceil(Int, rmass / qpp)) : Int(_get(r, "particles", 10000))
        qpp > 0.0 && rpart > 20_000_000 &&
            @warn "mass_per_parcel yields a very large parcel count" release=join(specs, ",") parcels=rpart
        # [release.chem_init]: release-scoped chemistry IC (species/key → value);
        # the boolean `h2o_from_plume` opts into the volcanic-plume water bridge.
        cinit = Dict{String,Float64}()
        h2o_plume = false
        for (k, v) in _get(r, "chem_init", Dict{String,Any}())
            if String(k) == "h2o_from_plume"
                h2o_plume = Bool(v)
            else
                cinit[String(k)] = Float64(v)
            end
        end
        for sid in specs
            push!(s.release, ReleaseBox(
                species   = sid,
                kind      = kind,
                geometry  = String(_get(r, "geometry", "box")),
                lon       = aslist(_get(r, "lon", nothing)),
                lat       = aslist(_get(r, "lat", nothing)),
                z         = haskey(r, "z") ? Float64[parse_length_m(z) for z in r["z"] if !_blank(z)] : [0.0, 100.0],
                zkind     = String(_get(r, "zkind", kind == "volcanic" ? "asl" : "agl")),
                time      = haskey(r, "time") ? DateTime[_asdt(t) for t in r["time"] if !_blank(t)] : DateTime[],
                particles = rpart,
                mass_per_parcel_kg = qpp,
                mass_kg   = rmass,
                npoints   = Int(_get(r, "npoints", 100)),
                comment   = String(_get(r, "comment", "")),
                from_file   = String(_get(r, "from_file", "")),
                file_format = fmt,
                volc_phases     = vphases,
                volc_pulses     = vpulses,
                volc_from_file  = String(_get(vo, "from_file", "")),
                volc_profile    = lowercase(String(_get(vo, "profile", "suzuki"))),
                volc_suzuki_a   = Float64(_get(vo, "suzuki_a", 4.0)),
                volc_suzuki_l   = Float64(_get(vo, "suzuki_lambda", 1.0)),
                plume_scheme    = lowercase(String(_get(pr, "scheme", "none"))),
                plume_frp_w     = Float64(_get(pr, "frp", 0.0)),
                plume_fire_area = Float64(_get(pr, "fire_area", 0.0)),
                plume_mer       = Float64(_get(pr, "mer", 0.0)),
                plume_temp_k    = Float64(_get(pr, "temperature", 800.0)),
                plume_water_frac= Float64(_get(pr, "water_frac", 0.05)),
                plume_profile   = String(_get(pr, "profile", "uniform")),
                plume_bottom    = String(_get(pr, "bottom", "pbl")),
                plume_mode      = lowercase(String(_get(pr, "mode", "layer"))),
                plume_washout   = Bool(_get(pr, "washout", true)),
                plume_entrainment = lowercase(String(_get(pr, "entrainment", "meanfield"))),
                plume_umbrella  = Bool(_get(pr, "umbrella_spread", false)),
                plume_spread_lambda = Float64(_get(pr, "spread_lambda", 0.2)),
                plume_umbrella_rsource = _lenget(pr, "umbrella_r_source", 10000.0),
                plume_umbrella_windshield = Bool(_get(pr, "umbrella_wind_shield", false)),
                plume_spread_umin = Float64(_get(pr, "spread_umin", 5.0)),
                plume_spread_tmax = _durget(pr, "spread_tmax", 43200.0),
                plume_spread_hmin = _lenget(pr, "spread_hmin", 0.0),
                plume_spread_zgain = Float64(_get(pr, "spread_zgain", 0.0)),
                plume_spread_zg_z0 = _lenget(pr, "spread_zgain_z0", 26.0e3),
                plume_spread_zg_z1 = _lenget(pr, "spread_zgain_z1", 30.0e3),
                plume_spread_azgain_n = Float64(_get(pr, "spread_azgain_n",
                                                _get(pr, "spread_azgain", 0.0))),
                plume_spread_azgain_s = Float64(_get(pr, "spread_azgain_s",
                                                _get(pr, "spread_azgain", 0.0))),
                plume_spread_az_z0 = _lenget(pr, "spread_azgain_ramp_z0", 0.0),
                plume_spread_az_z1 = _lenget(pr, "spread_azgain_ramp_z1", 0.0),
                plume_qc_crit   = Float64(_get(pr, "qc_crit", 5.0e-3)),
                plume_vapor_cap = Float64(_get(pr, "vapor_cap", Inf)),
                plume_umb_r_med = _lenget(pr, "umbrella_r_ice", 0.0),
                plume_umb_r_sg  = Float64(_get(pr, "umbrella_r_ice_sigma", 2.0)),
                plume_umb_h     = _lenget(pr, "umbrella_shed_depth", 2000.0),
                plume_umb_t     = _durget(pr, "umbrella_shed_time", 3600.0),
                plume_solid_r   = _lenget(pr, "tephra_r", 0.0),
                plume_vent_mix  = haskey(pr, "vent_mix") &&
                                  (_get(pr, "vent_mix", Dict{String,Any}()) !== false),
                plume_magma_t   = Float64(_get(_ventmix(pr), "magma_t", 1323.0)),
                plume_magma_water = Float64(_get(_ventmix(pr), "magma_water", 0.03)),
                plume_sw_core   = Float64(_get(_ventmix(pr), "seawater_core", 0.15)),
                plume_sw_edge   = Float64(_get(_ventmix(pr), "seawater_edge", 0.70)),
                plume_mix_exp   = Float64(_get(_ventmix(pr), "mix_exponent", 2.0)),
                spray_sst       = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "sst", -1.0)),
                spray_salinity  = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "salinity", 35.0)),
                spray_hs        = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "hs", -1.0)),
                spray_cd        = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "cd", -1.0)),
                spray_gust_sigma = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "gust_sigma", 0.0)),
                spray_rho       = Float64(_get(_get(r, "seaspray", Dict{String,Any}()), "rho_dry", 2017.0)),
                chem_init       = cinit,
                chem_h2o_plume  = h2o_plume,
            ))
        end
    end

    # physics
    ph = _get(d, "physics", Dict{String,Any}())
    s.physics.ctl   = Float64(_get(ph, "ctl", 5.0))
    s.physics.ifine = Int(_get(ph, "ifine", 4))
    s.physics.polar_stereo_lat = Float64(_get(ph, "polar_stereo_lat", 75.0))
    tb = _get(ph, "turbulence", Dict{String,Any}())
    tb isa AbstractString && (tb = Dict("scheme" => tb))      # allow turbulence = "hanna"
    s.physics.turbulence = Turbulence(
        scheme = String(_get(tb, "scheme", "hanna")),
        hanna_onthefly = Bool(_get(tb, "onthefly", false)),
        sigw_ft = Float64(_get(tb, "sigw_ft", 0.1)),
        kh_ft = Float64(_get(tb, "kh_ft", 0.0)),
        horizontal_tke = Bool(_get(tb, "horizontal_tke", false)),
        horizontal_hanna = Bool(_get(tb, "horizontal_hanna", false)),
        mesoscale = Bool(_get(tb, "mesoscale", false)),
        fturbmeso = Float64(_get(tb, "fturbmeso", 0.16)),
    )
    cv = _get(ph, "convection", Dict{String,Any}())
    cv isa AbstractString && (cv = Dict("scheme" => cv))
    s.physics.convection = Convection(
        scheme = String(_get(cv, "scheme", "off")),
        cape_min = Float64(_get(cv, "cape_min", 100.0)),
        cin_max  = Float64(_get(cv, "cin_max", 50.0)),
        tau_mix_s = _durget(cv, "tau_mix", 0.0),
        timestep_s = _durget(cv, "timestep", 0.0),
    )
    vt = _get(ph, "vertical", Dict{String,Any}())
    vt isa AbstractString && (vt = Dict("mode" => vt))
    s.physics.vertical = VerticalSpec(
        mode = String(_get(vt, "mode", "kinematic")),
        sigma_r = Float64(_get(vt, "sigma_r", 0.3)),
        sigma_t = Float64(_get(vt, "sigma_t", 0.1)),
        diabatic = Bool(_get(vt, "diabatic", true)),
        pt_winds = String(_get(vt, "pt_winds", "auto")),
        pt_pattern = String(_get(vt, "pt_pattern", "IFS_PT_MET")),
    )
    s.physics.surface_layer_wind = Bool(_get(ph, "surface_layer_wind", true))
    dp = _get(ph, "deposition", Dict{String,Any}())
    s.physics.deposition = Deposition(
        dry = Bool(_get(dp, "dry", false)),
        wet = Bool(_get(dp, "wet", false)),
    )
    mx = _get(ph, "mixing", Dict{String,Any}())
    mx isa AbstractString && (mx = Dict("scheme" => mx))      # allow mixing = "pairwise"
    s.physics.mixing = TransportMixing(
        scheme = String(_get(mx, "scheme", "off")),
        K_h = Float64(_get(mx, "K_h", 1000.0)),
        sigma0_h = _lenget(mx, "sigma0_h", 0.0),
        sigma_v = _lenget(mx, "sigma_v", 0.0),
        nsig = Float64(_get(mx, "nsig", 2.0)),
        max_neighbors = Int(_get(mx, "max_neighbors", 30)),
        background = Bool(_get(mx, "background", true)),
        every = Int(_get(mx, "every", 1)),
    )
    # observed-track ensemble nudging + wind-bias drift (HT thread)
    prd = _get(ph, "plume_radiation", Dict{String,Any}())
    let sps = _get(prd, "species", Any[])
        s.physics.prad_species = sps isa AbstractVector ? String.(sps) : [String(sps)]
    end
    pps = _get(prd, "partners", String[])
    s.physics.prad_partners = pps isa AbstractVector ? String.(pps) : [String(pps)]
    s.physics.prad_c_rad    = Float64(_get(prd, "c_rad", 0.002))
    s.physics.prad_lw_kappa = Float64(_get(prd, "lw_kappa", 0.0))
    s.physics.prad_h_ref    = _lenget(prd, "h_ref", 1000.0)
    s.physics.prad_q_bg     = Float64(_get(prd, "q_bg", 5.0))
    s.physics.prad_feedback = Bool(_get(prd, "feedback", true))
    s.physics.prad_tprime_tau = _durget(prd, "tprime_tau", 0.0)
    s.physics.prognostic_airm       = Bool(_get(ph, "prognostic_airm", false))
    s.physics.prognostic_airm_every = Int(_get(ph, "prognostic_airm_every", 12))
    # [tracers]: passive transport-history air-mass tags (top-level table)
    tr = _get(d, "tracers", Dict{String,Any}())
    s.physics.tags_enabled     = Bool(_get(tr, "enabled", false))
    s.physics.tags_strat_theta = Float64(_get(tr, "strat_theta", 380.0))
    s.physics.tags_tropo_origin = Bool(_get(tr, "tropo_origin", false))
    ng = _get(ph, "nudge", Dict{String,Any}())
    s.physics.nudge_file  = String(_get(ng, "file", ""))
    s.physics.nudge_tau_s = _durget(ng, "tau", 10800.0)
    s.physics.nudge_z     = Float64[_lenget(ng, "z_bottom", 18000.0),
                                    _lenget(ng, "z_top",    40000.0)]
    if haskey(ng, "z")
        zz = _get(ng, "z", Any[])
        length(zz) == 2 && (s.physics.nudge_z = Float64[Float64(zz[1]), Float64(zz[2])])
    end
    wb = _get(ph, "wind_bias", Dict{String,Any}())
    s.physics.windbias_u = Float64(_get(wb, "u", 0.0))
    s.physics.windbias_v = Float64(_get(wb, "v", 0.0))
    s.physics.windbias_z = Float64[_lenget(wb, "z_bottom", 18000.0),
                                   _lenget(wb, "z_top",    40000.0)]
    if haskey(wb, "z")
        zz = _get(wb, "z", Any[])
        length(zz) == 2 && (s.physics.windbias_z = Float64[Float64(zz[1]), Float64(zz[2])])
    end

    # output grid — accept BOTH field styles:
    #   origin = [x0,y0] / spacing = [dx,dy] / size = [nx,ny]   (vector form)
    #   lon0,lat0 / dx,dy / nx,ny                                (scalar form)
    out = _get(d, "output", Dict{String,Any}())
    g   = _get(out, "grid", Dict{String,Any}())
    origin  = _get(g, "origin", nothing)
    spacing = _get(g, "spacing", nothing)
    gsize   = _get(g, "size", nothing)
    gx0 = origin  !== nothing ? Float64(origin[1])  : Float64(_get(g, "lon0", -180.0))
    gy0 = origin  !== nothing ? Float64(origin[2])  : Float64(_get(g, "lat0", -90.0))
    gdx = spacing !== nothing ? Float64(spacing[1]) : Float64(_get(g, "dx", 1.0))
    gdy = spacing !== nothing ? Float64(spacing[2]) : Float64(_get(g, "dy", 1.0))
    gnx = gsize   !== nothing ? Int(gsize[1])       : Int(_get(g, "nx", 360))
    gny = gsize   !== nothing ? Int(gsize[2])       : Int(_get(g, "ny", 180))
    # vertical: [output.grid.vertical] {kind, levels}  OR legacy `heights` (= height_agl)
    vt = _get(g, "vertical", nothing)
    if vt !== nothing
        vkind  = String(_get(vt, "kind", "height_agl"))
        rawlev = filter(!_blank, _get(vt, "levels", Any[]))
        # empty levels → legacy `heights` (the GUI always writes the vertical
        # table; without this, saving would drop a heights-only grid)
        isempty(rawlev) &&
            (rawlev = filter(!_blank, _get(g, "heights", [0,100,500,1000,3000,10000])))
        levels = vkind == "height_agl" ? Float64[parse_length_m(h) for h in rawlev] :
                                         Float64[Float64(h) for h in rawlev]
    else
        vkind  = "height_agl"
        levels = Float64[parse_length_m(h) for h in filter(!_blank, _get(g, "heights", [0,100,500,1000,3000,10000]))]
    end
    s.output.grid = OutputGrid(
        crs = String(_get(g, "crs", "lonlat")), nest = Int(_get(g, "nest", 1)),
        lon0 = gx0, lat0 = gy0, nx = gnx, ny = gny, dx = gdx, dy = gdy,
        vkind = vkind, heights = levels,
    )
    s.output.quantity    = String(_get(out, "quantity", "concentration"))
    s.output.name        = String(_get(out, "name", ""))
    s.output.kernel      = let k = _get(out, "kernel", "uniform")
        k isa Bool ? (k ? "uniform" : "none") : lowercase(String(k))   # legacy Bool accepted
    end
    s.output.resolution  = Bool(_get(out, "resolution", false))
    s.output.bandwidth_floor = lowercase(String(_get(out, "bandwidth_floor", "off")))
    s.output.floor_sigmav    = Float64(_get(out, "floor_sigma_v", 0.5))
    s.output.floor_tl_s      = _durget(out, "floor_tl", 300.0)
    # [output.covariance] — flow-evolved covariance kernel parameters
    let cv = _get(out, "covariance", Dict{String,Any}())
        s.output.cov_mode      = lowercase(String(_get(cv, "mode", "physical")))
        s.output.cov_sigma0    = _lenget(cv, "sigma0", 50.0)
        s.output.cov_sigma_v   = Float64(_get(cv, "sigma_v", 0.5))
        s.output.cov_tl_s      = _durget(cv, "tl", 300.0)
        s.output.cov_sigma_max = _lenget(cv, "sigma_max", 5.0e5)
        s.output.cov_stride    = Int(_get(cv, "stride", 1))
    end
    # gridless receptors: table [output.receptors] {file, layer, [[station]]}
    # or a plain array [[output.receptors]] of inline stations
    let rc = _get(out, "receptors", nothing)
        if rc isa AbstractDict
            s.output.receptors_file  = String(_get(rc, "file", ""))
            s.output.receptors_layer = Float64[_lenget(rc, "layer_bottom", 0.0),
                                               _lenget(rc, "layer_top", 100.0)]
            s.output.receptors_hmin  = haskey(rc, "min_bandwidth") ?
                                       _lenget(rc, "min_bandwidth", -1.0) : -1.0
            rc = _get(rc, "station", nothing)
        end
        if rc isa AbstractVector
            s.output.receptors = [(id = String(_get(r, "id", "S$(k)")),
                                   lon = Float64(_get(r, "lon", 0.0)),
                                   lat = Float64(_get(r, "lat", 0.0)))
                                  for (k, r) in enumerate(rc)]
        end
    end
    s.output.per_release = Bool(_get(out, "per_release", true))
    s.output.per_species = Bool(_get(out, "per_species", true))
    s.output.group = lowercase(String(_get(out, "group", "")))
    s.output.group in ("", "auto", "species", "release", "box") ||
        error("[output].group = \"$(s.output.group)\" — use \"species\", \"release\"/\"box\", or omit (auto: species forward, release backward)")
    s.output.aod         = Bool(_get(out, "aod", false))
    # [output.reference]: string source name, or a table {source = "..."} /
    # {compare = "..."} — the [chemistry.init] source to compare the reacting
    # chemistry against (passive twin + X/X_pass/X_ref writer)
    s.output.reference   = let r = _get(out, "reference", "")
        r isa AbstractDict ? String(_get(r, "source", _get(r, "compare", ""))) : String(r)
    end
    # [output.tendencies]: list of species, or a table {species = [...]}
    s.output.tendencies  = let t = _get(out, "tendencies", String[])
        raw = t isa AbstractDict ? _get(t, "species", String[]) : t
        String[String(x) for x in raw]
    end
    s.output.path        = String(_get(out, "path", "juliaout/output/"))
    s.output.interval_s  = _durget(out, "interval", 3600.0)
    s.output.average_s   = _durget(out, "average",  s.output.interval_s)
    s.output.sample_s    = _durget(out, "sample",   900.0)
    tj = _get(out, "trajectory", Dict{String,Any}())
    s.output.trajectory.enabled     = Bool(_get(tj, "enabled", false))
    s.output.trajectory.interval_s  = _durget(tj, "interval", 0.0)
    s.output.trajectory.max_parcels = Int(_get(tj, "max_parcels", 10000))
    s.output.trajectory.species     = String[String(x) for x in _get(tj, "species", String[])]
    s.output.trajectory.format      = String(_get(tj, "format", "full"))

    # chemistry (optional [chemistry] table with nested [chemistry.cams] /
    # [chemistry.emissions])
    ch = _get(d, "chemistry", Dict{String,Any}())
    cams = _get(ch, "cams", Dict{String,Any}())
    emis = _get(ch, "emissions", Dict{String,Any}())
    light = _get(ch, "lightning", Dict{String,Any}())
    mixg = _get(ch, "mixing", Dict{String,Any}())
    hetg = _get(ch, "het", Dict{String,Any}())
    aerg = _get(ch, "aerosols", Dict{String,Any}())
    aerg isa Bool && (aerg = Dict{String,Any}("enabled" => aerg))  # aerosols = true
    # [chemistry.init] declarative block: init = "sources" +
    # [chemistry.init.sources.<name>] tables + [[chemistry.init.species]].
    # `init` may be a plain string ("cams"/"background") or the table itself.
    initg = _get(ch, "init", Dict{String,Any}())
    init_mode = initg isa AbstractString ? lowercase(String(initg)) :
                (initg isa AbstractDict && !isempty(initg)) ? "sources" : "background"
    initg isa AbstractDict || (initg = Dict{String,Any}())
    isources = InitSource[]
    for (nm, sd) in _get(initg, "sources", Dict{String,Any}())
        sd isa AbstractDict || continue
        push!(isources, InitSource(
            name  = String(nm),
            kind  = lowercase(String(_get(sd, "kind", "grib"))),
            dir   = String(_get(sd, "dir", "")),
            units = lowercase(String(_get(sd, "units", "mmr"))),
            vars  = Dict{String,String}(String(k) => String(v)
                                        for (k, v) in _get(sd, "vars", Dict{String,Any}())),
        ))
    end
    # legacy sugar: init = "cams" + [chemistry.cams] expands to the sources
    # form (init stays "cams" so existing consumers are untouched; the new
    # init layer sees a declared grib source + default)
    init_default_v = String(_get(initg, "default", ""))
    if init_mode == "cams" && isempty(isources)
        push!(isources, InitSource(name="cams", kind="grib",
                                   dir=String(_get(cams, "dir", "")), units="mmr"))
        isempty(init_default_v) && (init_default_v = "cams")
    end
    ispecies = InitSpecies[]
    for sd in _get(initg, "species", Any[])
        push!(ispecies, InitSpecies(
            name      = String(_get(sd, "name", "")),
            from      = lowercase(String(_get(sd, "from", ""))),
            predictor = String(_get(sd, "predictor", "")),
            curve     = String(_get(sd, "curve", "")),
            total     = String(_get(sd, "total", "")),
            regime    = lowercase(String(_get(sd, "regime", "any"))),
            value     = Float64(_get(sd, "value", 0.0)),
            profile   = String(_get(sd, "profile", "")),
        ))
    end
    s.chemistry = Chemistry(
        enabled     = Bool(_get(ch, "enabled", false)),
        mechanism   = String(_get(ch, "mechanism", "reducedqss")),
        timestep_s  = _durget(ch, "timestep", 3600.0),
        photolysis_timestep_s = _durget(ch, "photolysis_timestep", 0.0),
        init        = init_mode,
        init_default = init_default_v,
        init_sources = isources,
        init_species = ispecies,
        boundary    = Bool(_get(ch, "boundary", true)),
        lat_deg     = Float64(_get(ch, "lat_deg", 45.0)),
        h2o_ppb     = Float64(_get(ch, "h2o_ppb", 1.0e4)),
        h2o_prognostic = Bool(_get(ch, "h2o_prognostic", false)),
        ch4_oxidation  = Bool(_get(ch, "ch4_oxidation",  false)),
        sulfate_sad    = Bool(_get(ch, "sulfate_sad", false)),
        ice_reff_m     = _lenget(ch, "ice_reff", 1.0e-5),
        ice_shield_m   = _lenget(ch, "ice_shield", 0.0),
        sulfate_gmr_m  = _lenget(ch, "sulfate_gmr", 0.35e-6),
        sulfate_sigma_g = Float64(_get(ch, "sulfate_sigma_g", 1.3)),
        sulfate_kappa  = Float64(_get(ch, "sulfate_kappa", 0.5)),
        sulfate_so2_tau_s = _durget(ch, "sulfate_so2_tau", 17.1*86400.0),
        sulfate_so2_oh = Bool(_get(ch, "sulfate_so2_oh", false)),
        sulfate_so2_oh_feedback = Bool(_get(ch, "sulfate_so2_oh_feedback", true)),
        so4_init_kgm3  = Float64(_get(ch, "so4_init", 0.0)),
        so2_init_kgm3  = Float64(_get(ch, "so2_init", 0.0)),
        psc_dehydration = Bool(_get(ch, "psc_dehydration", false)),
        psc_rhi   = Float64(_get(ch, "psc_rhi", 1.0)),
        psc_tau_s = _durget(ch, "psc_tau", 6.0 * 3600.0),
        psc_denitrification = Bool(_get(ch, "psc_denitrification", false)),
        psc_nat_tau_s = _durget(ch, "psc_nat_tau", 5.0 * 86400.0),
        doy         = Int(_get(ch, "doy", 0)),
        # defaults must match the Chemistry struct (runspec.jl): the old loose
        # 1e-5/1e-3 TOML fallbacks silently degraded radical accuracy for any
        # config that omitted the keys
        reltol      = Float64(_get(ch, "reltol", 1.0e-6)),
        abstol      = Float64(_get(ch, "abstol", 1.0e-6)),
        precision   = lowercase(String(_get(ch, "precision", "fp64"))),
        solver      = lowercase(String(_get(ch, "solver", "rb23"))),
        photolysis  = String(_get(ch, "photolysis", "analytic")),
        lut_path    = String(_get(ch, "lut_path", "")),
        o3col_DU    = Float64(_get(ch, "o3col_DU", 300.0)),
        albedo      = Float64(_get(ch, "albedo", 0.05)),
        cloud       = Bool(_get(ch, "cloud", true)),
        cams_dir    = String(_get(cams, "dir", "")),
        cams_species= String[String(x) for x in _get(cams, "species", String[])],
        emissions_dir      = String(_get(emis, "dir", "")),
        emissions_files    = String[String(x) for x in _get(emis, "files", String[])],
        emissions_res      = Float64(_get(emis, "res", 0.4)),
        emissions_inject_m = _lenget(emis, "inject_depth", 150.0),
        emissions_nox_no   = Float64(_get(emis, "nox_no", 0.9)),
        emissions_diurnal       = String(_get(emis, "diurnal", "none")),
        emissions_diurnal_peak  = Float64(_get(emis, "diurnal_peak", 13.5)),
        emissions_diurnal_sigma = Float64(_get(emis, "diurnal_sigma", 2.0)),
        lightning_enabled    = Bool(_get(light, "enabled", false)),
        lightning_yield      = Float64(_get(light, "yield", 250.0)),
        lightning_umf_thresh = Float64(_get(light, "umf_thresh", 0.01)),
        lightning_h_min      = Float64(_get(light, "h_min", 5.0)),
        lightning_ut_lo      = _lenget(light, "ut_lo", 6000.0),
        lightning_ut_hi      = _lenget(light, "ut_hi", 14000.0),
        het_awet      = Float64(_get(hetg, "awet", 0.0)),
        mixing        = lowercase(String(_get(mixg, "scheme", "off"))),
        mixing_tau_s  = _durget(mixg, "tau", 1800.0),
        mixing_cells  = Int(_get(mixg, "cells", 1)),
        mixing_levels = Int(_get(mixg, "levels", 1)),
        aerosols      = Bool(_get(aerg, "enabled", !isempty(aerg))),
        aerosols_het  = Bool(_get(aerg, "het", true)),
        aerosols_aod  = Bool(_get(aerg, "aod", true)),
        aerosols_rh_growth = Bool(_get(aerg, "rh_growth", true)),
    )

    # task="chemistry" sugar (CHEM_INIT phase C4): the retired
    # chem_cams_pacific.jl harness always (i) domain-filled the [[release]]
    # box with air-mass parcels and (ii) compared the reacting run against its
    # init source.  The engine form is explicit — expand here so existing
    # task="chemistry" configs run unchanged through run_spec!:
    #   release.kind "release" → "domainfill";  output.reference → the
    #   default init source (when unset and one is declared).
    if s.task == "chemistry"
        for r in s.release
            r.kind == "release" && (r.kind = "domainfill")
        end
        if isempty(s.output.reference) && !isempty(s.chemistry.init_default)
            s.output.reference = s.chemistry.init_default
        end
    end
    return s
end

# ── RunSpec → dict (for writing) ──────────────────────────────────────────────
function runspec_to_dict(s::RunSpec)::Dict{String,Any}
    d = Dict{String,Any}()
    d["name"] = s.name
    d["task"] = s.task
    d["run"] = Dict{String,Any}(
        "direction" => s.direction,
        "wind" => s.wind,
        (isempty(s.checkpoint) ? () : ("checkpoint" => s.checkpoint,))...,
        (isempty(s.restart)    ? () : ("restart"    => s.restart,))...,
        (isempty(s.parcel_dump) ? () : ("parcel_dump" => s.parcel_dump,))...,
        "start" => s.start, "stop" => s.stop,
        "timestep" => format_duration(s.timestep_s),
        "benchmark" => s.benchmark,
    )
    d["meteorology"] = [begin
        md = Dict{String,Any}("model"=>m.model, "path"=>m.path, "variables"=>m.variables,
            "available"=>m.available, "crop"=>(m.crop == "explicit" ? m.crop_box : m.crop))
        m.max_met_gap_s != 6 * 3600.0 && (md["max_gap"] = format_duration(m.max_met_gap_s))
        md
    end for m in s.meteorology]
    d["species"] = map(s.species) do sp
        e = Dict{String,Any}("id"=>sp.id, "kind"=>sp.kind)
        sp.specnum != 0 && (e["specnum"] = sp.specnum)
        sp.decay_halflife_s > 0 && (e["decay_halflife"] = format_duration(sp.decay_halflife_s))
        if sp.kind == "aerosol"
            e["size"] = sp.size_distribution == "bins" ?
                Dict("distribution"=>"bins", "edges"=>sp.bin_edges_m) :
                Dict("distribution"=>"lognormal", "dg"=>sp.dg_m, "sigma_g"=>sp.sigma_g, "nbins"=>sp.nbins)
            e["density"] = sp.density
            e["aerosol_type"] = sp.aerosol_type
            sp.kappa >= 0 && (e["kappa"] = sp.kappa)
            sp.refractive_index != "auto" && (e["refractive_index"] = sp.refractive_index)
            e["settling"] = sp.settling
        end
        sp.molar_mass > 0 && (e["molar_mass"] = sp.molar_mass)   # review 2026-09-04 A7 item 3
        if sp.settle_r_eff > 0   # review 2026-09-04 A7 item 2: [species.settling] table (was dropped by the writer)
            st = Dict{String,Any}("r_eff" => sp.settle_r_eff, "rho" => sp.settle_rho)
            sp.kind == "aerosol" && (st["enabled"] = sp.settling)
            e["settling"] = st
        end
        e["dry_deposition"] = sp.dry_deposition
        e["wet_deposition"] = sp.wet_deposition
        if sp.ice
            ic = Dict{String,Any}("enabled"=>true, "rhi_crit"=>sp.ice_rhi_crit,
                                  "r_eff"=>sp.ice_r_eff_m, "h_layer"=>sp.ice_h_layer_m)
            sp.ice_q0 > 0 && (ic["q0"] = sp.ice_q0)
            e["ice"] = ic
        end
        sp.ice_scav_uptake > 0 &&
            (e["ice_scav"] = Dict{String,Any}("of"=>sp.ice_scav_of,
                                              "uptake"=>sp.ice_scav_uptake))
        if !isempty(sp.convert_to)
            cv = Dict{String,Any}("to"=>sp.convert_to,
                                  "tau"=>format_duration(sp.convert_tau_s))
            sp.convert_mass_ratio > 0 && (cv["mass_ratio"] = sp.convert_mass_ratio)
            # the thermodynamic gate and the ice closure were DROPPED here, so
            # a runspec → dict → runspec cycle silently returned a plain-τ
            # conversion (no freeze/sublimate gate, no lamina, no sectional)
            isempty(sp.convert_mode) || (cv["mode"] = sp.convert_mode)
            sp.convert_mode == "freeze" && (cv["rhi_crit"] = sp.convert_rhi_crit)
            sp.convert_lam_sigma > 0 && (cv["lamina_sigma"] = sp.convert_lam_sigma)
            sp.convert_sectional && (cv["sectional"] = true)
            e["convert"] = cv
        end
        e
    end
    d["release"] = map(s.release) do r
        e = Dict{String,Any}("species"=>r.species, "z"=>r.z, "zkind"=>r.zkind,
            "time"=>r.time, "particles"=>r.particles, "mass"=>string(r.mass_kg)*"kg")
        r.mass_per_parcel_kg > 0 &&
            (e["mass_per_parcel"] = string(r.mass_per_parcel_kg) * "kg")
        r.kind != "release" && (e["kind"] = r.kind)
        if !isempty(r.chem_init) || r.chem_h2o_plume
            ci = Dict{String,Any}(r.chem_init)
            r.chem_h2o_plume && (ci["h2o_from_plume"] = true)
            e["chem_init"] = ci
        end
        if r.kind == "volcanic"
            v = Dict{String,Any}("profile" => r.volc_profile)
            r.volc_profile == "suzuki" &&
                (v["suzuki_a"] = r.volc_suzuki_a; v["suzuki_lambda"] = r.volc_suzuki_l)
            isempty(r.volc_from_file) || (v["from_file"] = r.volc_from_file)
            if !isempty(r.volc_phases)
                v["phases"] = isempty(r.volc_pulses) ?
                    [Any[p[1], p[2], p[3], p[4], p[5]] for p in r.volc_phases] :
                    [Any[p[1], p[2], p[3], p[4], p[5], q[1], q[2], q[3]]
                     for (p, q) in zip(r.volc_phases, r.volc_pulses)]
            end
            e["volcanic"] = v
        end
        if r.plume_scheme != "none"
            p = Dict{String,Any}("scheme" => r.plume_scheme)
            r.plume_frp_w > 0 && (p["frp"] = r.plume_frp_w)
            r.plume_fire_area > 0 && (p["fire_area"] = r.plume_fire_area)
            if r.plume_scheme == "volcanic"
                p["mer"] = r.plume_mer
                p["temperature"] = r.plume_temp_k
                p["water_frac"] = r.plume_water_frac
            end
            r.plume_profile != "uniform" && (p["profile"] = r.plume_profile)
            r.plume_bottom != "pbl" && (p["bottom"] = r.plume_bottom)
            r.plume_mode != "layer" && (p["mode"] = r.plume_mode)
            r.plume_washout || (p["washout"] = false)
            isfinite(r.plume_vapor_cap) && (p["vapor_cap"] = r.plume_vapor_cap)
            r.plume_entrainment != "meanfield" && (p["entrainment"] = r.plume_entrainment)
            if r.plume_umbrella
                p["umbrella_spread"] = true
                r.plume_spread_lambda != 0.2 && (p["spread_lambda"] = r.plume_spread_lambda)
                r.plume_umbrella_rsource != 10000.0 && (p["umbrella_r_source"] = r.plume_umbrella_rsource)
                r.plume_umbrella_windshield && (p["umbrella_wind_shield"] = true)
                r.plume_spread_umin != 5.0 && (p["spread_umin"] = r.plume_spread_umin)
                r.plume_spread_tmax != 43200.0 && (p["spread_tmax"] = r.plume_spread_tmax)
                r.plume_spread_hmin > 0.0 && (p["spread_hmin"] = r.plume_spread_hmin)
                if r.plume_spread_zgain > 0.0
                    p["spread_zgain"] = r.plume_spread_zgain
                    p["spread_zgain_z0"] = r.plume_spread_zg_z0
                    p["spread_zgain_z1"] = r.plume_spread_zg_z1
                end
                if r.plume_spread_azgain_n > 0.0 || r.plume_spread_azgain_s > 0.0
                    p["spread_azgain_n"] = r.plume_spread_azgain_n
                    p["spread_azgain_s"] = r.plume_spread_azgain_s
                    if r.plume_spread_az_z1 > r.plume_spread_az_z0
                        p["spread_azgain_ramp_z0"] = r.plume_spread_az_z0
                        p["spread_azgain_ramp_z1"] = r.plume_spread_az_z1
                    end
                end
                r.plume_qc_crit != 5.0e-3 && (p["qc_crit"] = r.plume_qc_crit)
                if r.plume_umb_r_med > 0.0
                    p["umbrella_r_ice"] = r.plume_umb_r_med
                    p["umbrella_r_ice_sigma"] = r.plume_umb_r_sg
                    p["umbrella_shed_depth"] = r.plume_umb_h
                    p["umbrella_shed_time"] = r.plume_umb_t
                end
                r.plume_solid_r > 0.0 && (p["tephra_r"] = r.plume_solid_r)
            end
            if r.plume_vent_mix
                p["vent_mix"] = Dict{String,Any}(
                    "magma_t" => r.plume_magma_t,
                    "magma_water" => r.plume_magma_water,
                    "seawater_core" => r.plume_sw_core,
                    "seawater_edge" => r.plume_sw_edge,
                    "mix_exponent" => r.plume_mix_exp)
            end
            e["plume_rise"] = p
        end
        if r.kind == "seaspray"
            sp = Dict{String,Any}("salinity" => r.spray_salinity, "rho_dry" => r.spray_rho)
            r.spray_sst >= 0 && (sp["sst"] = r.spray_sst)
            r.spray_hs  >= 0 && (sp["hs"]  = r.spray_hs)
            r.spray_cd  >= 0 && (sp["cd"]  = r.spray_cd)
            r.spray_gust_sigma > 0 && (sp["gust_sigma"] = r.spray_gust_sigma)
            e["seaspray"] = sp
        end
        if !isempty(r.from_file)
            e["from_file"] = r.from_file; e["file_format"] = r.file_format
            isempty(r.lon) || (e["lon"] = r.lon; e["lat"] = r.lat)   # crop box
        else
            r.geometry == "box" || (e["geometry"] = r.geometry)
            e["lon"] = r.lon; e["lat"] = r.lat
            r.geometry in ("line","polygon") && (e["npoints"] = r.npoints)
        end
        isempty(r.comment) || (e["comment"] = r.comment)
        e
    end
    ph = s.physics
    d["physics"] = Dict{String,Any}(
        "ctl"=>ph.ctl, "ifine"=>ph.ifine,
        (ph.polar_stereo_lat != 75.0 ? ("polar_stereo_lat"=>ph.polar_stereo_lat,) : ())...,
        (ph.surface_layer_wind ? () : ("surface_layer_wind"=>false,))...,
        "turbulence"=>let td=Dict{String,Any}("scheme"=>ph.turbulence.scheme, "onthefly"=>ph.turbulence.hanna_onthefly, "sigw_ft"=>ph.turbulence.sigw_ft)
                          ph.turbulence.mesoscale && (td["mesoscale"]=true; td["fturbmeso"]=ph.turbulence.fturbmeso)
                          ph.turbulence.kh_ft > 0 && (td["kh_ft"]=ph.turbulence.kh_ft); ph.turbulence.horizontal_tke && (td["horizontal_tke"]=true); ph.turbulence.horizontal_hanna && (td["horizontal_hanna"]=true); td end,
        "convection"=>let cd2=Dict{String,Any}("scheme"=>ph.convection.scheme, "cape_min"=>ph.convection.cape_min, "cin_max"=>ph.convection.cin_max)
                          ph.convection.timestep_s > 0 && (cd2["timestep"]=format_duration(ph.convection.timestep_s))
                          ph.convection.tau_mix_s > 0 && (cd2["tau_mix"]=format_duration(ph.convection.tau_mix_s)); cd2 end,
        "vertical"=>Dict("mode"=>ph.vertical.mode, "sigma_r"=>ph.vertical.sigma_r, "sigma_t"=>ph.vertical.sigma_t, "diabatic"=>ph.vertical.diabatic, "pt_winds"=>ph.vertical.pt_winds, "pt_pattern"=>ph.vertical.pt_pattern),
        "deposition"=>Dict("dry"=>ph.deposition.dry, "wet"=>ph.deposition.wet),
        (isempty(ph.prad_species) ? () :
            ("plume_radiation"=>Dict{String,Any}(
                "species"=>ph.prad_species, "c_rad"=>ph.prad_c_rad,
                (isempty(ph.prad_partners) ? () : ("partners"=>ph.prad_partners,))...,
                (ph.prad_lw_kappa > 0 ? ("lw_kappa"=>ph.prad_lw_kappa,) : ())...,
                (ph.prad_h_ref != 1000.0 ? ("h_ref"=>ph.prad_h_ref,) : ())...,
                (ph.prad_tprime_tau > 0 ?
                    ("tprime_tau"=>format_duration(ph.prad_tprime_tau),) : ())...,
                "q_bg"=>ph.prad_q_bg, "feedback"=>ph.prad_feedback),))...,
        (ph.mixing.scheme != "off" ?
            ("mixing"=>Dict{String,Any}(
                "scheme"=>ph.mixing.scheme, "K_h"=>ph.mixing.K_h,
                "sigma0_h"=>ph.mixing.sigma0_h, "sigma_v"=>ph.mixing.sigma_v,
                "nsig"=>ph.mixing.nsig, "max_neighbors"=>ph.mixing.max_neighbors,
                "background"=>ph.mixing.background, "every"=>ph.mixing.every),) : ())...,
    )
    ph.tags_enabled && (d["tracers"] = Dict{String,Any}(
        "enabled"=>true, "strat_theta"=>ph.tags_strat_theta,
        "tropo_origin"=>ph.tags_tropo_origin))
    o = s.output
    d["output"] = Dict{String,Any}(
        "quantity"=>o.quantity, "name"=>o.name, "kernel"=>o.kernel,
        (o.resolution ? ("resolution"=>true,) : ())...,
        (o.bandwidth_floor != "off" ?
            ("bandwidth_floor"=>o.bandwidth_floor, "floor_sigma_v"=>o.floor_sigmav,
             "floor_tl"=>format_duration(o.floor_tl_s)) : ())...,
        (o.kernel == "covariance" ?
            ("covariance" => Dict{String,Any}(
                "mode"=>o.cov_mode,
                "sigma0"=>o.cov_sigma0, "sigma_v"=>o.cov_sigma_v,
                "tl"=>format_duration(o.cov_tl_s), "sigma_max"=>o.cov_sigma_max,
                "stride"=>o.cov_stride),) : ())...,
        ((!isempty(o.receptors) || !isempty(o.receptors_file)) ?
            ("receptors" => let rd = Dict{String,Any}(
                    "layer_bottom"=>o.receptors_layer[1], "layer_top"=>o.receptors_layer[2])
                o.receptors_hmin >= 0 && (rd["min_bandwidth"] = o.receptors_hmin)
                isempty(o.receptors_file) || (rd["file"] = o.receptors_file)
                isempty(o.receptors) ||
                    (rd["station"] = [Dict("id"=>r.id, "lon"=>r.lon, "lat"=>r.lat) for r in o.receptors])
                rd end,) : ())...,
        "per_release"=>o.per_release, "per_species"=>o.per_species,
        (isempty(o.group) ? () : ("group"=>o.group,))...,
        "aod"=>o.aod, "path"=>o.path,
        (isempty(o.reference) ? () : ("reference"=>o.reference,))...,
        (isempty(o.tendencies) ? () : ("tendencies"=>o.tendencies,))...,
        "interval"=>format_duration(o.interval_s), "average"=>format_duration(o.average_s),
        "sample"=>format_duration(o.sample_s),
        "grid"=>let gd=Dict{String,Any}("lon0"=>o.grid.lon0, "lat0"=>o.grid.lat0,
                    "nx"=>o.grid.nx, "ny"=>o.grid.ny, "dx"=>o.grid.dx, "dy"=>o.grid.dy)
                  o.grid.crs  != "lonlat"     && (gd["crs"] = o.grid.crs)
                  o.grid.crs  == "native"     && (gd["nest"] = o.grid.nest)
                  if o.grid.vkind == "height_agl"
                      gd["heights"] = o.grid.heights
                  else
                      gd["vertical"] = Dict("kind"=>o.grid.vkind, "levels"=>o.grid.heights)
                  end
                  gd end,
    )
    c = s.chemistry
    if c.enabled
        cd = Dict{String,Any}(
            "enabled"=>c.enabled, "mechanism"=>c.mechanism,
            "timestep"=>format_duration(c.timestep_s),
            "init"=>c.init, "boundary"=>c.boundary, "lat_deg"=>c.lat_deg,
            "h2o_ppb"=>c.h2o_ppb, "doy"=>c.doy, "reltol"=>c.reltol, "abstol"=>c.abstol, "precision"=>c.precision,
            (c.h2o_prognostic ? ("h2o_prognostic"=>true,) : ())...,
            (c.ch4_oxidation  ? ("ch4_oxidation"=>true,)  : ())...,
            (c.sulfate_sad ? ("sulfate_sad"=>true, "so4_init"=>c.so4_init_kgm3,
                              "so2_init"=>c.so2_init_kgm3) : ())...,
            (c.psc_dehydration ? ("psc_dehydration"=>true,) : ())...,
            (c.psc_denitrification ? ("psc_denitrification"=>true,) : ())...,
            (c.ice_reff_m != 1.0e-5 ? ("ice_reff"=>c.ice_reff_m,) : ())...,
            (c.ice_shield_m > 0.0   ? ("ice_shield"=>c.ice_shield_m,) : ())...,
            "solver"=>c.solver, "photolysis"=>c.photolysis,
        )
        c.photolysis_timestep_s > 0 && (cd["photolysis_timestep"]=format_duration(c.photolysis_timestep_s))
        c.photolysis == "lut" && (isempty(c.lut_path) || (cd["lut_path"] = c.lut_path))
        if c.photolysis == "lut"
            cd["o3col_DU"] = c.o3col_DU; cd["albedo"] = c.albedo
        end
        cd["cloud"] = c.cloud
        c.het_awet > 0 && (cd["het"] = Dict{String,Any}("awet"=>c.het_awet))
        c.aerosols && (cd["aerosols"] = Dict{String,Any}(
            "enabled"=>true, "het"=>c.aerosols_het, "aod"=>c.aerosols_aod,
            "rh_growth"=>c.aerosols_rh_growth))
        if c.mixing != "off"
            cd["mixing"] = Dict{String,Any}("scheme"=>c.mixing,
                "tau"=>format_duration(c.mixing_tau_s),
                "cells"=>c.mixing_cells, "levels"=>c.mixing_levels)
        end
        if c.init == "cams"
            cams = Dict{String,Any}("dir"=>c.cams_dir)
            isempty(c.cams_species) || (cams["species"] = c.cams_species)
            cd["cams"] = cams
        elseif c.init == "sources"
            ib = Dict{String,Any}()
            isempty(c.init_default) || (ib["default"] = c.init_default)
            if !isempty(c.init_sources)
                sb = Dict{String,Any}()
                for src in c.init_sources
                    e = Dict{String,Any}("kind"=>src.kind, "units"=>src.units)
                    isempty(src.dir)  || (e["dir"] = src.dir)
                    isempty(src.vars) || (e["vars"] = Dict{String,Any}(src.vars))
                    sb[src.name] = e
                end
                ib["sources"] = sb
            end
            if !isempty(c.init_species)
                ib["species"] = map(c.init_species) do sp
                    e = Dict{String,Any}("name"=>sp.name, "from"=>sp.from)
                    isempty(sp.predictor) || (e["predictor"] = sp.predictor)
                    isempty(sp.curve)     || (e["curve"] = sp.curve)
                    isempty(sp.total)     || (e["total"] = sp.total)
                    sp.regime == "any"    || (e["regime"] = sp.regime)
                    sp.value == 0.0       || (e["value"] = sp.value)
                    isempty(sp.profile)   || (e["profile"] = sp.profile)
                    e
                end
            end
            cd["init"] = ib
        end
        if !isempty(c.emissions_dir)
            em = Dict{String,Any}("dir"=>c.emissions_dir, "res"=>c.emissions_res,
                                  "inject_depth"=>c.emissions_inject_m,
                                  "nox_no"=>c.emissions_nox_no)
            if c.emissions_diurnal != "none"
                em["diurnal"] = c.emissions_diurnal
                em["diurnal_peak"] = c.emissions_diurnal_peak
                em["diurnal_sigma"] = c.emissions_diurnal_sigma
            end
            isempty(c.emissions_files) || (em["files"] = c.emissions_files)
            cd["emissions"] = em
        end
        d["chemistry"] = cd
    end
    return d
end

# Fill each species' material properties from the species DB (by `id`) wherever the
# TOML left them unset (molar mass 0, dep params < 0).  User-set values win.
function apply_species_db!(s::RunSpec)
    for sp in s.species
        db = get_species(sp.id)
        db === nothing && continue
        sp.molar_mass == 0.0     && db.molar_mass  > 0    && (sp.molar_mass = db.molar_mass)
        sp.decay_halflife_s == 0 && db.decay_halflife_s > 0 && (sp.decay_halflife_s = db.decay_halflife_s)
        sp.reldiff  < 0 && db.reldiff  > -9 && (sp.reldiff  = db.reldiff)
        sp.henry    < 0 && db.henry    > -9 && (sp.henry    = db.henry)
        sp.f0       < 0 && db.f0       > -9 && (sp.f0       = db.f0)
        sp.dryvel   < 0 && db.dryvel   > -9 && (sp.dryvel   = db.dryvel)
        sp.weta_gas < 0 && db.weta_gas > -9 && (sp.weta_gas = db.weta_gas)
        sp.wetb_gas < 0 && db.wetb_gas > -9 && (sp.wetb_gas = db.wetb_gas)
        if sp.kind == "aerosol" && db.density > 0        # named aerosol → its size/ρ
            sp.density == 2000.0 && (sp.density = db.density)
            sp.dg_m == 1.0e-6 && db.dquer > 0 && (sp.dg_m = db.dquer)
        end
    end
    return s
end

# ── Unknown-key check ─────────────────────────────────────────────────────────
# review 2026-09-04 A7 item 1: `_get(d, k, def)` is the only access path, so a
# misspelled or unknown key ([physics.mixing] Kh = 500) was silently ignored.
# The parsed TOML is flattened to dotted paths (array-of-tables → path without
# index) and compared with the schema (Field paths + schema_extra_paths).
# Schema is included AFTER this file in EasleeConfig.jl → resolve it at call time.
_schema_mod() = getfield(parentmodule(@__MODULE__), :Schema)

function _flatten_keys!(acc::Vector{String}, d::AbstractDict, prefix::String)
    for (k, v) in d
        p = isempty(prefix) ? String(k) : prefix * "." * String(k)
        if v isa AbstractDict
            isempty(v) ? push!(acc, p) : _flatten_keys!(acc, v, p)
        elseif v isa AbstractVector && !isempty(v) && all(x -> x isa AbstractDict, v)
            for e in v; _flatten_keys!(acc, e, p); end        # [[table]] → no index
        else
            push!(acc, p)
        end
    end
    return acc
end

function _levenshtein(a::AbstractString, b::AbstractString)
    ca, cb = collect(a), collect(b)
    prev = collect(0:length(cb)); cur = similar(prev)
    for i in eachindex(ca)
        cur[1] = i
        for j in eachindex(cb)
            cur[j+1] = min(prev[j+1] + 1, cur[j] + 1, prev[j] + (ca[i] == cb[j] ? 0 : 1))
        end
        prev, cur = cur, prev
    end
    return prev[end]
end

# A path is known when it is a schema path, a PREFIX of one (scalar sugar such
# as `turbulence = "hanna"`, `ice = true`, an empty table), or matches a
# wildcard pattern ("*" = any one segment) up to its own length.
function _path_known(p::String, exact::Vector{String}, patterns::Vector{Vector{SubString{String}}})
    p in exact && return true
    pd = p * "."
    any(k -> startswith(k, pd), exact) && return true
    segs = split(p, '.')
    for pat in patterns
        length(segs) <= length(pat) || continue
        all(i -> pat[i] == "*" || pat[i] == segs[i], eachindex(segs)) && return true
    end
    return false
end

"""
    unknown_config_keys(d) → Vector{Pair{String,String}}

Dotted paths present in the parsed TOML dict `d` that the schema does not
know, each paired with the nearest schema path (Levenshtein) as a suggestion.
"""
function unknown_config_keys(d::AbstractDict)
    known = String[String(k) for k in _schema_mod().schema_known_paths()]
    exact = filter(k -> !occursin("*", k), known)
    patterns = [split(k, '.') for k in known if occursin("*", k)]
    out = Pair{String,String}[]
    for p in sort!(_flatten_keys!(String[], d, ""))
        _path_known(p, exact, patterns) && continue
        # candidates: exact paths + wildcard patterns instantiated with this path's own segments
        segs = split(p, '.')
        cands = vcat(exact, [join([pat[i] == "*" ? segs[i] : pat[i] for i in eachindex(pat)], ".")
                             for pat in patterns if length(pat) == length(segs)])
        best = argmin(k -> _levenshtein(p, k), cands)
        push!(out, p => best)
    end
    return out
end

"""
    check_config_keys(d; source="")

Error (default) or warn (top-level `strict_keys = false`) listing every key of
the parsed TOML dict `d` that is not in the schema, with a suggestion each.
"""
function check_config_keys(d::AbstractDict; source::AbstractString="")
    unknown = unknown_config_keys(d)
    isempty(unknown) && return nothing
    strict = Bool(_get(d, "strict_keys", true))
    where = isempty(source) ? "" : " in $(source)"
    msg = "unknown configuration key" * (length(unknown) > 1 ? "s" : "") * where *
          " — not in the schema (easlee/src/config/schema.jl), the model would silently ignore " *
          (length(unknown) > 1 ? "them" : "it") * ":\n  - " *
          join(["$(k)  (did you mean \"$(sug)\"?)" for (k, sug) in unknown], "\n  - ")
    strict || (@warn msg; return nothing)
    error(msg * "\n  (strict_keys = false at the top level downgrades this to a warning)")
end

"""
    read_runspec(path) → RunSpec

Parse a TOML run configuration into a validated `RunSpec`: TOML → dict →
unknown-key check against the schema (`check_config_keys`; `strict_keys = false`
downgrades it to a warning) → `dict_to_runspec` (unit parsing + defaults),
species-DB fill-in of unset material properties, then `validate_runspec`
(throws listing every problem).  The entry point of every run:
`spec = read_runspec("configs/X.toml")`, then `build_run(spec)` / `run_spec!(spec)`.
"""
function read_runspec(path::AbstractString; check_paths::Bool=true)
    d = TOML.parsefile(path)
    check_config_keys(d; source=path)          # review 2026-09-04 A7 item 1
    return validate_runspec(apply_species_db!(dict_to_runspec(d)); check_paths)
end

"""
    write_runspec(spec, path) → path

Write a `RunSpec` back out as TOML (sorted keys, human-readable durations via
`format_duration`).  Inverse of `read_runspec` up to formatting/defaults.
"""
function write_runspec(s::RunSpec, path::AbstractString)
    open(path, "w") do io
        println(io, "# EASLEE run configuration")
        TOML.print(io, runspec_to_dict(s); sorted=true)
    end
    return path
end

end # module ConfigTOML
