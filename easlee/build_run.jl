
# ── build_run — RunSpec → runtime + physics objects ───────────────────────────
# Bridge layer (Main scope), included by scripts/CLI after runtime.jl and
# EasleeConfig.jl.  Strategy: emit temporary FLEXPART option files from the
# RunSpec and reuse the tested readflexpartinput + build_runtime ingestion,
# then assemble the physics objects (aerosol / hybrid vertical / convection /
# wet deposition) the new schemes need.  The user only ever sees TOML.

using Dates
using Printf
import Random

# Automatic AVAILABLE generation (scan_met_directory, write_available!).
include(joinpath(@__DIR__, "src", "io", "met_scan.jl"))

using Serialization: serialize, deserialize   # warm-start checkpoints

# Everything step_model! needs beyond (state, ctx, date), assembled from the spec.
mutable struct RunObjects
    convfields::Any        # Vector{ConvFields|nothing} or nothing
    diagconv::Any          # Vector{DiagConvFields|nothing} or nothing
    diagconvparams::Any
    wetdepparams::Any      # WetDepParams or nothing  (bulk gas scheme)
    wetdepstate::Any
    aerosol::Any           # AerosolModel or nothing
    hybridvert::Any        # HybridVert or nothing
    bench::Any             # BenchTimer or nothing
    chemstate::Any         # ChemState or nothing
    chemcfg::Any           # ChemConfig or nothing
    chembg::Any            # CamsBackground or nothing (init = "cams")
    cheminit::Any          # ChemInit ([chemistry.init] declarative layer) or nothing
    chemtl::Any            # BackgroundTimeline (bracketing slots for cheminit) or nothing
    chemcc::Any            # 3-D cloud fraction (met dict "cc") for photolysis cloud τ, or nothing
    chemqv::Any            # 3-D specific humidity (met dict "qv") for per-parcel H2O, or nothing
    spphys::Any            # SpeciesPhysics (decay / gas dry dep / wet A,B) or nothing
    icephys::Any           # IcePhysics ([species.ice] freezing/sedimentation) or nothing
    ssprays::Any           # Vector{SeaSprayBox} (kind="seaspray" sources; empty = none)
    convphys::Any          # SpeciesConvert ([species.convert]) or nothing
    pradphys::Any          # PlumeRadiation ([physics.plume_radiation]) or nothing
    tagphys::Any           # TransportTags ([tracers]) or nothing
    lightning::Any         # LightningNOx ([chemistry.lightning]) or nothing
    relchem::Any           # Dict{relbox=>Dict{String,Float64}} ([release.chem_init]) or nothing
    relh2oplume::Any       # Set{relbox} seeding cs.h2o from the plume qdet, or nothing
    aero_drydep::Bool      # review 2026-09-04 A6: aerosol dry deposition armed by the config
    aero_wetdep::Bool      # review 2026-09-04 A6: aerosol wet deposition armed by the config
end

# step_model! keyword tuple for a given RunObjects (only set what's active).
function step_kwargs(obj::RunObjects)
    kw = Dict{Symbol,Any}()
    obj.convfields  !== nothing && (kw[:convfields]  = obj.convfields)
    obj.diagconv    !== nothing && (kw[:diagconv]    = obj.diagconv;
                                    kw[:diagconvparams] = obj.diagconvparams)
    obj.wetdepparams!== nothing && (kw[:wetdepparams]= obj.wetdepparams;
                                    kw[:wetdepstate] = obj.wetdepstate)
    obj.aerosol     !== nothing && (kw[:aerosol]     = obj.aerosol;
                                    kw[:aerosol_drydep] = obj.aero_drydep;
                                    kw[:aerosol_wetdep] = obj.aero_wetdep)
    obj.hybridvert  !== nothing && (kw[:hybridvert]  = obj.hybridvert)
    obj.spphys      !== nothing && (kw[:species_physics] = obj.spphys)
    obj.icephys     !== nothing && (kw[:ice_physics] = obj.icephys)
    obj.convphys    !== nothing && (kw[:convert_physics] = obj.convphys)
    obj.pradphys    !== nothing && (kw[:plume_radiation] = obj.pradphys)
    obj.tagphys     !== nothing && (kw[:transport_tags] = obj.tagphys)
    obj.bench       !== nothing && (kw[:bench]       = obj.bench)
    return NamedTuple(kw)
end

# ── Chemistry from the [chemistry] block ──────────────────────────────────────
# Map the mechanism name → ChemMechanism singleton (chemistry/coupler.jl).
function _chem_mechanism(name::AbstractString)
    name == "reduced"    ? EASLEE.ReducedTrop() :
    name == "reducedqss" ? EASLEE.ReducedOHQSS() :
    name == "bb12"       ? EASLEE.BB12Mech() :
    name == "cb05"       ? EASLEE.CB05Mech() :
    name == "strat"      ? EASLEE.StratClOxFull() :
    name == "strat_brox" ? EASLEE.StratBrClOx() :
        error("unknown chemistry.mechanism \"$name\"")
end

# CAMS shortName → candidate chem-species names (matched against the mechanism's
# own species list).  Lets one map serve reduced / bb12 / cb05.
const _CAMS_TO_CHEMNAMES = Dict{String,Vector{String}}(
    # long-lived + NOx/Ox/HOx
    "go3"=>["O3"], "no"=>["NO"], "no2"=>["NO2"], "no3"=>["NO3"], "n2o5"=>["N2O5"],
    "hno3"=>["HNO3"], "ho2no2"=>["HNO4"], "ho2"=>["HO2"], "oh"=>["OH"],
    "h2o2"=>["H2O2"], "co"=>["CO"], "ch4_c"=>["CH4"], "h2"=>["H2"],
    # C1 / carbonyls / peroxides
    "ch3o2"=>["CH3O2"], "ch3ooh"=>["CH3OOH"], "hcho"=>["HCHO","CH2O"],
    "c2o3"=>["C2O3"], "pan"=>["PAN"], "ald2"=>["ALD2"], "ch3cocho"=>["CH3COCHO"],
    "ch3oh"=>["CH3OH"], "hcooh"=>["HCOOH"], "rooh"=>["ROOH"], "onit"=>["ORGNTR"],
    # alkane / alkene operators
    "par"=>["PAR"], "ror"=>["ROR"], "ole"=>["OLE","ISOP"], "xo2"=>["XO2"],
    "xo2n"=>["XO2N"], "rxpar"=>["RXPAR"],
    # NMVOC
    "c2h6"=>["C2H6"], "c2h5oh"=>["C2H5OH"], "c2h4"=>["C2H4"], "c3h8"=>["C3H8"],
    "ic3h7o2"=>["IC3H7O2"], "c3h6"=>["C3H6"], "hypropo2"=>["HYPROPO2"],
    "c5h8"=>["ISOP","OLE"], "ispd"=>["ISPD"], "c10h16"=>["C10H16"],
    "ch3coch3"=>["CH3COCH3"], "aco2"=>["ACO2"],
    # sulfur / ammonia
    "so2"=>["SO2"], "so4"=>["SO4"], "dms"=>["DMS"], "msa"=>["MSA"],
    "nh3"=>["NH3"], "nh2"=>["NH2"],
)

# Default variable map for the declarative init layer: mechanism species →
# background variable name.  CB05 uses the canonical CB05_FROM_CAMS bridge;
# other mechanisms derive from _CAMS_TO_CHEMNAMES (sorted shortnames, first
# mechanism-present candidate wins — deterministic).
function _default_chem_varmap(mech)
    if mech isa EASLEE.CB05Mech
        return Dict(v => k for (k, v) in sort(collect(EASLEE.CB05_FROM_CAMS)))
    end
    names = Set(EASLEE.chem_species(mech))
    vm = Dict{String,String}()
    for sn in sort!(collect(keys(_CAMS_TO_CHEMNAMES)))
        for c in _CAMS_TO_CHEMNAMES[sn]
            if c in names
                haskey(vm, c) || (vm[c] = sn)
                break
            end
        end
    end
    return vm
end

# (cams_shortname, chem_column) pairs for the species both the CAMS file and the
# mechanism have.  `bg` is a CamsBackground; `mech` the chosen mechanism.
function _cams_chem_pairs(bg, mech)
    names = EASLEE.chem_species(mech)
    col   = Dict(n => i for (i, n) in enumerate(names))
    pairs = Tuple{String,Int}[]
    for (sn, cands) in _CAMS_TO_CHEMNAMES
        EASLEE.species_index(bg, sn) == 0 && continue
        for c in cands
            if haskey(col, c); push!(pairs, (sn, col[c])); break; end
        end
    end
    return pairs
end

"""
    build_chemistry(spec, n_max) → (cs::ChemState, cfg::ChemConfig, bg, ci, tl)

Construct the chemistry state/config from `spec.chemistry`.  Rows are filled with
clean-troposphere background defaults; seeding happens per-step in the driver loop
(it needs parcel positions/time).  `bg` is the current background slot (or
`nothing`), `ci` the declarative ChemInit ([chemistry.init], incl. the legacy
init="cams" sugar), `tl` its bracketing BackgroundTimeline.  All-nothing when
disabled.
"""
function build_chemistry(spec::EasleeConfig.RunSpec, n_max::Int)
    c = spec.chemistry
    c.enabled || return (nothing, nothing, nothing, nothing, nothing)
    # This distribution ships transport + aerosol only: the reactive-chemistry
    # module (mechanisms, solvers, photolysis, heterogeneous chemistry,
    # emissions, adjoint) is not part of it.  Fail with that sentence rather
    # than with an UndefVarError from the first chemistry symbol below.
    error("[chemistry] enabled = true, but the reactive-chemistry module is not " *
          "included in this EASLEE distribution (transport + aerosol only). " *
          "Remove the [chemistry] block or set enabled = false.")
    mech = _chem_mechanism(c.mechanism)
    doy  = c.doy > 0 ? c.doy : Dates.dayofyear(spec.start)

    # Photolysis: analytic cos-SZA, or the TUV lookup table (CB05 only).
    # The analytic CB05 fallback underestimates J badly vs the LUT (J_O1D ~5.5×,
    # J_NO2 ~1.8× low at overhead sun — see PHOTOLYSIS_TUV.md §Audit): warn hard.
    if c.photolysis == "analytic" && mech isa EASLEE.CB05Mech
        msg = "chemistry: cb05 with photolysis=\"analytic\" — the analytic J values are " *
              "far below the TUV LUT (J_O1D ~5.5× low). Set photolysis=\"lut\"" *
              (c.init == "cams" ? "; analytic J invalidates CAMS comparisons." : ".")
        @warn msg
    end
    photolut = nothing
    if c.photolysis == "lut"
        mech isa EASLEE.CB05Mech ||
            @warn "chemistry.photolysis=\"lut\" only affects the cb05 mechanism; ignored for $(c.mechanism)"
        lpath = isempty(c.lut_path) ?
                joinpath(@__DIR__, "..", "TUV", "cb05_photolysis_lut.nc") : c.lut_path
        photolut = EASLEE.load_photolut(lpath)
        @info "chemistry: TUV photolysis LUT loaded" file=basename(lpath)
    end

    # GPU chemistry gate: the engine feeds per-parcel H2O/lat/lon into the
    # solver (scientifically required — never run with the scalar h2o_ppb
    # fallback).  The bb12, reduced/reducedqss, strat/strat-qss and cb05-rb23
    # GPU kernels all take the per-parcel inputs (they enter the CPU-side
    # parameter build — for cb05 the whole keff incl. LUT J/cloud-τ/khet is
    # built per parcel on CPU); mechanisms without a per-parcel-capable GPU
    # path keep the CPU force (typemax): strat_brox and cb05-fbdf.
    gth  = mech isa Union{EASLEE.BB12Mech, EASLEE.ReducedTrop, EASLEE.ReducedOHQSS,
                          EASLEE.StratClOxFull, EASLEE.StratClOxQSS} ? 256 :
           (mech isa EASLEE.CB05Mech && c.solver == "rb23") ? 256 : typemax(Int)
    c.precision in ("fp64", "fp32") ||
        error("[chemistry].precision must be \"fp64\" or \"fp32\" (got \"$(c.precision)\")")
    prec = c.precision == "fp32" ? Float32 : Float64
    prec === Float32 && !(mech isa EASLEE.BB12Mech || mech isa EASLEE.CB05Mech) &&
        @warn "[chemistry].precision=\"fp32\" only affects the bb12/cb05 fast paths; ignored for $(c.mechanism)"
    c.solver in ("rb23", "fbdf") ||
        error("[chemistry].solver must be \"rb23\" or \"fbdf\" (got \"$(c.solver)\")")
    prec === Float32 && mech isa EASLEE.CB05Mech && c.solver == "fbdf" &&
        @warn "[chemistry].precision=\"fp32\" is honoured by solver=\"rb23\" only — the CB05 FBDF path runs Float64"
    cfg  = EASLEE.ChemConfig(mech; lat_deg=c.lat_deg, doy=doy, H2O_ppb=c.h2o_ppb,
                             awet=c.het_awet,
                             h2o_prognostic=c.h2o_prognostic, ch4_oxidation=c.ch4_oxidation,
                             sulfate_sad=c.sulfate_sad, sulfate_gmr_m=c.sulfate_gmr_m,
                             sulfate_sigma_g=c.sulfate_sigma_g, sulfate_kappa=c.sulfate_kappa,
                             sulfate_so2_tau_s=c.sulfate_so2_tau_s,
                             sulfate_so2_oh=c.sulfate_so2_oh,
                             sulfate_so2_oh_feedback=c.sulfate_so2_oh_feedback,
                             so4_init=c.so4_init_kgm3, so2_init=c.so2_init_kgm3,
                             psc_dehydration=c.psc_dehydration, psc_rhi=c.psc_rhi, psc_tau_s=c.psc_tau_s,
                             psc_denitrification=c.psc_denitrification, psc_nat_tau_s=c.psc_nat_tau_s,
                             reltol=c.reltol, abstol=c.abstol, gpu_thresh=gth,
                             precision=prec, solver=c.solver,
                             ice_reff_m=c.ice_reff_m, ice_shield_m=c.ice_shield_m,
                             photolut=photolut, o3col_DU=c.o3col_DU, albedo=c.albedo)
    cs   = EASLEE.ChemState(n_max, cfg; passive = !isempty(spec.output.reference),
                            aerosols = c.aerosols,
                            tendency_species = spec.output.tendencies,
                            tropo = spec.physics.tags_tropo_origin)
    c.aerosols &&
        @info "chemistry.aerosols: carried CAMS aermr → per-parcel A_wet (het) + overhead AOD550 (photolysis)" het=c.aerosols_het

    # Background fill (every mechanism column that has a default).  For CB05 reuse
    # cb05_background; for others fill the species the background dict knows.
    if mech isa EASLEE.CB05Mech
        EASLEE.init_chemistry!(cs, EASLEE.cb05_background(), cfg)
    end

    # ── Declarative init layer (CHEM_INIT.md phases A–C) ─────────────────────
    # Any declared [chemistry.init] sources/species — including the legacy
    # init="cams" sugar, which toml_io expands to a grib source — build a
    # ChemInit + a bracketing BackgroundTimeline: initial condition and inflow
    # reseeding are then TIME-INTERPOLATED between the available slots (the
    # old engine path froze one file for the whole run).
    bg = nothing; ci = nothing; tl = nothing
    if !isempty(c.init_sources) || !isempty(c.init_species)
        bgrow = mech isa EASLEE.CB05Mech ? EASLEE.cb05_background() : nothing
        ci = EASLEE.build_chem_init(c, EASLEE.chem_species(mech);
                                    varmap=_default_chem_varmap(mech), background=bgrow)
        dname = !isempty(c.init_default) ? c.init_default :
                (length(c.init_sources) == 1 ? c.init_sources[1].name : "")
        for st in ci.steps
            st.kind === :source && st.source != dname &&
                error("chemistry.init: species \"$(st.name)\" reads source \"$(st.source)\" — " *
                      "one file-backed source per run for now (the default \"$dname\"); " *
                      "put a per-species variable override under the default source instead")
        end
        dsrc = isempty(dname) ? nothing : ci.sources[dname]
        if dsrc !== nothing && dsrc.kind in ("grib", "netcdf")
            species = dsrc.kind == "grib" && !isempty(c.cams_species) ? c.cams_species : nothing
            tl = EASLEE.build_background_timeline(Symbol(dsrc.kind), dsrc.dir;
                                                  species=species, units=dsrc.units,
                                                  at=spec.start,
                                                  window=_crop_window(spec))
            EASLEE.set_init_time!(ci, tl.bg0, tl.bg1, EASLEE.timeline_weight(tl, spec.start))
            bg = tl.bg0
            # FAIL LOUDLY rather than silently background-filling absent species:
            # if the config requested CAMS species the file does not contain, the
            # chosen mechanism cannot be initialised from CAMS as asked.
            if dsrc.kind == "grib" && !isempty(c.cams_species)
                missing = [s for s in c.cams_species if EASLEE.species_index(bg, s) == 0]
                if !isempty(missing)
                    have = sort(collect(keys(bg.species)))
                    error("""
                    chemistry: cannot run mechanism="$(c.mechanism)" with the "$dname" source as configured.
                      $(length(missing)) of $(length(c.cams_species)) requested CAMS species are NOT in
                      $(basename(tl.slots[tl.i0].path)):
                        missing : $(join(missing, ", "))
                        present : $(join(have, ", "))
                    Fix one of:
                      • regenerate the CAMS file(s) keeping all requested species
                        (scripts/split_cams.sh --species ...), or
                      • trim chemistry.cams.species to what the file actually has, or
                      • use a lighter mechanism, e.g.  mechanism = "bb12"  with
                        photolysis = "analytic"  (needs only the basic CAMS species).
                    """)
                end
            end
            @info "chemistry: [chemistry.init] source timeline" source=dname kind=dsrc.kind slots=length(tl.slots) first=tl.slots[1].valid last=tl.slots[end].valid
        end
        print(EASLEE.chem_init_summary(ci))
    end
    @info "chemistry enabled" mechanism=c.mechanism species=EASLEE.n_chem_cols(mech) init=c.init chem_timestep=c.timestep_s doy=doy
    return (cs, cfg, bg, ci, tl)
end

# Resolve the [meteorology].crop policy → (lon0, lat0, lon1, lat1) or nothing.
# "auto" crops only for a (single) domainfill release — exact there, since the
# per-step cull keeps every parcel inside the box; "off" never; an explicit
# crop_box = [lon0, lon1, lat0, lat1] is honoured on any run (warned: parcels
# reaching the crop edge are terminated — regional-domain semantics).
function _crop_window(spec; warn::Bool = false)
    cropcfg = lowercase(spec.meteorology[1].crop)
    if cropcfg == "auto"
        dfr = filter(r -> r.kind == "domainfill", spec.release)
        length(dfr) == 1 || return nothing
        r = dfr[1]
        return (Float64(r.lon[1]), Float64(minimum(r.lat)),
                Float64(r.lon[end]), Float64(maximum(r.lat)))
    elseif cropcfg == "off"
        return nothing
    else
        cb = spec.meteorology[1].crop_box
        length(cb) == 4 || error("[meteorology].crop must be \"auto\", \"off\", or [lon0, lon1, lat0, lat1]")
        warn && !any(r -> r.kind == "domainfill", spec.release) &&
            @warn "[meteorology].crop on a non-domainfill run: parcels reaching the crop edge are TERMINATED there (regional-domain semantics) — make sure the box covers everywhere your trajectories go"
        return (cb[1], cb[3], cb[2], cb[4])
    end
end

# ── FLEXPART file emission from a RunSpec ─────────────────────────────────────
_fdate(d::DateTime) = Dates.format(d, "yyyymmdd")
_ftime(d::DateTime) = Dates.format(d, "HHMMSS")
_iout(q) = q == "mixing_ratio" ? 2 : q == "both" ? 3 : 1

function _write_command(io, s)
    println(io, "&COMMAND")
    println(io, "LDIRECT= ", s.direction == "backward" ? -1 : 1, ",")
    println(io, "IBDATE= ", _fdate(s.start), ",");  println(io, "IBTIME= ", _ftime(s.start), ",")
    println(io, "IEDATE= ", _fdate(s.stop),  ",");  println(io, "IETIME= ", _ftime(s.stop),  ",")
    println(io, "LOUTSTEP= ",   round(Int, s.output.interval_s), ",")
    println(io, "LOUTAVER= ",   round(Int, s.output.average_s),  ",")
    println(io, "LOUTSAMPLE= ", round(Int, s.output.sample_s),   ",")
    println(io, "LSYNCTIME= ",  round(Int, s.timestep_s), ",")
    println(io, "CTL= ", s.physics.ctl, ",")
    println(io, "IFINE= ", s.physics.ifine, ",")
    println(io, "IOUT= ", _iout(s.output.quantity), ",")
    println(io, "LCONVECTION= ", s.physics.convection.scheme == "massflux" ? 1 : 0, ",")
    println(io, "IOUTPUTFOREACHRELEASE= ", s.output.per_release ? 1 : 0, ",")
    # [output].group → initoutgrid column grouping (auto: species fwd, release bwd)
    println(io, "OUTPUTGROUP= ", isempty(s.output.group) ? "auto" : s.output.group, ",")
    println(io, "MDOMAINFILL= 0,")
    println(io, "IND_SOURCE= 1,");  println(io, "IND_RECEPTOR= 1,")
    println(io, "IPIN= 0,");  println(io, "IFLUX= 0,");  println(io, "SURF_ONLY= 0,")
    println(io, "/")
end

# Map species ids → FLEXPART numbers (specnum if set, else a stable index).
function _species_numbers(s)
    num = Dict{String,Int}()
    next = 900
    for sp in s.species
        num[sp.id] = sp.specnum > 0 ? sp.specnum : (next += 1)
    end
    return num
end

function _write_releases(io, s, specnum)
    # Only plain box releases go to the temp file; file/inventory and other
    # geometries are expanded later in _expand_releases! (which overwrites
    # input["RELEASE"]), and domainfill parcels are placed by the engine
    # (domainfill_init!/domainfill_inject! in run_spec!).
    boxes = filter(r -> r.kind == "release" && isempty(r.from_file) && r.geometry == "box",
                   s.release)
    boxnums = [specnum[r.species] for r in boxes]
    println(io, "&RELEASES_CTRL")
    println(io, " NSPEC= ", max(1, length(unique(boxnums))), ",")
    println(io, " SPECNUM_REL= ", isempty(boxnums) ? "1" : join(boxnums, ", "), ",")
    println(io, " /")
    if isempty(boxes)
        # all releases are file-based or domainfill: emit one placeholder box
        # (1 parcel, negligible mass, centred on the first release's box) so the
        # namelist parses; _expand_releases! replaces file releases before the
        # run is built, and the domainfill population is placed by the engine.
        r0 = s.release[1]
        lonc = length(r0.lon) >= 2 ? 0.5 * (r0.lon[1] + r0.lon[2]) :
               isempty(r0.lon) ? 0.0 : r0.lon[1]
        latc = length(r0.lat) >= 2 ? 0.5 * (r0.lat[1] + r0.lat[2]) :
               isempty(r0.lat) ? 0.0 : r0.lat[1]
        println(io, "&RELEASE")
        println(io, " IDATE1= ", _fdate(r0.time[1]), ","); println(io, " ITIME1= ", _ftime(r0.time[1]), ",")
        println(io, " IDATE2= ", _fdate(r0.time[2]), ","); println(io, " ITIME2= ", _ftime(r0.time[2]), ",")
        println(io, " LON1= ", lonc, ",\n LON2= ", lonc, ",\n LAT1= ", latc, ",\n LAT2= ", latc, ",\n Z1= 0.0,\n Z2= 100.0,")
        println(io, " ZKIND= 1,\n MASS= 1.0e-30,\n PARTS= 1,\n COMMENT= \"placeholder\",\n /")
        return
    end
    for r in boxes
        println(io, "&RELEASE")
        println(io, " IDATE1= ", _fdate(r.time[1]), ",");  println(io, " ITIME1= ", _ftime(r.time[1]), ",")
        println(io, " IDATE2= ", _fdate(r.time[2]), ",");  println(io, " ITIME2= ", _ftime(r.time[2]), ",")
        println(io, " LON1= ", r.lon[1], ",");  println(io, " LON2= ", r.lon[2], ",")
        println(io, " LAT1= ", r.lat[1], ",");  println(io, " LAT2= ", r.lat[2], ",")
        println(io, " Z1= ", r.z[1], ",");      println(io, " Z2= ", r.z[2], ",")
        println(io, " ZKIND= ", r.zkind == "asl" ? 2 : 1, ",")
        println(io, " MASS= ", r.mass_kg, ",")
        println(io, " PARTS= ", r.particles, ",")
        println(io, " COMMENT= \"", isempty(r.comment) ? r.species : r.comment, "\",")
        println(io, " /")
    end
end

function _write_outgrid(io, g)
    # Heads-up when the config requests an output-grid capability whose runtime
    # binning is not wired yet (config is parsed/validated, but conc.jl/initoutgrid
    # still bin by height-AGL and infer the horizontal CRS from the met projection).
    g.vkind == "height_agl" ||
        @warn "output.grid.vertical.kind=\"$(g.vkind)\" is parsed but the runtime still bins by height_agl (TODO)"
    # review 2026-09-04 A7 item 6: crs/nest are parsed and validated but the runtime ignores them
    g.crs == "lonlat" ||
        @warn "output.grid.crs=\"$(g.crs)\" is parsed but IGNORED: the runtime output grid is always lon/lat (lon0/lat0/dx/dy in degrees)"
    g.nest == 1 ||
        @warn "output.grid.nest=$(g.nest) is parsed but IGNORED: the output grid is not tied to a nest projection"
    println(io, "&OUTGRID")
    println(io, " OUTLON0= ", g.lon0, ",");  println(io, " OUTLAT0= ", g.lat0, ",")
    println(io, " NUMXGRID= ", g.nx, ",");   println(io, " NUMYGRID= ", g.ny, ",")
    println(io, " DXOUT= ", g.dx, ",");      println(io, " DYOUT= ", g.dy, ",")
    println(io, " OUTHEIGHTS= ", join(g.heights, ","), ",")
    # forward-compat keys (ignored by the current initoutgrid; for the staged
    # native/lonlat-override + vertical-kind binning)
    println(io, " CRS= ", g.crs, ",")
    println(io, " OUTNEST= ", g.nest, ",")
    println(io, " OUTVKIND= ", g.vkind, ",")
    println(io, " /")
end

# Parse a loose-format locations / inventory file → Vector{(lat, lon, lev, weight)}.
# Accepts comma- or whitespace-separated columns; skips blank and comment
# (# or //) lines and any header line whose first field isn't numeric.
function _parse_release_file(path::AbstractString, fmt::AbstractString)
    locs = Tuple{Float64,Float64,Float64,Float64}[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#") || startswith(line, "//")) && continue
        toks = split(line, r"[,\s]+"; keepempty=false)
        length(toks) >= 3 || continue
        nums = tryparse.(Float64, toks)
        (nums[1] === nothing || nums[2] === nothing || nums[3] === nothing) && continue   # header
        if fmt == "lonlatlev"
            lon, lat, lev = nums[1], nums[2], nums[3]
        else                       # latlonlev (default) and inventory share lat,lon first
            lat, lon, lev = nums[1], nums[2], nums[3]
        end
        # inventory: 3rd column is the EMISSION value at a single level; treat
        # level as 4th col if present, weight = emission.
        local w::Float64, z::Float64
        if fmt == "inventory"
            z = (length(nums) >= 4 && nums[4] !== nothing) ? nums[4] : 0.0
            w = lev                # the value column
        else
            z = lev
            w = (length(nums) >= 4 && nums[4] !== nothing) ? nums[4] : 1.0
        end
        (isfinite(w) && w > 0) || (w = 1.0)
        push!(locs, (lat, lon, z, w))
    end
    return locs
end

# Gridded NetCDF emission inventory (CAMS-GLOB-* style: lat × lon × time, one
# variable per sector in kg m⁻² s⁻¹) → point sources.  The sector fluxes are
# summed at the time slice matching the release month, area-weighted to an
# emission RATE [kg/s] per cell, cropped to the release lon/lat box when given,
# and thresholded to the dominant cells (cumulative 99.9 % of the total, capped
# at `maxpoints`) so a 6.5-million-cell global grid yields a tractable number
# of point sources.
function _parse_inventory_nc(path::AbstractString, r; maxpoints::Int=50_000)
    v = NetCDF.open(path)
    lat = Float64.(NetCDF.ncread(path, "lat"))
    lon = Float64.(NetCDF.ncread(path, "lon"))
    nlon, nlat = length(lon), length(lat)

    # time slice: months-since axis → match the release start month
    tidx = 1
    if haskey(v.vars, "time")
        tv = v.vars["time"]; u = get(tv.atts, "units", "")
        m = match(r"months since (\d{4})-(\d{2})", u)
        if m !== nothing
            base_y, base_m = parse(Int, m.captures[1]), parse(Int, m.captures[2])
            want = (year(r.time[1]) - base_y) * 12 + (month(r.time[1]) - base_m)
            tvals = Int.(tv[:])
            j = findfirst(==(want), tvals)
            tidx = j === nothing ? clamp(want + 1, 1, length(tvals)) : j
        end
    end

    # sum the sector flux variables (3-D, kg m-2 s-1) at the time slice.
    # CAMS-GLOB-ANT files carry the per-sector vars AND their `sum` — use `sum`
    # alone when present (adding everything double-counts), and skip `x_*`
    # sub-sectors whose parent `x` is also a variable (fef_coal/gas/oil vs fef).
    names3 = [name for (name, var) in v.vars if ndims(var) == 3 &&
              occursin("m-2 s-1", lowercase(get(var.atts, "units", "")))]
    use = "sum" in names3 ? ["sum"] :
          [n for n in names3 if !(occursin("_", n) && String(split(n, "_")[1]) in names3)]
    flux = zeros(Float64, nlon, nlat)
    nsect = 0
    for name in use
        sl = v.vars[name][:, :, tidx]              # (lon, lat)
        @inbounds for j in 1:nlat, i in 1:nlon
            f = Float64(sl[i, j]); isfinite(f) && f > 0 && (flux[i, j] += f)
        end
        nsect += 1
    end
    nsect == 0 && error("inventory $(basename(path)): no kg m-2 s-1 sector variables found")

    # optional PER-CELL injection height "inject_z" [m] (BB plume rise, e.g. GFAS
    # mami/apt): overrides the config z box per cell when present.
    injz = nothing
    if haskey(v.vars, "inject_z")
        iz = v.vars["inject_z"]
        injz = ndims(iz) == 3 ? iz[:, :, tidx] : iz[:, :]     # (lon, lat)
    end

    # injection height: mid-point of the release z box (zkind honoured at add!).
    # Default release z = [0, 100] ⇒ 50 m ≈ surface for anthropogenic inventories.
    zrel = length(r.z) >= 2 ? 0.5 * (Float64(r.z[1]) + Float64(r.z[2])) :
           (length(r.z) >= 1 ? Float64(r.z[1]) : 0.0)

    # optional crop to the release lon/lat box
    loncrop = length(r.lon) == 2 ? extrema(r.lon) : (-Inf, Inf)
    latcrop = length(r.lat) == 2 ? extrema(r.lat) : (-Inf, Inf)
    dlon = abs(lon[2] - lon[1]); dlat = abs(lat[2] - lat[1])
    R = 6.371e6; deg = π / 180

    cells = Tuple{Float64,Float64,Float64,Float64}[]  # (lat, lon, emission_rate kg/s, z_m)
    @inbounds for j in 1:nlat
        (latcrop[1] <= lat[j] <= latcrop[2]) || continue
        area = (dlon * deg * R * cosd(lat[j])) * (dlat * deg * R)
        for i in 1:nlon
            flux[i, j] > 0 || continue
            (loncrop[1] <= lon[i] <= loncrop[2]) || continue
            zc = injz === nothing ? zrel :
                 (isfinite(injz[i, j]) && injz[i, j] > 0 ? Float64(injz[i, j]) : zrel)
            push!(cells, (lat[j], lon[i], flux[i, j] * area, zc))
        end
    end
    isempty(cells) && error("inventory $(basename(path)): no emission in the release region")

    sort!(cells; by = c -> -c[3])                  # strongest first
    total = sum(c[3] for c in cells)
    keep = Tuple{Float64,Float64,Float64,Float64}[]
    cum = 0.0
    for c in cells
        push!(keep, (c[1], c[2], c[4], c[3]))      # (lat, lon, per-cell z, weight)
        cum += c[3]
        (length(keep) >= maxpoints || cum >= 0.999 * total) && break
    end
    zinfo = injz === nothing ? "z=$(round(zrel))m (config)" :
            "per-cell z $(round(minimum(k[3] for k in keep)))–$(round(maximum(k[3] for k in keep)))m"
    @info "inventory parsed" file=basename(path) sectors=nsect month=tidx cells=length(cells) kept=length(keep) total_kg_per_s=round(total, sigdigits=4) inject=zinfo
    return keep
end

# ── Volcanic column source (kind = "volcanic") ────────────────────────────────
# Normalised vertical mass weights over `n` sublayers of [zb, zt]:
#   suzuki   — Suzuki (1983): k(z) ∝ [(1 − z/H)·exp(A·(z/H − 1))]^λ with H = zt
#              (the standard eruption-column shape; A≈4 puts the mass near the
#              umbrella, larger A → top-heavier)
#   umbrella — Gaussian centred at zb + 0.85·(zt − zb), σ = (zt − zb)/8
#              (satellite-style umbrella cloud with a thin stem)
#   uniform  — equal mass per sublayer
function _volc_profile_weights(profile::AbstractString, zb::Real, zt::Real,
                               n::Int, A::Real, λ::Real)
    zb = Float64(zb); zt = Float64(zt); A = Float64(A); λ = Float64(λ)
    zmid = [zb + (l - 0.5) * (zt - zb) / n for l in 1:n]
    w = if profile == "uniform"
        ones(n)
    elseif profile == "umbrella"
        zc = zb + 0.85 * (zt - zb); σ = (zt - zb) / 8
        [exp(-((z - zc) / σ)^2) for z in zmid]
    else                                        # suzuki
        H = zt
        [max((1 - z / H) * exp(A * (z / H - 1)), 0.0)^λ for z in zmid]
    end
    s = sum(w)
    s > 0 || error("volcanic: profile \"$profile\" has zero weight over [$zb, $zt] m")
    return w ./ s
end

# CSV injection table: t0, t1, z_bottom [m], z_top [m], frac  (comments/# ok).
function _parse_volcanic_csv(path::AbstractString)
    phases = Tuple{DateTime,DateTime,Float64,Float64,Float64}[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        f = split(line, r"[,\s]+"; keepempty=false)
        length(f) >= 5 || continue
        push!(phases, (DateTime(f[1]), DateTime(f[2]),
                       parse(Float64, f[3]), parse(Float64, f[4]), parse(Float64, f[5])))
    end
    isempty(phases) && error("volcanic from_file has no usable phases: $path")
    return phases
end

# ── Geometry sampling → point locations (lat, lon) ────────────────────────────
# Sample npoints along a polyline (lon/lat vertices), evenly by arc length.
function _sample_line(lon, lat, npoints)
    seg = [sqrt((lon[i+1]-lon[i])^2 + (lat[i+1]-lat[i])^2) for i in 1:length(lon)-1]
    total = sum(seg); total <= 0 && return [(lat[1], lon[1])]
    cum = cumsum(seg); pts = Tuple{Float64,Float64}[]
    for k in 0:npoints-1
        d = total * k / max(npoints - 1, 1)
        i = findfirst(c -> c >= d, cum); i === nothing && (i = length(seg))
        d0 = i == 1 ? 0.0 : cum[i-1]
        f = seg[i] > 0 ? (d - d0) / seg[i] : 0.0
        push!(pts, (lat[i] + f*(lat[i+1]-lat[i]), lon[i] + f*(lon[i+1]-lon[i])))
    end
    return pts
end

# Ray-casting point-in-polygon test.
function _in_polygon(x, y, vx, vy)
    inside = false; n = length(vx); j = n
    @inbounds for i in 1:n
        if ((vy[i] > y) != (vy[j] > y)) &&
           (x < (vx[j]-vx[i]) * (y-vy[i]) / (vy[j]-vy[i]) + vx[i])
            inside = !inside
        end
        j = i
    end
    return inside
end

# Sample npoints uniformly inside a polygon (rejection within the bbox).
function _sample_polygon(lon, lat, npoints)
    lo1, lo2 = extrema(lon); la1, la2 = extrema(lat)
    pts = Tuple{Float64,Float64}[]; tries = 0
    while length(pts) < npoints && tries < 200 * npoints
        x = lo1 + rand() * (lo2 - lo1); y = la1 + rand() * (la2 - la1)
        _in_polygon(x, y, lon, lat) && push!(pts, (y, x))
        tries += 1
    end
    isempty(pts) && push!(pts, (sum(lat)/length(lat), sum(lon)/length(lon)))  # degenerate
    return pts
end

# Rebuild input["RELEASE"] + SPECNUM_REL from the spec: box releases stay one
# entry; point/line/polygon/profile and file releases expand to point sources,
# distributing the release's total particles and mass.
function _expand_releases!(input, spec, specnum)
    K = ("LON1","LON2","LAT1","LAT2","Z1","Z2","ZKIND","MASS","PARTS",
         "IDATE1","ITIME1","IDATE2","ITIME2","COMMENT")
    R = Dict{String,Vector{String}}(k => String[] for k in K)
    specnums = Int[]
    bspecies = String[]                      # species id per expanded box (relbox order)
    bspecidx = Int[]                         # spec.release index per expanded box
    bphase   = Int[]                         # volcanic phase index per box (0 = n/a)
    _cur_species = Ref("")                   # set per release before its add! calls
    _cur_idx = Ref(0)
    _cur_phase = Ref(0)
    add!(lon1,lon2,lat1,lat2,z1,z2,zk,mass,parts,t1,t2,com,snum) = begin
        push!(R["LON1"], string(lon1)); push!(R["LON2"], string(lon2))
        push!(R["LAT1"], string(lat1)); push!(R["LAT2"], string(lat2))
        push!(R["Z1"], string(z1));     push!(R["Z2"], string(z2))
        push!(R["ZKIND"], string(zk));  push!(R["MASS"], string(mass)); push!(R["PARTS"], string(parts))
        push!(R["IDATE1"], _fdate(t1)); push!(R["ITIME1"], _ftime(t1))
        push!(R["IDATE2"], _fdate(t2)); push!(R["ITIME2"], _ftime(t2))
        push!(R["COMMENT"], com); push!(specnums, snum)
        push!(bspecies, _cur_species[]); push!(bspecidx, _cur_idx[])
        push!(bphase, _cur_phase[])
    end
    for (kspec, r) in enumerate(spec.release)
        r.kind == "domainfill" && continue   # engine-placed (domainfill_init!/_inject!)
        zk = r.zkind == "asl" ? 2 : 1
        snum = specnum[r.species]
        _cur_species[] = r.species
        _cur_idx[] = kspec
        if r.kind == "volcanic"
            # Time-resolved eruption column at the vent: each phase becomes NL
            # Suzuki/umbrella-weighted sublayer releases over its own window.
            phases = isempty(r.volc_from_file) ? r.volc_phases :
                     _parse_volcanic_csv(r.volc_from_file)
            fsum = sum(p[5] for p in phases)
            NL = 24
            kept_w = 0.0; nsub = 0
            plan = Tuple{Float64,Float64,DateTime,DateTime,Float64,Int,Int}[]  # z1,z2,t0,t1,wgt,parts,phase
            for (iph, (t0v, t1v, zb, zt, frac)) in enumerate(phases)
                w = _volc_profile_weights(r.volc_profile, zb, zt, NL,
                                          r.volc_suzuki_a, r.volc_suzuki_l)
                for l in 1:NL
                    wgt = (frac / fsum) * w[l]
                    parts = round(Int, r.particles * wgt)
                    parts >= 1 || continue          # unresolved tail (mass renormalised below)
                    z1 = zb + (l - 1) * (zt - zb) / NL
                    z2 = zb + l * (zt - zb) / NL
                    push!(plan, (z1, z2, t0v, t1v, wgt, parts, iph))
                    kept_w += wgt; nsub += 1
                end
            end
            isempty(plan) && error("release \"$(r.species)\": volcanic column unresolved — increase particles")
            for (z1, z2, t0v, t1v, wgt, parts) in plan
                add!(r.lon[1], r.lon[1], r.lat[1], r.lat[1], z1, z2, zk,
                     r.mass_kg * wgt / kept_w, parts, t0v, t1v,
                     isempty(r.comment) ? "volcanic $(r.species)" : r.comment, snum)
            end
            @info "expanded volcanic release" species=r.species phases=length(phases) sublayers=nsub profile=r.volc_profile mass_kept_pct=round(100kept_w, digits=1)
        elseif r.kind == "seaspray"
            # met-driven OSSA source: the box is reserved for species/relbox
            # attribution only — the engine spawns the parcels from the flux
            # (build_seaspray/apply_seaspray! in run_spec!): 0 scheduled
            # particles, 0 scheduled mass.
            add!(r.lon[1], r.lon[2], r.lat[1], r.lat[2], r.z[1], r.z[2], zk,
                 0.0, 0, r.time[1], r.time[2],
                 isempty(r.comment) ? "seaspray $(r.species)" : r.comment, snum)
        elseif isempty(r.from_file) && r.geometry == "box"
            add!(r.lon[1], r.lon[2], r.lat[1], r.lat[2], r.z[1], r.z[2],
                 zk, r.mass_kg, r.particles, r.time[1], r.time[2],
                 isempty(r.comment) ? r.species : r.comment, snum)
        elseif isempty(r.from_file)
            # geometry → point sources, sharing the total particles/mass equally
            pts = r.geometry == "point"   ? [(r.lat[1], r.lon[1])] :
                  r.geometry == "profile" ? [(r.lat[1], r.lon[1])] :   # levels handled below
                  r.geometry == "line"    ? _sample_line(r.lon, r.lat, r.npoints) :
                                            _sample_polygon(r.lon, r.lat, r.npoints)
            if r.geometry == "profile"
                # one point at each z level
                n = length(r.z); pp = max(1, r.particles ÷ n)
                for z in r.z
                    add!(r.lon[1], r.lon[1], r.lat[1], r.lat[1], z, z, zk,
                         r.mass_kg / n, pp, r.time[1], r.time[2], r.species, snum)
                end
            else
                z1 = r.z[1]; z2 = length(r.z) >= 2 ? r.z[2] : r.z[1]
                n = length(pts); pp = max(1, r.particles ÷ n)
                for (lat, lon) in pts
                    add!(lon, lon, lat, lat, z1, z2, zk, r.mass_kg / n, pp,
                         r.time[1], r.time[2], r.species, snum)
                end
            end
            @info "expanded $(r.geometry) release" species=r.species points=(r.geometry=="profile" ? length(r.z) : length(pts))
        else
            # gridded NetCDF inventory vs loose-format text file
            locs = (r.file_format == "inventory" && occursin(r"\.nc4?$"i, r.from_file)) ?
                   _parse_inventory_nc(r.from_file, r) :
                   _parse_release_file(r.from_file, r.file_format)
            isempty(locs) && error("release file has no usable locations: $(r.from_file)")
            sumw = sum(l[4] for l in locs)
            # Distribute the total particles by weight.  With more locations than
            # particles the sub-particle tail rounds to zero and is dropped; the
            # released MASS is then renormalised over the kept (resolvable) cells
            # so the inventory's total emission is conserved, just concentrated
            # in the dominant cells.
            kept = Tuple{Float64,Float64,Float64,Float64,Int}[]   # lat,lon,z,frac,parts
            for (lat, lon, z, w) in locs
                frac = w / sumw
                parts = round(Int, r.particles * frac)
                parts >= 1 && push!(kept, (lat, lon, z, frac, parts))
            end
            isempty(kept) && error("release \"$(r.species)\": particles=$(r.particles) too few " *
                                   "for $(length(locs)) locations — increase particles")
            fkept = sum(k[4] for k in kept)
            for (lat, lon, z, frac, parts) in kept
                add!(lon, lon, lat, lat, z, z, zk, r.mass_kg * frac / fkept, parts,
                     r.time[1], r.time[2], r.species, snum)
            end
            @info "expanded file release" file=basename(r.from_file) locations=length(locs) kept=length(kept) species=r.species particles=r.particles emission_kept=round(100*fkept, digits=1)
        end
    end
    # [species.convert]: every box whose species converts gets a SHADOW box
    # for the daughter species (0 scheduled particles/mass) — converted parcels
    # are re-boxed there by apply_convert!, so species physics / aerosol /
    # output attribution follow the daughter from the moment of conversion.
    byid_cv = Dict(sp.id => sp for sp in spec.species)
    bconvert = zeros(Int, length(bspecies))
    # CHAINS.  A shadow box may itself convert (A -> B -> C, e.g.
    # vapour -> ice -> tagged-sublimated-water), so newly appended boxes MUST be
    # revisited.  Iterating over a FIXED 1:nb0 silently dropped every link after
    # the first: the config parsed, the first conversion armed, and the rest of
    # the chain was never wired (diagnosed 2026-08-01 on the HT water-tagging
    # run — `sublimate_boxes = 0` while `freeze_boxes = 42`, and the ice, unable
    # to return, simply sedimented 36 Tg of water out of the stratosphere).
    # `lineage` carries the species already visited on the way to a box, so a
    # cycle (A -> B -> A) closes onto the existing box instead of regressing
    # forever; that reproduces the old hard-coded 2-cycle exactly.
    lineage = [Set{String}([bspecies[b]]) for b in eachindex(bspecies)]
    b = 0
    while b < length(bspecies)
        b += 1
        sp = get(byid_cv, bspecies[b], nothing)
        (sp === nothing || isempty(sp.convert_to)) && continue
        # already wired (cycle closed onto an existing box)?
        bconvert[b] != 0 && continue
        # would this link re-enter a species already on this chain? then it is a
        # cycle: point back at the ancestor box carrying that species.
        if sp.convert_to in lineage[b]
            anc = findfirst(a -> bspecies[a] == sp.convert_to &&
                            issubset(lineage[a], lineage[b]), eachindex(bspecies))
            anc === nothing || (bconvert[b] = anc)
            bconvert[b] != 0 && continue
        end
        _cur_species[] = sp.convert_to
        _cur_idx[]     = bspecidx[b]
        _cur_phase[]   = bphase[b]
        add!(R["LON1"][b], R["LON2"][b], R["LAT1"][b], R["LAT2"][b],
             R["Z1"][b], R["Z2"][b], R["ZKIND"][b], 0.0, 0,
             DateTime(R["IDATE1"][b] * R["ITIME1"][b], dateformat"yyyymmddHHMMSS"),
             DateTime(R["IDATE2"][b] * R["ITIME2"][b], dateformat"yyyymmddHHMMSS"),
             "convert " * bspecies[b] * "->" * sp.convert_to,
             get(specnum, sp.convert_to, 0))
        push!(bconvert, 0)                      # the shadow may convert onward
        push!(lineage, union(lineage[b], Set([sp.convert_to])))
        bconvert[b] = length(bspecies)          # parent box → its shadow index
        # REVERSIBLE pair (vapour ⇄ ice): a phase change must be able to come
        # BACK.  When the daughter declares a conversion to the parent species,
        # close the 2-cycle immediately (identical to the pre-2026-08-01
        # behaviour); longer chains are handled by the loop revisiting the
        # appended box.
        let dsp = get(byid_cv, sp.convert_to, nothing)
            if dsp !== nothing && dsp.convert_to == bspecies[b]
                bconvert[end] = b
            end
        end
    end
    input["_box_convert"] = bconvert
    input["RELEASE"] = R
    haskey(input, "RELEASES_CTRL") || (input["RELEASES_CTRL"] = Dict{String,Any}())
    input["RELEASES_CTRL"]["SPECNUM_REL"] = string.(specnums)
    input["nbreleasebox"] = length(specnums)
    input["_box_species"] = bspecies         # relbox → species id (SpeciesPhysics)
    input["_box_specidx"] = bspecidx         # relbox → spec.release index (plume cfg)
    input["_box_phase"]   = bphase           # relbox → volcanic phase (0 = n/a)
    return nothing
end

# Resolve the AVAILABLE file:
#   • explicit path → use it;
#   • "auto" → SCAN the met directory (header valid time, filename fallback)
#     and write a fresh AVAILABLE into the work dir — the user never hand-writes
#     one.  An existing AVAILABLE in the met dir is ignored on "auto" so the
#     listing always matches what is actually on disk.
function _resolve_available(m::EasleeConfig.MetSource, workdir::AbstractString, idx::Int)
    if m.available != "auto" && isfile(m.available)
        return m.available
    end
    varsfile = _variables_path(m)
    multistep = _is_multistep(varsfile)
    include_only, exclude = _scan_filters(varsfile)
    @info "scanning meteorology directory for $(m.model) (auto AVAILABLE)…" path=m.path
    entries = scan_met_directory(m.path; model=m.model, multistep=multistep,
                                 include_only=include_only, exclude=exclude)
    isempty(entries) &&
        error("meteorology source \"$(m.model)\": no readable met files in $(m.path)")
    avail = joinpath(workdir, "AVAILABLE_$(idx)")
    write_available!(avail, entries, varsfile, m.path)
    @info "  found $(length(entries)) time steps: $(first(entries)[1]) → $(last(entries)[1])"
    return avail
end

# Guard the run window against the met archive.  Only the START side is fatal:
# a start before the first met time makes the time-interpolation index 0 (an
# out-of-bounds read deep in the advance kernel — the bug this catches).  A stop
# AFTER the last met time is benign — readmetfiles! holds the last frame, so the
# final steps run on stale met; that is the historical behaviour (e.g. the CAPTEX
# validation, whose releases legitimately overrun a gapped archive), so WARN only
# rather than aborting the run.  Forecast-step archives (base + steps 1…N) often
# start their day at base+1h, not 00:00.
function _check_met_window(avail::AbstractString, start::DateTime, stop::DateTime, model;
                           max_gap_s::Float64 = 6 * 3600.0)
    ts = DateTime[]
    for line in eachline(avail)
        m = match(r"^(\d{8}) (\d{6})\s", line)
        m === nothing && continue
        push!(ts, DateTime(m.captures[1] * m.captures[2], "yyyymmddHHMMSS"))
    end
    isempty(ts) && return
    sort!(ts)
    t0, t1 = ts[1], ts[end]
    lo, hi = min(start, stop), max(start, stop)   # backward runs swap
    lo < t0 && error("run.start $lo is before the first $model met time $t0 — " *
                     "move run.start to ≥ $t0 (this archive's day begins there)")
    gaph = max_gap_s / 3600.0
    _dh(a, b) = Dates.value(b - a) / 3.6e6        # DateTime diff → hours
    # trailing overrun: run extends past the archive.  A short hold (≤ max_gap)
    # is the intended forecast-tail behaviour; a big one means a missing file.
    if hi > t1
        over = _dh(t1, hi)
        max_gap_s > 0 && over > gaph &&
            error("run.stop $hi is $(round(over, digits=1)) h past the last $model met " *
                  "time $t1 — a gap > $(round(gaph, digits=1)) h means a MISSING met file, " *
                  "not an intended overrun; add the file(s) or set [meteorology].max_gap to override")
        @warn "run.stop $hi is after the last $model met time $t1 — the final " *
              "steps will run on the last available met (held constant)" maxlog=3
    end
    # internal gaps: consecutive met times, whose interval overlaps [lo,hi], more
    # than max_gap apart → the run would interpolate/hold across missing met.
    if max_gap_s > 0
        for i in 1:length(ts)-1
            a, b = ts[i], ts[i+1]
            (a < hi && b > lo) || continue        # pair straddles the run window
            dh = _dh(a, b)
            dh > gaph &&
                error("$model met has a $(round(dh, digits=1)) h gap between $a and $b " *
                      "inside the run window [$lo … $hi] — a MISSING met file (the run " *
                      "would interpolate across it); add the file(s) or set " *
                      "[meteorology].max_gap to override")
        end
    end
    return
end

# Resolve the variables namelist path for a met source (listvariables/<name>.vars).
function _variables_path(m::EasleeConfig.MetSource)
    isempty(m.variables) && error("meteorology \"$(m.model)\": set `variables` to a namelist name")
    isfile(m.variables) && return m.variables                      # absolute/explicit
    cand = joinpath(@__DIR__, "..", "listvariables", m.variables * ".vars")
    isfile(cand) && return abspath(cand)
    error("meteorology \"$(m.model)\": variables namelist not found: $(m.variables)")
end

# Does the variables namelist declare multistep = true (several times per file)?
function _is_multistep(varsfile::AbstractString)
    isfile(varsfile) || return false
    for line in eachline(varsfile)
        occursin(r"^\s*multistep\s*=\s*true"i, line) && return true
    end
    return false
end

# Scan filters from the namelist: include only files matching the companion
# SOURCE token (the main files), and always exclude the isentropic PT winds.
function _scan_filters(varsfile::AbstractString)
    include_only = Regex[]; exclude = Regex[r"IFS_PT_MET", r"FULL_PT"]
    isfile(varsfile) || return include_only, exclude
    esc(t) = Regex(replace(String(t), r"([.\-\[\]()+*?^$\\])" => s"\\\1"))
    for line in eachline(varsfile)
        m = match(r"^\s*companions\s*=\s*(.+?)\s*(?:#.*)?$"i, line)
        m === nothing && continue
        for pat in split(m.captures[1], ",")
            parts = split(strip(pat), ">")
            length(parts) == 2 && push!(include_only, esc(strip(parts[1])))
        end
    end
    return include_only, exclude
end

# Wind dimensionality → which parcel indices the advance may update.
# (update_indx, update_indy, update_indl).  A frozen index keeps each parcel's
# own value, so 2D/1D runs need no met modification — see updatepos!.
function _wind_mask(wind::AbstractString)::NTuple{3,Bool}
    w = lowercase(replace(wind, "_"=>"-"))
    w in ("3d", "")        && return (true,  true,  true)
    w == "2d-lonlat"       && return (true,  true,  false)  # freeze level
    w == "2d-latalt"       && return (false, true,  true)   # freeze lon
    w == "2d-lonalt"       && return (true,  false, true)   # freeze lat
    w in ("1d-vertical","1d","column") && return (false, false, true)  # column only
    @warn "unknown run.wind=\"$wind\" — using 3D"; return (true, true, true)
end

"""
    build_run(spec; prefetch=true, extra_slots=0, workdir=mktempdir())
        → (rt, obj::RunObjects)

Build the runtime and the physics objects for a TOML/RunSpec run.  Emits
temporary FLEXPART option files into `workdir` and ingests them through the
tested `readflexpartinput` + `build_runtime` path — an internal mechanism, the
user only ever sees TOML — then assembles the aerosol / hybrid-vertical /
convection / wet-deposition / chemistry objects.  `rt` is the runtime tuple
from `build_runtime`; pass `step_kwargs(obj)` on to `step_model!`.
"""
function build_run(spec::EasleeConfig.RunSpec; prefetch::Bool=true,
                   extra_slots::Int=0, workdir::AbstractString=mktempdir())
    EasleeConfig.validate_runspec(spec)
    optdir = joinpath(workdir, "options"); mkpath(optdir)
    mkpath(spec.output.path)

    specnum = _species_numbers(spec)
    open(io -> _write_command(io, spec),             joinpath(optdir, "COMMAND"), "w")
    open(io -> _write_releases(io, spec, specnum),   joinpath(optdir, "RELEASES"), "w")
    open(io -> _write_outgrid(io, spec.output.grid), joinpath(optdir, "OUTGRID"), "w")

    # pathnames (one nest block per meteorology source)
    pnfile = joinpath(workdir, "pathnames")
    open(pnfile, "w") do io
        for (k, m) in enumerate(spec.meteorology)
            avail = _resolve_available(m, workdir, k)
            _check_met_window(avail, spec.start, spec.stop, m.model; max_gap_s = m.max_met_gap_s)
            if k == 1
                println(io, optdir, "/")
                println(io, endswith(spec.output.path, "/") ? spec.output.path : spec.output.path * "/")
                println(io, endswith(m.path, "/") ? m.path : m.path * "/")
                println(io, avail)
            else
                println(io, endswith(m.path, "/") ? m.path : m.path * "/")
                println(io, avail)
            end
        end
    end

    input = readflexpartinput(pnfile)
    # Expand file/inventory and non-box geometry releases into point sources
    # directly in the RELEASE vectors (fast for >10⁴ locations).  Domainfill
    # releases are engine-placed, not expanded.
    any(r -> r.kind == "volcanic" ||
             (r.kind != "domainfill" && (!isempty(r.from_file) || r.geometry != "box")),
        spec.release) && _expand_releases!(input, spec, specnum)
    # Domainfill needs a slot pool ≈ 2× the fill population: the fill itself plus
    # headroom for boundary inflow between inject and cull (matches the harness).
    df_parts = sum(r -> r.kind == "domainfill" ? r.particles : 0, spec.release; init=0)
    # ── Domain crop ([meteorology].crop, default "auto") ──────────────────────
    # Load only the release box + 2 halo cells of each met field instead of the
    # full (often global) grid.  SAFE ONLY for domain-fill runs: the per-step
    # cull keeps every parcel inside the box, so met outside box+halo is
    # provably never touched.  Ordinary releases roam arbitrarily far → never
    # auto-cropped.  An EXPLICIT crop = [lon0, lon1, lat0, lat1] is honoured on
    # any run (power users who know their trajectories stay inside), with a loud
    # warning: parcels reaching the crop edge are terminated there, exactly like
    # a regional-met run's domain edge.
    crop = _crop_window(spec; warn=true)
    # seaspray releases schedule 0 particles — reserve their parcel-pool slots
    # from the TOML `particles` field (the spawn CAP for the run)
    ss_parts = sum((r.kind == "seaspray" ? r.particles : 0) for r in spec.release; init=0)
    rt = build_runtime(input; prefetch=prefetch,
                       extra_slots=extra_slots + 2 * df_parts + ss_parts,
                       crop=crop)
    ctx = rt.ctx
    if crop !== nothing
        # fail loudly if the output grid pokes outside the cropped met window
        g1 = ctx.nests[1].gridconfig
        lonspan = (Float64(g1.lonmin), Float64(g1.lonmin) + (Int(g1.Ni) - 1) * Float64(g1.dlon))
        latspan = (Float64(g1.latmin), Float64(g1.latmin) + (Int(g1.Nj) - 1) * Float64(g1.dlat))
        oc0 = ctx.output
        olon0 = Float64(oc0.lon0); olon1 = olon0 + Float64(oc0.dlon) * Int(oc0.Niout)
        olat0 = Float64(oc0.lat0); olat1 = olat0 + Float64(oc0.dlat) * Int(oc0.Njout)
        # map into the crop's (possibly seam-crossing, monotone) lon frame
        _inlon(l) = (lonspan[1] + mod(l - lonspan[1], 360.0)) <= lonspan[2] + 1e-6
        (_inlon(olon0) && _inlon(olon1 - 1e-9) &&
         latspan[1] - 1e-6 <= olat0 && olat1 <= latspan[2] + 1e-6) ||
            error("output grid ($(olon0)…$(olon1)E, $(olat0)…$(olat1)N) lies outside the cropped met window ($(lonspan[1])…$(lonspan[2])E, $(latspan[1])…$(latspan[2])N) — shrink the OUTGRID or set [meteorology].crop=\"off\"")
    end

    # ── turbulence tag per nest ───────────────────────────────────────────────
    tag = let t = spec.physics.turbulence
        t.scheme == "none"    ? EASLEE.Types.NoTurb() :
        t.scheme == "thomson" ? EASLEE.Types.Thomson(true, false) :
                                EASLEE.Types.Hanna(t.hanna_onthefly, Float32(t.sigw_ft))
    end
    # ── wind dimensionality per nest (freeze suppressed parcel indices) ───────
    wmask = _wind_mask(spec.wind)
    for n in 1:length(ctx.nests)
        ctx.nests[n].nestrunconfig.turbulence = tag
        ctx.nests[n].nestrunconfig.wind_mask  = wmask
        ctx.nests[n].nestrunconfig.mesoscale  = spec.physics.turbulence.mesoscale
        ctx.nests[n].nestrunconfig.fturbmeso  = Float32(spec.physics.turbulence.fturbmeso)
        ctx.nests[n].nestrunconfig.kh_ft      = Float32(spec.physics.turbulence.kh_ft)
        # Only native-e nests can supply a horizontal variance; on a Hanna-fill
        # nest (ERA5) the archived tke is 0.5*sigma_w^2 and carries none, so the
        # term stays off there instead of fabricating one.
        ctx.nests[n].nestrunconfig.horiz_tke  = spec.physics.turbulence.horizontal_tke &&
                                                ctx.nests[n].nestrunconfig.vmix_native
        # Hanna horizontal Langevin on the Hanna-fill nests only: a native-e nest
        # already has the real horizontal variance (horiz_tke); Hanna's PBL
        # profiles are the observation-based substitute where the archive has
        # no horizontal information.
        ctx.nests[n].nestrunconfig.horiz_hanna = spec.physics.turbulence.horizontal_hanna &&
                                                 !ctx.nests[n].nestrunconfig.vmix_native
        ctx.nests[n].nestrunconfig.mo_surface = spec.physics.surface_layer_wind
    end
    @info(spec.physics.surface_layer_wind ?
          "Monin-Obukhov surface-layer wind ON ([physics].surface_layer_wind)" :
          "Monin-Obukhov surface-layer wind OFF - wind constant below the lowest model level (legacy)")
    if spec.physics.turbulence.horizontal_tke
        spec.physics.mixing.scheme != "off" &&
            error("[physics.turbulence].horizontal_tke and [physics.mixing] are both armed: one K either moves parcels or exchanges mass, never both (double counting)")
        spec.physics.turbulence.scheme == "none" &&
            error("[physics.turbulence].horizontal_tke needs a turbulence scheme (scheme = \"thomson\"); it rides on the same per-sub-step TKE the vertical Langevin step reads")
        # native-e is a property of the DRIVER namelist (listvariables/*.vars
        # #info tke_kind), published by the reader via set_vertical_mixing!.
        # A Hanna-fill driver archives 0.5*sigma_w^2 only — no horizontal
        # information — so deriving sigma_u from it would be fabricating data.
        let armed = [n for n in 1:length(ctx.nests)
                       if ctx.nests[n].nestrunconfig.vmix_native]
            isempty(armed) &&
                error("[physics.turbulence].horizontal_tke requires a native-e driver (listvariables/*.vars #info tke_kind = e, e.g. AROME/MesoNH): this run's met archives 0.5*sigma_w^2 (Hanna fill), which carries no horizontal variance")
            for n in armed
                nrc = ctx.nests[n].nestrunconfig
                @info "TKE-driven horizontal turbulence ON ([physics.turbulence].horizontal_tke)" nest=n r_w=nrc.vmix_rw r_h=0.5*(1-nrc.vmix_rw) c_k=nrc.vmix_ck
            end
            length(armed) < length(ctx.nests) &&
                @info "  (off on the Hanna-fill nest(s): no horizontal variance in 0.5*sigma_w^2)" nests_off=setdiff(1:length(ctx.nests), armed)
        end
        # 2026-09-21: the KA kernel now carries the horizontal Langevin term
        # (advance_ka.jl, use_htke), so this refusal is retired.  The guard that
        # matters lives in advance_dispatch!: an explicitly requested backend
        # refuses to downgrade, so an unported switch can no longer be dropped
        # silently on ANY path.
    end
    if spec.physics.turbulence.horizontal_hanna
        spec.physics.mixing.scheme != "off" &&
            error("[physics.turbulence].horizontal_hanna and [physics.mixing] are both armed: one K either moves parcels or exchanges mass, never both (double counting)")
        spec.physics.turbulence.scheme == "none" &&
            error("[physics.turbulence].horizontal_hanna needs a turbulence scheme (scheme = \"thomson\" or \"hanna\"); it rides on the same per-sub-step PBL parameters")
        let armed = [n for n in 1:length(ctx.nests) if ctx.nests[n].nestrunconfig.horiz_hanna]
            isempty(armed) &&
                error("[physics.turbulence].horizontal_hanna: every nest is a native-e driver — use horizontal_tke there (Hanna's PBL profiles are the substitute for archives WITHOUT a horizontal variance)")
            @info "Hanna horizontal Langevin ON ([physics.turbulence].horizontal_hanna): σ_u/σ_v/T_Lu/T_Lv in the PBL, FLEXPART convention (grid components)" nests=armed
            length(armed) < length(ctx.nests) &&
                @info "  (off on the native-e nest(s): horizontal_tke is the term there)" nests_off=setdiff(1:length(ctx.nests), armed)
        end
        # 2026-09-21: ported to the KA kernel (use_hhanna) — refusal retired, see
        # the horizontal_tke note above.
    end
    if spec.physics.turbulence.kh_ft > 0
        spec.physics.mixing.scheme != "off" &&
            error("[physics.turbulence].kh_ft and [physics.mixing] are both armed: one K either moves parcels or exchanges mass, never both (double counting)")
        @info "free-atmosphere horizontal diffusivity ON ([physics.turbulence].kh_ft)" Kh_m2s=spec.physics.turbulence.kh_ft sigma_7d_km=round(sqrt(2*spec.physics.turbulence.kh_ft*7*86400)/1000, digits=1)
    end
    wmask == (true,true,true) || @info "wind dimensionality" wind=spec.wind update_xyl=wmask
    if spec.physics.turbulence.mesoscale
        @info "mesoscale velocity fluctuations ON" fturbmeso=spec.physics.turbulence.fturbmeso windint_s=ctx.nests[1].nestrunconfig.windint_s
        get(ENV, "EASLEE_ADVANCE", "") == "ka_gpu" &&
            @warn "mesoscale fluctuations are CPU-advance only — the KA/GPU path ignores them"
    end

    # ── polar-stereographic caps ([physics].polar_stereo_lat) ────────────────
    # Poleward of ±polar_stereo_lat the horizontal position update runs in a
    # pole-tangent stereographic plane (FLEXPART switchnorth/switchsouth).
    # Enabled per pole and per nest only where it is meaningful: a global
    # Rectilinear grid whose last/first latitude row reaches that pole.
    if spec.physics.polar_stereo_lat < 90.0
        plat = Float32(spec.physics.polar_stereo_lat)
        for n in 1:length(ctx.nests)
            g = ctx.nests[n].gridconfig
            (g isa EASLEE.Types.Rectilinear && g.globalx && g.dlat > 0.0f0) || continue
            nrc = ctx.nests[n].nestrunconfig
            latmax = g.latmin + Float32(g.Nj - 1) * g.dlat
            if latmax >= 90.0f0 - g.dlat / 2
                nrc.polar_n = true
                nrc.switchnorth_j = (plat - g.latmin) / g.dlat + 1.0f0
            end
            if g.latmin <= -90.0f0 + g.dlat / 2
                nrc.polar_s = true
                nrc.switchsouth_j = (-plat - g.latmin) / g.dlat + 1.0f0
            end
            (nrc.polar_n || nrc.polar_s) &&
                @info "polar-stereographic transport ON" nest=n switch_lat=spec.physics.polar_stereo_lat north=nrc.polar_n south=nrc.polar_s
        end
    end

    # :none | :uniform | :adaptive | :covariance | :covshape
    # (kernel="covariance" splits on [output.covariance].mode: physical → P is
    #  the width; shape → kNN area × P anisotropy)
    ctx.output.kernel = spec.output.kernel == "covariance" && spec.output.cov_mode == "shape" ?
        :covshape : Symbol(spec.output.kernel)
    ctx.output.resolution = spec.output.resolution
    # adaptive kNN bandwidth refresh cadence: once per output frame (the kNN
    # pass is the dominant output cost).  samples-per-frame = interval/sample.
    ctx.output.bw_stride = spec.output.sample_s > 0 ?
        max(1, round(Int, spec.output.interval_s / spec.output.sample_s)) : 1
    # physical (relative-dispersion) bandwidth floor — grid deposit, kNN pass
    # and receptors all read these (0 = off)
    ctx.output.bwfloor_sv = spec.output.bandwidth_floor == "dispersion" ?
        spec.output.floor_sigmav : 0.0
    ctx.output.bwfloor_tl = spec.output.floor_tl_s
    ctx.output.bwfloor_sv > 0.0 &&
        @info "adaptive-kernel bandwidth floor: Taylor dispersion" sigma_v=ctx.output.bwfloor_sv T_L_s=ctx.output.bwfloor_tl
    spec.output.kernel in ("adaptive", "covariance") && get(ENV, "EASLEE_ADVANCE", "") == "ka_gpu" &&
        @warn "output.kernel=\"$(spec.output.kernel)\" uses the CPU deposit path (GPU deposit stays uniform)"

    # ── plume rise per release box (box index → cfg; applied at release time) ──
    # Expansion multiplies release boxes (volcanic sublayers, file/geometry
    # points), so key the plume config by the EXPANDED box index (relbox), via
    # the _box_specidx map; without expansion, box index == spec index.
    specidx = haskey(input, "_box_specidx") ? Vector{Int}(input["_box_specidx"]) :
              collect(1:length(spec.release))
    # per-population adaptive-kernel bandwidths (task #35): population = SPEC
    # release (all expanded boxes of one release share a KDE sample; domainfill
    # and boundary parcels with relbox=0 form population 0)
    (spec.output.kernel == "adaptive" || ctx.output.kernel === :covshape) &&
        EASLEE.Operators.set_output_popmap!(specidx)
    bphase  = haskey(input, "_box_phase") ? Vector{Int}(input["_box_phase"]) :
              zeros(Int, length(specidx))
    # OSSA sea-spray sources (kind = "seaspray"): precompute the ocean cells
    # per reserved box; parcels are spawned by apply_seaspray! in run_spec!
    ssprays = EASLEE.Operators.SeaSprayBox[]
    for (b, kr) in enumerate(specidx)
        (1 <= kr <= length(spec.release) && spec.release[kr].kind == "seaspray") || continue
        ssb = EASLEE.Operators.build_seaspray(ctx, rt.metfields2D[1], b, spec.release[kr])
        ssb === nothing && continue
        push!(ssprays, ssb)
        @info "sea-spray source armed (OSSA)" species=spec.release[kr].species ocean_cells=length(ssb.ii) quantum_kg=ssb.quantum gust_sigma_ms=ssb.gust_sigma
    end
    pcfg = Dict{Int,NamedTuple}()
    for (b, k) in enumerate(specidx)
        r = spec.release[k]
        lowercase(r.plume_scheme) == "none" && continue
        ph = b <= length(bphase) ? bphase[b] : 0
        # per-phase source pulse (multi-blast scenario) overrides the release-level
        # MER/T0/water_frac — e.g. Hunga Tonga climactic vs sustained phases
        mer, tk, nw = r.plume_mer, r.plume_temp_k, r.plume_water_frac
        if ph >= 1 && ph <= length(r.volc_pulses)
            mer, tk, nw = r.volc_pulses[ph]
        end
        vm = r.plume_vent_mix ?
             EASLEE.Operators.VentMix(T_mag = r.plume_magma_t,
                                      n_mag = r.plume_magma_water,
                                      x_core = r.plume_sw_core,
                                      x_edge = r.plume_sw_edge,
                                      p_mix = r.plume_mix_exp) : nothing
        pcfg[b] = (scheme = r.plume_scheme, frp = r.plume_frp_w, area = r.plume_fire_area,
                   profile = r.plume_profile, bottom = r.plume_bottom,
                   mer = mer, temp = tk, waterfrac = nw,
                   mode = r.plume_mode, washout = r.plume_washout,
                   entrainment = r.plume_entrainment, species = r.species,
                   umbrella = r.plume_umbrella, spread_lambda = r.plume_spread_lambda,
                   umbrella_rsource = r.plume_umbrella_rsource,
                   umbrella_windshield = r.plume_umbrella_windshield,
                   spread_umin = r.plume_spread_umin,
                   spread_tmax = r.plume_spread_tmax,
                   spread_hmin = r.plume_spread_hmin,
                   spread_zgain = r.plume_spread_zgain,
                   spread_zg_z0 = r.plume_spread_zg_z0,
                   spread_zg_z1 = r.plume_spread_zg_z1,
                   spread_azgain_n = r.plume_spread_azgain_n,
                   spread_azgain_s = r.plume_spread_azgain_s,
                   spread_az_z0 = r.plume_spread_az_z0,
                   spread_az_z1 = r.plume_spread_az_z1,
                   qc_crit = r.plume_qc_crit, vapor_cap = r.plume_vapor_cap,
                   umb_r_med = r.plume_umb_r_med,
                   umb_r_sg = r.plume_umb_r_sg, umb_h = r.plume_umb_h,
                   umb_t = r.plume_umb_t, solid_r = r.plume_solid_r,
                   ventmix = vm, specidx = k, phase = ph)
    end
    EASLEE.Operators.set_plume_config!(pcfg)
    isempty(pcfg) || @info "plume rise enabled" boxes=sort(collect(keys(pcfg))) schemes=[pcfg[k].scheme for k in sort(collect(keys(pcfg)))]

    # ── [release.chem_init]: release-scoped chemistry IC, keyed by EXPANDED box ──
    # (relbox); a parcel's relbox → obj.relchem → apply_release_chem! at seeding.
    relchem = Dict{Int,Dict{String,Float64}}()
    relh2oplume = Set{Int}()
    for (b, k) in enumerate(specidx)
        (1 <= k <= length(spec.release)) || continue
        isempty(spec.release[k].chem_init) || (relchem[b] = spec.release[k].chem_init)
        spec.release[k].chem_h2o_plume && push!(relh2oplume, b)
    end
    isempty(relchem) || @info "release-scoped chemistry IC" boxes=sort(collect(keys(relchem)))
    # ICE shadow boxes of the bridged releases ([species.convert] targets): the
    # per-step water/ice sync must follow a parcel when it freezes into the
    # shadow box and comes back — the box SET is (vapour boxes, their shadows).
    relh2oice = Set{Int}()
    if !isempty(relh2oplume) && haskey(input, "_box_convert")
        bcv = Vector{Int}(input["_box_convert"])
        for b in relh2oplume
            1 <= b <= length(bcv) && bcv[b] > 0 && push!(relh2oice, bcv[b])
        end
    end
    isempty(relh2oplume) || @info "plume-water bridge (physics → cs.h2o/cs.ice, per step)" vapor=sort(collect(relh2oplume)) ice=sort(collect(relh2oice))

    # ── physics objects ───────────────────────────────────────────────────────
    has_aerosol = any(sp -> sp.kind == "aerosol", spec.species)
    conv = spec.physics.convection
    wet_any = spec.physics.deposition.wet || any(sp -> sp.wet_deposition, spec.species)

    convfields = nothing; diagconv = nothing; diagconvparams = nothing
    if conv.scheme == "massflux" || wet_any
        convfields = [EASLEE.ConvFields(rt.metfields3D[n]) for n in 1:rt.nest]
        all(c -> c === nothing, convfields) && (convfields = nothing)
    end
    if conv.scheme == "cape"
        # 3-D dict passed too: archived qv (when mapped) feeds the source-parcel θe
        # — required on ERA5, whose reader stores dry ρ (Tv reconstruction → q≡0).
        diagconv = [EASLEE.DiagConvFields(rt.metfields2D[n], rt.metfields3D[n]) for n in 1:rt.nest]
        diagconvparams = EASLEE.DiagConvParams(cape_min=Float32(conv.cape_min),
                                               cin_max=Float32(conv.cin_max),
                                               tau_mix=Float32(conv.tau_mix_s))
    end

    # wet/dry deposition: aerosol species use the size-dependent path inside the
    # AerosolModel; gas species use the bulk WetDepParams.
    wetdepparams = nothing; wetdepstate = nothing
    gas_wet = spec.physics.deposition.wet || any(sp -> sp.wet_deposition && sp.kind == "gas", spec.species)
    if gas_wet && convfields !== nothing
        # review 2026-09-04 A6: the bulk gas scheme used to be switched off whenever an aerosol
        # species was present (it had no per-box mask and would have scavenged the aerosol boxes a
        # second time).  WetDepParams.boxon now restricts it to the gas boxes that asked for it.
        gasmask = Bool[]
        if haskey(input, "_box_species")
            spk = Dict(sp.id => sp for sp in spec.species)
            gasmask = [haskey(spk, sid) && spk[sid].kind == "gas" &&
                       (spec.physics.deposition.wet || spk[sid].wet_deposition)
                       for sid in Vector{String}(input["_box_species"])]
        elseif has_aerosol
            @warn "gas wet_deposition SKIPPED: no per-box species map, cannot separate gas from aerosol boxes"
            gasmask = nothing
        end
        if gasmask !== nothing
            wetdepparams = EASLEE.WetDepParams(boxon = gasmask)
            wetdepstate  = EASLEE.WetDepState(ctx)
            isempty(gasmask) || @info "gas wet deposition armed" boxes=count(gasmask) of=length(gasmask)
        end
    end

    aerosol = nothing
    aero_dry = true; aero_wet = true        # review 2026-09-04 A6: gates passed to step_model!
    if has_aerosol
        sp = first(filter(s -> s.kind == "aerosol", spec.species))   # ONE bins set for all aerosol species (validate_runspec enforces agreement — review A7 item 4)
        # review 2026-09-04 A7 item 3: refractive_index was parsed and never passed (m_ri)
        m_ri = EasleeConfig.parse_refractive_index(sp.refractive_index)   # nothing = "auto" → type default
        bins = if sp.size_distribution == "bins"
            EASLEE.AerosolBins(sp.bin_edges_m; rho_p=sp.density,
                species=Symbol(sp.aerosol_type), m_ri=m_ri,
                kappa = sp.kappa >= 0 ? Float32(sp.kappa) : nothing)
        else
            EASLEE.AerosolBins(dg=sp.dg_m, sigma_g=sp.sigma_g, nbins=sp.nbins,
                rho_p=sp.density, species=Symbol(sp.aerosol_type), m_ri=m_ri,
                kappa = sp.kappa >= 0 ? Float32(sp.kappa) : nothing)
        end
        m_ri === nothing || @info "aerosol refractive index from [species].refractive_index" species=sp.id m_ri=m_ri
        # mixed gas+aerosol runs: restrict the aerosol physics to the boxes
        # whose species is an aerosol (legacy single-species: all boxes)
        bxa = Bool[]
        if haskey(input, "_box_species") && any(sp -> sp.kind != "aerosol", spec.species)
            byid_a = Dict(spp.id => spp for spp in spec.species)
            bxa = [get(byid_a, sid, nothing) !== nothing &&
                   byid_a[sid].kind == "aerosol"
                   for sid in Vector{String}(input["_box_species"])]
        end
        # emission-driven aerosol = an aerosol-species DOMAINFILL + [chemistry.emissions]:
        # the parcels' state.mass is the domainfill AIR mass, so aerosol mass rides in pmass.
        aero_emit = !isempty(spec.chemistry.emissions_dir) &&
                    any(r -> r.kind == "domainfill" && r.species == sp.id, spec.release)
        # per-box settling mask from [species].settling (all-true → empty = legacy)
        sxa = Bool[]
        if haskey(input, "_box_species")
            byid_s = Dict(spp.id => spp for spp in spec.species)
            sxa = [get(byid_s, sid, nothing) !== nothing && byid_s[sid].settling
                   for sid in Vector{String}(input["_box_species"])]
            all(sxa) && (sxa = Bool[])
        end
        # review 2026-09-04 A6: aerosol dry/wet deposition used to run regardless of
        # [species].dry_deposition/wet_deposition and [physics.deposition]; both are now gated
        # in step_model! through aerosol_drydep / aerosol_wetdep (RunObjects.aero_*).
        aer_sp   = filter(s -> s.kind == "aerosol", spec.species)
        aero_dry = spec.physics.deposition.dry || any(s -> s.dry_deposition, aer_sp)
        aero_wet = spec.physics.deposition.wet || any(s -> s.wet_deposition, aer_sp)
        aero_dry || @info "aerosol dry deposition OFF (no aerosol species sets dry_deposition, [physics.deposition].dry = false)"
        aero_dry && !spec.physics.deposition.dry && !all(s -> s.dry_deposition, aer_sp) &&
            @warn "aerosol dry deposition applies to ALL aerosol boxes: the per-species dry_deposition switch cannot be honoured (AerosolModel has one href, no per-box deposition mask)" on=[s.id for s in aer_sp if s.dry_deposition] off=[s.id for s in aer_sp if !s.dry_deposition]
        aero_wet || @info "aerosol wet deposition OFF (no aerosol species sets wet_deposition, [physics.deposition].wet = false)"
        aerosol = EASLEE.AerosolModel(length(rt.state.status), ctx, bins;
                                      href = 30.0,   # dry deposition gated in step_model! (aerosol_drydep)
                                      boxon=bxa, settleon=sxa, emit=aero_emit)
    end

    hybridvert = nothing
    if spec.physics.vertical.mode in ("hybrid", "isentropic")
        v = spec.physics.vertical
        ptdir = v.pt_winds == "auto" ? spec.meteorology[1].path :
                v.pt_winds == "off"  ? "" : v.pt_winds
        hybridvert = EASLEE.HybridVert(length(rt.state.status), rt;
            sigma_r=v.sigma_r, sigma_t=v.sigma_t, diabatic=v.diabatic,
            ptdir=ptdir, ptpattern=v.pt_pattern)
    end

    bench = spec.benchmark ? EASLEE.BenchTimer() : nothing

    # ── chemistry (optional) ──────────────────────────────────────────────────
    chemstate, chemcfg, chembg, cheminit, chemtl = build_chemistry(spec, length(rt.state.status))

    # ── per-species process physics (decay / gas dry dep / below-cloud A,B) ──
    # relbox → species id: recorded by _expand_releases!, else the plain-box order.
    box_ids = haskey(input, "_box_species") ? Vector{String}(input["_box_species"]) :
              String[r.species for r in spec.release
                     if r.kind == "release" && isempty(r.from_file) && r.geometry == "box"]
    spphys = nothing
    icephys = nothing
    convphys = nothing
    pradphys = nothing
    # ── Per-parcel DIAGNOSED layer thickness h_v ([physics].prognostic_airm) ──
    # update_airm! fills this with the local vertical spread of the parcel
    # cloud (kNN).  It used to be allocated ONLY inside the [species.ice]
    # branch, so a run using the [species.convert]/[species.fall] ice path
    # (no [species.ice] block — e.g. ht_full3) left it empty and
    # plume_radiation silently fell back to the CONSTANT h_ref = 1 km.  Since
    # the lw_kappa cooling is an edge effect ∝ 1/h, that made the descent
    # magnitude untrustworthy (handoff §2.3).  Allocate it once, here, and
    # share the SAME vector with every consumer.
    hvshared = spec.physics.prognostic_airm ?
               zeros(Float32, length(rt.state.status)) : Float32[]
    # [tracers]: passive transport-history tags (independent of release boxes).
    # Capture the STATIC land-sea mask once (raw archived 2-D dict, nest 1) for the
    # marine-vs-continental tag; the typed metfields2D struct carries only PBL fields.
    tag_lsm = nothing
    if spec.physics.tags_enabled
        m2raw = isdefined(rt, :metfields2D) && !isempty(rt.metfields2D) ? rt.metfields2D[1] : nothing
        if m2raw !== nothing && haskey(m2raw, "lsm")
            cand = Matrix{Float32}(Float32.(@view m2raw["lsm"][:, :, 1]))
            # Some readers register an "lsm" key as a zero PLACEHOLDER (no real
            # land-sea mask archived).  An all-ocean mask is not usable — treat
            # it as absent so cont_frac isn't a silent, misleading 0.
            if any(>=(0.5f0), cand)
                tag_lsm = cand
                @info "transport tags: land-sea mask captured ([tracers])" frac_land=round(count(>=(0.5f0), cand) / length(cand), digits=3)
            else
                @warn "transport tags: archived land-sea mask is all-ocean (a placeholder) — the continental (cont_frac) tag will be 0; use ERA5 with a real _SFC_ lsm companion to enable marine/continental attribution" maxlog=1
            end
        else
            @warn "transport tags: no land-sea mask (lsm) in the archive — the continental (cont_frac) tag will be 0; add ERA5 _SFC_ companion fields to enable it" maxlog=1
        end
    end
    tagphys = spec.physics.tags_enabled ?
        EASLEE.TransportTags(length(rt.state.status);
                             strat_theta=spec.physics.tags_strat_theta, lsm=tag_lsm) : nothing
    tagphys === nothing ||
        @info "transport tags armed ([tracers])" strat_theta_K=spec.physics.tags_strat_theta lsm=(tag_lsm !== nothing)
    if !isempty(box_ids)
        nb   = length(box_ids)
        lam  = zeros(Float32, nb)
        weta = fill(-1.0f0, nb); wetb = fill(-1.0f0, nb)
        ddep = falses(nb)
        vd   = fill(-1.0f0, nb); rel = fill(1.0f0, nb)
        hen  = zeros(Float32, nb); f0v = zeros(Float32, nb)
        setr = zeros(Float32, nb); setrho = fill(917.0f0, nb)   # [species.settling]
        byid = Dict(sp.id => sp for sp in spec.species)
        for (b, sid) in enumerate(box_ids)
            sp = get(byid, sid, nothing)
            sp === nothing && continue
            sp.decay_halflife_s > 0 && (lam[b] = Float32(log(2.0) / sp.decay_halflife_s))
            if sp.settle_r_eff > 0                      # espèce condensée qui tombe
                setr[b]   = Float32(sp.settle_r_eff)
                setrho[b] = Float32(sp.settle_rho)
            end
            if sp.kind == "gas" && sp.wet_deposition && sp.weta_gas >= 0
                weta[b] = Float32(sp.weta_gas)
                wetb[b] = Float32(max(sp.wetb_gas, 0.0))
            end
            if sp.kind == "gas" && (sp.dry_deposition || spec.physics.deposition.dry)
                if sp.dryvel >= 0 || (sp.reldiff > 0 && sp.henry >= 0)
                    ddep[b] = true
                    vd[b]   = Float32(sp.dryvel)          # <0 → resistance model
                    rel[b]  = Float32(max(sp.reldiff, 1.0))
                    hen[b]  = Float32(max(sp.henry, 0.0))
                    f0v[b]  = Float32(clamp(sp.f0 < 0 ? 0.0 : sp.f0, 0.0, 1.0))
                else
                    @warn "species \"$sid\": dry_deposition requested but no dryvel/reldiff/henry (species DB) — gas dry deposition skipped"
                end
            end
        end
        if any(>(0f0), lam) || any(ddep) || any(>=(0f0), weta) || any(>(0f0), setr)
            depgrid = [zeros(Float64, Int(n.gridconfig.Ni), Int(n.gridconfig.Nj))
                       for n in ctx.nests]
            spphys = EASLEE.SpeciesPhysics(lam, weta, wetb, ddep, vd, rel, hen, f0v,
                                           30.0f0, depgrid, setr, setrho)
            any(>(0f0), setr) && @info "species settling armed ([species.settling])" boxes=count(>(0f0), setr) r_eff_um=round(maximum(setr)*1e6, digits=1)
            @info "species physics" decay_boxes=count(>(0f0), lam) gas_drydep_boxes=count(ddep) species_wet_ab=count(>=(0f0), weta)
        end
        # [species.ice]: transport-phase freezing/sedimentation for tagged water
        if any(sp -> sp.ice, spec.species)
            ion = falses(nb); irhi = fill(1.5f0, nb); ireff = zeros(Float32, nb)
            ihl = fill(500.0f0, nb); iq0 = zeros(Float32, nb)
            iscav = zeros(Float32, nb)
            for (b, sid) in enumerate(box_ids)
                sp = get(byid, sid, nothing)
                sp === nothing && continue
                sp.ice_scav_uptake > 0 && (iscav[b] = Float32(sp.ice_scav_uptake))
                sp.ice || continue
                ion[b]  = true
                irhi[b] = Float32(sp.ice_rhi_crit)
                ireff[b] = Float32(sp.ice_r_eff_m)      # v_t(T,p) computed per parcel
                ihl[b]  = Float32(sp.ice_h_layer_m)
                iq0[b]  = Float32(sp.ice_q0)
            end
            if any(ion)
                # [physics].prognostic_airm → per-parcel DIAGNOSED escape scale
                # (update_airm! fills it with the local plume thickness) —
                # the SHARED vector, so ice and plume_radiation see the same h
                icephys = EASLEE.IcePhysics(ion, irhi, ireff, ihl, iq0,
                                            zeros(Float64, nb), iscav,
                                            zeros(Float64, nb), hvshared)
                # indicative v_t at 25 km (220 K / 25 hPa) — the actual fall speed
                # is evaluated at each parcel's own T, p (Stokes–Cunningham)
                vt25 = EASLEE.Operators.settling_velocity(2f0 * maximum(ireff), 917f0, 220f0, 2500f0)
                @info "ice physics armed ([species.ice])" boxes=count(ion) partner_boxes=count(>(0f0), iscav) r_eff_um=round(maximum(ireff)*1e6, digits=1) v_t_25km_cms=round(vt25*100, digits=2) rhi_crit=maximum(irhi)
            end
        end
        # [physics.plume_radiation]: LW cooling from excess plume water
        if !isempty(spec.physics.prad_species)
            pon  = [sid in spec.physics.prad_species  for sid in box_ids]
            ppar = [sid in spec.physics.prad_partners for sid in box_ids]
            if any(pon)
                # per-parcel layer thickness: update_airm! fills the SHARED hv,
                # which is the same diagnostic the ice module uses.  No longer
                # conditional on [species.ice] being armed — the convert/fall
                # ice path gets the trajectory-derived thickness too
                # (empty only when prognostic_airm is off -> h_ref fallback).
                hvshare = hvshared
                # [physics.plume_radiation].tprime_tau: prognostic cold anomaly
                # T′ — allocated per slot only when armed (0 = legacy path)
                tpt = spec.physics.prad_tprime_tau
                pradphys = EASLEE.PlumeRadiation(pon, ppar,
                    spec.physics.prad_c_rad / 86400.0,
                    spec.physics.prad_lw_kappa, spec.physics.prad_q_bg,
                    spec.physics.prad_feedback,
                    zeros(Float32, length(rt.state.status)),
                    hvshare, spec.physics.prad_h_ref,
                    tpt,
                    tpt > 0 ? zeros(Float32, length(rt.state.status)) : Float32[])
                @info "plume radiation armed ([physics.plume_radiation])" boxes=count(pon) partner_boxes=count(ppar) c_rad_K_day_ppmv=spec.physics.prad_c_rad lw_kappa=spec.physics.prad_lw_kappa thickness=(isempty(hvshare) ? "h_ref $(spec.physics.prad_h_ref) m" : "per-parcel (ice.hv)") q_bg_ppmv=spec.physics.prad_q_bg feedback=spec.physics.prad_feedback tprime_tau_h=(tpt > 0 ? round(tpt / 3600, digits=2) : "OFF (legacy: instant equilibrium)")
            end
        end
        # [species.convert]: stochastic in-flight conversion to the shadow boxes
        if any(sp -> !isempty(sp.convert_to), spec.species) && haskey(input, "_box_convert")
            bcv = Vector{Int}(input["_box_convert"])
            ncb = length(bcv)
            clam = zeros(Float32, ncb); ctgt = zeros(Int32, ncb); crat = ones(Float32, ncb)
            cmode = zeros(Int8, ncb); crhi = fill(1.5f0, ncb)   # gate thermodynamique
            clsg  = zeros(Float32, ncb)                        # largeur lamina sous-maille
            csbr  = zeros(Float32, ncb)                        # rayon du cristal (sublimation diffusive)
            bids2 = Vector{String}(input["_box_species"])
            for b in 1:ncb
                bcv[b] > 0 || continue
                sp = byid[bids2[b]]
                clam[b] = Float32(1.0 / sp.convert_tau_s)
                cmode[b] = sp.convert_mode == "freeze"    ? Int8(1) :
                           sp.convert_mode == "sublimate" ? Int8(2) : Int8(0)
                crhi[b]  = Float32(sp.convert_rhi_crit)
                clsg[b]  = Float32(sp.convert_lam_sigma)
                # the sublimating species IS the ice: its own settling radius
                sp.convert_mode == "sublimate" && sp.settle_r_eff > 0 &&
                    (csbr[b] = Float32(sp.settle_r_eff))
                ctgt[b] = Int32(bcv[b])
                if sp.convert_mass_ratio > 0
                    crat[b] = Float32(sp.convert_mass_ratio)
                else
                    dsp = byid[sp.convert_to]
                    if sp.molar_mass > 0 && dsp.molar_mass > 0
                        crat[b] = Float32(dsp.molar_mass / sp.molar_mass)
                    else
                        @warn "species convert $(sp.id)→$(sp.convert_to): molar masses unknown — mass_ratio = 1" maxlog=1
                    end
                end
            end
            # ── SECTIONAL ladder ([species.convert].sectional) ───────────────
            # Order the sublimating (ICE) SPECIES by crystal radius and link
            # each to the next smaller one.  Boxes are expanded per release, so
            # a demotion must land on a box of the target species that shares
            # this box's release identity (`bspecidx`) — otherwise the parcel
            # would jump releases and corrupt the output attribution.
            crice = Float32[]; cbdn = zeros(Int32, ncb); cbre = zeros(Float32, ncb)
            if any(sp -> sp.convert_sectional, spec.species)
                bidx = haskey(input, "_box_specidx") ?
                       Vector{Int}(input["_box_specidx"]) : zeros(Int, ncb)
                # radius per sublimating species id (the ICE bins)
                rsp = Dict{String,Float64}()
                for sp in spec.species
                    sp.convert_mode == "sublimate" && sp.settle_r_eff > 0 &&
                        (rsp[sp.id] = sp.settle_r_eff)
                end
                ladder = sort(collect(keys(rsp)); by = s -> rsp[s])   # small → large
                pos = Dict(s => k for (k, s) in enumerate(ladder))
                nlad = length(ladder)
                for b in 1:ncb
                    sid = bids2[b]
                    k = get(pos, sid, 0)
                    k > 0 || continue
                    byid[sid].convert_sectional || continue
                    rk = rsp[ladder[k]]
                    # lower edge = geometric midpoint with the bin below; for
                    # the smallest bin, one ratio-step below it (crystal gone)
                    rlo = k > 1 ? sqrt(rsp[ladder[k-1]] * rk) :
                          (nlad > 1 ? rk * sqrt(rk / rsp[ladder[2]]) : 0.3 * rk)
                    cbre[b] = Float32(rlo)
                    if k > 1
                        tgt = ladder[k-1]
                        # same release first, then any box of that species
                        j = findfirst(a -> bids2[a] == tgt && bidx[a] == bidx[b],
                                      1:ncb)
                        j === nothing && (j = findfirst(a -> bids2[a] == tgt, 1:ncb))
                        j === nothing || (cbdn[b] = Int32(j))
                    end                       # k == 1 → 0 = off the ladder
                end
                if any(>(0f0), cbre)
                    crice = zeros(Float32, length(rt.state.status))
                    @info "SECTIONAL ice armed ([species.convert].sectional)" bins=nlad ladder_um=round.([rsp[s]*1e6 for s in ladder], digits=2) boxes=count(>(0f0), cbre) demotable=count(>(0), cbdn) bottom_bin_boxes=count(i -> cbre[i] > 0f0 && cbdn[i] == 0, 1:ncb)
                else
                    @warn "[species.convert].sectional set but no ICE bin ladder was built — is `mode = \"sublimate\"` + [species.fall].r_eff present on the ice species?"
                end
            end
            if any(>(0), ctgt)
                convphys = EASLEE.SpeciesConvert(clam, ctgt, crat, zeros(Float64, ncb),
                                                 cmode, crhi, clsg, csbr,
                                                 crice, cbdn, cbre;
                                                 diag_slots = isempty(spec.parcel_dump) ?
                                                     0 : length(rt.state.status))
                isempty(spec.parcel_dump) ||
                    @info "convert-gate diagnostics armed ([run].parcel_dump)" period=spec.parcel_dump slots=length(rt.state.status)
                any(>(0f0), csbr) && @info "sublimation diffusion-limited ([species.fall].r_eff)" boxes=count(>(0f0), csbr) r_um=round.(sort(unique(csbr[csbr.>0]))*1e6, digits=2)
                any(>(0f0), clsg) && @info "sub-grid lamina armed ([species.convert].lamina_sigma)" boxes=count(>(0f0), clsg) sigma=maximum(clsg)
                any(!=(Int8(0)), cmode) && @info "conversion thermodynamique ([species.convert].mode)" freeze_boxes=count(==(Int8(1)), cmode) sublimate_boxes=count(==(Int8(2)), cmode)
                @info "species conversion armed ([species.convert])" converting_boxes=count(>(0), ctgt) tau_d=round(1 / (86400 * maximum(clam)), digits=1) mass_ratio=maximum(crat)
            end
        end
    end

    # cloud correction of photolysis: archived 3-D cloud fraction (CLD companions)
    chemcc = nothing
    if chemstate !== nothing && spec.chemistry.cloud
        if haskey(rt.metfields3D[1], "cc")
            chemcc = rt.metfields3D[1]["cc"]
            @info "chemistry: photolysis cloud correction ON (archived cloud fraction cc)"
        else
            @warn "chemistry.cloud=true but no 3-D cloud fraction (cc) in the met — photolysis stays clear-sky. Use a variables namelist with CLD companions (e.g. ERA5_CONV)."
        end
    end
    # per-parcel water vapour from the archived specific humidity (drives the
    # O(1D)+H2O primary OH source — a constant H2O is wrong at both ends)
    chemqv = nothing
    if chemstate !== nothing
        if haskey(rt.metfields3D[1], "qv")
            chemqv = rt.metfields3D[1]["qv"]
            @info "chemistry: per-parcel H2O from archived qv"
        else
            @warn "chemistry: no archived qv in the met — H2O fixed at $(chemcfg.H2O_ppb) ppb everywhere (radicals will be biased). Map qv in the variables namelist."
        end
    end

    # lightning NOx ([chemistry.lightning]): needs a mechanism with NO + archived
    # convective mass flux (massflux convection).
    lightning = nothing
    if chemcfg !== nothing && spec.chemistry.lightning_enabled
        no_col = findfirst(==("NO"), EASLEE.chem_species(chemcfg.mechanism))
        if no_col === nothing
            @warn "chemistry.lightning: mechanism has no NO species — disabled"
        elseif convfields === nothing
            @warn "chemistry.lightning: needs archived convective mass flux (physics.convection=massflux) — disabled"
        else
            lightning = EASLEE.LightningNOx(spec.chemistry.lightning_yield,
                            Float32(spec.chemistry.lightning_umf_thresh),
                            spec.chemistry.lightning_h_min,
                            Float32(spec.chemistry.lightning_ut_lo),
                            Float32(spec.chemistry.lightning_ut_hi), no_col)
            @info "lightning NOx ON" yield_mol_per_flash=spec.chemistry.lightning_yield h_min_km=spec.chemistry.lightning_h_min
        end
    end

    obj = RunObjects(convfields, diagconv, diagconvparams,
                     wetdepparams, wetdepstate, aerosol, hybridvert, bench,
                     chemstate, chemcfg, chembg, cheminit, chemtl, chemcc, chemqv,
                     spphys, icephys, ssprays, convphys, pradphys, tagphys, lightning,
                     isempty(relchem) ? nothing : relchem,
                     isempty(relh2oplume) ? nothing :
                         (vapor = relh2oplume, ice = relh2oice),
                     aero_dry, aero_wet)
    return rt, obj
end

# Per-step chemistry: refresh per-parcel T/P + lat/lon, seed newly released
# parcels from CAMS (init="cams"), then integrate the mechanism over
# `apply_every` transport steps.  Returns nothing.  `seeded` tracks which slots
# already got their CAMS initial composition (slots are reused, so reset on
# release is the caller's concern via the active mask).
function apply_chemistry_step!(rt, obj::RunObjects, spec, itime::Int,
                               currentdate::DateTime, indt, t_abs0::Float64,
                               dt_s::Float64, seeded::BitVector;
                               emis = nothing, dfill = nothing, kinject::Int = 0)
    obj.chemstate === nothing && return nothing
    cs = obj.chemstate; cfg = obj.chemcfg
    state = rt.state
    nestctx = rt.ctx.nests[1]; gc = nestctx.gridconfig
    EASLEE.timeit!(obj.bench, :chem_extract) do
        EASLEE.extract_tp!(cs, state, nestctx.metfields3D, indt; qv=obj.chemqv)
        # per-parcel cloud optical depth above/below for the photolysis cloud factor
        obj.chemcc === nothing ||
            EASLEE.extract_cloud_tau!(cs, state, nestctx.metfields3D, obj.chemcc, indt)
    end

    # PHYSICS → CHEMISTRY : per-step water/ice sync (masters the one-shot seed;
    # see sync_plume_water! — the parcel's water budget lives in state.mass)
    obj.relh2oplume === nothing ||
        EASLEE.timeit!(obj.bench, :chem_sync) do
            EASLEE.sync_plume_water!(cs, state, obj.relh2oplume.vapor, obj.relh2oplume.ice)
        end

    # declarative init: rotate the background bracket to `currentdate` first,
    # so the o3col/seeding below read the current slots
    if obj.cheminit !== nothing && obj.chemtl !== nothing
        EASLEE.timeline_advance!(obj.chemtl, currentdate)
        obj.chembg = obj.chemtl.bg0
        EASLEE.set_init_time!(obj.cheminit, obj.chemtl.bg0, obj.chemtl.bg1,
                              EASLEE.timeline_weight(obj.chemtl, currentdate))
    end

    # per-parcel lat/lon for spatially-resolved photolysis (solar zenith angle);
    # and, when the LUT photolysis is active with a CAMS field, a per-parcel total
    # O3 column [DU] integrated from the CAMS go3 profile (data-driven o3col axis).
    do_o3col = cfg.photolut !== nothing && obj.chembg !== nothing
    # per-parcel overhead aerosol AOD (photolysis aerosol_factor) — same
    # per-step column pattern as the O3 column, from the same background
    do_aod = !isempty(cs.aod_above) && obj.chembg !== nothing
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        la, lo = EASLEE.Operators.indyindxtolatlon(gc, Float64(state.indy[i]), Float64(state.indx[i]))
        cs.lat[i] = la; cs.lon[i] = lo
        if do_o3col
            duc = EASLEE.cams_o3_column_DU(obj.chembg, lo, la, 0.0)
            isnan(duc) || (cs.o3col[i] = duc)
        end
        if do_aod
            av = EASLEE.cams_aod550_above(obj.chembg, lo, la, cs.P_Pa[i] / 100.0, 0.0;
                                          rh_growth=spec.chemistry.aerosols_rh_growth)
            isnan(av) || (cs.aod_above[i] = av)
        end
    end

    # carried CAMS aerosol map (aero slot → background var index), refreshed
    # against the current bracket slot; empty when [chemistry.aerosols] is off
    # or the file has no aermr fields
    amap = Tuple{Int,Int}[]
    if size(cs.aero, 2) > 0 && obj.chembg !== nothing
        for (t, n) in enumerate(EASLEE.Operators.CAMS_AERO_NAMES)
            vi = EASLEE.species_index(obj.chembg, n)
            vi > 0 && push!(amap, (t, vi))
        end
    end
    @inline function _seed_aero!(i, vals)
        isempty(amap) && return
        @inbounds for t in axes(cs.aero, 2); cs.aero[i, t] = 0.0; end   # slot reuse
        @inbounds for (t, vi) in amap
            cs.aero[i, t] = max(Float64(vals[vi]), 0.0)
        end
    end

    # Seeding of newly active parcels (initial condition AND boundary/domainfill
    # reseeding — slots reused after a cull get `seeded` cleared by the caller).
    EASLEE.timeit!(obj.bench, :cams_reseed) do
    if obj.cheminit !== nothing
        # declarative path: time-interpolated sources + derivations (init.jl);
        # PSS/family steps get the run's own LUT J at the parcel when available
        ci  = obj.cheminit
        doy = spec.chemistry.doy > 0 ? spec.chemistry.doy : Dates.dayofyear(currentdate)
        t_abs = t_abs0 + itime * dt_s
        haspss = any(st -> st.kind === :steady_state || st.kind === :family, ci.steps)
        @inbounds for i in eachindex(state.status)
            (state.status[i] == 0x01 && !seeded[i]) || continue
            lat, lon = cs.lat[i], cs.lon[i]
            p_hPa = cs.P_Pa[i] / 100.0
            Jv = nothing
            if haspss && cfg.photolut !== nothing
                czr = EASLEE.cos_sza_raw(t_abs; lat_deg=lat, lon_deg=lon, doy=doy)
                Jv = czr <= 0.0 ? ntuple(_ -> 0.0, EASLEE.NJ_CB05) :
                     EASLEE.photo_j(cfg.photolut, acosd(clamp(czr, 0.0, 1.0)),
                                    7.64 * log(101325.0 / max(cs.P_Pa[i], 1.0)),
                                    cs.o3col[i], cs.albedo[i]; T_K=cs.T_K[i])
            end
            EASLEE.init_parcel!(view(cs.chem, i, :), ci, lon, lat, p_hPa;
                                T=cs.T_K[i], h2o_ppb=cs.h2o[i],
                                t_abs=t_abs, doy=doy, J=Jv)
            # carried aerosol load from the current bracket slot (aermr vary
            # slowly over one 3-h slot; seeding happens once per parcel)
            isempty(amap) ||
                _seed_aero!(i, EASLEE.interp_cams(obj.chembg, lon, lat, p_hPa, 0.0))
            # release-scoped plume IC overlays the background for eruption parcels
            if obj.relchem !== nothing
                ov = get(obj.relchem, Int(state.relbox[i]), nothing)
                ov === nothing || EASLEE.apply_release_chem!(cs, ov, i)
            end
            if obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor
                EASLEE.seed_h2o_from_plume!(cs, state.mass[i], state.airm[i], i)
            end
            # freeze the passive twin == the fully initialised composition
            isempty(cs.chem_pass) || (@views cs.chem_pass[i, :] .= cs.chem[i, :])
            # a recycled slot must start its tendency accumulation from zero
            isempty(cs.tend_cols) || (@views cs.tend_chem[i, :] .= 0.0)
            # plume-origin (tropo) air tag = 1.0 on any volcanic-release parcel
            if !isempty(cs.tropo)
                isplume = (obj.relchem !== nothing && haskey(obj.relchem, Int(state.relbox[i]))) ||
                          (obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor)
                isplume && (cs.tropo[i] = 1.0)
            end
            seeded[i] = true
        end
    elseif obj.chembg !== nothing
        # legacy path (no [chemistry.init] block): one frozen slot, direct pairs
        pairs = _cams_chem_pairs(obj.chembg, cfg.mechanism)
        tfrac = 0.0
        @inbounds for i in eachindex(state.status)
            (state.status[i] == 0x01 && !seeded[i]) || continue
            lat, lon = cs.lat[i], cs.lon[i]
            p_hPa = cs.P_Pa[i] / 100.0
            vals = EASLEE.interp_cams(obj.chembg, lon, lat, p_hPa, tfrac)
            for (sn, col) in pairs
                cs.chem[i, col] = max(Float64(vals[EASLEE.species_index(obj.chembg, sn)]), 0.0)
            end
            _seed_aero!(i, vals)
            if obj.relchem !== nothing
                ov = get(obj.relchem, Int(state.relbox[i]), nothing)
                ov === nothing || EASLEE.apply_release_chem!(cs, ov, i)
            end
            if obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor
                EASLEE.seed_h2o_from_plume!(cs, state.mass[i], state.airm[i], i)
            end
            isempty(cs.chem_pass) || (@views cs.chem_pass[i, :] .= cs.chem[i, :])
            isempty(cs.tend_cols) || (@views cs.tend_chem[i, :] .= 0.0)
            # plume-origin (tropo) air tag = 1.0 on any volcanic-release parcel
            if !isempty(cs.tropo)
                isplume = (obj.relchem !== nothing && haskey(obj.relchem, Int(state.relbox[i]))) ||
                          (obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor)
                isplume && (cs.tropo[i] = 1.0)
            end
            seeded[i] = true
        end
    elseif obj.relchem !== nothing || obj.relh2oplume !== nothing
        # no CAMS / declarative background — plume IC on a bare mechanism state
        @inbounds for i in eachindex(state.status)
            (state.status[i] == 0x01 && !seeded[i]) || continue
            if obj.relchem !== nothing
                ov = get(obj.relchem, Int(state.relbox[i]), nothing)
                ov === nothing || EASLEE.apply_release_chem!(cs, ov, i)
            end
            if obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor
                EASLEE.seed_h2o_from_plume!(cs, state.mass[i], state.airm[i], i)
            end
            isempty(cs.chem_pass) || (@views cs.chem_pass[i, :] .= cs.chem[i, :])
            isempty(cs.tend_cols) || (@views cs.tend_chem[i, :] .= 0.0)
            # plume-origin (tropo) air tag = 1.0 on any volcanic-release parcel
            if !isempty(cs.tropo)
                isplume = (obj.relchem !== nothing && haskey(obj.relchem, Int(state.relbox[i]))) ||
                          (obj.relh2oplume !== nothing && Int(state.relbox[i]) in obj.relh2oplume.vapor)
                isplume && (cs.tropo[i] = 1.0)
            end
            seeded[i] = true
        end
    end
    end  # timeit! :cams_reseed

    # integrate over the chemistry operator-split window (chemistry.timestep,
    # expressed as an integer number of transport steps)
    ae = max(1, round(Int, spec.chemistry.timestep_s / dt_s))
    if itime % ae == 0
        # surface emissions into the domain-fill parcels (AFTER CAMS seeding so
        # the increment is not overwritten, BEFORE the mechanism integration)
        if emis !== nothing && dfill !== nothing
            EASLEE.timeit!(obj.bench, :emissions) do
                EASLEE.apply_emissions!(cs, state, emis, gc, dfill.mass_per_parcel,
                                        ae * dt_s, kinject; when=currentdate)
            end
        end
        # lightning NOx into the UT parcels of convective columns ([chemistry.lightning])
        if obj.lightning !== nothing && dfill !== nothing && obj.convfields !== nothing
            EASLEE.apply_lightning_nox!(obj.lightning, cs, state, rt.ctx, obj.convfields,
                dfill.mass_per_parcel, ae * dt_s,
                EASLEE.MetAdapter.met_time_index_array(currentdate, rt.ctx.nests))
        end
        # parcel–parcel composition exchange ([chemistry.mixing] = "pairwise"):
        # mix BEFORE the mechanism so freshly diluted plume parcels react at
        # their mixed state.  Only seeded parcels take part — an unseeded slot
        # holds zeros and would drain its partner.
        if spec.chemistry.mixing == "pairwise"
            EASLEE.timeit!(obj.bench, :chem_mixing) do
                EASLEE.apply_chem_mixing!(cs, state, ae * dt_s;
                                          tau    = spec.chemistry.mixing_tau_s,
                                          cells  = spec.chemistry.mixing_cells,
                                          levels = spec.chemistry.mixing_levels,
                                          mask   = seeded)
            end
        end
        EASLEE.timeit!(obj.bench, :chemistry) do
            EASLEE.apply_chemistry!(cs, ae * dt_s, t_abs0 + itime * dt_s, cfg)
        end
    end
    return nothing
end

# ── Inter-trajectory mixing during transport ([physics.mixing] = "pairwise") ─
# Pairwise Fickian exchange of tracer MASS between co-located parcels of the
# SAME species (chemistry/mixing.jl old scheme).  Within one release every
# parcel represents an equal air mass, so exchanging mass ≡ exchanging mixing
# ratio (washout-scaled per-parcel masses ARE genuine q differences).  The
# coverage-based background term (q_bg = 0) entrains zero-tracer ambient air:
# parcels fully surrounded by plume neighbours (cov ≈ 1) only exchange among
# themselves; edge / isolated parcels dilute at rate (1−cov)·K_h/σ₀².
# Σmass is conserved by the exchange; the background term is the physical
# leak of tracer into untracked ambient air.
function apply_transport_mixing!(state, ctx, oc, mx, indt, dt::Float64; tropo_vec=nothing)
    N   = length(state.status)
    gc  = ctx.nests[1].gridconfig
    m3  = ctx.nests[1].metfields3D
    dx_m, dy_m = EASLEE._grid_spacing_m(gc, state)
    σ0  = mx.sigma0_h > 0 ? mx.sigma0_h : max(dx_m, dy_m)
    mcfg = EASLEE.MixingConfig(sigma0_h = σ0, K_h = mx.K_h, nsig = mx.nsig,
                               max_neighbors = mx.max_neighbors)
    x  = Vector{Float64}(undef, N); y  = similar(x)
    z  = Vector{Float64}(undef, N); yl = similar(x)
    q  = Matrix{Float64}(undef, N, 1)
    on = falses(N)
    @inbounds for i in 1:N
        state.status[i] == 0x01 || continue
        state.relbox[i] > 0     || continue      # boundary/domainfill parcels: skip
        on[i] = true
        x[i]  = Float64(state.indx[i]) * dx_m
        y[i]  = Float64(state.indy[i]) * dy_m
        z[i]  = Float64(EASLEE.Operators.interpolation4Dlevsingle(
                    m3.zlev, state.indx[i], state.indy[i],
                    state.indl[i], Float32(indt)))
        yl[i] = mx.sigma_v > 0 ? mx.sigma_v : max(Float64(state.miss3[i]), 10.0)
    end
    q_bg = mx.background ? [0.0] : nothing
    npairs = 0
    act = falses(N)
    for sn in sort(unique(Int.(oc.boxspecies)))
        fill!(act, false)
        n_s = 0
        @inbounds for i in 1:N
            on[i] || continue
            if Int(oc.boxspecies[Int(state.relbox[i])]) == sn
                act[i] = true; q[i, 1] = Float64(state.mass[i]); n_s += 1
            end
        end
        n_s < 2 && continue
        # merge_packed=false: point-release parcels are legitimately dense and
        # the cluster merger would relocate the LOCAL x/y copies mid-loop,
        # corrupting the next species' pair search
        npairs += EASLEE.apply_mixing_xyz!(q, x, y, z, yl, act, dt, mcfg;
                                           q_bg = q_bg, merge_packed = false)
        F = eltype(state.mass)
        @inbounds for i in 1:N
            act[i] && (state.mass[i] = F(q[i, 1]))
        end
    end
    # plume-origin (tropo) air-fraction tag: mix as a mixing ratio across ALL
    # active plume parcels together (one air mass, both releases), always
    # diluting toward 0 (ambient stratosphere carries no tropo air) — this is
    # the entrainment signal that makes f_tropo decay as the plume mixes out.
    if tropo_vec !== nothing && !isempty(tropo_vec)
        fill!(act, false); n_t = 0
        @inbounds for i in 1:N
            on[i] || continue
            act[i] = true; q[i, 1] = tropo_vec[i]; n_t += 1
        end
        if n_t >= 2
            EASLEE.apply_mixing_xyz!(q, x, y, z, yl, act, dt, mcfg;
                                     q_bg = [0.0], merge_packed = false)
            @inbounds for i in 1:N
                act[i] && (tropo_vec[i] = q[i, 1])
            end
        end
    end
    return npairs
end

# Bin the per-parcel transport-history tags onto the output (lon,lat,height)
# grid → cell-mean fields (NaN where empty), a snapshot at the frame boundary
# (like AOD).  Same cell mapping as chemconc_sample!/conccalc!.  One pass fills
# all four fields: strat_frac / pbl_frac / cont_frac ∈ [0,1] and conv_count ≥ 0.
function _grid_tag_fields(tags, state, ctx, oc, t_array; nest::Int = 1)
    Ni, Nj, Nk = Int(oc.Niout), Int(oc.Njout), Int(oc.Nkout)
    ssum = zeros(Float64, Ni, Nj, Nk); psum = zeros(Float64, Ni, Nj, Nk)
    csum = zeros(Float64, Ni, Nj, Nk); vsum = zeros(Float64, Ni, Nj, Nk)
    cnt  = zeros(Int, Ni, Nj, Nk)
    nestctx = ctx.nests[nest]; g = nestctx.gridconfig; lc = nestctx.levelconfig
    zlev = nestctx.metfields3D.zlev; zi = oc.outinterfaces
    indt = Float32(t_array[nest])
    @inbounds for i in eachindex(state.status)
        (state.status[i] == 0x01 && tags.n_steps[i] > 0) || continue
        indx = state.indx[i]; indy = state.indy[i]; indl = state.indl[i]
        ii = Base.unsafe_trunc(Int, EASLEE.Operators.indxtoxout(g, indx, indy))
        jj = Base.unsafe_trunc(Int, EASLEE.Operators.indytoyout(g, indx, indy))
        (1 <= ii <= Ni && 1 <= jj <= Nj) || continue
        zout = EASLEE.Operators.indltozlev(lc, zlev, indx, indy, indl, indt)
        zout < 0f0 && continue
        l = EASLEE.Operators.find_out_layer(zout, zi); l == 0 && continue
        ssum[ii, jj, l] += Float64(EASLEE.Operators.strat_fraction(tags, i))
        psum[ii, jj, l] += Float64(EASLEE.Operators.pbl_fraction(tags, i))
        csum[ii, jj, l] += Float64(EASLEE.Operators.cont_fraction(tags, i))
        vsum[ii, jj, l] += Float64(tags.conv_count[i])
        cnt[ii, jj, l]  += 1
    end
    strat = fill(NaN32, Ni, Nj, Nk); pbl = fill(NaN32, Ni, Nj, Nk)
    cont  = fill(NaN32, Ni, Nj, Nk); conv = fill(NaN32, Ni, Nj, Nk)
    @inbounds for idx in eachindex(strat)
        if cnt[idx] > 0
            c = cnt[idx]
            strat[idx] = Float32(ssum[idx] / c); pbl[idx] = Float32(psum[idx] / c)
            cont[idx]  = Float32(csum[idx] / c); conv[idx] = Float32(vsum[idx] / c)
        end
    end
    return (strat = strat, pbl = pbl, cont = cont, conv = conv)
end

# ── Time-framed grid output (FLEXPART LOUTSTEP/LOUTAVER/LOUTSAMPLE) ──────────
# Streamed per frame: each frame is frozen from the running accumulators and
# IMMEDIATELY appended as one time slice to the NetCDF (unlimited time dim), so
# no full time series is ever resident — only one per-frame scratch buffer.
mutable struct GridWriter
    path      :: String
    lon       :: Vector{Float64}
    lat       :: Vector{Float64}
    hgt       :: Vector{Float64}
    start     :: DateTime
    spec      :: Any
    has_conc  :: Bool
    has_mr    :: Bool
    has_res   :: Bool
    has_aod   :: Bool                       # sulfate/aerosol column AOD550 (output.aod)
    has_tags  :: Bool                       # transport-history tag profile ([tracers])
    gspec_ids :: Vector{Int}                # species numbers for per-species frames
    spnames   :: Dict{Int,String}
    cbuf      :: Array{Float32,3}           # conc scratch (Ni,Nj,Nk); 0 where empty
    mbuf      :: Array{Float32,3}           # mr scratch
    rbuf      :: Array{Float32,2}           # resolution scratch (Ni,Nj); NaN where empty
    aodbuf    :: Array{Float32,2}           # AOD550 scratch (Ni,Nj); 0-size when off
    tagbuf    :: Array{Float32,3}           # strat_frac scratch (Ni,Nj,Nk); 0-size when off
    tagpbl    :: Array{Float32,3}           # pbl_frac scratch  (Ni,Nj,Nk); 0-size when off
    tagcont   :: Array{Float32,3}           # cont_frac scratch (Ni,Nj,Nk); 0-size when off
    tagconv   :: Array{Float32,3}           # conv_count scratch (Ni,Nj,Nk); 0-size when off
    sbuf      :: Dict{Int,Array{Float32,3}} # per-species scratch
    created   :: Bool
end

function GridWriter(oc, path::AbstractString, start::DateTime, spec;
                    has_conc, has_mr, has_res, has_aod = false, has_tags = false,
                    gspec_ids, spnames)
    Ni, Nj, Nk = Int(oc.Niout), Int(oc.Njout), Int(oc.Nkout)
    lon = Float64[oc.lon0 + (i - 0.5) * oc.dlon for i in 1:Ni]
    lat = Float64[oc.lat0 + (j - 0.5) * oc.dlat for j in 1:Nj]
    zi  = oc.outinterfaces
    hgt = Float64[k == 1 ? Float64(zi[1]) : 0.5 * (Float64(zi[k-1]) + Float64(zi[k])) for k in 1:Nk]
    ids = collect(Int, gspec_ids)
    GridWriter(String(path), lon, lat, hgt, start, spec,
               has_conc, has_mr, has_res, has_aod, has_tags, ids, spnames,
               has_conc ? zeros(Float32, Ni, Nj, Nk) : zeros(Float32, 0,0,0),
               has_mr   ? zeros(Float32, Ni, Nj, Nk) : zeros(Float32, 0,0,0),
               has_res  ? fill(NaN32, Ni, Nj)        : fill(NaN32, 0,0),
               has_aod  ? zeros(Float32, Ni, Nj)     : zeros(Float32, 0,0),
               has_tags ? fill(NaN32, Ni, Nj, Nk)    : fill(NaN32, 0,0,0),
               has_tags ? fill(NaN32, Ni, Nj, Nk)    : fill(NaN32, 0,0,0),
               has_tags ? fill(NaN32, Ni, Nj, Nk)    : fill(NaN32, 0,0,0),
               has_tags ? fill(NaN32, Ni, Nj, Nk)    : fill(NaN32, 0,0,0),
               Dict(n => zeros(Float32, Ni, Nj, Nk) for n in ids),
               false)
end

# Create the grid NetCDF (once, when the first frame is frozen): conc/mr/<id>_conc
# on (lon,lat,height,time) + resolution(lon,lat,time), all with an UNLIMITED time
# dimension so frames append one slice at a time.  Float64 on disk (matches the
# former write-once output byte-for-byte).
function _grid_create!(gw::GridWriter)
    path = gw.path
    isfile(path) && rm(path)
    lon, lat, hgt = gw.lon, gw.lat, gw.hgt
    gw.has_conc && nccreate(path, "conc", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "kg m-3",
                    "long_name" => "concentration, mean over the trailing output.average window"))
    gw.has_mr && nccreate(path, "mr", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "kg kg-1",
                    "long_name" => "mass mixing ratio, mean over the trailing output.average window"))
    for n in gw.gspec_ids
        id = get(gw.spnames, n, "spec" * string(n))
        nccreate(path, id * "_conc", "lon", lon, "lat", lat, "height", hgt, "time", 0;
            atts = Dict("units" => "kg m-3",
                        "long_name" => id * " concentration, mean over the trailing output.average window"))
    end
    gw.has_res && nccreate(path, "resolution", "lon", lon, "lat", lat, "time", 0;
        atts = Dict("units" => "m",
                    "long_name" => "R(x): mass-weighted mean adaptive-kernel bandwidth — the true resolution of the gridded field (R < grid spacing: resolved; R > grid spacing: raise N)"))
    gw.has_aod && nccreate(path, "AOD550", "lon", lon, "lat", lat, "time", 0;
        atts = Dict("units" => "1",
                    "long_name" => "aerosol optical depth at 550 nm (column, hygroscopic-grown; snapshot at each frame boundary)"))
    gw.has_tags && nccreate(path, "strat_frac", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "1", "_FillValue" => NaN,
                    "long_name" => "stratospheric-residence fraction ([tracers]): cell-mean, over parcels, of the fraction of each trajectory spent above θ=strat_theta"))
    gw.has_tags && nccreate(path, "pbl_frac", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "1", "_FillValue" => NaN,
                    "long_name" => "PBL-contact fraction ([tracers]): cell-mean fraction of each trajectory spent with z_agl ≤ pblh"))
    gw.has_tags && nccreate(path, "cont_frac", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "1", "_FillValue" => NaN,
                    "long_name" => "continental fraction ([tracers]): cell-mean share of PBL-contact time over land (lsm ≥ 0.5); marine = 1 − cont_frac"))
    gw.has_tags && nccreate(path, "conv_count", "lon", lon, "lat", lat, "height", hgt, "time", 0;
        atts = Dict("units" => "1", "_FillValue" => NaN,
                    "long_name" => "convective displacements ([tracers]): cell-mean count of deep-convective redistributions per parcel"))
    nccreate(path, "time", "time", 0; atts = Dict("units" => "seconds since $(gw.start)",
                                                  "long_name" => "end of each averaging window"))
    NetCDF.ncputatt(path, "global", Dict(
        "start" => string(gw.start), "interval_s" => gw.spec.output.interval_s,
        "average_s" => gw.spec.output.average_s, "sample_s" => gw.spec.output.sample_s))
    gw.created = true
    return nothing
end

# Freeze the running accumulators into the per-frame scratch — conc = Σm/(nsamp·V)
# [kg m⁻³] (massarray, summed over boxes) and mr = Σ(m/ρ)/(nsamp·V) [kg kg⁻¹]
# (concarray) — then append that one frame as time slice `t` and reset the
# resolution accumulator.  Empty conc/mr cells are 0; empty resolution is NaN.
function _freeze_grid_frame!(gw::GridWriter, oc, nsamp::Int, t::Int, t_sec::Float64;
                             aod_field = nothing, tag_field = nothing)
    ns  = max(nsamp, 1)
    vol = oc.volume
    Ni, Nj, Nk = size(vol)
    gw.has_conc && fill!(gw.cbuf, 0.0f0)
    gw.has_mr   && fill!(gw.mbuf, 0.0f0)
    for g in values(gw.sbuf); fill!(g, 0.0f0); end
    @inbounds for k in 1:Nk, j in 1:Nj, i in 1:Ni
        v = vol[i, j, k]
        (isfinite(v) && v > 0) || continue
        if gw.has_conc
            s = 0.0
            for b in axes(oc.massarray, 4); s += oc.massarray[i, j, k, b, 1]; end
            gw.cbuf[i, j, k] = Float32(s / (ns * v))
        end
        # per-species frames (output.per_species with >1 species): group the
        # output COLUMNS by their species number (oc.classspecies — one column
        # per box under "release" grouping, per species under "species")
        if !isempty(gw.sbuf)
            for b in axes(oc.massarray, 4)
                m = oc.massarray[i, j, k, b, 1]
                m == 0.0 && continue
                g = get(gw.sbuf, Int(oc.classspecies[b]), nothing)
                g === nothing && continue
                g[i, j, k] += Float32(m / (ns * v))
            end
        end
        if gw.has_mr
            s = 0.0
            for b in axes(oc.concarray, 4); s += Float64(oc.concarray[i, j, k, b, 1]); end
            gw.mbuf[i, j, k] = Float32(s / (ns * v))
        end
    end
    # R(x): mass-weighted mean kernel width over the window (NaN where empty)
    if gw.has_res
        fill!(gw.rbuf, NaN32)
        @inbounds for j in 1:Nj, i in 1:Ni
            gw.rbuf[i, j] = oc.resden[i, j] > 0.0 ?
                Float32(oc.resnum[i, j] / oc.resden[i, j]) : NaN32
        end
        fill!(oc.resnum, 0.0); fill!(oc.resden, 0.0)
    end

    gw.created || _grid_create!(gw)
    # append this frame as one time slice (grows the unlimited time dim to `t`)
    gw.has_conc && NetCDF.ncwrite(reshape(Float64.(gw.cbuf), Ni,Nj,Nk,1), gw.path, "conc";
                       start = [1,1,1,t], count = [Ni,Nj,Nk,1])
    gw.has_mr   && NetCDF.ncwrite(reshape(Float64.(gw.mbuf), Ni,Nj,Nk,1), gw.path, "mr";
                       start = [1,1,1,t], count = [Ni,Nj,Nk,1])
    for (n, g) in gw.sbuf
        id = get(gw.spnames, n, "spec" * string(n))
        NetCDF.ncwrite(reshape(Float64.(g), Ni,Nj,Nk,1), gw.path, id * "_conc";
                       start = [1,1,1,t], count = [Ni,Nj,Nk,1])
    end
    gw.has_res && NetCDF.ncwrite(reshape(Float64.(gw.rbuf), Ni,Nj,1), gw.path, "resolution";
                       start = [1,1,t], count = [Ni,Nj,1])
    if gw.has_aod
        if aod_field !== nothing
            @inbounds for j in 1:Nj, i in 1:Ni; gw.aodbuf[i, j] = Float32(aod_field[i, j]); end
        else
            fill!(gw.aodbuf, 0.0f0)
        end
        NetCDF.ncwrite(reshape(Float64.(gw.aodbuf), Ni,Nj,1), gw.path, "AOD550";
                       start = [1,1,t], count = [Ni,Nj,1])
    end
    if gw.has_tags
        Nk = size(gw.tagbuf, 3)
        # tag_field is the NamedTuple (strat,pbl,cont,conv) from _grid_tag_fields
        for (buf, key, var) in ((gw.tagbuf, :strat, "strat_frac"),
                                (gw.tagpbl, :pbl,  "pbl_frac"),
                                (gw.tagcont, :cont, "cont_frac"),
                                (gw.tagconv, :conv, "conv_count"))
            if tag_field !== nothing
                src = getproperty(tag_field, key)
                @inbounds for k in 1:Nk, j in 1:Nj, i in 1:Ni; buf[i,j,k] = Float32(src[i,j,k]); end
            else
                fill!(buf, NaN32)
            end
            NetCDF.ncwrite(reshape(Float64.(buf), Ni,Nj,Nk,1), gw.path, var;
                           start = [1,1,1,t], count = [Ni,Nj,Nk,1])
        end
    end
    NetCDF.ncwrite([t_sec], gw.path, "time"; start = [t], count = [1])
    return nothing
end

_grid_close!(gw::GridWriter) = (gw.created && NetCDF.ncclose(gw.path); nothing)

"""
    run_spec!(spec; prefetch=true, verbose=true) → (rt, obj, fname)

Run the whole simulation described by `spec`: `build_run`, then the time loop
(transport + physics, plus whatever the config enables — domainfill,
chemistry/emissions/mixing, receptors, covariance kernel, nudging, trajectory
CSV).  Writes the time-framed grid NetCDF `grid_<name>.nc` (one frame per
`output.interval`, FLEXPART LOUTSTEP/LOUTAVER semantics; single whole-run mean
for backward or sub-interval runs) and, with chemistry on, `chem_<name>.nc`.
Returns the runtime tuple, the `RunObjects`, and the grid file path.
"""
# ── Warm start ([run] checkpoint / restart) ──────────────────────────────────
# The checkpoint carries everything the loop cannot recreate: the full parcel
# state, per-parcel chemistry, the physics accumulators, and the date.  The
# release quota is per-step, so resuming inside/after a release window never
# double-releases; met/output are re-derived from the config.
const _CHK_VERSION = 1

_chk_fingerprint(spec, rt) = (spec.name, length(spec.species), length(spec.release),
                              length(rt.state.status))

function _checkpoint_save(path::AbstractString, rt, obj, spec, d::DateTime, seeded)
    st = rt.state
    chk = Dict{Symbol,Any}(
        :version => _CHK_VERSION, :date => d,
        :fingerprint => _chk_fingerprint(spec, rt),
        :seeded => copy(seeded))
    for f in fieldnames(typeof(st))
        f in (:nest_buffers,) && continue
        chk[Symbol("state_", f)] = copy(getfield(st, f))
    end
    if obj.chemstate !== nothing
        chk[:chem] = copy(obj.chemstate.chem)
        isempty(obj.chemstate.chem_pass) || (chk[:chem_pass] = copy(obj.chemstate.chem_pass))
        isempty(obj.chemstate.tend_cols) || (chk[:tend_chem] = copy(obj.chemstate.tend_chem))
        # SWV-prognostic mode: the parcel-carried water budget must survive
        # a warm restart (mock chemstates in tests lack these fields → guard)
        if hasproperty(obj.chemstate, :swv_init) && !isempty(obj.chemstate.swv_init)
            chk[:swv_init] = copy(obj.chemstate.swv_init)
            chk[:swv_h2o]  = copy(obj.chemstate.h2o)
            isempty(obj.chemstate.ch4) || (chk[:swv_ch4] = copy(obj.chemstate.ch4))
        end
        # sulfate SAD state (so4/so2/awet), plume-origin tag, dehyd/denit, CAMS aero —
        # NOT recreatable on restart (the SO2 release fired at t0) → must be saved.
        # (same guard as swv_init above: mock chemstates in tests carry only
        # `chem`, so every optional field must be probed before it is read)
        _cs = obj.chemstate
        _has(f) = hasproperty(_cs, f) && !isempty(getproperty(_cs, f))
        _has(:so4)   && (chk[:so4]   = copy(_cs.so4))
        _has(:so2)   && (chk[:so2]   = copy(_cs.so2))
        _has(:awet)  && (chk[:awet]  = copy(_cs.awet))
        _has(:tropo) && (chk[:tropo] = copy(_cs.tropo))
        _has(:dehyd) && (chk[:dehyd] = copy(_cs.dehyd))
        _has(:denit) && (chk[:denit] = copy(_cs.denit))
        hasproperty(_cs, :aero) && size(_cs.aero, 2) > 0 && (chk[:aero] = copy(_cs.aero))
        _has(:aod_above) && (chk[:aod_above] = copy(_cs.aod_above))
    end
    obj.icephys === nothing ||
        (chk[:ice_frozen] = copy(obj.icephys.frozen);
         chk[:ice_scavenged] = copy(obj.icephys.scavenged))
    obj.convphys === nothing ||
        (chk[:converted] = copy(obj.convphys.converted);
         chk[:cv_demoted] = copy(obj.convphys.demoted);
         chk[:cv_gone]    = copy(obj.convphys.gone);
         # per-parcel crystal radius: without it a warm start re-seeds every
         # crystal at its bin's NOMINAL radius and silently resets the cascade
         isempty(obj.convphys.rice) || (chk[:cv_rice] = copy(obj.convphys.rice)))
    obj.pradphys === nothing ||
        (chk[:trad] = copy(obj.pradphys.trad);
         # per-parcel T′: without it a warm start re-spins the cold anomaly
         # from zero (τ ≈ 10.8 h of wrong saturation gates)
         isempty(obj.pradphys.tanom) || (chk[:pr_tanom] = copy(obj.pradphys.tanom)))
    obj.tagphys === nothing ||
        (chk[:tag_ss] = copy(obj.tagphys.strat_steps);
         chk[:tag_ns] = copy(obj.tagphys.n_steps);
         chk[:tag_age] = copy(obj.tagphys.age_strat);
         chk[:tag_ps] = copy(obj.tagphys.pbl_steps);
         chk[:tag_ls] = copy(obj.tagphys.land_steps);
         chk[:tag_agp] = copy(obj.tagphys.age_pbl);
         chk[:tag_cv] = copy(obj.tagphys.conv_count))
    isempty(obj.ssprays) ||
        (chk[:spray] = [(copy(b.acc), b.emitted[], b.emitted_mean[], b.emitted_gusty[], b.nspawn[])
                        for b in obj.ssprays])
    tmp = path * ".tmp"
    open(io -> serialize(io, chk), tmp, "w")
    mv(tmp, path; force=true)
    @info "checkpoint saved" file=path date=d active=count(==(0x01), rt.state.status)
    return nothing
end

# ── Per-parcel diagnostic snapshot ([run] parcel_dump) ───────────────────────
# 2026-08-05 pre-registration, steps 1–3: the freeze/demotion attribution and
# the sectional-τ instrument both need per-parcel state at a time CADENCE the
# checkpoint (one file, overwritten) cannot provide.  Each snapshot is
# self-contained: parcel state + gate state + the cumulative flux accumulators
# (deltas between snapshots = fluxes per interval) + the box metadata needed
# to map box → species/bin offline.  Diagnostic only — reading it back is the
# analysis scripts' job, the model never loads these.
function _parcel_snapshot_save(path::AbstractString, rt, obj, d::DateTime)
    st = rt.state
    snap = Dict{Symbol,Any}(
        :date => d,
        :status => copy(st.status), :relbox => copy(st.relbox),
        :mass => copy(st.mass), :airm => copy(st.airm),
        :indl => copy(st.indl),
        # horizontal position (grid index) → cell-level sectional τ offline;
        # id → parcel tracking across snapshots (slots can be recycled)
        :indx => copy(st.indx), :indy => copy(st.indy), :id => copy(st.id))
    obj.hybridvert === nothing ||
        (snap[:theta_p] = copy(obj.hybridvert.theta_p))
    if obj.convphys !== nothing
        cv = obj.convphys
        isempty(cv.rice)   || (snap[:rice] = copy(cv.rice))
        if !isempty(cv.gate_T)
            snap[:gate_T] = copy(cv.gate_T); snap[:gate_qsi] = copy(cv.gate_qsi)
            snap[:gate_p] = copy(cv.gate_p)
        end
        snap[:converted] = copy(cv.converted)
        snap[:demoted]   = copy(cv.demoted)
        snap[:gone]      = copy(cv.gone)
        snap[:fz] = hcat(cv.fz_mT, cv.fz_mTp, cv.fz_mlq,
                         cv.fz_mlqsi, cv.fz_mlp, cv.fz_mdT)
        # box metadata: mode (1 freeze / 2 sublimate), the gate threshold, the
        # bin's nominal radius, ladder links — everything offline needs
        snap[:box_mode] = copy(cv.mode); snap[:box_rhi] = copy(cv.rhi_crit)
        snap[:box_r] = copy(cv.subl_radius); snap[:box_target] = copy(cv.target)
        snap[:box_down] = copy(cv.bin_down); snap[:box_redge] = copy(cv.bin_redge)
    end
    obj.pradphys === nothing || isempty(obj.pradphys.tanom) ||
        (snap[:tanom] = copy(obj.pradphys.tanom))
    # per-parcel kNN vertical spread h_v (update_airm!'s `hv_out`).  With airm
    # and the gate's own p/T this decomposes the closure:
    #   airm = rho * 4*pi*sqrt(3) * h_h^2 * h_v / k   ->   h_h is recoverable.
    # That is the test of whether a freezing parcel's airm contraction is a
    # real concentration (h_h and h_v shrink together) or the settling
    # artefact (h_v collapses alone, neighbours being ranked horizontally).
    if obj.icephys !== nothing && !isempty(obj.icephys.hv)
        snap[:hv] = copy(obj.icephys.hv)
    elseif obj.pradphys !== nothing && !isempty(obj.pradphys.hv)
        snap[:hv] = copy(obj.pradphys.hv)
    end
    tmp = path * ".tmp"
    open(io -> serialize(io, snap), tmp, "w")
    mv(tmp, path; force=true)
    return nothing
end

function _checkpoint_load!(rt, obj, spec, path::AbstractString, seeded)::DateTime
    chk = open(deserialize, path)
    chk[:version] == _CHK_VERSION ||
        error("checkpoint version $(chk[:version]) ≠ $_CHK_VERSION")
    fp = _chk_fingerprint(spec, rt)
    saved = chk[:fingerprint]
    (saved[1] == fp[1] && saved[2] == fp[2] && saved[3] == fp[3]) ||
        @warn "checkpoint config fingerprint differs — same-config continuation expected" saved=saved current=fp
    st = rt.state
    nsave = length(chk[:state_status])
    nsave <= length(st.status) ||
        error("checkpoint pool ($nsave) larger than this run's pool ($(length(st.status))) — keep the release/particle setup identical")
    for f in fieldnames(typeof(st))
        f in (:nest_buffers,) && continue
        key = Symbol("state_", f)
        haskey(chk, key) || continue
        copyto!(view(getfield(st, f), 1:nsave), chk[key])
    end
    copyto!(view(seeded, 1:nsave), chk[:seeded])
    if obj.chemstate !== nothing && haskey(chk, :chem)
        n = min(size(obj.chemstate.chem, 1), size(chk[:chem], 1))
        obj.chemstate.chem[1:n, :] .= chk[:chem][1:n, :]
        haskey(chk, :chem_pass) && !isempty(obj.chemstate.chem_pass) &&
            (obj.chemstate.chem_pass[1:n, :] .= chk[:chem_pass][1:n, :])
        haskey(chk, :tend_chem) && !isempty(obj.chemstate.tend_cols) &&
            (obj.chemstate.tend_chem[1:n, :] .= chk[:tend_chem][1:n, :])
    end
    if obj.chemstate !== nothing && haskey(chk, :swv_init) &&
       hasproperty(obj.chemstate, :swv_init) && !isempty(obj.chemstate.swv_init)
        n = length(chk[:swv_init])
        copyto!(view(obj.chemstate.swv_init, 1:n), chk[:swv_init])
        copyto!(view(obj.chemstate.h2o, 1:n), chk[:swv_h2o])
        haskey(chk, :swv_ch4) && !isempty(obj.chemstate.ch4) &&
            copyto!(view(obj.chemstate.ch4, 1:n), chk[:swv_ch4])
        # ch4_oxidation activé PAR-DESSUS un checkpoint qui ne le portait pas :
        # les cs.ch4 alloués restent à zéro et les parcelles restaurées (swv_init
        # pinné) ne seront jamais re-semées → Cl+CH4 à taux NUL, l'inverse du
        # but.  Backfill au CH4 initial pour tout slot restauré resté à zéro.
        if !haskey(chk, :swv_ch4) && !isempty(obj.chemstate.ch4)
            cs_ = obj.chemstate
            nfil = 0
            @inbounds for i in 1:n
                if cs_.swv_init[i] && cs_.ch4[i] == 0.0
                    cs_.ch4[i] = cs_.ch4_init; nfil += 1
                end
            end
            nfil > 0 && @info "restart: cs.ch4 backfilled to ch4_init (checkpoint predates ch4_oxidation)" parcels=nfil
        end
    end
    if obj.chemstate !== nothing            # sulfate SAD / tropo tag / dehyd-denit / aero
        _cs = obj.chemstate
        # probe by symbol: mock chemstates (tests) carry only `chem`, and old
        # checkpoints predate these keys — both sides must be optional
        function _rload!(k)
            (haskey(chk, k) && hasproperty(_cs, k)) || return
            v = getproperty(_cs, k)
            isempty(v) || copyto!(view(v, 1:length(chk[k])), chk[k])
        end
        _rload!(:so4);   _rload!(:so2);   _rload!(:awet);  _rload!(:tropo)
        _rload!(:dehyd); _rload!(:denit); _rload!(:aod_above)
        haskey(chk, :aero) && hasproperty(_cs, :aero) && size(_cs.aero, 2) > 0 &&
            (m = size(chk[:aero], 1); _cs.aero[1:m, :] .= chk[:aero][1:m, :])
    end
    obj.icephys === nothing || !haskey(chk, :ice_frozen) ||
        (obj.icephys.frozen .= chk[:ice_frozen];
         obj.icephys.scavenged .= chk[:ice_scavenged])
    if obj.convphys !== nothing
        haskey(chk, :converted)  && (obj.convphys.converted .= chk[:converted])
        haskey(chk, :cv_demoted) && (obj.convphys.demoted .= chk[:cv_demoted])
        haskey(chk, :cv_gone)    && (obj.convphys.gone .= chk[:cv_gone])
        haskey(chk, :cv_rice) && !isempty(obj.convphys.rice) &&
            copyto!(view(obj.convphys.rice, 1:length(chk[:cv_rice])), chk[:cv_rice])
        haskey(chk, :cv_rice) || isempty(obj.convphys.rice) ||
            @warn "restart: checkpoint predates cv_rice — sectional radii re-seed at NOMINAL bin radii (cascade state lost)"
    end
    if obj.pradphys !== nothing
        haskey(chk, :trad) &&
            copyto!(view(obj.pradphys.trad, 1:length(chk[:trad])), chk[:trad])
        haskey(chk, :pr_tanom) && !isempty(obj.pradphys.tanom) &&
            copyto!(view(obj.pradphys.tanom, 1:length(chk[:pr_tanom])), chk[:pr_tanom])
        haskey(chk, :pr_tanom) || isempty(obj.pradphys.tanom) ||
            @warn "restart: checkpoint predates pr_tanom — T′ re-spins from zero (τ of wrong gates)"
    end
    obj.tagphys === nothing || !haskey(chk, :tag_ss) || begin
        n = length(chk[:tag_ss])
        copyto!(view(obj.tagphys.strat_steps, 1:n), chk[:tag_ss])
        copyto!(view(obj.tagphys.n_steps, 1:n), chk[:tag_ns])
        copyto!(view(obj.tagphys.age_strat, 1:n), chk[:tag_age])
        if haskey(chk, :tag_ps)     # increment-2 tags (older checkpoints omit them)
            copyto!(view(obj.tagphys.pbl_steps, 1:n), chk[:tag_ps])
            copyto!(view(obj.tagphys.land_steps, 1:n), chk[:tag_ls])
            copyto!(view(obj.tagphys.age_pbl, 1:n), chk[:tag_agp])
            copyto!(view(obj.tagphys.conv_count, 1:n), chk[:tag_cv])
        end
    end
    if haskey(chk, :spray)
        for (b, sv) in zip(obj.ssprays, chk[:spray])
            length(b.acc) == length(sv[1]) && (b.acc .= sv[1])
            b.emitted[] = sv[2]; b.emitted_mean[] = sv[3]
            b.emitted_gusty[] = sv[4]; b.nspawn[] = sv[5]
        end
    end
    EASLEE.State.sync_dtime_ns!(st)
    @info "checkpoint restored" file=path date=chk[:date] active=count(==(0x01), st.status)
    return chk[:date]::DateTime
end

function run_spec!(spec::EasleeConfig.RunSpec; prefetch::Bool=true, verbose::Bool=true)
    rt, obj = build_run(spec; prefetch=prefetch)
    ctx = rt.ctx; rc = ctx.runconfig
    kw  = step_kwargs(obj)
    # [physics].prognostic_airm: only the populations whose q = mass/airm is
    # actually consumed (ice closure / plume radiation) — the domain-fill
    # parcels keep the air mass they were released with.
    airm_boxes = nothing
    if spec.physics.prognostic_airm
        bx = Set{Int}()
        obj.icephys === nothing || (for b in eachindex(obj.icephys.on)
            obj.icephys.on[b] && push!(bx, b) end)
        obj.pradphys === nothing || (for b in eachindex(obj.pradphys.on)
            (obj.pradphys.on[b] || obj.pradphys.partner[b]) && push!(bx, b) end)
        isempty(bx) || (airm_boxes = bx)
        @info "prognostic air mass armed ([physics].prognostic_airm)" boxes=length(bx) every_steps=spec.physics.prognostic_airm_every
    end
    ts  = rc.timestep
    dt_s = Float64(abs(Dates.value(ts))) / 1000.0
    loopstart = rc.direction == -1 ? rc.enddate : rc.startdate
    nsteps = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(ts)))
    # chemistry bookkeeping (no-op when obj.chemstate === nothing)
    do_chem = obj.chemstate !== nothing
    seeded  = falses(length(rt.state.status))
    # ── warm start ([run] restart): resume at the saved step ────────────────
    it0 = 0
    if !isempty(spec.restart)
        rc.direction == -1 && error("[run].restart: backward runs not supported")
        isfile(spec.restart) || error("[run].restart: file not found: $(spec.restart)")
        chkdate = _checkpoint_load!(rt, obj, spec, spec.restart, seeded)
        it0 = round(Int, Dates.value(chkdate - loopstart) / abs(Dates.value(ts))) + 1
        (0 < it0 <= nsteps) ||
            error("[run].restart: checkpoint date $chkdate outside this run's window")
        # never overwrite the original run's output files
        base = isempty(spec.output.name) ? spec.name : spec.output.name
        spec.output.name = base * Dates.format(chkdate, "_rmmddHH")
        @info "warm start: resuming" from=chkdate step=it0 of=nsteps output_suffix=spec.output.name
    end
    chk_period_s = spec.checkpoint in ("", "end") ? 0.0 :
                   EasleeConfig.ConfigTOML.parse_duration_s(spec.checkpoint)
    chk_path = joinpath(spec.output.path,
                        "checkpoint_" * (isempty(spec.name) ? "run" : spec.name) * ".jls")
    snap_period_s = isempty(spec.parcel_dump) ? 0.0 :
                    EasleeConfig.ConfigTOML.parse_duration_s(spec.parcel_dump)
    snap_dir = joinpath(spec.output.path, "snapshots")
    snap_period_s > 0 && (mkpath(snap_dir);
        @info "parcel snapshots armed ([run].parcel_dump)" dir=snap_dir period=spec.parcel_dump)
    t_abs0  = Dates.value(Dates.Time(loopstart)) / 1e9   # s since UTC midnight
    # domainfill release (kind="domainfill"): place the fill population at t=0,
    # then per-step box-edge inflow injection + cull of parcels leaving the box.
    dfill = nothing
    # Fixed-seed RNG for domain-fill parcel placement so a run is reproducible
    # run-to-run (and across thread counts — placement is serial).  Julia's
    # default RNG is entropy-seeded per process, which otherwise reshuffles every
    # parcel each launch.  Created once; the SAME stream feeds init + every inject.
    df_rng = Random.Xoshiro(0x0EA51EE5)
    dfrel = filter(r -> r.kind == "domainfill", spec.release)
    if !isempty(dfrel)
        length(dfrel) == 1 || error("run_spec!: only one domainfill release per run")
        r = dfrel[1]
        read_met_step!(rt, loopstart)
        indt0 = EASLEE.MetAdapter.met_time_index_array(loopstart, ctx.nests)[1]
        dfill = EASLEE.build_domainfill(ctx.nests[1], indt0, r.lon, r.lat, r.z[end],
                                        r.particles)
        placed = it0 > 0 ? 0 :   # warm start: the fill population is restored
                 EASLEE.domainfill_init!(rt.state, ctx, dfill, loopstart, indt0;
                                         rng=df_rng, seeded=seeded)
        @info "domainfill" placed levels=length(dfill.keep_lev) mass_per_parcel=dfill.mass_per_parcel box="$(dfill.lon0)–$(dfill.lon1)E $(dfill.lat0)–$(dfill.lat1)N ≤$(round(Int, dfill.ztop)) m"
    end
    # surface/plume emissions ([chemistry.emissions]) into the domain-fill parcels —
    # gases (needs chem) and/or aerosol mass (needs an emission-driven AerosolModel)
    emis = nothing; kinject = 0
    aero_emit = obj.aerosol !== nothing && !isempty(obj.aerosol.pmass)
    if (do_chem || aero_emit) && !isempty(spec.chemistry.emissions_dir)
        dfill === nothing &&
            error("chemistry.emissions requires a [[release]] kind=\"domainfill\"")
        dcyc = EASLEE.DiurnalCycle(Symbol(spec.chemistry.emissions_diurnal);
                   peak=spec.chemistry.emissions_diurnal_peak,
                   sigma=spec.chemistry.emissions_diurnal_sigma)
        emis = EASLEE.load_emissions(spec.chemistry.emissions_dir,
                   spec.chemistry.emissions_files,
                   do_chem ? obj.chemcfg.mechanism : nothing, spec.start;
                   lon0=dfill.lon0, lon1=dfill.lon1, lat0=dfill.lat0, lat1=dfill.lat1,
                   res=spec.chemistry.emissions_res,
                   nox_no_frac=spec.chemistry.emissions_nox_no, diurnal=dcyc)
        if emis !== nothing
            spec.chemistry.emissions_diurnal == "none" ||
                @info "emissions diurnal cycle" scheme=spec.chemistry.emissions_diurnal peak_LST=spec.chemistry.emissions_diurnal_peak
            ks = [k for k in eachindex(dfill.lev_h) if dfill.lev_h[k] <= spec.chemistry.emissions_inject_m]
            kinject = isempty(ks) ? 1 : maximum(ks)
            @info "emissions: injecting below $(spec.chemistry.emissions_inject_m) m (model levels 1–$kinject)"
        end
    end
    # gridded chemistry output: one frame per output interval (uses the sync-step
    # sampling cadence — chemconc_sample! each step, averaged per interval)
    cc = nothing; out_n = 1; nframe = 0
    # [output.reference]: compare the reacting chemistry to its init source.
    # ref_on needs the passive twin (allocated in ChemState) + a live source
    # bracket (obj.chemtl) to co-sample X_ref per parcel.  ref_cmap maps each
    # chem column to the source's vector index (built once from the bracket).
    ref_on   = !isempty(spec.output.reference) && do_chem &&
               obj.chemtl !== nothing && !isempty(obj.chemstate.chem_pass)
    ref_cmap = Tuple{Int,Int}[]
    if ref_on
        # mechanism species → source variable index.  Species names (O3, NO2)
        # are NOT the source's variable names (CAMS shortNames go3, no2) — resolve
        # via the source's own `vars` override, else the mechanism default map
        # (_default_chem_varmap: CB05_FROM_CAMS for cb05), exactly as
        # build_chem_init does for the default fill.
        mech = obj.chemcfg.mechanism
        src  = get(obj.cheminit.sources, spec.output.reference, nothing)
        vm   = _default_chem_varmap(mech)
        specn = collect(EASLEE.chem_species(mech))
        for (col, name) in enumerate(specn)
            v  = src === nothing ? get(vm, name, name) : get(src.vars, name, get(vm, name, name))
            vi = EASLEE.species_index(obj.chemtl.bg0, v)
            vi > 0 && push!(ref_cmap, (col, vi))
        end
        @info "output.reference: passive twin + $(length(ref_cmap))-species source comparison" source=spec.output.reference
    end
    if do_chem
        out_n  = max(1, round(Int, spec.output.interval_s / dt_s))
        NT     = max(1, nsteps ÷ out_n)
        cc_name = isempty(spec.output.name) ? spec.name : spec.output.name
        cc     = EASLEE.ChemConc(ctx.output,
                     collect(EASLEE.chem_species(obj.chemcfg.mechanism)), NT;
                     passive = ref_on, reference = ref_on,
                     aod = spec.chemistry.aerosols && spec.chemistry.aerosols_aod,
                     h2o = spec.chemistry.h2o_prognostic,
                     sad = spec.chemistry.sulfate_sad,
                     tropo = spec.physics.tags_tropo_origin && !isempty(obj.chemstate.tropo),
                     tendency = copy(obj.chemstate.tend_names),
                     path = joinpath(spec.output.path, "chem_$(cc_name).nc"))
    end
    # per-trajectory CSV (no-op unless [output.trajectory].enabled)
    traj = EASLEE.traj_open(spec, rt.state, dt_s)
    # ── time-framed grid output (forward runs with ≥ 1 full interval) ─────────
    # FLEXPART semantics: a frame is frozen at every output.interval, holding the
    # mean over the trailing output.average window sampled every output.sample.
    # The legacy single whole-run average is kept for backward runs (the SRR grid
    # is a time integral) and for runs shorter than one interval.
    oc      = ctx.output
    out_ng  = max(1, round(Int, spec.output.interval_s / dt_s))
    avg_ng  = max(1, round(Int, spec.output.average_s  / dt_s))
    samp_ng = max(1, round(Int, spec.output.sample_s   / dt_s))
    NTg     = spec.direction == "forward" ? nsteps ÷ out_ng : 0
    framed  = NTg >= 1
    # per-species framed grids (output.per_species with >1 species): the total
    # `conc` sums the release boxes, which drowns a 1e9-kg species under a
    # 3.3e12-kg one (HT HCl vs H2O) — write <id>_conc alongside
    spnum   = _species_numbers(spec)
    multisp = framed && oc.iout in (1, 3) && spec.output.per_species &&
              length(spec.species) > 1
    gspec_ids = multisp ? Int[Int(n) for (id, n) in spnum] : Int[]
    spnames   = multisp ? Dict(Int(n) => String(id) for (id, n) in spnum) : Dict{Int,String}()
    # streamed grid writer (created lazily on the first frame; frames appended
    # one time slice at a time — no full time series held in RAM)
    grid_out  = isempty(spec.output.name) ? spec.name : spec.output.name
    gw = framed ? GridWriter(oc, joinpath(spec.output.path, "grid_$(grid_out).nc"),
                             loopstart, spec;
                             has_conc = oc.iout in (1, 3), has_mr = oc.iout in (2, 3),
                             has_res = oc.resolution,
                             has_aod = spec.output.aod && obj.aerosol !== nothing,
                             has_tags = obj.tagphys !== nothing,
                             gspec_ids = gspec_ids,
                             spnames = spnames) : nothing
    # gridless receptors (KDE point evaluation, any grid kernel)
    recs = nothing
    if !isempty(spec.output.receptors) || !isempty(spec.output.receptors_file)
        rl = spec.output.receptors_layer
        # Optional bandwidth floor (receptors.jl).  DEFAULT 0 — the kernel
        # width must come from the trajectories themselves (mean distance to
        # the surrounding parcels = the dispersion the random walk carries in
        # centroid spread instead of per-parcel covariance).  min_bandwidth
        # is an explicit representativeness experiment, never silent.
        hmin = max(spec.output.receptors_hmin, 0.0)
        # NOTE: the dispersion bandwidth_floor is DELIBERATELY NOT wired into the
        # receptor (fsv=0).  The grid deposit renormalises each parcel's weights
        # to conserve its mass, so widening h just redistributes fixed mass over
        # more cells (neutral).  The receptor point-KDE does NOT renormalise
        # (C = Σ mᵢ·G₂ integrates to mᵢ over ALL space), so flooring h wider near
        # a source — where parcel density is high — makes every station sum over
        # many more parcels' tails and OVER-COUNTS (ETEX etex_floor A/B, 2026-07-07:
        # receptor FB +0.29→+1.63, 8021 conc ×350).  The receptor's
        # representativeness knob is the fixed `min_bandwidth`; the age-growing
        # per-parcel footprint fix is σ(age)/√N_local (TODO §F), not this floor.
        fsv  = 0.0
        ftl  = spec.output.floor_tl_s
        # kernel="covariance": the receptor evaluates the flow-evolved P
        # (full bivariate Gaussian) instead of the isotropic hband KDE
        ucov = spec.output.kernel == "covariance"
        ushp = ucov && spec.output.cov_mode == "shape"
        recs = if !isempty(spec.output.receptors_file)
            EASLEE.Operators.receptorset_from_file(spec.output.receptors_file,
                                                   rl[1], rl[2]; hmin=hmin, fsv=fsv, ftl=ftl,
                                                   usecov=ucov, shape=ushp)
        else
            EASLEE.Operators.ReceptorSet([r.id for r in spec.output.receptors],
                                         [r.lon for r in spec.output.receptors],
                                         [r.lat for r in spec.output.receptors],
                                         rl[1], rl[2]; hmin=hmin, fsv=fsv, ftl=ftl,
                                         usecov=ucov, shape=ushp)
        end
        @info "gridless receptors" stations=length(recs.ids) layer=(rl[1], rl[2]) hmin_m=round(hmin, digits=0) floor_sigma_v=fsv covariance=ucov
    end
    nframe_g = 0; nsamp_win = 0
    # Receptors read per-parcel bandwidths (state.hband).  When the GRID kernel
    # is adaptive, conccalc! already refreshes them; otherwise (uniform/none
    # grid + receptors) refresh here once per output frame so the gridless KDE
    # still has real bandwidths instead of the h→1 m floor.
    recs_need_bw = recs !== nothing &&
                   !(ctx.output.kernel in (:adaptive, :covariance, :covshape))
    last_bw_frame = Ref(-1)
    # observed-track nudging + wind-bias drift (HT thread; TOML-driven)
    tnudge = isempty(spec.physics.nudge_file) ? nothing :
        EASLEE.Operators.load_track_nudge(spec.physics.nudge_file,
                                          spec.physics.nudge_tau_s, spec.physics.nudge_z)
    tnudge === nothing ||
        @info "track nudging armed" file=basename(spec.physics.nudge_file) tau_s=spec.physics.nudge_tau_s z=spec.physics.nudge_z rows=length(tnudge.t_ms)
    wbias = (spec.physics.windbias_u == 0.0 && spec.physics.windbias_v == 0.0) ? nothing :
        EASLEE.Operators.WindBias(spec.physics.windbias_u, spec.physics.windbias_v,
                                  spec.physics.windbias_z[1], spec.physics.windbias_z[end])
    wbias === nothing ||
        @info "wind-bias drift armed" u_ms=spec.physics.windbias_u v_ms=spec.physics.windbias_v z=spec.physics.windbias_z
    # inter-trajectory mixing during transport ([physics.mixing]; TOML-driven)
    pmix = spec.physics.mixing.scheme == "pairwise" ? spec.physics.mixing : nothing
    pmix === nothing ||
        @info "inter-trajectory mixing armed" K_h=pmix.K_h sigma0_h=(pmix.sigma0_h > 0 ? pmix.sigma0_h : "auto (met grid)") sigma_v=(pmix.sigma_v > 0 ? pmix.sigma_v : "turbulence yl") background=pmix.background every=pmix.every
    # OSSA sea-spray sources (built in build_run; empty vector = none)
    ssprays = obj.ssprays
    # flow-evolved covariance kernel ([output.covariance]): evolve every
    # parcel's P once per cov_stride sync steps by the resolved deformation +
    # turbulent growth; the deposit/receptors read P at evaluation time.
    do_cov = ctx.output.kernel in (:covariance, :covshape)
    do_cov && @info "covariance kernel armed" mode=spec.output.cov_mode sigma0_m=spec.output.cov_sigma0 sigma_v=spec.output.cov_sigma_v tl_s=spec.output.cov_tl_s stride=spec.output.cov_stride
    framed && @info "grid output: $(NTg) frames every $(spec.output.interval_s) s (avg $(spec.output.average_s) s, sample $(spec.output.sample_s) s)"
    for itime in it0:nsteps
        d = loopstart + itime * ts
        # periodic warm-start checkpoint ([run] checkpoint = duration)
        chk_period_s > 0 && itime > it0 &&
            itime % max(1, round(Int, chk_period_s / dt_s)) == 0 &&
            _checkpoint_save(chk_path, rt, obj, spec, d - ts, seeded)
        # periodic per-parcel diagnostic snapshot ([run] parcel_dump = duration)
        snap_period_s > 0 && itime > it0 &&
            itime % max(1, round(Int, snap_period_s / dt_s)) == 0 &&
            _parcel_snapshot_save(joinpath(snap_dir,
                "snap_" * Dates.format(d - ts, "yyyymmdd_HHMM") * ".jls"),
                rt, obj, d - ts)
        EASLEE.timeit!(obj.bench, :metread) do
            read_met_step!(rt, d)
        end
        ninj = 0
        if dfill !== nothing && itime > 0
            indt_df = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            ninj = EASLEE.timeit!(obj.bench, :dfill) do
                EASLEE.domainfill_inject!(rt.state, ctx, dfill, d, indt_df, dt_s;
                                          rng=df_rng, seeded=seeded)
            end
        end
        # OSSA sea-spray sources: met-driven surface parcel spawning
        if !isempty(ssprays) && itime > 0
            indt_ss = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            EASLEE.timeit!(obj.bench, :seaspray) do
                for ssb in ssprays
                    EASLEE.Operators.apply_seaspray!(rt.state, ctx, rt.metfields2D[1], ssb, d, indt_ss, dt_s)
                end
            end
        end
        # framed output: deposit only within the averaging window, at the sample cadence
        r_out  = framed ? mod(itime, out_ng) : -1
        in_win = framed && (r_out == 0 || r_out > out_ng - avg_ng)
        do_out = framed ? (itime > 0 && in_win && itime % samp_ng == 0) : true
        # [physics].prognostic_airm: refresh the represented air mass from the
        # current parcel spacing every `prognostic_airm_every` transport steps
        # (the dispersion timescale is hours — no need to redo the kNN each step)
        do_airm = spec.physics.prognostic_airm && itime > 0 &&
                  itime % max(spec.physics.prognostic_airm_every, 1) == 0
        EASLEE.step_model!(rt.state, ctx, d; do_output=do_out, do_airm=do_airm,
                           airm_boxes=airm_boxes, kw...)
        do_out && (nsamp_win += 1)
        # engine-driven aerosol BB emission (no chemistry step): inject GFAS OC/BC
        # MASS at the per-cell plume height with the SEVIRI-tuned diurnal + local time
        if emis !== nothing && dfill !== nothing && aero_emit && !do_chem && itime > 0
            EASLEE.timeit!(obj.bench, :emissions) do
                EASLEE.apply_emissions!(nothing, rt.state, emis, ctx.nests[1].gridconfig,
                    dfill.mass_per_parcel, dt_s, kinject;
                    when=d, aerosol=obj.aerosol, lev_h=Float64.(dfill.lev_h))
            end
        end
        if do_cov && itime % spec.output.cov_stride == 0
            EASLEE.timeit!(obj.bench, :covkernel) do
                indt_c = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)
                dt_cov = dt_s * spec.output.cov_stride
                for nn in 1:length(ctx.nests)
                    EASLEE.Operators.update_parcel_cov!(rt.state,
                        ctx.nests[nn].metfields3D, ctx.nests[nn].gridconfig,
                        indt_c[nn], nn, dt_cov;
                        sigma0=spec.output.cov_sigma0, sigma_v=spec.output.cov_sigma_v,
                        tl=spec.output.cov_tl_s, sigma_max=spec.output.cov_sigma_max)
                end
            end
        end
        if recs !== nothing && do_out
            # refresh bandwidths for the gridless KDE (once per output frame)
            # when the grid kernel isn't already doing it
            if recs_need_bw
                fr = out_ng > 0 ? itime ÷ out_ng : itime
                if fr != last_bw_frame[]
                    EASLEE.Operators.update_bandwidths!(rt.state, ctx.nests[1].gridconfig, 1;
                        floor_sv=ctx.output.bwfloor_sv, floor_tl=ctx.output.bwfloor_tl)
                    last_bw_frame[] = fr
                end
            end
            indt_r = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            EASLEE.Operators.receptor_sample!(recs, rt.state, ctx.nests[1].gridconfig,
                ctx.nests[1].levelconfig, ctx.nests[1].metfields3D.zlev, indt_r, 1)
        end
        # umbrella-intrusion spreading (volcanic parcels mode, umbrella_spread=true)
        if !isempty(EASLEE.Operators.UMBRELLA_REG)
            gcn = ctx.nests[1].gridconfig
            glx = ctx.nests[1].metconfig.globalx
            indt_u = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            for (_, (umb, bset)) in EASLEE.Operators.UMBRELLA_REG
                EASLEE.Operators.apply_umbrella_spread!(rt.state, gcn, glx,
                                                        Int(gcn.Ni), umb, bset, d, dt_s;
                                                        met=ctx.nests[1].metfields3D,
                                                        indt=indt_u)
            end
        end
        # inter-trajectory mixing ([physics.mixing] = "pairwise"): Fickian
        # tracer-mass exchange between co-located same-species parcels +
        # coverage-based edge entrainment of ambient air
        if pmix !== nothing && itime > 0 && itime % pmix.every == 0
            indt_m = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            np_mix = EASLEE.timeit!(obj.bench, :mixing) do
                apply_transport_mixing!(rt.state, ctx, oc, pmix, indt_m,
                                        dt_s * pmix.every;
                                        tropo_vec = (obj.chemstate === nothing ? nothing :
                                                     obj.chemstate.tropo))
            end
            verbose && itime % max(1, round(Int, 3600.0 / dt_s)) == 0 &&
                @info "mixing" step=itime pairs=np_mix
        end
        # observed-track ensemble nudging + wind-bias drift ([physics.nudge] /
        # [physics.wind_bias]) — post-advection position corrections in a band
        if itime > 0 && (tnudge !== nothing || wbias !== nothing)
            gcn  = ctx.nests[1].gridconfig
            lcn  = ctx.nests[1].levelconfig
            zln  = ctx.nests[1].metfields3D.zlev
            indt_n = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            if tnudge !== nothing
                nmv, dlo, dla = EASLEE.Operators.apply_track_nudge!(
                    rt.state, tnudge, gcn, lcn, zln, indt_n, d, dt_s)
                nmv > 0 && itime % max(1, round(Int, 3600.0 / dt_s)) == 0 &&
                    @info "track nudge" step=itime parcels=nmv dlon=round(dlo, digits=4) dlat=round(dla, digits=4)
            end
            wbias === nothing ||
                EASLEE.Operators.apply_wind_bias!(rt.state, wbias, gcn, lcn, zln, indt_n, dt_s)
        end
        EASLEE.timeit!(() -> EASLEE.State.sync_dtime_ns!(rt.state), obj.bench, :sync)
        if dfill !== nothing
            indt_df = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            ncull = EASLEE.timeit!(obj.bench, :dfill) do
                EASLEE.domainfill_cull!(rt.state, ctx, dfill, indt_df; seeded=seeded)
            end
            verbose && itime % 4 == 0 &&
                @info "domainfill" step=itime active=count(==(0x01), rt.state.status) injected=ninj culled=ncull
        end
        if framed && itime > 0 && r_out == 0 && nframe_g < NTg
            nframe_g += 1
            recs !== nothing && EASLEE.Operators.receptor_flush!(recs, itime * dt_s)
            aodf = (gw.has_aod && obj.aerosol !== nothing) ?
                   EASLEE.compute_aod(obj.aerosol, rt.state, ctx,
                       EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)) : nothing
            tagf = (gw.has_tags && obj.tagphys !== nothing) ?
                   _grid_tag_fields(obj.tagphys, rt.state, ctx, oc,
                       EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)) : nothing
            _freeze_grid_frame!(gw, oc, nsamp_win, nframe_g, itime * dt_s;
                                aod_field = aodf, tag_field = tagf)
            fill!(oc.concarray, 0.0f0); fill!(oc.massarray, 0.0)
            nsamp_win = 0
        end
        if do_chem
            indt = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1]
            apply_chemistry_step!(rt, obj, spec, itime, d, indt, t_abs0, dt_s, seeded;
                                  emis=emis, dfill=dfill, kinject=kinject)
            # reference co-sampling: obj.chemtl was rotated to `d` inside
            # apply_chemistry_step!, so its bracket is current for this step
            refbg = ref_on ?
                (obj.chemtl.bg0, obj.chemtl.bg1,
                 EASLEE.timeline_weight(obj.chemtl, d), ref_cmap) : nothing
            EASLEE.timeit!(() -> EASLEE.chemconc_sample!(cc, obj.chemstate, rt.state, ctx;
                                                          refbg=refbg), obj.bench, :chemgrid)
            if itime > 0 && itime % out_n == 0 && nframe < length(cc.times)
                nframe += 1
                EASLEE.chemconc_freeze!(cc, nframe, itime * dt_s)
                ref_on && (dg = EASLEE.chem_frame_diag(cc, nframe);
                           isempty(dg) || @info "chem frame $nframe (X vs $(spec.output.reference)):$dg")
                if spec.chemistry.aerosols
                    css = obj.chemstate
                    act = findall(css.active_mask)
                    if !isempty(act)
                        med(v) = round(sort(v[act])[cld(length(act), 2)], sigdigits=3)
                        @info "aerosols frame $nframe (per-parcel medians)" A_wet_m2m3=med(css.awet) AOD550_above=med(css.aod_above) cloud_tau_above=med(css.cld_above)
                    end
                end
            end
        end
        if traj !== nothing
            indt_t = EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)   # per-nest array
            EASLEE.traj_write!(traj, rt.state, ctx, d, indt_t, itime)
        end
        verbose && itime % 10 == 0 &&
            @info "step $itime/$nsteps  $d  active=$(count(==(0x01), rt.state.status))"
    end
    EASLEE.traj_close(traj)
    outname = isempty(spec.output.name) ? spec.name : spec.output.name
    fname = joinpath(spec.output.path, "grid_$(outname).nc")
    if framed && nframe_g > 0
        _grid_close!(gw)            # frames were streamed at each interval
        if recs !== nothing
            rpath = joinpath(dirname(fname), "receptors_" * spec.name * ".csv")
            EASLEE.Operators.receptor_write(recs, rpath, loopstart)
            @info "receptors → $rpath" rows=length(recs.rows)
        end
        @info "grid → $fname  ($(nframe_g) time frames)"
    else
        initnetcdf!(ctx.output, fname)
        writenetcdf(ctx.output, fname, nsteps + 1)
    end
    if do_chem && nframe > 0
        EASLEE.chemconc_close!(cc)
        @info "chemistry grid → $(cc.path)  ($(nframe) frames, $(length(cc.species)) species)"
    end
    isempty(spec.checkpoint) ||
        _checkpoint_save(chk_path, rt, obj, spec, loopstart + nsteps * ts, seeded)
    for ssb in ssprays
        enh = ssb.emitted_mean[] > 0 ? ssb.emitted_gusty[] / ssb.emitted_mean[] : 1.0
        @info "sea spray (OSSA) budget" spawned_parcels=ssb.nspawn[] emitted_kg=round(ssb.emitted[], sigdigits=4) mean_field_kg=round(ssb.emitted_mean[], sigdigits=4) gust_integrated_kg=round(ssb.emitted_gusty[], sigdigits=4) transient_enhancement=round(enh, digits=3)
    end
    if obj.pradphys !== nothing
        pr = obj.pradphys
        act = [i for i in eachindex(rt.state.status)
               if rt.state.status[i] == 0x01 &&
                  1 <= Int(rt.state.relbox[i]) <= length(pr.on) &&
                  pr.on[Int(rt.state.relbox[i])]]
        if !isempty(act)
            w  = Float64.(rt.state.mass[act])
            dT = Float64.(pr.trad[act])
            mw = sum(w .* dT) / max(sum(w), 1e-300)
            @info "plume radiation: accumulated in-plume ΔT_rad" parcels=length(act) mass_weighted_mean_K=round(mw, digits=3) min_K=round(minimum(dT), digits=2) p10_K=round(sort(dT)[max(1, length(dT) ÷ 10)], digits=2)
            let d = EASLEE.Operators.PRAD_DBG
                @info "plume radiation: descent debug" ok=d.ok[] dthdz0=d.dthdz0[] p_low=d.plow[] nonfinite=d.nonfin[] mean_dz_m=(d.ok[] > 0 ? round(d.dzsum[]/d.ok[], digits=2) : NaN)
            end
            # point-of-use verification (trap #1): the layer thickness the
            # 1/h edge-cooling law actually used.  h_ref means the diagnostic
            # never reached the physics — the descent magnitude is then not
            # trustworthy (handoff §2.3).
            if !isempty(pr.hv)
                hva = Float64.(pr.hv[act]); hvp = hva[hva .> 0]
                @info "plume radiation: per-parcel layer thickness h_v (prognostic_airm)" parcels_with_hv=length(hvp) frac=round(length(hvp) / max(length(hva), 1), digits=3) median_m=(isempty(hvp) ? NaN : round(sort(hvp)[cld(length(hvp), 2)], digits=1)) min_m=(isempty(hvp) ? NaN : round(minimum(hvp), digits=1)) max_m=(isempty(hvp) ? NaN : round(maximum(hvp), digits=1)) h_ref_fallback_m=pr.h_ref
            else
                @warn "plume radiation: h_v EMPTY — the 1/h cooling used the constant h_ref; enable [physics].prognostic_airm" h_ref_m=pr.h_ref
            end
            # the STANDING anomaly the saturation gates actually saw —
            # expect ≈ −Q_rad·min(t, τ); compare to the in-plume sonde−ERA5
            # residual (−4/−5 K)
            if pr.tprime_tau > 0 && !isempty(pr.tanom)
                ta  = Float64.(pr.tanom[act])
                mwt = sum(w .* ta) / max(sum(w), 1e-300)
                @info "plume radiation: standing cold anomaly T′ ([physics.plume_radiation].tprime_tau)" tau_h=round(pr.tprime_tau / 3600, digits=1) mass_weighted_mean_K=round(mwt, digits=2) min_K=round(minimum(ta), digits=2) frac_below_1K=round(count(<(-1.0), ta) / length(ta), digits=3)
            end
        end
    end
    if obj.tagphys !== nothing
        tg = obj.tagphys
        act = [i for i in eachindex(rt.state.status) if rt.state.status[i] == 0x01 &&
               tg.n_steps[i] > 0]
        if !isempty(act)
            sf = [Float64(EASLEE.strat_fraction(tg, i)) for i in act]
            nowstrat = count(i -> tg.age_strat[i] == 0.0f0, act)
            @info "transport tags: stratospheric residence ([tracers])" parcels=length(act) mean_strat_frac=round(sum(sf) / length(sf), digits=3) frac_currently_stratospheric=round(nowstrat / length(act), digits=3)
            pf = [Float64(EASLEE.Operators.pbl_fraction(tg, i)) for i in act]
            cf = [Float64(EASLEE.Operators.cont_fraction(tg, i)) for i in act]
            everpbl = count(i -> tg.pbl_steps[i] > 0, act)
            mean_cont = everpbl > 0 ? sum(cf) / everpbl : 0.0
            @info "transport tags: PBL / marine-continental / convection ([tracers])" mean_pbl_frac=round(sum(pf) / length(pf), digits=3) frac_ever_in_pbl=round(everpbl / length(act), digits=3) mean_cont_frac_of_pbl=round(mean_cont, digits=3) mean_conv_count=round(sum(tg.conv_count[i] for i in act) / length(act), digits=2)
        end
    end
    if obj.convphys !== nothing
        ctot = sum(obj.convphys.converted)
        ctot > 0 &&
            @info "species conversion budget ([species.convert])" total_kg=round(ctot, sigdigits=4) per_box=round.(filter(>(0.0), obj.convphys.converted), sigdigits=3)
        # point-of-use check (trap #1) for the sectional path: the SURVIVING
        # crystal-size distribution.  If the radii have not moved off their
        # seeded bin values, the r² law never acted and size sorting is absent
        # — exactly the null that ht_full3 produced with the exponential.
        if !isempty(obj.convphys.rice)
            rv = Float64[]
            wv = Float64[]
            for i in eachindex(rt.state.status)
                rt.state.status[i] == 0x01 || continue
                obj.convphys.rice[i] > 0f0 || continue
                push!(rv, Float64(obj.convphys.rice[i]) * 1e6)   # µm
                push!(wv, Float64(rt.state.mass[i]))
            end
            if !isempty(rv)
                sw = sum(wv)
                rmw = sum(rv .* wv) / max(sw, 1e-300)
                sp_ = sortperm(rv)
                cw = cumsum(wv[sp_]) ./ max(sw, 1e-300)
                q(f) = rv[sp_][max(1, findfirst(>=(f), cw))]
                @info "SECTIONAL ice: surviving crystal radii ([species.convert].sectional)" parcels=length(rv) mass_weighted_mean_um=round(rmw, digits=3) p10_um=round(q(0.10), digits=3) median_um=round(q(0.50), digits=3) p90_um=round(q(0.90), digits=3) min_um=round(minimum(rv), digits=4) max_um=round(maximum(rv), digits=3)
            end
        end
    end
    if obj.icephys !== nothing
        tot = sum(obj.icephys.frozen)
        @info "ice physics: sedimented out of the tagged budget" total_kg=round(tot, sigdigits=4) total_Tg=round(tot / 1e9, sigdigits=4) per_box=round.(obj.icephys.frozen ./ 1e9, sigdigits=3)
        stot = sum(obj.icephys.scavenged)
        stot > 0 &&
            @info "ice physics: co-scavenged by the ice ([species.ice_scav])" total_kg=round(stot, sigdigits=4) total_Tg=round(stot / 1e9, sigdigits=4) per_box=round.(obj.icephys.scavenged ./ 1e9, sigdigits=3)
    end
    obj.bench !== nothing && EASLEE.bench_report(obj.bench)
    return rt, obj, fname
end
