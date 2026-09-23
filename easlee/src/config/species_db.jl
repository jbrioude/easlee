# ─────────────────────────────────────────────────────────────────────────────
# species_db.jl — EASLEE species database.
#
# One record per species with the numbers needed for transport + removal:
#   molar weight, radioactive decay, DRY deposition (gas: reldiff/Henry/f0 or a
#   prescribed vd; aerosol: density/diameter/σ), WET deposition (below-cloud gas
#   A,B; in-cloud aerosol), and the OH reaction k(T)=C·T^N·exp(−D/T) + photolysis,
#   from which a chemical lifetime is derived (needs an [OH] estimate — below).
#
# Field names + units mirror the FLEXPART SPECIES_NNN.namelist so the FLEXPART-11
# set imports verbatim via `parse_flexpart_species` / `import_flexpart_dir!`
# (point it at a FLEXPART options/SPECIES directory).  -9.9/-9e-10 = "unset".
#
# The seeded core below uses values I can source confidently (exact molar weights,
# physical half-lives, standard JPL/IUPAC OH kinetics).  FLEXPART-specific removal
# params (Henry/reldiff/f0, scavenging A,B) are best imported from the SPECIES
# files rather than hand-entered — left `unset` here where not certain.
# ─────────────────────────────────────────────────────────────────────────────
module SpeciesDB

export SpeciesProps, SPECIES_DB, get_species, parse_flexpart_species,
       import_flexpart_dir!, estimate_oh, chemical_lifetime_s, load_species_toml!

import TOML
const UNSET = -9.99e0

Base.@kwdef mutable struct SpeciesProps
    name::String
    molar_mass::Float64        = 29.0      # g/mol
    decay_halflife_s::Float64  = UNSET     # <0 → stable
    # dry deposition — gas
    reldiff::Float64           = UNSET     # diffusivity ratio D(H2O)/D(gas)
    henry::Float64             = UNSET     # Henry's-law constant [M atm⁻¹]
    f0::Float64                = UNSET     # reactivity factor (Wesely) 0–1
    dryvel::Float64            = UNSET     # prescribed dry-dep velocity [m s⁻¹]
    # dry deposition / settling — aerosol
    density::Float64           = UNSET     # [kg m⁻³]  (<0 → gas)
    dquer::Float64             = 0.0       # count-median diameter [m]
    dsigma::Float64            = 0.0       # geometric σ
    # wet deposition
    weta_gas::Float64          = UNSET     # below-cloud gas A
    wetb_gas::Float64          = UNSET     # below-cloud gas B
    crain_aero::Float64        = UNSET     # in-cloud aerosol (rain)
    csnow_aero::Float64        = UNSET     # in-cloud aerosol (snow)
    ccn_aero::Float64          = UNSET     # CCN efficiency
    in_aero::Float64           = UNSET     # IN efficiency
    # gas-phase chemistry:  k_OH(T) = oh_c · T^oh_n · exp(−oh_d/T)  [cm³ molec⁻¹ s⁻¹]
    oh_c::Float64              = UNSET
    oh_d::Float64              = UNSET
    oh_n::Float64              = UNSET
    jrate::Float64             = 0.0       # reference photolysis rate [s⁻¹] (noon, clear)
    note::String               = ""
end

isgas(s::SpeciesProps) = s.density < 0
isset(x) = x > -9.0

# ── FLEXPART SPECIES_NNN.namelist parser ─────────────────────────────────────
function parse_flexpart_species(path::AbstractString)::SpeciesProps
    txt = read(path, String)
    val(key) = (m = match(Regex("P$key\\s*=\\s*([-0-9.Ee+]+)"), txt); m === nothing ? UNSET : parse(Float64, m.captures[1]))
    val2(a, b) = (v = val(a); isset(v) ? v : val(b))   # first key, else alias
    nm = match(r"PSPECIES\s*=\s*\"?\s*([A-Za-z0-9_\-]+)", txt)
    # Radionuclide files leave PDECAY=-9.9 and put the half-life in the comment,
    # e.g. "half life (948886704 sec)".  Recover it when the value is unset.
    decay = val("DECAY")
    if !isset(decay)
        hm = match(r"half\s*life\s*\(([0-9.eE+]+)\s*sec"i, txt)
        hm !== nothing && (decay = parse(Float64, hm.captures[1]))
    end
    SpeciesProps(
        name = nm === nothing ? basename(path) : String(nm.captures[1]),
        molar_mass = val("WEIGHTMOLAR"), decay_halflife_s = decay,
        reldiff = val("RELDIFF"), henry = val("HENRY"), f0 = val("F0"), dryvel = val("DRYVEL"),
        density = val("DENSITY"), dquer = val2("DQUER", "DIA"), dsigma = val("DSIGMA"),
        weta_gas = val("WETA_GAS"), wetb_gas = val("WETB_GAS"),
        crain_aero = val("CRAIN_AERO"), csnow_aero = val("CSNOW_AERO"),
        ccn_aero = val("CCN_AERO"), in_aero = val("IN_AERO"),
        oh_c = val("OHCCONST"), oh_d = val("OHDCONST"), oh_n = val("OHNCONST"),
        note = "imported from $(basename(path))",
    )
end

"Populate the DB from a FLEXPART SPECIES directory (named files + Radionuclides/),
recursively. Skips READMEs / source files."
function import_flexpart_dir!(dir::AbstractString; db = SPECIES_DB)
    n = 0
    for (root, _, files) in walkdir(dir), fn in files
        startswith(fn, "SPECIES_") || continue
        any(endswith(fn, e) for e in (".README", ".f90", ".namelist~")) && continue
        sp = parse_flexpart_species(joinpath(root, fn)); db[uppercase(sp.name)] = sp; n += 1
    end
    @info "imported $n FLEXPART species from $dir"; return n
end

# ── OH estimate for chemical-lifetime diagnostics ────────────────────────────
# A simple diurnally-averaged tropospheric [OH] [molec cm⁻³]: a global reference
# (~1e6) scaled by insolation (cos of solar declination × latitude) and season.
# For a real run, OH comes from a climatology or an external field;
# this is only for the lifetime column.
function estimate_oh(; lat_deg = 0.0, doy = 80, oh_ref = 1.0e6)
    decl = 23.44 * cosd(360.0 * (doy - 172) / 365.0)         # solar declination
    s    = max(cosd(lat_deg - decl), 0.05)                   # insolation proxy
    return oh_ref * (0.4 + 1.2 * s)                          # ~0.5–1.6 ×1e6
end

# ── Chemical lifetime [s] from OH reaction + photolysis ──────────────────────
function chemical_lifetime_s(s::SpeciesProps; T = 273.15, oh = 1.0e6)
    k_oh = isset(s.oh_c) ? s.oh_c * T^(isset(s.oh_n) ? s.oh_n : 0.0) *
                           exp(-(isset(s.oh_d) ? s.oh_d : 0.0) / T) : 0.0
    loss = k_oh * oh + s.jrate
    return loss > 0 ? 1.0 / loss : Inf
end

get_species(name) = get(SPECIES_DB, uppercase(name), nothing)

# ── The database: loaded from the friendly TOML (easlee/data/species.toml,
#    regenerated by scripts/gen_species_db.jl from FLEXPART + EASLEE chemistry). ─
const SPECIES_DB = Dict{String,SpeciesProps}()
_default_db_path() = joinpath(@__DIR__, "..", "..", "data", "species.toml")

"Load the species DB from a friendly TOML (dry/wet/oh are inline tables or \"?\")."
function load_species_toml!(path::AbstractString = _default_db_path(); db = SPECIES_DB)
    isfile(path) || return 0
    num(x)    = x isa Number ? Float64(x) : UNSET          # "?" → UNSET
    gd(d, k)  = (d isa AbstractDict && haskey(d, k) && d[k] isa Number) ? Float64(d[k]) : UNSET
    n = 0
    for e in get(TOML.parsefile(path), "species", Any[])
        nm = String(e["name"]); dry = get(e, "dry", nothing); wet = get(e, "wet", nothing); oh = get(e, "oh", nothing)
        dec = String(get(e, "decay", "stable"))
        decay = dec == "stable" ? UNSET : parse(Float64, replace(dec, "s" => ""))
        db[uppercase(nm)] = SpeciesProps(name = nm,
            molar_mass = num(get(e, "molar_mass", "?")), decay_halflife_s = decay,
            reldiff = gd(dry, "reldiff"), henry = gd(dry, "henry"), f0 = gd(dry, "f0"),
            density = gd(dry, "density"), dquer = max(gd(dry, "dquer"), 0.0), dsigma = max(gd(dry, "dsigma"), 0.0),
            weta_gas = gd(wet, "below_a"), wetb_gas = gd(wet, "below_b"),
            crain_aero = gd(wet, "crain"), csnow_aero = gd(wet, "csnow"),
            oh_c = gd(oh, "c"), oh_d = gd(oh, "d"), oh_n = gd(oh, "n"))
        n += 1
    end
    return n
end

function __init__()
    try
        load_species_toml!()
    catch e
        @warn "species DB load failed (run scripts/gen_species_db.jl)" exception = e
    end
end

end # module SpeciesDB
