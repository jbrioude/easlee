
# ── src/io/output.jl ──────────────────────────────────────────────────────────
# NetCDF output: file initialisation and grid writing.
#
# The two accumulators filled by conccalc! give both output kinds directly:
#   massarray = Σ m_p              → concentration  Σm/(V·nsamples)   [kg m⁻³]
#   concarray = Σ m_p/ρ_p          → mass mixing r. Σ(m/ρ)/(V·nsamp.) [kg kg⁻¹]
# In BACKWARD runs each sample is weighted Δt·(m/M_released) in conccalc!
# (Seibert & Frank 2004), so the same arrays yield source-receptor
# sensitivities:
#   massarray → sens    = ΣΔt·(m/M)/V        [s m⁻³]
#   concarray → sens_mr = ΣΔt·(m/(ρM))/V     [s kg⁻¹]
# (time-integrated — no nsamples normalisation).
#
# Variables written (selected by COMMAND IOUT — 1/9 conc, 2/10 mr, 3/11 both):
#   spec{NNN}_conc / spec{NNN}_mr        per species, aggregated over its boxes
#   spec{NNN}_sens / spec{NNN}_sens_mr   backward equivalents
#   conc_box / mr_box / sens_box / …     per release box (numpoint dim, Nc > 1)
#   release_species(numpoint)            box → species map
#   conc                                 LEGACY variable (unchanged semantics:
#                                        concarray·1e12/nsamples/V, "ng/m3"
#                                        label) kept for the validate_* scripts.
#
# Per-species variables need IOUTPUTFOREACHRELEASE=1 when NSPEC > 1 (with a
# single output class the species cannot be separated after accumulation —
# the aggregate is then written as "specall_*").

using NetCDF
using Printf

_spec_suffix(s::Int) = s == 0 ? "specall" : @sprintf("spec%03d", s)

# (species, columns) groups for the per-species aggregation.  Columns are
# release boxes under "release" grouping, species under "species" grouping
# ([output].group) — oc.classspecies is the per-COLUMN species number either
# way.  Falls back to boxspecies for legacy dicts whose classspecies length
# does not match the accumulator (identity mapping there anyway).
function _spec_groups(oc)
    Nc = Int(oc.Nc)
    cs = length(oc.classspecies) == Nc ? oc.classspecies : oc.boxspecies
    if Nc == 1
        sp = unique(cs)
        return length(sp) == 1 && Int(sp[1]) != 0 ? [(Int(sp[1]), [1])] : [(0, [1])]
    end
    return [(Int(s), findall(==(Int32(s)), cs))
            for s in sort(unique(Int.(cs)))]
end

# per-column species vector actually written/grouped (legacy-dict safe)
_class_species(oc) =
    length(oc.classspecies) == Int(oc.Nc) ? oc.classspecies : oc.boxspecies

# ── NetCDF file initialisation ────────────────────────────────────────────────

# Tolerant field read: `oc` is an OutputConfig here, but legacy checkpoints and
# a few tests still pass the raw Dict, and older serialised OutputConfigs
# predate the CRS fields.
_ocget(oc, k::Symbol, dflt) = oc isa AbstractDict ? get(oc, String(k), dflt) :
                              (hasproperty(oc, k) ? getproperty(oc, k) : dflt)

function initnetcdf!(oc, filename)
    Niout  = Int(oc.Niout)
    Njout  = Int(oc.Njout)
    relbox = Int(oc.Nc)
    backward = oc.direction == -1

    latitude  = [oc.lat0 + oc.dlat*(j-1) for j in 1:Njout]
    longitude = [oc.lon0 + oc.dlon*(i-1) for i in 1:Niout]

    # ── Horizontal axes ──────────────────────────────────────────────────────
    # On a native (Lambert) output CRS lon0/dlon are PROJECTION METRES, so the
    # historic hard-coded "degrees_east" put the grid at 0N 0E with a 1300-degree
    # spacing for every reader that trusted it.  Label the axes for what they
    # actually are and ship the true geographic coordinates as CF 2-D auxiliary
    # coordinate variables.  A lonlat output grid is unchanged (bit-for-bit).
    proj_m = _ocget(oc, :proj_m, false)
    # CF: on a projected grid the data variables point at the 2-D aux coords.
    _cf(a::Dict) = proj_m ? merge(a, Dict("coordinates" => "longitude latitude")) : a
    lat_dim = proj_m ?
        NcDim("y", Njout, atts=Dict("units"=>"m", "standard_name"=>"projection_y_coordinate",
                                    "long_name"=>"y coordinate of projection"); values=latitude) :
        NcDim("latitude",  Njout, atts=Dict("units"=>"degrees_north"); values=latitude)
    lon_dim = proj_m ?
        NcDim("x", Niout, atts=Dict("units"=>"m", "standard_name"=>"projection_x_coordinate",
                                    "long_name"=>"x coordinate of projection"); values=longitude) :
        NcDim("longitude", Niout, atts=Dict("units"=>"degrees_east");  values=longitude)
    alt_dim = NcDim("height",    length(oc.outheight), atts=Dict("units"=>"meters"); values=oc.outheight)
    box_dim = NcDim("numpoint",  relbox, atts=Dict("units"=>"#"); values=collect(1:relbox))
    dims3   = [lon_dim, lat_dim, alt_dim]
    dims4   = [lon_dim, lat_dim, alt_dim, box_dim]

    vars = NcVar[]

    # True geographic coordinates for a projected output grid (CF: the data
    # variables carry `coordinates = "longitude latitude"`).
    if proj_m && !isempty(_ocget(oc, :lon2d, Matrix{Float32}(undef,0,0)))
        push!(vars, NcVar("longitude", [lon_dim, lat_dim];
                          atts=Dict("units"=>"degrees_east", "standard_name"=>"longitude",
                                    "long_name"=>"longitude of grid cell centre"),
                          t=Float32, compress=4))
        push!(vars, NcVar("latitude", [lon_dim, lat_dim];
                          atts=Dict("units"=>"degrees_north", "standard_name"=>"latitude",
                                    "long_name"=>"latitude of grid cell centre"),
                          t=Float32, compress=4))
    end

    # Legacy variable — semantics unchanged (validate_* scripts read it).
    push!(vars, NcVar("conc", relbox == 1 ? dims3 : dims4;
                      atts=_cf(Dict("units"=>"ng/m3")), t=Float32, compress=4))

    do_conc = oc.iout in (1, 3)
    do_mr   = oc.iout in (2, 3)

    for (s, _) in _spec_groups(oc)
        sfx = _spec_suffix(s)
        if backward
            push!(vars, NcVar("$(sfx)_sens", dims3;
                atts=_cf(Dict("units"=>"s m-3",
                          "long_name"=>"source-receptor sensitivity (concentration receptor)",
                          "formula"=>"sum_t dt*(m_parcel/M_released)/V_cell")),
                t=Float32, compress=4))
            push!(vars, NcVar("$(sfx)_sens_mr", dims3;
                atts=_cf(Dict("units"=>"s kg-1",
                          "long_name"=>"source-receptor sensitivity (mixing-ratio receptor)",
                          "formula"=>"sum_t dt*(m_parcel/(rho*M_released))/V_cell")),
                t=Float32, compress=4))
        else
            do_conc && push!(vars, NcVar("$(sfx)_conc", dims3;
                atts=_cf(Dict("units"=>"kg m-3", "long_name"=>"mass concentration")),
                t=Float32, compress=4))
            do_mr && push!(vars, NcVar("$(sfx)_mr", dims3;
                atts=_cf(Dict("units"=>"kg kg-1", "long_name"=>"mass mixing ratio")),
                t=Float32, compress=4))
        end
    end

    # per-column variables: skipped under "species" grouping, where each column
    # IS a species and conc_box/mr_box would duplicate the spec%03d_* variables
    perbox = relbox > 1 && oc.outgroup !== :species
    if perbox
        if backward
            push!(vars, NcVar("sens_box",    dims4; atts=_cf(Dict("units"=>"s m-3")),  t=Float32, compress=4))
            push!(vars, NcVar("sens_mr_box", dims4; atts=_cf(Dict("units"=>"s kg-1")), t=Float32, compress=4))
        else
            do_conc && push!(vars, NcVar("conc_box", dims4;
                atts=_cf(Dict("units"=>"kg m-3")), t=Float32, compress=4))
            do_mr && push!(vars, NcVar("mr_box", dims4;
                atts=_cf(Dict("units"=>"kg kg-1")), t=Float32, compress=4))
        end
    end
    relbox > 1 && push!(vars, NcVar("release_species", [box_dim];
        atts=Dict("long_name"=>"FLEXPART species number of each output column ([output].group: release box or species)"), t=Int32))

    isfile(filename) && rm(filename)
    nc = NetCDF.create(filename, vars...)
    oc.ncid = nc
    relbox > 1 && ncwrite(Vector{Int32}(_class_species(oc)), filename, "release_species")
    if proj_m && !isempty(_ocget(oc, :lon2d, Matrix{Float32}(undef,0,0)))
        ncwrite(Matrix{Float32}(_ocget(oc, :lon2d, Matrix{Float32}(undef,0,0))), filename, "longitude")
        ncwrite(Matrix{Float32}(_ocget(oc, :lat2d, Matrix{Float32}(undef,0,0))), filename, "latitude")
    end
    NetCDF.sync(nc)
end

# ── Grid writer ───────────────────────────────────────────────────────────────
# Normalise the accumulators and write every variable declared above.

function writenetcdf(oc, filename, nbsamples)
    nc     = oc.ncid
    relbox = Int(oc.Nc)
    vol    = oc.volume
    backward = oc.direction == -1

    @assert nbsamples > 0 "nbsamples must be > 0"

    good = isfinite.(vol) .& (vol .> 0)
    # Backward grids are time-integrals: no nbsamples averaging.
    fnorm = backward ? 1.0f0 : 1.0f0 / Float32(nbsamples)

    # cell normalisation (i,j,k) → field/V (zero where the volume is degenerate)
    function _norm(src::AbstractArray{<:Real,3})
        out = zeros(Float32, size(src))
        @inbounds for k in axes(src,3), j in axes(src,2), i in axes(src,1)
            good[i,j,k] && (out[i,j,k] = Float32(src[i,j,k]) * fnorm / vol[i,j,k])
        end
        return out
    end

    # Legacy "conc" (unchanged semantics: concarray·1e12/nsamples/V).
    legacy_norm(src) = begin
        out = zeros(Float32, size(src))
        @inbounds for idx in CartesianIndices(src)
            i, j, k = Tuple(idx)[1], Tuple(idx)[2], Tuple(idx)[3]
            good[i,j,k] && (out[idx] = Float32(src[idx]) * (1.0f12 / nbsamples) / vol[i,j,k])
        end
        out
    end
    if relbox == 1
        ncwrite(legacy_norm(oc.concarray[:, :, :, 1, 1]), filename, "conc")
    else
        ncwrite(legacy_norm(oc.concarray[:, :, :, :, 1]), filename, "conc")
    end

    do_conc = oc.iout in (1, 3)
    do_mr   = oc.iout in (2, 3)

    for (s, boxes) in _spec_groups(oc)
        sfx = _spec_suffix(s)
        magg = dropdims(sum(view(oc.massarray, :, :, :, boxes, 1); dims=4); dims=(4,))
        cagg = dropdims(sum(view(oc.concarray, :, :, :, boxes, 1); dims=4); dims=(4,))
        if backward
            ncwrite(_norm(magg), filename, "$(sfx)_sens")
            ncwrite(_norm(cagg), filename, "$(sfx)_sens_mr")
        else
            do_conc && ncwrite(_norm(magg), filename, "$(sfx)_conc")
            do_mr   && ncwrite(_norm(cagg), filename, "$(sfx)_mr")
        end
    end

    if relbox > 1 && oc.outgroup !== :species
        mbox = zeros(Float32, Int(oc.Niout), Int(oc.Njout), Int(oc.Nkout), relbox)
        cbox = zeros(Float32, size(mbox))
        for b in 1:relbox
            mbox[:, :, :, b] = _norm(view(oc.massarray, :, :, :, b, 1))
            cbox[:, :, :, b] = _norm(view(oc.concarray, :, :, :, b, 1))
        end
        if backward
            ncwrite(mbox, filename, "sens_box")
            ncwrite(cbox, filename, "sens_mr_box")
        else
            do_conc && ncwrite(mbox, filename, "conc_box")
            do_mr   && ncwrite(cbox, filename, "mr_box")
        end
    end

    NetCDF.sync(nc)
end
