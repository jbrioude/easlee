# ──────────────────────────────────────────────────────────────────────────────
# core/config.jl — typed runtime contexts (Config module)
#
# Bridges the legacy Dict-based initialisation (io/initialize.jl) to the typed
# runtime: converts Main-scope grid/level objects to src/types structs, builds
# the ReleaseSchedule, and assembles NestContext / ModelContext.
# Entry point: build_context (exported from EASLEE).
# ──────────────────────────────────────────────────────────────────────────────

module Config

using Dates
using ..State
using ..Types: Lambert, Rectilinear, Mercator, AbstractGrid
using ..Types: SigmaPressureLevel, AbstractLevel
using ..Types: MetFields3D, MetFields2D, MetConfig
using ..Types: RunConfig, NestRunConfig, OutputConfig
using ..Types: NoTurb

# ── Legacy-type converters ────────────────────────────────────────────────────
# The legacy initialize.jl produces Main.Lambert / Main.Rectilinear /
# Main.sigmapressurelevel objects.  These functions convert them to the typed
# structs from src/types/grid.jl and src/types/level.jl so the rest of the
# code never sees untyped legacy objects.

function to_typed_grid(g, globalx::Bool)
    T = typeof(g)
    name = string(nameof(T))

    if name == "Lambert"
        return Lambert(
            g.Ni, g.Nj, g.Nk,
            g.dx, g.dy,
            g.proj_stdlon, g.proj_truelat1, g.proj_truelat2,
            g.proj_hemi, g.proj_cone, g.proj_polei, g.proj_polej, g.proj_rebydx,
            g.coeffx, g.offsetx, g.coeffy, g.offsety,
            g.Uedge, g.Vedge, g.Wedge,
            globalx,
        )
    elseif name == "Rectilinear"
        return Rectilinear(
            g.Ni, g.Nj, g.Nk,
            g.dx, g.dy,
            g.lonmin, g.latmin,
            g.dxconst, g.dyconst,
            g.coeffx, g.offsetx, g.coeffy, g.offsety,
            g.Uedge, g.Vedge, g.Wedge,
            g.dlon, g.dlat,
            globalx,
        )
    elseif name == "Mercator"
        return Mercator(
            g.Ni, g.Nj, g.Nk,
            g.dx, g.dy,
            g.coslat0, g.lonor, g.lator, g.y_or,
            g.lon0out, g.lat0out, g.dxout, g.dyout,
            g.Uedge, g.Vedge, g.Wedge,
            globalx,
        )
    else
        error("Unknown gridconfig type: $T")
    end
end

function to_typed_level(l)
    T = typeof(l)
    name = string(nameof(T))

    if name == "sigmapressurelevel"
        Nk = Int32(l.Nk)
        # WRF has no ak/bk (eta coords): initlevelconfig stores scalar dummies (1.0, 0.0).
        # ztoindl never reads ak/bk — it uses the zlev met array — so fill is safe.
        ak = l.ak isa AbstractVector ? Vector{Float32}(l.ak) : fill(Float32(l.ak), Nk)
        bk = l.bk isa AbstractVector ? Vector{Float32}(l.bk) : fill(Float32(l.bk), Nk)
        return SigmaPressureLevel(Nk, ak, bk)
    else
        error("Unknown levelconfig type: $T")
    end
end

# NestContext is fully parametric — all fields are concretely typed.
# No Any-dispatch barriers remain; the compiler can specialise every method
# that receives a NestContext on the exact grid/level/fields types.
struct NestContext{G<:AbstractGrid, L<:AbstractLevel,
                  M3<:MetFields3D, M2<:MetFields2D}
    id::Int
    metconfig::MetConfig      # typed met-grid metadata
    gridconfig::G
    levelconfig::L
    metfields3D::M3
    metfields2D::M2
    nestrunconfig::NestRunConfig   # per-nest run state (t_index, turbulence, ctl, direction)
end

struct ModelContext
    nests::Vector{NestContext}
    releases::State.ReleaseSchedule
    output::OutputConfig     # typed output-grid config + accumulators
    runconfig::RunConfig     # typed global run parameters + mutable step state
end

@inline function _getreleaseblock(runconfig, input)
    if haskey(runconfig, "RELEASES")
        return runconfig["RELEASES"]
    elseif haskey(input, "RELEASE")
        return input["RELEASE"]
    else
        error("No release block found. Expected runconfig[\"RELEASES\"] or input[\"RELEASE\"].")
    end
end

function build_release_schedule(runconfig, input;
                                F::Type{<:AbstractFloat}=Float32,
                                I::Type{<:Integer}=Int32)
    r = _getreleaseblock(runconfig, input)

    required = ("LON1","LON2","LAT1","LAT2","Z1","Z2",
                "IDATE1","ITIME1","IDATE2","ITIME2",
                "PARTS","MASS","ZKIND")

    for k in required
        haskey(r, k) || error("Release block is missing key \"$k\".")
    end

    nrel = length(r["LON1"])

    for k in required
        length(r[k]) == nrel || error("Release key \"$k\" has length $(length(r[k])) but expected $nrel.")
    end

    timerelease1 = Vector{DateTime}(undef, nrel)
    timerelease2 = Vector{DateTime}(undef, nrel)

    for i in 1:nrel
        t1 = r["IDATE1"][i] * r["ITIME1"][i]
        t2 = r["IDATE2"][i] * r["ITIME2"][i]
        timerelease1[i] = DateTime(t1, dateformat"yyyymmddHHMMSS")
        timerelease2[i] = DateTime(t2, dateformat"yyyymmddHHMMSS")
    end

    # Species per box from RELEASES_CTRL SPECNUM_REL (single value applies to
    # all boxes; a list must match the number of boxes).
    species = fill(I(1), nrel)
    ctrl = haskey(input, "RELEASES_CTRL") ? input["RELEASES_CTRL"] :
           haskey(runconfig, "RELEASES_CTRL") ? runconfig["RELEASES_CTRL"] : nothing
    if ctrl !== nothing && haskey(ctrl, "SPECNUM_REL")
        vals = Int[]
        for entry in ctrl["SPECNUM_REL"], tok in split(String(entry), ",")
            isempty(strip(tok)) || push!(vals, parse(Int, strip(tok)))
        end
        if length(vals) == nrel
            species = I.(vals)
        elseif length(vals) >= 1
            species = fill(I(vals[1]), nrel)
        end
    end

    return State.ReleaseSchedule{F, I}(
        parse.(F, r["LON1"]),
        parse.(F, r["LON2"]),
        parse.(F, r["LAT1"]),
        parse.(F, r["LAT2"]),
        parse.(F, r["Z1"]),
        parse.(F, r["Z2"]),
        timerelease1,
        timerelease2,
        parse.(I, r["PARTS"]),
        parse.(F, r["MASS"]),
        string.(r["ZKIND"]),
        I.(collect(1:nrel)),
        species,
    )
end

"""
    build_context(input, runconfig, outputconfig, metconfig, gridconfig,
                  levelconfig, metfields3D, metfields2D; F=Float32, I=Int32)
        -> ModelContext

Assemble the typed ModelContext from the legacy Dict-based initialisation
products: one NestContext per met nest (typed grid/level/met-field structs that
alias the same arrays the readers fill), the ReleaseSchedule, and the typed
OutputConfig / RunConfig.
"""
function build_context(input, runconfig, outputconfig, metconfig, gridconfig, levelconfig, metfields3D, metfields2D;
                       F::Type{<:AbstractFloat}=Float32, I::Type{<:Integer}=Int32)
    nnest = length(metconfig)

    @assert length(gridconfig)  == nnest
    @assert length(levelconfig) == nnest
    @assert length(metfields3D) == nnest
    @assert length(metfields2D) == nnest

    nests = Vector{NestContext}(undef, nnest)
    for n in 1:nnest
        # MetConfig is built first so globalx can be read from it directly,
        # eliminating the need for the per-nest runconfig[n] Dict.
        typed_mc    = MetConfig(metconfig[n])
        typed_grid  = to_typed_grid(gridconfig[n], typed_mc.globalx)
        typed_level = to_typed_level(levelconfig[n])
        # Convert the legacy Dict arrays to typed MetFields structs.
        # The struct fields are references to the SAME underlying arrays as the
        # Dict, so readmetfiles!(…, dict, …) writes remain visible through the
        # struct automatically — no extra synchronisation step needed.
        typed_m3  = MetFields3D(metfields3D[n])
        typed_m2  = MetFields2D(metfields2D[n])
        typed_nrc = NestRunConfig(
            (0, 0),
            NoTurb(),
            Float32(runconfig["ctl"]),
            Int8(runconfig["direction"]),
            (true, true, true),          # 3D by default; build_run overrides from spec.wind
            false, 0.16f0,               # mesoscale fluctuations off by default
            # OU decorrelation = met archive interval (from the AVAILABLE scan)
            Float32(length(typed_mc.filedate) >= 2 ?
                    Dates.value(typed_mc.filedate[2] - typed_mc.filedate[1]) / 1000 : 3600.0),
            0.0f0,                       # kh_ft off by default (build_run sets from spec)
            false,                       # horiz_tke off by default (build_run sets from spec)
            false,                       # horiz_hanna off by default (build_run sets from spec)
            # this nest's own (native_e, r_w, c_k) from vmix_from_namelist
            let vm = get(metconfig[n], "vmix", (false, NaN32, NaN32))
                (Bool(vm[1]), Float32(vm[2]), Float32(vm[3]))
            end...,
            false, false, Inf32, -Inf32, # polar caps off (build_run enables from spec + grid)
            false,                       # MO surface-layer wind (build_run sets from spec)
        )
        nests[n] = NestContext(
            n,
            typed_mc,
            typed_grid,
            typed_level,
            typed_m3,
            typed_m2,
            typed_nrc,
        )
    end

    releases = build_release_schedule(runconfig, input; F=F, I=I)

    # Output semantics consumed by initnetcdf!/writenetcdf — recorded on the
    # OutputConfig so the writers need no extra context.
    outputconfig["iout"]       = get(runconfig, "iout", 1)
    outputconfig["direction"]  = runconfig["direction"]
    outputconfig["boxspecies"] = collect(Int32,  releases.species)
    # backward norm is per OUTPUT COLUMN: aggregate released mass over the
    # boxes mapped to each column (identity under "release" grouping)
    let mb = collect(Float32, releases.mass_per_box),
        ocl = Vector{Int32}(get(outputconfig, "outclass", Int32[]))
        if isempty(ocl) || isempty(mb)
            outputconfig["relmass"] = mb
        else
            rm = zeros(Float32, Int(maximum(ocl)))
            for (b, c) in enumerate(ocl)
                b <= length(mb) && (rm[Int(c)] += mb[b])
            end
            outputconfig["relmass"] = rm
        end
    end
    typed_output = OutputConfig(outputconfig)

    Nc = get(runconfig, "Nc", 1)
    typed_rc = RunConfig(runconfig, Nc)

    return ModelContext(
        nests,
        releases,
        typed_output,
        typed_rc,
    )
end

end
