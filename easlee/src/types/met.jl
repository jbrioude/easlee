
# ── Typed met-field containers ────────────────────────────────────────────────
#
# MetFields3D / MetFields2D replace the legacy Dict{Any,Any} that was passed
# to every physics kernel.  Named struct fields give the compiler concrete types
# throughout the hot path (no boxing, no dict lookup overhead).
#
# A = concrete 4-D array type (CPU: Array{Float32,4}, GPU: CuArray{Float32,4})
# B = corresponding 3-D slice type for 2-D surface fields
#
# Construction from the legacy Dict is done once in runtime.jl.  The
# struct fields are references to the *same* underlying arrays as the Dict, so
# readmetfiles!(…, dict, …) that writes into those arrays in-place is
# automatically visible through the struct — no synchronisation needed.

struct MetFields3D{A<:AbstractArray{Float32,4}}
    Ux::A        # zonal wind (m/s)
    Vy::A        # meridional wind (m/s)
    levdot::A    # vertical velocity in level-index space (1/s)
    rho::A       # air density (kg/m³)
    theta::A     # potential temperature (K)
    tke::A       # turbulent kinetic energy (m²/s²)
    zlev::A      # geopotential height of model levels (m)
    drhodz::A    # ∂ρ/∂z (kg/m⁴)
    dsigwdz::A   # ∂(σ_w)/∂z
    tls::A       # turbulent length scale (m)
    puv::A       # pressure on UV levels (Pa)  — zeroed for AROME
    T::A         # temperature (K)             — zeroed for AROME
end

struct MetFields2D{B<:AbstractArray{Float32,3}}
    pblh::B      # planetary boundary-layer height (m)
    ustar::B     # friction velocity (m/s)
    wstar::B     # convective velocity scale (m/s)
    ol::B        # Obukhov length (m)
end

# ── Convective / cloud fields (ERA5 IFS_EA_CLD companion files) ───────────────
# Grid-mean time-mean convective diagnostics on the same (Ni, Nj, Nk, 2) layout
# as MetFields3D.  Sign conventions as archived by the IFS:
#   umf  ≥ 0  updraught mass flux        [kg m⁻² s⁻¹]
#   dmf  ≤ 0  downdraught mass flux      [kg m⁻² s⁻¹]  (negative = downward)
#   udr  ≥ 0  updraught detrainment rate [kg m⁻³ s⁻¹]
#   ddr  ≥ 0  downdraught detrainment rate [kg m⁻³ s⁻¹]
#   tpf  ≥ 0  total precipitation flux   [kg m⁻² s⁻¹]  (3-D: flux through level)
#   cc   0–1  cloud fraction
struct ConvFields{A<:AbstractArray{Float32,4}}
    umf::A
    dmf::A
    udr::A
    ddr::A
    tpf::A
    cc::A
end

# ── Diagnostic-convection fields (CAPE-gated scheme) ─────────────────────────
# For datasets without archived convective mass fluxes.  cape/cin are 2-D
# surface diagnostics (Ni, Nj, 2); cin may be a zero array when not archived.
# qv (3-D specific humidity, when the reader archives it) feeds the source
# parcel's θe: reconstruction from Tv = p/(R_d·ρ) FAILS on datasets whose
# reader stores DRY density (ERA5: rho = p/(R_d·T) → q ≡ 0 → no EL is ever
# found despite high CAPE).
struct DiagConvFields{B<:AbstractArray{Float32,3},Q}
    cape::B
    cin::B
    qv::Q          # AbstractArray{Float32,4} or Nothing
end

function DiagConvFields(d2::AbstractDict, d3::AbstractDict = Dict{String,Any}())
    haskey(d2, "cape") || return nothing
    B = typeof(d2["cape"])
    cin = haskey(d2, "cin") ? d2["cin"] : zeros(Float32, size(d2["cape"])...)
    qv  = haskey(d3, "qv") ? d3["qv"] : nothing
    DiagConvFields{B,typeof(qv)}(d2["cape"], cin, qv)
end

# Build from the legacy metfields3D Dict (shares array references, like
# MetFields3D).  Returns `nothing` when the fields the caller needs were not
# listed in the run's variables file.
#   require = :convection (default) — needs the updraft mass flux umf;
#   require = :wetdep — needs only the precipitation flux tpf (cc optional).
# Wet deposition reads tpf and cc and nothing else, so it must not depend on
# umf: until 2026-10-13 it did (one constructor for both), which switched wet
# deposition off silently on every dataset without a convective mass flux —
# GFS among them, whose reader builds tpf from the surface precipitation rate
# (reported by Yue Huang).
function ConvFields(d::AbstractDict; require::Symbol = :convection)
    require in (:convection, :wetdep) || throw(ArgumentError("require = :convection | :wetdep"))
    key = require === :convection ? "umf" : "tpf"
    haskey(d, key) || return nothing
    A = typeof(d[key])
    zero_arr() = zeros(Float32, size(d[key])...)
    # wet deposition reads tpf and cc only (and takes the level count from cc):
    # the absent convective fields share one minimal array instead of four
    # full-size ones (~1.4 GB on a global 0.25° grid)
    tiny = require === :wetdep ? zeros(Float32, ntuple(_ -> 1, ndims(d[key]))...) : nothing
    conv_or_zero(k) = haskey(d, k) ? d[k] : (tiny === nothing ? zero_arr() : tiny)
    ConvFields{A}(
        conv_or_zero("umf"), conv_or_zero("dmf"), conv_or_zero("udr"), conv_or_zero("ddr"),
        haskey(d, "tpf") ? d["tpf"] : zero_arr(),
        haskey(d, "cc")  ? d["cc"]  : zero_arr(),
    )
end

# ── Convenience constructors from legacy Dict ─────────────────────────────────
# Shares array references — O(1), no allocation.

function MetFields3D(d::AbstractDict)
    A = typeof(d["Ux"])
    # Required fields: Ux Vy levdot rho theta tke zlev
    # Optional fields: drhodz dsigwdz tls puv T
    # (absent for WRF and other models that don't produce them directly)
    ref = d["Ux"]
    zero_arr() = zeros(Float32, size(ref)...)
    MetFields3D{A}(
        d["Ux"], d["Vy"], d["levdot"],
        d["rho"], d["theta"], d["tke"], d["zlev"],
        haskey(d, "drhodz")  ? d["drhodz"]  : zero_arr(),
        haskey(d, "dsigwdz") ? d["dsigwdz"] : zero_arr(),
        haskey(d, "tls")     ? d["tls"]     : zero_arr(),
        haskey(d, "puv")     ? d["puv"]     : zero_arr(),
        haskey(d, "T")       ? d["T"]       : zero_arr(),
    )
end

function MetFields2D(d::AbstractDict)
    # pblh/ustar/wstar/ol are not yet populated by all model readers; fall back
    # to a zero array of the right size when the key is absent.
    # All existing 2D arrays share (Ni, Nj, 2); grab any to infer dimensions.
    local ref_arr
    if isempty(d)
        # No 2D fields at all — return a 1×1×2 zero placeholder.
        ref_arr = zeros(Float32, 1, 1, 2)
    else
        ref_arr = first(values(d))
    end
    zero_f() = zeros(Float32, size(ref_arr)...)
    B = Array{Float32,3}
    MetFields2D{B}(
        haskey(d, "pblh")  ? d["pblh"]  : zero_f(),
        haskey(d, "ustar") ? d["ustar"] : zero_f(),
        haskey(d, "wstar") ? d["wstar"] : zero_f(),
        haskey(d, "ol")    ? d["ol"]    : zero_f(),
    )
end

# ── MetConfig — typed replacement for the legacy metconfig Dict ───────────────
#
# Constructed once from the Dict produced by gridcheck_onefile! in initialize.jl.
# The legacy Dict is still passed to readmetfiles! (unchanged) so no reader code
# needs to change.  NestContext.metconfig holds this struct instead of the Dict,
# eliminating the last Any-typed field from NestContext.
#
# filedate is shared by reference with the Dict entry; it is populated once at
# startup and never mutated again, so sharing is safe.

struct MetConfig
    model::String            # "AROME", "ERA5", …
    metfileformat::String    # "NETCDF", "GRIB", …
    Ni::Int32
    Nj::Int32
    Nk::Int32
    dx::Float32              # grid spacing in x (metres for Lambert, degrees for rectilinear)
    dy::Float32              # grid spacing in y
    globalx::Bool            # true for longitude-wrapping domains
    filedate::Vector{DateTime}       # available met file timestamps (shared ref with Dict)
    filename::Vector{String}         # available met file paths    (shared ref with Dict)
    # Coordinate arrays: 2-D (Ni×Nj) for Lambert / curvilinear grids,
    # 1-D (Ni,) or (Nj,) for regular lat-lon / GRIB grids.
    # AbstractArray covers both cases without boxing.
    longitudes::AbstractArray{Float32}
    latitudes::AbstractArray{Float32}
end

function MetConfig(d::AbstractDict)
    ni = Int32(get(d, "Ni", 0))
    nj = Int32(get(d, "Nj", 0))
    # Convert to Float32 regardless of the source element type (GRIB gives Float64).
    lons = Float32.(get(d, "longitudes", zeros(Float32, ni)))
    lats = Float32.(get(d, "latitudes",  zeros(Float32, nj)))
    MetConfig(
        get(d, "model",         ""),
        get(d, "metfileformat", ""),
        ni, nj,
        Int32(get(d, "Nk", 0)),
        Float32(get(d, "dx", 0.0)),
        Float32(get(d, "dy", 0.0)),
        get(d, "globalx", false),
        get(d, "filedate",  DateTime[]),
        get(d, "filename",  String[]),
        lons, lats,
    )
end

# ── MetSnapshot (kept for future GPU snapshot path) ───────────────────────────
#
# Holds two bracketing time slices per 3D/2D field with the time-interpolation
# weight pre-computed.  Not yet wired into the main CPU path.
#
# F = float precision (Float32 or Float64)
# A3 = 4D array type for 3D fields (Ni×Nj×Nk×2)
# A2 = 3D array type for 2D fields (Ni×Nj×2)
struct MetSnapshot{F<:AbstractFloat,
                   A3<:AbstractArray{F,4},
                   A2<:AbstractArray{F,3}}
    # Velocity
    u::A3
    v::A3
    w::A3
    # Thermodynamics
    rho::A3
    theta::A3
    tke::A3
    # Geometry
    zlev::A3    # geopotential height of model levels
    # 2D surface fields
    pblh::A2    # planetary boundary layer height
    ustar::A2   # friction velocity
    wstar::A2   # convective velocity scale
    ol::A2      # Obukhov length
    # Time interpolation weight
    indt::F
end
