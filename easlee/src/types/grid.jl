
# Horizontal grid types.
# All projections share Ni, Nj, Nk and stagger flags; projection-specific
# parameters live in the concrete subtypes.

abstract type AbstractGrid end

struct Rectilinear <: AbstractGrid
    Ni::Int32
    Nj::Int32
    Nk::Int32
    dx::Float32
    dy::Float32
    lonmin::Float32
    latmin::Float32
    dxconst::Float32
    dyconst::Float32
    coeffx::Float32
    offsetx::Float32
    coeffy::Float32
    offsety::Float32
    Uedge::Bool
    Vedge::Bool
    Wedge::Bool
    dlon::Float32
    dlat::Float32
    globalx::Bool
end

struct Lambert <: AbstractGrid
    Ni::Int32
    Nj::Int32
    Nk::Int32
    dx::Float32
    dy::Float32
    proj_stdlon::Float32
    proj_truelat1::Float32
    proj_truelat2::Float32
    proj_hemi::Float32
    proj_cone::Float32
    proj_polei::Float32
    proj_polej::Float32
    proj_rebydx::Float32
    coeffx::Float32
    offsetx::Float32
    coeffy::Float32
    offsety::Float32
    Uedge::Bool
    Vedge::Bool
    Wedge::Bool
    globalx::Bool
end

# Mercator projection (MesoNH conformal with RPK = 0).
# Grid spacing dx/dy is uniform in the MAP plane [m]; the scale factor
# k = cos(lat0)/cos(lat) converts physical to map distances.
struct Mercator <: AbstractGrid
    Ni::Int32
    Nj::Int32
    Nk::Int32
    dx::Float32        # map-plane spacing [m]
    dy::Float32
    coslat0::Float32   # cos of reference latitude LAT0
    lonor::Float32     # lon of grid point (1,1) [deg]
    lator::Float32     # lat of grid point (1,1) [deg]
    y_or::Float32      # Mercator y of grid point (1,1) [m]
    lon0out::Float32   # output-grid origin [deg]
    lat0out::Float32
    dxout::Float32     # output-grid spacing [deg]
    dyout::Float32
    Uedge::Bool
    Vedge::Bool
    Wedge::Bool
    globalx::Bool
end

# Turbulence scheme tags — used for dispatch on turbulence!()
struct Thomson
    BL89::Bool
    read_tls::Bool
end

struct Hanna
    turbswitch::Bool    # true = real Hanna on-the-fly (option b)
    sigw_ft::Float32    # free-troposphere background σ_w (m/s); default 0.1
end
Hanna(turbswitch::Bool) = Hanna(turbswitch, 0.1f0)

# No-turbulence tag — pure advection, validation / diagnostic runs.
struct NoTurb end
