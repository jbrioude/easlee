# ──────────────────────────────────────────────────────────────────────────────
# io/tools.jl — grid/projection types and interpolation helpers (Main scope)
#
# Defines the horizontal grid structs (Rectilinear, Lambert, Mercator), the
# vertical-level structs (sigmapressurelevel, ...), the turbulence marker types
# (Thomson, Hanna), scalar multi-linear interpolation kernels (1D/3D/4D), and
# lat/lon ↔ grid-index conversions incl. the WRF Lambert-conformal projection
# math (lc_cone / set_lc / latlontoindyindx — WPS map_utils lineage).
# The canonical batch interpolation routines live in src/met/interp.jl; these
# scalar variants serve Main-scope callers (initialize.jl, runtime.jl).
# ──────────────────────────────────────────────────────────────────────────────

using GRIB
using Dates

const FM = Float32
const FH = Float64
const FL = Float16
const IM = Int32
const IH = Int64
const IL = Int16

struct sigmapressurelevel #in ECMWF and AROME for instance
   Nk
   ak
   bk
end
mutable struct sigmageopotentiallevel #in MESONH for instance
   Nk
   wname
   zname
   Zlevels
end
mutable struct pressurelevel #in ECMWF or NCEP/GFS
   Nk
   wname
   zname
   plevels
end
mutable struct hhllevel #in ICON
   Nk
   wname
   zname
   hhl
end
mutable struct potentialtemperaturelevel #for diabatic motion
   Nk
   wname
   zname
   thetalevels
end

#there are different z level types, and horizontal coordinate types. the z level type
#is defined within the horizontal coordinate type to make a single struct

struct Rectilinear
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
end
struct Lambert
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
   geoout::Bool   # output grid geographic (true) or native metres (false)
end
struct Mercator   # MesoNH conformal projection with RPK = 0
   Ni::Int32
   Nj::Int32
   Nk::Int32
   dx::Float32
   dy::Float32
   coslat0::Float32
   lonor::Float32
   lator::Float32
   y_or::Float32
   lon0out::Float32
   lat0out::Float32
   dxout::Float32
   dyout::Float32
   Uedge::Bool
   Vedge::Bool
   Wedge::Bool
end
struct Thomson
   BL89::Bool
   read_tls::Bool   
end
struct Hanna
   turbswitch::Bool
end
mutable struct Lambertconformal #it is assumed that one can get the exact location based on the projection information provided
   zleveltype
   latcen
   loncen
end

#linear interpolation
#1D
@inline c0(di) = (1.f0 - di)
@inline c1(di) = di
#3D
@inline c000(di, dj, dk) = (1.f0 - di) * (1.f0 - dj) * (1.f0 - dk)
@inline c001(di, dj, dk) = (1.f0 - di) * (1.f0 - dj) *      dk
@inline c010(di, dj, dk) = (1.f0 - di) *      dj  * (1.f0 - dk)
@inline c011(di, dj, dk) = (1.f0 - di) *      dj  *      dk
@inline c100(di, dj, dk) =      di  * (1.f0 - dj) * (1.f0 - dk)
@inline c101(di, dj, dk) =      di  * (1.f0 - dj) *      dk
@inline c110(di, dj, dk) =      di  *      dj  * (1.f0 - dk)
@inline c111(di, dj, dk) =      di  *      dj  *      dk
#4D — time orientation fixed 2026-08-11 (same mirror as met/interp.jl: slot l
# must carry (1 − dl); see test/test_time_interp_orientation.jl)
@inline c0000(di, dj, dk, dl) = (1.f0- di) * (1.f0- dj) * (1.f0- dk) * (1.f0- dl)
@inline c0010(di, dj, dk, dl) = (1.f0- di) * (1.f0- dj) *      dk  * (1.f0- dl)
@inline c0100(di, dj, dk, dl) = (1.f0- di) *      dj  * (1.f0- dk) * (1.f0- dl)
@inline c0110(di, dj, dk, dl) = (1.f0- di) *      dj  *      dk  * (1.f0- dl)
@inline c1000(di, dj, dk, dl) =      di  * (1.f0- dj) * (1.f0- dk) * (1.f0- dl)
@inline c1010(di, dj, dk, dl) =      di  * (1.f0- dj) *      dk  * (1.f0- dl)
@inline c1100(di, dj, dk, dl) =      di  *      dj  * (1.f0- dk) * (1.f0- dl)
@inline c1110(di, dj, dk, dl) =      di  *      dj  *      dk  * (1.f0- dl)
@inline c0001(di, dj, dk, dl) = (1.f0- di) * (1.f0- dj) * (1.f0- dk) *      dl
@inline c0011(di, dj, dk, dl) = (1.f0- di) * (1.f0- dj) *      dk  *      dl
@inline c0101(di, dj, dk, dl) = (1.f0- di) *      dj  * (1.f0- dk) *      dl
@inline c0111(di, dj, dk, dl) = (1.f0- di) *      dj  *      dk  *      dl
@inline c1001(di, dj, dk, dl) =      di  * (1.f0- dj) * (1.f0- dk) *      dl
@inline c1011(di, dj, dk, dl) =      di  * (1.f0- dj) *      dk  *      dl
@inline c1101(di, dj, dk, dl) =      di  *      dj  * (1.f0- dk) *      dl
@inline c1111(di, dj, dk, dl) =      di  *      dj  *      dk  *      dl

#1D
@inline interpolate1D(field, di, i) =
    @inbounds (c0(di)*field[i] + c1(di)*field[i+1])

#3D
@inline interpolate3D(field, di, dj, dk, i, j, k) =
    @inbounds (  c000(di, dj, dk) * field[i,   j,   k  ]
               + c001(di, dj, dk) * field[i,   j,   k+1]
               + c010(di, dj, dk) * field[i,   j+1, k  ]
               + c011(di, dj, dk) * field[i,   j+1, k+1]
               + c100(di, dj, dk) * field[i+1, j,   k  ]
               + c101(di, dj, dk) * field[i+1, j,   k+1]
               + c110(di, dj, dk) * field[i+1, j+1, k  ]
               + c111(di, dj, dk) * field[i+1, j+1, k+1])

@inline @views interpolate3Dirreg(field, di, dj, dk, i, j, k,l) =
    @inbounds (  c000(di, dj, dk) * field[i,   j,   l, k  ]
               + c001(di, dj, dk) * field[i,   j,   l, k+1]
               + c010(di, dj, dk) * field[i,   j+1, l, k  ]
               + c011(di, dj, dk) * field[i,   j+1, l, k+1]
               + c100(di, dj, dk) * field[i+1, j,   l, k  ]
               + c101(di, dj, dk) * field[i+1, j,   l, k+1]
               + c110(di, dj, dk) * field[i+1, j+1, l, k  ]
               + c111(di, dj, dk) * field[i+1, j+1, l, k+1])
#4D
@inline interpolate4D(field, di, dj, dk, dl, i, j, k, l) =
               (  c0000(di, dj, dk, dl) * field[i,   j,   k  , l  ]
               + c0010(di, dj, dk, dl) * field[i,   j,   k+1, l  ]
               + c0100(di, dj, dk, dl) * field[i,   j+1, k  , l  ]
               + c0110(di, dj, dk, dl) * field[i,   j+1, k+1, l  ]
               + c1000(di, dj, dk, dl) * field[i+1, j,   k  , l  ]
               + c1010(di, dj, dk, dl) * field[i+1, j,   k+1, l  ]
               + c1100(di, dj, dk, dl) * field[i+1, j+1, k  , l  ]
               + c1110(di, dj, dk, dl) * field[i+1, j+1, k+1, l  ]
               + c0001(di, dj, dk, dl) * field[i,   j,   k  , l+1]
               + c0011(di, dj, dk, dl) * field[i,   j,   k+1, l+1]
               + c0101(di, dj, dk, dl) * field[i,   j+1, k  , l+1]
               + c0111(di, dj, dk, dl) * field[i,   j+1, k+1, l+1]
               + c1001(di, dj, dk, dl) * field[i+1, j,   k  , l+1]
               + c1011(di, dj, dk, dl) * field[i+1, j,   k+1, l+1]
               + c1101(di, dj, dk, dl) * field[i+1, j+1, k  , l+1]
               + c1111(di, dj, dk, dl) * field[i+1, j+1, k+1, l+1])

#interpolation functions — scalar single-point variants used by Main-scope callers.
# Batch variants (interpolation4Dall, interpolation4D, interpolation4Dlev) and the
# irregular-z variants (interpolation4D_irregz, interpolation4D_irregzsingle) have
# been removed — the canonical versions live in src/met/interp.jl (EASLEE module).
function interpolation1Dsingle(field,i)
   dt = mod(i,1)
   l = Base.unsafe_trunc(Int,i)
   return interpolate1D(field,dt,l)
end
@inline function get_vertical_profile(field::Array{Float32}, indx::Float32, indy::Float32, indt::Float32, nz::Int32)
    # This function extracts a vertical profile at a given x,y location by performing 
    # trilinear interpolation at each vertical level
    
    # Initialize the profile array
    profile = Vector{Float32}(undef, nz)
    
    # Get the integer and fractional parts of the indices
    i1 = Base.unsafe_trunc(Int32, indx)
    i2 = i1 + 1
    j1 = Base.unsafe_trunc(Int32, indy)
    j2 = j1 + 1
    t1 = Base.unsafe_trunc(Int32, indt)
    t2 = min(t1 + 1, size(field, 4))  # Ensure we don't go beyond array bounds
    
    di = mod(indx, 1.0f0)
    dj = mod(indy, 1.0f0)
    dt = mod(indt, 1.0f0)
    
    # For each vertical level, interpolate horizontally and in time
    @inbounds for k in 1:nz
        # Interpolate at time t1 using bilinear interpolation in x,y
        v1 = (1.0f0 - di) * (1.0f0 - dj) * field[i1, j1, k, t1] +
             (1.0f0 - di) * dj * field[i1, j2, k, t1] +
             di * (1.0f0 - dj) * field[i2, j1, k, t1] +
             di * dj * field[i2, j2, k, t1]
        
        # Interpolate at time t2 using bilinear interpolation in x,y
        v2 = (1.0f0 - di) * (1.0f0 - dj) * field[i1, j1, k, t2] +
             (1.0f0 - di) * dj * field[i1, j2, k, t2] +
             di * (1.0f0 - dj) * field[i2, j1, k, t2] +
             di * dj * field[i2, j2, k, t2]
        
        # Final linear interpolation in time
        profile[k] = (1.0f0 - dt) * v1 + dt * v2
    end
    
    return profile
end
# Stacked read of a full-level (class F) field at level index z — the same map
# as src/met/interp.jl (_split_t, _split_z_stagger, _clamp_top_k): index z + ½
# for z ≥ 1.5, 2z − 1 below (indl = 1 → the surface slot 1, 1.5 → slot 2), and
# the top stencil shifted down so the k+1 point stays in bounds.  Until
# 2026-10-22 this Main-scope copy kept the old near-surface map 2z + ½ (slots 2–3
# at the ground) and no top clamp; Main-scope callers (test_freeze_diag,
# test_sectional_ice, where it shadows the imported module function) now get
# the engine's read.
function interpolation4Dsingle(field,x,y,z,t)
   if t >= 2f0
     dl = 1.0f0; l = 1
   else
     dl = Float32(mod(t,1)); l = Base.unsafe_trunc(Int,t)
   end
   zz = z >= 1.5f0 ? Float32(z) + 0.5f0 : max(1.0f0, 2f0 * Float32(z) - 1.0f0)
   dz = Float32(mod(zz,1)); kk = Base.unsafe_trunc(Int,zz)
   if kk >= size(field,3)
     dz = 1.0f0; kk = size(field,3) - 1
   end
   return interpolate4D(field, Float32(mod(x,1)), Float32(mod(y,1)), dz, dl,
                        Base.unsafe_trunc(Int,x), Base.unsafe_trunc(Int,y), kk, l)
end

function indytolat(gridconfig::Rectilinear,indy)
   return indy*gridconfig.dy + gridconfig.latmin -1.f0
end
function dxtodindx(gridconfig::Rectilinear,dx,indx,indy)
   return dx*gridconfig.dxconst / cosd(indytolat(gridconfig,indy))
end
function dytodindy(gridconfig::Rectilinear,dy,indx,indy)
   return dy*gridconfig.dyconst
end
# (Lambert dx→index conversions live in met/coords.jl, with the map factor)
function lontoindx(gridconfig::Rectilinear,lon)
 # wrap into [lonmin, lonmin+360): a −180..180 release/receptor lon must land
 # on a 0–360-convention met crop (CAPTEX: lonmin=260, release −84.05 → 275.95);
 # in-frame lons are unchanged, genuinely-outside points map beyond Ni.
 return  mod.(lon .- gridconfig.lonmin, 360.) ./ gridconfig.dlon .+ 1.
end
function lattoindy(gridconfig::Rectilinear,lat)
 return  (lat .- gridconfig.latmin) ./ gridconfig.dlat .+ 1.
end
function latlontoindyindx(gridconfig::Rectilinear,lat,lon)
 return lattoindy(gridconfig::Rectilinear,lat),lontoindx(gridconfig::Rectilinear,lon)
end

# Constants used in the calculations
const RAD_PER_DEG = π / 180
const DEG_PER_RAD = 180 / π

function lc_cone(truelat1, truelat2)
    if abs(truelat1 - truelat2) > 0.1
        # Secant projection cone factor calculation
        cone = (log(cos(truelat1 * RAD_PER_DEG)) - log(cos(truelat2 * RAD_PER_DEG))) /
               (log(tan((90 - abs(truelat1)) * RAD_PER_DEG * 0.5)) -
                log(tan((90 - abs(truelat2)) * RAD_PER_DEG * 0.5)))
    else
        # Tangent projection cone factor calculation
        cone = sin(abs(truelat1) * RAD_PER_DEG)
    end
    return cone
end

function set_lc(proj_truelat1, proj_truelat2, proj_lon1, proj_stdlon, proj_rebydx, proj_hemi, proj_lat1)
    # Compute cone factor
    proj_cone = lc_cone(proj_truelat1, proj_truelat2)

    # Adjust longitude differences to stay out of the "cut zone"
    deltalon1 = mod(proj_lon1 - proj_stdlon + 180, 360) - 180

    # Convert truelat1 to radians and compute cosine for later use
    tl1r = proj_truelat1 * RAD_PER_DEG
    ctl1r = cos(tl1r)

    # Compute the radius to the lower-left (south-west) corner
    proj_rsw = proj_rebydx * ctl1r / proj_cone *
               (tan((90 * proj_hemi - proj_lat1) * RAD_PER_DEG / 2) /
                tan((90 * proj_hemi - proj_truelat1) * RAD_PER_DEG / 2)) ^ proj_cone

    # Compute the pole point
    arg = proj_cone * (deltalon1 * RAD_PER_DEG)
    proj_polei = 1.0 - proj_hemi * proj_rsw * sin(arg)
    proj_polej = 1.0 + proj_rsw * cos(arg)

    # Return the computed values
    return proj_cone, proj_rsw, proj_polei, proj_polej
end

function latlontoindyindx(gridconfig::Lambert,lat,lon)
    # Normalize the longitude difference
    deltalon = ((lon .- gridconfig.proj_stdlon .+ 180.) .% 360.) .- 180.

    # Convert truelat1 to radians and compute cos for later use
    tl1r = gridconfig.proj_truelat1 * RAD_PER_DEG
    ctl1r = cos(tl1r)

    # Compute the tangent operations in a single line to minimize function calls
    tan_half_diff = tan.((90. * gridconfig.proj_hemi .- lat) .* RAD_PER_DEG ./ 2.)
    tan_half_true_lat = tan((90. * gridconfig.proj_hemi - gridconfig.proj_truelat1) * RAD_PER_DEG / 2.)

    # Calculate the radius to the desired point
    rm = gridconfig.proj_rebydx * ctl1r / gridconfig.proj_cone .* (tan_half_diff ./ tan_half_true_lat).^gridconfig.proj_cone

    # Argument for sinusoidal projection adjustment
    arg = gridconfig.proj_cone .* (deltalon .* RAD_PER_DEG)

    # Calculate Cartesian coordinates
    i = gridconfig.proj_polei .+ gridconfig.proj_hemi .* rm .* sin.(arg)
    j = gridconfig.proj_polej .- rm .* cos.(arg)

    # Adjust coordinates for the southern hemisphere
    if gridconfig.proj_hemi == -1
        i = 2. .- i
        j = 2. .- j
    end
    return j, i
end
function indyindxtolatlon(gridconfig::Lambert,j,i)
    # Check hemisphere and adjust indices if necessary
    inew = gridconfig.proj_hemi == -1 ? 2 - i : i
    jnew = gridconfig.proj_hemi == -1 ? 2 - j : j

    # Compute radius squared to the (i, j) location
    xx = inew - gridconfig.proj_polei
    yy = gridconfig.proj_polej - jnew
    r2 = xx^2 + yy^2
    r = sqrt(r2) / gridconfig.proj_rebydx

    if r2 == 0
        lat = gridconfig.proj_hemi * 90
        lon = gridconfig.proj_stdlon
    else
        # Compute longitude (two-argument atan = Fortran atan2)
        lon = gridconfig.proj_stdlon + DEG_PER_RAD * atan(gridconfig.proj_hemi * xx, yy) / gridconfig.proj_cone
        lon = (lon + 360) % 360

        # Compute latitude
        chi1 = (90 - gridconfig.proj_hemi * gridconfig.proj_truelat1) * RAD_PER_DEG
        chi2 = (90 - gridconfig.proj_hemi * gridconfig.proj_truelat2) * RAD_PER_DEG

        if chi1 == chi2
            chi = 2 * atan((r / tan(chi1))^(1 / gridconfig.proj_cone) * tan(chi1 / 2))
        else
            chi = 2 * atan((r * gridconfig.proj_cone / sin(chi1))^(1 / gridconfig.proj_cone) * tan(chi1 / 2))
        end

        lat = (90 - chi * DEG_PER_RAD) * gridconfig.proj_hemi
    end

    return lat, lon
end

# Met-grid index → output-grid coordinate (scalar; the canonical versions live
# in src/met/coords.jl inside the EASLEE module).
@inline indxtoxout(gridconfig::Rectilinear,indx,indy) = indx*gridconfig.coeffx + gridconfig.offsetx
@inline indytoyout(gridconfig::Rectilinear,indx,indy) = indy*gridconfig.coeffy + gridconfig.offsety
# (Lambert output mapping: met/coords.jl, through lon/lat for a geographic grid)
@inline test_z(z,outheight)::Int = findfirst((x) -> x >= z,outheight)
