
# ── Horizontal coordinate transforms ─────────────────────────────────────────
# All functions return (indy, indx) to match the FLEXPART convention
# (j=latitude index, i=longitude index).

const _RAD_PER_DEG = π / 180.0
const _DEG_PER_RAD = 180.0 / π

# Rectilinear (regular lat/lon) ───────────────────────────────────────────────

@inline function latlontoindyindx(g::Rectilinear, lat, lon)
    # Normalize longitude into the grid's [lonmin, lonmin+360) range — for ANY
    # rectilinear grid, not just global ones: a 0–360-convention REGIONAL crop
    # (CAPTEX ERA5: 260–300°E) must accept a −180..180 release/receptor lon
    # (−84.05 ≡ 275.95).  In-frame lons are unchanged; genuinely-outside lons
    # map beyond Ni and are rejected downstream exactly as before.
    lon_n = mod(lon - g.lonmin, 360.0f0) + g.lonmin
    indy = (lat   - g.latmin) / g.dlat + 1.0f0
    indx = (lon_n - g.lonmin) / g.dlon + 1.0f0
    return indy, indx
end

@inline function dxtodindx(g::Rectilinear, dx, indx, indy)
    lat = (indy - 1.0f0) * g.dy + g.latmin
    return dx * g.dxconst / cosd(lat)
end

@inline dytodindy(g::Rectilinear, dy, indx, indy) = dy * g.dyconst

# On a zonally-periodic (global) met grid the output grid may cross the grid's
# lon seam (e.g. OUTGRID 150–230°E on a −180..180 met frame): the linear map
# then puts parcels east of the seam one full 360° period outside the output
# window.  Wrap xout by the met-grid period toward [0.5, ∞) so a
# dateline-crossing OUTGRID fills on both sides of the seam.
@inline function indxtoxout(g::Rectilinear, indx, indy)
    xout = indx * g.coeffx + g.offsetx
    if g.globalx
        xper = 360.0f0 * g.coeffx / g.dlon        # one lon wrap in xout units
        xout < 1.0f0 && (xout += xper)              # cells start at xout = 1 (2026-10-13)
        xout - xper >= 1.0f0 && (xout -= xper)
    end
    return xout
end
@inline indytoyout(g::Rectilinear, indx, indy) = indy * g.coeffy + g.offsety

# Inverse: grid index → geographic coordinates (Rectilinear is trivial)
@inline function indyindxtolatlon(g::Rectilinear, indy, indx)
    lat = g.latmin + (indy - 1.0) * g.dlat
    lon = g.lonmin + (indx - 1.0) * g.dlon
    return lat, lon
end

# Polar-stereographic caps (FLEXPART switchnorth/switchsouth analogue) ────────
# Poleward of the configured switch latitude the lat/lon metric update breaks
# down: dxtodindx carries a 1/cos(lat) factor (→ ∞ at the pole) and a parcel
# crossing the pole cannot be represented by a monotone indy step at all.
# Following FLEXPART (advance_mod update_xy, ngrid = -1/-2), the horizontal
# position update is done in a pole-tangent stereographic plane instead:
#
#   ρ = 2R·tan((90° − hemi·lat)/2)     x = ρ·sin(lon)     y = −hemi·ρ·cos(lon)
#
# with hemi = +1 (north cap) / −1 (south cap).  The projection is conformal:
# local (east, north) → map (x, y) is a pure rotation by hemi·lon, and map
# distance = k · ground distance with k = 2/(1 + sin(hemi·lat)) (k = 1 at the
# pole).  The accumulated displacement (dxsave east, dysave north) is sampled
# in the parcel's start-of-step local frame — which is FIXED during the
# sub-step loop — so one rotation of the accumulated total is exactly
# equivalent to FLEXPART's per-grid-point uupol/vvpol rotation.  k is
# evaluated at the start latitude (FLEXPART's cgszll convention).
# Float64 internals: near the pole ρ is the difference of large/small terms
# and Float32 loses the longitude.

const _PS_2R = 2.0 * 6.371e6    # stereo plane: ρ = 2R·tan(colat/2)

@inline function polar_stereo_update(g::Rectilinear, hemi::Float32,
                                     indy, indx, dxsave, dysave)
    lat, lon = indyindxtolatlon(g, indy, indx)
    h    = Float64(hemi)
    slon = sind(lon)
    clon = cosd(lon)
    ρ = _PS_2R * tand((90.0 - h * lat) / 2.0)
    x = ρ * slon
    y = -h * ρ * clon
    k = 2.0 / (1.0 + sind(h * lat))          # map = k × ground distance
    x += k * (Float64(dxsave) * clon - h * Float64(dysave) * slon)
    y += k * (h * Float64(dxsave) * slon + Float64(dysave) * clon)
    ρn   = sqrt(x * x + y * y)
    latn = h * (90.0 - 2.0 * atand(ρn / _PS_2R))
    lonn = ρn > 0.0 ? atand(x, -h * y) : lon   # exactly at the pole: keep lon
    return latlontoindyindx(g, Float32(latn), Float32(lonn))
end

# Non-rectilinear grids never activate the polar caps (build_run only enables
# them on global Rectilinear grids that reach the pole); this fallback exists
# so the polar branch in updatepos! typechecks for every grid type.
@inline polar_stereo_update(g::AbstractGrid, hemi::Float32, indy, indx, dxsave, dysave) =
    (indy + dytodindy(g, dysave, indx, indy),
     indx + dxtodindx(g, dxsave, indx, indy))

# Lambert conformal ───────────────────────────────────────────────────────────

function latlontoindyindx(g::Lambert, lat, lon)
    deltalon = mod(lon - g.proj_stdlon + 180.0, 360.0) - 180.0

    tl1r  = g.proj_truelat1 * _RAD_PER_DEG
    ctl1r = cos(tl1r)

    tan_half_diff     = tan.((90.0 * g.proj_hemi .- lat) .* _RAD_PER_DEG ./ 2.0)
    tan_half_truelat1 = tan((90.0 * g.proj_hemi - g.proj_truelat1) * _RAD_PER_DEG / 2.0)

    rm  = g.proj_rebydx * ctl1r / g.proj_cone .*
          (tan_half_diff ./ tan_half_truelat1) .^ g.proj_cone

    arg = g.proj_cone .* (deltalon .* _RAD_PER_DEG)

    i = g.proj_polei .+ g.proj_hemi .* rm .* sin.(arg)
    j = g.proj_polej .-                rm .* cos.(arg)

    if g.proj_hemi == -1
        i = 2.0 .- i
        j = 2.0 .- j
    end
    return j, i   # (indy, indx)
end

function indyindxtolatlon(g::Lambert, j, i)
    inew = g.proj_hemi == -1 ? 2.0 - i : i
    jnew = g.proj_hemi == -1 ? 2.0 - j : j

    xx = inew - g.proj_polei
    yy = g.proj_polej - jnew
    r2 = xx^2 + yy^2
    r  = sqrt(r2) / g.proj_rebydx

    if r2 == 0
        lat = g.proj_hemi * 90.0
        lon = g.proj_stdlon
    else
        lon = g.proj_stdlon + _DEG_PER_RAD * atan(g.proj_hemi * xx, yy) / g.proj_cone
        lon = mod(lon + 360.0, 360.0)

        chi1 = (90.0 - g.proj_hemi * g.proj_truelat1) * _RAD_PER_DEG
        chi2 = (90.0 - g.proj_hemi * g.proj_truelat2) * _RAD_PER_DEG

        if chi1 == chi2
            chi = 2.0 * atan((r / tan(chi1))^(1.0/g.proj_cone) * tan(chi1/2.0))
        else
            chi = 2.0 * atan((r * g.proj_cone / sin(chi1))^(1.0/g.proj_cone) * tan(chi1/2.0))
        end
        lat = (90.0 - chi * _DEG_PER_RAD) * g.proj_hemi
    end
    return lat, lon
end

# Float32, type-stable twin of indyindxtolatlon(::Lambert) for the output
# mapping, which runs inside the KernelAbstractions deposit: the Float64 version
# returns a Float32 lon on the pole branch and Float64 otherwise (a Union in a
# GPU kernel) and computes in FP64 (slow or unavailable on GPU backends).
# Same construction; ~1e-5 deg precision, far below any output cell.
@inline function _lambert_latlon32(g::Lambert, j, i)
    hemi = g.proj_hemi
    inew = hemi == -1.0f0 ? 2.0f0 - Float32(i) : Float32(i)
    jnew = hemi == -1.0f0 ? 2.0f0 - Float32(j) : Float32(j)
    xx = inew - g.proj_polei
    yy = g.proj_polej - jnew
    r2 = xx * xx + yy * yy
    r2 == 0.0f0 && return hemi * 90.0f0, g.proj_stdlon
    r  = sqrt(r2) / g.proj_rebydx
    n  = g.proj_cone
    lon = g.proj_stdlon + (180.0f0 / Float32(π)) * atan(hemi * xx, yy) / n
    lon = mod(lon + 360.0f0, 360.0f0)
    chi1 = (90.0f0 - hemi * g.proj_truelat1) * (Float32(π) / 180.0f0)
    chi2 = (90.0f0 - hemi * g.proj_truelat2) * (Float32(π) / 180.0f0)
    chi  = chi1 == chi2 ?
        2.0f0 * atan((r / tan(chi1))^(1.0f0 / n) * tan(chi1 / 2.0f0)) :
        2.0f0 * atan((r * n / sin(chi1))^(1.0f0 / n) * tan(chi1 / 2.0f0))
    return (90.0f0 - chi * (180.0f0 / Float32(π))) * hemi, lon
end

"""
    lambert_mapfac(g::Lambert, indx, indy) -> m

Map scale factor of the Lambert conformal grid at a grid position: a grid step
covers DX/m on the ground.  With r the projected distance from the pole in
Earth radii and chi the colatitude (same construction as indyindxtolatlon),
m = n·r / sin(chi), n = proj_cone; m = 1 on the true latitude(s), < 1 between
two true latitudes, > 1 outside.  GPU-safe (Float32 arithmetic, no allocation).

Displacements used dx / DX until 2026-10-13, i.e. m = 1 everywhere: WRF winds
are physical m/s, so parcels moved by the factor m too slow or too fast
(-4 % at 45 N on a 30/60 N secant grid, and more towards the domain edges).
"""
@inline function lambert_mapfac(g::Lambert, indx, indy)
    i = Float32(indx); j = Float32(indy)
    inew = g.proj_hemi == -1 ? 2.0f0 - i : i
    jnew = g.proj_hemi == -1 ? 2.0f0 - j : j
    xx = inew - g.proj_polei
    yy = g.proj_polej - jnew
    r  = sqrt(xx * xx + yy * yy) / g.proj_rebydx
    r > 0.0f0 || return 1.0f0
    n    = g.proj_cone
    chi1 = (90.0f0 - g.proj_hemi * g.proj_truelat1) * Float32(π) / 180.0f0
    chi2 = (90.0f0 - g.proj_hemi * g.proj_truelat2) * Float32(π) / 180.0f0
    chi  = chi1 == chi2 ?
        2.0f0 * atan((r / tan(chi1))^(1.0f0 / n) * tan(chi1 / 2.0f0)) :
        2.0f0 * atan((r * n / sin(chi1))^(1.0f0 / n) * tan(chi1 / 2.0f0))
    return n * r / sin(chi)
end

@inline dxtodindx(g::Lambert, dx, indx, indy) = dx * lambert_mapfac(g, indx, indy) / g.dx
@inline dytodindy(g::Lambert, dy, indx, indy) = dy * lambert_mapfac(g, indx, indy) / g.dy

# Output grid on a Lambert nest ([output.grid].crs): geographic → through the
# parcel's lon/lat, coeffx = 1/dxout [1/deg], offsetx = lon0 [deg]; native →
# linear in projection metres.  Cell ii spans xout ∈ [ii, ii+1) either way.
@inline function indxtoxout(g::Lambert, indx, indy)
    g.geoout || return indx * g.coeffx + g.offsetx          # crs = "native" (metres)
    _, lon = _lambert_latlon32(g, indy, indx)
    return mod(lon - g.offsetx, 360.0f0) * g.coeffx + 1.0f0
end
@inline function indytoyout(g::Lambert, indx, indy)
    g.geoout || return indy * g.coeffy + g.offsety
    lat, _ = _lambert_latlon32(g, indy, indx)
    return (lat - g.offsety) * g.coeffy + 1.0f0
end

# Mercator (MesoNH conformal RPK = 0) ─────────────────────────────────────────
# Map y = R·cos(lat0)·ln(tan(π/4 + lat/2)); lon is linear in x.
# Scale factor k = cos(lat0)/cos(lat): map distance = k × ground distance.

const _MERC_R = 6.371f6

@inline _merc_y(coslat0::Float32, lat) =
    _MERC_R * coslat0 * log(tan(Float32(π)/4.0f0 + Float32(lat) * Float32(π)/360.0f0))

@inline function latlontoindyindx(g::Mercator, lat, lon)
    dl   = mod(Float32(lon) - g.lonor + 180.0f0, 360.0f0) - 180.0f0
    indx = dl * (Float32(π)/180.0f0) * _MERC_R * g.coslat0 / g.dx + 1.0f0
    indy = (_merc_y(g.coslat0, lat) - g.y_or) / g.dy + 1.0f0
    return indy, indx
end

@inline function indyindxtolatlon(g::Mercator, indy, indx)
    lon = g.lonor + (Float32(indx) - 1.0f0) * g.dx / (_MERC_R * g.coslat0) * (180.0f0/Float32(π))
    y   = g.y_or + (Float32(indy) - 1.0f0) * g.dy
    lat = 2.0f0 * atan(exp(y / (_MERC_R * g.coslat0))) * (180.0f0/Float32(π)) - 90.0f0
    return lat, lon
end

@inline function dxtodindx(g::Mercator, dx, indx, indy)
    lat, _ = indyindxtolatlon(g, indy, indx)
    return dx * (g.coslat0 / cosd(lat)) / g.dx
end

@inline function dytodindy(g::Mercator, dy, indx, indy)
    lat, _ = indyindxtolatlon(g, indy, indx)
    return dy * (g.coslat0 / cosd(lat)) / g.dy
end

# Output-grid mapping: exact lon/lat → fractional OUTGRID coordinate, cell ii
# spanning [ii, ii+1).  The "+ one met grid spacing" this used to add copied
# the Rectilinear offset error (fixed 2026-10-13): it shifted deposits by
# (met spacing / output spacing − 1) cells.
@inline function indxtoxout(g::Mercator, indx, indy)
    _, lon = indyindxtolatlon(g, indy, indx)
    return (lon - g.lon0out) / g.dxout + 1.0f0
end

@inline function indytoyout(g::Mercator, indx, indy)
    lat, _ = indyindxtolatlon(g, indy, indx)
    return (lat - g.lat0out) / g.dyout + 1.0f0
end

# ── Vertical coordinate transform ────────────────────────────────────────────

function ztoindl(lev::SigmaPressureLevel, indx, indy, z, indt, met3D)
    # met3D may be a MetFields3D struct (normal path) or a bare Array (legacy/test path).
    zlev = met3D isa MetFields3D ? met3D.zlev :
           met3D isa AbstractArray ? met3D :
           met3D["zlev"]
    indl = interpolation4D_irregzsingle(zlev, indx, indy, z, indt,
                                        Int32(1), lev.Nk)
    return indl === nothing ? 1.0f0 : indl
end

# ── Level-index identity transforms (sigma-pressure coords) ──────────────────
# For sigma-pressure levels the "level" coordinate IS the index, so these
# are identity functions. Other level types (HHL, pressure, theta) would
# override them with real inversions.

@inline levtoindl(lev::SigmaPressureLevel, indx, indy, indl, indt, met3D) = Float32(indl)
@inline indltolev(lev::SigmaPressureLevel, indx, indy, indl, indt, met3D) = Float32(indl)

# ── Output-grid helpers ────────────────────────────────────────────────────────

@inline indltozlev(lev::SigmaPressureLevel, matzlev, indx, indy, indl, indt) =
    interpolation4Dlevsingle(matzlev, indx, indy, indl, indt)
