# physics/track_nudge.jl — observed-track ensemble nudging + wind-bias drift
# (Hunga Tonga thread: put the stratospheric cloud where the observations say
# while preserving its internal structure — prerequisite for the Vömel et al.
# radiosonde validation over Australia).
#
# Two corrections, both applied to PARCEL MOTION (the archived met is never
# modified), both confined to a height band and both TOML-driven:
#
#   • TrackNudge ([physics.nudge]): each sync step, compute the horizontal
#     CENTROID of the active in-band parcels and relax it toward the
#     time-interpolated observed track centroid with timescale τ —
#     Δx = (x_obs − x_ens)·(1 − e^(−dt/τ)) — applied as ONE RIGID SHIFT to
#     every in-band parcel.  Internal structure and relative dispersion are
#     untouched; only the advection bias is corrected.  Outside the track's
#     time span the run is free (no clamping, no extrapolation).
#     Track file = the ht_cloud_envelope.csv format:
#         time, lon_west, lon_east, lat_south, lat_north[, …]
#     (ISO times; lons in any convention — unwrapped internally).
#
#   • WindBias ([physics.wind_bias]): constant per-parcel drift (u, v) [m/s]
#     inside the band — the manual stratospheric-wind tuning knob for the
#     free-running leg (e.g. ERA5 meridional bias sensitivity).
#
# Longitude arithmetic is done in an UNWRAPPED frame anchored on the track
# (dateline-safe: the HT cloud crosses ±180 within hours).

struct TrackNudge
    t_ms  :: Vector{Int64}       # absolute DateTime value [ms], ascending
    clon  :: Vector{Float64}     # observed centroid lon, UNWRAPPED continuous
    clat  :: Vector{Float64}
    tau   :: Float64             # relaxation time [s]
    z0    :: Float64             # height band [m]
    z1    :: Float64
end

"""
    load_track_nudge(file, tau_s, z) → TrackNudge

Parse the observed-track CSV (`time, lon_west, lon_east, lat_south,
lat_north`).  Envelope centres become the track centroids; longitudes are
unwrapped so a dateline-crossing track stays continuous.
"""
function load_track_nudge(file::AbstractString, tau_s::Real, z::AbstractVector)
    t = Int64[]; lo = Float64[]; la = Float64[]
    for ln in eachline(file)
        s = strip(ln)
        (isempty(s) || startswith(s, "#") || startswith(lowercase(s), "time")) && continue
        f = split(s, ",")
        length(f) >= 5 || continue
        push!(t, Dates.value(DateTime(strip(f[1]))))
        w = parse(Float64, f[2]); e = parse(Float64, f[3])
        # centre of a possibly wrap-straddling box: take the shorter arc
        de = mod(e - w + 180.0, 360.0) - 180.0
        push!(lo, w + 0.5 * de)
        push!(la, 0.5 * (parse(Float64, f[4]) + parse(Float64, f[5])))
    end
    length(t) >= 2 || error("[physics.nudge] track $file: need ≥ 2 rows (time, lon_west, lon_east, lat_south, lat_north)")
    p = sortperm(t)
    t = t[p]; lo = lo[p]; la = la[p]
    # unwrap lon so the track is continuous across ±180
    for k in 2:length(lo)
        lo[k] = lo[k-1] + (mod(lo[k] - lo[k-1] + 180.0, 360.0) - 180.0)
    end
    return TrackNudge(t, lo, la, Float64(tau_s), Float64(z[1]), Float64(z[end]))
end

# observed centroid at date (nothing outside the track span)
function _track_at(tn::TrackNudge, date::DateTime)
    tm = Dates.value(date)
    (tm < tn.t_ms[1] || tm > tn.t_ms[end]) && return nothing
    k = searchsortedlast(tn.t_ms, tm)
    k >= length(tn.t_ms) && return (tn.clon[end], tn.clat[end])
    f = (tm - tn.t_ms[k]) / max(tn.t_ms[k+1] - tn.t_ms[k], 1)
    return ((1 - f) * tn.clon[k] + f * tn.clon[k+1],
            (1 - f) * tn.clat[k] + f * tn.clat[k+1])
end

"""
    apply_track_nudge!(state, tn, gconfig, lconfig, zlev, indt, date, dt_s)
        → (n_moved, dlon, dlat)

One nudge step: rigid-shift the in-band ensemble toward the observed track
centroid.  Returns what was applied (degrees) for logging; (0, 0, 0) when the
date is outside the track or no parcel is in the band.
"""
function apply_track_nudge!(state::State.ParcelState, tn::TrackNudge,
                            gconfig, lconfig, zlev, indt, date::DateTime,
                            dt_s::Float64)
    obs = _track_at(tn, date)
    obs === nothing && return (0, 0.0, 0.0)
    olon, olat = obs

    # in-band ensemble centroid, lon unwrapped around the OBSERVED centroid
    idx = Int[]
    slon = 0.0; slat = 0.0
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        z = indltozlev(lconfig, zlev, state.indx[i], state.indy[i],
                       state.indl[i], Float32(indt))
        (tn.z0 <= z <= tn.z1) || continue
        la, lo = indyindxtolatlon(gconfig, Float64(state.indy[i]), Float64(state.indx[i]))
        push!(idx, i)
        slon += olon + (mod(lo - olon + 180.0, 360.0) - 180.0)
        slat += la
    end
    isempty(idx) && return (0, 0.0, 0.0)
    elon = slon / length(idx); elat = slat / length(idx)

    f = 1.0 - exp(-dt_s / max(tn.tau, 1.0))
    dlon = (olon - elon) * f
    dlat = (olat - elat) * f
    _shift_parcels!(state, idx, gconfig, dlon, dlat)
    return (length(idx), dlon, dlat)
end

# rigid horizontal shift in grid-index space (Rectilinear: affine lon/lat map)
function _shift_parcels!(state, idx::Vector{Int}, g::Rectilinear,
                         dlon::Float64, dlat::Float64)
    dix = Float32(dlon / Float64(g.dlon))
    diy = Float32(dlat / Float64(g.dlat))
    nx  = Float32(g.Ni)
    @inbounds for i in idx
        x = state.indx[i] + dix
        g.globalx && (x = mod(x - 1.0f0, nx) + 1.0f0)   # wrap on global met
        state.indx[i] = x
        state.indy[i] = state.indy[i] + diy
    end
    return nothing
end

# ── Wind-bias drift ───────────────────────────────────────────────────────────

struct WindBias
    u  :: Float64        # [m/s] eastward drift added inside the band
    v  :: Float64        # [m/s] northward
    z0 :: Float64
    z1 :: Float64
end

"""
    apply_wind_bias!(state, wb, gconfig, lconfig, zlev, indt, dt_s) → n_moved

Constant per-parcel drift inside the height band: Δlon = u·dt/(R_E·cosφ),
Δlat = v·dt/R_E (per parcel — each at its own latitude).
"""
function apply_wind_bias!(state::State.ParcelState, wb::WindBias,
                          gconfig::Rectilinear, lconfig, zlev, indt,
                          dt_s::Float64)
    (wb.u == 0.0 && wb.v == 0.0) && return 0
    dlat = wb.v * dt_s / 111320.0
    diy  = Float32(dlat / Float64(gconfig.dlat))
    nx   = Float32(gconfig.Ni)
    n = 0
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        z = indltozlev(lconfig, zlev, state.indx[i], state.indy[i],
                       state.indl[i], Float32(indt))
        (wb.z0 <= z <= wb.z1) || continue
        la, _ = indyindxtolatlon(gconfig, Float64(state.indy[i]), Float64(state.indx[i]))
        dlon = wb.u * dt_s / (111320.0 * max(cosd(la), 0.05))
        x = state.indx[i] + Float32(dlon / Float64(gconfig.dlon))
        gconfig.globalx && (x = mod(x - 1.0f0, nx) + 1.0f0)
        state.indx[i] = x
        state.indy[i] = state.indy[i] + diy
        n += 1
    end
    return n
end
