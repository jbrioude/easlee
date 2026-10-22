
# ── Box-restricted domain fill (release kind = "domainfill") ──────────────────
#
# Engine form of the FLEXPART DOMAINFILL mode over a lon/lat/z sub-box (the
# [[release]] box), generalised from the chem_cams_pacific.jl harness:
#
#   1. domainfill_init!   — place nparts parcels at t=0: uniform in lon/lat over
#      the box, model layers drawn ∝ air mass (ρ·Δz), each parcel carrying
#      mass_per_parcel = box air mass / nparts.
#   2. domainfill_inject! — each step, inject parcels along the 4 box edges
#      proportional to the incoming mass flux ρ·v_in·width·Δz·Δt.  The library
#      inject_boundary_parcels! sweeps the WHOLE met boundary — on global met
#      that injects around the entire globe (then most are culled), the dominant
#      per-step cost; this sweeps only the box faces.
#   3. domainfill_cull!   — free parcels that left the box (horizontally, or
#      above the fill-layer top), keeping the population ≈ nparts and the slot
#      pool bounded.
#
# All domain-fill parcels carry relbox = 0 (the boundary-injection sentinel):
# conc.jl clamps 0 → box 1 for gridded output, and the cull uses it to leave
# ordinary release parcels untouched.  One domainfill box per run.
#
# The optional `seeded` BitVector (chemistry bookkeeping in run_spec!) has the
# bit cleared for every slot this code writes or frees, so init="cams" re-seeds
# inflow parcels with the boundary composition.
#
# Vertical layout (contract 2026-10-22): layer k is [zlev[k], zlev[k+1]) — level
# indices [k, k+1), where its parcels are placed — and holds full level k+1.
# Its air mass is ρ[k+1]·(zlev[k+1] − zlev[k]) and its inflow uses the wind of
# full level k+1, both read with the stacked interpolation at index k + 1/2.
# Until 2026-10-22 the weight was ρ[k]·Δz centred on half level k (ρ and Δz
# half a level below the layer the parcels went into) and the face reads were
# the plain interpolation at k + 1/2 (the mean of slots k and k+1).

struct DomainFillBox
    lon0::Float64; lon1::Float64        # continuous frame, lon1 > lon0 (dateline-crossing OK)
    lat0::Float64; lat1::Float64
    ztop::Float64                       # fill-layer top [m]
    keep_lev::Vector{Int}               # layers [k, k+1) whose full level (mid height) is ≤ ztop
    lev_w::Vector{Float64}              # air-mass weight (ρ·Δz) per keep_lev entry
    lev_h::Vector{Float64}              # height of EVERY half level k = zlev[k] (box-centre column)
    mass_per_parcel::Float32
    nparts::Int
end

# Wrap a met-frame longitude into the box's continuous [lon0, lon0+360) frame.
@inline _df_wraplon(df::DomainFillBox, lon) = mod(lon - df.lon0, 360.0) + df.lon0

# Met-grid spacing in degrees for the box sweeps (edge faces + mass integral).
@inline _df_dlon(gc) = gc isa Rectilinear ? Float64(gc.dlon) : 0.5
@inline _df_dlat(gc) = gc isa Rectilinear ? Float64(gc.dlat) : 0.5

const _DF_MDEG = 111_320.0              # metres per degree (meridional)

# Air density of layer k = full level k+1: the stacked read at index k + 1/2.
@inline _df_rho(rhoarr, fx::Float32, fy::Float32, k::Int, indt::Float32) =
    interpolation4Dsingle(rhoarr, fx, fy, Float32(k) + 0.5f0, indt)

# First free slot after `after` (monotone scan — callers thread the cursor).
@inline function _df_next_free(state, after::Int)
    i = after
    @inbounds while i < length(state.status)
        i += 1
        state.status[i] == 0x00 && return i
    end
    return 0
end

"""
    build_domainfill(nestctx, indt, lon, lat, ztop, nparts) → DomainFillBox

Layer weights (box-centre column) + box air mass → mass_per_parcel.  Call after
the first `read_met_step!` so ρ/zlev are loaded.
"""
function build_domainfill(nestctx, indt::Float32, lon::AbstractVector, lat::AbstractVector,
                          ztop::Real, nparts::Int)
    gc  = nestctx.gridconfig
    met = nestctx.metfields3D
    Nk  = Int(nestctx.levelconfig.Nk)
    lon0, lon1 = Float64(lon[1]), Float64(lon[2])
    lon1 <= lon0 && (lon1 += 360.0)               # dateline-crossing box
    lat0, lat1 = min(lat[1], lat[2]), max(lat[1], lat[2])

    # per-layer ρ·Δz weight + half-level heights at the box-centre column
    indy, indx = latlontoindyindx(gc, 0.5 * (lat0 + lat1), 0.5 * (lon0 + lon1))
    fx, fy = Float32(indx), Float32(indy)
    w = zeros(Float64, Nk - 1); h = zeros(Float64, Nk)
    for k in 1:Nk
        zk = interpolation4Dlevsingle(met.zlev, fx, fy, Float32(k), indt)
        h[k] = isfinite(zk) ? zk : 1.0e9
    end
    for k in 1:(Nk - 1)
        ρ  = _df_rho(met.rho, fx, fy, k, indt)
        dz = _dz_lev(met.zlev, fx, fy, k, indt)
        (isfinite(ρ) && ρ > 0 && isfinite(dz)) && (w[k] = ρ * dz)
    end
    # a layer belongs to the fill when its full level (mid height) is ≤ ztop
    keep = [k for k in 1:(Nk - 1) if 0.5 * (h[k] + h[k+1]) <= ztop && w[k] > 0]
    isempty(keep) && error("domainfill: no model layers below ztop = $ztop m")

    # box air mass: ρ·area·Δz over met-resolution cells × kept layers
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    tot = 0.0
    latq = lat0 + dlat / 2
    while latq < lat1
        area = (dlat * _DF_MDEG) * (dlon * _DF_MDEG * cosd(latq))
        lonq = lon0 + dlon / 2
        while lonq < lon1
            iy, ix = latlontoindyindx(gc, latq, lonq)
            qx, qy = Float32(ix), Float32(iy)
            for k in keep
                ρ  = _df_rho(met.rho, qx, qy, k, indt)
                dz = _dz_lev(met.zlev, qx, qy, k, indt)
                (isfinite(ρ) && ρ > 0 && isfinite(dz)) && (tot += ρ * area * dz)
            end
            lonq += dlon
        end
        latq += dlat
    end
    tot > 0 || error("domainfill: box air mass is zero — box outside the met domain?")
    return DomainFillBox(lon0, lon1, lat0, lat1, Float64(ztop), keep, w[keep], h,
                         Float32(tot / nparts), nparts)
end

"""
    domainfill_init!(state, ctx, df, date, indt; nest=1, rng, seeded=nothing) → n placed

Place the initial fill: uniform in lon/lat over the box, layers ∝ air mass
(a parcel of layer k at level index k + U[0, 1)).
Parcels have `is_new = false` and `Dtime = date`, so the advance picks them up
in the step at `date` and `initialize_new_parcels!` leaves them alone.
"""
function domainfill_init!(state::State.ParcelState, ctx, df::DomainFillBox,
                          date::DateTime, indt::Float32;
                          nest::Int = 1, rng = Random.default_rng(), seeded = nothing)
    gc = ctx.nests[nest].gridconfig
    # integer parcel count per layer ∝ weight; rounding residual → heaviest layer
    wsum = sum(df.lev_w)
    counts = [round(Int, wk / wsum * df.nparts) for wk in df.lev_w]
    counts[argmax(df.lev_w)] += df.nparts - sum(counts)
    nplaced = 0; slot = 0
    for (li, k) in enumerate(df.keep_lev)
        for _ in 1:counts[li]
            slot = _df_next_free(state, slot)
            slot == 0 && return nplaced                  # pool exhausted
            lon = df.lon0 + rand(rng) * (df.lon1 - df.lon0)
            lat = df.lat0 + rand(rng) * (df.lat1 - df.lat0)
            indy, indx = latlontoindyindx(gc, lat, lon)
            _write_boundary_parcel!(state, slot, nest, Float32(indx), Float32(indy),
                                    Float32(k) + rand(rng, Float32),
                                    df.mass_per_parcel, date)
            seeded === nothing || (seeded[slot] = false)
            nplaced += 1
        end
    end
    return nplaced
end

"""
    domainfill_inject!(state, ctx, df, currentdate, indt, dt_s; nest=1, rng, seeded)
        → n injected

Inject parcels along the 4 box edges ∝ the incoming mass flux
ρ·v_in·width·Δz·Δt (stochastic rounding of the fractional part).  Call BEFORE
the advance so injected parcels (Dtime = currentdate) move in the same step.
W/E faces are skipped for a zonally-periodic box; S/N faces at the poles.
"""
function domainfill_inject!(state::State.ParcelState, ctx, df::DomainFillBox,
                            currentdate::DateTime, indt::Float32, dt_s::Real;
                            nest::Int = 1, rng = Random.default_rng(), seeded = nothing)
    nestctx = ctx.nests[nest]
    gc  = nestctx.gridconfig
    met = nestctx.metfields3D
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    dtf  = Float32(dt_s)
    ninj = 0; slot = 0

    place! = (lon, lat, k) -> begin
        slot = _df_next_free(state, slot)
        slot == 0 && return false
        indy, indx = latlontoindyindx(gc, lat, lon)
        _write_boundary_parcel!(state, slot, nest, Float32(indx), Float32(indy),
                                Float32(k) + rand(rng, Float32),
                                df.mass_per_parcel, currentdate)
        seeded === nothing || (seeded[slot] = false)
        ninj += 1
        return true
    end

    # West (inflow U>0) / East (inflow U<0) — skipped for a zonally-periodic box
    if (df.lon1 - df.lon0) < 360.0 - dlon
        dy_m = Float32(dlat * _DF_MDEG)
        for (lonf, usign) in ((df.lon0, 1.0f0), (df.lon1, -1.0f0))
            latq = df.lat0 + dlat / 2
            while latq < df.lat1
                indy, indx = latlontoindyindx(gc, latq, lonf)
                fx, fy = Float32(indx), Float32(indy)
                for k in df.keep_lev
                    kf = Float32(k) + 0.5f0          # layer mid index → full level k+1
                    u = interpolateU(gc.Uedge, met.Ux, fx, fy, kf, indt) * usign
                    (isfinite(u) && u > 0.0f0) || continue
                    ρ = _df_rho(met.rho, fx, fy, k, indt)
                    (isfinite(ρ) && ρ > 0.0f0) || continue
                    nr = (ρ * u * dy_m * _dz_lev(met.zlev, fx, fy, k, indt) * dtf) /
                         df.mass_per_parcel
                    isfinite(nr) || continue
                    n = floor(Int, nr)
                    rand(rng, Float32) < (nr - Float32(n)) && (n += 1)
                    for _ in 1:n
                        jlat = clamp(latq + (rand(rng) - 0.5) * dlat,
                                     df.lat0 + 1e-3, df.lat1 - 1e-3)
                        place!(lonf + usign * 0.05 * dlon, jlat, k) || return ninj
                    end
                end
                latq += dlat
            end
        end
    end
    # South (inflow V>0) / North (inflow V<0)
    for (latf, vsign) in ((df.lat0, 1.0f0), (df.lat1, -1.0f0))
        abs(latf) >= 89.9 && continue
        dx_m = Float32(dlon * _DF_MDEG * cosd(latf))
        lonq = df.lon0 + dlon / 2
        while lonq < df.lon1
            indy, indx = latlontoindyindx(gc, latf, lonq)
            fx, fy = Float32(indx), Float32(indy)
            for k in df.keep_lev
                kf = Float32(k) + 0.5f0              # layer mid index → full level k+1
                v = interpolateV(gc.Vedge, met.Vy, fx, fy, kf, indt) * vsign
                (isfinite(v) && v > 0.0f0) || continue
                ρ = _df_rho(met.rho, fx, fy, k, indt)
                (isfinite(ρ) && ρ > 0.0f0) || continue
                nr = (ρ * v * dx_m * _dz_lev(met.zlev, fx, fy, k, indt) * dtf) /
                     df.mass_per_parcel
                isfinite(nr) || continue
                n = floor(Int, nr)
                rand(rng, Float32) < (nr - Float32(n)) && (n += 1)
                for _ in 1:n
                    jlon = clamp(lonq + (rand(rng) - 0.5) * dlon,
                                 df.lon0 + 1e-3, df.lon1 - 1e-3)
                    place!(jlon, latf + vsign * 0.05 * dlat, k) || return ninj
                end
            end
            lonq += dlon
        end
    end
    return ninj
end

"""
    domainfill_cull!(state, ctx, df, indt; nest=1, seeded=nothing) → n freed

Free domain-fill parcels (relbox = 0) that left the box horizontally or rose
above the fill-layer top (5 % + 100 m margin).  Ordinary release parcels
(relbox ≥ 1) are never touched.
"""
function domainfill_cull!(state::State.ParcelState, ctx, df::DomainFillBox,
                          indt::Float32; nest::Int = 1, seeded = nothing)
    nestctx = ctx.nests[nest]
    gc  = nestctx.gridconfig
    met = nestctx.metfields3D
    ztop_m = df.ztop * 1.05 + 100.0
    nfreed = 0
    @inbounds for i in eachindex(state.status)
        (state.status[i] == 0x01 && state.relbox[i] == 0) || continue
        lat, lon = indyindxtolatlon(gc, Float64(state.indy[i]), Float64(state.indx[i]))
        lonw = _df_wraplon(df, lon)
        inbox = (df.lon0 <= lonw <= df.lon1) && (df.lat0 <= lat <= df.lat1)
        if inbox
            z = interpolation4Dlevsingle(met.zlev, state.indx[i], state.indy[i],
                                         state.indl[i], indt)
            inbox = !(isfinite(z) && z > ztop_m)
        end
        if !inbox
            state.status[i] = 0x00
            seeded === nothing || (seeded[i] = false)
            nfreed += 1
        end
    end
    return nfreed
end
