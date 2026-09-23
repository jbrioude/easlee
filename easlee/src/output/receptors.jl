
# ── Gridless receptor output (ADAPTIVE_KERNEL companion) ─────────────────────
#
# Pointwise evaluation of the SAME adaptive-KDE estimator the deposit uses —
# no intermediate Eulerian grid, no cell-average representativeness error:
#
#   C(x_s) = (1/Δz) · Σ_{i : z_i ∈ layer} m_i · G₂(‖x_s − X_i‖; h_i)
#   G₂(r; h) = exp(−r²/2h²) / (2πh²)                       [kg m⁻³ directly]
#
# Two simplifications vs the grid path: NO renormalisation (a point evaluation
# estimates the field's value, it does not partition mass) and NO bandwidth
# ambiguity (sample-point estimator: h belongs to the CONTRIBUTING parcel —
# the same choice that makes the grid deposit conservative; one meaning of h
# everywhere).  R(x_s) = Σ m·w·h / Σ m·w comes out of the same sum: a station
# report reads "C = 340 pg m⁻³, resolved at 15 km".
#
# Time semantics follow FLEXPART unchanged: evaluate at the sample cadence,
# average over the output window, one series per station.  Support truncated
# at 3h_i like the deposit.  Stations are few and static — a plain per-parcel
# loop over stations with a bounding-box reject is O(N·S_bbox) and negligible
# next to the deposit.

mutable struct ReceptorSet
    ids  :: Vector{String}
    lon  :: Vector{Float64}
    lat  :: Vector{Float64}
    z1   :: Float64                  # sampling layer [m], same z convention as
    z2   :: Float64                  # the output grid (find_out_layer heights)
    # Physical bandwidth floor [m].  The kNN hband measures STATISTICAL
    # spacing; near a fresh point release the parcels sit in a filament far
    # thinner than any parcel's true subgrid footprint after hours of
    # turbulent transport, and a point evaluation there reports the
    # unphysically thin filament (ETEX station 8021, 20 km from the source:
    # h ≈ 60 m → 20× the observed 3-h mean while the timing was exactly
    # right).  Flooring h at a footprint scale (default: half the output
    # grid spacing) bounds the point-evaluation variance; set 0 for the
    # pure sample-point estimator.
    hmin :: Float64
    # Physical (relative-dispersion) bandwidth floor: σv [m/s] and T_L [s] of
    # the Taylor growth (_hfloor_taylor in conc.jl), evaluated per parcel at
    # its AGE.  fsv = 0 disables.  Unlike hmin (a fixed representativeness
    # scale), this is the width of atmosphere the parcel physically stands
    # for — it also covers parcels born after the last kNN refresh (hband 0).
    fsv  :: Float64
    ftl  :: Float64
    # true → read the flow-evolved covariance (state cxx/cxy/cyy, kernel =
    # "covariance") and evaluate the FULL bivariate Gaussian
    # exp(−mahal/2)/(2π√detP) — the receptor honours the cxy tilt the
    # separable grid deposit drops.  false → isotropic hband KDE.
    usecov :: Bool
    # shape mode (usecov && shape): rescale the tilted P to det = h⁴ (h = the
    # kNN hband) — the statistical area with the physical orientation; the
    # kernel size then vanishes with N (no double counting of dispersion).
    shape  :: Bool
    # within-window accumulators (mean over nsamp at flush)
    csum :: Vector{Float64}          # Σ_samples C(x_s)
    rnum :: Vector{Float64}          # Σ m·w·h   (across samples — mass-weighted R)
    rden :: Vector{Float64}          # Σ m·w
    nsamp:: Int
    # flushed rows: (t_s, station index, conc, R, nsamp)
    rows :: Vector{NTuple{5,Float64}}
end

function ReceptorSet(ids::Vector{String}, lon::Vector{Float64}, lat::Vector{Float64},
                     z1::Real, z2::Real; hmin::Real = 0.0,
                     fsv::Real = 0.0, ftl::Real = 300.0, usecov::Bool = false,
                     shape::Bool = false)
    n = length(ids)
    ReceptorSet(ids, lon, lat, Float64(z1), Float64(z2), Float64(hmin),
                Float64(fsv), Float64(ftl), usecov, shape,
                zeros(n), zeros(n), zeros(n), 0, NTuple{5,Float64}[])
end

"Load a DATEM station list (`id lat lon` per line, comments skipped)."
function receptorset_from_file(path::AbstractString, z1::Real, z2::Real; hmin::Real = 0.0,
                               fsv::Real = 0.0, ftl::Real = 300.0, usecov::Bool = false,
                               shape::Bool = false)
    ids = String[]; lon = Float64[]; lat = Float64[]
    for ln in eachline(path)
        occursin(r"^\s*\d", ln) || continue
        f = split(ln); length(f) >= 3 || continue
        push!(ids, String(f[1]))
        push!(lat, parse(Float64, f[2])); push!(lon, parse(Float64, f[3]))
    end
    isempty(ids) && error("no stations parsed from $path (expected: id lat lon)")
    ReceptorSet(ids, lon, lat, z1, z2; hmin=hmin, fsv=fsv, ftl=ftl, usecov=usecov,
                shape=shape)
end

"""
    receptor_sample!(rs, state, gconfig, lconfig, zlev_mat, indt, nn) → n_hits

One KDE evaluation of every station (call at the output SAMPLE cadence, after
the deposit pass so `state.hband` is fresh).  Parcel loop threaded with
per-thread station buffers (stations are few).
"""
function receptor_sample!(rs::ReceptorSet, state::State.ParcelState,
                          gconfig, lconfig, zlev_mat, indt, nn::Int)
    ns = length(rs.ids)
    ns == 0 && return 0
    nn8 = Int8(nn)
    dz = rs.z2 - rs.z1
    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    cbuf = [zeros(Float64, ns) for _ in 1:nthr]
    nbuf = [zeros(Float64, ns) for _ in 1:nthr]
    dbuf = [zeros(Float64, ns) for _ in 1:nthr]
    hits = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in eachindex(state.status)
        (state.status[i] == 0x01 && state.nest[i] == nn8) || continue
        z = indltozlev(lconfig, zlev_mat, state.indx[i], state.indy[i],
                       state.indl[i], Float32(indt))
        (rs.z1 <= z < rs.z2) || continue
        la, lo = indyindxtolatlon(gconfig, Float64(state.indy[i]), Float64(state.indx[i]))
        m = Float64(state.mass[i])
        tid = Threads.threadid()
        cb = cbuf[tid]; nb = nbuf[tid]; db = dbuf[tid]
        if rs.usecov
            # flow-evolved covariance: full bivariate Gaussian (with the cxy
            # tilt — a sheared filament passes a station AT ITS WIDTH along
            # the station's direction, not the isotropic mean).  hmin acts as
            # a per-axis variance floor.
            pa = max(Float64(state.cxx[i]), 1.0, rs.hmin^2)
            pd = max(Float64(state.cyy[i]), 1.0, rs.hmin^2)
            pb = Float64(state.cxy[i])
            blim = 0.99 * sqrt(pa * pd)
            pb = clamp(pb, -blim, blim)
            pdet = pa * pd - pb * pb
            if rs.shape
                # statistical area, physical orientation: P → P·h²/√det
                h = max(Float64(state.hband[i]), 1.0, rs.hmin)
                sc = h * h / sqrt(pdet)
                pa *= sc; pd *= sc; pb *= sc
                pdet = pa * pd - pb * pb          # = h⁴ (up to roundoff)
            end
            r3x = 3.0 * sqrt(pa); r3y = 3.0 * sqrt(pd)
            heff = pdet^0.25                       # area-equivalent width for R
            @inbounds for s in 1:ns
                dy = (la - rs.lat[s]) * 111_320.0
                abs(dy) > r3y && continue
                dlon = mod(lo - rs.lon[s] + 180.0, 360.0) - 180.0
                dx = dlon * 111_320.0 * cosd(rs.lat[s])
                abs(dx) > r3x && continue
                mahal = (pd * dx * dx - 2.0 * pb * dx * dy + pa * dy * dy) / pdet
                mahal > 9.0 && continue            # 3σ ellipse (matches 3h cut)
                w = exp(-0.5 * mahal) / (2π * sqrt(pdet))   # G₂ [m⁻²]
                cb[s] += m * w / dz
                nb[s] += m * w * heff
                db[s] += m * w
                Threads.atomic_add!(hits, 1)
            end
            continue
        end
        h = max(Float64(state.hband[i]), 1.0, rs.hmin)
        if rs.fsv > 0.0
            hf = _hfloor_taylor(_parcel_age_s(state, i), rs.fsv, rs.ftl)
            h < hf && (h = hf)
        end
        r3 = 3.0 * h
        @inbounds for s in 1:ns
            dy = (la - rs.lat[s]) * 111_320.0
            abs(dy) > r3 && continue
            # wrap the lon difference into (−180, 180]: parcel lons come back
            # in the MET frame (0–360 on the CAPTEX crop: −84.05 ≡ 275.95)
            dlon = mod(lo - rs.lon[s] + 180.0, 360.0) - 180.0
            dx = dlon * 111_320.0 * cosd(rs.lat[s])
            abs(dx) > r3 && continue
            r2 = dx * dx + dy * dy
            r2 > r3 * r3 && continue
            w = exp(-r2 / (2.0 * h * h)) / (2π * h * h)     # G₂ [m⁻²]
            cb[s] += m * w / dz
            nb[s] += m * w * h
            db[s] += m * w
            Threads.atomic_add!(hits, 1)
        end
    end
    for t in 1:nthr
        rs.csum .+= cbuf[t]
        rs.rnum .+= nbuf[t]
        rs.rden .+= dbuf[t]
    end
    rs.nsamp += 1
    return hits[]
end

"Freeze the window means into rows (called at each output-interval boundary)."
function receptor_flush!(rs::ReceptorSet, t_s::Real)
    ns = max(rs.nsamp, 1)
    for s in eachindex(rs.ids)
        c = rs.csum[s] / ns
        R = rs.rden[s] > 0.0 ? rs.rnum[s] / rs.rden[s] : NaN
        push!(rs.rows, (Float64(t_s), Float64(s), c, R, Float64(rs.nsamp)))
    end
    fill!(rs.csum, 0.0); fill!(rs.rnum, 0.0); fill!(rs.rden, 0.0)
    rs.nsamp = 0
    return nothing
end

"Write the accumulated station time series (CSV: time, id, conc, R, nsamp)."
function receptor_write(rs::ReceptorSet, path::AbstractString, start)
    open(path, "w") do io
        println(io, "# gridless receptor time series — C at the station coordinate")
        println(io, "# at the model's native resolution R (mass-weighted mean kernel")
        println(io, "# width of the contributing parcels); start = $start")
        println(io, "time_s,station,conc_kgm3,resolution_m,nsamp")
        for (t, s, c, R, n) in rs.rows
            println(io, Int(t), ",", rs.ids[Int(s)], ",", c, ",",
                    isnan(R) ? "NaN" : round(R, digits=1), ",", Int(n))
        end
    end
    return path
end
