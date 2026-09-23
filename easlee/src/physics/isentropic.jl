
# ── Hybrid kinematic/isentropic vertical transport ────────────────────────────
#
# In the stratosphere, kinematic vertical velocities (etadot/ω) are noisy and
# produce spurious vertical dispersion; transport there is far more accurate
# on isentropes with an explicit diabatic drift (CLaMS, MPTRAC).  This module
# implements the MPTRAC-style hybrid (after Mahowald et al. 2002):
#
#   • each parcel carries a PROGNOSTIC potential temperature θ_p;
#   • the diabatic vertical motion is  θ̇ = (θ/T)·(∂T/∂t)_rad  from the ERA5
#     archived radiative tendencies (avg_ttlwr + avg_ttswr, K/s — the leading
#     diabatic term above the tropopause);
#   • after the kinematic advection step, the parcel's vertical position is
#     blended between the kinematic result and the isentropic target
#     (θ_p + θ̇·Δt, inverted through the local θ ladder):
#
#         indl ← indl_kin + w(σ)·(indl_θ − indl_kin)
#
#     with a smooth weight w(σ): 0 for σ = p/p_s > σ_r (troposphere, pure
#     kinematic), 1 for σ < σ_t (stratosphere, pure isentropic-diabatic),
#     cosine ramp in between.  Both limits are mass-consistent and the blend
#     is smooth in space and time, so the well-mixed criterion holds through
#     the tropopause region — the same argument as MPTRAC's ζ coordinate.
#
# In the w = 1 region θ is conserved EXACTLY between diabatic increments
# (the prognostic θ_p is not re-interpolated from the grid), eliminating the
# numerical vertical diffusion of repeated grid↔parcel conversions — the
# CLaMS property that makes isentropic transport attractive for the
# stratosphere.
#
# Requirements: metfields3D theta/T/puv (always present) + the radiative
# tendency fields mapped in the namelist (ERA5_CONV.vars):
#     ttlwr = avg_ttlwr,  ttswr = avg_ttswr

# Diabatic-tendency arrays for one nest (shared refs with the met dict).
# Preferred: the TOTAL parameterised temperature tendency (namelist
# `ttpm = mttpm` — radiation + latent + diffusion); fallback: LW + SW
# radiative tendencies (`ttlwr = avg_ttlwr`, `ttswr = avg_ttswr`).
"""
Diabatic temperature-tendency fields [K/s] of one nest, feeding the θ̇ drift
in `hybrid_vertical!` (tt2 === nothing when the total tendency `ttpm` is
archived; otherwise tt1+tt2 = LW+SW radiative tendencies).
"""
struct IsentropicFields{A<:AbstractArray{Float32,4}}
    tt1::A
    tt2::Union{A,Nothing}
end

function IsentropicFields(d::AbstractDict)
    if haskey(d, "ttpm")
        A = typeof(d["ttpm"])
        return IsentropicFields{A}(d["ttpm"], nothing)
    end
    (haskey(d, "ttlwr") && haskey(d, "ttswr")) || return nothing
    A = typeof(d["ttlwr"])
    IsentropicFields{A}(d["ttlwr"], d["ttswr"])
end

# ── Isentropic-level winds (ERA5 IFS_PT_MET_* files) ──────────────────────────
# u/v archived ON θ-surfaces.  These — not model-level winds vertically
# interpolated to the parcel's isentrope — must drive the horizontal motion in
# the isentropic regime: the quadrilinear model-level interpolation does its
# horizontal part at constant LEVEL, blending winds across isentropes wherever
# θ has a horizontal gradient on the level (vortex edge, jets).  The PT files
# interpolate along the isentrope on the native IFS grid instead.
mutable struct PTWinds
    u::Array{Float32,4}             # (Ni, Nj, Nθ, 2) — lat ascending
    v::Array{Float32,4}
    thetalev::Vector{Float32}       # ascending [K]
    lon0::Float32; dlon::Float32; Ni::Int
    lat0::Float32; dlat::Float32; Nj::Int
    filedates::Vector{DateTime}
    filenames::Vector{String}
    t1::Int; t2::Int                # loaded bracket (0 = none)
end

function PTWinds(directory::AbstractString; pattern::AbstractString="IFS_PT_MET")
    files = sort(filter(f -> occursin(pattern, basename(f)),
                        readdir(directory, join=true)))
    isempty(files) && return nothing
    dates = DateTime[]
    for f in files
        g = GribFile(f); m = read(g)
        push!(dates, DateTime(string(m["dataDate"]), dateformat"yyyymmdd") +
                     Dates.Hour(div(Int(m["dataTime"]), 100)))
        GRIB.destroy(g)
    end
    ord   = sortperm(dates)
    files = files[ord]; dates = dates[ord]
    # geometry + θ ladder from the first file
    levs = Float64[]
    local Ni, Nj, lon0, lat1, dlon, dlat
    GribFile(files[1]) do g
        for msg in g
            m = Message(msg)
            m["shortName"] == "u" || continue
            if isempty(levs)
                Ni   = Int(m["Ni"]); Nj = Int(m["Nj"])
                lons = reshape(m["longitudes"], Ni, Nj)
                lats = reshape(m["latitudes"],  Ni, Nj)
                lon0 = Float32(lons[1, 1])
                dlon = Float32(lons[2, 1] - lons[1, 1])
                # stored reversed in j (ascending), like the main reader
                lat1 = Float32(lats[1, Nj]); dlat = Float32(lats[1, Nj-1] - lats[1, Nj])
            end
            lv = Float64(m["level"]); lv in levs || push!(levs, lv)
        end
    end
    sort!(levs)
    pt = PTWinds(zeros(Float32, Ni, Nj, length(levs), 2),
                 zeros(Float32, Ni, Nj, length(levs), 2),
                 Float32.(levs), lon0, dlon, Ni, lat1, dlat, Nj,
                 dates, files, 0, 0)
    return pt
end

function _load_pt_file!(pt::PTWinds, fname::AbstractString, slot::Int)
    lmap = Dict{Float64,Int}(Float64(v) => k for (k, v) in enumerate(pt.thetalev))
    GribFile(fname) do g
        for msg in g
            m = Message(msg)
            nm = m["shortName"]
            (nm == "u" || nm == "v") || continue
            k = get(lmap, Float64(m["level"]), 0)
            k == 0 && continue
            dst = nm == "u" ? pt.u : pt.v
            dst[:, :, k, slot] = reverse(m["values"], dims=2)
        end
    end
    return nothing
end

# Ensure the loaded 3-hourly bracket covers `date` (cheap full reloads).
function read_pt_step!(pt::PTWinds, date::DateTime)
    n  = length(pt.filedates)
    i1 = something(findlast(d -> d <= date, pt.filedates), 1)
    i2 = min(i1 + 1, n)
    if (pt.t1, pt.t2) != (i1, i2)
        _load_pt_file!(pt, pt.filenames[i1], 1)
        i2 != i1 && _load_pt_file!(pt, pt.filenames[i2], 2)
        i2 == i1 && (pt.u[:, :, :, 2] = pt.u[:, :, :, 1]; pt.v[:, :, :, 2] = pt.v[:, :, :, 1])
        pt.t1, pt.t2 = i1, i2
    end
    return nothing
end

# u, v on the isentrope θ at (lat, lon, date) — quadrilinear, vertical in θ.
@inline function _pt_uv(pt::PTWinds, lat::Float32, lon::Float32, theta::Float32,
                        date::DateTime)
    nl = length(pt.thetalev)
    (theta < pt.thetalev[1] || theta > pt.thetalev[nl]) && return (NaN32, NaN32)
    ix = mod(lon - pt.lon0, 360.0f0) / pt.dlon + 1.0f0
    iy = (lat - pt.lat0) / pt.dlat + 1.0f0
    (1.0f0 <= ix < Float32(pt.Ni) && 1.0f0 <= iy < Float32(pt.Nj)) || return (NaN32, NaN32)
    k = 1
    @inbounds while k < nl - 1 && pt.thetalev[k+1] < theta
        k += 1
    end
    fk = (theta - pt.thetalev[k]) / (pt.thetalev[k+1] - pt.thetalev[k])
    d1 = pt.filedates[pt.t1]; d2 = pt.filedates[pt.t2]
    ft = d2 == d1 ? 0.0f0 :
         clamp(Float32(Dates.value(date - d1)) / Float32(Dates.value(d2 - d1)), 0.0f0, 1.0f0)
    i = Base.unsafe_trunc(Int, ix); fi = ix - i
    j = Base.unsafe_trunc(Int, iy); fj = iy - j
    @inline tri(A, t) = @inbounds begin
        c00 = (1-fi)*A[i,j,  k,  t] + fi*A[i+1,j,  k,  t]
        c10 = (1-fi)*A[i,j+1,k,  t] + fi*A[i+1,j+1,k,  t]
        c01 = (1-fi)*A[i,j,  k+1,t] + fi*A[i+1,j,  k+1,t]
        c11 = (1-fi)*A[i,j+1,k+1,t] + fi*A[i+1,j+1,k+1,t]
        (1-fk)*((1-fj)*c00 + fj*c10) + fk*((1-fj)*c01 + fj*c11)
    end
    u = (1-ft)*tri(pt.u, 1) + ft*tri(pt.u, 2)
    v = (1-ft)*tri(pt.v, 1) + ft*tri(pt.v, 2)
    return (u, v)
end

mutable struct HybridVert
    fields::Vector{Union{IsentropicFields,Nothing}}  # per nest
    pt::Union{PTWinds,Nothing}  # isentropic-level winds (nest 1)
    theta_p::Vector{Float32}    # prognostic parcel θ [K]; 0 = unset
    sigma_r::Float32            # blend bottom: pure kinematic below (σ > σ_r)
    sigma_t::Float32            # blend top: pure isentropic above (σ < σ_t)
    diabatic::Bool              # false → θ̇ forced to 0 (validation runs)
end

"""
    HybridVert(n_max, rt; sigma_r=0.3, sigma_t=0.1, diabatic=true,
               ptdir=met directory, ptpattern="IFS_PT_MET")

Build the hybrid-vertical handler from the runtime: collects the diabatic
tendency fields of every nest, and the isentropic-level wind files (PT) from
the met directory when present (without them the horizontal motion falls back
to model-level winds sampled on the isentrope).
"""
function HybridVert(n_max::Integer, rt; sigma_r::Real=0.3, sigma_t::Real=0.1,
                    diabatic::Bool=true,
                    ptdir::AbstractString=rt.input[1]["pathmetfields"],
                    ptpattern::AbstractString="IFS_PT_MET")
    flds = Union{IsentropicFields,Nothing}[IsentropicFields(rt.metfields3D[n])
                                           for n in 1:rt.nest]
    pt = try
        PTWinds(ptdir; pattern=ptpattern)
    catch err
        @warn "PT (isentropic-level wind) files unusable; falling back to model-level winds" err
        nothing
    end
    pt === nothing && @info "no $(ptpattern)* files in $ptdir — isentropic horizontal motion uses model-level winds"
    HybridVert(flds, pt, zeros(Float32, Int(n_max)),
               Float32(sigma_r), Float32(sigma_t), diabatic)
end

# Smooth blend weight from σ = p/p_s: 0 in the troposphere, 1 in the strat.
@inline function _hyb_weight(sigma::Float32, sr::Float32, st::Float32)
    sigma >= sr && return 0.0f0
    sigma <= st && return 1.0f0
    x = (sr - sigma) / (sr - st)
    return 0.5f0 * (1.0f0 - cos(Float32(π) * x))
end

# Invert the local θ ladder: θ → continuous level index.  Restricted to the
# index range above kmin (θ monotonic in the upper troposphere/stratosphere).
@inline function _theta_to_indl(th::AbstractVector{Float32}, tq::Float32,
                                kmin::Int, nz::Int)
    tq <= th[kmin] && return Float32(kmin)
    @inbounds for k in (kmin+1):nz
        if tq <= th[k]
            return Float32(k - 1) + (tq - th[k-1]) / max(th[k] - th[k-1], 1.0f-3)
        end
    end
    return Float32(nz - 1)
end

"""
    hybrid_vertical!(hyb, state, ctx, dt_s, t_array)

Apply the hybrid kinematic/isentropic vertical correction to all active
parcels (call after advection).  Parcels in the blend/stratospheric region
are pulled onto their prognostic isentrope θ_p (+ diabatic drift); θ_p is
(re)initialised from the grid for new parcels and parcels in the kinematic
region.
"""
function hybrid_vertical!(hyb::HybridVert, state, ctx, dt_s::Real, t_array)
    dt = Float32(dt_s) * Float32(Int(ctx.runconfig.direction))  # signed: backward runs reverse θ̇
    currentdate = ctx.runconfig.currentdate
    hyb.pt !== nothing && read_pt_step!(hyb.pt, currentdate)
    for n in 1:length(ctx.nests)
        isf = hyb.fields[n]
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8,
                      eachindex(state.status))
        isempty(idx) && continue

        m3   = nestctx.metfields3D
        nz   = size(m3.zlev, 3)
        Ni   = Int(nestctx.gridconfig.Ni)
        Nj   = Int(nestctx.gridconfig.Nj)
        indt = t_array[n]

        nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
        thbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]

        Threads.@threads :static for i in idx
            tid  = Threads.threadid()
            indx = state.indx[i]
            indy = state.indy[i]
            (indx < 1.0f0 || indx >= Float32(Ni) ||
             indy < 1.0f0 || indy >= Float32(Nj)) && continue

            indl = clamp(state.indl[i], 1.0f0, Float32(nz - 1))
            p  = interpolation4Dsingle(m3.puv, indx, indy, indl, indt)
            ps = interpolation4Dsingle(m3.puv, indx, indy, 1.0f0, indt)
            (isfinite(p) && isfinite(ps) && ps > 1.0f4) || continue
            w = _hyb_weight(p / ps, hyb.sigma_r, hyb.sigma_t)

            th = thbufs[tid]
            get_vertical_profile!(th, m3.theta, indx, indy, indt, nz)

            if w <= 0.0f0
                # Pure kinematic region: keep θ_p synchronised with the grid
                # so the parcel enters the blend zone with the right value.
                k0 = clamp(Base.unsafe_trunc(Int, indl), 1, nz - 1)
                hyb.theta_p[i] = th[k0] + (indl - k0) * (th[k0+1] - th[k0])
                continue
            end

            # Initialise θ_p on first contact with the blend zone.
            k0 = clamp(Base.unsafe_trunc(Int, indl), 1, nz - 1)
            th_grid = th[k0] + (indl - k0) * (th[k0+1] - th[k0])
            hyb.theta_p[i] <= 0.0f0 && (hyb.theta_p[i] = th_grid)

            # Diabatic drift θ̇ = (θ/T)·(dT/dt)_diab
            thdot = 0.0f0
            if hyb.diabatic && isf !== nothing
                T  = interpolation4Dsingle(m3.T, indx, indy, indl, indt)
                dT = interpolation4Dsingle(isf.tt1, indx, indy, indl, indt)
                isf.tt2 !== nothing &&
                    (dT += interpolation4Dsingle(isf.tt2, indx, indy, indl, indt))
                if isfinite(T) && isfinite(dT) && T > 150.0f0
                    thdot = hyb.theta_p[i] / T * dT
                end
            end
            th_tgt = hyb.theta_p[i] + thdot * dt

            # Horizontal correction with the ISENTROPIC-LEVEL winds: the
            # model-level advection step interpolated u,v horizontally at
            # constant level, mixing winds across isentropes wherever θ has a
            # gradient on the level.  Replace the fraction w of the
            # displacement with PT-file winds sampled on the parcel's θ.
            if hyb.pt !== nothing
                lat, lon = indyindxtolatlon(nestctx.gridconfig, indy, indx)
                uth, vth = _pt_uv(hyb.pt, Float32(lat), Float32(lon),
                                  hyb.theta_p[i], currentdate)
                if isfinite(uth) && isfinite(vth)
                    uml = interpolateU(nestctx.gridconfig.Uedge, m3.Ux, indx, indy, indl, indt)
                    vml = interpolateV(nestctx.gridconfig.Vedge, m3.Vy, indx, indy, indl, indt)
                    if isfinite(uml) && isfinite(vml)
                        dxc = w * (uth - uml) * dt
                        dyc = w * (vth - vml) * dt
                        state.indx[i] += dxtodindx(nestctx.gridconfig, dxc, indx, indy)
                        state.indy[i] += dytodindy(nestctx.gridconfig, dyc, indx, indy)
                    end
                end
            end

            # Invert the ladder over the full column.  θ is strictly monotonic
            # above the boundary layer, and this branch only runs at w > 0
            # (σ < σ_r), so the search range is safe.  A bracket tied to the
            # post-kinematic k0 must NOT be used: a single large etadot-noise
            # excursion would clamp the inversion at the bracket edge and the
            # parcel would ratchet off its isentrope.
            indl_th = _theta_to_indl(th, th_tgt, 2, nz)

            indl_new = indl + w * (indl_th - indl)
            state.indl[i] = clamp(indl_new, 1.0f0, Float32(nz - 1))

            # Prognostic θ update: exact in the isentropic limit (no grid
            # re-interpolation when w = 1), grid-relaxed inside the blend.
            k1 = clamp(Base.unsafe_trunc(Int, state.indl[i]), 1, nz - 1)
            th_new = th[k1] + (state.indl[i] - k1) * (th[k1+1] - th[k1])
            hyb.theta_p[i] = w * th_tgt + (1.0f0 - w) * th_new
        end
    end
    return nothing
end
