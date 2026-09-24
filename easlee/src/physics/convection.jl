
# ── Convective redistribution from archived IFS mass fluxes ───────────────────
#
# Monte-Carlo redistribution of parcels in convective columns, driven by the
# ERA5/IFS archived convection diagnostics (CLD companion files):
#
#   umf  M_u(z) ≥ 0   updraught mass flux            [kg m⁻² s⁻¹]
#   dmf  M_d(z) ≤ 0   downdraught mass flux          [kg m⁻² s⁻¹]
#   udr  D_u(z) ≥ 0   updraught detrainment rate     [kg m⁻³ s⁻¹]   (unused — see below)
#   ddr  D_d(z) ≥ 0   downdraught detrainment rate   [kg m⁻³ s⁻¹]   (unused — see below)
#
# WELL-MIXED BY CONSTRUCTION.  The defining test of a convection scheme is that
# an initially uniform (in air mass) tracer stays uniform.  An earlier version
# entrained parcels with E=∂M/∂z+D and let the rest subside at a velocity
# w_env=−(M_u+M_d)/ρ; because the non-local updraught conveyor's discrete flux
# through each interface did not exactly equal the archived M_u, the up- and
# down-transport did not cancel and the sub-cloud drained to χ≈0.5 (a structural,
# step-size-independent bias — confirmed by a Δt = 15/5/1 min sweep all
# plateauing at the same violation).  This scheme removes that bias:
#
# 1. Entrainment and detrainment are taken DIRECTLY from the flux divergence of
#    the interface mass fluxes  M_u^if[k]=½(M_u[k]+M_u[k+1]):
#        entrain_up[k] = max(M_u^if[k] − M_u^if[k−1], 0)·Δt    (out of layer k)
#        detrain_up[k] = max(M_u^if[k−1] − M_u^if[k], 0)·Δt    (into  layer k)
#    so the parcel conveyor's flux through every interface equals the archived
#    M_u exactly (mirror for the downdraught).  The archived udr/ddr are not
#    needed — the flux divergence already is entrainment minus detrainment.
#
# 2. A parcel entrained into the updraught rides the flux and detrains in layer
#    l with probability 1 − M_u^if[l]/M_u^if[l−1] (the flux fraction lost there),
#    detraining unconditionally where the flux dies (cloud top).  This keeps the
#    physical NON-LOCAL transport (sub-cloud air lofted to cloud top in one
#    step).  Downdraught capture is the mirror, riding downward.
#
# 3. Compensating subsidence is the EXACT mass-conserving residual, not a
#    velocity.  With m1[k] the air-mass image of the conveyor and
#    S[k]=Σ_{l≤k}(dm[l]−m1[l]) the cumulative deficit below interface k, every
#    parcel then descends across interface c−1/2 with probability
#    S[c−1]/(dm[c]+S[c−1]).  This composite maps the air-mass vector dm onto
#    itself in expectation ⇒ a uniform tracer stays uniform, while the conveyor
#    above carries the asymmetric transport.  (Verified: BL std(χ) 0.155→0.021.)
#
# All probabilities use the grid-mean, time-mean archived fluxes, which is
# consistent with parcels representing grid-mean air masses.  The scheme is
# applied with dt = |model timestep| in both directions.  BACKWARD runs use the
# ADJOINT of the forward kernel (2026-09-24; until then the forward kernel was
# applied unchanged, with the comment "its own adjoint" — it is not: the conveyor
# is strongly non-symmetric).  For a mass-conserving kernel P(k→l) the adjoint is
# the mass-weighted transpose P*(l→k) = P(k→l)·dm[k]/dm[l]; because the forward
# step is the composition  subsidence ∘ conveyor, the backward step is
# conveyor* ∘ subsidence*, each reversed exactly as a walk over the pooled masses:
#   subsidence*: from l, ascend across l+1/2 with probability S[l]/(dm[l]+S[l−1])
#                (the share of the pool at l that had arrived from above), repeat;
#   conveyor*:   at c, came from the resident air with probability (dm−mo)/m1,
#                from the updraught with up_in/m1, from the downdraught with
#                dn_in/m1; the draught origin is found walking DOWN the updraught
#                (up the downdraught): at layer j the entrained share of the
#                draught flux is (F_{j+1/2} − F_{j−1/2})/F_{j+1/2}.
# test/test_convection_adjoint.jl checks Q(l→k)·dm[l] = P(k→l)·dm[k] by Monte-Carlo
# and that a uniform ensemble stays uniform under both kernels.  The `cape`
# scheme (diagnostic_convection!) mixes a column toward uniform: its kernel is
# symmetric in air mass and needs no transpose.
#
# Vertical indexing: profiles are collocated with zlev (k = 1 surface … Nk top)
# exactly like every other 3-D met field after the GRIB reader's surface-layer
# shift.  The parcel's continuous level index indl lives in the same space.

const _CONV_MF_EPS = 1.0f-5   # [kg m⁻² s⁻¹] below this the draught has died out

"""
    convective_redistribution!(state, ctx, convfields, dt_s, t_array; rng)

Apply one convective redistribution step to all active parcels.
`convfields` is a Vector (one entry per nest) of `ConvFields` or `nothing`;
nests with `nothing` are skipped.  `dt_s` is the step length in seconds.
"""
function convective_redistribution!(state, ctx, convfields, dt_s::Real, t_array;
                                    rng = nothing, conv_count = nothing)
    for n in 1:length(ctx.nests)
        cf = convfields[n]
        cf === nothing && continue
        nestctx = ctx.nests[n]

        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8,
                      eachindex(state.status))
        isempty(idx) && continue

        _convect_nest!(state, nestctx, cf, Float32(dt_s), t_array[n], idx;
                       conv_count = conv_count,
                       backward = Int(ctx.runconfig.direction) == -1)   # adjoint kernel (2026-09-24)
    end
    return nothing
end

function _convect_nest!(state, nestctx, cf, dt::Float32, indt, idx;
                        conv_count = nothing, backward::Bool = false)
    zlev_arr = nestctx.metfields3D.zlev
    rho_arr  = nestctx.metfields3D.rho
    nz       = size(zlev_arr, 3)
    Ni       = Int(nestctx.gridconfig.Ni)
    Nj       = Int(nestctx.gridconfig.Nj)

    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    mkbuf() = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    zbufs, rhobufs        = mkbuf(), mkbuf()
    mubufs, mdbufs        = mkbuf(), mkbuf()      # archived up/down mass flux (levels)
    dmbufs, m1bufs, Sbufs = mkbuf(), mkbuf(), mkbuf()  # layer air mass, post-conveyor image, cum. deficit

    Threads.@threads :static for i in idx
        tid = Threads.threadid()

        indx = state.indx[i]
        indy = state.indy[i]
        (indx < 1.0f0 || indx >= Float32(Ni) ||
         indy < 1.0f0 || indy >= Float32(Nj)) && continue

        mu = mubufs[tid]
        get_vertical_profile!(mu, cf.umf, indx, indy, indt, nz)
        md = mdbufs[tid]
        get_vertical_profile!(md, cf.dmf, indx, indy, indt, nz)

        # Fast exit: no convection in this column.
        mumax = 0.0f0
        mdmax = 0.0f0
        @inbounds for k in 1:nz
            mu[k] = max(mu[k], 0.0f0)
            md[k] = max(-md[k], 0.0f0)       # md now holds the MAGNITUDE m_d ≥ 0
            mu[k] > mumax && (mumax = mu[k])
            md[k] > mdmax && (mdmax = md[k])
        end
        (mumax < _CONV_MF_EPS && mdmax < _CONV_MF_EPS) && continue

        z   = zbufs[tid]
        get_vertical_profile!(z,   zlev_arr, indx, indy, indt, nz)
        rho = rhobufs[tid]
        get_vertical_profile!(rho, rho_arr,  indx, indy, indt, nz)

        # ── Build the column's mass-conserving redistribution from the archived
        #    fluxes.  Entrainment/detrainment are taken DIRECTLY from the flux
        #    divergence of the interface mass fluxes (Mu_if, Md_if), so the
        #    parcel conveyor's flux through every interface equals the archived
        #    flux exactly — the prerequisite for the up- and down-transport to
        #    cancel for a uniform tracer.  Fills dm (layer air mass), m1 (the
        #    air-mass image after the conveyor) and S (cumulative deficit below
        #    each level) that drive the compensating-subsidence descent.
        dm = dmbufs[tid]; m1 = m1bufs[tid]; S = Sbufs[tid]
        _build_conv_column!(dm, m1, S, z, rho, mu, md, dt, nz)

        indl = state.indl[i]
        k0   = clamp(Base.unsafe_trunc(Int, indl), 1, nz - 1)
        kc   = _convect_column(_ParcelRNG(state, i), k0, dm, m1, S, mu, md, dt, nz, backward)
        # 3) Land at a uniformly random height inside the destination layer.
        @inbounds z_new = z[kc] + _conv_rand(state, i) * (z[min(kc + 1, nz)] - z[kc])
        # Met fields can carry NaN at masked cells; never write a non-finite
        # level index (same policy as the advection kernel).
        isfinite(z_new) || continue
        z_new < 0.0f0 && (z_new = -z_new)         # reflect at ground
        z_new > z[nz - 1] && (z_new = z[nz - 1])  # keep below model top
        state.indl[i] = _z_to_indl(z, z_new, nz)

        # [tracers]: count a convective displacement (destination layer ≠ origin).
        # Per-parcel write under @threads :static (unique i) → no data race.
        (conv_count !== nothing && kc != k0) && (@inbounds conv_count[i] += Int32(1))
    end
    return nothing
end

# ── Column kernel (forward and adjoint) ─────────────────────────────────────
# Per-parcel random numbers from the parcel's own xorshift stream, as a callable.
struct _ParcelRNG{S}
    state::S
    i::Int
end
@inline (r::_ParcelRNG)() = _conv_rand(r.state, r.i)

"""
    _convect_column(rng, k0, dm, m1, S, mu, md, dt, nz, backward) → destination layer

One convective redistribution of a parcel in layer `k0` on a column described by
the layer air masses `dm`, the post-conveyor image `m1`, the cumulative deficit `S`
and the interface up/down mass fluxes (`mu`, `md`, level-collocated, magnitudes).
`backward=false`: conveyor then compensating subsidence (the forward kernel);
`backward=true`: the mass-weighted transpose — reversed subsidence then reversed
conveyor.  `rng()` returns a uniform number in [0,1).
"""
@inline function _convect_column(rng, k0::Int, dm, m1, S, mu, md, dt::Float32, nz::Int, backward::Bool)
    if !backward
        # 1) Conveyor: entrain into the up/downdraught (probabilities = the flux
        #    fraction entrained out of this layer) and walk to the detrainment
        #    level; otherwise the parcel stays (it is "environment" air).
        @inbounds eu_out = dm[k0] > 0f0 ? max(_uflux(mu, k0) - _uflux(mu, k0 - 1), 0f0) * dt / dm[k0] : 0f0
        @inbounds ed_out = dm[k0] > 0f0 ? max(_dflux(md, k0 - 1) - _dflux(md, k0), 0f0) * dt / dm[k0] : 0f0
        tot = eu_out + ed_out
        if tot > 1.0f0; eu_out /= tot; ed_out /= tot; end
        r = rng()
        kc = if r < eu_out
            _flux_walk(rng, mu, k0, nz, +1)        # ascend in updraught
        elseif r < eu_out + ed_out
            _flux_walk(rng, md, k0, nz, -1)        # descend in downdraught
        else
            k0                                     # environment (resident)
        end
        # 2) Compensating subsidence: descend across interfaces to fill the
        #    exact cumulative air-mass deficit S left by the conveyor.  The
        #    descent probability at interface c−1/2 is S[c-1]/(dm[c]+S[c-1]),
        #    which maps the post-conveyor image m1 back onto dm in expectation,
        #    so an initially uniform (in air mass) ensemble stays uniform.
        @inbounds while kc > 1
            denom = dm[kc] + S[kc - 1]
            (denom > 0f0 && rng() < S[kc - 1] / denom) || break
            kc -= 1
        end
        return kc
    else
        # 1*) Reversed subsidence: the pool at level c after the forward descent
        #     holds dm[c] + S[c-1] of which S[c] arrived from above — ascend with
        #     that share, and repeat at the level above.
        kc = k0
        @inbounds while kc < nz - 1
            denom = dm[kc] + (kc > 1 ? S[kc - 1] : 0f0)      # the pool at kc after the forward descent
            (denom > 0f0 && rng() < S[kc] / denom) || break
            kc += 1
        end
        # 2*) Reversed conveyor: the post-conveyor mass m1[c] is resident air
        #     (dm − mo), updraught detrainment (up_in) or downdraught detrainment
        #     (dn_in); the draught's origin is found by walking back along it.
        @inbounds begin
            up_in  = max(_uflux(mu, kc - 1) - _uflux(mu, kc), 0f0) * dt
            dn_in  = max(_dflux(md, kc) - _dflux(md, kc - 1), 0f0) * dt
            mtot   = max(m1[kc], 1f-30)
            r = rng()
            if r < up_in / mtot
                return _flux_walk_adjoint(rng, mu, kc, nz, +1)
            elseif r < (up_in + dn_in) / mtot
                return _flux_walk_adjoint(rng, md, kc, nz, -1)
            else
                return kc
            end
        end
    end
end

# Reverse walk along a draught: a parcel that detrained at layer `l` from the
# updraught (dir=+1) came from a layer below.  Going down from l−1, the draught
# flux at interface j+1/2 is made of the air entrained at j, share
# (F_{j+1/2} − F_{j−1/2})/F_{j+1/2}, and of air from further below; stop where
# the draught flux vanishes.  Mirror for the downdraught (dir=−1, walking up).
@inline function _flux_walk_adjoint(rng, m::AbstractVector{Float32}, l::Int, nz::Int, dir::Int)
    j = l - dir
    @inbounds while 1 <= j <= nz - 1
        f_here  = dir > 0 ? _uflux(m, j)     : _dflux(m, j - 1)   # flux leaving j toward l
        f_below = dir > 0 ? _uflux(m, j - 1) : _dflux(m, j)       # flux entering j from further away
        if f_here < _CONV_MF_EPS                                  # no draught leaves j: it started at j+dir
            j += dir; break
        end
        pent = clamp((f_here - f_below) / f_here, 0.0f0, 1.0f0)   # entrained-at-j share
        rng() < pent && break
        j -= dir
    end
    return clamp(j, 1, nz - 1)
end

# Interface mass fluxes (collocated halfway between levels k and k+1); 0 at the
# surface (k=0) and model top (k=nz).  Magnitudes ≥ 0.
@inline _uflux(mu::AbstractVector{Float32}, k::Int) =
    (k < 1 || k >= length(mu)) ? 0.0f0 : 0.5f0 * (mu[k] + mu[k + 1])
@inline _dflux(md::AbstractVector{Float32}, k::Int) =
    (k < 1 || k >= length(md)) ? 0.0f0 : 0.5f0 * (md[k] + md[k + 1])

# Build the per-column air-mass redistribution bookkeeping.
#   dm[k] = ρ̄·Δz   layer air mass (k = 1 … nz−1; dm[nz] = 0)
#   m1[k]          air-mass image after the up/downdraught conveyor
#   S[k]  = Σ_{l≤k}(dm[l] − m1[l])   cumulative deficit below interface k
# By flux telescoping ∑(dm−m1)=0, so S[nz-1]≈0; and since the conveyor moves
# air upward, lower layers are in deficit ⇒ S[k] ≥ 0 (a valid downward current).
@inline function _build_conv_column!(dm, m1, S, z, rho, mu, md, dt::Float32, nz::Int)
    @inbounds for k in 1:nz-1
        dz   = max(z[k + 1] - z[k], 1.0f0)
        rmid = 0.5f0 * (max(rho[k], 1.0f-4) + max(rho[k + 1], 1.0f-4))
        dm[k] = rmid * dz
        # net flux convergence → detrainment into k; divergence → entrainment out
        up_in  = max(_uflux(mu, k - 1) - _uflux(mu, k), 0.0f0) * dt   # updraught detrains here
        up_out = max(_uflux(mu, k) - _uflux(mu, k - 1), 0.0f0) * dt   # entrained up
        dn_in  = max(_dflux(md, k) - _dflux(md, k - 1), 0.0f0) * dt   # downdraught detrains here
        dn_out = max(_dflux(md, k - 1) - _dflux(md, k), 0.0f0) * dt   # entrained down
        mo = up_out + dn_out
        mo > dm[k] && (mo *= dm[k] / mo)                 # cannot remove more than the layer holds
        m1[k] = dm[k] - mo + up_in + dn_in
    end
    dm[nz] = 0.0f0; m1[nz] = 0.0f0
    acc = 0.0f0
    @inbounds for k in 1:nz-1
        acc += dm[k] - m1[k]
        S[k] = max(acc, 0.0f0)
    end
    S[nz] = 0.0f0
    return nothing
end

# Walk a parcel entrained into an up- (dir=+1) or downdraught (dir=−1) at layer
# k0, returning the DESTINATION LAYER.  The entrained air joins the draught
# flux and rides one layer at a time; on entering layer l it detrains with
# probability 1 − F_out/F_in, the fraction of the draught flux lost across that
# layer (so the per-parcel walk reproduces the flux-divergence detrainment used
# by `_build_conv_column!`).  It always detrains where the flux dies (cloud top
# / downdraught base) or at the domain edge.
@inline function _flux_walk(rng, m::AbstractVector{Float32},
                            k0::Int, nz::Int, dir::Int)
    l = k0
    @inbounds while true
        l += dir
        if l < 1 || l > nz - 1
            l -= dir                      # rode off the end → detrain at last valid layer
            break
        end
        # flux entering layer l (upstream side) vs leaving it (downstream side)
        f_in  = dir > 0 ? _uflux(m, l - 1) : _dflux(m, l)
        f_out = dir > 0 ? _uflux(m, l)     : _dflux(m, l - 1)
        f_in < _CONV_MF_EPS && break       # draught has died here → detrain
        pdet = 1.0f0 - clamp(f_out / f_in, 0.0f0, 1.0f0)
        rng() < pdet && break
    end
    return clamp(l, 1, nz - 1)
end

# Height of a parcel at continuous level index indl on profile z.
@inline function _z_at(z::AbstractVector{Float32}, indl::Float32, nz::Int)
    k = clamp(Base.unsafe_trunc(Int, indl), 1, nz - 1)
    f = indl - Float32(k)
    return @inbounds z[k] + f * (z[k + 1] - z[k])
end

# Inverse: height → continuous level index (linear search, profiles are short).
@inline function _z_to_indl(z::AbstractVector{Float32}, zq::Float32, nz::Int)
    zq <= z[1] && return 1.0f0
    @inbounds for k in 2:nz
        if zq <= z[k]
            return Float32(k - 1) + (zq - z[k - 1]) / max(z[k] - z[k - 1], 1.0f-3)
        end
    end
    return Float32(nz - 1)
end

# Per-parcel uniform random number drawn from the parcel's xorshift64 stream
# (thread-safe: each parcel is visited by exactly one thread).  Self-seeds
# parcels whose stream was never initialised.
@inline function _conv_rand(state, i::Int)
    s = state.rng_state[i]
    s == UInt64(0) && (s = rng_seed(i))
    s, u = xorshift_uniform01(s)
    state.rng_state[i] = s
    return u
end

# ════════════════════════════════════════════════════════════════════════════
# Diagnostic convection for datasets WITHOUT archived mass fluxes (GFS, …),
# following the MPTRAC "extreme convection parameterization" (Hoffmann et al.):
#
#   1. A convective event is diagnosed where archived CAPE ≥ cape_min
#      (and CIN ≤ cin_max when archived).
#   2. The equilibrium level (EL, level of neutral buoyancy) is found from the
#      θe of a near-surface parcel vs the environment's saturated θe profile
#      (Bolton 1980 formulas; q reconstructed from Tv = p/(R_d·ρ)).
#   3. Parcels below the EL are redistributed UNIFORMLY IN AIR MASS over
#      [0, z_EL] with probability  p = 1 − exp(−Δt/τ_mix)  (τ_mix → 0 gives
#      the MPTRAC "extreme" limit: the convective column mixes completely
#      every time step).
#
# Uniform-in-mass redistribution preserves the well-mixed criterion by
# construction.  Compared with the explicit mass-flux scheme above, this is an
# upper bound on convective vertical transport — use the archived-flux scheme
# whenever umf/detrainment fields are available.

"""
Tunables for the diagnostic (CAPE-gated) convection scheme; see
`diagnostic_convection!`.  Defaults follow the MPTRAC extreme-convection
setup (τ_mix = 0 → complete column mixing each triggered step).
"""
Base.@kwdef struct DiagConvParams
    cape_min::Float32  = 100.0f0    # [J/kg] event threshold
    cin_max::Float32   = 50.0f0     # [J/kg] skip strongly capped columns
    tau_mix::Float32   = 0.0f0      # [s] mixing timescale; 0 = extreme limit
    zel_min::Float32   = 2000.0f0   # [m] ignore events shallower than this
end

"""
    diagnostic_convection!(state, ctx, diagconv, dt_s, t_array; params)

CAPE-gated convective redistribution.  `diagconv` is a Vector (one per nest)
of `DiagConvFields` or `nothing`.
"""
function diagnostic_convection!(state, ctx, diagconv, dt_s::Real, t_array;
                                params::DiagConvParams = DiagConvParams())
    for n in 1:length(ctx.nests)
        df = diagconv[n]
        df === nothing && continue
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8,
                      eachindex(state.status))
        isempty(idx) && continue
        _diag_convect_nest!(state, nestctx, df, Float32(dt_s), t_array[n], idx, params)
    end
    return nothing
end

# Bolton (1980) saturation specific humidity and equivalent potential temperature.
@inline function _qsat_bolton(T::Float32, p::Float32)
    tc = T - 273.15f0
    es = 611.2f0 * exp(17.67f0 * tc / (tc + 243.5f0))
    es = min(es, 0.5f0 * p)
    return 0.622f0 * es / (p - 0.378f0 * es)
end

@inline function _theta_e_bolton(T::Float32, p::Float32, q::Float32)
    q  = clamp(q, 1.0f-8, 0.04f0)
    rh = clamp(q / max(_qsat_bolton(T, p), 1.0f-8), 0.01f0, 1.0f0)
    TL = 1.0f0 / (1.0f0 / (T - 55.0f0) - log(rh) / 2840.0f0) + 55.0f0
    th = T * (100000.0f0 / p)^(0.2854f0 * (1.0f0 - 0.28f0 * q))
    return th * exp((3376.0f0 / TL - 2.54f0) * q * (1.0f0 + 0.81f0 * q))
end

function _diag_convect_nest!(state, nestctx, df, dt::Float32, indt, idx,
                             p::DiagConvParams)
    m3   = nestctx.metfields3D
    nz   = size(m3.zlev, 3)
    Ni   = Int(nestctx.gridconfig.Ni)
    Nj   = Int(nestctx.gridconfig.Nj)
    pmix = p.tau_mix <= 0.0f0 ? 1.0f0 : 1.0f0 - exp(-dt / p.tau_mix)

    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    zbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    rbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    tbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    pbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    cbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]   # cumulative air mass
    qbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]   # archived qv (when present)

    Threads.@threads :static for i in idx
        tid  = Threads.threadid()
        indx = state.indx[i]
        indy = state.indy[i]
        (indx < 1.0f0 || indx >= Float32(Ni) ||
         indy < 1.0f0 || indy >= Float32(Nj)) && continue

        cape = interpolate2Dsingle(df.cape, indx, indy, indt)
        (isfinite(cape) && cape >= p.cape_min) || continue
        cin = interpolate2Dsingle(df.cin, indx, indy, indt)
        (isfinite(cin) && cin > p.cin_max) && continue

        pmix < 1.0f0 && _conv_rand(state, i) >= pmix && continue

        z = zbufs[tid];  get_vertical_profile!(z,  m3.zlev, indx, indy, indt, nz)
        r = rbufs[tid];  get_vertical_profile!(r,  m3.rho,  indx, indy, indt, nz)
        T = tbufs[tid];  get_vertical_profile!(T,  m3.T,    indx, indy, indt, nz)
        pp = pbufs[tid]; get_vertical_profile!(pp, m3.puv,  indx, indy, indt, nz)

        # Source parcel: lowest model level.  q from the ARCHIVED humidity when the
        # reader provides it; the Tv = p/(R_d·ρ) reconstruction is only a fallback —
        # it yields q ≡ 0 on readers that store dry density (ERA5), killing the EL.
        q2 = if df.qv !== nothing
            qb = qbufs[tid]
            get_vertical_profile!(qb, df.qv, indx, indy, indt, nz)
            qb[2]
        else
            (pp[2] / (287.05f0 * max(r[2], 1.0f-4) * max(T[2], 150.0f0)) - 1.0f0) / 0.608f0
        end
        the_p = _theta_e_bolton(T[2], pp[2], max(q2, 0.0f0))

        # Equilibrium level: topmost level (below 17 km) where the lifted
        # parcel is still buoyant w.r.t. the saturated environment.
        kel = 0
        @inbounds for k in nz:-1:3
            (isfinite(z[k]) && z[k] <= 17000.0f0) || continue
            thes = _theta_e_bolton(T[k], max(pp[k], 100.0f0), _qsat_bolton(T[k], max(pp[k], 100.0f0)))
            if the_p > thes
                kel = k
                break
            end
        end
        kel == 0 && continue
        z[kel] < p.zel_min && continue

        # Only parcels below the EL take part in the convective mixing.
        zp = _z_at(z, state.indl[i], nz)
        zp > z[kel] && continue

        # Redistribute uniformly in air mass over [0, z_EL].
        cm = cbufs[tid]
        cm[1] = 0.0f0
        @inbounds for k in 2:kel
            dm = 0.5f0 * (r[k-1] + r[k]) * (z[k] - z[k-1])
            cm[k] = cm[k-1] + (isfinite(dm) ? max(dm, 0.0f0) : 0.0f0)
        end
        cm[kel] <= 0.0f0 && continue
        u = _conv_rand(state, i) * cm[kel]
        kdst = kel - 1
        @inbounds for k in 2:kel
            if cm[k] >= u
                kdst = k - 1
                break
            end
        end
        f = (u - cm[kdst]) / max(cm[kdst+1] - cm[kdst], 1.0f-6)
        z_new = z[kdst] + f * (z[kdst+1] - z[kdst])
        isfinite(z_new) || continue
        state.indl[i] = _z_to_indl(z, clamp(z_new, 0.0f0, z[nz-1]), nz)
    end
    return nothing
end
