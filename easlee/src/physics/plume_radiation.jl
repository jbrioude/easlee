# physics/plume_radiation.jl
# ─────────────────────────────────────────────────────────────────────────────
# In-plume longwave radiative cooling from EXCESS water vapour
# ([physics.plume_radiation]) — the observed few-K cold anomaly inside the
# Hunga Tonga water plume (radiosonde analyses; radiative-transfer estimates
# give heating rates of a few tenths of K/day per ~100 ppmv excess) and the
# slow radiative DESCENT of the layer.
#
# Closure: a LINEARISED cooling-to-space term per parcel,
#
#     dT/dt = −c_rad · max(q_ppmv − q_bg, 0)
#
# with q from the parcel's ice-closure air mass (q = mass/airm — the same
# per-parcel mixing ratio the freezing module uses).  c_rad is an EXPLICIT
# calibratable coefficient (configured in K/day per ppmv); this is not a
# radiation code — validate it against published heating rates and the sonde
# anomalies before citing.
#
# Two effects:
#   • diagnostic: each parcel accumulates its ΔT_rad along the trajectory
#     (`trad`) — mass-weighted stats at run end; compare to the in-plume
#     sonde temperature anomalies;
#   • optional feedback (`feedback = true`): the cooled parcel sinks to its
#     new buoyancy-equilibrium level through the ambient θ profile,
#     dz = dθ_parcel/(∂θ/∂z) with dθ = dT·(p0/p)^κ — the radiative descent
#     of the water layer.
# ─────────────────────────────────────────────────────────────────────────────

const _PRAD_MW = 28.9647 / 18.0153

struct PlumeRadiation
    on::Vector{Bool}          # release boxes whose WATER produces the cooling
    partner::Vector{Bool}     # boxes that FEEL the co-located cooling without
                              # producing it (e.g. the sulfate aerosol riding the
                              # moist plume): they get the met-cell mean excess-q
                              # cooling of the water parcels around them + the
                              # same buoyancy descent.  As the aerosol settles
                              # out of the moist layer the cell excess → 0 and
                              # its radiative descent stops — the aerosol/water
                              # SEPARATION emerges from the co-location itself.
    c_rad::Float64            # K/s per ppmv excess (OPTICALLY THIN limit)
    # ── LW OPACITY ([physics.plume_radiation].lw_kappa) ──────────────────────
    # The linear c_rad·excess is the OPTICALLY THIN limit, where the layer
    # thickness cancels: eps ~ kappa·q·h and the layer mass ~ rho·h, so
    # dT/dt ~ c_rad·q independent of h.  The HT plume is NOT in that limit —
    # at 1000-3000 ppmv over 1-3 km the emissivity saturates (eps -> 1), and a
    # saturated layer cools at its EDGES, so the rate falls as 1/h.
    # With lw_kappa > 0 [1/(ppmv·m)]:
    #     dT/dt = -(c_rad/(kappa·h)) · (1 - exp(-kappa·h·q_excess))
    # which reduces EXACTLY to -c_rad·q_excess when kappa·h·q << 1, so
    # lw_kappa = 0 keeps the legacy behaviour (hard cap) bit-for-bit.
    # h is the per-parcel layer thickness from update_airm! (ice.hv) — so both
    # the concentration AND the thickness come from the trajectories.
    lw_kappa::Float64
    q_bg::Float64             # background ppmv subtracted before cooling
    feedback::Bool            # buoyancy-equilibrium descent
    trad::Vector{Float32}     # per-parcel accumulated ΔT_rad [K] (slot-indexed;
                              # zeroed for inactive slots each call)
    hv::Vector{Float32}       # per-parcel layer thickness [m] from update_airm!
                              # (empty = fall back to h_ref)
    h_ref::Float64            # fallback thickness [m] when hv is unavailable
    # ── PROGNOSTIC COLD ANOMALY T′ ([physics.plume_radiation].tprime_tau) ────
    # With the hybrid vertical the parcel is in buoyancy equilibrium BY
    # CONSTRUCTION (indl is re-derived from θ_p, so ambient θ at the parcel's
    # level ≡ θ_p): radiative cooling becomes descent instantly and the parcel
    # is never colder than its surroundings.  The observed plume is: sondes
    # show a −4/−5 K in-plume anomaly, and at 220 K (d ln e_si/dT ≈ 0.13/K)
    # those 5 K are a factor ~1.9 on RHi — the difference between 2.6 K and
    # 7–8 K above the frost point (HT_ICE_FINDINGS §5).  A model without this
    # state can never saturate for the right reason.
    #
    # Closure: T′ is the DISEQUILIBRIUM reservoir.  Cooling dT accumulates in
    # T′ (the parcel gets colder than ambient AT ITS LEVEL — θ_p untouched);
    # T′ relaxes on tprime_tau by driving the buoyancy descent (the parcel
    # sinks, adiabatic recompression supplies the heat).  Steady state
    # T′ ≈ −Q_rad·τ; total descent over long times equals the legacy scheme's,
    # only the phase lag and the standing anomaly differ.  τ is the mesoscale
    # buoyancy-adjustment time of a wide thin layer, ~(L/h)/N ≈ 0.5–1 day for
    # L~1000 km, h~1 km, N~0.02 s⁻¹ — an estimate, not a tuning knob; validate
    # the resulting T′ against the in-plume sonde−ERA5 residual.
    # The saturation gates (apply_convert_thermo!, apply_ice!) consume `tanom`
    # as e_si(T_ambient + T′) at ambient p.  tprime_tau = 0 → legacy path
    # bit-for-bit (cooling straight to descent, no anomaly state).
    tprime_tau::Float64       # relaxation time [s]; 0 = off (legacy)
    tanom::Vector{Float32}    # per-parcel T′ [K] ≤ 0 (slot-indexed; empty = off)
end

# |T′| safety cap [K]: ~3× the observed −4/−5 K anomaly.  Guards against a
# runaway if q or the cooling rate misbehaves; a real anomaly beyond this
# would collapse by convective overturning anyway.
const _TANOM_CAP = 15.0f0

# legacy 5-arg form (no partners)
PlumeRadiation(on::Vector{Bool}, c_rad::Real, q_bg::Real, feedback::Bool,
               trad::Vector{Float32}) =
    PlumeRadiation(on, falses(length(on)), Float64(c_rad), 0.0, Float64(q_bg),
                   feedback, trad, Float32[], 1000.0, 0.0, Float32[])
# 6-arg form used by build_run (partners, no LW opacity)
PlumeRadiation(on::Vector{Bool}, partner::Vector{Bool}, c_rad::Real, q_bg::Real,
               feedback::Bool, trad::Vector{Float32}) =
    PlumeRadiation(on, partner, Float64(c_rad), 0.0, Float64(q_bg),
                   feedback, trad, Float32[], 1000.0, 0.0, Float32[])
# 9-arg form (LW opacity + thickness, no T′ state)
PlumeRadiation(on::Vector{Bool}, partner::Vector{Bool}, c_rad::Real, lw_kappa::Real,
               q_bg::Real, feedback::Bool, trad::Vector{Float32},
               hv::Vector{Float32}, h_ref::Real) =
    PlumeRadiation(on, partner, Float64(c_rad), Float64(lw_kappa), Float64(q_bg),
                   feedback, trad, hv, Float64(h_ref), 0.0, Float32[])

"""
Radiative cooling rate [K/s] for excess `ex` ppmv in a layer of thickness `h` m.
kappa = 0 -> legacy linear law with the hard cap.
"""
@inline function _prad_rate(pr::PlumeRadiation, ex::Float64, h::Float64)
    if pr.lw_kappa <= 0.0
        return min(pr.c_rad * ex, _PRAD_CAP_K_DAY / 86400.0)
    end
    τ = pr.lw_kappa * h * ex
    return (pr.c_rad / (pr.lw_kappa * h)) * (τ > 1e-4 ? -expm1(-τ) : τ)
end

any_prad(pr::PlumeRadiation) = any(pr.on) || any(pr.partner)

"""
    apply_plume_radiation!(state, ctx, pr, dt_s, t_array) → (n, ΣmΔT)

One radiative-cooling step for the parcels of the enabled boxes (needs the
ice-closure `airm`; parcels without it are skipped).  Returns the touched
count and the mass-weighted ΔT sum of this step (for budget logging).
"""
# DEBUG counters — why the descent no-ops (printed in the run-end budget log)
const PRAD_DBG = (dthdz0 = Threads.Atomic{Int}(0), plow = Threads.Atomic{Int}(0),
                  nonfin = Threads.Atomic{Int}(0), ok = Threads.Atomic{Int}(0),
                  dzsum = Threads.Atomic{Float64}(0.0))

# Saturation of the cooling-to-space closure [K/day]: the linear c_rad·excess
# diverges in the near-source plume (q ~ 10⁴-10⁶ ppmv, and the ice-closure q
# never dilutes) while a real optically-thick layer cools at a few K/day at
# most (Sellitto 2022; observed layer subsidence ≈ −2 to −3 K/day in θ·T/θ).
const _PRAD_CAP_K_DAY = 5.0

# buoyancy-equilibrium descent of parcel i after a radiative ΔT (shared by the
# water pass and the partner pass).  With the HYBRID/ISENTROPIC vertical mode
# the parcel's z is re-imposed from its prognostic θ_p every step — moving indl
# alone is UNDONE next step.  The radiative ΔT is a genuine diabatic θ̇, so when
# `thetap` (HybridVert.theta_p) is provided the descent is applied THERE
# (dθ = dT·(p₀/p)^κ); the indl move is kept for the pure-kinematic case.
@inline function _prad_descend!(state, m3, i::Int, dT::Float64, nz::Int,
                                indt::Float32, z::Vector{Float32},
                                thetap::Union{Nothing,Vector{Float32}} = nothing)
    indx = state.indx[i]; indy = state.indy[i]
    indl = clamp(state.indl[i], 1.0f0, Float32(nz - 1))
    θ1 = interpolation4Dsingle(m3.theta, indx, indy,
                               max(indl - 0.5f0, 1.0f0), indt)
    θ2 = interpolation4Dsingle(m3.theta, indx, indy,
                               min(indl + 0.5f0, Float32(nz)), indt)
    get_vertical_profile!(z, m3.zlev, indx, indy, indt, nz)
    z1 = _z_at(z, max(indl - 0.5f0, 1.0f0), nz)
    z2 = _z_at(z, min(indl + 0.5f0, Float32(nz)), nz)
    dθdz = (Float64(θ2) - Float64(θ1)) / max(abs(Float64(z2) - Float64(z1)), 1.0)
    abs(dθdz) > 1e-5 || (Threads.atomic_add!(PRAD_DBG.dthdz0, 1); return nothing)
    p = Float64(interpolation4Dsingle(m3.puv, indx, indy, indl, indt))
    p > 10.0 || (Threads.atomic_add!(PRAD_DBG.plow, 1); return nothing)
    dθ = dT * (1.0e5 / p)^0.2857
    dz = dθ / dθdz                                         # < 0: descent
    znew = _z_at(z, indl, nz) + Float32(dz)
    isfinite(znew) || (Threads.atomic_add!(PRAD_DBG.nonfin, 1); return nothing)
    state.indl[i] = _z_to_indl(z, max(znew, 0.1f0), nz)
    # hybrid/isentropic mode: the diabatic increment must live on θ_p, or the
    # next hybrid step re-imposes the old isentrope and undoes the descent
    if thetap !== nothing && thetap[i] > 0.0f0
        thetap[i] += Float32(dθ)
    end
    Threads.atomic_add!(PRAD_DBG.ok, 1)
    Threads.atomic_add!(PRAD_DBG.dzsum, dz)
    return nothing
end

function apply_plume_radiation!(state::State.ParcelState, ctx, pr::PlumeRadiation,
                                dt_s::Real, t_array; nest::Int = 1,
                                thetap::Union{Nothing,Vector{Float32}} = nothing)
    nb = length(pr.on)
    nestctx = ctx.nests[nest]
    m3 = nestctx.metfields3D
    nz = Int(nestctx.levelconfig.Nk)
    indt = Float32(t_array[nest])
    n8 = eltype(state.nest)(nest)
    ntouch = Threads.Atomic{Int}(0)
    mdt = Threads.Atomic{Float64}(0.0)
    nt = Threads.maxthreadid()
    zbufs = [zeros(Float32, nz) for _ in 1:nt]
    haspartner = any(pr.partner)
    # T′ reservoir active?  relfrac = the fraction of the anomaly relaxed (into
    # descent) this step; the remainder stands as the parcel's cold anomaly.
    hastp = pr.tprime_tau > 0.0 && !isempty(pr.tanom)
    relfrac = hastp ? -expm1(-Float64(dt_s) / pr.tprime_tau) : 0.0
    # DIRECTION (fixed 2026-08-07).  The caller passes |dt|, so without this the
    # plume cooled — and therefore SANK — even in a backward run, which is the
    # wrong sign: a layer that descends radiatively forward in time was HIGHER
    # in the past.  Only the TENDENCY flips; `relfrac` above is a relaxation
    # timescale, not a tendency, and stays positive.
    #
    # NOTE: the T′ reservoir is a MEMORY process and is NOT time-reversible —
    # charging it backward is not the inverse of relaxing it forward.  For
    # backward runs use `tprime_tau = ""` (cooling goes straight into the
    # descent, which is what the ~7-day integrated displacement depends on).
    sgn = Float64(Int(ctx.runconfig.direction))
    hastp && sgn < 0 && @warn "plume_radiation: T′ reservoir is not time-reversible in a backward run — set tprime_tau=\"\" for trajectory work" tprime_tau=pr.tprime_tau
    # per-thread met-cell accumulators for the partner pass: Σ(m·ex) and Σm of
    # the WATER parcels per cell (same cell key as the ice co-scavenging)
    exd = [Dict{Int64,Float64}() for _ in 1:nt]
    wtd = [Dict{Int64,Float64}() for _ in 1:nt]
    Threads.@threads :static for i in eachindex(state.status)
        if state.status[i] != 0x01
            pr.trad[i] = 0.0f0                 # slot hygiene: reset freed slots
            hastp && (pr.tanom[i] = 0.0f0)
            continue
        end
        state.nest[i] == n8 || continue
        # ── T′ relaxation — EVERY active parcel carrying an anomaly, whatever
        # its box (a parcel that converted or left the plume still
        # re-equilibrates): the relaxed portion drives the buoyancy descent,
        # the rest stands as the cold anomaly the saturation gates see.
        if hastp && pr.tanom[i] < -1.0f-6
            rel = Float64(pr.tanom[i]) * relfrac           # K, < 0
            pr.tanom[i] -= Float32(rel)
            pr.feedback &&
                _prad_descend!(state, m3, i, rel, nz, indt,
                               zbufs[Threads.threadid()], thetap)
        end
        b = Int(state.relbox[i])
        (1 <= b <= nb && pr.on[b]) || continue
        airm = Float64(state.airm[i])
        airm > 0.0 || continue
        q  = Float64(state.mass[i]) / airm * _PRAD_MW * 1e6    # ppmv
        ex = q - pr.q_bg
        ex > 0.0 || continue
        if haspartner
            tid = Threads.threadid()
            _, _, _, _, ii, jj, kk, _ = _parcel_indices(
                state.indx[i], state.indy[i], state.indl[i], indt)
            key = _ice_cellkey(ii, jj, kk)
            m = Float64(state.mass[i])
            exd[tid][key] = get(exd[tid], key, 0.0) + m * ex
            wtd[tid][key] = get(wtd[tid], key, 0.0) + m
        end
        hpar = (!isempty(pr.hv) && pr.hv[i] > 0f0) ? Float64(pr.hv[i]) : pr.h_ref
        dT = -_prad_rate(pr, ex, hpar) * Float64(dt_s) * sgn   # K (< 0 forward, > 0 backward)
        pr.trad[i] += Float32(dT)
        Threads.atomic_add!(ntouch, 1)
        Threads.atomic_add!(mdt, Float64(state.mass[i]) * dT)
        if hastp
            # cooling charges the anomaly reservoir; descent comes from its
            # relaxation above, NOT from dT directly (no double counting)
            pr.tanom[i] = max(pr.tanom[i] + Float32(dT), -_TANOM_CAP)
        else
            pr.feedback &&
                _prad_descend!(state, m3, i, dT, nz, indt, zbufs[Threads.threadid()], thetap)
        end
    end
    # ── partner pass: boxes riding the moist plume feel the CELL's mean excess
    # (mass-weighted over the co-located water parcels).  No airm needed — the
    # cooling is a property of the air the partner sits in, not of its own mass.
    if haspartner
        exm = exd[1]
        for t in 2:nt, (k, v) in exd[t]; exm[k] = get(exm, k, 0.0) + v; end
        wtm = wtd[1]
        for t in 2:nt, (k, v) in wtd[t]; wtm[k] = get(wtm, k, 0.0) + v; end
        if !isempty(exm)
            Threads.@threads :static for i in eachindex(state.status)
                state.status[i] == 0x01 || continue
                state.nest[i] == n8     || continue
                b = Int(state.relbox[i])
                (1 <= b <= nb && pr.partner[b]) || continue
                _, _, _, _, ii, jj, kk, _ = _parcel_indices(
                    state.indx[i], state.indy[i], state.indl[i], indt)
                key = _ice_cellkey(ii, jj, kk)
                w = get(wtm, key, 0.0)
                w > 0.0 || continue
                ex = get(exm, key, 0.0) / w                    # cell-mean excess ppmv
                ex > 0.0 || continue
                hpar = (!isempty(pr.hv) && pr.hv[i] > 0f0) ? Float64(pr.hv[i]) : pr.h_ref
                dT = -_prad_rate(pr, ex, hpar) * Float64(dt_s) * sgn
                pr.trad[i] += Float32(dT)
                Threads.atomic_add!(ntouch, 1)
                if hastp
                    # partner (e.g. the ICE bins) charges its own T′ — the
                    # relaxation in the main loop already handled this parcel
                    pr.tanom[i] = max(pr.tanom[i] + Float32(dT), -_TANOM_CAP)
                else
                    pr.feedback &&
                        _prad_descend!(state, m3, i, dT, nz, indt, zbufs[Threads.threadid()], thetap)
                end
            end
        end
    end
    return ntouch[], mdt[]
end
