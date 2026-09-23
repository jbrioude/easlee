# physics/ice_physics.jl
# ─────────────────────────────────────────────────────────────────────────────
# Transport-phase ice freezing + sedimentation for tagged water species
# ([[species]] [species.ice]) — the first-48-hour stratospheric water loss
# that the eruption-column washout cannot see (it ends at detrainment).
#
# Physics: a freshly injected volcanic water cloud at 19–25 km is massively
# supersaturated with respect to ICE at < 200 K.  Above a homogeneous-
# nucleation threshold RHi_crit the excess vapour condenses and the ice
# sediments out of the layer on the fall time scale h_layer/v_t — permanently
# removed from the stratospheric budget (re-evaporation happens below the
# tropopause, outside the tagged budget).  Fingerprints in the HT
# observations: the satellite footprint collapse at +9 h (ice evaporation)
# and the 200 Tg emergent vs ~150 Tg MLS burden gap.
#
# The Lagrangian closure: a parcel carries tagged MASS, not a mixing ratio.
# Each parcel is assigned the AIR MASS it represents at release,
#     M_air = m₀ / q₀,
# with q₀ the water mass fraction of its detrained plume member
# (plume_rise.mode="parcels": the ensemble member's saturation-adjusted qv at
# detrainment — physically the right number, per parcel) or the prescribed
# [species.ice] q0 otherwise.  The parcel's mixing ratio at any later time is
# then q = m(t)/M_air — consistent with the pairwise-mixing engine, where
# exchanging mass between equal-air-mass parcels IS exchanging mixing ratio.
#
# Saturation over ice: Murphy & Koop (2005), QJRMS 131, Eq. (7) — valid
# 110–273 K.  Terminal velocity: Stokes–Cunningham at the parcel's LOCAL T, p
# (aerosol.jl settling_velocity, ice density 917 kg/m³) at the prescribed
# effective radius — pressure-dependent slip matters at stratospheric p.
# (The former closure v_t = 12·√(r/1 mm) m/s is a LARGE-hydrometeor plume-
# washout law: at r_eff = 20 µm it gave 1.7 m/s ≈ 145 km/day, draining the
# stratospheric layer ~30× too fast.  Stokes–Cunningham at 25 km gives
# ~6 cm/s ≈ 5 km/day for the same radius.)
# ─────────────────────────────────────────────────────────────────────────────

"""
    esat_ice(T_K) → Pa

Saturation vapour pressure over hexagonal ice, Murphy & Koop (2005) Eq. (7).
"""
@inline esat_ice(T::Real) =
    exp(9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T)

# v1.1 review A5: ONE saturation closure for the plume column / plume-phase /
# umbrella adjustments (plume_rise.jl `_esat`, plume_phase.jl `_esat_ice`),
# replacing the Bolton liquid fit (valid −35…+35 °C, 1.7–1.9× e_si at
# 200–215 K).  Murphy & Koop (2005) Eq. (10) over LIQUID water (123–332 K)
# above the homogeneous-freezing temperature, Eq. (7) over ICE below it.
"""
    esat_liq_mk(T_K) → Pa

Saturation vapour pressure over liquid water, Murphy & Koop (2005) Eq. (10).
"""
@inline esat_liq_mk(T::Real) =
    exp(54.842763 - 6763.22 / T - 4.210 * log(T) + 0.000367 * T +
        tanh(0.0415 * (T - 218.8)) *
        (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T))

"Homogeneous-freezing switch [K] between the ice and liquid saturation laws."
const _ESAT_T_ICE = 235.0

"""
    esat_mk(T_K) → Pa

Saturation vapour pressure of the condensate phase a plume parcel holds at
`T`: over ice (Murphy & Koop Eq. 7) below 235 K, over liquid (Eq. 10) above.
"""
@inline esat_mk(T::Real) = T < _ESAT_T_ICE ? esat_ice(T) : esat_liq_mk(T)

"""
    IcePhysics

Per-release-box ice microphysics parameters (built from the [[species]]
[species.ice] blocks, like SpeciesPhysics).  `on[b] == false` boxes cost one
lookup.  `q0[b] > 0` lazily assigns M_air = m/q0 to parcels the plume
ensemble did not initialise (state.airm == 0).
"""
struct IcePhysics
    on::Vector{Bool}
    rhi_crit::Vector{Float32}     # nucleation threshold on RHi (1.5 ≈ Koop)
    r_eff::Vector{Float32}        # effective ice radius [m] → v_t(T,p) per parcel
    h_layer::Vector{Float32}      # settling escape scale [m]
    q0::Vector{Float32}           # prescribed release mixing ratio [kg/kg]; 0 = plume/off
    frozen::Vector{Float64}       # per-box removed-mass accumulator [kg]
    # ── ice co-scavenging ([species.ice_scav]) ───────────────────────────────
    # Partner boxes (soluble gases riding the ice: SO2, HCl, HNO3 …) lose mass
    # in proportion to the FRACTION of the co-located water that froze out
    # this step: dm = m · uptake · f_shed(cell).  uptake ∈ [0,1] is the
    # retained-in-ice fraction (1 = fully retained dissolved tracer).
    scav_uptake::Vector{Float32}  # 0 = not a partner box
    scavenged::Vector{Float64}    # per-box co-scavenged-mass accumulator [kg]
    # ── DIAGNOSED settling escape scale ([physics].prognostic_airm) ───────────
    # Per-PARCEL local vertical spread of the parcel cloud [m], filled by
    # update_airm! (the h_v of its k nearest neighbours).  When > 0 it REPLACES
    # the configured per-box `h_layer`: the ice has to fall out of the layer the
    # plume actually occupies, not out of a prescribed 500 m.  Empty / 0 → the
    # configured value is used.  (HT: the model diagnoses an ~8 km thick plume,
    # so a 500 m escape scale sheds 6-7× too fast.)
    hv::Vector{Float32}
end

IcePhysics(on, rhi, reff, hl, q0, frozen) =
    IcePhysics(on, rhi, reff, hl, q0, frozen,
               zeros(Float32, length(on)), zeros(Float64, length(on)), Float32[])
# 8-arg form (with scavenging) kept working; `hv` defaults to empty = off
IcePhysics(on, rhi, reff, hl, q0, frozen, scav, scavd) =
    IcePhysics(on, rhi, reff, hl, q0, frozen, scav, scavd, Float32[])

any_ice(ice::IcePhysics) = any(ice.on) || any(>(0f0), ice.scav_uptake)

# met cell key for the co-scavenging accumulation (cell × level bin)
@inline _ice_cellkey(ii::Integer, jj::Integer, kk::Integer) =
    (Int64(ii) << 40) | (Int64(jj) << 20) | Int64(kk)

const _RHO_ICE = 917.0f0          # hexagonal ice [kg/m³]

# LEGACY plume-washout closure (large hydrometeors in the eruption column) —
# kept for the plume-phase code paths; NOT used for stratospheric ice, where
# it over-settles µm crystals ~30×.  Use settling_velocity(2r, 917, T, p).
ice_vterm(r_eff_m::Real) = 12.0 * sqrt(max(r_eff_m, 1e-7) / 1e-3)

"""
    apply_ice!(state, ctx, ice, dt_s, t_array) → removed [kg]

One freezing/sedimentation step for every active parcel of an ice-enabled
release box: T, p from the met at the parcel (θ/ρ interpolation, as the
chemistry extract_tp!); q = m/M_air; above `rhi_crit`·q_si the excess vapour
sheds at rate 1 − exp(−dt·v_t/h_layer), with v_t the Stokes–Cunningham fall
speed of an r_eff crystal at the PARCEL's T, p (slip correction grows with
altitude).  Removed mass accumulates in `ice.frozen` (the stratosphere→
troposphere ice flux of the tagged budget).

`tanom`: per-parcel radiative cold anomaly T′ ≤ 0 from
`[physics.plume_radiation].tprime_tau` — when given, e_si is evaluated at
T_ambient + T′ (ambient p): the plume's self-cooling reaches the freezing
threshold, which the ambient-θ inversion alone can never see (buoyancy
equilibrium under the hybrid vertical).
"""
function apply_ice!(state::State.ParcelState, ctx, ice::IcePhysics,
                    dt_s::Real, t_array;
                    tanom::Union{Nothing,AbstractVector} = nothing)
    nb = length(ice.on)
    total = Threads.Atomic{Float64}(0.0)
    boxtot  = [Threads.Atomic{Float64}(0.0) for _ in 1:nb]
    scavtot = [Threads.Atomic{Float64}(0.0) for _ in 1:nb]
    haspartner = any(>(0f0), ice.scav_uptake)
    nt = Threads.maxthreadid()
    for n in 1:length(ctx.nests)
        nestctx = ctx.nests[n]
        m3 = nestctx.metfields3D
        indt = Float32(t_array[n])
        n8 = eltype(state.nest)(n)
        # per-thread cell accumulators (merged below): water mass before the
        # step and ice mass shed, per met cell × level bin
        wat_d  = [Dict{Int64,Float64}() for _ in 1:nt]
        shed_d = [Dict{Int64,Float64}() for _ in 1:nt]
        Threads.@threads :static for i in eachindex(state.status)
            state.status[i] == 0x01 || continue
            state.nest[i] == n8     || continue
            b = Int(state.relbox[i])
            (1 <= b <= nb && ice.on[b]) || continue
            airm = Float64(state.airm[i])
            if airm <= 0.0
                ice.q0[b] > 0f0 || continue
                airm = Float64(state.mass[i]) / Float64(ice.q0[b])
                airm > 0.0 || continue
                state.airm[i] = eltype(state.airm)(airm)
            end
            tid = Threads.threadid()
            # θ/ρ → T, p: identical inversion to the chemistry extract_tp!
            # (θ = T·(P0/p)^(Rd/cp), p = ρ·Rd·T ⇒ T = θ^(cp/cv)·(ρ·Rd/P0)^(Rd/cv))
            di, dj, dz, dl, ii, jj, kk, l = _parcel_indices(
                state.indx[i], state.indy[i], state.indl[i], indt)
            key = _ice_cellkey(ii, jj, kk)
            if haspartner
                wat_d[tid][key] = get(wat_d[tid], key, 0.0) + Float64(state.mass[i])
            end
            θ = Float64(interp4(m3.theta, di, dj, dz, dl, ii, jj, kk, l))
            ρ = max(Float64(interp4(m3.rho, di, dj, dz, dl, ii, jj, kk, l)), 1e-6)
            T = max(θ^(1004.0 / 716.95) * (ρ * 287.05 / 1.0e5)^(287.05 / 716.95), 100.0)
            p = max(ρ * 287.05 * T, 1.0)
            # radiative cold anomaly T′: e_si at the PLUME temperature,
            # p stays ambient (constant-pressure anomaly)
            if tanom !== nothing && tanom[i] < 0.0f0
                T = max(T + Float64(tanom[i]), 100.0)
            end
            esi  = esat_ice(T)
            esi < 0.9 * p || continue                    # too warm/low: no ice
            qsi  = 0.622 * esi / (p - esi)               # kg/kg
            q    = Float64(state.mass[i]) / airm
            q > Float64(ice.rhi_crit[b]) * qsi || continue
            vt  = settling_velocity(2.0f0 * ice.r_eff[b], _RHO_ICE,
                                    Float32(T), Float32(p))
            # escape scale: the DIAGNOSED local plume thickness when available
            # (update_airm! fills ice.hv), else the configured per-box value.
            hlay = (!isempty(ice.hv) && ice.hv[i] > 0f0) ?
                   Float64(ice.hv[i]) : Float64(ice.h_layer[b])
            fsh = -expm1(-Float64(dt_s) * Float64(vt) / max(hlay, 1.0))
            dm  = (q - qsi) * airm * fsh
            dm  = min(dm, Float64(state.mass[i]))
            state.mass[i] -= eltype(state.mass)(dm)
            Threads.atomic_add!(total, dm)
            Threads.atomic_add!(boxtot[b], dm)
            haspartner &&
                (shed_d[tid][key] = get(shed_d[tid], key, 0.0) + dm)
        end
        # ── co-scavenging pass: partner parcels in cells where ice shed ──────
        if haspartner
            shed = shed_d[1]
            for t in 2:nt, (k, v) in shed_d[t]
                shed[k] = get(shed, k, 0.0) + v
            end
            isempty(shed) && continue
            wat = wat_d[1]
            for t in 2:nt, (k, v) in wat_d[t]
                wat[k] = get(wat, k, 0.0) + v
            end
            Threads.@threads :static for i in eachindex(state.status)
                state.status[i] == 0x01 || continue
                state.nest[i] == n8     || continue
                b = Int(state.relbox[i])
                (1 <= b <= nb && ice.scav_uptake[b] > 0f0) || continue
                di, dj, dz, dl, ii, jj, kk, l = _parcel_indices(
                    state.indx[i], state.indy[i], state.indl[i], indt)
                s = get(shed, _ice_cellkey(ii, jj, kk), 0.0)
                s > 0.0 || continue
                f  = min(s / max(get(wat, _ice_cellkey(ii, jj, kk), 0.0), 1e-300), 1.0)
                dm = min(Float64(state.mass[i]) * Float64(ice.scav_uptake[b]) * f,
                         Float64(state.mass[i]))
                state.mass[i] -= eltype(state.mass)(dm)
                Threads.atomic_add!(scavtot[b], dm)
            end
        end
    end
    @inbounds for b in 1:nb
        ice.frozen[b]    += boxtot[b][]
        ice.scavenged[b] += scavtot[b][]
    end
    return total[]
end
