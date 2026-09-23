# physics/transport_tags.jl
# ─────────────────────────────────────────────────────────────────────────────
# Passive per-parcel TRANSPORT-HISTORY tags ([tracers]) for full-profile air-mass
# attribution — the "how did this air get here" axis, complementing chemistry +
# the adjoint (the "what is it / where did the precursors come from" axis).
#
# Each parcel carries its own trajectory, so these are cheap linear diagnostics
# written each step from the met the parcel sits in, then gridded at output as
# cell-mean fields — e.g. the fraction of stratospheric-origin air at each
# altitude over a receptor.  They never perturb the transport or chemistry.
#
# Increment 1 — STRATOSPHERIC residence / STE:
#   strat_frac = (steps with θ ≥ strat_theta) / (active steps)     ∈ [0,1]
#     a proxy for time spent in stratospheric air (θ ≈ 380 K is the tropical
#     "overworld"/tropopause); dominates the 20–30 km profile, whose O3 comes
#     in from the boundary field and descends, not from in-situ production.
#   age_strat = time since the parcel was last stratospheric [s] (huge if never)
#     — recency of stratospheric contact, i.e. how freshly descended the air is.
#
# Increment 2 — PBL contact, MARINE/CONTINENTAL origin, CONVECTIVE history:
#   pbl_frac  = (steps with z_agl ≤ pblh) / (active steps)          ∈ [0,1]
#     fraction of the trajectory spent in boundary-layer contact — where surface
#     emissions, dry deposition and the diurnal cycle act on the air mass.
#   cont_frac = (PBL-contact steps over land) / (PBL-contact steps)  ∈ [0,1]
#     of the time in surface contact, the share over CONTINENT (lsm ≥ 0.5) vs
#     ocean — the marine-vs-continental origin axis.  0 if never in the PBL.
#   age_pbl   = time since the parcel was last in the PBL [s] (huge if never)
#     — recency of surface contact, i.e. how long the air has been aloft/aged.
#   conv_count = number of convective redistributions the parcel underwent
#     — a count of deep-convective vertical displacements (incremented by the
#     convection operator whenever a parcel's destination layer ≠ its origin).
#
# The PBL/marine tags need the 2-D pblh (live, typed nestctx.metfields2D.pblh)
# and the land-sea mask lsm (STATIC — captured once at construction from the raw
# archived 2-D dict, since the typed MetFields2D struct carries only PBL fields).
# conv_count is filled in-place by convective_redistribution! via the optional
# `conv_count` vector it is handed (no perturbation of the transport itself).
# ─────────────────────────────────────────────────────────────────────────────

mutable struct TransportTags
    strat_theta :: Float64                    # θ threshold for "stratospheric" [K]
    lsm         :: Union{Nothing,Matrix{Float32}}  # static land-sea mask (nest 1) or nothing
    # increment 1 (stratospheric)
    strat_steps :: Vector{Int32}    # steps with θ ≥ strat_theta (slot-indexed)
    n_steps     :: Vector{Int32}    # active steps (denominator for the fractions)
    age_strat   :: Vector{Float32}  # s since last stratospheric (AGE_NEVER if never)
    # increment 2 (PBL / marine-continental / convection)
    pbl_steps   :: Vector{Int32}    # steps with z_agl ≤ pblh (PBL contact)
    land_steps  :: Vector{Int32}    # PBL-contact steps over land (lsm ≥ 0.5)
    age_pbl     :: Vector{Float32}  # s since last PBL contact (AGE_NEVER if never)
    conv_count  :: Vector{Int32}    # # convective redistributions undergone
end

const _TAG_AGE_NEVER = 3.0f9        # ≫ any run length [s]: "never"
const _TAG_PBLH_FLOOR = 50.0f0      # [m] minimal usable pblh (below → treat as unset)

TransportTags(n_max::Int; strat_theta::Float64 = 380.0,
              lsm::Union{Nothing,Matrix{Float32}} = nothing) =
    TransportTags(strat_theta, lsm,
                  zeros(Int32, n_max), zeros(Int32, n_max), fill(_TAG_AGE_NEVER, n_max),
                  zeros(Int32, n_max), zeros(Int32, n_max), fill(_TAG_AGE_NEVER, n_max),
                  zeros(Int32, n_max))

any_tags(t::TransportTags) = true

"reset a (re)used slot's tags — call at release init, like airm/trad."
@inline function reset_tags!(t::TransportTags, i::Int)
    @inbounds begin
        t.strat_steps[i] = 0; t.n_steps[i] = 0; t.age_strat[i] = _TAG_AGE_NEVER
        t.pbl_steps[i]  = 0; t.land_steps[i] = 0; t.age_pbl[i] = _TAG_AGE_NEVER
        t.conv_count[i] = 0
    end
    return nothing
end

"""
    apply_transport_tags!(state, ctx, tags, dt_s, t_array; nest=1) → nactive

One tagging step: sample θ, the PBL height and the land-sea mask at each active
parcel and update the stratospheric / PBL-contact / marine-continental counters
and clocks.  Freed slots are reset (slot hygiene).  Returns the active count.
(conv_count is updated separately, in the convection operator.)
"""
function apply_transport_tags!(state::State.ParcelState, ctx, tags::TransportTags,
                               dt_s::Real, t_array; nest::Int = 1)
    nestctx = ctx.nests[nest]
    m3   = nestctx.metfields3D
    zlev = m3.zlev
    pblh = nestctx.metfields2D.pblh          # typed MetFields2D field (live via ctx)
    lsm  = tags.lsm                          # static mask captured at construction
    g    = nestctx.gridconfig
    Ni   = Int(g.Ni); Nj = Int(g.Nj)
    nz   = Int(nestctx.levelconfig.Nk)
    indt = Float32(t_array[nest])
    it   = clamp(round(Int, indt), 1, size(pblh, 3))
    n8   = eltype(state.nest)(nest)
    θc   = Float32(tags.strat_theta)
    dts  = Float32(dt_s)
    nact = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in eachindex(state.status)
        if state.status[i] != 0x01
            tags.strat_steps[i] = 0; tags.n_steps[i] = 0        # freed-slot hygiene
            tags.age_strat[i] = _TAG_AGE_NEVER
            tags.pbl_steps[i] = 0; tags.land_steps[i] = 0
            tags.age_pbl[i] = _TAG_AGE_NEVER; tags.conv_count[i] = 0
            continue
        end
        state.nest[i] == n8 || continue
        indx = state.indx[i]; indy = state.indy[i]
        indl = clamp(state.indl[i], 1.0f0, Float32(nz))
        @inbounds begin
            tags.n_steps[i] += Int32(1)
            # ── stratospheric (θ) ────────────────────────────────────────────
            θ = interpolation4Dsingle(m3.theta, indx, indy, indl, indt)
            if θ >= θc
                tags.strat_steps[i] += Int32(1)
                tags.age_strat[i] = 0.0f0
            else
                tags.age_strat[i] += dts
            end
            # ── PBL contact + marine/continental ─────────────────────────────
            hpbl = Float32(pblh[clamp(round(Int, indx), 1, Ni),
                                clamp(round(Int, indy), 1, Nj), it])
            if hpbl >= _TAG_PBLH_FLOOR
                z_amsl = interpolation4Dsingle(zlev, indx, indy, indl, indt)
                z_sfc  = interpolation4Dsingle(zlev, indx, indy, 1.0f0, indt)
                z_agl  = z_amsl - z_sfc
                if z_agl <= hpbl
                    tags.pbl_steps[i] += Int32(1)
                    tags.age_pbl[i] = 0.0f0
                    if lsm !== nothing
                        ix = clamp(round(Int, indx), 1, size(lsm, 1))
                        iy = clamp(round(Int, indy), 1, size(lsm, 2))
                        lsm[ix, iy] >= 0.5f0 && (tags.land_steps[i] += Int32(1))
                    end
                else
                    tags.age_pbl[i] += dts
                end
            else
                tags.age_pbl[i] += dts
            end
        end
        Threads.atomic_add!(nact, 1)
    end
    return nact[]
end

"per-parcel stratospheric fraction ∈ [0,1] (0 if the parcel has no steps yet)."
@inline strat_fraction(t::TransportTags, i::Int) =
    @inbounds t.n_steps[i] > 0 ? Float32(t.strat_steps[i]) / Float32(t.n_steps[i]) : 0.0f0

"per-parcel PBL-contact fraction ∈ [0,1] (0 if no steps yet)."
@inline pbl_fraction(t::TransportTags, i::Int) =
    @inbounds t.n_steps[i] > 0 ? Float32(t.pbl_steps[i]) / Float32(t.n_steps[i]) : 0.0f0

"per-parcel continental fraction of PBL-contact time ∈ [0,1] (0 if never in PBL)."
@inline cont_fraction(t::TransportTags, i::Int) =
    @inbounds t.pbl_steps[i] > 0 ? Float32(t.land_steps[i]) / Float32(t.pbl_steps[i]) : 0.0f0
