
# ── Per-species process physics (decay, gas dry deposition, wet-dep coeffs) ──
#
# The species DB / TOML carry per-species parameters (decay half-life, gas
# dry-deposition reldiff/henry/f0 or a prescribed vd, below-cloud scavenging
# A/B) that used to be parsed and dropped.  SpeciesPhysics maps them onto the
# RELEASE-BOX index (state.relbox — each box releases one species), which is
# what a parcel actually carries, so the per-parcel loops stay integer lookups.
#
#   apply_decay!          m ← m·exp(−ln2·Δt/T½)             (radioactive decay)
#   dry_deposition_gas!   resistance model vd = 1/(Ra+Rb+Rc) in the lowest href
#                         metres, or a prescribed vd; removed mass accumulated
#                         on a per-nest grid (like the aerosol path).
#   wet-dep A/B           consumed by _wetdep_nest! (wetdep.jl) for the
#                         below-cloud pathway; the in-cloud (nucleation)
#                         pathway keeps the bulk coefficients until CCN/IN
#                         scavenging is wired.
#
# Rc closure (no land-use database yet): Wesely (1989)-style bulk surface
# resistance calibrated for WATER surfaces,
#     Rc = 1 / (H*/3000 + f0/2000)      [s m⁻¹]
# (SO2: H*≈1e5 → Rc≈0.03·10³/H*… ≈30 s/m, vd≈0.8 cm/s; O3: f0=1 → Rc≈2000 s/m,
# vd≈0.05 cm/s — both match open-ocean observations).  Over land this is only
# indicative; use the species' prescribed `dryvel` or extend with a land-use
# lookup when land deposition matters.

"""
Per-release-box species process parameters (decay λ, below-cloud wet A/B, gas
dry-deposition inputs) indexed by `state.relbox`, plus the per-nest
dry-deposition accumulation grids.  Built by the runtime from the species
DB / TOML; consumed by `apply_decay!`, `dry_deposition_gas!` and
`_wetdep_nest!`.
"""
struct SpeciesPhysics
    lambda ::Vector{Float32}    # decay constant ln2/T½ [1/s] per relbox (0 = stable)
    weta   ::Vector{Float32}    # below-cloud scavenging A per relbox (<0 → bulk default)
    wetb   ::Vector{Float32}    # below-cloud scavenging B per relbox
    drydep ::Vector{Bool}       # gas dry deposition active for this box
    vd     ::Vector{Float32}    # prescribed vd [m/s] per relbox (<0 → resistance model)
    reldiff::Vector{Float32}    # D(H2O)/D(gas) per relbox (Rb)
    henry  ::Vector{Float32}    # Henry's-law constant [M/atm] (Rc)
    f0     ::Vector{Float32}    # Wesely reactivity 0–1 (Rc)
    href   ::Float32            # dry-dep application depth [m]
    depgrid::Vector{Matrix{Float64}}   # per-nest dry-deposited mass [kg per met cell]
    # ── [species.settling]: GRAVITATIONAL FALL of a condensed-phase species ──
    # A tagged condensed species (volcanic ice, ash …) does not follow the air:
    # it falls at its Stokes–Cunningham terminal velocity.  Needed to represent
    # sedimenting ice as its OWN trajectory population — ice that condenses at
    # 28 km re-evaporates several km LOWER, not where it formed, so a reservoir
    # attached to the parent parcel would put the water back at the wrong level.
    settle_r  ::Vector{Float32}  # effective particle RADIUS [m] (0 = no settling)
    settle_rho::Vector{Float32}  # particle density [kg/m³] (ice: 917)
end

"True when at least one release box has a nonzero decay constant."
any_decay(sp::SpeciesPhysics)  = any(>(0f0), sp.lambda)
"True when at least one release box has gas dry deposition enabled."
any_drydep(sp::SpeciesPhysics) = any(sp.drydep)

"""
    apply_decay!(state, sp, dt_s) → total mass decayed

Radioactive decay of parcel mass, per release box.  Skipped boxes (λ=0) cost
one lookup.  relbox 0 (domain-fill / boundary parcels) never decays.
"""
function apply_decay!(state::State.ParcelState, sp::SpeciesPhysics, dt_s::Real)
    nb = length(sp.lambda)
    dt = Float32(dt_s)
    total = 0.0
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        b = Int(state.relbox[i])
        (1 <= b <= nb && sp.lambda[b] > 0f0) || continue
        dm = state.mass[i] * (-expm1(-sp.lambda[b] * dt))
        state.mass[i] -= dm
        total += Float64(dm)
    end
    return total
end

# Gas quasi-laminar resistance: Rb = (2/κu*)·(Sc/Pr)^(2/3), Sc = reldiff·ν/D_H2O.
@inline function _rb_gas(reldiff::Float32, nu::Float32, ustar::Float32)
    D_h2o = 2.4f-5                       # H2O diffusivity in air [m²/s]
    sc    = reldiff * nu / D_h2o
    return clamp((2.0f0 / (0.4f0 * ustar)) * (sc / 0.72f0)^(2.0f0/3.0f0), 0.1f0, 1.0f4)
end

# Bulk (water-surface) Wesely-style surface resistance from H* and f0.
@inline _rc_gas(henry::Float32, f0::Float32) =
    clamp(1.0f0 / (max(henry, 0f0) / 3000.0f0 + max(f0, 0f0) / 2000.0f0 + 1f-8),
          1.0f0, 1.0f5)

"""
    dry_deposition_gas!(state, ctx, sp, dt_s, t_array) → total mass removed

Resistance-model (or prescribed-vd) gas dry deposition for parcels in the
lowest `sp.href` metres: m ← m·exp(−vd·Δt/href); removed mass accumulates on
`sp.depgrid[n]`.  Skipped in backward runs.
"""
function dry_deposition_gas!(state::State.ParcelState, ctx, sp::SpeciesPhysics,
                             dt_s::Real, t_array)
    Int(ctx.runconfig.direction) == -1 && return 0.0
    dt = Float32(dt_s)
    nb = length(sp.lambda)
    total = Threads.Atomic{Float64}(0.0)
    for n in 1:length(ctx.nests)
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8 &&
                           1 <= Int(state.relbox[i]) <= nb && sp.drydep[Int(state.relbox[i])],
                      eachindex(state.status))
        isempty(idx) && continue
        m3 = nestctx.metfields3D; m2 = nestctx.metfields2D
        nz = size(m3.zlev, 3)
        Ni = Int(nestctx.gridconfig.Ni); Nj = Int(nestctx.gridconfig.Nj)
        indt = t_array[n]
        dep  = sp.depgrid[n]
        nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
        zbufs   = [Vector{Float32}(undef, nz) for _ in 1:nthr]
        removed = zeros(Float64, nthr)
        Threads.@threads :static for i in idx
            tid  = Threads.threadid()
            indx = state.indx[i]; indy = state.indy[i]
            (indx < 1.0f0 || indx >= Float32(Ni) ||
             indy < 1.0f0 || indy >= Float32(Nj)) && continue
            z = zbufs[tid]
            get_vertical_profile!(z, m3.zlev, indx, indy, indt, nz)
            zp = _z_at(z, clamp(state.indl[i], 1.0f0, Float32(nz - 1)), nz)
            zp > sp.href && continue
            b = Int(state.relbox[i])
            vd = sp.vd[b]
            if vd < 0f0                              # resistance model
                # lowest full level = slot 2, read at indl 1.5 (layout contract)
                T  = interpolation4Dsingle(m3.T,   indx, indy, 1.5f0, indt)
                p  = interpolation4Dsingle(m3.puv, indx, indy, 1.5f0, indt)
                us = interpolate2Dsingle(m2.ustar, indx, indy, indt)
                ol = interpolate2Dsingle(m2.ol,    indx, indy, indt)
                (isfinite(T) && isfinite(p) && T > 150.0f0) || continue
                isfinite(us) || (us = 0.2f0); us = max(us, 0.05f0)
                isfinite(ol) || (ol = 1.0f6)
                mu   = 1.458f-6 * T * sqrt(T) / (T + 110.4f0)
                nu   = mu / (p / (287.05f0 * T))
                aol  = abs(ol) < 1.0f0 ? (ol < 0.0f0 ? -1.0f0 : 1.0f0) : ol
                ra   = (log(sp.href / 0.03f0) - _psi_h(sp.href / aol) + _psi_h(0.03f0 / aol)) /
                       (0.4f0 * us)
                isfinite(ra) || (ra = 100.0f0)
                ra   = clamp(ra, 1.0f0, 1.0f4)
                vd   = 1.0f0 / (ra + _rb_gas(sp.reldiff[b], nu, us) +
                                _rc_gas(sp.henry[b], sp.f0[b]))
            end
            vd > 0f0 || continue
            frem = -expm1(-vd * dt / sp.href)
            dm = state.mass[i] * frem
            state.mass[i] -= dm
            removed[tid] += Float64(dm)
            ii = clamp(round(Int, indx), 1, Ni)
            jj = clamp(round(Int, indy), 1, Nj)
            Atomix.@atomic dep[ii, jj] += Float64(dm)
        end
        Threads.atomic_add!(total, sum(removed))
    end
    return total[]
end

# ── [species.convert]: in-flight stochastic species conversion ───────────────
# Each parcel of a converting box switches to the daughter's SHADOW box with
# probability 1 − exp(−dt/τ) per step (per-parcel xorshift stream — thread-
# safe and reproducible), rescaling its mass by M_daughter/M_parent.  Mass is
# conserved in expectation; the daughter inherits aerosol physics/deposition/
# output attribution through its own relbox.  Ice `airm` is reset (the
# closure belongs to the water species).
struct SpeciesConvert
    lambda::Vector{Float32}     # 1/τ per box (0 = no conversion)
    target::Vector{Int32}       # parent box → shadow box (0 = none)
    ratio::Vector{Float32}      # mass rescale on conversion
    converted::Vector{Float64}  # per-parent-box converted-mass accumulator [kg]
    # ── THERMODYNAMIC gating ([species.convert].mode) ────────────────────────
    # 0 = plain first-order τ (legacy).
    # 1 = FREEZE   : convert only where q > rhi_crit·q_si (vapour → ice).
    # 2 = SUBLIMATE: convert only where q < q_si          (ice → vapour).
    # The pair makes the phase change REVERSIBLE: an ice parcel that falls into
    # warmer, subsaturated air gives its water back — instead of the one-way
    # ratchet of the removal-only closure.  q = mass/airm, q_si over ice at the
    # parcel's own (T, p) — identical inversion to the ice module.
    mode::Vector{Int8}
    rhi_crit::Vector{Float32}   # freeze threshold on RHi (mode 1)
    # ── SUB-GRID LAMINA ([species.convert].lamina_sigma) ─────────────────────
    # q = mass/airm is the parcel MEAN.  The real air the parcel represents is
    # not uniform: the plume is laminated, and ice lives in the saturated
    # laminae.  Gating on the mean makes every parcel flash-sublimate, because
    # at 20 hPa/220 K a 5 µm crystal sublimates in 0.65 s once exposed — ~460
    # e-foldings per 5-minute step.  Ice and the saturated lamina are the SAME
    # OBJECT; a model that ignores the sub-grid structure returns "no ice" for
    # the same reason a 3 km satellite average does.
    #
    # Demonstrated 2026-08-02: with five ice size bins spanning a 25x range in
    # fall speed, ALL bins retained 2.4-3.2 % — i.e. the loss was entirely
    # size-INDEPENDENT sublimation, not sedimentation.  Size resolution alone
    # cannot fix it.
    #
    # Treatment: assume a lognormal sub-grid distribution of q about the parcel
    # mean with relative width `lam_sigma`.  The ELIGIBLE FRACTION of the parcel
    # (supersaturated for freezing, subsaturated for sublimation) then scales
    # the conversion probability.  lam_sigma = 0 reproduces the old step
    # function exactly, so existing configs are bit-identical.
    lam_sigma::Vector{Float32}
    # ── DIFFUSION-LIMITED SUBLIMATION ([species.convert].diffusion_limited) ──
    # The prescribed `tau` is a single number for every crystal size and every
    # humidity.  The real timescale for a sublimating crystal is
    #     t = r^2 * rho_ice / (2 * D_v * (1-S) * rho_sat)
    # i.e. it goes as r^2 AND as 1/(1-S).  At 20 hPa/220 K in subsaturated air
    # that is 0.02-8 s across our bins — everything vanishes inside one timestep
    # and SIZE CANNOT MATTER.  Only near saturation (S -> 1) does r^2 separate
    # the bins, which is exactly where the lamina puts the ice.  So the size
    # dependence and the lamina are the SAME mechanism: both act through (1-S),
    # and a fixed tau destroys both.  With this on, small crystals evaporate
    # first, the surviving r_eff RISES, and tau_optical = 1.5*IWP/(rho*r) falls
    # from BOTH ends — the mechanism that makes the plume go subvisible.
    # radius[b] > 0 enables it for that box (taken from [species.fall].r_eff).
    subl_radius::Vector{Float32}
    # ── SECTIONAL: per-parcel radius + BIN-TO-BIN transfer ([species.convert]
    # .sectional) ────────────────────────────────────────────────────────────
    # The exponential p = 1 − exp(−dt/t_sub) SATURATES: with t_sub = 0.02–8 s
    # against a 300 s step it is 1 for every bin, which ERASES the r² law it
    # was built to express (8-day evidence, ht_full3: r_eff flat 3.10–3.44 µm
    # over the whole run, bin fractions constant).  Neither the lamina nor the
    # T′ anomaly can fix that: both act on (1−S), which multiplies every bin
    # ALIKE.  Size selectivity has to come from the SIZE COORDINATE itself.
    #
    # Sectional treatment (ICE_UPGRADE_DESIGN §2), integrated exactly rather
    # than through a rate constant.  For diffusional sublimation
    #     r dr/dt = −K,   K = D_v(1−S)ρ_sat/ρ_ice   ⇒   r'² = r² − 2K·dt
    # so over one step a crystal keeps the mass fraction
    #     m'/m = (r'/r)³ = (1 − dt/t_sub)^{3/2},    t_sub = r²/(2K)
    # and loses the rest to vapour.  This is the SAME physics the exponential
    # approximated, but it stays size-dependent for every dt: a bin with
    # dt ≥ t_sub is annihilated while a bin with t_sub ≫ dt barely shrinks,
    # and t_sub ∝ r² spans a factor ~400 across our five bins.
    #
    # `rice` is the parcel's OWN crystal radius [m] (slot-indexed, seeded from
    # the box r_eff on first use), so the radius EVOLVES instead of being
    # pinned to the box.  When it shrinks past the next bin edge the parcel is
    # demoted to the smaller ICE box (`bin_down`); below the smallest bin the
    # crystal is gone and the parcel converts to vapour.  Settling then reads
    # the same radius, so sedimentation becomes size-resolved as well.
    # Empty `rice` (sectional off) = the legacy exponential path, unchanged.
    rice::Vector{Float32}        # per-PARCEL crystal radius [m]; empty = off
    bin_down::Vector{Int32}      # ICE box → next SMALLER ICE box (0 = smallest)
    bin_redge::Vector{Float32}   # ICE box → lower radius edge [m] of that bin
    # ── FLUX / GATE DIAGNOSTICS (2026-08-05 pre-registration) ────────────────
    # Per-box accumulators [kg], zero-cost when no event fires.  `demoted` is
    # the mass RELABELED one bin down (a pure relabel — the parcel keeps its
    # mass); `gone` is the crystal-vanished subset of `converted` (r'² ≤ 0
    # within a step, or shrinkage off the bottom of the ladder).
    demoted::Vector{Float64}
    gone::Vector{Float64}
    # Freeze-event moments, per SOURCE (vapour) box: Σ over mode-1 conversions
    # of m·X.  X = ambient T (no T′), the T′ actually applied, ln q, ln q_si
    # (at the gate's T), ln p, and the STAGGERED−PLAIN θ/ρ-read temperature
    # offset at the gate — the §0c.2-item-7 vertical-convention split.  Since
    # the 2026-10-22 layout contract settling and the hybrid ladder read θ
    # stacked like the gate, so fz_mdT no longer measures a live split: it is
    # the error the pre-contract plain read made, kept for comparability.
    fz_mT::Vector{Float64}
    fz_mTp::Vector{Float64}
    fz_mlq::Vector{Float64}
    fz_mlqsi::Vector{Float64}
    fz_mlp::Vector{Float64}
    fz_mdT::Vector{Float64}
    # Per-parcel gate state at the LAST convert step that evaluated this
    # parcel (slot-indexed; empty = diagnostics off).  gate_T includes T′.
    gate_T::Vector{Float32}
    gate_qsi::Vector{Float32}
    gate_p::Vector{Float32}
end
# 11-arg form (the pre-diagnostics signature): allocates the accumulators;
# per-parcel gate vectors only when `diag_slots` (the parcel-pool size) is set.
function SpeciesConvert(lambda, target, ratio, converted, mode, rhi, lsg, sbr,
                        rice, bdn, bre; diag_slots::Integer = 0)
    nb = length(lambda)
    gv() = diag_slots > 0 ? zeros(Float32, diag_slots) : Float32[]
    SpeciesConvert(lambda, target, ratio, converted, mode, rhi, lsg, sbr,
                   rice, bdn, bre,
                   zeros(Float64, nb), zeros(Float64, nb),
                   zeros(Float64, nb), zeros(Float64, nb), zeros(Float64, nb),
                   zeros(Float64, nb), zeros(Float64, nb), zeros(Float64, nb),
                   gv(), gv(), gv())
end
# legacy 4-arg form: plain τ conversion, no thermodynamic gate
SpeciesConvert(lambda, target, ratio, converted) =
    SpeciesConvert(lambda, target, ratio, converted,
                   zeros(Int8, length(lambda)), fill(1.5f0, length(lambda)),
                   zeros(Float32, length(lambda)), zeros(Float32, length(lambda)))
SpeciesConvert(lambda, target, ratio, converted, mode, rhi) =
    SpeciesConvert(lambda, target, ratio, converted, mode, rhi,
                   zeros(Float32, length(lambda)), zeros(Float32, length(lambda)))
# 8-arg form (lamina + diffusion-limited radius, no sectional state)
SpeciesConvert(lambda, target, ratio, converted, mode, rhi, lsg, sbr) =
    SpeciesConvert(lambda, target, ratio, converted, mode, rhi, lsg, sbr,
                   Float32[], zeros(Int32, length(lambda)),
                   zeros(Float32, length(lambda)))

"True when the sectional (per-parcel radius + bin transfer) path is armed."
any_sectional(cv::SpeciesConvert) = !isempty(cv.rice)

"""
Fraction of a lognormal sub-grid distribution (mean `q`, relative width `σ`)
lying above the threshold `qt`.  σ = 0 → step function.
"""
@inline function _lamina_frac(q::Float64, qt::Float64, σ::Float64)
    σ > 1e-6 || return q > qt ? 1.0 : 0.0
    (q > 0.0 && qt > 0.0) || return 0.0
    # median of the lognormal is q·exp(-σ²/2) so that its MEAN is q
    z = (log(qt / q) + 0.5 * σ * σ) / (σ * sqrt(2.0))
    return 0.5 * _erfc(z)
end

"Abramowitz & Stegun 7.1.26 complementary error function (|ε| < 1.5e-7)."
@inline function _erfc(x::Float64)
    s = sign(x); a = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * a)
    y = 1.0 - (((((1.061405429t - 1.453152027) * t) + 1.421413741) * t
                - 0.284496736) * t + 0.254829592) * t * exp(-a * a)
    return s >= 0 ? 1.0 - y : 1.0 + y
end
any_thermo_convert(cv::SpeciesConvert) = any(!=(Int8(0)), cv.mode)

any_convert(cv::SpeciesConvert) = any(>(0), cv.target)

function apply_convert!(state::State.ParcelState, cv::SpeciesConvert, dt_s::Real)
    nb = length(cv.target)
    tot = 0.0
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        b = Int(state.relbox[i])
        (1 <= b <= nb && cv.target[b] > 0) || continue
        p = -expm1(-Float64(dt_s) * Float64(cv.lambda[b]))
        _conv_rand(state, i) <= p || continue
        m0 = Float64(state.mass[i])
        state.relbox[i] = eltype(state.relbox)(cv.target[b])
        state.mass[i]   = eltype(state.mass)(m0 * Float64(cv.ratio[b]))
        state.airm[i]   = eltype(state.airm)(0)
        cv.converted[b] += m0
        tot += m0
    end
    return tot
end

# ─────────────────────────────────────────────────────────────────────────────
# [species.convert].mode : THERMODYNAMICALLY GATED phase change
# ─────────────────────────────────────────────────────────────────────────────
"""
    apply_convert_thermo!(state, ctx, cv, dt_s, t_array; nest=1) → converted_kg

Same stochastic box switch as `apply_convert!`, but the conversion is allowed
only where the local thermodynamics permits it:

  * mode 1 (FREEZE)    : q > rhi_crit·q_si — supersaturated over ice
  * mode 2 (SUBLIMATE) : q < q_si          — subsaturated, the ice gives back

with q = mass/airm and q_si over ice at the parcel's own (T, p).  Boxes with
mode 0 fall back to the plain τ conversion.

`tanom`: the per-parcel radiative cold anomaly T′ ≤ 0 from
`[physics.plume_radiation].tprime_tau` (PlumeRadiation.tanom).  When given,
the saturation is evaluated at the PLUME temperature, e_si(T_ambient + T′) at
ambient p — without it the plume's self-cooling never reaches the gate: the
hybrid vertical keeps the parcel in buoyancy equilibrium, so the ambient-θ
inversion below always returns the ERA5 temperature (a −5 K anomaly at 220 K
is a factor ~1.9 on RHi — HT_ICE_FINDINGS §5).

Together with `[species.settling]` on the ice box this makes the ice a real
transported population: it condenses aloft, FALLS, and sublimates lower down
where the air is subsaturated — the water returns to the atmosphere at the
level the ice actually reached, and the budget is reversible.
"""
function apply_convert_thermo!(state::State.ParcelState, ctx, cv::SpeciesConvert,
                               dt_s::Real, t_array; nest::Int = 1,
                               tanom::Union{Nothing,AbstractVector} = nothing)
    nb = length(cv.target)
    nestctx = ctx.nests[nest]
    m3 = nestctx.metfields3D
    indt = Float32(t_array[nest])
    n8 = eltype(state.nest)(nest)
    tot = 0.0
    sect = any_sectional(cv)
    # buffers for the PLAIN-convention θ/ρ read at freeze events (fz_mdT) —
    # this loop is serial, one pair suffices
    nzc = Int(nestctx.levelconfig.Nk)
    thbuf = zeros(Float32, nzc); rhobuf = zeros(Float32, nzc)
    diagp = !isempty(cv.gate_T)
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        state.nest[i] == n8     || continue
        b = Int(state.relbox[i])
        (1 <= b <= nb && cv.target[b] > 0) || continue
        md = cv.mode[b]
        # sectional: seed this parcel's crystal radius from its box the first
        # time it is seen as ice (also re-seeds a parcel that just re-froze,
        # whose radius was zeroed when its crystal vanished)
        if sect && md == Int8(2) && cv.rice[i] <= 0f0 && cv.subl_radius[b] > 0f0
            cv.rice[i] = cv.subl_radius[b]
        end
        lam_scale = 1.0
        q_save = qsi_save = T_save = p_save = esi_save = Ta_save = 0.0
        if md != Int8(0)
            airm = Float64(state.airm[i])
            airm > 0.0 || continue                    # no closure yet: skip
            θ = Float64(interpolation4Dsingle(m3.theta, state.indx[i], state.indy[i],
                                              state.indl[i], indt))
            ρ = max(Float64(interpolation4Dsingle(m3.rho, state.indx[i], state.indy[i],
                                                  state.indl[i], indt)), 1e-6)
            T = max(θ^(1004.0 / 716.95) * (ρ * 287.05 / 1.0e5)^(287.05 / 716.95), 100.0)
            p = max(ρ * 287.05 * T, 1.0)
            Tamb = T                                  # ambient, before T′
            # radiative cold anomaly T′: e_si at the PLUME temperature, p stays
            # ambient (the anomaly is at constant pressure)
            if tanom !== nothing && tanom[i] < 0.0f0
                T = max(T + Float64(tanom[i]), 100.0)
            end
            esi = esat_ice(T)
            esi < 0.9 * p || continue
            qsi = 0.622 * esi / (p - esi)
            q   = Float64(state.mass[i]) / airm
            q_save, qsi_save, T_save, p_save, esi_save = q, qsi, T, p, esi
            Ta_save = Tamb
            if diagp
                cv.gate_T[i] = Float32(T); cv.gate_qsi[i] = Float32(qsi)
                cv.gate_p[i] = Float32(p)
            end
            σ   = Float64(cv.lam_sigma[b])
            if md == Int8(1)
                # FREEZE: the eligible fraction is the SUPERSATURATED tail
                fel = _lamina_frac(q, Float64(cv.rhi_crit[b]) * qsi, σ)
            else
                # SUBLIMATE: only the ice OUTSIDE the saturated lamina is
                # exposed, i.e. the SUBSATURATED fraction 1 − f_lam
                fel = 1.0 - _lamina_frac(q, qsi, σ)
            end
            fel > 1e-6 || continue
            lam_scale = fel
        end
        lam_eff = Float64(cv.lambda[b])
        # ── SECTIONAL path: integrate r² = r₀² − 2K·dt exactly ───────────────
        # Size selectivity survives here because it lives in the RADIUS, not in
        # a rate constant that saturates.  The parcel loses the mass fraction
        # 1 − (r'/r)³ to vapour and KEEPS the rest as smaller crystals; when r'
        # crosses the bin's lower edge the parcel is demoted a bin.
        if sect && md == Int8(2) && qsi_save > 0.0
            r0 = Float64(cv.rice[i])
            r0 > 0.0 || (r0 = Float64(cv.subl_radius[b]))
            if r0 > 0.0
                S  = min(q_save / qsi_save, 0.999999)
                Dv = 2.11e-5 * (1.01325e5 / p_save) * (T_save / 273.15)^1.94
                rs = esi_save / (461.5 * T_save)
                K  = Dv * max(1.0 - S, 0.0) * rs / 917.0        # m²/s
                # only the EXPOSED fraction sublimates (lamina; σ=0 → 0 or 1)
                r2 = r0 * r0 - 2.0 * K * Float64(dt_s) * lam_scale
                if r2 <= 0.0
                    # crystal gone within the step → all mass to vapour
                    m0 = Float64(state.mass[i])
                    state.relbox[i] = eltype(state.relbox)(cv.target[b])
                    state.mass[i]   = eltype(state.mass)(m0 * Float64(cv.ratio[b]))
                    state.airm[i]   = eltype(state.airm)(0)
                    cv.rice[i]      = 0.0f0
                    cv.converted[b] += m0
                    cv.gone[b]      += m0
                    tot += m0
                    continue
                end
                rn = sqrt(r2)
                # mass lost to vapour this step: m·(1 − (r'/r)³).  Applied
                # DETERMINISTICALLY (no stochastic all-or-nothing), so the
                # r² law is not washed out by sampling noise.  The vapour is
                # credited to the daughter box's budget; the parcel keeps its
                # ice identity and its remaining mass.
                fkeep = (rn / r0)^3
                m0    = Float64(state.mass[i])
                dm    = m0 * (1.0 - fkeep)
                if dm > 0.0
                    state.mass[i] = eltype(state.mass)(m0 - dm)
                    cv.converted[b] += dm
                    tot += dm
                end
                cv.rice[i] = Float32(rn)
                # bin demotion: r' below this bin's lower edge → smaller box
                if cv.bin_redge[b] > 0f0 && rn < Float64(cv.bin_redge[b])
                    nb2 = Int(cv.bin_down[b])
                    if nb2 > 0
                        state.relbox[i] = eltype(state.relbox)(nb2)
                        cv.demoted[b]  += Float64(state.mass[i])
                    else
                        # off the bottom of the ladder: the crystal is gone
                        mr = Float64(state.mass[i])
                        state.relbox[i] = eltype(state.relbox)(cv.target[b])
                        state.mass[i]   = eltype(state.mass)(mr * Float64(cv.ratio[b]))
                        state.airm[i]   = eltype(state.airm)(0)
                        cv.rice[i]      = 0.0f0
                        cv.converted[b] += mr
                        cv.gone[b]      += mr
                        tot += mr
                    end
                end
                continue
            end
        end
        if md == Int8(2) && cv.subl_radius[b] > 0f0 && qsi_save > 0.0
            # diffusion-limited: t = r^2*rho_ice/(2*Dv*(1-S)*rho_sat)
            S  = min(q_save / qsi_save, 0.999999)
            Dv = 2.11e-5 * (1.01325e5 / p_save) * (T_save / 273.15)^1.94
            rs = esi_save / (461.5 * T_save)
            tsub = Float64(cv.subl_radius[b])^2 * 917.0 /
                   (2.0 * Dv * max(1.0 - S, 1e-6) * rs)
            lam_eff = 1.0 / max(tsub, 1.0e-3)
        end
        pconv = -expm1(-Float64(dt_s) * lam_eff)
        md == Int8(0) || (pconv *= lam_scale)
        _conv_rand(state, i) <= pconv || continue
        m0 = Float64(state.mass[i])
        # FREEZE-event moments (pre-registration step 1): the thermodynamic
        # environment of every vapour→ice conversion, mass-weighted.  fz_mdT
        # measures the staggered−plain read offset at the gate: the gate's
        # `interpolation4Dsingle` samples θ/ρ half a level ABOVE the plain
        # read that settling and the hybrid ladder made before the 2026-10-22
        # layout contract (§0c.2 item 7); both now read stacked like the gate.
        if md == Int8(1) && qsi_save > 0.0
            cv.fz_mT[b]    += m0 * Ta_save
            cv.fz_mTp[b]   += m0 * (T_save - Ta_save)
            cv.fz_mlq[b]   += m0 * log(q_save)
            cv.fz_mlqsi[b] += m0 * log(qsi_save)
            cv.fz_mlp[b]   += m0 * log(p_save)
            get_vertical_profile!(thbuf, m3.theta, state.indx[i], state.indy[i],
                                  indt, nzc)
            get_vertical_profile!(rhobuf, m3.rho, state.indx[i], state.indy[i],
                                  indt, nzc)
            il  = clamp(state.indl[i], 1.0f0, Float32(nzc - 1))
            θpl = Float64(_z_at(thbuf, il, nzc))
            ρpl = max(Float64(_z_at(rhobuf, il, nzc)), 1e-6)
            Tpl = max(θpl^(1004.0 / 716.95) *
                      (ρpl * 287.05 / 1.0e5)^(287.05 / 716.95), 100.0)
            cv.fz_mdT[b]   += m0 * (Ta_save - Tpl)
        end
        state.relbox[i] = eltype(state.relbox)(cv.target[b])
        state.mass[i]   = eltype(state.mass)(m0 * Float64(cv.ratio[b]))
        state.airm[i]   = eltype(state.airm)(0)
        cv.converted[b] += m0
        tot += m0
    end
    return tot
end

# ─────────────────────────────────────────────────────────────────────────────
# [species.settling] : gravitational fall of a condensed-phase tagged species
# ─────────────────────────────────────────────────────────────────────────────
any_settling(sp::SpeciesPhysics) = any(>(0f0), sp.settle_r)

"""
    apply_species_settling!(state, ctx, sp, dt_s, t_array; nest=1) → n_moved

Move every parcel of a settling box DOWN by v_t·dt, with v_t the
Stokes–Cunningham terminal velocity at the parcel's own (T, p) — the same
`settling_velocity` the ice module uses.  This is what turns the ice shadow
species into its own falling trajectory population.

`thetap`: the hybrid/isentropic scheme's PROGNOSTIC parcel θ.  It MUST be
passed whenever `[physics.vertical] mode = "hybrid"/"isentropic"` is active.
Sedimentation is a CROSS-ISENTROPIC motion — the crystal falls THROUGH the air,
it does not carry its isentrope with it — so the parcel's θ has to be reset to
the ambient θ of the level it lands on.  Without this the next hybrid step
re-inverts θ_p (isentropic.jl: indl ← indl + w·(indl_θ − indl)) and snaps the
parcel straight back up: with w → 1 in the stratosphere the fall is cancelled
EXACTLY, settling becomes a silent no-op, and r_eff has no effect at all.
(Diagnosed 2026-07-31 on the HT r_eff bracket: 10/30/50 µm gave an identical
ice centroid, 26.8→26.2 km in 24 h, where 50 µm owes ~30 km/day.)
"""
function apply_species_settling!(state::State.ParcelState, ctx, sp::SpeciesPhysics,
                                 dt_s::Real, t_array; nest::Int = 1,
                                 thetap::Union{Nothing,AbstractVector} = nothing,
                                 rice::Union{Nothing,AbstractVector} = nothing)
    nb = length(sp.settle_r)
    nestctx = ctx.nests[nest]
    m3 = nestctx.metfields3D
    nz = Int(nestctx.levelconfig.Nk)
    indt = Float32(t_array[nest])
    n8 = eltype(state.nest)(nest)
    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    zbufs = [zeros(Float32, nz) for _ in 1:nthr]
    thbufs = [zeros(Float32, nz) for _ in 1:nthr]
    nmov = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        state.nest[i] == n8     || continue
        b = Int(state.relbox[i])
        (1 <= b <= nb && sp.settle_r[b] > 0f0) || continue
        θ = Float64(interpolation4Dsingle(m3.theta, state.indx[i], state.indy[i],
                                          state.indl[i], indt))
        ρ = max(Float64(interpolation4Dsingle(m3.rho, state.indx[i], state.indy[i],
                                              state.indl[i], indt)), 1e-6)
        T = max(θ^(1004.0 / 716.95) * (ρ * 287.05 / 1.0e5)^(287.05 / 716.95), 100.0)
        p = max(ρ * 287.05 * T, 1.0)
        # sectional: the parcel's OWN (shrinking) crystal radius, so a crystal
        # that has sublimated down falls slower — sedimentation is then
        # size-resolved rather than pinned to the box's nominal r_eff
        rp = (rice !== nothing && rice[i] > 0f0) ? rice[i] : sp.settle_r[b]
        vt = settling_velocity(2f0 * rp, sp.settle_rho[b],
                               Float32(T), Float32(p))
        vt > 0f0 || continue
        z = zbufs[Threads.threadid()]
        get_vertical_profile!(z, m3.zlev, state.indx[i], state.indy[i], indt, nz)
        indl = clamp(state.indl[i], 1.0f0, Float32(nz))
        znew = _z_at(z, indl, nz) - Float32(Float64(vt) * Float64(dt_s))
        isfinite(znew) || continue
        state.indl[i] = _z_to_indl(z, max(znew, 0.1f0), nz)
        # cross-isentropic: rebase the prognostic θ on the level just reached,
        # else the hybrid step undoes this fall (see docstring).
        #
        # CONVENTION: θ_p must be read with the SAME map the hybrid ladder
        # inverts (isentropic.jl `_theta_at_indl`/`_theta_to_indl`).  A
        # half-level mismatch between the two (2026-08-05, the ASCENT PUMP:
        # 13–17 K on a stratospheric ladder) makes the hybrid relocate the
        # parcel after every settling step, ~0.4 km upward per step on the
        # freeze→settle→sublimate population (ht_asc4/ht_asc5).  θ is a
        # full-level field, so since the 2026-10-22 layout contract both sides
        # use the stacked read (slot k at indl k − ½); zlev above stays plain.
        if thetap !== nothing
            th = thbufs[Threads.threadid()]
            get_vertical_profile!(th, m3.theta, state.indx[i], state.indy[i], indt, nz)
            θn = _theta_at_indl(th, clamp(state.indl[i], 1.0f0, Float32(nz - 1)), nz)
            isfinite(θn) && θn > 0f0 && (thetap[i] = θn)
        end
        Threads.atomic_add!(nmov, 1)
    end
    return nmov[]
end
