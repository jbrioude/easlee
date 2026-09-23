
# ── Lagrangian volcanic plume phase (VOLCANIC_PLUME.md Level 2, steps 1+2) ────
#
# Individual trajectories RISE from the vent on their own buoyancy instead of
# being placed at a prescribed layer.  Each plume parcel carries a moist
# thermodynamic state (T, w, q_vapour, q_condensate, solids fraction) plus
# ORIGIN-FRACTION accumulators, and is sub-stepped (Δt ≈ 1 s — the sync step is
# far too coarse for a 100–250 m/s column) until it detrains (w → ~0 above the
# vent), where its final height and provenance are recorded.
#
# Entrainment (v1, "continuous dilution"): a parcel rising dz mixes with
# ambient air at the fractional rate λ(z)·dz, with λ(z) = (1/Q)dQ/dz taken from
# the 1-D column (`volcanic_column_profiles`).  Mixing is mass-weighted:
#   x ← (1−δ)x + δ·x_ambient,  δ = λ|w|Δt      (T, q_v, w; q_c/solids dilute)
# and the origin bookkeeping tags δ by the altitude band it was swallowed in
# (lower troposphere < z_lt, free troposphere < z_tp, stratosphere above), so
# the DETRAINED ensemble's provenance is an emergent quantity directly
# comparable to `volcanic_column_provenance` (same λ ⇒ the ensemble MEAN should
# agree; the heights, spread, and saturation state are the new content).
# v2 (planned) replaces the prescribed λ with pairwise exchange against
# ambient-labelled parcels through the mixing engine.
#
# Saturation adjustment (same _esat/latent-heat closure as the column model)
# keeps rising parcels at 100 % RH once condensing — the RS41 constraint.
# A Rayleigh damping τ_d kills Brunt–Väisälä ringing around the neutral level.
#
# Washout (step 3): condensate above a holding capacity q_crit is shed as
# precipitation on a timescale τ_fall (irreversible: mass leaves the parcel,
# which also UNLOADS it — a real buoyancy source).  Soluble gases partition
# into the condensate by effective Henry equilibrium
#   f_diss = HRTL/(1+HRTL),  HRTL = H_eff(T)·R·T·L,  L = q_c·ρ/ρ_w  [vol/vol]
# (an ice-retention factor applies below 250 K) and the dissolved share leaves
# with each shed increment.  Per parcel we track the SURVIVING FRACTION of the
# initial SO₂/HCl load and a full water budget (initial + entrained =
# carried + fallen) in units of the initial parcel mass.
# Washout happens in TWO phases: partially during the minutes-long column
# ascent (the parcel arrives warm — its neutral density is set by the solids
# load — so much of the water is still vapour), then dominantly in the
# umbrella as the plume cools to the ~210 K ambient and dehydrates
# (`plume_umbrella_equilibrate!`).

"""
Per-parcel plume-phase state (allocated for the full slot pool; only slots
marked `active` are integrated; `done` slots have detrained).
Origin fractions satisfy fsrc+flt+fft+fst = 1 at all times.
"""
mutable struct PlumeParcels
    active::BitVector
    done::BitVector
    z::Vector{Float64}          # height ASL [m]
    w::Vector{Float64}          # vertical velocity [m/s]
    T::Vector{Float64}          # temperature [K]
    qv::Vector{Float64}         # vapour mass fraction
    qc::Vector{Float64}         # condensate (liquid+ice) mass fraction
    msol::Vector{Float64}       # solids mass fraction
    fsrc::Vector{Float64}       # provenance fractions
    flt::Vector{Float64}
    fft::Vector{Float64}
    fst::Vector{Float64}
    ztgt::Vector{Float64}       # fountain-collapse detrainment target (NaN until drawn)
    m::Vector{Float64}          # parcel mass / initial mass (grows by entrainment,
                                # shrinks by precipitation fallout)
    went::Vector{Float64}       # entrained water   [per initial parcel mass]
    wfell::Vector{Float64}      # precipitated water [per initial parcel mass]
    fso2::Vector{Float64}       # surviving fraction of the initial SO₂ load
    fhcl::Vector{Float64}       # surviving fraction of the initial HCl load
    r::Vector{Float64}          # normalised radial position in the column
                                # (0 = axis, 1 = edge; used by v2 pairwise mode)
    fmagw::Vector{Float64}      # magmatic fraction of the parcel's SOURCE water
                                # (init-time label; 0 without a VentMix)
    # COARSE-MODE budget ([release.plume_rise].vapor_cap, 2026-08-07): vapour
    # retained above cap·q_si(T_ambient) at the umbrella equilibration —
    # physically, water that could only exist as CONDENSATE and should leave
    # as large, fast-sedimenting crystals rather than be delivered as
    # (mislabelled) supersaturated vapour.  [per initial parcel mass]
    wcoarse::Vector{Float64}
    # umbrella-equilibration shed ([per initial parcel mass]): the share of
    # `wfell` that fell AT the umbrella (not during column ascent) — the
    # first-day coarse ice cloud (Aeolus: opaque, cleared within a day),
    # emitted by the release as a fast-falling ICE_C species.
    wumb::Vector{Float64}
end

PlumeParcels(n::Integer) = PlumeParcels(
    falses(n), falses(n),
    zeros(n), zeros(n), fill(280.0, n), zeros(n), zeros(n), zeros(n),
    zeros(n), zeros(n), zeros(n), zeros(n), fill(NaN, n),
    ones(n), zeros(n), zeros(n), ones(n), ones(n), zeros(n), zeros(n),
    zeros(n), zeros(n))

"Murphy & Koop (2005) saturation vapour pressure over ICE [Pa] — the cap is
an ice-saturation criterion (the umbrella is at ~210-240 K).
v1.1 review A5/A6: single definition lives in ice_physics.jl (`esat_ice`)."
_esat_ice(T) = esat_ice(T)

"""
Ambient + column profiles the plume phase integrates through (built once per
release from the met column at the vent + the 1-D column model).
"""
struct PlumeProfiles
    z::Vector{Float64}          # ambient grid [m ASL]
    Ta::Vector{Float64}
    pa::Vector{Float64}
    qa::Vector{Float64}
    lamz::Vector{Float64}       # λ(z) on the SAME grid (0 above column top)
    ua::Vector{Float64}         # 1-D column speed u(z) (0 outside the span) —
                                # v2 edge engulfment is λ·u·dt (collective rate)
    z_vent::Float64
    z_top::Float64              # 1-D column top (diagnostic)
    z_nbh::Float64
    z_lt::Float64               # provenance band edges
    z_tp::Float64
end

"""
    build_plume_profiles(mer, T0, nw0, z, T, p, qv;
                         alpha=0.10, u0=200.0, z_vent=z[1], n_g0=0.03,
                         z_lt=5000.0, z_tp=17000.0) → PlumeProfiles
"""
function build_plume_profiles(mer::Real, T0::Real, nw0::Real,
                              z::AbstractVector, T::AbstractVector,
                              p::AbstractVector, qv::AbstractVector;
                              alpha::Real=0.10, u0::Real=200.0,
                              z_vent::Real=z[1], n_g0::Real=0.03,
                              z_lt::Real=5000.0, z_tp::Real=17000.0)
    zg, lam, up, ztop, znbh, ok =
        volcanic_column_profiles(mer, T0, nw0, z, T, p, qv;
                                 alpha=alpha, u0=u0, z_vent=z_vent, n_g0=n_g0)
    ok || @warn "plume phase: 1-D column not buoyant — λ profile may be degenerate" maxlog=1
    # resample λ and the column speed onto the ambient grid (0 outside the span)
    lamz = zeros(length(z)); ua = zeros(length(z))
    for (i, zz) in enumerate(z)
        if !isempty(zg) && zg[1] <= zz <= zg[end]
            k = clamp(searchsortedlast(zg, zz), 1, length(zg) - 1)
            f = clamp((zz - zg[k]) / max(zg[k+1] - zg[k], 1e-3), 0.0, 1.0)
            lamz[i] = (1 - f) * lam[k] + f * lam[min(k + 1, length(lam))]
            ua[i]   = (1 - f) * up[k]  + f * up[min(k + 1, length(up))]
        end
    end
    return PlumeProfiles(collect(Float64, z), collect(Float64, T),
                         collect(Float64, p), collect(Float64, qv),
                         lamz, ua, float(z_vent), ztop, znbh, float(z_lt), float(z_tp))
end

"""
In-plume washout parameters (all per-substance values are EFFECTIVE Henry
coefficients at 298 K [M/atm] with a Clausius–Clapeyron-type T dependence
`H(T) = h·exp(dhr·(1/T − 1/298))`, and an ice-retention factor applied below
250 K).  Defaults are set for a water-rich phreatomagmatic column:
- `h_so2 = 3` — SO₂ at plume pH ≈ 2 (volcanic condensate is strongly acidified
  by the co-erupted HCl/H₂SO₄, which suppresses S(IV) dissociation; the
  effective solubility collapses toward the physical Henry constant),
- `h_hcl = 2e6` — HCl dissociates completely; uptake is solubility-unlimited
  wherever condensate exists (also efficiently co-deposited on ice → r_ice=1),
- `qc_crit = 5 g/kg` — hydrometeor holding capacity of the updraft;
  `tau_fall = 240 s` — shedding (autoconversion+fallout) timescale.
"""
Base.@kwdef struct PlumeWashout
    qc_crit::Float64  = 5.0e-3
    tau_fall::Float64 = 240.0
    h_so2::Float64    = 3.0
    dhr_so2::Float64  = 3100.0
    rice_so2::Float64 = 0.5
    h_hcl::Float64    = 2.0e6
    dhr_hcl::Float64  = 9000.0
    rice_hcl::Float64 = 1.0
    # settling=true replaces the constant τ_fall with SIZE-DEPENDENT shedding:
    # an effective hydrometeor radius from the condensate load (n_ice particles
    # per kg of parcel air), a graupel-like terminal velocity
    # v_t = 12·(r/1 mm)^0.5 m/s, and shedding time τ = h_parcel/v_t.  Inside a
    # vigorous updraft hydrometeors fall into plume air below and are RECYCLED
    # — net loss only where v_t > w_gate·|w| (slow top, fountain, umbrella),
    # which is why the real column lofted its ice before dropping it.
    settling::Bool    = false
    n_ice::Float64    = 2.0e5    # hydrometeors per kg of parcel air
    h_parcel::Float64 = 500.0    # depth a hydrometeor falls to leave the parcel [m]
    w_gate::Float64   = 0.2      # shed only when v_t > w_gate·|w|
    # UMBRELLA retention (post-detrainment).  `qc_crit` above is the holding
    # capacity of an UPDRAFT; the umbrella has w ≈ 0, so what it keeps is set by
    # SEDIMENTATION over its own thickness and lifetime, not by an updraft.  With
    # a lognormal crystal spectrum (mass-weighted) the retained fraction is
    #     f = ∫ n(r)r³ exp(−v_t(r)·t/h) dr / ∫ n(r)r³ dr
    # A single radius is knife-edge (50 µm → 52 % retained, 200 µm → 0.0 %); the
    # spectrum smooths it.  `umb_r_med = 0` ⇒ disabled, fall back to `qc_crit`.
    # For HT the legacy 5 g/kg is EQUIVALENT to r ≈ 105 µm at h = 2 km, t = 1 h.
    umb_r_med::Float64 = 0.0     # mass-median crystal radius [m]
    umb_r_sg::Float64  = 2.0     # geometric standard deviation of the spectrum
    umb_h::Float64     = 2000.0  # umbrella thickness the crystals fall through [m]
    umb_t::Float64     = 3600.0  # time over which the umbrella sheds [s]
    # TEPHRA fallout during ascent.  v1 never sheds `msol`, so every parcel
    # carries its full ash load to detrainment and must arrive HOT (density
    # compensation: +40 to +185 K over ambient) — which places the detrainment
    # high and lets the water survive as vapour.  With `solid_r > 0` the coarse
    # ash sheds where its terminal velocity beats the updraft gate
    # (v_t > w_gate·|w|, falling h_parcel to leave — the same gate/geometry as
    # the hydrometeor `settling` branch).  0 ⇒ legacy (no solids fallout).
    solid_r::Float64   = 0.0     # mass-median coarse-ash radius [m]
end

"""
Mass-weighted fraction of the condensate an umbrella RETAINS: crystals of a
lognormal spectrum settle at Stokes–Cunningham speed through `umb_h` for
`umb_t`.  Evaluated at the ambient (T, p) of the detrainment level.
"""
function umbrella_retained_fraction(wash::PlumeWashout, T::Float64, p::Float64)
    wash.umb_r_med > 0.0 || return NaN
    rho = p / (_RA * T)
    mu  = 1.458e-6 * T^1.5 / (T + 110.4)
    lam = 2.0 * mu / (rho * sqrt(8.0 * _RA * T / π))
    lsg = log(wash.umb_r_sg)
    num = 0.0; den = 0.0
    lo = log(wash.umb_r_med) - 4.0 * lsg
    hi = log(wash.umb_r_med) + 4.0 * lsg
    n = 200; dl = (hi - lo) / n
    for i in 0:n                                   # trapèze en ln(r)
        lr = lo + i * dl
        r  = exp(lr)
        w  = exp(-(lr - log(wash.umb_r_med))^2 / (2.0 * lsg^2)) * r^3
        (i == 0 || i == n) && (w *= 0.5)
        kn = lam / r
        cc = 1.0 + kn * (1.257 + 0.4 * exp(-1.1 / kn))
        vt = 2.0 * 917.0 * 9.81 * r * r / (9.0 * mu) * cc
        num += w * exp(-vt * wash.umb_t / wash.umb_h)
        den += w
    end
    return den > 0.0 ? num / den : NaN
end

# effective radius [m] and terminal velocity [m/s] of the shed hydrometeors
@inline function _pp_vterm(qc::Float64, wash::PlumeWashout)
    vol = qc / (wash.n_ice * 917.0)                # ice volume per particle [m³]
    ri = cbrt(3.0 * vol / (4π))
    return min(12.0 * sqrt(ri / 1.0e-3), 30.0)
end

# dissolved fraction of a gas in the parcel condensate (effective Henry
# equilibrium; L = condensate volume per volume of plume [vol/vol])
@inline function _pp_fdiss(h298::Float64, dhr::Float64, rice::Float64,
                           T::Float64, L::Float64)
    H = h298 * exp(dhr * (1.0 / T - 1.0 / 298.15))
    hrtl = H * 0.08206 * T * L
    f = hrtl / (1.0 + hrtl)
    return T < 250.0 ? f * rice : f
end

@inline function _pp_interp(pr::PlumeProfiles, a::Vector{Float64}, zq::Float64)
    k = clamp(searchsortedlast(pr.z, zq), 1, length(pr.z) - 1)
    f = clamp((zq - pr.z[k]) / max(pr.z[k+1] - pr.z[k], 1e-3), 0.0, 1.0)
    return (1 - f) * a[k] + f * a[k+1]
end

# Deterministic per-parcel uniform in [0,1): splitmix64 finalizer of the parcel
# index.  The fountain-collapse target level ztgt is drawn AT MOST ONCE per
# parcel, so a pure hash of i gives an i.i.d.-quality draw that is reproducible
# regardless of thread count (a bare rand() inside the @threads loop made the
# plume ensemble depend on Julia's thread count).
@inline function _pp_hash_u01(i::Integer)
    x = UInt64(i) * 0x9E3779B97F4A7C15
    x ⊻= x >> 30;  x *= 0xBF58476D1CE4E5B9
    x ⊻= x >> 27;  x *= 0x94D049BB133111EB
    x ⊻= x >> 31
    return Float64(x >>> 11) * 0x1p-53
end

"""
Radial magma–seawater mixing model at the vent: the eruptive jet is
magmatic-dominated on the axis and ocean-water-dominated at the edge
(`x(r) = x_core + (x_edge − x_core)·r^p_mix`, x = seawater mass fraction of
the erupted mixture).  Per-parcel initial state comes from the
phreatomagmatic enthalpy balance

    melt·c_m·(T_mag − T) = x·[c_w(T_b − T_sea) + L_v + c_pv(T − T_b)]

solved in two regimes: full vaporisation (common T > 373 K) or, where the
water load exceeds the vaporisation cap f_max = c_m·ΔT/(c_m·ΔT + c_w·ΔT_w +
L_v) (~0.29 at T_mag = 1323 K — VOLCANIC_THEORY.md §3), a boiling-point
mixture (T = 373 K) carrying the unvaporised remainder as LIQUID condensate
qc₀ — the heavy, wet edge skirt.  Magmatic volatiles (`n_mag` of the melt)
always exsolve as vapour; their exsolution enthalpy is neglected (small).
"""
Base.@kwdef struct VentMix
    T_mag::Float64  = 1323.0     # magma temperature [K] (~1050 °C arc andesite)
    T_sea::Float64  = 298.0      # ocean surface [K]
    n_mag::Float64  = 0.03       # magmatic volatile (water) wt fraction of melt
    x_core::Float64 = 0.15       # seawater mass fraction on the axis
    x_edge::Float64 = 0.70       # … at the column edge
    p_mix::Float64  = 2.0        # radial profile exponent
    cp_magma::Float64 = 1100.0   # silicate melt heat capacity [J/kg/K]
end

"""
    vent_mix_state(vm, r) → (T0, qv0, qc0, fmag, melt)

Initial thermodynamic state of a parcel erupting at normalised radius `r`.
`fmag` = magmatic fraction of the parcel's source water, `melt` = melt mass
fraction (the SO₂/HCl carrier — gas loads scale with it).
"""
function vent_mix_state(vm::VentMix, r::Real)
    x = clamp(vm.x_core + (vm.x_edge - vm.x_core) * clamp(r, 0.0, 1.0)^vm.p_mix,
              0.0, 0.95)
    melt = 1.0 - x
    cm = vm.cp_magma; cw = _CPL; Lv = 2.26e6; Tb = 373.15
    # full-vaporisation branch: common T from the enthalpy balance
    Tfull = (melt * cm * vm.T_mag + x * (_CPV * Tb - cw * (Tb - vm.T_sea) - Lv)) /
            (melt * cm + x * _CPV)
    if Tfull >= Tb || x <= 0.0
        T0 = Tfull; fv = 1.0
    else
        # boiling-point mixture: only fv of the seawater vaporises
        T0 = Tb
        fv = clamp((melt * cm * (vm.T_mag - Tb) / x - cw * (Tb - vm.T_sea)) / Lv,
                   0.0, 1.0)
    end
    qv0 = fv * x + vm.n_mag * melt
    qc0 = (1.0 - fv) * x
    wmag = vm.n_mag * melt
    fmag = wmag / max(wmag + x, 1.0e-12)
    return (T0 = T0, qv0 = qv0, qc0 = qc0, fmag = fmag, melt = melt)
end

# area-weighted vent means (r drawn area-uniform ⇒ a = r² uniform on [0,1])
"""
    vent_mix_means(vm; nq=64) → (T0, nw0)

Area-weighted vent-mean mixture temperature [K] and total water mass fraction
of a `VentMix` radial profile (nq-point quadrature, r drawn area-uniform) —
the bulk (T0, nw0) equivalent of the radially resolved vent state.
"""
function vent_mix_means(vm::VentMix; nq::Int=64)
    T̄ = 0.0; n̄w = 0.0
    for k in 1:nq
        s = vent_mix_state(vm, sqrt((k - 0.5) / nq))
        T̄ += s.T0; n̄w += s.qv0 + s.qc0
    end
    return (T0 = T̄ / nq, nw0 = n̄w / nq)
end

"""
    plume_parcel_init!(pp, i, prof; T0, nw0, n_g0=0.03, w0=200.0,
                       ventmix=nothing, rng)

Arm slot `i` as a freshly erupted parcel at the vent (small jitter on w and T
so the ensemble develops a spread).  With `ventmix::VentMix` the initial
(T, qv, qc) come from the parcel's radial position via `vent_mix_state`
(hot dry-ish magmatic core → boiling wet edge) instead of the uniform
(T0, nw0).
"""
function plume_parcel_init!(pp::PlumeParcels, i::Int, prof::PlumeProfiles;
                            T0::Real, nw0::Real, n_g0::Real=0.03,
                            w0::Real=200.0,
                            ventmix::Union{Nothing,VentMix}=nothing,
                            rng=Random.default_rng())
    pp.active[i] = true;  pp.done[i] = false
    pp.z[i]  = prof.z_vent + 50.0 * rand(rng)
    pp.w[i]  = w0 * (0.8 + 0.4 * rand(rng))
    pp.r[i] = sqrt(rand(rng))            # area-uniform over the vent cross-section
    if ventmix === nothing
        pp.T[i]  = T0 * (0.95 + 0.10 * rand(rng))
        pp.qv[i] = nw0
        pp.qc[i] = 0.0
        pp.fmagw[i] = 0.0
    else
        s = vent_mix_state(ventmix, pp.r[i])
        pp.T[i]  = s.T0 * (0.97 + 0.06 * rand(rng))
        pp.qv[i] = s.qv0
        pp.qc[i] = s.qc0
        pp.fmagw[i] = s.fmag
    end
    pp.msol[i] = max(1.0 - pp.qv[i] - pp.qc[i] - n_g0, 0.0)
    pp.fsrc[i] = 1.0; pp.flt[i] = 0.0; pp.fft[i] = 0.0; pp.fst[i] = 0.0
    pp.ztgt[i] = NaN
    pp.m[i] = 1.0; pp.went[i] = 0.0; pp.wfell[i] = 0.0; pp.wcoarse[i] = 0.0
    pp.wumb[i] = 0.0
    pp.fso2[i] = 1.0; pp.fhcl[i] = 1.0
    return nothing
end

# mixture density (gas = air + vapour ideal mix; solids/condensate incompressible)
@inline function _pp_rho(T::Float64, qv::Float64, qc::Float64, ms::Float64, p::Float64)
    m_air = max(1.0 - qv - qc - ms, 0.0)
    m_gas = m_air + qv
    Rg = m_gas > 0 ? (m_air * _RA + qv * _RV) / m_gas : _RA
    return 1.0 / (m_gas * Rg * T / p + ms / 2500.0 + qc / 1000.0)
end
@inline _pp_cp(qv, qc, ms) = (1 - qv - qc - ms) * _CPA + qv * _CPV + qc * _CPL + ms * _CPS

"""
    plume_parcels_substep!(pp, prof, dt_total; dt_sub=1.0, tau_d=300.0,
                           w_detrain=0.5, washout=nothing) → n_detrained_this_call

Integrate every active parcel over `dt_total` seconds in `dt_sub` sub-steps:
entrain (origin-tagged), buoyancy + damping, rise, saturation-adjust, and —
when `washout::PlumeWashout` is given — shed precipitation and scavenge the
dissolved share of the parcel's SO₂/HCl load with it.  A parcel detrains when
|w| < `w_detrain` at least 2 km above the vent (it has risen and run out of
buoyancy); its slot flips to `done` with the final state preserved.
Thread-parallel over parcels (fully independent in v1).
"""
function plume_parcels_substep!(pp::PlumeParcels, prof::PlumeProfiles,
                                dt_total::Real; dt_sub::Real=1.0,
                                tau_d::Real=300.0, w_detrain::Real=0.5,
                                washout::Union{Nothing,PlumeWashout}=nothing)
    nsub = max(1, round(Int, dt_total / dt_sub))
    dt = dt_total / nsub
    fshed = washout === nothing ? 0.0 : 1.0 - exp(-dt / washout.tau_fall)
    ndet = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in eachindex(pp.active)
        pp.active[i] || continue
        detrained = false
        for _ in 1:nsub
            if _pp_step_one!(pp, i, prof, dt, tau_d, w_detrain, washout, fshed, 1.0)
                detrained = true
                break
            end
        end
        if detrained
            pp.active[i] = false; pp.done[i] = true
            Threads.atomic_add!(ndet, 1)
        end
    end
    return ndet[]
end

# One physics substep for parcel `i` (reads/writes the pp arrays in place):
# entrainment (λ scaled by `entrain_scale` — v2 gates it to the edge shell),
# buoyancy + damping + rise, saturation adjustment, washout, detrainment.
# Returns true when the parcel detrains this substep.
@inline function _pp_step_one!(pp::PlumeParcels, i::Int, prof::PlumeProfiles,
                               dt::Float64, tau_d::Real, w_detrain::Real,
                               washout::Union{Nothing,PlumeWashout},
                               fshed::Float64, entrain_scale::Float64)
    z = pp.z[i]; w = pp.w[i]; T = pp.T[i]
    qv = pp.qv[i]; qc = pp.qc[i]; ms = pp.msol[i]
    fs = pp.fsrc[i]; fl = pp.flt[i]; ff = pp.fft[i]; fx = pp.fst[i]
    ztgt = pp.ztgt[i]
    m = pp.m[i]; went = pp.went[i]; wfell = pp.wfell[i]
    fso2 = pp.fso2[i]; fhcl = pp.fhcl[i]
    detrained = false
    begin
            # ── entrain: dilution fraction δ = λ(z)·|w|·dt, origin-tagged ─────
            lam = _pp_interp(prof, prof.lamz, z) * entrain_scale
            δ = min(lam * abs(w) * dt, 0.25)
            if δ > 0
                Ta = _pp_interp(prof, prof.Ta, z)
                qa = max(_pp_interp(prof, prof.qa, z), 0.0)
                T  = (1 - δ) * T + δ * Ta
                qv = (1 - δ) * qv + δ * qa
                qc *= (1 - δ);  ms *= (1 - δ);  w *= (1 - δ)
                fs *= (1 - δ);  fl *= (1 - δ);  ff *= (1 - δ);  fx *= (1 - δ)
                zrel = z - prof.z_vent
                if zrel < prof.z_lt;      fl += δ
                elseif zrel < prof.z_tp;  ff += δ
                else;                     fx += δ
                end
                # mass bookkeeping: δ is the AMBIENT share of the final mass
                dm = m * δ / (1 - δ)
                went += qa * dm
                m += dm
            end
            # ── buoyancy + damping, rise ──────────────────────────────────────
            pz = _pp_interp(prof, prof.pa, z)
            Ta = _pp_interp(prof, prof.Ta, z)
            ρa = pz / (_RA * Ta)
            ρp = _pp_rho(T, qv, qc, ms, pz)
            w += dt * (_G * (ρa - ρp) / ρp - w / tau_d)
            znew = z + w * dt
            # adiabatic T change along the move (dry until saturation handles latent)
            T -= _G / _pp_cp(qv, qc, ms) * (znew - z)
            z = clamp(znew, prof.z_vent, prof.z[end] - 200.0)
            # ── saturation adjustment: Newton-damped condensation/evaporation ─
            # All-at-once adjustment is parity-unstable under deep
            # supersaturation (condense-all → +300 K → re-evaporate → nets to
            # nothing); the Clausius–Clapeyron-damped Newton step converges.
            pz = _pp_interp(prof, prof.pa, z)
            for _ in 1:8
                es = _esat(T)
                es >= 0.9 * pz && break              # too hot to condense at this p
                m_gas = max(1.0 - qc - ms, 0.0)
                qsat = 0.622 * es / (pz - es) * m_gas
                dq = qv - qsat
                abs(dq) < 1.0e-6 && break
                cp = _pp_cp(qv, qc, ms)
                dqs_dT = qsat * _LV / (_RV * T^2)
                dc = clamp(dq / (1 + _LV * dqs_dT / cp), -qc, qv)
                qv -= dc; qc += dc
                T += _LV * dc / cp
            end
            T = max(T, 150.0)
            # ── washout: shed condensate above the holding capacity ───────────
            # The shed increment carries the Henry-dissolved share of the
            # remaining SO₂/HCl.  Removing mass renormalises all mass fractions
            # (the parcel UNLOADS — extra buoyancy, a real effect).
            if washout !== nothing && qc > washout.qc_crit
                fsh = fshed
                if washout.settling
                    vt = _pp_vterm(qc, washout)
                    fsh = vt > washout.w_gate * abs(w) ?
                          1.0 - exp(-dt * vt / washout.h_parcel) : 0.0
                end
                dqf = (qc - washout.qc_crit) * fsh
                if dqf > 1.0e-12
                    ρp = _pp_rho(T, qv, qc, ms, pz)
                    L = qc * ρp / 1000.0                     # condensate vol/vol
                    frem = dqf / qc                          # share of condensate shed
                    fso2 *= 1 - frem * _pp_fdiss(washout.h_so2, washout.dhr_so2,
                                                 washout.rice_so2, T, L)
                    fhcl *= 1 - frem * _pp_fdiss(washout.h_hcl, washout.dhr_hcl,
                                                 washout.rice_hcl, T, L)
                    wfell += dqf * m
                    m *= 1 - dqf
                    rn = 1 / (1 - dqf)
                    qc = (qc - dqf) * rn
                    qv *= rn; ms *= rn
                end
            end
            # ── tephra fallout: coarse ash sheds where v_t beats the updraft ──
            # Same gate/geometry as the hydrometeor settling branch.  Shedding
            # solids UNLOADS the parcel (mass renormalisation), so its density
            # compensation — the +40..+185 K hot arrival — relaxes with it.
            if washout !== nothing && washout.solid_r > 0.0 && ms > 1.0e-6
                vts = min(12.0 * sqrt(washout.solid_r / 1.0e-3), 30.0)
                if vts > washout.w_gate * abs(w)
                    dms = ms * (1.0 - exp(-dt * vts / washout.h_parcel))
                    if dms > 1.0e-12
                        m *= 1 - dms
                        rn = 1 / (1 - dms)
                        ms = (ms - dms) * rn
                        qv *= rn; qc *= rn
                    end
                end
            end
            # ── detrainment ───────────────────────────────────────────────────
            # (a) SLOW and NO LONGER BUOYANT above the vent (a bare |w|≈0 check
            #     fires at every ballistic apex; a buoyant parcel at w≈0 must
            #     re-accelerate, not detrain);
            # (b) fountain collapse: a parcel DESCENDING through the umbrella
            #     band [NBH, top] joins the intrusion there — individual parcels
            #     otherwise overshoot to their personal (undiluted) neutral
            #     level far above the collective column top.
            if abs(w) < w_detrain && z > prof.z_vent + 2000.0
                pz2 = _pp_interp(prof, prof.pa, z)
                ρa2 = pz2 / (_RA * _pp_interp(prof, prof.Ta, z))
                if _pp_rho(T, qv, qc, ms, pz2) >= 0.998 * ρa2
                    if z <= prof.z_top * 1.02
                        detrained = true
                    else
                        # low-dilution parcel equilibrated ABOVE the collective
                        # top: parameterised fountain mixing folds it into the
                        # umbrella band (draw a target level once) — warm it
                        # adiabatically along the fold and re-equilibrate water
                        isnan(ztgt) && (ztgt = prof.z_nbh + _pp_hash_u01(i) * (prof.z_top - prof.z_nbh))
                        T += _G / _pp_cp(qv, qc, ms) * (z - ztgt)
                        z = ztgt
                        pzf = _pp_interp(prof, prof.pa, z)
                        for _ in 1:8
                            esf = _esat(T)
                            esf >= 0.9 * pzf && break
                            m_gasf = max(1.0 - qc - ms, 0.0)
                            qsf = 0.622 * esf / (pzf - esf) * m_gasf
                            dqf = qv - qsf
                            abs(dqf) < 1.0e-6 && break
                            cpf = _pp_cp(qv, qc, ms)
                            dcf = clamp(dqf / (1 + _LV * qsf * _LV / (_RV * T^2) / cpf), -qc, qv)
                            qv -= dcf; qc += dcf
                            T += _LV * dcf / cpf
                        end
                        detrained = true
                    end
                end
            end
            # fountain collapse: descending material spreads THROUGH the
            # umbrella band — draw a per-parcel detrainment level in
            # [NBH, top] on the first downward crossing, detrain on reaching it
            if !detrained && w < -w_detrain && z <= prof.z_top
                isnan(ztgt) && (ztgt = prof.z_nbh + _pp_hash_u01(i) * (prof.z_top - prof.z_nbh))
                z <= ztgt && (detrained = true)
            end
    end
    pp.z[i] = z; pp.w[i] = w; pp.T[i] = T
    pp.qv[i] = qv; pp.qc[i] = qc; pp.msol[i] = ms
    pp.fsrc[i] = fs; pp.flt[i] = fl; pp.fft[i] = ff; pp.fst[i] = fx
    pp.ztgt[i] = ztgt
    pp.m[i] = m; pp.went[i] = went; pp.wfell[i] = wfell
    pp.fso2[i] = fso2; pp.fhcl[i] = fhcl
    return detrained
end

# ── v2: pairwise (θ, q, w) exchange between parcels ──────────────────────────
#
# Replaces the mean-field entrainment with the "Dirac in a crowd"
# configuration: parcels carry a normalised radial position r ∈ [0,1] in the
# column; only the OUTER SHELL (r > redge) engulfs raw ambient air — at a rate
# λ(z)/f_shell so the column-mean dilution still equals λ(z) — and interior
# parcels feel the environment only secondhand, through mass-weighted
# conservative pairwise exchange with r-adjacent neighbours in their height
# bin (the same F_ij ∝ (x_j − x_i) form as the puff-mixing engine).  This
# produces the radial structure the mean-field model cannot: a hot,
# source-rich core and a cool, diluted edge.

# Conservative pairwise exchange over all active parcels: bin by height, sort
# each bin by r, pair r-ADJACENT parcels (random parity so partners rotate),
# and relax each pair toward its mass-weighted mean by μ = 1 − exp(−dt/τ_mix).
# Exchanges the intensive state (T, q_v, q_c, solids, w, provenance) and the
# per-initial-load gas trackers via their per-mass concentration (conserves
# the total gas exactly).
function _plume_pairwise_pass!(pp::PlumeParcels, dt::Real;
                               tau_mix::Real=60.0, dzbin::Real=1000.0,
                               rng=Random.default_rng())
    μ = 1.0 - exp(-dt / tau_mix)
    idx = findall(pp.active)
    length(idx) < 2 && return nothing
    # bin by height
    zmin = minimum(view(pp.z, idx))
    bins = Dict{Int,Vector{Int}}()
    for i in idx
        push!(get!(bins, floor(Int, (pp.z[i] - zmin) / dzbin), Int[]), i)
    end
    for (_, b) in bins
        nb = length(b)
        nb < 2 && continue
        sort!(b; by=i -> pp.r[i])
        k0 = 1 + (rand(rng) < 0.5 ? 0 : 1)          # random parity
        for k in k0:2:(nb - 1)
            i = b[k]; j = b[k + 1]
            mi = pp.m[i]; mj = pp.m[j]; s = mi + mj
            wi = mj / s; wj = mi / s                # push i by mj/s, j by mi/s
            @inline function _xch!(a::Vector{Float64})
                Δ = μ * (a[j] - a[i])
                a[i] += Δ * wi; a[j] -= Δ * wj
            end
            _xch!(pp.T); _xch!(pp.qv); _xch!(pp.qc); _xch!(pp.msol); _xch!(pp.w)
            _xch!(pp.fsrc); _xch!(pp.flt); _xch!(pp.fft); _xch!(pp.fst)
            # gas loads: exchange the per-mass concentration c = f/m
            for a in (pp.fso2, pp.fhcl)
                ci = a[i] / mi; cj = a[j] / mj
                Δ = μ * (cj - ci)
                a[i] = (ci + Δ * wi) * mi; a[j] = (cj - Δ * wj) * mj
            end
        end
    end
    return nothing
end

# Turbulent radial wander with homogenisation time τ_r.  The walk is in
# NORMALISED AREA a = r² (reflecting at axis and edge): its stationary density
# is uniform in a, i.e. area-uniform in the cross-section — the invariant of
# incompressible mixing, and the density the shell fraction f_shell = 1−r_e²
# normalisation assumes.  (A walk in r itself drifts the population toward
# uniform-in-r, draining the shell and under-delivering entrainment.)
function _plume_radial_walk!(pp::PlumeParcels, dt::Real; tau_r::Real=240.0,
                             rng=Random.default_rng())
    σ = sqrt(dt / tau_r)
    for i in eachindex(pp.active)
        pp.active[i] || continue
        a = pp.r[i]^2 + σ * randn(rng)
        a = abs(a); a > 1.0 && (a = 2.0 - a)
        pp.r[i] = sqrt(clamp(a, 0.0, 1.0))
    end
    return nothing
end

"""
    plume_parcels_substep_v2!(pp, prof, dt_total; dt_sub=1.0, tau_d=300.0,
                              w_detrain=0.5, washout=nothing, redge=0.75,
                              tau_mix=60.0, tau_r=240.0, dzbin=1000.0)
        → n_detrained_this_call

v2 of the plume phase: pairwise-exchange entrainment (see above).  Same
per-parcel physics as `plume_parcels_substep!`, but ambient engulfment is
confined to the outer shell (r > `redge`, rate λ/f_shell) and parcels mix
with each other pairwise (τ_mix) while wandering radially (τ_r).  The
column-MEAN dilution matches the 1-D λ(z), so ensemble provenance and
detrainment heights should agree with v1 — the radial structure is the new
content.
"""
function plume_parcels_substep_v2!(pp::PlumeParcels, prof::PlumeProfiles,
                                   dt_total::Real; dt_sub::Real=1.0,
                                   tau_d::Real=300.0, w_detrain::Real=0.5,
                                   washout::Union{Nothing,PlumeWashout}=nothing,
                                   redge::Real=0.75, tau_mix::Real=60.0,
                                   tau_r::Real=240.0, dzbin::Real=1000.0)
    nsub = max(1, round(Int, dt_total / dt_sub))
    dt = dt_total / nsub
    fshed = washout === nothing ? 0.0 : 1.0 - exp(-dt / washout.tau_fall)
    fshell = 1.0 - redge^2
    ndet = 0
    for _ in 1:nsub
        _plume_pairwise_pass!(pp, dt; tau_mix=tau_mix, dzbin=dzbin)
        _plume_radial_walk!(pp, dt; tau_r=tau_r)
        ndetA = Threads.Atomic{Int}(0)
        Threads.@threads :static for i in eachindex(pp.active)
            pp.active[i] || continue
            # Edge engulfment at the COLLECTIVE rate λ(z)·u_col(z)/f_shell:
            # the 1-D λ is defined per unit rise OF THE COLUMN, and a cooled,
            # slowed edge parcel keeps being scraped by the shear interface —
            # using its own |w| would under-dilute the ensemble.
            # (_pp_step_one! computes δ = λ·scale·|w|·dt, so scale carries u/|w|.)
            esc = 0.0
            if pp.r[i] > redge
                ucol = _pp_interp(prof, prof.ua, pp.z[i])
                esc = ucol / (max(abs(pp.w[i]), 0.5) * fshell)
            end
            if _pp_step_one!(pp, i, prof, dt, tau_d, w_detrain, washout, fshed, esc)
                pp.active[i] = false; pp.done[i] = true
                Threads.atomic_add!(ndetA, 1)
            end
        end
        ndet += ndetA[]
    end
    return ndet
end

"""
    plume_umbrella_equilibrate!(pp, prof; washout, q_ret=washout.qc_crit) → n

Post-detrainment umbrella step (the first hours after injection): each DONE
parcel relaxes to the AMBIENT temperature at its detrainment level (the
umbrella mixes with stratospheric air and radiates — the environment is the
heat bath, so the condensation enthalpy is exported rather than warming the
parcel), condenses down to saturation, and sheds condensate above the retained
amount `q_ret` — carrying the Henry-dissolved share of the remaining SO₂/HCl
with it (ice retention applies: the umbrella is at ~210 K).

This is where MOST of the vaporised ocean water is removed: the column ascent
is only minutes and arrives warm (its neutral density is set by the solids
load), so in-column washout is partial; the RS41-observed saturated, ice-rich
umbrella then dehydrates toward `q_sat(T_ambient)` + retained ice.  Splitting
the fallout between the two phases is exactly the "washed out right at the
beginning vs surviving in the stratosphere" budget.
"""
function plume_umbrella_equilibrate!(pp::PlumeParcels, prof::PlumeProfiles;
                                     washout::Union{Nothing,PlumeWashout}=nothing,
                                     q_ret::Real=washout === nothing ? 0.0 : washout.qc_crit,
                                     vapor_cap::Real=Inf)
    n = 0
    for i in eachindex(pp.done)
        pp.done[i] || continue
        z = pp.z[i]
        Ta = _pp_interp(prof, prof.Ta, z)
        pz = _pp_interp(prof, prof.pa, z)
        T = Ta
        qv = pp.qv[i]; qc = pp.qc[i]; ms = pp.msol[i]; m = pp.m[i]
        # isothermal saturation adjustment at the ambient temperature
        es = _esat(T)
        if es < 0.9 * pz
            m_gas = max(1.0 - qc - ms, 0.0)
            qsat = 0.622 * es / (pz - es) * m_gas
            dc = clamp(qv - qsat, -qc, qv)
            qv -= dc; qc += dc
        end
        # shed: either down to a fixed holding capacity (legacy `q_ret`) or by
        # SEDIMENTATION of the crystal spectrum (physical — see PlumeWashout)
        qkeep = q_ret
        if washout !== nothing && washout.umb_r_med > 0.0
            fr = umbrella_retained_fraction(washout, T, pz)
            isnan(fr) || (qkeep = qc * fr)
        end
        if washout !== nothing && qc > qkeep
            dqf = qc - qkeep
            ρg = pz / (_RA * T)
            L = qc * ρg / 1000.0
            frem = dqf / qc
            pp.fso2[i] *= 1 - frem * _pp_fdiss(washout.h_so2, washout.dhr_so2,
                                               washout.rice_so2, T, L)
            pp.fhcl[i] *= 1 - frem * _pp_fdiss(washout.h_hcl, washout.dhr_hcl,
                                               washout.rice_hcl, T, L)
            pp.wfell[i] += dqf * m
            pp.wumb[i]  += dqf * m       # umbrella share -> ICE_C emission
            m *= 1 - dqf
            r = 1 / (1 - dqf)
            qc = (qc - dqf) * r
            qv *= r; ms *= r
        end
        # ── vapour cap ([release.plume_rise].vapor_cap, 2026-08-07) ──────────
        # The saturation adjustment above uses LIQUID e_s at the ambient T, so
        # at 228–240 K a parcel legally RETAINS vapour far above what it can
        # hold against the transport model's ICE gates once placed — that
        # vapour flash-freezes after handoff (41–53 Tg in the first hours of
        # every run) instead of leaving as the observed fast-sedimenting
        # coarse ice (Aeolus: opaque cloud removed within day 1).  Cap the
        # retained vapour at vapor_cap·q_si(T_ambient); the excess goes to
        # `wcoarse`, the coarse-mode budget the release emits as a large-
        # crystal species.  Inf = off (bit-identical legacy).  Bookkeeping
        # identical to the condensate shed above.
        if isfinite(vapor_cap)
            esi = _esat_ice(T)
            if esi < 0.9 * pz
                m_gas = max(1.0 - qc - ms, 0.0)
                qsi = 0.622 * esi / (pz - esi) * m_gas
                cap = vapor_cap * qsi
                if qv > cap && cap < 1.0
                    # shed enough that the POST-renormalisation vapour fraction
                    # equals the cap: (qv - dqv)/(1 - dqv) = cap
                    dqv = (qv - cap) / (1 - cap)
                    pp.wcoarse[i] += dqv * m
                    m *= 1 - dqv
                    r = 1 / (1 - dqv)
                    qv = (qv - dqv) * r          # = cap exactly
                    qc *= r; ms *= r
                end
            end
        end
        pp.T[i] = T; pp.qv[i] = qv; pp.qc[i] = qc; pp.msol[i] = ms; pp.m[i] = m
        n += 1
    end
    return n
end

"""
    plume_release_ensemble(mer, T0, nw0, z, T, p, qv; alpha=0.10, n=2000,
                           washout=PlumeWashout(), max_time=3600.0)
        → (z_det, fso2, fhcl, fh2o, prof)

Run a representative plume-phase ensemble through the given ambient column and
return the DETRAINMENT distributions the release machinery samples from:
`z_det` [m ASL on the input z grid], and per-member surviving fractions of the
initial SO₂ / HCl / water load (`fh2o` = carried water / lofted water).
`washout=nothing` disables scavenging (all fractions 1).  Stragglers that have
not detrained by `max_time` are folded at their current height.
"""
function plume_release_ensemble(mer::Real, T0::Real, nw0::Real,
                                z::AbstractVector, T::AbstractVector,
                                p::AbstractVector, qv::AbstractVector;
                                alpha::Real=0.10, n::Int=2000,
                                washout::Union{Nothing,PlumeWashout}=PlumeWashout(),
                                entrainment::Symbol=:meanfield,   # :meanfield | :pairwise
                                ventmix::Union{Nothing,VentMix}=nothing,
                                max_time::Real=3600.0,
                                vapor_cap::Real=Inf,
                                rng=Random.default_rng())
    # With a vent mix the 1-D column (λ, u, z_top/NBH references) runs on the
    # AXIS (core) state: the ascending jet IS the magmatic core — the wet
    # boiling edge (beyond the vaporisation cap) is negatively buoyant and
    # collapses into the low skirt, which the 1-D model does not describe.
    # (The area-MEAN mixture sits at 373 K, water-loaded, and cannot rise —
    # driving the profiles with it degenerates λ and the detrainment bands.)
    if ventmix !== nothing
        s0 = vent_mix_state(ventmix, 0.0)
        T0 = s0.T0; nw0 = s0.qv0
    end
    prof = build_plume_profiles(mer, T0, nw0, z, T, p, qv; alpha=alpha)
    pp = PlumeParcels(n)
    for i in 1:n
        plume_parcel_init!(pp, i, prof; T0=T0, nw0=nw0, ventmix=ventmix, rng=rng)
    end
    # per-parcel lofted source water + melt (gas-load) weight, captured at init
    w0src = [pp.qv[i] + pp.qc[i] for i in 1:n]
    gasw  = ones(n)
    if ventmix !== nothing
        melts = [vent_mix_state(ventmix, pp.r[i]).melt for i in 1:n]
        m̄ = sum(melts) / n
        m̄ > 0 && (gasw = melts ./ m̄)
    end
    t = 0.0
    while any(pp.active) && t < max_time
        entrainment === :pairwise ?
            plume_parcels_substep_v2!(pp, prof, 30.0; washout=washout) :
            plume_parcels_substep!(pp, prof, 30.0; washout=washout)
        t += 30.0
    end
    for i in 1:n                     # stragglers: freeze where they are
        pp.active[i] && (pp.active[i] = false; pp.done[i] = true)
    end
    washout === nothing ||
        plume_umbrella_equilibrate!(pp, prof; washout=washout,
                                    vapor_cap=vapor_cap)
    fh2o = [(pp.qv[i] + pp.qc[i]) * pp.m[i] / max(w0src[i] + pp.went[i], 1e-12)
            for i in 1:n]
    # coarse-mode budget: same normalisation as fh2o (fraction of the member's
    # lofted+entrained water), so fh2o + fcoarse + wfell-share ≈ 1
    fcoarse = [pp.wcoarse[i] / max(w0src[i] + pp.went[i], 1e-12) for i in 1:n]
    fumb    = [pp.wumb[i]    / max(w0src[i] + pp.went[i], 1e-12) for i in 1:n]
    return (z_det = copy(pp.z), fso2 = copy(pp.fso2), fhcl = copy(pp.fhcl),
            fh2o = fh2o, fcoarse = fcoarse, fumb = fumb,
            fmagw = copy(pp.fmagw), gasw = gasw, w0src = w0src,
            qdet = copy(pp.qv),   # vapour fraction at detrainment [kg/kg] — ice closure q₀
            prof = prof)
end

"""
    plume_phase_stats(pp) → NamedTuple

Ensemble diagnostics over the DETRAINED parcels: detrainment-height
distribution and mean provenance — the quantities to compare against the
1-D `volcanic_column_provenance` and the observed cloud layer.
"""
function plume_phase_stats(pp::PlumeParcels)
    idx = findall(pp.done)
    isempty(idx) && return (n=0,)
    zs = pp.z[idx]
    q(p) = partialsort(copy(zs), max(1, round(Int, p * length(zs))))
    (n = length(idx),
     z_mean = sum(zs) / length(zs),
     z_p10 = q(0.10), z_p50 = q(0.50), z_p90 = q(0.90),
     fsrc = sum(pp.fsrc[idx]) / length(idx),
     flt  = sum(pp.flt[idx])  / length(idx),
     fft  = sum(pp.fft[idx])  / length(idx),
     fst  = sum(pp.fst[idx])  / length(idx),
     qv_mean = sum(pp.qv[idx]) / length(idx),
     qc_mean = sum(pp.qc[idx]) / length(idx),
     # ── washout budget (units of the initial parcel mass, ensemble means) ──
     w_carried = sum((pp.qv[i] + pp.qc[i]) * pp.m[i] for i in idx) / length(idx),
     w_entrained = sum(pp.went[idx]) / length(idx),
     w_fallen  = sum(pp.wfell[idx]) / length(idx),
     fso2 = sum(pp.fso2[idx]) / length(idx),
     fhcl = sum(pp.fhcl[idx]) / length(idx),
     fmagw = sum(pp.fmagw[idx]) / length(idx))
end
