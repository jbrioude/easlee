# ─────────────────────────────────────────────────────────────────────────────
# physics/plume_rise.jl — biomass-burning injection height.
#
# Two schemes, selected per release by [release.plume_rise].scheme:
#   "sofiev"  (default) — Sofiev et al. (2012) FRP-based semi-empirical plume top.
#   "freitas" (option)  — simplified Freitas-style 1-D entraining buoyant plume.
#
# Both return an injection LAYER [z_bottom, z_top] (m AGL) over which the release
# code distributes the parcels.  Pure functions — no EASLEE deps — so they unit-test
# standalone and slot into the release seeding once wired.
# ─────────────────────────────────────────────────────────────────────────────

# ── Sofiev et al. (2012) plume-top height [m AGL] ────────────────────────────
#   FRP_W  fire radiative power [W];  H_abl boundary-layer height [m];
#   N2_FT  Brunt–Väisälä² in the free troposphere [s⁻²].
# Coefficients fitted to MISR (Sofiev et al. 2012, ACP, Eqs. 10+13): α=0.24,
# β=170 m, γ=0.35, δ=0.6, Pf0=1e6 W, N0²=2.5e-4 s⁻² — verified against the paper
# (2026-07-02).  This is the one-step form; the paper's two-step variant (Eq. 19/21
# parameters for FT plumes) scores slightly worse overall and is not implemented.
function plume_top_sofiev(FRP_W::Real, H_abl::Real, N2_FT::Real;
                          α=0.24, β=170.0, γ=0.35, δ=0.6,
                          Pf0=1.0e6, N2_0=2.5e-4)
    N2 = max(float(N2_FT), 0.0)                      # only stable FT suppresses rise
    return α * H_abl + β * (max(float(FRP_W), 0.0) / Pf0)^γ * exp(-δ * N2 / N2_0)
end

# ── Mastin et al. (2009) eruption-column calibration ──────────────────────────
# H = 2.00 · V^0.241 with H the plume height [km above vent] and V the
# dense-rock-equivalent volumetric flow [m³ s⁻¹] (ρ_DRE default 2500 kg m⁻³).
# Converts an observed column height ↔ mass eruption rate when building a
# `kind="volcanic"` source term (see VOLCANIC_PLUME.md).
"""
    mastin_height_km(mer_kg_s; rho_dre=2500.0) → column height [km above vent]

Mastin et al. (2009) calibration H = 2.00·V^0.241 from the mass eruption
rate [kg/s] (V = MER/ρ_DRE, dense-rock-equivalent flow [m³/s]).
"""
mastin_height_km(mer_kg_s::Real; rho_dre::Real=2500.0) =
    2.00 * (mer_kg_s / rho_dre)^0.241
"""
    mastin_mer_kg_s(height_km; rho_dre=2500.0) → mass eruption rate [kg/s]

Inverse of `mastin_height_km`: MER from an observed column height [km].
"""
mastin_mer_kg_s(height_km::Real; rho_dre::Real=2500.0) =
    rho_dre * (height_km / 2.00)^(1 / 0.241)

# ── Volcanic eruption column (1-D steady moist plume, Woods-1988 family) ──────
#
# Top-hat entraining plume driven by the SOURCE energetics — mass eruption rate
# MER [kg/s] (total mixture: tephra + magmatic gas + vaporised external water),
# mixture temperature T0 [K] (after magma–seawater interaction) and the water
# mass fraction n_w0 (vapour at the vent).  Integrates mass/momentum/energy
# fluxes upward through the AMBIENT column (z [m ASL], T [K], p [Pa], qv),
# entraining at u_e = α·u, condensing water when saturated (latent-heat release
# — the phreatomagmatic boost that lofted Hunga Tonga far beyond a dry column
# of the same MER), until the vertical velocity vanishes:
#
#   Q = r²ρu   (mass flux /π)          dQ/dz = 2 r α u ρ_a
#   P = Q u    (momentum flux /π)      dP/dz = g r² (ρ_a − ρ)
#   F = Q(c_p T + g z + u²/2)          dF/dz = dQ/dz (c_pa T_a + g z) + L dC/dz
#
# Bulk c_p over (dry gas, vapour, condensate, solids); mixture density with
# ideal gas + incompressible solids/condensate.  Returns the TOP height (u→0),
# the neutral-buoyancy height (ρ = ρ_a crossing from below) and the suggested
# injection layer [NBH, top] — the umbrella spreads at/above NBH.
#
# Deliberately bulk: no wind bending (irrelevant for a >10⁸ kg/s column), no
# fallout, single condensate phase.  α default 0.10 (Woods).
const _CPA = 1005.0;  const _CPV = 1996.0;  const _CPL = 4186.0; const _CPS = 1100.0
const _RA  = 287.05;  const _RV  = 461.5;   const _LV  = 2.5e6;  const _G   = 9.81
# v1.1 review A5: Murphy–Koop (ice < 235 K, liquid above; ice_physics.jl
# `esat_mk`) instead of the Bolton liquid fit, which sat 1.7–1.9× above e_si
# at the 200–215 K the column, plume-phase and umbrella adjustments reach.
_esat(T) = esat_mk(T)                                            # Pa

"""
    volcanic_column_rise(mer, T0, n_w0, z_amb, T_amb, p_amb, qv_amb;
                         alpha=0.10, u0=200.0, z_vent=0.0, n_g0=0.03)
        → (z_top, z_nbh, ok)

`n_g0` is the magmatic-gas mass fraction (adds to the vent gas phase).
Heights in m ASL.  `ok=false` when the column collapses (never buoyant).
"""
function volcanic_column_rise(mer::Real, T0::Real, n_w0::Real,
                              z_amb::AbstractVector, T_amb::AbstractVector,
                              p_amb::AbstractVector, qv_amb::AbstractVector;
                              alpha::Real=0.10, u0::Real=200.0,
                              z_vent::Real=0.0, n_g0::Real=0.03)
    # ambient interpolators (linear in z)
    interp(zq, a) = begin
        k = clamp(searchsortedlast(z_amb, zq), 1, length(z_amb) - 1)
        f = clamp((zq - z_amb[k]) / max(z_amb[k+1] - z_amb[k], 1e-3), 0.0, 1.0)
        (1 - f) * a[k] + f * a[k+1]
    end
    ρ_amb(zq) = interp(zq, p_amb) / (_RA * interp(zq, T_amb))

    # vent state: gas = magmatic gas (as air) + water vapour; solids the rest
    m_v = float(n_w0); m_g = float(n_g0); m_s = max(1 - m_v - m_g, 0.0); m_l = 0.0
    T = float(T0); z = float(z_vent); u = float(u0)
    cpmix() = m_g * _CPA + m_v * _CPV + m_l * _CPL + m_s * _CPS
    function ρmix(pz)
        m_gas = m_g + m_v
        Rg = m_gas > 0 ? (m_g * _RA + m_v * _RV) / m_gas : _RA
        vol_gas = m_gas * Rg * T / pz
        return 1.0 / (vol_gas + m_s / 2500.0 + m_l / 1000.0)
    end
    p = interp(z, p_amb)
    ρ = ρmix(p)
    Q = mer / (π * 1.0)                     # flux per π: Q = r²ρu
    Qs = Q                                  # source mass flux (per π), for mixing bookkeeping
    P = Q * u
    F = Q * (cpmix() * T + _G * z + 0.5 * u^2)
    was_buoyant = false; z_nbh = NaN
    dz = 50.0
    zmax = z_amb[end] - 100
    while z < zmax
        pz  = interp(z, p_amb)
        Ta  = interp(z, T_amb)
        qa  = max(interp(z, qv_amb), 0.0)
        ρa  = pz / (_RA * Ta)
        ρ   = ρmix(pz)
        r2  = Q / (ρ * u)
        r2 <= 0 && break
        r   = sqrt(r2)
        # buoyancy bookkeeping: NBH = last upward crossing ρ < ρa → ρ ≥ ρa
        if ρ < ρa
            was_buoyant = true
        elseif was_buoyant && isnan(z_nbh)
            z_nbh = z
        end
        dQ = 2 * r * alpha * u * ρa
        dP = _G * r2 * (ρa - ρ)
        dF = dQ * (_CPA * (1 - qa) * Ta + _CPV * qa * Ta + _G * z)
        # entrained composition
        Qn = Q + dQ * dz
        w  = dQ * dz / Qn                       # entrained mass fraction this step
        m_g = (1 - w) * m_g + w * (1 - qa)
        m_v = (1 - w) * m_v + w * qa
        m_l = (1 - w) * m_l
        m_s = (1 - w) * m_s
        P += dP * dz
        F += dF * dz
        Q  = Qn
        u  = P / Q
        u <= 1.0 && break
        z += dz
        # temperature from the energy flux, then condensation adjustment
        T = (F / Q - _G * z - 0.5 * u^2) / cpmix()
        T = max(T, 150.0)
        pz = interp(z, p_amb)
        for _ in 1:3                             # saturation adjustment
            m_gas = m_g + m_v
            qsat = 0.622 * _esat(T) / max(pz - _esat(T), 1.0) * m_gas
            if m_v > qsat && T < 320.0
                dc = min(m_v - qsat, m_v)
                m_v -= dc; m_l += dc
                T += _LV * dc / cpmix()
            else
                break
            end
        end
    end
    ok = was_buoyant && z > z_vent + 2000
    isnan(z_nbh) && (z_nbh = 0.8 * z)            # top reached without ρ recrossing
    return (z, z_nbh, ok)
end

"""
    volcanic_column_provenance(mer, T0, n_w0, z_amb, T_amb, p_amb, qv_amb;
                               alpha=0.10, u0=200.0, z_vent=0.0, n_g0=0.03,
                               z_lt=5_000.0, z_tp=17_000.0)
        → (z_top, z_nbh, ok, amp, f_source, f_lower_trop, f_free_trop, f_strat)

Same column as `volcanic_column_rise`, additionally tagging every entrained
mass increment by its origin altitude: below `z_lt` (lower troposphere),
`z_lt`–`z_tp` (free troposphere), above `z_tp` (stratosphere).  `amp` is the
total mass amplification Q_top/Q_vent; the `f_*` are mass fractions OF THE AIR
DELIVERED TO THE UMBRELLA (source mixture + entrained, summing to 1).
The mass detrained at the top is the whole flux (top-hat, no lateral
detrainment), so these are the umbrella's provenance fractions.
"""
function volcanic_column_provenance(mer::Real, T0::Real, n_w0::Real,
                                    z_amb::AbstractVector, T_amb::AbstractVector,
                                    p_amb::AbstractVector, qv_amb::AbstractVector;
                                    alpha::Real=0.10, u0::Real=200.0,
                                    z_vent::Real=0.0, n_g0::Real=0.03,
                                    z_lt::Real=5_000.0, z_tp::Real=17_000.0)
    interp(zq, a) = begin
        k = clamp(searchsortedlast(z_amb, zq), 1, length(z_amb) - 1)
        f = clamp((zq - z_amb[k]) / max(z_amb[k+1] - z_amb[k], 1e-3), 0.0, 1.0)
        (1 - f) * a[k] + f * a[k+1]
    end
    m_v = float(n_w0); m_g = float(n_g0); m_s = max(1 - m_v - m_g, 0.0); m_l = 0.0
    T = float(T0); z = float(z_vent); u = float(u0)
    cpmix() = m_g * _CPA + m_v * _CPV + m_l * _CPL + m_s * _CPS
    function ρmix(pz)
        m_gas = m_g + m_v
        Rg = m_gas > 0 ? (m_g * _RA + m_v * _RV) / m_gas : _RA
        return 1.0 / (m_gas * Rg * T / pz + m_s / 2500.0 + m_l / 1000.0)
    end
    Q = mer / π
    Q0 = Q
    M_lt = 0.0; M_ft = 0.0; M_st = 0.0        # entrained mass (per π) by origin
    P = Q * u
    F = Q * (cpmix() * T + _G * z + 0.5 * u^2)
    was_buoyant = false; z_nbh = NaN
    dz = 50.0
    while z < z_amb[end] - 100
        pz = interp(z, p_amb); Ta = interp(z, T_amb); qa = max(interp(z, qv_amb), 0.0)
        ρa = pz / (_RA * Ta)
        ρ  = ρmix(pz)
        r2 = Q / (ρ * u); r2 <= 0 && break
        r  = sqrt(r2)
        if ρ < ρa
            was_buoyant = true
        elseif was_buoyant && isnan(z_nbh)
            z_nbh = z
        end
        dQ = 2 * r * alpha * u * ρa
        dm = dQ * dz
        if z - z_vent < z_lt;      M_lt += dm
        elseif z - z_vent < z_tp;  M_ft += dm
        else;                      M_st += dm
        end
        dP = _G * r2 * (ρa - ρ)
        dF = dQ * (_CPA * (1 - qa) * Ta + _CPV * qa * Ta + _G * z)
        Qn = Q + dm
        w  = dm / Qn
        m_g = (1 - w) * m_g + w * (1 - qa); m_v = (1 - w) * m_v + w * qa
        m_l = (1 - w) * m_l;                m_s = (1 - w) * m_s
        P += dP * dz; F += dF * dz; Q = Qn
        u = P / Q
        u <= 1.0 && break
        z += dz
        T = max((F / Q - _G * z - 0.5 * u^2) / cpmix(), 150.0)
        pz = interp(z, p_amb)
        for _ in 1:3
            m_gas = m_g + m_v
            qsat = 0.622 * _esat(T) / max(pz - _esat(T), 1.0) * m_gas
            if m_v > qsat && T < 320.0
                dc = min(m_v - qsat, m_v); m_v -= dc; m_l += dc
                T += _LV * dc / cpmix()
            else
                break
            end
        end
    end
    ok = was_buoyant && z > z_vent + 2000
    isnan(z_nbh) && (z_nbh = 0.8 * z)
    amp = Q / Q0
    return (z, z_nbh, ok, amp, Q0 / Q, M_lt / Q, M_ft / Q, M_st / Q)
end

"""
    volcanic_column_profiles(mer, T0, n_w0, z_amb, T_amb, p_amb, qv_amb;
                             alpha=0.10, u0=200.0, z_vent=0.0, n_g0=0.03)
        → (zgrid, lam, uprof, z_top, z_nbh, ok)

Same integration as `volcanic_column_rise`, returning the height profiles the
Lagrangian plume phase needs: `lam(z) = (1/Q)·dQ/dz` [1/m] (fractional
entrainment rate — a parcel rising dz dilutes by lam·dz) and the column rise
speed `uprof(z)`.  `zgrid` spans vent → top at the integration step.
"""
function volcanic_column_profiles(mer::Real, T0::Real, n_w0::Real,
                                  z_amb::AbstractVector, T_amb::AbstractVector,
                                  p_amb::AbstractVector, qv_amb::AbstractVector;
                                  alpha::Real=0.10, u0::Real=200.0,
                                  z_vent::Real=0.0, n_g0::Real=0.03)
    interp(zq, a) = begin
        k = clamp(searchsortedlast(z_amb, zq), 1, length(z_amb) - 1)
        f = clamp((zq - z_amb[k]) / max(z_amb[k+1] - z_amb[k], 1e-3), 0.0, 1.0)
        (1 - f) * a[k] + f * a[k+1]
    end
    m_v = float(n_w0); m_g = float(n_g0); m_s = max(1 - m_v - m_g, 0.0); m_l = 0.0
    T = float(T0); z = float(z_vent); u = float(u0)
    cpmix() = m_g * _CPA + m_v * _CPV + m_l * _CPL + m_s * _CPS
    function ρmix(pz)
        m_gas = m_g + m_v
        Rg = m_gas > 0 ? (m_g * _RA + m_v * _RV) / m_gas : _RA
        return 1.0 / (m_gas * Rg * T / pz + m_s / 2500.0 + m_l / 1000.0)
    end
    Q = mer / π
    P = Q * u
    F = Q * (cpmix() * T + _G * z + 0.5 * u^2)
    was_buoyant = false; z_nbh = NaN
    dz = 50.0
    zgrid = Float64[]; lam = Float64[]; uprof = Float64[]
    while z < z_amb[end] - 100
        pz = interp(z, p_amb); Ta = interp(z, T_amb); qa = max(interp(z, qv_amb), 0.0)
        ρa = pz / (_RA * Ta)
        ρ  = ρmix(pz)
        r2 = Q / (ρ * u); r2 <= 0 && break
        r  = sqrt(r2)
        if ρ < ρa
            was_buoyant = true
        elseif was_buoyant && isnan(z_nbh)
            z_nbh = z
        end
        dQ = 2 * r * alpha * u * ρa
        push!(zgrid, z); push!(lam, dQ / Q); push!(uprof, u)
        dP = _G * r2 * (ρa - ρ)
        dF = dQ * (_CPA * (1 - qa) * Ta + _CPV * qa * Ta + _G * z)
        Qn = Q + dQ * dz
        w  = dQ * dz / Qn
        m_g = (1 - w) * m_g + w * (1 - qa); m_v = (1 - w) * m_v + w * qa
        m_l = (1 - w) * m_l;                m_s = (1 - w) * m_s
        P += dP * dz; F += dF * dz; Q = Qn
        u = P / Q
        u <= 1.0 && break
        z += dz
        T = max((F / Q - _G * z - 0.5 * u^2) / cpmix(), 150.0)
        pz = interp(z, p_amb)
        for _ in 1:3
            m_gas = m_g + m_v
            qsat = 0.622 * _esat(T) / max(pz - _esat(T), 1.0) * m_gas
            if m_v > qsat && T < 320.0
                dc = min(m_v - qsat, m_v); m_v -= dc; m_l += dc
                T += _LV * dc / cpmix()
            else
                break
            end
        end
    end
    ok = was_buoyant && z > z_vent + 2000
    isnan(z_nbh) && (z_nbh = 0.8 * z)
    return (zgrid, lam, uprof, z, z_nbh, ok)
end

# ── PBL-top diagnostic from a θ profile (parcel method) ──────────────────────
# Lowest height where θ exceeds the near-surface θ by dθ — the mixed-layer top.
# Used when the met reader did not archive a usable boundary-layer height.
function estimate_pbl_height(zagl::AbstractVector, θcol::AbstractVector;
                             dθ = 0.5, zmin = 100.0)
    θ0 = θcol[1]
    @inbounds for k in 2:length(zagl)
        θcol[k] >= θ0 + dθ && return max(float(zagl[k]), zmin)
    end
    return max(0.3 * float(zagl[end]), zmin)
end

# ── Brunt–Väisälä² from a potential-temperature profile ──────────────────────
@inline brunt_vaisala2(θ::Real, dθdz::Real; g = 9.81) = (g / θ) * dθdz

# Free-troposphere N²: mean of N² over a slab just above the PBL top.
function n2_free_trop(zs::AbstractVector, θa::AbstractVector, H_abl::Real;
                      slab = 2000.0)
    g = 9.81; n = length(zs); acc = 0.0; cnt = 0
    for k in 1:n-1
        zc = 0.5 * (zs[k] + zs[k+1])
        (zc >= H_abl && zc <= H_abl + slab) || continue
        dθdz = (θa[k+1] - θa[k]) / (zs[k+1] - zs[k])
        acc += (g / θa[k]) * max(dθdz, 0.0); cnt += 1
    end
    return cnt > 0 ? acc / cnt : 2.5e-4              # fallback ≈ N2_0
end

# ── FRP → convective heat flux [W] (for the Freitas plume source) ────────────
# A fraction f_conv of the radiative power feeds the convective updraught.
# Total fire power ≈ FRP / f_rad (f_rad ≈ 0.10–0.15); convective ≈ f_conv of total.
frp_to_convective_heat(FRP_W::Real; f_rad = 0.13, f_conv = 0.5) =
    max(float(FRP_W), 0.0) / f_rad * f_conv

# ── Simplified Freitas-style 1-D entraining buoyant plume ────────────────────
# Steady top-hat plume (Morton–Taylor–Turner) integrated upward in the ambient
# θ(z) until the vertical velocity vanishes.  Flux state (Q,M,B):
#   Q = R²w (volume),  M = R²w² (momentum),  B = R²w·b (buoyancy), b = g(θ_p−θ_a)/θ_a
#   dQ/dz = 2 αe √M ;  dM/dz = Q·B/M ;  dB/dz = −Q·N²(z)
# Returns (z_neutral, z_top) [m AGL].  Dry — latent-heat (pyroCb) boost is a v2.
#   Q_heat fire convective heat flux [W]; fire_area [m²]; ambient zs[m] θa[K] (AGL).
function plume_layer_freitas(Q_heat::Real, fire_area::Real,
                             zs::AbstractVector, θa::AbstractVector;
                             αe = 0.1, w0 = 1.0, g = 9.81,
                             ρ = 1.0, cp = 1004.0)
    n = length(zs); n >= 2 || return (0.0, 0.0)
    R0 = sqrt(max(float(fire_area), 1.0) / π)        # source radius [m]
    θ0 = θa[1]
    # initial buoyancy from the convective heat flux: F0 = g·Q_h/(π ρ cp θ0)
    b0 = g * max(float(Q_heat), 0.0) / (π * ρ * cp * θ0) / max(R0^2 * w0, 1e-6)
    Q = R0^2 * w0
    M = R0^2 * w0^2
    B = R0^2 * w0 * b0
    θ_at(z) = begin
        z <= zs[1] && return θa[1]; z >= zs[end] && return θa[end]
        k = searchsortedlast(zs, z); t = (z - zs[k]) / (zs[k+1] - zs[k])
        θa[k] * (1 - t) + θa[k+1] * t
    end
    N2_at(z) = begin
        k = clamp(searchsortedlast(zs, z), 1, n - 1)
        brunt_vaisala2(θa[k], (θa[k+1] - θa[k]) / (zs[k+1] - zs[k]))
    end
    dz = 25.0; z = zs[1]; ztop = zs[1]; zneu = zs[1]; prevb = b0
    while z < zs[end] && M > 0.0 && Q > 0.0
        N2 = N2_at(z)
        dQ = 2 * αe * sqrt(max(M, 0.0))
        dM = Q * B / M
        dB = -Q * N2
        Q += dQ * dz; M += dM * dz; B += dB * dz
        z += dz
        b = Q > 0 ? B / Q : 0.0
        (prevb > 0 && b <= 0) && (zneu = z)          # neutral-buoyancy level
        prevb = b
        ztop = z
    end
    zneu <= zs[1] && (zneu = 0.5 * ztop)             # fallback if never crossed
    return (zneu, ztop)
end

# ── Dispatch: injection layer [z_bottom, z_top] for a release ────────────────
#   scheme "sofiev"|"freitas"; bottom "pbl"|"ground"|"fraction:<f>".
function plume_injection_layer(scheme::AbstractString, FRP_W::Real, fire_area::Real,
                               H_abl::Real, zs::AbstractVector, θa::AbstractVector;
                               bottom = "pbl")
    s = lowercase(scheme)
    if s == "freitas"
        zneu, ztop = plume_layer_freitas(frp_to_convective_heat(FRP_W), fire_area, zs, θa)
        return (zneu, ztop)
    else  # sofiev (default)
        ztop = plume_top_sofiev(FRP_W, H_abl, n2_free_trop(zs, θa, H_abl))
        # "pbl": if the plume tops out above the PBL it detrains in [PBL, top];
        #        if it stays within the PBL, the smoke fills the column [0, top].
        zb = bottom == "ground" ? 0.0 :
             startswith(bottom, "fraction:") ? parse(Float64, split(bottom, ":")[2]) * ztop :
             (ztop > H_abl ? H_abl : 0.0)
        return (min(zb, ztop), ztop)
    end
end
