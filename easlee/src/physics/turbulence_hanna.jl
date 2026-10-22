
# ── Hanna (1982) / FLEXPART boundary-layer turbulence parameterization ────────
#
# Used when 3D TKE is not available (e.g. ERA5).  Two operating modes:
#
#   Option a (pre-compute): hanna_fill_tke! populates the TKE and TLS fields
#     during met initialisation; the standard Thomson/Langevin path then uses
#     those fields without any code-path change in the physics kernel.
#
#   Option b (on-the-fly): turbulence_hanna_v1 computes σ_w, dσ_w/dz and T_L
#     at every Langevin sub-step from 2D surface parameters (ustar, pblh,
#     wstar, ol) interpolated at the parcel's horizontal position.
#
# Stability regimes (following FLEXPART / Hanna 1982):
#   Neutral    : pblh / |L| < 1      — exponential decay with height
#   Convective : L < 0, pblh/|L| ≥ 1 — Lenschow / Hanna CBL profile
#   Stable     : L > 0, pblh/|L| ≥ 1 — linear decrease with height
#
# References:
#   Hanna S.R. (1982) Applications in Air Pollution Modelling, in
#     Nieuwstadt & van Dop (eds), Atmos. Turbulence and Air Pollution Modelling.
#   Stohl et al. (2005) Atmos. Chem. Phys. 5, 2461-2474 (FLEXPART v6 paper).

# ── Wstar from Monin-Obukhov scaling ─────────────────────────────────────────
# w* = u* · (-h/L)^(1/3)  for L < 0 (unstable); 0 otherwise.
@inline function hanna_wstar(ustar::Float32, pblh::Float32, ol::Float32)::Float32
    ol < 0.0f0 || return 0.0f0
    return ustar * cbrt(max(-pblh / ol, 0.0f0))
end

# ── Core Hanna (1982) / FLEXPART BL parameterization ─────────────────────────
#
# Returns all turbulence quantities needed by the Langevin equation:
#   sigu, sigv : horizontal velocity std-devs (m/s)  — for future use
#   sigw       : vertical velocity std-dev (m/s)
#   dsigwdz    : ∂σ_w/∂z (m/s per m) — analytically computed, used in drift term
#   tlu, tlv   : horizontal Lagrangian T_L (s)       — for future use
#   tlw        : vertical  Lagrangian T_L (s)
#
# Input z is height AGL.  Below pblh the BL parameterization applies;
# outside this function's scope (z ≥ pblh) the caller applies free-tropo values.
#
@inline function hanna_bl(
    z    ::Float32,
    pblh ::Float32,
    ol   ::Float32,
    ust  ::Float32,
    wst  ::Float32,
)::NTuple{7,Float32}
    sigu::Float32    = 0.0f0
    sigv::Float32    = 0.0f0
    sigw::Float32    = 0.0f0
    dsigwdz::Float32 = 0.0f0
    tlu::Float32     = 0.0f0
    tlv::Float32     = 0.0f0
    tlw::Float32     = 0.0f0

    # Non-finite Obukhov length (failed heat-flux derivation) must classify as
    # NEUTRAL, not stable: NaN compares false everywhere, so `ol < 0` below
    # would silently route NaN-L columns into the stable branch (σw decays with
    # height, convective turbulence killed).  aol = Inf → h/aol = 0 → neutral.
    aol = isfinite(ol) ? abs(ol) : Inf32
    h   = max(pblh, 50.0f0)

    # ── Neutral: pblh / |L| < 1 ───────────────────────────────────────────────
    if aol < 1.0f-6 || h / aol < 1.0f0
        ust   = max(ust, 1.0f-4)
        corr  = z / ust           # z / u* [s]  (FLEXPART convention)
        sigu  = 1.0f-2 + 2.0f0 * ust * exp(-3.0f-4 * corr)
        sigw0 = 1.3f0 * ust * exp(-2.0f-4 * corr)
        # dsigwdz before the floor term — derivative of 1.3·u*·exp(-2e-4·z/u*)
        dsigwdz = -2.0f-4 / ust * sigw0         # ∂/∂z [1.3·u*·exp(-2e-4·z/u*)]
        sigw    = sigw0 + 1.0f-2
        sigv    = sigw
        tlu_raw = 0.5f0 * z / sigw / max(1.0f0 + 1.5f-3 * corr, 1.0f-6)
        tlu  = clamp(tlu_raw, 0.5f0, 300.0f0)
        tlv  = tlu
        tlw  = tlu

    # ── Unstable / convective: L < 0, pblh/|L| ≥ 1 ───────────────────────────
    elseif ol < 0.0f0
        zeta = max(z / h, 1.0f-6)

        sigu = 1.0f-2 + ust * cbrt(12.0f0 - 0.5f0 * h / ol)
        sigv = sigu

        # Lenschow / Hanna (1982) σ_w formula
        inner = max(1.2f0 * wst^2 * (1.0f0 - 0.9f0 * zeta) * zeta^(2.0f0/3.0f0) +
                    (1.8f0 - 1.4f0 * zeta) * ust^2, 0.0f0)
        sigw_bl = sqrt(inner)

        # Analytical ∂σ_w/∂z (chain rule on inner)
        # d(inner)/dz = (1/h) · d(inner)/d(zeta)
        # d/d(zeta) [1.2·w*²·(1-0.9z)·z^(2/3) + (1.8-1.4z)·u*²]
        #  = 1.2·w*²·(-0.9·zeta^(2/3) + (1-0.9z)·(2/3)·zeta^(-1/3)) - 1.4·u*²
        d_inner_dzeta = 1.2f0 * wst^2 *
                         (-0.9f0 * zeta^(2.0f0/3.0f0) +
                          (1.0f0 - 0.9f0 * zeta) * (2.0f0/3.0f0) * max(zeta, 1.0f-9)^(-1.0f0/3.0f0)) -
                        1.4f0 * ust^2
        dsigwdz = sigw_bl > 1.0f-6 ? 0.5f0 * d_inner_dzeta / (h * sigw_bl) : 0.0f0

        sigw = sigw_bl + 1.0f-2

        # Horizontal TL (Hanna 1982)
        tlu = clamp(0.15f0 * h / sigu, 0.5f0, 300.0f0)
        tlv = tlu

        # Vertical TL — three sub-regimes (Hanna 1982 / FLEXPART)
        if z < aol
            tlw = 0.1f0 * z / max(sigw * (0.55f0 - 0.38f0 * z / aol), 1.0f-6)
        elseif zeta < 0.1f0
            tlw = 0.59f0 * z / max(sigw, 1.0f-6)
        else
            tlw = 0.15f0 * h / max(sigw, 1.0f-6) * (1.0f0 - exp(-5.0f0 * zeta))
        end
        tlw = clamp(tlw, 0.5f0, 300.0f0)

    # ── Stable: L > 0, pblh/|L| ≥ 1 ─────────────────────────────────────────
    else
        zoh  = min(z / h, 1.0f0)
        ust  = max(ust, 1.0f-4)
        sigw = 1.3f0 * ust * max(1.0f0 - zoh, 0.05f0) + 1.0f-2
        dsigwdz = -1.3f0 * ust / h   # constant slope
        sigu = 2.2f0 * ust * max(1.0f0 - zoh, 0.05f0) + 1.0f-2
        sigv = sigu
        tlu  = clamp(0.07f0 * h / sigu, 0.5f0, 300.0f0)
        tlv  = tlu
        tlw  = clamp(0.1f0 * h / max(sigw, 1.0f-6), 0.5f0, 300.0f0)
    end

    # σ_w floor — and where the floor binds, the TRUE local gradient is zero:
    # keeping the analytic (negative) dσ_w/dz there injects a spurious
    # downward drift into the Langevin equation and piles parcels up at the
    # ground (well-mixed violation, worst in the stable regime).
    if sigw < 0.2f0
        sigw    = 0.2f0
        dsigwdz = 0.0f0
    end
    return sigu, sigv, sigw, dsigwdz, tlu, tlv, tlw
end

# ── Full Hanna profile at height z (BL + entrainment zone + free tropo) ──────
#
# The BL and free-troposphere regimes are joined by a finite ENTRAINMENT ZONE
# [EZ_LO·h, EZ_HI·h] instead of a hard switch at z = h:
#   • a discontinuous σ_w² at h is exactly where Langevin schemes violate the
#     well-mixed criterion (the drift term ∂σ_w²/∂z becomes a delta function
#     the integrator cannot see) — FLEXPART avoids this by reflecting parcels
#     at h, closing the PBL;
#   • the blend resolves the gradient over a finite inversion-layer depth, so
#     trajectories cross the PBL top with a controlled, physical exchange rate
#     and the drift correction stays consistent with the σ_w profile.
# Variances (not σ_w) are blended — the energetically consistent choice — with
# a cosine weight; dσ_w/dz includes the analytic weight-gradient term.
const EZ_LO = 0.85f0    # entrainment zone bottom, fraction of pblh
const EZ_HI = 1.15f0    # entrainment zone top

@inline function hanna_sigw_tlw(
    z      ::Float32,
    pblh   ::Float32,
    ustar  ::Float32,
    wstar  ::Float32,
    ol     ::Float32,
    sigw_ft::Float32,
)::Tuple{Float32,Float32,Float32}   # (sigw, dsigwdz, tlw)
    h  = max(pblh, 50.0f0)
    z1 = EZ_LO * h
    z2 = EZ_HI * h
    sigw_f = max(sigw_ft, 0.05f0)
    tlw_f  = clamp(100.0f0 * sigw_f, 0.5f0, 300.0f0)
    if z >= z2
        return sigw_f, 0.0f0, tlw_f
    elseif z <= z1
        _, _, sigw, dsigwdz, _, _, tlw = hanna_bl(z, pblh, ol, ustar, wstar)
        return sigw, dsigwdz, tlw
    end
    # Entrainment zone: Hanna formulas remain well-defined up to EZ_HI·h
    # (the convective inner term is floored at 0), so evaluate them at z.
    _, _, sigw_b, dsig_b, _, _, tlw_b = hanna_bl(z, pblh, ol, ustar, wstar)
    x    = (z - z1) / (z2 - z1)
    w    = 0.5f0 * (1.0f0 + cos(Float32(π) * x))                 # 1 → 0
    dwdz = -0.5f0 * Float32(π) / (z2 - z1) * sin(Float32(π) * x)
    var  = w * sigw_b^2 + (1.0f0 - w) * sigw_f^2
    sigw = sqrt(max(var, 0.25f0 * sigw_f^2))
    dsigwdz = (w * 2.0f0 * sigw_b * dsig_b + dwdz * (sigw_b^2 - sigw_f^2)) /
              (2.0f0 * max(sigw, 1.0f-3))
    tlw  = clamp(w * tlw_b + (1.0f0 - w) * tlw_f, 0.5f0, 300.0f0)
    return max(sigw, 0.05f0), dsigwdz, tlw
end

# ── Option a: pre-compute TKE field from Hanna ───────────────────────────────
#
# Called from metgribfieldcalculation! once per time slice.
# Overwrites metfields3D["tke"][:,:,:,indt] with 0.5·σ_w²  and
#             metfields3D["tls"][:,:,:,indt] with ℓ = 0.5·T_L·σ_w.
#
function hanna_fill_tke!(metfields3D, metfields2D, indt::Int,
                         Ni::Int, Nj::Int, Nk::Int;
                         sigw_ft::Float32 = 0.1f0)
    for j in 1:Nj
        for i in 1:Ni
            ustar_ij = Float32(metfields2D["ustar"][i, j, indt])
            pblh_ij  = Float32(metfields2D["pblh"][i, j, indt])
            ol_ij    = Float32(metfields2D["ol"][i, j, indt])
            wstar_ij = hanna_wstar(ustar_ij, pblh_ij, ol_ij)

            for k in 1:Nk
                z = Float32(metfields3D["zlev"][i, j, k, indt])
                sigw, _, tlw = hanna_sigw_tlw(z, pblh_ij, ustar_ij,
                                               wstar_ij, ol_ij, sigw_ft)
                metfields3D["tke"][i, j, k, indt] = 0.5f0 * sigw^2
                metfields3D["tls"][i, j, k, indt] = 0.5f0 * tlw * sigw
            end
        end
    end
    return nothing
end

# ── Archived-K_z fill (ERA5 avg_tdch) ────────────────────────────────────────
#
# Some drivers archive the host model's OWN vertical diffusion coefficient.
# ERA5 has `avg_tdch` (paramId 235014, "time-mean turbulent diffusion
# coefficient for heat", m²/s) on every model level, already inside the
# IFS_EA_MET files — it was simply never read.  Heat is the right coefficient
# for a passive tracer (turbulent Prandtl number ≈ 1); momentum would not be.
#
# SPLITTING K INTO (σ_w, T_L).  A diffusivity is one number and the Langevin
# step needs two, K = σ_w²·T_L, so a second relation is required.  Three were
# compared on a neutral column (u* = 0.227 m/s, ERA5 K):
#
#   z [m]   K      ℓ=BL89 geometric      σ_w from similarity     Hanna alone
#                  σ_w      T_L          σ_w     T_L=K/σ_w²      σ_w    T_L
#     100  14.4   0.288     174          0.280      183          0.280   107
#     325  13.9   0.086    1900          0.232      259          0.232   223
#     600   8.0   0.027   11250          0.184      236          0.184   329
#
# Taking T_L = 2ℓ/σ_w from the BL89 mixing length gives σ_w = K/(2ℓ), which
# collapses with height: ℓ_geom = 0.25·z carries no stability correction, so by
# 600 m σ_w is an order below the observed neutral 1.3·u* ≈ 0.30 m/s and T_L is
# three hours.  The variance, not the diffusivity, is the quantity similarity
# theory constrains well.
#
# So σ_w is kept at its Hanna (1982) value — unchanged, which makes `tke` here
# BIT-IDENTICAL to the Hanna fill, including the Thomson drift term that reads
# ∂σ_w/∂z from it — and the ARCHIVED K sets the Lagrangian timescale instead:
#
#     T_L = K_z / σ_w²        (so σ_w²·T_L = K_z exactly)
#
# The model's effective vertical diffusivity therefore becomes the host model's
# own, while the velocity variance stays on firm observational ground.  This is
# a clean single factor against a Hanna run: only `tls` changes.
#
# Above the boundary layer the archived K falls to ~0, so T_L reverts to the
# Hanna value across the SAME entrainment zone [EZ_LO·h, EZ_HI·h] the σ_w fill
# uses.  `metfields3D["tke"]` holds the archived K on entry (the namelist maps
# `tke = avg_tdch`) and 0.5·σ_w² on exit; the column is copied out first.
function kz_fill_tke!(metfields3D, metfields2D, indt::Int,
                      Ni::Int, Nj::Int, Nk::Int;
                      sigw_ft::Float32 = 0.1f0, kmin::Float32 = 1.0f-3,
                      tl_ft::Float32 = 100.0f0, kmin_ft::Float32 = 0.5f0,
                      sigw_ft_max::Float32 = 2.0f0)
    Kcol = Vector{Float32}(undef, Nk)
    for j in 1:Nj, i in 1:Ni
        @inbounds for k in 1:Nk
            Kcol[k] = max(Float32(metfields3D["tke"][i, j, k, indt]), 0.0f0)
        end
        ust = Float32(metfields2D["ustar"][i, j, indt])
        pbl = Float32(metfields2D["pblh"][i, j, indt])
        ol  = Float32(metfields2D["ol"][i, j, indt])
        wst = hanna_wstar(ust, pbl, ol)
        h   = max(pbl, 50.0f0)
        z1  = EZ_LO * h
        z2  = EZ_HI * h
        @inbounds for k in 1:Nk
            z = Float32(metfields3D["zlev"][i, j, k, indt])
            sigw, _, tlw_h = hanna_sigw_tlw(z, pbl, ust, wst, ol, sigw_ft)
            K   = Kcol[k]
            tlw = tlw_h
            if z >= z2
                # ── FREE TROPOSPHERE ───────────────────────────────────────
                # Hanna has nothing here: it returns a CONSTANT sigw_ft with no
                # spatial structure.  The archived field does carry structure —
                # on 20220115 18Z, 31 % of points 500 m above the BL top exceed
                # 1 m²/s and 6.5 % exceed 10, clustered at |lat| 25-32° (convective
                # outflow), with a second population at 8-12 km near |lat| 38-48°
                # (jet shear / clear-air turbulence).  Discarding it was the
                # first version's main limitation.
                #
                # No similarity theory gives sigma_w up here, so the inversion
                # runs the other way from the PBL branch: hold a free-troposphere
                # Lagrangian timescale and let K set the VARIANCE,
                #     sigma_w = sqrt(K/T_L,ft),
                # capped so a single extreme cell cannot produce a runaway
                # velocity.  Below kmin_ft the layer is quiescent and the Hanna
                # free-troposphere pair is kept, so clear air is unchanged.
                if K > kmin_ft
                    sigw = clamp(sqrt(K / tl_ft), sigw_ft, sigw_ft_max)
                    tlw  = clamp(K / max(sigw * sigw, 1.0f-6), 0.5f0, 300.0f0)
                end
            elseif K > kmin
                # ── BOUNDARY LAYER ─────────────────────────────────────────
                # sigma_w stays at its Hanna value (so the stored tke, and the
                # Thomson drift read from its gradient, match the Hanna fill) and
                # the archived K sets the timescale: sigma_w^2*T_L = K exactly.
                tl_k = clamp(K / max(sigw * sigw, 1.0f-6), 0.5f0, 300.0f0)
                if z > z1
                    x = (z - z1) / (z2 - z1)
                    w = 0.5f0 * (1.0f0 + cos(Float32(pi) * x))      # 1 → 0
                    tlw = w * tl_k + (1.0f0 - w) * tlw_h
                else
                    tlw = tl_k
                end
            end
            metfields3D["tke"][i, j, k, indt] = 0.5f0 * sigw * sigw
            metfields3D["tls"][i, j, k, indt] = 0.5f0 * tlw * sigw
        end
    end
    return nothing
end

# ── 2D surface-field bilinear + time interpolation ───────────────────────────
# field: (Ni, Nj, 2) — two bracketing time slices (same layout as 3D arrays).
@inline function interpolate2Dsingle(
    field ::AbstractArray{Float32,3},
    indx  ::Float32,
    indy  ::Float32,
    indt  ::Float32,
)::Float32
    dl, l = _split_t(indt)
    i1 = Base.unsafe_trunc(Int32, indx)
    j1 = Base.unsafe_trunc(Int32, indy)
    i2 = i1 + Int32(1)
    j2 = j1 + Int32(1)
    l2 = min(l + 1, size(field, 3))
    di = Float32(mod(indx, 1.0f0))
    dj = Float32(mod(indy, 1.0f0))
    w00 = (1.0f0 - di) * (1.0f0 - dj)
    w01 = (1.0f0 - di) * dj
    w10 = di * (1.0f0 - dj)
    w11 = di * dj
    @inbounds begin
        v1 = w00*field[i1,j1,l]  + w01*field[i1,j2,l]  +
             w10*field[i2,j1,l]  + w11*field[i2,j2,l]
        v2 = w00*field[i1,j1,l2] + w01*field[i1,j2,l2] +
             w10*field[i2,j1,l2] + w11*field[i2,j2,l2]
    end
    return v1 * (1.0f0 - dl) + v2 * dl
end

# ── Option b: on-the-fly Hanna turbulence step ───────────────────────────────
#
# Langevin structure mirrors turbulence_v3(::Thomson,...) but derives σ_w,
# dσ_w/dz, and T_L analytically from Hanna (1982) / FLEXPART instead of from
# a stored 3D TKE field.
#
# The analytically computed dσ_w/dz eliminates the finite-difference
# approximation used in the pre-compute option, giving exact consistency
# between the drift term and the σ_w profile.
#
# `rhoprof` / `rhogradprof` are ρ and ∂ρ/∂z at the half levels `height`
# (get_density_profiles_half!, met/interp.jl), not the raw class-F columns.
#
# Returns: (z_new, delz, wp_new, tlw, tke_z, yl)
# v1.1: count of guarded w' resets in turbulence_hanna_v1 (trap-#1 visibility)
const _HANNA_WP_RESETS = Threads.Atomic{Int}(0)

@inline function turbulence_hanna_v1(
    hanna       ::Hanna,
    indl        ::Float32,
    dt          ::Float32,
    ldt         ::Float32,
    ifine       ::Float32,
    wp          ::Float32,
    height      ::AbstractVector{Float32},
    rhoprof     ::AbstractVector{Float32},
    rhogradprof ::AbstractVector{Float32},
    ustar       ::Float32,
    pblh        ::Float32,
    wstar       ::Float32,
    ol          ::Float32,
)
    dt_left ::Float32 = dt
    yl      ::Float32 = 0.0f0
    tlw     ::Float32 = 1.0f0
    sigw    ::Float32 = 0.2f0
    tke_z   ::Float32 = 0.0f0
    z       ::Float32 = interpolation1Dsingle(height, indl)
    delz    ::Float32 = 0.0f0
    zwall   ::Float32 = 0.0f0

    while dt_left > eps(Float32)
        z       = interpolation1Dsingle(height, indl)
        rhoa    = interpolation1Dsingle(rhoprof, indl)
        rhograd = interpolation1Dsingle(rhogradprof, indl)
        rhoaux  = rhograd / max(rhoa, 1.0f-12)

        # Hanna analytical σ_w, dσ_w/dz, T_L with the blended entrainment
        # zone across the PBL top (see hanna_sigw_tlw).
        sigw, dsigwdz, tlw = hanna_sigw_tlw(z, pblh, ustar, wstar, ol, hanna.sigw_ft)

        yl    = 0.5f0 * tlw * sigw    # mixing length ℓ = T_L·σ_w/2
        tke_z = 0.5f0 * sigw^2

        dt_sub ::Float32 = max(min(tlw / 5.0f0 / ifine, dt / ifine), 0.5f0)
        dt_sub = min(dt_sub, dt_left)
        dt_left -= dt_sub

        dtftlw ::Float32 = dt_sub / max(tlw, 0.5f0)
        rw     ::Float32 = exp(-dtftlw)

        term1  ::Float32 = rw * wp
        term2  ::Float32 = clamp(randn(Float32), -3.0f0, 3.0f0) *
                           sqrt(max(1.0f0 - rw * rw, 0.0f0))
        term3  ::Float32 = tlw * (1.0f0 - rw) * (dsigwdz + rhoaux * sigw)

        wp_try ::Float32 = term1 + term2 + term3
        # v1.1: guard against a rare non-finite / runaway normalised velocity
        # (test_wellmixed_pbl produced a NaN height in 2 of 9 threaded runs;
        # the state was not captured).  |w'| > 20 sigma is unphysical: reset
        # to a fresh draw, count it, warn once so the event stays visible.
        if !(isfinite(wp_try) && abs(wp_try) <= 20.0f0)
            Threads.atomic_add!(_HANNA_WP_RESETS, 1)
            @warn "hanna: non-finite or runaway w' reset (see _HANNA_WP_RESETS)" maxlog=1 wp_try term1 term2 term3 sigw dsigwdz tlw z
            wp_try = clamp(randn(Float32), -3.0f0, 3.0f0)
        end
        delz            = wp_try * sigw * dt_sub
        z_try  ::Float32 = z + delz

        if z_try < zwall
            z    = 2.0f0 * zwall - z_try
            wp   = -wp_try
            z    = max(z, zwall + 1.0f-4)
            indl = _indl_below_top(height, z)
        else
            z  = z_try
            wp = wp_try
        end

        z    = clamp(z, height[1] + 1.0f-4, height[end] - 1.0f-4)
        indl = _indl_below_top(height, z)
    end

    return z, delz, wp, tlw, tke_z, yl
end

# ── Horizontal Hanna closure ([physics.turbulence].horizontal_hanna) ─────────
#
# Hanna (1982) gives σ_u, σ_v, T_Lu, T_Lv in the PBL from the SAME (u*, w*, h, L)
# the vertical fill uses; hanna_bl already computes them, they were never wired.
# This is FLEXPART's in-PBL horizontal Langevin (advance.f90:
#   up = ru·up + N(0,1)·sigu·√(1−ru²), ru = exp(−dt/tlu), likewise vp with
#   sigv/tlv), applied to the GRID components with no rotation into the wind
# frame — the FLEXPART convention, kept here so the two models agree.
#
# Above the PBL Hanna is undefined.  The same cosine entrainment-zone weight
# hanna_sigw_tlw uses for σ_w fades the horizontal VARIANCES to zero across
# [EZ_LO·h, EZ_HI·h]; the free-troposphere horizontal term, if wanted, remains
# [physics.turbulence] kh_ft (FLEXPART d_trop).  T_L is held at its last BL
# value above so a parcel's u′/v′ decays smoothly instead of being cut.
#
# Returns (sigu, tlu, sigv, tlv); sigu = sigv = 0 above EZ_HI·h.
@inline function hanna_h(z::Float32, pblh::Float32, ol::Float32,
                         ust::Float32, wst::Float32)::NTuple{4,Float32}
    h  = max(pblh, 50.0f0)
    z1 = EZ_LO * h
    z2 = EZ_HI * h
    z >= z2 && return (0.0f0, 300.0f0, 0.0f0, 300.0f0)
    sigu, sigv, _, _, tlu, tlv, _ = hanna_bl(z, pblh, ol, ust, wst)
    z <= z1 && return (sigu, tlu, sigv, tlv)
    x = (z - z1) / (z2 - z1)
    w = 0.5f0 * (1.0f0 + cos(Float32(π) * x))                 # 1 → 0
    return (sqrt(w) * sigu, tlu, sqrt(w) * sigv, tlv)         # variances × w
end

