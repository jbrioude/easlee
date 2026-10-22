
# ── Vertical turbulence — Thomson Langevin step on TKE profiles (BL89 ℓ) ──────
# Stochastic vertical displacement of one parcel per model sub-step from a
# 1-D TKE/θ/ρ profile: turbulence_v3 is the active scheme (dispatch tags
# Thomson / Hanna / NoTurb); turbulent_length_scale_v2 (Bougeault–Lacarrère
# 1989-style mixing length) and tke_partition (per-driver σ_w/T_L closure,
# Fix C) are its building blocks.  The Langevin drift matches FLEXPART's
# well-mixed formulation (Thomson 1987); density correction included.
# 4-D/GPU variant: turbulence_v3_4d in turbulence_4d.jl.

# Turbulent mixing length from TKE and local stability.
# Returns ℓmix [m].
function turbulent_length_scale_v2(indz, indzp,
                                    tke_here::Float32,
                                    z_here::Float32,
                                    z_prof::AbstractVector{Float32},
                                    theta_prof::AbstractVector{Float32};
                                    c0::Float32  = 0.25f0,
                                    ℓmax::Float32 = 500.0f0)
    dist_bot = max(z_here - z_prof[1],   1f-3)
    dist_top = max(z_prof[end] - z_here, 1f-3)
    ℓ_geom   = c0 * min(dist_bot, dist_top)

    g    = 9.81f0
    θ0   = theta_prof[indz]
    θ1   = theta_prof[indzp]
    z0   = z_prof[indz]
    z1   = z_prof[indzp]
    dθdz = (θ1 - θ0) / (z1 - z0 + 1f-12)
    θ_bar = 0.5f0 * (θ0 + θ1)
    N2   = (g / θ_bar) * dθdz

    if N2 > 0f0
        w_tke = sqrt(max(2f0 * tke_here, 0f0))
        ℓ_b   = w_tke / sqrt(max(N2, 1f-12))
    else
        ℓ_b = ℓmax
    end

    return min(ℓ_geom, ℓ_b, ℓmax)
end

# Floors for σ_w and T_L (env-tunable for well-mixed sensitivity tests).  The defaults
# (0.2 m/s, 0.5 s) put a kink in σ_w(z) and clamp T_L near the surface, both of which can
# break the well-mixed balance at the floor→active transition.
const SIGW_FLOOR = parse(Float32, get(ENV, "SIGW_FLOOR", "0.2"))
const TLW_FLOOR  = parse(Float32, get(ENV, "TLW_FLOOR",  "0.5"))

# ── Fix C: driver-consistent vertical mixing (ADAPTIVE_KERNEL §8) ─────────────
# The met field named "tke" carries DIFFERENT quantities per driver:
#   • Hanna-fill drivers (ERA5 / GFS / MERRA-2 / AROME-HDF): the reader writes
#     0.5·σ_w,H² (w-variance only, hanna_fill_tke!).  √(2·tke) is exactly
#     Hanna's σ_w — these runs are correct and MUST stay bit-for-bit identical.
#   • Native-e drivers (MesoNH, AROME native e): the full e = ½(σ_u²+σ_v²+σ_w²).
#     There √(2e) assigns all three components' energy to w, and the implied
#     diffusivity K = σ_w²·τ_Lw = 2√2·ℓ·√e is over an order of magnitude above
#     the CBR scalar closure C_K·ℓ·√e.
# The closure that maps the archived field to the two quantities transport
# actually needs — σ_w² and K_z — is a per-DRIVER property.  It is selected by
# the reader from the listvariables/*.vars namelist (#info keys tke_kind /
# vmix_rw / vmix_ck — NO closure constants live in this file) and published
# once per run via set_vertical_mixing!.  τ_Lw = K_z/σ_w² everywhere.
"""
Per-driver vertical-mixing closure mapping the archived "tke" field to the
(σ_w², K_z) transport needs (see the Fix C note above).  `native_e=false` =
Hanna-fill semantics (field is 0.5·σ_w²); `native_e=true` = full TKE e with
CBR constants `r_w`/`c_k` from the variables namelist.
"""
struct VMixClosure
    native_e :: Bool      # true → archived field is the full TKE e (CBR closure)
    r_w      :: Float32   # CBR vertical variance share: σ_w² = 2·e·r_w   (namelist)
    c_k      :: Float32   # CBR scalar diffusivity const: K_z = C_K·ℓ·√e  (namelist)
end

# Default = Hanna-fill semantics (all current ERA5/GFS/MERRA-2 runs).  The CBR
# constants are deliberately NaN here: they are namelist data, never code data,
# and the Hanna branch never reads them.
const VMIX = Ref(VMixClosure(false, NaN32, NaN32))

"""
    set_vertical_mixing!(native_e; r_w=NaN32, c_k=NaN32)

Publish the per-driver vertical-mixing closure (called by the runtime after the
reader has parsed the variables namelist; see vmix_from_namelist).
"""
function set_vertical_mixing!(native_e::Bool; r_w::Float32=NaN32, c_k::Float32=NaN32)
    if native_e
        (isfinite(r_w) && r_w > 0.0f0) ||
            error("set_vertical_mixing!: native-e closure needs r_w > 0 (namelist vmix_rw)")
        (isfinite(c_k) && c_k > 0.0f0) ||
            error("set_vertical_mixing!: native-e closure needs c_k > 0 (namelist vmix_ck)")
    end
    VMIX[] = VMixClosure(native_e, r_w, c_k)
    return VMIX[]
end
"""
    vertical_mixing() → VMixClosure

The run-global vertical-mixing closure published by `set_vertical_mixing!`.
"""
vertical_mixing() = VMIX[]

# w-share of the archived field's vertical gradient: σ_w² = 2·wshare·field, so
# ∂σ_w/∂z = wshare·∂(field)/∂z / σ_w.  Hanna fill: field = 0.5σ_w² → wshare = 1
# (multiplying by the literal 1.0f0 is exact in IEEE — bit-for-bit invariant).
@inline _vmix_wshare(vm::VMixClosure) = vm.native_e ? vm.r_w : 1.0f0

# Partition the archived turbulence field into (τ_Lw [s], σ_w [m/s]).
# Dissolved form of the old tke_partition: τ_Lw = K_z/σ_w² with σ_w²/K_z from
# the per-driver closure.
@inline function tke_partition(vm::VMixClosure, tke_z::Float32, yl::Float32)
    if vm.native_e
        # Native-e (CBR): σ_w² = 2·e·r_w (vertical variance share) and
        # K_z = C_K·ℓ·√e (CBR SCALAR constant; ℓ = yl from
        # turbulent_length_scale_v2).  τ_Lw = K_z/σ_w², floors as before.
        e_z  = max(tke_z, 0.0f0)
        sigw = max(sqrt(2.0f0 * vm.r_w * e_z), SIGW_FLOOR)
        kz   = vm.c_k * yl * sqrt(e_z)
        tlw  = clamp(kz / (sigw * sigw), TLW_FLOOR, 60.0f0)
        return tlw, sigw
    else
        # Hanna fill: field = 0.5·σ_w,H² ⇒ σ_w = √(2·field) = σ_w,H exactly.
        # K_z is DEFINED as σ_w²·τ_Lw,H with τ_Lw,H = 2·yl/σ_w (the historic
        # partition), so
        #     τ_Lw = K_z/σ_w² = (σ_w²·2·yl/σ_w)/σ_w² = 2·yl/σ_w
        # — the σ_w² cancels ALGEBRAICALLY (not numerically), so the expression
        # below is the same Float32 op sequence as the old tke_partition:
        # bit-for-bit identical on ERA5/GFS/MERRA-2 by construction.
        sigw = max(sqrt(2.0f0 * tke_z), SIGW_FLOOR)
        tlw  = clamp(2.0f0 * yl / sigw, TLW_FLOOR, 60.0f0)
        return tlw, sigw
    end
end

# Back-compat 2-arg form: reads the run-global closure.  CPU callers only —
# GPU kernels must take the closure as an argument (a Ref load is a host
# pointer dereference, illegal in device code); see turbulence_v3_4d.
tke_partition(tke_z::Float32, yl::Float32) = tke_partition(VMIX[], tke_z, yl)

# ── Horizontal partition of the SAME archived field (native-e drivers only) ───
# A native-e driver (AROME native e, MesoNH) archives the full
# e = ½(σ_u²+σ_v²+σ_w²), so the horizontal variance is already in the field the
# vertical step reads — no new met input is needed to drive a horizontal
# Langevin step.  Horizontal isotropy:
#     σ_u² = σ_v² = 2·e·r_h  with  r_h = (1 − r_w)/2
# so r_w + 2·r_h = 1 and the three variances sum back to 2e EXACTLY: the
# horizontal share is whatever the namelist's vertical share (vmix_rw) leaves
# over, never an independent tunable.  K_h reuses the CBR scalar closure of
# K_z (C_K·ℓ·√e), hence τ_Lh = K_h/σ_u².  τ_Lh is capped at 600 s (horizontal
# Lagrangian timescales are an order above the vertical 60 s cap).
#
# Hanna-fill drivers (ERA5/GFS/MERRA-2/AROME-HDF) archive 0.5·σ_w,H² — the
# w-variance ONLY, with no horizontal information in it — so this returns
# (0, 0) there and build_run refuses [physics.turbulence].horizontal_tke.
@inline function tke_partition_h(vm::VMixClosure, tke_z::Float32, yl::Float32)
    vm.native_e || return (0.0f0, 0.0f0)
    e_z  = max(tke_z, 0.0f0)
    r_h  = 0.5f0 * (1.0f0 - vm.r_w)
    sigh = max(sqrt(2.0f0 * r_h * e_z), SIGW_FLOOR)
    kh   = vm.c_k * yl * sqrt(e_z)
    tlh  = clamp(kh / (sigh * sigh), TLW_FLOOR, 600.0f0)
    return tlh, sigh
end

# Level index of height z on the profile, kept below the top level: at a 16 km
# top, height[end] − 1e-4 rounds to height[end] in Float32 (spacing 0.002 m), the
# search returned exactly Nk and the next interpolation1Dsingle read element
# Nk+1 (AROME-HDF with turbulence, BoundsError 43 of 42, 2026-10-22).  The KA
# bodies are protected by their samplers' top clamp.
@inline _indl_below_top(height::AbstractVector{Float32}, z::Float32) =
    min(interpolation1Dsingleirreg(height, z), Float32(length(height)) - 1.0f-3)

# ── NoTurb pass-through ───────────────────────────────────────────────────────
# Pure advection — no stochastic vertical displacement.
# Returns the same signature as turbulence_v3 but with zero perturbation.
@inline function turbulence_v3(
    ::NoTurb,
    indl::Float32, dt::Float32, ldt::Float32, ifine::Float32, wp::Float32,
    height::AbstractVector{Float32},
    pttprof::AbstractVector{Float32},
    tkeprof::AbstractVector{Float32},
    rhoprof::AbstractVector{Float32},
    rhogradprof::AbstractVector{Float32},
    sigw_col::AbstractVector{Float32}=Float32[],
    tlw_col::AbstractVector{Float32}=Float32[],
    yl_col::AbstractVector{Float32}=Float32[],
    ; vm::VMixClosure=VMIX[],   # v1.1: accepted and ignored (pure advection)
)
    z = interpolation1Dsingle(height, indl)
    return z, 0.0f0, 0.0f0, ldt, 0.0f0, 0.0f0, indl   # (z, delz, wp, tlw, tke_z, yl, indl)
end

# ── Thomson (BL89) full turbulence step ───────────────────────────────────────
# Advances a single parcel through one model sub-step using the stochastic
# Langevin equation.  Returns (z_new, delz, wp_new, tlw, tke_z, yl).
#
# Every profile is indexed like `height` (the zlev column, half levels): tke is
# class H as archived; pttprof, rhoprof and rhogradprof must be the class-F
# θ, ρ and ∂ρ/∂z AT the half levels (get_vertical_profile_half! /
# get_density_profiles_half!, met/interp.jl — layout contract 2026-10-22).
#
# dsigwdz is computed online as a forward difference of clamped sigw over the
# layer [indz, indz+1] containing the parcel.  This is consistent with
# tke_partition and avoids the spurious drift that a precomputed centred
# difference produces in the TKE-floor region above the PBL top.
function turbulence_v3(
    ::Thomson,
    indl::Float32,
    dt::Float32,
    ldt::Float32,
    ifine::Float32,
    wp::Float32,
    height::AbstractVector{Float32},
    pttprof::AbstractVector{Float32},
    tkeprof::AbstractVector{Float32},
    rhoprof::AbstractVector{Float32},
    rhogradprof::AbstractVector{Float32},
    sigw_col::AbstractVector{Float32}=Float32[],   # unused (column memo reverted — broke well-mixed)
    tlw_col::AbstractVector{Float32}=Float32[],
    yl_col::AbstractVector{Float32}=Float32[];
    # This nest's closure.  Defaults to the run-global one so every existing
    # single-met call site is unchanged; a nested run passes the nest's own.
    vm::VMixClosure=VMIX[],
)
    dt_left::Float32 = dt
    yl::Float32      = 0.0f0
    tlw::Float32     = 1.0f0
    sigw::Float32    = 0.2f0
    tke_z::Float32   = 0.0f0
    z::Float32       = interpolation1Dsingle(height, indl)
    delz::Float32    = 0.0f0
    zwall::Float32   = 0.0f0
    # Per-driver σ_w²/K_z closure (Fix C) — hoisted out of the sub-step loop.
    wshare  = _vmix_wshare(vm)

    while dt_left > eps(Float32)
        z        = interpolation1Dsingle(height, indl)
        tke_z    = interpolation1Dsingle(tkeprof, indl)
        rhoa     = interpolation1Dsingle(rhoprof, indl)
        rhograd  = interpolation1Dsingle(rhogradprof, indl)
        rhoaux   = rhograd / max(rhoa, 1.0f-12)

        # Use the layer *containing* the parcel (trunc(indl) gives the lower bound).
        indz  = clamp(Base.unsafe_trunc(Int, indl), 1, length(height) - 1)
        indzp = indz + 1

        yl      = turbulent_length_scale_v2(indz, indzp, tke_z, z, height, pttprof)
        tlw, sigw = tke_partition(vm, tke_z, yl)
        # ∂σ_w/∂z CONSISTENT with the value σ_w²=2·wshare·field: ∂σ_w/∂z =
        # (wshare/σ_w)·∂(field)/∂z using the SAME layer-linear field slope
        # (√(linear)≠linear(√), so differencing √ at the endpoints is
        # inconsistent with the value).  Zero where σ_w sits at its floor.
        # Hanna fill: wshare = 1 (exact multiply — unchanged numbers).
        dz_layer  = height[indzp] - height[indz]
        dtkedz    = dz_layer > 0.1f0 ? (tkeprof[indzp] - tkeprof[indz]) / dz_layer : 0.0f0
        dsigwdz   = sigw > SIGW_FLOOR ? wshare * dtkedz / sigw : 0.0f0

        dt_sub::Float32 = max(min(tlw / 5.0f0 / ifine, dt / ifine), 0.5f0)
        dt_sub = min(dt_sub, dt_left)
        dt_left -= dt_sub

        dtftlw::Float32 = dt_sub / max(tlw, 0.5f0)
        rw::Float32     = exp(-dtftlw)

        term1::Float32 = rw * wp
        term2::Float32 = clamp(randn(Float32), -3.0f0, 3.0f0) * sqrt(max(1.0f0 - rw * rw, 0.0f0))
        # Well-mixed drift for the NORMALISED velocity ŵ=wp.  NB: transforming the physical
        # Thomson drift ½(1+ŵ²)∂σ_w²/∂z to ŵ=w/σ(z) makes the velocity-dependent ŵ² piece
        # CANCEL against the Itô coordinate term (∂ŵ/∂z·dz), leaving exactly ∂σ_w/∂z +
        # σ_w·ρ⁻¹∂ρ/∂z.  (Adding a (1+ŵ²) factor here double-counts it → BL runaway.)
        term3::Float32 = tlw * (1.0f0 - rw) * (dsigwdz + rhoaux * sigw)

        wp_try::Float32 = term1 + term2 + term3
        # a NaN met value in the column (masked cell, or the 2026-10-22 GFS
        # surface slot) made wp, z and then the level index NaN, and the next
        # profile read went out of bounds: reset to a fresh draw, as the Hanna
        # bodies do (turbulence_hanna.jl), warned once
        if !isfinite(wp_try)
            @warn "thomson: non-finite w' reset (NaN in the turbulence profiles)" maxlog=1 term1 term2 term3 sigw dsigwdz rhoaux tlw z
            wp_try = clamp(randn(Float32), -3.0f0, 3.0f0)
        end
        delz = wp_try * sigw * dt_sub
        z_try::Float32  = z + delz

        if z_try < zwall
            z  = 2.0f0 * zwall - z_try
            wp = -wp_try
            z  = max(z, zwall + 1.0f-4)
            indl = _indl_below_top(height, z)
        else
            z  = z_try
            wp = wp_try
        end

        z    = clamp(z, height[1] + 1.0f-4, height[end] - 1.0f-4)
        indl = _indl_below_top(height, z)
    end

    return z, delz, wp, tlw, tke_z, yl, _indl_below_top(height, z)
end

# ── Hanna option a (pre-computed TKE) ─────────────────────────────────────────
# When turbulence = Hanna(false), hanna_fill_tke! has already written σ_w/T_L
# into the tke/tls fields at met-read time, so the per-parcel step IS the
# Thomson Langevin update reading those profiles.  (Hanna(true) / option b is
# the on-the-fly path handled in advance.jl before turbulence_v3 is reached.)
# The Thomson body uses its tag only for dispatch, so the dummy tag is safe.
@inline turbulence_v3(::Hanna,
    indl::Float32, dt::Float32, ldt::Float32, ifine::Float32, wp::Float32,
    height::AbstractVector{Float32}, pttprof::AbstractVector{Float32},
    tkeprof::AbstractVector{Float32}, rhoprof::AbstractVector{Float32},
    rhogradprof::AbstractVector{Float32},
    sigw_col::AbstractVector{Float32}=Float32[],
    tlw_col::AbstractVector{Float32}=Float32[],
    yl_col::AbstractVector{Float32}=Float32[];
    vm::VMixClosure=VMIX[]) =                       # v1.1: this nest's closure, forwarded
    turbulence_v3(Thomson(true, false), indl, dt, ldt, ifine, wp, height,
                  pttprof, tkeprof, rhoprof, rhogradprof, sigw_col, tlw_col, yl_col;
                  vm = vm)
