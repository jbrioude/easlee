
# ── Buffer-free turbulence kernel ─────────────────────────────────────────────
# turbulence_v3_4d is equivalent to turbulence_v3 but operates directly on the
# 4-D met arrays instead of pre-extracted 1-D profiles.
#
# Why: on GPU there is no cheap way to allocate per-thread Vector{Float32}(nz).
# Replacing get_vertical_profile! + interpolation1Dsingle with direct 4D
# interpolation removes all heap allocation from the hot loop and makes the
# function GPU-kernel-safe.
#
# Performance on CPU: one extra bilinear interpolation per field lookup vs. the
# profile approach, but the profile extraction is eliminated entirely → net
# wash or small win at typical nz=60.

# ── Turbulence length scale (4D version) ─────────────────────────────────────
# Accesses the theta/z profiles at exactly two integer level indices.  theta is
# class F: its value at the half levels indz/indzp (heights zlev[indz]/zlev[indzp])
# comes from half_level_value_4d, the same stencil as the CPU profile
# (get_vertical_profile_half!, layout contract 2026-10-22).
@inline function turbulent_length_scale_v2_4d(
    indz::Int, indzp::Int,
    tke_here::Float32, z_here::Float32,
    indx::Float32, indy::Float32, indt::Float32, nk::Int32,
    zlev_arr  ::AbstractArray{Float32,4},
    theta_arr ::AbstractArray{Float32,4},
)
    z_bot  = interpolation4Dlevsingle(zlev_arr,  indx, indy, 1.0f0,         indt)
    z_top  = interpolation4Dlevsingle(zlev_arr,  indx, indy, Float32(nk),   indt)
    dist_bot = max(z_here - z_bot, 1.0f-3)
    dist_top = max(z_top  - z_here, 1.0f-3)
    ℓ_geom   = 0.25f0 * min(dist_bot, dist_top)

    θ0 = half_level_value_4d(theta_arr, indx, indy, indz,  indt)
    θ1 = half_level_value_4d(theta_arr, indx, indy, indzp, indt)
    z0 = interpolation4Dlevsingle(zlev_arr,  indx, indy, Float32(indz),  indt)
    z1 = interpolation4Dlevsingle(zlev_arr,  indx, indy, Float32(indzp), indt)

    dθdz  = (θ1 - θ0) / (z1 - z0 + 1.0f-12)
    θ_bar = 0.5f0 * (θ0 + θ1)
    N2    = (9.81f0 / θ_bar) * dθdz
    ℓmax  = 500.0f0

    ℓ_b = N2 > 0.0f0 ?
        sqrt(max(2.0f0 * tke_here, 0.0f0)) / sqrt(max(N2, 1.0f-12)) :
        ℓmax

    return min(ℓ_geom, ℓ_b, ℓmax)
end

# ── NoTurb pass-through (4D, GPU-safe) ───────────────────────────────────────
@inline function turbulence_v3_4d(
    ::NoTurb,
    indl   ::Float32, indx::Float32, indy::Float32,
    dt     ::Float32, ldt::Float32,  ifine::Float32, wp::Float32,
    indt   ::Float32, nk::Int32,     rng_state::UInt64,
    zlev_arr   ::AbstractArray{Float32,4},
    theta_arr  ::AbstractArray{Float32,4},
    tke_arr    ::AbstractArray{Float32,4},
    rho_arr    ::AbstractArray{Float32,4},
    vmix       ::VMixClosure,          # unused; keeps the call sites uniform
)
    z = interpolation4Dlevsingle(zlev_arr, indx, indy, indl, indt)
    return z, 0.0f0, 0.0f0, ldt, 0.0f0, 0.0f0, rng_state
end

# ── Main turbulence step (4D, GPU-safe) ───────────────────────────────────────
# Replaces turbulence_v3.  Arguments:
#   ::Thomson          — dispatch tag (same as before)
#   indl               — current vertical index (fractional)
#   indx, indy         — horizontal parcel position (for 4D lookup)
#   dt, ldt, ifine, wp — same as turbulence_v3
#   indt               — met time index
#   nk                 — number of vertical levels
#   zlev, theta, tke, rho — full 4D met arrays (dsigwdz computed online; ρ, θ
#                        and ∂ρ/∂z taken at the half levels, see met/interp.jl)
#   rng_state          — UInt64 xorshift state (updated on return)
#
# Returns: (z_new, delz, wp_new, ldt_new, tke_z, yl, rng_state_new)
@inline function turbulence_v3_4d(
    ::Thomson,
    indl   ::Float32,
    indx   ::Float32,
    indy   ::Float32,
    dt     ::Float32,
    ldt    ::Float32,
    ifine  ::Float32,
    wp     ::Float32,
    indt   ::Float32,
    nk     ::Int32,
    rng_state::UInt64,
    zlev_arr   ::AbstractArray{Float32,4},
    theta_arr  ::AbstractArray{Float32,4},
    tke_arr    ::AbstractArray{Float32,4},
    rho_arr    ::AbstractArray{Float32,4},
    vmix       ::VMixClosure,   # per-driver σ_w²/K_z closure (Fix C) — passed as a
                                # kernel argument: GPU code cannot read the VMIX Ref
)
    dt_left ::Float32 = dt
    yl      ::Float32 = 0.0f0
    tlw     ::Float32 = 1.0f0
    sigw    ::Float32 = 0.2f0
    tke_z   ::Float32 = 0.0f0
    z       ::Float32 = interpolation4Dlevsingle(zlev_arr, indx, indy, indl, indt)
    delz    ::Float32 = 0.0f0
    zwall   ::Float32 = 0.0f0
    wshare  ::Float32 = _vmix_wshare(vmix)

    while dt_left > eps(Float32)
        z       = interpolation4Dlevsingle(zlev_arr,   indx, indy, indl, indt)
        tke_z   = interpolation4Dlevsingle(tke_arr,    indx, indy, indl, indt)
        # ρ is class F: half-level ρ and ∂ρ/∂z, as the CPU profiles (contract)
        rhoa, rhograd = half_level_density_4d(rho_arr, zlev_arr, indx, indy, indl, indt)
        rhoaux  = rhograd / max(rhoa, 1.0f-12)

        indz  = clamp(Base.unsafe_trunc(Int, indl), 1, Int(nk) - 1)
        indzp = indz + 1

        yl       = turbulent_length_scale_v2_4d(indz, indzp, tke_z, z,
                                                indx, indy, indt, nk,
                                                zlev_arr, theta_arr)
        tlw, sigw = tke_partition(vmix, tke_z, yl)

        # ∂σ_w/∂z CONSISTENT with the value σ_w² = 2·wshare·field — the same form as
        # the CPU body (turbulence.jl:243): (wshare/σ_w)·∂(field)/∂z on the layer's
        # linear field slope, zero where σ_w sits at its floor.
        # 2026-09-22: this body still differenced the CLAMPED σ_w at the two levels
        # ((σ_k+1 − σ_k)/Δz), the form the CPU comment warns is inconsistent
        # (√(linear) ≠ linear(√)).  The CPU had been corrected, the KA twin had not,
        # and the two advances carried systematically different drift terms — a
        # +1.7 % mass offset on the KA side in an all-switches-off control
        # (docs/BENCHMARKS.md).  Aligned with the reference.
        tke_k    = interpolation4Dlevsingle(tke_arr,  indx, indy, Float32(indz),  indt)
        tke_kp1  = interpolation4Dlevsingle(tke_arr,  indx, indy, Float32(indzp), indt)
        z_k      = interpolation4Dlevsingle(zlev_arr, indx, indy, Float32(indz),  indt)
        z_kp1    = interpolation4Dlevsingle(zlev_arr, indx, indy, Float32(indzp), indt)
        dz_layer = z_kp1 - z_k
        dtkedz   = dz_layer > 0.1f0 ? (tke_kp1 - tke_k) / dz_layer : 0.0f0
        dsigwdz  = sigw > SIGW_FLOOR ? wshare * dtkedz / sigw : 0.0f0   # same floor constant as the CPU body

        dt_sub ::Float32 = max(min(tlw / 5.0f0 / ifine, dt / ifine), 0.5f0)
        dt_sub = min(dt_sub, dt_left)
        dt_left -= dt_sub

        dtftlw ::Float32 = dt_sub / max(tlw, 0.5f0)
        rw     ::Float32 = exp(-dtftlw)

        rng_state, zrand = xorshift_randn_clamped(rng_state)

        term1  ::Float32 = rw * wp
        term2  ::Float32 = zrand * sqrt(max(1.0f0 - rw * rw, 0.0f0))
        term3  ::Float32 = tlw * (1.0f0 - rw) * (dsigwdz + rhoaux * sigw)

        wp_try ::Float32 = term1 + term2 + term3
        if !isfinite(wp_try)            # NaN met value in the column: as turbulence_v3
            rng_state, zr2 = xorshift_randn_clamped(rng_state)
            wp_try = zr2
        end
        delz            = wp_try * sigw * dt_sub
        z_try  ::Float32 = z + delz

        if z_try < zwall
            z    = 2.0f0 * zwall - z_try
            wp   = -wp_try
            z    = max(z, zwall + 1.0f-4)
            indl = interpolation4D_irregzsingle(zlev_arr, indx, indy, z, indt, Int32(1), nk)
        else
            z  = z_try
            wp = wp_try
        end

        z_bot = interpolation4Dlevsingle(zlev_arr, indx, indy, 1.0f0,       indt)
        z_top = interpolation4Dlevsingle(zlev_arr, indx, indy, Float32(nk), indt)
        z    = clamp(z, z_bot + 1.0f-4, z_top - 1.0f-4)
        indl = interpolation4D_irregzsingle(zlev_arr, indx, indy, z, indt, Int32(1), nk)
    end

    return z, delz, wp, tlw, tke_z, yl, rng_state
end

# ── Hanna ON-THE-FLY (4D, GPU-safe) ──────────────────────────────────────────
# 2026-09-21: the device-side twin of `turbulence_hanna_v1` (turbulence_hanna.jl).
# Same Langevin body as `turbulence_v3_4d` above — identical sub-stepping, the
# same reflection and clamping — with two substitutions:
#
#   • the column is sampled from the 4-D met arrays at (indx, indy, indl, indt)
#     instead of the hoisted 1-D profiles the CPU path builds once per sync step
#     (`height`, `rhoprof`, `rhogradprof`), which a kernel cannot carry; ρ and
#     ∂ρ/∂z at the half levels from half_level_density_4d (the CPU stencil);
#   • σ_w, dσ_w/dz and T_L come from the Hanna (1982) analytic profile
#     `hanna_sigw_tlw(z, pblh, ustar, wstar, ol, sigw_ft)` instead of the TKE
#     field, so the surface parameters are needed per parcel: they arrive as
#     three extra 2-D met arrays interpolated at the parcel's position.
#
# Differences from the CPU function, both deliberate:
#   • the RNG is the per-parcel xorshift stream (`xorshift_randn_clamped`), as
#     for Thomson on this path — the backends are statistically equivalent, not
#     bitwise (see BENCHMARKS.md §2);
#   • the runaway-w' guard keeps the reset but drops the atomic counter and the
#     `@warn`: neither is legal in a GPU kernel.  The event is silent here.
@inline function turbulence_hanna_4d(
    indl   ::Float32,
    indx   ::Float32,
    indy   ::Float32,
    dt     ::Float32,
    ldt    ::Float32,
    ifine  ::Float32,
    wp     ::Float32,
    indt   ::Float32,
    nk     ::Int32,
    rng_state::UInt64,
    zlev_arr   ::AbstractArray{Float32,4},
    rho_arr    ::AbstractArray{Float32,4},
    ustar  ::Float32,
    pblh   ::Float32,
    wstar  ::Float32,
    ol     ::Float32,
    sigw_ft::Float32,
)
    dt_left ::Float32 = dt
    yl      ::Float32 = 0.0f0
    tlw     ::Float32 = 1.0f0
    sigw    ::Float32 = 0.2f0
    tke_z   ::Float32 = 0.0f0
    z       ::Float32 = interpolation4Dlevsingle(zlev_arr, indx, indy, indl, indt)
    delz    ::Float32 = 0.0f0
    zwall   ::Float32 = 0.0f0

    while dt_left > eps(Float32)
        z       = interpolation4Dlevsingle(zlev_arr,   indx, indy, indl, indt)
        rhoa, rhograd = half_level_density_4d(rho_arr, zlev_arr, indx, indy, indl, indt)
        rhoaux  = rhograd / max(rhoa, 1.0f-12)

        sigw, dsigwdz, tlw = hanna_sigw_tlw(z, pblh, ustar, wstar, ol, sigw_ft)
        yl    = 0.5f0 * tlw * sigw
        tke_z = 0.5f0 * sigw^2

        dt_sub ::Float32 = max(min(tlw / 5.0f0 / ifine, dt / ifine), 0.5f0)
        dt_sub = min(dt_sub, dt_left)
        dt_left -= dt_sub

        dtftlw ::Float32 = dt_sub / max(tlw, 0.5f0)
        rw     ::Float32 = exp(-dtftlw)

        rng_state, zrand = xorshift_randn_clamped(rng_state)

        term1  ::Float32 = rw * wp
        term2  ::Float32 = zrand * sqrt(max(1.0f0 - rw * rw, 0.0f0))
        term3  ::Float32 = tlw * (1.0f0 - rw) * (dsigwdz + rhoaux * sigw)

        wp_try ::Float32 = term1 + term2 + term3
        if !(isfinite(wp_try) && abs(wp_try) <= 20.0f0)
            rng_state, zr2 = xorshift_randn_clamped(rng_state)
            wp_try = zr2
        end
        delz            = wp_try * sigw * dt_sub
        z_try  ::Float32 = z + delz

        if z_try < zwall
            z    = 2.0f0 * zwall - z_try
            wp   = -wp_try
            z    = max(z, zwall + 1.0f-4)
            indl = interpolation4D_irregzsingle(zlev_arr, indx, indy, z, indt, Int32(1), nk)
        else
            z  = z_try
            wp = wp_try
        end

        z_bot = interpolation4Dlevsingle(zlev_arr, indx, indy, 1.0f0,       indt)
        z_top = interpolation4Dlevsingle(zlev_arr, indx, indy, Float32(nk), indt)
        z    = clamp(z, z_bot + 1.0f-4, z_top - 1.0f-4)
        indl = interpolation4D_irregzsingle(zlev_arr, indx, indy, z, indt, Int32(1), nk)
    end

    return z, delz, wp, tlw, tke_z, yl, rng_state
end
