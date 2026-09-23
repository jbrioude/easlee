
# ── KernelAbstractions advance kernel ─────────────────────────────────────────
# advance_parcel_kernel! processes one parcel per work-item (GPU thread / CPU
# task).  It is identical on CPU (CPU() backend) and CUDA (CUDABackend()).
#
# All Dict/struct lookups are done BEFORE the kernel launch in
# advance_active_parcels_ka! so that only primitive types and plain arrays enter
# the kernel — the only types that are GPU-safe.
#
# Gridconfig is passed as a Rectilinear/Lambert value (isbitstype → GPU-safe).
# Levelconfig only needs Nk::Int32 (SigmaPressureLevel.ak/bk are not used here).

using KernelAbstractions
using KernelAbstractions.Extras

# ── Boundary-condition helpers (GPU-safe, no dispatch on Any) ─────────────────

@inline function _dxtodindx_rect(g::Rectilinear, dx::Float32,
                                  indx::Float32, indy::Float32)
    lat = (indy - 1.0f0) * g.dy + g.latmin
    return dx * g.dxconst / cos(lat * Float32(π) / 180.0f0)
end

@inline _dytodindy_rect(g::Rectilinear, dy::Float32) = dy * g.dyconst

@inline function _dxtodindx_lambert(g::Lambert, dx::Float32)
    return dx / g.dx
end

@inline _dytodindy_lambert(g::Lambert, dy::Float32) = dy / g.dy

# ── mesoscale sigma over the 2x2x2x2 met neighbourhood (GPU-safe twin of
# _std16 / _mesoscale_sigmas in advance.jl) ──────────────────────────────────
@inline function _std16_4d(a, i0::Int, j0::Int, k0::Int)
    s = 0.0f0; s2 = 0.0f0
    @inbounds for l in 1:2, k in k0:(k0 + 1), j in j0:(j0 + 1), i in i0:(i0 + 1)
        v = a[i, j, k, l]
        isfinite(v) || return 0.0f0
        s += v; s2 += v * v
    end
    m = s / 16.0f0
    var = s2 / 16.0f0 - m * m
    return var > 0.0f0 ? sqrt(var) : 0.0f0
end

@inline function _mesoscale_sigmas_4d(uarr, varr, warr, xi::Float32, yi::Float32,
                                      zi::Float32, nz_int::Int)
    i0 = clamp(Base.unsafe_trunc(Int, xi), 1, size(uarr, 1) - 1)
    j0 = clamp(Base.unsafe_trunc(Int, yi), 1, size(uarr, 2) - 1)
    k0 = clamp(Base.unsafe_trunc(Int, zi), 1, nz_int - 1)
    return (_std16_4d(uarr, i0, j0, k0), _std16_4d(varr, i0, j0, k0), _std16_4d(warr, i0, j0, k0))
end

# ── Inner parcel body (generic, inlined into the kernel) ──────────────────────

@inline function _advance_one!(
    # parcel state (scalars for the current parcel)
    indx_r::Float32, indy_r::Float32, indl_r::Float32,
    up::Float32, vp::Float32, wp::Float32, ldt::Float32,
    rng_st::UInt64,
    # met arrays (concrete typed)
    uarray      ::AbstractArray{Float32,4},
    varray      ::AbstractArray{Float32,4},
    levdotarray ::AbstractArray{Float32,4},
    zlev_arr    ::AbstractArray{Float32,4},
    theta_arr   ::AbstractArray{Float32,4},
    tke_arr     ::AbstractArray{Float32,4},
    rho_arr    ::AbstractArray{Float32,4},
    drhodz_arr ::AbstractArray{Float32,4},
    # 2026-09-21: surface parameters for on-the-fly Hanna (length-0 when off)
    ustar_arr  ::AbstractArray{Float32,3},
    pblh_arr   ::AbstractArray{Float32,3},
    ol_arr     ::AbstractArray{Float32,3},
    use_hanna  ::Bool,
    sigw_ft    ::Float32,
    # 2026-09-21: the last five threads-only switches (Jerome: a GPU run must be
    # GPU everywhere).  mesou/mesov/mesow are per-parcel OU state, in and out.
    mesou      ::Float32,
    mesov      ::Float32,
    mesow      ::Float32,
    use_meso   ::Bool,
    fturbmeso  ::Float32,
    r_meso     ::Float32,
    rs_meso    ::Float32,
    kh_ft      ::Float32,      # 0 = off
    use_htke   ::Bool,
    use_hhanna ::Bool,
    use_mo     ::Bool,
    # grid/physics config
    gridconfig,          # Rectilinear or Lambert (isbitstype, dispatches statically)
    turb,                # Thomson, Hanna, or NoTurb — no annotation so both dispatch statically
    vmix       ::VMixClosure,   # per-driver σ_w²/K_z closure (Fix C, isbits)
    indt_n     ::Float32,
    nk         ::Int32,
    nsub       ::Int32,     # sub-steps per sync step (host: _nsubsteps, advance.jl)
    dt_inc     ::Float32,
    ifine      ::Float32,
    t_direction::Float32,
    swn_j      ::Float32,   # polar caps (Inf32 / -Inf32 = off): indy beyond
    sws_j      ::Float32,   # these rows → stereo-plane position update
)
    # Polar-stereographic cap membership — fixed for the whole step (winds are
    # sampled at the entry position; see the CPU kernel for the rationale).
    hemi ::Float32 = indy_r > swn_j ? 1.0f0 : (indy_r < sws_j ? -1.0f0 : 0.0f0)

    # 2026-09-22 double-check: the CPU samples the horizontal-Hanna surface
    # parameters ONCE per sync step at the ENTRY position (advance.jl:187-195),
    # hoisted above the sub-step loop.  Do the same, so the two paths sample the
    # same air for sigma_u/sigma_v.
    hh_ust ::Float32 = 0.0f0; hh_pbl ::Float32 = 0.0f0
    hh_ol  ::Float32 = 0.0f0; hh_wst ::Float32 = 0.0f0
    if use_hhanna
        ix_h ::Float32 = clamp(indx_r, 1.0f0, Float32(size(zlev_arr, 1)) - 1.0f0)
        iy_h ::Float32 = clamp(indy_r, 1.0f0, Float32(size(zlev_arr, 2)) - 1.0f0)
        hh_ust = interpolate2Dsingle(ustar_arr, ix_h, iy_h, indt_n)
        hh_pbl = interpolate2Dsingle(pblh_arr,  ix_h, iy_h, indt_n)
        hh_ol  = interpolate2Dsingle(ol_arr,    ix_h, iy_h, indt_n)
        hh_wst = hanna_wstar(hh_ust, hh_pbl, hh_ol)
    end

    dxsave    ::Float32 = 0.0f0
    dysave    ::Float32 = 0.0f0
    dindl     ::Float32 = 0.0f0
    uout      ::Float32 = 0.0f0
    vout      ::Float32 = 0.0f0
    wout      ::Float32 = 0.0f0
    tke_z     ::Float32 = 0.0f0
    yl        ::Float32 = 0.0f0
    totdelz   ::Float32 = 0.0f0

    # review 2026-09-04 A5: integer sub-step count from the host (_nsubsteps,
    # advance.jl) — the Float32 `while inc_time < float_timestep; inc_time += dt_inc`
    # accumulation ran ctl+1 sub-steps for ~20 % of (Δt, ctl) pairs (Δt=3 s,
    # ctl=9 → 10).  dt_inc = Δt/nsub exactly; always the full sub-step duration.
    for _ in 1:nsub
        dt ::Float32 = dt_inc

        ux     ::Float32 = interpolateU(gridconfig.Uedge, uarray,      indx_r, indy_r, indl_r, indt_n)
        vy     ::Float32 = interpolateV(gridconfig.Vedge, varray,      indx_r, indy_r, indl_r, indt_n)
        levdot ::Float32 = interpolateW(gridconfig.Wedge, levdotarray, indx_r, indy_r, indl_r, indt_n)

        # [physics].surface_layer_wind — below the lowest model level the archived
        # wind is a slab; put the parcel back on the Monin-Obukhov similarity
        # profile anchored at that level (mirrors advance.jl:261-277).
        if use_mo && indl_r < 2.0f0
            z_sl  ::Float32 = interpolation4Dlevsingle(zlev_arr, indx_r, indy_r, indl_r, indt_n)
            z1_sl ::Float32 = interpolation4Dlevsingle(zlev_arr, indx_r, indy_r, 2.0f0,  indt_n)
            if isfinite(z_sl) && z_sl < z1_sl
                ust_sl ::Float32 = interpolate2Dsingle(ustar_arr, indx_r, indy_r, indt_n)
                ol_sl  ::Float32 = interpolate2Dsingle(ol_arr,    indx_r, indy_r, indt_n)
                fac_sl ::Float32 = mo_wind_factor(z_sl, z1_sl, ust_sl, ol_sl,
                                                  sqrt(ux * ux + vy * vy))
                ux *= fac_sl
                vy *= fac_sl
            end
        end

        # review 2026-09-04 A4: masked/terrain NaN in the met → zero displacement
        # (CPU kernel advance.jl:243-250), never a NaN position / unsafe_trunc crash.
        dxsave += (isfinite(ux) ? ux : 0.0f0) * dt * t_direction
        dysave += (isfinite(vy) ? vy : 0.0f0) * dt * t_direction

        # review 2026-09-04 A1: PER-SUB-STEP vertical displacement, not the
        # accumulated one.  indl_r is committed every sub-step below, so adding the
        # running total again applied the mean vertical wind (ctl+1)/2 times
        # (3x at ctl=5).  Mirrors advance.jl:245-254 (dlev_step).
        # For SigmaPressureLevel: indltolev and levtoindl are identity → dindl = dlev_step
        dlev_step ::Float32 = (isfinite(levdot) ? levdot : 0.0f0) * dt * t_direction
        dindl  = dlev_step

        prov = hemi == 0.0f0 ?
            (indy_r + _dytodindy(gridconfig, dysave, indx_r, indy_r),
             indx_r + _dxtodindx(gridconfig, dxsave, indx_r, indy_r)) :
            polar_stereo_update(gridconfig, hemi, indy_r, indx_r, dxsave, dysave)
        indy_prov ::Float32 = prov[1]
        indx_prov ::Float32 = prov[2]
        indl_prov ::Float32 = indl_r + dindl

        if turb isa NoTurb
            # review 2026-09-04 A4: pure advection mirrors advance.jl:268-281 —
            # no z↔indl round trip through an unclamped indl_prov (indl < 1 read
            # level 0 out of bounds); reflect at the ground, clamp at the top.
            new_indl = indl_prov
            if new_indl <= 1.0f0
                indl_r = abs(new_indl - 1.0f0) + 1.0f0
            else
                indl_r = min(new_indl, Float32(nk) - 1.0f-3)
            end
        else
            # review 2026-09-04 A4: clamp the provisional position to the valid
            # interpolation range before the turbulence step (advance.jl:285-286);
            # the 4-D samplers use unsafe_trunc and would read index 0 / Nk+1.
            indx_safe ::Float32 = clamp(indx_prov, 1.0f0, Float32(size(zlev_arr, 1)) - 1.0f0)
            indy_safe ::Float32 = clamp(indy_prov, 1.0f0, Float32(size(zlev_arr, 2)) - 1.0f0)
            indl_safe ::Float32 = clamp(indl_prov, 1.0f0, Float32(nk) - 1.0f-3)

            # 2026-09-21: on-the-fly Hanna samples the surface parameters at the
            # PROVISIONAL horizontal position, as the CPU path does
            # (advance.jl:317-338), then runs the same Langevin body with the
            # Hanna analytic σ_w / T_L instead of the TKE profile.
            z_new ::Float32, delz ::Float32, wp, ldt, tke_z, yl, rng_st =
                if Bool(use_hanna)
                    ust_h  ::Float32 = interpolate2Dsingle(ustar_arr, indx_safe, indy_safe, indt_n)
                    pbl_h  ::Float32 = interpolate2Dsingle(pblh_arr,  indx_safe, indy_safe, indt_n)
                    ol_h   ::Float32 = interpolate2Dsingle(ol_arr,    indx_safe, indy_safe, indt_n)
                    wst_h  ::Float32 = hanna_wstar(ust_h, pbl_h, ol_h)
                    turbulence_hanna_4d(indl_safe, indx_safe, indy_safe,
                                        dt, ldt, ifine, wp, indt_n, nk, rng_st,
                                        zlev_arr, rho_arr, drhodz_arr,
                                        ust_h, pbl_h, wst_h, ol_h, Float32(sigw_ft))
                else
                    turbulence_v3_4d(turb, indl_safe, indx_safe, indy_safe,
                                     dt, ldt, ifine, wp, indt_n, nk, rng_st,
                                     zlev_arr, theta_arr, tke_arr,
                                     rho_arr, drhodz_arr, vmix)
                end

            totdelz += delz

            # ── horizontal turbulence (advance.jl:375-400) ────────────────────
            # exact OU update over the sub-step; up/vp persist across sub-steps
            # and across sync steps through uturb/vturb, so they are returned.
            if use_htke
                tlh ::Float32, sigh ::Float32 = tke_partition_h(vmix, tke_z, yl)
                if sigh > 0.0f0
                    rr ::Float32 = exp(-dt / tlh)
                    sq ::Float32 = sqrt(max(1.0f0 - rr * rr, 0.0f0))
                    rng_st, zu = xorshift_randn_clamped(rng_st)
                    rng_st, zv = xorshift_randn_clamped(rng_st)
                    up = up * rr + sigh * sq * zu
                    vp = vp * rr + sigh * sq * zv
                    dxsave += up * dt
                    dysave += vp * dt
                end
            elseif use_hhanna
                # Hanna sigma_u/T_Lu on grid-x, sigma_v/T_Lv on grid-y, no rotation
                # into the wind frame (FLEXPART advance.f90:392-401); surface
                # parameters are the ones hoisted at entry, as on the CPU.
                sigu, tlu, sigv, tlv = hanna_h(max(z_new, 0.0f0), hh_pbl, hh_ol, hh_ust, hh_wst)
                ru ::Float32 = exp(-dt / tlu)
                rv ::Float32 = exp(-dt / tlv)
                rng_st, zu2 = xorshift_randn_clamped(rng_st)
                rng_st, zv2 = xorshift_randn_clamped(rng_st)
                up = up * ru + sigu * sqrt(max(1.0f0 - ru * ru, 0.0f0)) * zu2
                vp = vp * rv + sigv * sqrt(max(1.0f0 - rv * rv, 0.0f0)) * zv2
                dxsave += up * dt
                dysave += vp * dt
            end

            if z_new > 0.0f0
                dindl = interpolation4D_irregzsingle(zlev_arr, indx_r, indy_r,
                                                     z_new, indt_n, Int32(1), nk) - indl_r
            else
                dindl = 0.0f0
            end

            new_indl = indl_r + dindl
            if new_indl <= 1.0f0
                indl_r = abs(new_indl - 1.0f0) + 1.0f0
            else
                indl_r += dindl
            end
        end

        # match the CPU kernel (advance.jl): u/v turbulence is not carried across
        # steps UNLESS a horizontal scheme owns them (then uturb/vturb are the
        # OU memory and must survive — advance.jl:384/398 set uout = up)
        uout = (use_htke || use_hhanna) ? up : 0.0f0
        vout = (use_htke || use_hhanna) ? vp : 0.0f0
        wout = wp
    end

    # ── mesoscale wind fluctuations (advance.jl:415-425) ─────────────────────
    # OU velocity with amplitude fturbmeso*sigma over the 2x2x2x2 neighbourhood;
    # acts in ALL layers, vertical in level-index space with ground reflection.
    ts_full ::Float32 = dt_inc * Float32(nsub)
    if use_meso
        usig, vsig, wsig = _mesoscale_sigmas_4d(uarray, varray, levdotarray,
                                                indx_r, indy_r, indl_r, Int(nk))
        if mesou == 0.0f0 && mesov == 0.0f0 && mesow == 0.0f0
            rng_st, z1 = xorshift_randn_clamped(rng_st); mesou = z1 * usig * fturbmeso
            rng_st, z2 = xorshift_randn_clamped(rng_st); mesov = z2 * vsig * fturbmeso
            rng_st, z3 = xorshift_randn_clamped(rng_st); mesow = z3 * wsig * fturbmeso
        else
            rng_st, z1 = xorshift_randn_clamped(rng_st); mesou = r_meso * mesou + rs_meso * z1 * usig * fturbmeso
            rng_st, z2 = xorshift_randn_clamped(rng_st); mesov = r_meso * mesov + rs_meso * z2 * vsig * fturbmeso
            rng_st, z3 = xorshift_randn_clamped(rng_st); mesow = r_meso * mesow + rs_meso * z3 * wsig * fturbmeso
        end
        dxsave += mesou * ts_full * t_direction
        dysave += mesov * ts_full * t_direction
        ml_new ::Float32 = indl_r + mesow * ts_full * t_direction
        ml_new <= 1.0f0 && (ml_new = abs(ml_new - 1.0f0) + 1.0f0)
        indl_r = min(ml_new, Float32(nk) - 1.0f-3)
    end

    # ── free-atmosphere horizontal diffusivity (advance.jl:439-450) ──────────
    # memoryless random walk above the PBL top: delta = sqrt(2*K_h*dt)*N(0,1)
    if kh_ft > 0.0f0
        ix_s ::Float32 = clamp(indx_r, 1.0f0, Float32(size(zlev_arr, 1)) - 1.0f0)
        iy_s ::Float32 = clamp(indy_r, 1.0f0, Float32(size(zlev_arr, 2)) - 1.0f0)
        il_s ::Float32 = clamp(indl_r, 1.0f0, Float32(nk) - 1.0f-3)
        zpar ::Float32 = interpolation4Dlevsingle(zlev_arr, ix_s, iy_s, il_s, indt_n)
        pbl2 ::Float32 = interpolate2Dsingle(pblh_arr, ix_s, iy_s, indt_n)
        if isfinite(zpar) && zpar > max(pbl2, 50.0f0)
            sk ::Float32 = sqrt(2.0f0 * kh_ft * ts_full)
            rng_st, zk1 = xorshift_randn_clamped(rng_st); dxsave += sk * zk1
            rng_st, zk2 = xorshift_randn_clamped(rng_st); dysave += sk * zk2
        end
    end

    if hemi == 0.0f0
        indx_r += _dxtodindx(gridconfig, dxsave, indx_r, indy_r)
        indy_r += _dytodindy(gridconfig, dysave, indx_r, indy_r)
    else
        # Stereo-plane update — exact across the pole (see CPU kernel); nudge
        # off the exact pole row so the interpolation stencil stays in bounds.
        fin = polar_stereo_update(gridconfig, hemi, indy_r, indx_r, dxsave, dysave)
        indy_r = clamp(fin[1], 1.0f0, Float32(size(uarray, 2)) - 1.0f-3)
        indx_r = fin[2]
    end

    return indx_r, indy_r, indl_r, uout, vout, wout, ldt, tke_z, yl, totdelz, rng_st,
           mesou, mesov, mesow
end

# Dispatch wrappers so the compiler can monomorphise _advance_one! on grid type.
@inline _dxtodindx(g::Rectilinear, dx, indx, indy) =
    dx * g.dxconst / cos(((indy - 1.0f0) * g.dy + g.latmin) * Float32(π) / 180.0f0)
@inline _dxtodindx(g::Lambert, dx, indx, indy) = dx / g.dx
@inline _dytodindy(g::Rectilinear, dy, indx, indy) = dy * g.dyconst
@inline _dytodindy(g::Lambert, dy, indx, indy)     = dy / g.dy
# Mercator (Meso-NH): delegate to the GPU-safe coords formula (cosd/atan/exp/log).
@inline _dxtodindx(g::Mercator, dx, indx, indy) = dxtodindx(g, dx, indx, indy)
@inline _dytodindy(g::Mercator, dy, indx, indy) = dytodindy(g, dy, indx, indy)

# ── KernelAbstractions kernel ─────────────────────────────────────────────────

# 2026-09-23: physics switches reach the kernel as Val{Bool} so the GPU compiler
# specialises on their VALUE and dead-code-eliminates the branches that are off.
# Passing them as runtime Bools compiled every branch into one kernel and cost
# 44 % of the GPU advance with all five switches off (bisect 4345524 → 2bdf8f8:
# 116 → 167 s on the ETEX 1 M benchmark).  _kb unwraps either form.
@inline _kb(b::Bool) = b
@inline _kb(::Val{B}) where {B} = B

@kernel function advance_parcel_kernel!(
    # active parcel index list
    @Const(idx_arr),
    # parcel state arrays (read+write)
    indx_arr, indy_arr, indl_arr,
    uturb_arr, vturb_arr, wturb_arr, dt_arr,
    miss2_arr, miss3_arr, miss4_arr,
    Dtime_ns_arr,          # Int64 — nanoseconds since epoch
    rng_state_arr,         # UInt64 xorshift state per parcel
    # met arrays (read-only)
    @Const(uarray), @Const(varray), @Const(levdotarray),
    @Const(zlev_arr), @Const(theta_arr), @Const(tke_arr),
    @Const(rho_arr), @Const(drhodz_arr),
    # 2026-09-21: surface parameters for ON-THE-FLY Hanna (3-D: Ni x Nj x 2 slots).
    # Zero-length arrays when use_hanna is false — KA still needs the argument.
    @Const(ustar_arr), @Const(pblh_arr), @Const(ol_arr),
    # 2026-09-21: mesoscale OU state (read+write)
    mesou_arr, mesov_arr, mesow_arr,
    # grid / level config + turbulence tag (isbitstype structs, no type annotations
    # in the @kernel signature — KA 0.9 GPU compiler specialises on runtime type)
    gridconfig,
    BL89,            # Bool — Thomson.BL89
    read_tls,        # Bool — Thomson.read_tls
    use_turb,        # Bool — false → NoTurb (no turbulence), true → Thomson
    use_hanna,       # Bool — true → Hanna on-the-fly (overrides the Thomson body)
    sigw_ft,         # Float32 — free-troposphere σ_w floor (Hanna.sigw_ft)
    use_meso,        # Bool    — [physics.turbulence].mesoscale
    fturbmeso,       # Float32 — mesoscale amplitude
    r_meso,          # Float32 — OU memory exp(-2Δt/T_met)
    rs_meso,         # Float32 — sqrt(1-r²)
    kh_ft_v,         # Float32 — free-atmosphere K_h (0 = off)
    kh_on,           # Val{Bool} — kh_ft_v > 0, so the branch folds at compile time
    use_htke,        # Bool    — horizontal_tke
    use_hhanna,      # Bool    — horizontal_hanna
    use_mo,          # Bool    — surface_layer_wind (Monin-Obukhov)
    vmix,            # VMixClosure — per-driver σ_w²/K_z closure (Fix C, isbits)
    indt_n,
    nk,
    nsub,            # Int32 — sub-steps per sync step (review 2026-09-04 A5)
    dt_inc,
    ifine,
    t_direction,
    swn_j,           # polar-cap switch rows (Inf32 / -Inf32 = off)
    sws_j,
    timestep_ns,
    Ni,
    Nj,
    globalx,
    dlon,
)
    j = @index(Global, Linear)
    # KA kernels cannot have bare return — guard with if.
    # review 2026-09-04 A4: entry test of the CPU kernel (advance.jl:176-178).
    # A parcel released outside the domain has OOB indx/indy; the unsafe_trunc
    # samplers would read index 0 / > Ni.  Skip it and write NOTHING: Dtime_ns
    # stays at currentdate, so reassign_escaped_parcels! frees the slot exactly
    # as after the CPU `continue`.  On a global grid [Ni, Ni+1) is the valid
    # antimeridian seam (halo column) → allow up to Ni+1.
    if j <= length(idx_arr) &&
       !(indx_arr[Int(idx_arr[j])] < 1.0f0 ||
         indx_arr[Int(idx_arr[j])] >= (globalx ? Float32(Ni) + 1.0f0 : Float32(Ni)) ||
         indy_arr[Int(idx_arr[j])] < 1.0f0 ||
         indy_arr[Int(idx_arr[j])] >= Float32(Nj))
        i = Int(idx_arr[j])

        indx_r  = indx_arr[i]
        indy_r  = indy_arr[i]
        indl_r  = indl_arr[i]
        up      = uturb_arr[i]
        vp      = vturb_arr[i]
        wp      = wturb_arr[i]
        ldt     = dt_arr[i]
        rng_st  = rng_state_arr[i]

        # Dispatch on turbulence flag with concrete types in each branch so the
        # GPU compiler can resolve the turbulence_v3_4d method statically.
        kh_eff ::Float32 = _kb(kh_on) ? Float32(kh_ft_v) : 0.0f0
        if _kb(use_turb)
            indx_r, indy_r, indl_r, uout, vout, wout, ldt, tke_z, yl, totdelz, rng_st,
            mesou_n, mesov_n, mesow_n =
                _advance_one!(indx_r, indy_r, indl_r, up, vp, wp, ldt, rng_st,
                              uarray, varray, levdotarray,
                              zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
                              ustar_arr, pblh_arr, ol_arr, _kb(use_hanna), Float32(sigw_ft),
                              mesou_arr[i], mesov_arr[i], mesow_arr[i],
                              _kb(use_meso), Float32(fturbmeso), Float32(r_meso), Float32(rs_meso),
                              kh_eff, _kb(use_htke), _kb(use_hhanna), _kb(use_mo),
                              gridconfig, Thomson(_kb(BL89), _kb(read_tls)), vmix, indt_n, nk,
                              nsub, dt_inc, ifine, t_direction, swn_j, sws_j)
        else
            indx_r, indy_r, indl_r, uout, vout, wout, ldt, tke_z, yl, totdelz, rng_st,
            mesou_n, mesov_n, mesow_n =
                _advance_one!(indx_r, indy_r, indl_r, up, vp, wp, ldt, rng_st,
                              uarray, varray, levdotarray,
                              zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
                              ustar_arr, pblh_arr, ol_arr, false, 0.1f0,
                              mesou_arr[i], mesov_arr[i], mesow_arr[i],
                              _kb(use_meso), Float32(fturbmeso), Float32(r_meso), Float32(rs_meso),
                              kh_eff, _kb(use_htke), _kb(use_hhanna), _kb(use_mo),
                              gridconfig, NoTurb(), vmix, indt_n, nk,
                              nsub, dt_inc, ifine, t_direction, swn_j, sws_j)
        end

        # Write back parcel state
        indx_arr[i]      = indx_r
        indy_arr[i]      = indy_r
        indl_arr[i]      = indl_r
        uturb_arr[i]     = uout
        vturb_arr[i]     = vout
        wturb_arr[i]     = wout
        dt_arr[i]        = ldt
        miss2_arr[i]     = tke_z
        miss3_arr[i]     = yl
        miss4_arr[i]     = totdelz
        rng_state_arr[i] = rng_st
        mesou_arr[i]     = mesou_n
        mesov_arr[i]     = mesov_n
        mesow_arr[i]     = mesow_n

        # Pole overshoot from the lat/lon path (switch row close to the pole):
        # mirror across the pole instead of killing (FLEXPART HSO/AL fix).
        if swn_j < Inf32 && indy_r > Float32(Nj)
            indy_r  = 2.0f0 * Float32(Nj) - indy_r
            indx_r += 180.0f0 / dlon
            indy_arr[i] = indy_r
            indx_arr[i] = indx_r
        elseif sws_j > -Inf32 && indy_r < 1.0f0
            indy_r  = 2.0f0 - indy_r
            indx_r += 180.0f0 / dlon
            indy_arr[i] = indy_r
            indx_arr[i] = indx_r
        end

        # Boundary conditions
        new_Dtime = Dtime_ns_arr[i]
        if globalx
            # Periodic wrap into [1, Ni+1).  Upper threshold is Ni+1, not Ni: the
            # halo column at Ni+1 makes [Ni, Ni+1) a valid antimeridian seam, so
            # wrapping at >= Ni instead pushed crossing parcels into [0,1) where
            # the next interpolation reads column 0 (the 180° void / OOB read).
            indx_r < 1.0f0                && (indx_arr[i] += 360.0f0 / dlon)
            indx_r >= Float32(Ni) + 1.0f0 && (indx_arr[i] -= 360.0f0 / dlon)
        else
            if indx_r < 1.0f0 || indx_r >= Float32(Ni)
                new_Dtime = Int64(1)
            end
        end
        if indy_r < 1.0f0 || indy_r >= Float32(Nj)
            new_Dtime = Int64(2)
        end
        # Advance parcel time if still alive
        if new_Dtime > Int64(2)
            new_Dtime += timestep_ns
        end
        Dtime_ns_arr[i] = new_Dtime
    end
end

# ── New-parcel initialisation kernel ─────────────────────────────────────────
# GPU-safe equivalent of initialize_new_parcels!.
# Each work-item initialises one newly-released parcel:
#   • samples random release height via xorshift_uniform01
#   • converts z → indl via interpolation4D_irregzsingle (no alloc)
#   • computes tke/yl/tlw via turbulent_length_scale_v2_4d (no profile buffer)
#   • samples random wturb via xorshift_randn_clamped
#
# idx_arr  — indices of new parcels (CPU findall result, passed as CuArray on GPU)
# zrelease1/2_arr — per-box release height bounds (small, upload once at run start)

@kernel function initialize_parcels_kernel!(
    # parcel state (read+write)
    indl_arr,
    wturb_arr,
    dt_arr,
    miss2_arr,
    miss3_arr,
    rng_state_arr,
    # parcel state (read-only)
    @Const(idx_arr),
    @Const(indx_arr),
    @Const(indy_arr),
    @Const(relbox_arr),
    # release box heights (read-only, nboxes elements)
    @Const(zrelease1_arr),
    @Const(zrelease2_arr),
    # met arrays (read-only)
    @Const(zlev_arr),
    @Const(theta_arr),
    @Const(tke_arr),
    # config scalars
    vmix   ::VMixClosure,   # per-driver σ_w²/K_z closure (Fix C, isbits)
    indt_n ::Float32,
    nk     ::Int32,
    ctl    ::Float32,
    Ni     ::Int32,
    Nj     ::Int32,
)
    j = @index(Global, Linear)
    if j <= length(idx_arr)
        i = Int(idx_arr[j])

        indx = indx_arr[i]
        indy = indy_arr[i]
        ibox = Int(relbox_arr[i])

        z1 = zrelease1_arr[ibox]
        z2 = zrelease2_arr[ibox]

        rng_st = rng_state_arr[i]

        # Sample random release height in [z1, z2]
        rng_st, u = xorshift_uniform01(rng_st)
        z_release = z1 + u * (z2 - z1)

        # Clamp to valid grid range
        safe_indx = clamp(indx, 1.0f0, Float32(Ni - 1))
        safe_indy = clamp(indy, 1.0f0, Float32(Nj - 1))

        # z → fractional level index (GPU-safe, no alloc)
        indl = interpolation4D_irregzsingle(zlev_arr, safe_indx, safe_indy,
                                             z_release, indt_n, Int32(1), nk)
        indl = clamp(indl, 1.0f0, Float32(nk - 1))

        # Turbulence quantities at release level
        tke_z  = interpolation4Dlevsingle(tke_arr,  safe_indx, safe_indy, indl, indt_n)
        z_here = interpolation4Dlevsingle(zlev_arr, safe_indx, safe_indy, indl, indt_n)

        indz  = clamp(Base.unsafe_trunc(Int32, indl), Int32(1), nk - Int32(1))
        indzp = indz + Int32(1)

        yl = turbulent_length_scale_v2_4d(indz, indzp, tke_z, z_here,
                                           safe_indx, safe_indy, indt_n, nk,
                                           zlev_arr, theta_arr)
        tlw, _ = tke_partition(vmix, tke_z, yl)

        # Random initial vertical velocity (normalised)
        rng_st, wturb = xorshift_randn_clamped(rng_st)

        # Write back parcel state
        indl_arr[i]      = indl
        wturb_arr[i]     = wturb
        dt_arr[i]        = min(tlw / ctl / 5.0f0, 10.0f0)
        miss2_arr[i]     = tke_z
        miss3_arr[i]     = yl
        rng_state_arr[i] = rng_st
    end
end

# ── GPU-dispatch wrapper for new-parcel initialisation ────────────────────────
# Arguments:
#   backend         — KA backend (CUDABackend() or CPU())
#   gs              — device parcel-state NamedTuple
#   gpu_met         — MetFields3D with device arrays (CuArray on GPU, Array on CPU)
#   new_idx         — device array of new parcel indices (CuArray or Vector{Int})
#   zrelease1/2_dev — device vectors of release heights (upload once at run start)
#   indt_n          — Float32 met-time index
#   nk, ctl, Ni, Nj — Int32/Float32 config scalars

"""
    initialize_new_parcels_gpu!(backend, gs, gpu_met, new_idx,
                                zrelease1_dev, zrelease2_dev,
                                indt_n, nk, ctl, Ni, Nj)

GPU-safe equivalent of `initialize_new_parcels!`: samples each new parcel's
release height [m] in its box's [z1, z2], converts to a fractional level
index, and seeds the turbulence state (tke/ℓ/T_L, random wturb) from the
per-parcel xorshift RNG.  See the argument list in the comment above.
"""
function initialize_new_parcels_gpu!(
    backend,
    gs,
    gpu_met,
    new_idx,
    zrelease1_dev,
    zrelease2_dev,
    indt_n ::Float32,
    nk     ::Int32,
    ctl    ::Float32,
    Ni     ::Int32,
    Nj     ::Int32,
)
    isempty(new_idx) && return nothing

    kernel! = initialize_parcels_kernel!(backend, 256)
    kernel!(
        gs.indl, gs.wturb, gs.dt, gs.miss2, gs.miss3, gs.rng_state,
        new_idx,
        gs.indx, gs.indy, gs.relbox,
        zrelease1_dev, zrelease2_dev,
        gpu_met.zlev, gpu_met.theta, gpu_met.tke,
        VMIX[],                       # closure read host-side, passed as isbits
        indt_n, nk, ctl, Ni, Nj;
        ndrange = length(new_idx),
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end

# ── CPU/GPU dispatch wrapper ───────────────────────────────────────────────────

"""
    advance_active_parcels_ka!(state, ctx, currentdate, indt_array; backend=nothing)

KernelAbstractions advance step for all parcels due at `currentdate`
(EASLEE_ADVANCE=ka_cpu/ka_gpu path).  `backend=nothing` → threaded CPU();
pass `CUDABackend()` for GPU (met + state uploaded to resident device
buffers each step).  Turbulence: Thomson, Hanna(false) (pre-computed fields,
same Thomson body as the CPU) or NoTurb — Hanna on-the-fly (turbswitch=true)
is CPU-only (advance.jl); advance_dispatch! forces the threads path for it.
"""
function advance_active_parcels_ka!(state, ctx, currentdate, indt_array;
                                     backend=nothing)
    timestep = ctx.runconfig.timestep
    nnest    = length(ctx.nests)

    # Default: CPU() with nthreads block size
    ka_backend = isnothing(backend) ? CPU(; static=true) : backend
    block_size = ka_backend isa CPU ? min(256, Threads.nthreads() * 4) : 256

    for n in 1:nnest
        n8 = Int8(n)
        idx = findall(i -> state.Dtime[i] == currentdate && state.nest[i] == n8,
                      eachindex(state.Dtime))
        isempty(idx) && continue

        nestctx    = ctx.nests[n]
        nrc        = nestctx.nestrunconfig   # typed NestRunConfig

        # Extract all typed met arrays before kernel launch
        met         = nestctx.metfields3D
        met2        = nestctx.metfields2D     # 2026-09-21: surface params for on-the-fly Hanna
        uarray      = met.Ux
        varray      = met.Vy
        levdotarray = met.levdot
        zlev_arr    = met.zlev
        theta_arr   = met.theta
        tke_arr     = met.tke
        rho_arr     = met.rho
        drhodz_arr = met.drhodz

        gridconfig  = nestctx.gridconfig   # Rectilinear or Lambert (isbitstype)
        nk          = Int32(nestctx.levelconfig.Nk)

        float_timestep = Float32(abs(Dates.value(timestep)) / 1000.0f0)
        # review 2026-09-04 A5: integer sub-step count, same helper as the CPU
        nsub_i, dt_inc = _nsubsteps(float_timestep, nrc.ctl)
        nsub   = Int32(nsub_i)
        ifine  = 5.0f0
        t_dir  = Float32(nrc.direction)
        indt_n = indt_array[n]

        turb_raw = nrc.turbulence
        # KA kernel supports Thomson (full turbulence) and NoTurb (pure advection).
        # review 2026-09-04 A2: Hanna(false) (pre-computed σ_w/T_L written into the
        # tke field at met-read time — the ERA5/GFS/MERRA-2/WRF default) runs the
        # Thomson Langevin body on the CPU with the dummy tag Thomson(true,false)
        # (turbulence.jl:259-268); do the same here.  Only on-the-fly Hanna
        # (turbswitch=true) is CPU-only — advance_dispatch! forces :threads for it.
        # 2026-09-21: on-the-fly Hanna now runs on this path too (it used to force
        # :threads through _ka_unsupported_switch).  use_turb still gates the
        # Thomson body; use_hanna overrides it inside the kernel.
        use_hanna = turb_raw isa Hanna && turb_raw.turbswitch
        sigw_ft   = use_hanna ? Float32(turb_raw.sigw_ft) : 0.1f0
        # 2026-09-21: the remaining four switches + Monin-Obukhov, all on this path now
        use_meso   = Bool(nrc.mesoscale)
        fturbmeso  = Float32(nrc.fturbmeso)
        ts_f       = Float32(abs(ctx.runconfig.timestep.value) / 1000.0)
        r_meso     = exp(-2.0f0 * ts_f / max(Float32(nrc.windint_s), ts_f))
        rs_meso    = sqrt(max(1.0f0 - r_meso * r_meso, 0.0f0))
        kh_ft_v    = Float32(nrc.kh_ft)
        use_htke   = Bool(nrc.horiz_tke)
        use_hhanna = Bool(nrc.horiz_hanna)
        use_mo     = Bool(nrc.mo_surface)
        # the surface arrays are needed by Hanna, horizontal Hanna, K_h aloft and MO
        need_sfc   = use_hanna || use_hhanna || use_mo || kh_ft_v > 0.0f0
        use_turb = turb_raw isa Thomson || turb_raw isa Hanna
        turb_tag = turb_raw isa Thomson ? turb_raw : Thomson(true, false)
        turb_BL89     = use_turb ? Bool(turb_tag.BL89)     : false
        turb_read_tls = use_turb ? Bool(turb_tag.read_tls) : false

        Ni = Int32(nestctx.gridconfig.Ni)
        Nj = Int32(nestctx.gridconfig.Nj)
        globalx = Bool(nestctx.gridconfig.globalx)
        # dlon only exists on Rectilinear; Lambert domains are always regional (globalx=false)
        dlon    = gridconfig isa Rectilinear ? Float32(gridconfig.dlon) : 1.0f0
        timestep_ns = Int64(Dates.value(timestep))  # milliseconds as Int64

        # Ensure rng_state array exists and is initialised for new parcels
        _ensure_rng_states!(state, idx)

        # Inbound Dtime → Dtime_ns sync.  The kernel uses Dtime_ns as the parcel
        # clock (advances it by timestep_ns and uses the 1/2 markers for OOB
        # kills), but newly released parcels only have Dtime set — their Dtime_ns
        # is still 0, which the kernel reads as "dead" (≤2) and never advances,
        # collapsing Dtime to the epoch so the parcel is never deposited.  Mirror
        # the live Dtime into Dtime_ns for the active set before launching.
        @inbounds for i in idx
            state.Dtime_ns[i] = Int64(Dates.value(state.Dtime[i]))
        end

        idx_arr = idx  # Vector{Int}

        # Polar-cap switch rows (Inf32/-Inf32 sentinels = cap off)
        swn_j = nrc.polar_n ? nrc.switchnorth_j : Inf32
        sws_j = nrc.polar_s ? nrc.switchsouth_j : -Inf32

        # empty 3-D stubs when Hanna is off: KA needs the argument, the GPU path
        # then uploads 0 bytes and the kernel never reads them
        _empty2d = Array{Float32,3}(undef, 0, 0, 0)
        ustar_arr = need_sfc ? met2.ustar : _empty2d
        pblh_arr  = need_sfc ? met2.pblh  : _empty2d
        ol_arr    = need_sfc ? met2.ol    : _empty2d

        kernel! = advance_parcel_kernel!(ka_backend, block_size)
        scal = (gridconfig, Val(Bool(turb_BL89)), Val(Bool(turb_read_tls)), Val(Bool(use_turb)),
                Val(Bool(use_hanna)), sigw_ft,
                Val(Bool(use_meso)), fturbmeso, r_meso, rs_meso,
                kh_ft_v, Val(kh_ft_v > 0.0f0), Val(Bool(use_htke)), Val(Bool(use_hhanna)), Val(Bool(use_mo)),
                VMIX[],   # per-driver σ_w²/K_z closure — host Ref read, isbits into the kernel
                indt_n, nk, nsub, dt_inc, ifine, t_dir,
                swn_j, sws_j,
                timestep_ns, Ni, Nj, globalx, dlon)

        if ka_backend isa CPU
            kernel!(
                idx_arr,
                state.indx, state.indy, state.indl,
                state.uturb, state.vturb, state.wturb, state.dt,
                state.miss2, state.miss3, state.miss4,
                state.Dtime_ns, state.rng_state,
                uarray, varray, levdotarray,
                zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
                ustar_arr, pblh_arr, ol_arr,
                state.mesou, state.mesov, state.mesow,
                scal...;
                ndrange = length(idx_arr),
            )
            KernelAbstractions.synchronize(ka_backend)
        else
            # GPU: upload met + state into resident device buffers, launch, copy
            # the kernel-written state arrays back.  (Residency/gather is a future
            # optimisation; this is correctness-first.)
            # review 2026-09-04 A5: nest index + met bracket → per-nest resident
            # buffers, met re-uploaded only when the bracket/fields changed.
            _ka_gpu_advance!(kernel!, ka_backend, state, idx_arr,
                             uarray, varray, levdotarray,
                             zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
                             ustar_arr, pblh_arr, ol_arr,
                             scal, n, nrc.t_index)
        end

        # Sync Dtime_ns back to Dtime (Int64 ns → DateTime)
        _sync_dtime_from_ns!(state, idx)
    end
    return nothing
end

# ── GPU launch (_KAGpuBuffers, _KAGPU, _ka_gpu_advance!): moved to
# advance_ka_gpu.jl, included by core/operators.jl only when the GPU stack is
# loaded (EASLEE_ADVANCE=ka_gpu / EASLEE_GPU=1).  The :ka_gpu dispatch above
# only reaches _ka_gpu_advance! when CUDA is functional.

# ── Helper: initialise per-parcel RNG state on first use ─────────────────────
# review 2026-09-04 A7: the stream was seeded from the slot index alone and the
# slot is zeroed on reuse (boundary.jl:88), so successive occupants of a slot
# replayed identical noise and ensemble seeds were ignored.  Seed from
# (slot, first-advance time, run seed): Dtime == currentdate for every parcel in
# idx, so a reused slot (later injection) gets a different stream, and
# EASLEE_SEED (rng.jl RNG_RUN_SEED) shifts the whole ensemble.
function _ensure_rng_states!(state, idx)
    if !hasproperty(state, :rng_state)
        error("ParcelState missing rng_state field — add it to state.jl")
    end
    run_seed = RNG_RUN_SEED[]
    for i in idx
        state.rng_state[i] == UInt64(0) &&
            (state.rng_state[i] = rng_seed(i, Int64(Dates.value(state.Dtime[i])), run_seed))
    end
end

# ── DateTime ↔ Int64 sync ─────────────────────────────────────────────────────
# Dtime_ns stores Dates.value(dt) (milliseconds since epoch as Int64).
# We sync it back after the kernel advances it.
function _sync_dtime_from_ns!(state, idx)
    for i in idx
        ns = state.Dtime_ns[i]
        if ns == Int64(1)
            state.Dtime[i] = DateTime(1)
        elseif ns == Int64(2)
            state.Dtime[i] = DateTime(2)
        else
            state.Dtime[i] = DateTime(Dates.UTInstant(Millisecond(ns)))
        end
    end
end
