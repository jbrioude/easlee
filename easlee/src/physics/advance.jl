
# ── Parcel advection kernel (CPU) ─────────────────────────────────────────────
# FLEXPART lineage: ~ advance.f90 — mean-wind + Langevin-turbulence sub-stepped
# position update with ground reflection, mesoscale OU velocity fluctuations
# (interpol_mesoscale/turbulence_mesoscale) and polar-stereo caps (ngrid=-1/-2).
# advance! is the entry point called from Operators; updatepos! carries the
# threaded per-parcel loop.
#
# Threading: the parcel loop in updatepos! is parallelised with @threads :static.
# Each thread owns a private set of 6 profile buffers (zero heap allocation in
# the hot loop).  randn(Float32) is thread-safe in Julia ≥ 1.7 (per-thread PRNG).
#
# Type discipline: metfields3D_n is now a MetFields3D struct with concrete array
# fields, so the compiler fully specialises every interpolation call without any
# Dict-lookup overhead or dynamic-dispatch boxing.

# Vertical half-width (in model levels) for the windowed turbulence-profile
# extraction.  The coefficient profiles (theta/tke/rho/drhodz) are only read by
# the Langevin step at the parcel's own level ±1, so extracting the full Nk column
# every sub-step is wasted work; a window of ±W around the parcel covers every
# level it can reach in one sub-step.  PROFW=0 restores the full-column extraction
# (baseline); the default is on.  See get_vertical_profile_window! in met/interp.jl.
const TURB_PROFILE_WINDOW = parse(Int, get(ENV, "PROFW", "24"))

# Persistent per-(nz,nthr) scratch profile buffers for updatepos! — allocated
# once and reused across all model steps (the buffers are per-parcel scratch,
# nothing survives a parcel iteration).  8 buffer-vectors × nthr thread slots.
const _PROFILE_BUF_CACHE = Dict{Tuple{Int,Int}, NTuple{8, Vector{Vector{Float32}}}}()
function _profile_buffers(nz_int::Int, nthr::Int)
    get!(_PROFILE_BUF_CACHE, (nz_int, nthr)) do
        ntuple(_ -> [Vector{Float32}(undef, nz_int) for _ in 1:nthr], 8)
    end
end

function advance!(
    trajsin,
    timestep,
    rc_n::NestRunConfig,
    metfields3D_n::MetFields3D,
    metfields2D_n::MetFields2D,
    gridconfig_n,
    levelconfig_n,
    indt_n,
    idx,
)
    uarray      = metfields3D_n.Ux
    varray      = metfields3D_n.Vy
    levdotarray = metfields3D_n.levdot
    # Pull geometry from gridconfig — already typed, no Dict lookup.
    globalx::Bool = gridconfig_n.globalx
    Ni::Int32     = gridconfig_n.Ni
    Nj::Int32     = gridconfig_n.Nj

    updatepos!(trajsin, idx, indt_n, timestep, rc_n,
               gridconfig_n, levelconfig_n,
               uarray, varray, levdotarray, metfields3D_n, metfields2D_n, globalx, Ni, Nj)
    return nothing
end

# σ of a met field over the 2×2×2 grid points surrounding a parcel, both time
# slots (16 values) — FLEXPART's interpol_mesoscale: the spatial spread of the
# resolved winds near the grid cutoff, proxy for the unresolved sub-grid part.
@inline function _std16(a, i0::Int, j0::Int, k0::Int)
    s = 0.0f0; s2 = 0.0f0
    @inbounds for l in 1:2, k in k0:(k0 + 1), j in j0:(j0 + 1), i in i0:(i0 + 1)
        v = a[i, j, k, l]
        isfinite(v) || return 0.0f0          # masked/terrain cell → no fluctuation
        s += v; s2 += v * v
    end
    m = s / 16.0f0
    var = s2 / 16.0f0 - m * m
    return var > 0.0f0 ? sqrt(var) : 0.0f0
end

@inline function _mesoscale_sigmas(uarr, varr, warr,
                                   xi::Float32, yi::Float32, zi::Float32, nz_int::Int)
    i0 = clamp(Base.unsafe_trunc(Int, xi), 1, size(uarr, 1) - 1)
    j0 = clamp(Base.unsafe_trunc(Int, yi), 1, size(uarr, 2) - 1)
    k0 = clamp(Base.unsafe_trunc(Int, zi), 1, nz_int - 1)
    return (_std16(uarr, i0, j0, k0),
            _std16(varr, i0, j0, k0),
            _std16(warr, i0, j0, k0))
end

# OU update of the per-parcel mesoscale velocity (FLEXPART turbulence_mesoscale):
# v ← r·v + √(1−r²)·N(0,1)·fm·σ, stationary std = fm·σ, lag-1 autocorr = r.
# A fresh parcel (exact-zero triple) gets a full-variance draw — no spin-up bias.
@inline function _mesoscale_update!(trajsin, i::Int, usig::Float32, vsig::Float32,
                                    wsig::Float32, r::Float32, rs::Float32, fm::Float32)
    if trajsin.mesou[i] == 0.0f0 && trajsin.mesov[i] == 0.0f0 &&
       trajsin.mesow[i] == 0.0f0
        trajsin.mesou[i] = randn(Float32) * usig * fm
        trajsin.mesov[i] = randn(Float32) * vsig * fm
        trajsin.mesow[i] = randn(Float32) * wsig * fm
    else
        trajsin.mesou[i] = r * trajsin.mesou[i] + rs * randn(Float32) * usig * fm
        trajsin.mesov[i] = r * trajsin.mesov[i] + rs * randn(Float32) * vsig * fm
        trajsin.mesow[i] = r * trajsin.mesow[i] + rs * randn(Float32) * wsig * fm
    end
    return nothing
end

# review 2026-09-04 A2: number of CTL sub-steps and their length.  The former
# `while inc_time < float_timestep; inc_time += dt_inc` accumulated in Float32
# and ran ctl+1 sub-steps for ~20 % of (Δt, ctl) pairs (e.g. Δt=3 s, ctl=9 → 10).
# ctl is an integer by contract (Float32 in the config struct).
@inline function _nsubsteps(float_timestep::Float32, ctl::Real)
    nsub = max(1, round(Int, ctl))
    return nsub, float_timestep / Float32(nsub)
end

function updatepos!(
    trajsin, idx, indt_n, timestep, rc_n::NestRunConfig,
    gridconfig, levelconfig,
    uarray, varray, levdotarray, metfields3Din::MetFields3D,
    metfields2Din::MetFields2D,
    globalx, Ni, Nj,
)
    nz::Int32               = size(metfields3Din.zlev, 3)
    nz_int                  = Int(nz)
    float_timestep::Float32 = abs(timestep.value / 1000.0f0)
    nsub, dt_inc            = _nsubsteps(float_timestep, rc_n.ctl)   # review 2026-09-04 A2
    ifine::Float32          = 5.0f0
    # review 2026-09-04 A1: one shared top-of-model clamp for the parcel level
    # index (was written out three times).
    indl_top::Float32       = Float32(nz) - 1.0f-3
    t_direction::Float32    = Float32(rc_n.direction)
    # This nest's own sigma_w^2/K_z closure: with AROME nested inside ERA5 the
    # two drivers declare different ones, and applying the inner nest's CBR
    # constants to the outer nest's Hanna-fill field (or vice versa) would
    # misread the archived tke by a factor 2*r_w.
    vm_n = VMixClosure(rc_n.vmix_native, rc_n.vmix_rw, rc_n.vmix_ck)

    # NestRunConfig.turbulence is Union{Thomson,Hanna,NoTurb} — concrete typed, no runtime dispatch.
    _turb = rc_n.turbulence

    # MetFields3D is a concrete struct — all field accesses below are direct
    # memory reads with no dict-lookup or boxing overhead.
    zlev_arr    = metfields3Din.zlev
    theta_arr   = metfields3Din.theta
    tke_arr     = metfields3Din.tke
    rho_arr     = metfields3Din.rho
    drhodz_arr  = metfields3Din.drhodz
    # One set of 8 scratch profile buffers per thread slot.  These are pure
    # per-parcel scratch, so they persist across model steps — cached by
    # (nz, nthr) instead of re-allocated at every updatepos! call (the
    # per-step GenericMemory churn the profiler flagged).  Lookup runs on the
    # calling thread before the @threads region → no race.
    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    (height_bufs, pttprof_bufs, tkeprof_bufs, rhoprof_bufs,
     rhogradprof_bufs, sigwcol_bufs, tlwcol_bufs, ylcol_bufs) =
        _profile_buffers(nz_int, nthr)

    # review 2026-09-04 A4: mesoscale OU memory r = exp(−2Δt/T_met) and its
    # complement are per-run constants — hoisted out of the parcel loop.
    r_meso::Float32  = exp(-2.0f0 * float_timestep / max(rc_n.windint_s, float_timestep))
    rs_meso::Float32 = sqrt(max(1.0f0 - r_meso * r_meso, 0.0f0))

    # idx is already a Vector{Int} (from findall in advance_active_parcels!)
    # so @threads :static can iterate it without an extra collect().
    Threads.@threads :static for i in idx
        tid = Threads.threadid()

        height_buf       = height_bufs[tid]
        pttprof_buf      = pttprof_bufs[tid]
        tkeprof_buf      = tkeprof_bufs[tid]
        rhoprof_buf      = rhoprof_bufs[tid]
        rhogradprof_buf  = rhogradprof_bufs[tid]
        sigwcol_buf      = sigwcol_bufs[tid]
        tlwcol_buf       = tlwcol_bufs[tid]
        ylcol_buf        = ylcol_bufs[tid]

        dxsave::Float32   = 0.0f0
        dysave::Float32   = 0.0f0
        dindl::Float32    = 0.0f0

        ldt::Float32 = trajsin.dt[i]
        up::Float32  = trajsin.uturb[i]
        vp::Float32  = trajsin.vturb[i]
        wp::Float32  = trajsin.wturb[i]

        # Hanna horizontal Langevin: the PBL parameters are sampled ONCE per
        # sync step at the entry position (like the hoisted turbulence
        # profiles); hanna_h then evaluates σ_u/σ_v/T_L at each sub-step's z.
        hh_ust::Float32 = 0.0f0; hh_pbl::Float32 = 0.0f0
        hh_ol::Float32  = 0.0f0; hh_wst::Float32 = 0.0f0
        if rc_n.horiz_hanna
            ix_h = clamp(trajsin.indx[i], 1.0f0, Float32(size(zlev_arr, 1) - 1))
            iy_h = clamp(trajsin.indy[i], 1.0f0, Float32(size(zlev_arr, 2) - 1))
            hh_ust = interpolate2Dsingle(metfields2Din.ustar, ix_h, iy_h, indt_n)
            hh_pbl = interpolate2Dsingle(metfields2Din.pblh,  ix_h, iy_h, indt_n)
            hh_ol  = interpolate2Dsingle(metfields2Din.ol,    ix_h, iy_h, indt_n)
            hh_wst = hanna_wstar(hh_ust, hh_pbl, hh_ol)
        end

        uout::Float32    = 0.0f0
        vout::Float32    = 0.0f0
        wout::Float32    = 0.0f0
        tke_z::Float32   = 0.0f0
        yl::Float32      = 0.0f0
        totdelz::Float32 = 0.0f0

        # Parcels released outside the domain have OOB indx/indy.
        # All @inbounds interpolation calls below use unsafe_trunc and will
        # segfault on index 0 or > Ni/Nj.  Skip them here; Dtime stays at
        # currentdate so reassign_escaped_parcels! will free the slot.
        # On a global grid the met arrays carry a wrap (halo) column at Ni+1
        # (= copy of column 1), so indx in [Ni, Ni+1) is a VALID antimeridian
        # seam position that interpolates across 180°.  Killing it here (>= Ni)
        # punched a parcel-void at 180°; allow up to Ni+1 when globalx.
        indx_hi::Float32 = globalx ? Float32(Ni) + 1.0f0 : Float32(Ni)
        (trajsin.indx[i] < 1.0f0 || trajsin.indx[i] >= indx_hi ||
         trajsin.indy[i] < 1.0f0 || trajsin.indy[i] >= Float32(Nj)) && continue

        # Wind dimensionality: remember the entry indices so frozen dimensions can
        # be restored after the advance (2D/1D runs — see rc_n.wind_mask).
        ix0::Float32 = trajsin.indx[i]; iy0::Float32 = trajsin.indy[i]; il0::Float32 = trajsin.indl[i]

        # Polar-stereographic cap membership (FLEXPART ngrid = -1/-2) — decided
        # once per step from the entry row: wind sampling below stays at the
        # fixed entry position, so the cap and its rotation angle are constant
        # over the sub-steps.  hemi = 0 → ordinary lat/lon metric update.
        hemi::Float32 = rc_n.polar_n && iy0 > rc_n.switchnorth_j ?  1.0f0 :
                        rc_n.polar_s && iy0 < rc_n.switchsouth_j ? -1.0f0 : 0.0f0

        # ── Vertical turbulence profiles: built ONCE per sync step, hoisted out
        # of the CTL sub-step loop below.  Over one model step a parcel advects
        # and mixes only a short vertical/horizontal distance, so the column at
        # the entry position (ix0,iy0) is a good approximation for every sub-step.
        # This removes the per-sub-step 4D→column re-extraction that dominated the
        # turbulence cost.  Full column (not the windowed path) → no window drift
        # as the parcel moves vertically over the step.
        if !(_turb isa NoTurb)
            indx_safe0 = clamp(ix0, 1.0f0, Float32(size(zlev_arr, 1) - 1))
            indy_safe0 = clamp(iy0, 1.0f0, Float32(size(zlev_arr, 2) - 1))
            get_vertical_profile!(height_buf,              zlev_arr,   indx_safe0, indy_safe0, indt_n, nz_int)
            get_vertical_profile_nearest!(pttprof_buf,     theta_arr,  indx_safe0, indy_safe0, indt_n, nz_int)
            get_vertical_profile_nearest!(tkeprof_buf,     tke_arr,    indx_safe0, indy_safe0, indt_n, nz_int)
            get_vertical_profile_nearest!(rhoprof_buf,     rho_arr,    indx_safe0, indy_safe0, indt_n, nz_int)
            get_vertical_profile_nearest!(rhogradprof_buf, drhodz_arr, indx_safe0, indy_safe0, indt_n, nz_int)
            fill!(sigwcol_buf, NaN32)      # invalide le mémo colonne σw/τ_Lw/ℓ

        end

        for _ in 1:nsub   # review 2026-09-04 A2: exact integer sub-step count
            # dt_inc divides float_timestep evenly (ctl is integer), so there is
            # no fractional remainder — each sub-step is always a full dt_inc slice.
            # Using min() here gave dt=0 on the last sub-step, shortening displacement.
            dt::Float32 = dt_inc

            local ux::Float32 = interpolateU(gridconfig.Uedge, uarray,      trajsin.indx[i], trajsin.indy[i], trajsin.indl[i], indt_n)
            local vy::Float32 = interpolateV(gridconfig.Vedge, varray,      trajsin.indx[i], trajsin.indy[i], trajsin.indl[i], indt_n)
            levdot::Float32 = interpolateW(gridconfig.Wedge, levdotarray, trajsin.indx[i], trajsin.indy[i], trajsin.indl[i], indt_n)

            # Monin-Obukhov surface layer (surface_layer.jl): below the lowest
            # model level the archived wind is a slab — put the parcel back on
            # the similarity profile anchored at that level.
            # review 2026-09-04 A3: z_sl < z1_sl ⇔ indl < 2 (z1_sl is zlev at
            # indl = 2 and zlev is monotone in k) — skip the two 16-point
            # interpolations for every parcel above the lowest model level.
            if rc_n.mo_surface && trajsin.indl[i] < 2.0f0
                z_sl  = interpolation4Dlevsingle(zlev_arr, trajsin.indx[i],
                                                 trajsin.indy[i], trajsin.indl[i], indt_n)
                z1_sl = interpolation4Dlevsingle(zlev_arr, trajsin.indx[i],
                                                 trajsin.indy[i], 2.0f0, indt_n)
                if isfinite(z_sl) && z_sl < z1_sl
                    ust_sl = interpolate2Dsingle(metfields2Din.ustar, trajsin.indx[i],
                                                 trajsin.indy[i], indt_n)
                    ol_sl  = interpolate2Dsingle(metfields2Din.ol, trajsin.indx[i],
                                                 trajsin.indy[i], indt_n)
                    fac_sl = mo_wind_factor(z_sl, z1_sl, ust_sl, ol_sl,
                                            sqrt(ux * ux + vy * vy))
                    ux *= fac_sl
                    vy *= fac_sl
                end
            end

            # Met fields can have NaN at masked/terrain cells; treat as zero displacement
            # rather than propagating NaN into dxsave → indx_prov → unsafe_trunc crash.
            dxsave += (isfinite(ux) ? ux : 0.0f0) * dt * t_direction
            dysave += (isfinite(vy) ? vy : 0.0f0) * dt * t_direction
            # Per-sub-step vertical displacement only (NOT accumulated across sub-steps).
            # Unlike horizontal (dxsave/dysave accumulated, indx/indy updated once after loop),
            # indl is updated inside the loop each sub-step (turbulence changes vertical
            # position incrementally).  Using the accumulated dlev here would apply the total
            # displacement again on every sub-step → (ctl+1)/2 times too much vertical wind.
            dlev_step = (isfinite(levdot) ? levdot : 0.0f0) * dt * t_direction

            # Identity for SigmaPressureLevel; met3D arg is ignored so Any is fine here
            lev_here = indltolev(levelconfig, trajsin.indx[i], trajsin.indy[i], trajsin.indl[i], indt_n, metfields3Din)
            dindl    = levtoindl(levelconfig, trajsin.indx[i], trajsin.indy[i], lev_here + dlev_step, indt_n, metfields3Din) - trajsin.indl[i]

            indl_prov::Float32 = trajsin.indl[i] + dindl

            if _turb isa NoTurb
                # Pure advection: the vertical follows the mean (kinematic) wind only —
                # indl_prov already includes the levdot displacement.  Skip the costly
                # 138-level profile extraction + Langevin step (they only feed turbulence).
                new_indl = indl_prov
                if new_indl <= 1.0f0
                    trajsin.indl[i] = abs(new_indl - 1.0f0) + 1.0f0   # reflect at ground
                else
                    # Clamp at the model top: vertical advection (levdot) can push indl past
                    # Nk; without this the next sub-step's interpolateW reads field[…,Nk+1]
                    # out of bounds (segfault).  nz = size(zlev,3) = Nk.
                    trajsin.indl[i] = min(new_indl, indl_top)
                end
                uout = 0.0f0; vout = 0.0f0; wout = wp
            else
                # Turbulence profiles were built once per sync step at the entry
                # position (hoisted above the sub-step loop) — nothing to
                # re-extract here; the buffers are reused across all sub-steps.

                # Option b: real Hanna on-the-fly — surface params interpolated at
                # the provisional parcel position; Thomson path otherwise.
                local z_new::Float32, delz::Float32, indl_turb::Float32
                if _turb isa Hanna && _turb.turbswitch
                    # review 2026-09-04 A4: provisional horizontal position
                    # (mean wind) is consumed only here — computed only here.
                    # In a polar cap the lat/lon metric is unusable (1/cos(lat)
                    # lon blow-up) — use the stereo-plane update there too so the
                    # surface parameters are sampled at a sane location.
                    prov = hemi == 0.0f0 ?
                        (trajsin.indy[i] + dytodindy(gridconfig, dysave, trajsin.indx[i], trajsin.indy[i]),
                         trajsin.indx[i] + dxtodindx(gridconfig, dxsave, trajsin.indx[i], trajsin.indy[i])) :
                        polar_stereo_update(gridconfig, hemi, trajsin.indy[i], trajsin.indx[i], dxsave, dysave)
                    # Clamp provisional position: @inbounds in interpolate2Dsingle will segfault
                    # on index 0 if a parcel near the boundary drifts outside [1, Ni-1] × [1, Nj-1].
                    indy_safe = clamp(Float32(prov[1]), 1.0f0, Float32(size(zlev_arr, 2) - 1))
                    indx_safe = clamp(Float32(prov[2]), 1.0f0, Float32(size(zlev_arr, 1) - 1))
                    ustar_here = interpolate2Dsingle(metfields2Din.ustar, indx_safe, indy_safe, indt_n)
                    pblh_here  = interpolate2Dsingle(metfields2Din.pblh,  indx_safe, indy_safe, indt_n)
                    ol_here    = interpolate2Dsingle(metfields2Din.ol,    indx_safe, indy_safe, indt_n)
                    wstar_here = hanna_wstar(ustar_here, pblh_here, ol_here)
                    z_new, delz, wp, ldt, tke_z, yl =
                        turbulence_hanna_v1(_turb, indl_prov, dt, ldt, ifine, wp,
                                            height_buf, rhoprof_buf, rhogradprof_buf,
                                            ustar_here, pblh_here, wstar_here, ol_here)
                    indl_turb = z_new > 0.0f0 ?
                        ztoindl(levelconfig, trajsin.indx[i], trajsin.indy[i], z_new, indt_n, zlev_arr) :
                        trajsin.indl[i]
                else
                    # turbulence_v3 already maintains indl on the SAME hoisted height
                    # ladder — returning it skips the ztoindl re-search per sub-step.
                    z_new, delz, wp, ldt, tke_z, yl, indl_turb =
                        turbulence_v3(_turb, indl_prov, dt, ldt, ifine, wp,
                                      height_buf, pttprof_buf, tkeprof_buf,
                                      rhoprof_buf, rhogradprof_buf,
                                      sigwcol_buf, tlwcol_buf, ylcol_buf;
                                      vm = vm_n)
                    z_new <= 0.0f0 && (indl_turb = trajsin.indl[i])
                end

                totdelz += delz
                dindl = indl_turb - trajsin.indl[i]

                # Reflect at ground (indl < 1)
                new_indl = trajsin.indl[i] + dindl
                if new_indl <= 1.0f0
                    trajsin.indl[i] = abs(new_indl - 1.0f0) + 1.0f0
                else
                    trajsin.indl[i] += dindl
                end

                # ── Horizontal Langevin step on the SAME 3-D TKE ──────────
                # A native-e driver archives e = ½(σ_u²+σ_v²+σ_w²), so the
                # horizontal variance is already in the profile the vertical
                # step just read — tke_z and yl are its return values, sampled
                # at this sub-step's position.  OU velocity WITH memory:
                #     u' ← u'·r + σ_u·√(1−r²)·N(0,1),   r = exp(−dt/τ_Lh)
                # exact for constant coefficients over any dt, so no stability
                # constraint on the sub-step.  Direction-agnostic like the
                # vertical step (the noise is symmetric) — no t_direction.
                # up/vp persist across sync steps via trajsin.uturb/vturb.
                if rc_n.horiz_tke
                    tlh, sigh = tke_partition_h(vm_n, tke_z, yl)
                    if sigh > 0.0f0
                        rr = exp(-dt / tlh)
                        sq = sqrt(max(1.0f0 - rr * rr, 0.0f0))
                        up = up * rr + sigh * sq * randn(Float32)
                        vp = vp * rr + sigh * sq * randn(Float32)
                        dxsave += up * dt
                        dysave += vp * dt
                    end
                    uout = up; vout = vp
                elseif rc_n.horiz_hanna
                    # ── Hanna horizontal Langevin (FLEXPART in-PBL form) ───────
                    # σ_u/T_Lu on the grid-x component, σ_v/T_Lv on grid-y, no
                    # rotation into the wind frame (FLEXPART advance.f90:392-401).
                    # z_new is the sub-step height AGL the vertical step just
                    # produced; above the entrainment zone hanna_h returns σ = 0
                    # and u′/v′ decay on T_L instead of being cut.
                    sigu, tlu, sigv, tlv = hanna_h(max(z_new, 0.0f0), hh_pbl, hh_ol, hh_ust, hh_wst)
                    ru = exp(-dt / tlu); rv = exp(-dt / tlv)
                    up = up * ru + sigu * sqrt(max(1.0f0 - ru * ru, 0.0f0)) * randn(Float32)
                    vp = vp * rv + sigv * sqrt(max(1.0f0 - rv * rv, 0.0f0)) * randn(Float32)
                    dxsave += up * dt
                    dysave += vp * dt
                    uout = up; vout = vp
                else
                    uout = 0.0f0; vout = 0.0f0
                end
                wout = wp
            end
        end  # sub-step loop

        # ── Mesoscale velocity fluctuations (FLEXPART usig/vsig/wsig) ─────────
        # Per-parcel OU velocity, updated ONCE per sync step: memory
        # r = exp(−2Δt/T_met) (T_met = met archive interval — the timescale the
        # archived winds cannot resolve), amplitude fturbmeso·σ over the
        # 2×2×2×2 grid-point neighbourhood (sub-grid wind uncertainty).  This
        # is distinct from the PBL turbulence above: it acts in ALL layers.
        # Vertical is in level-index space (the FLEXPART η branch: our levdot
        # IS an index-space velocity) with ground reflection / top clamp.
        if rc_n.mesoscale
            usig, vsig, wsig = _mesoscale_sigmas(uarray, varray, levdotarray,
                                                 trajsin.indx[i], trajsin.indy[i],
                                                 trajsin.indl[i], nz_int)
            _mesoscale_update!(trajsin, i, usig, vsig, wsig, r_meso, rs_meso, rc_n.fturbmeso)
            dxsave += trajsin.mesou[i] * float_timestep * t_direction
            dysave += trajsin.mesov[i] * float_timestep * t_direction
            ml_new = trajsin.indl[i] + trajsin.mesow[i] * float_timestep * t_direction
            ml_new <= 1.0f0 && (ml_new = abs(ml_new - 1.0f0) + 1.0f0)  # ground reflection
            trajsin.indl[i] = min(ml_new, indl_top)
        end

        # ── FLEXPART free-atmosphere horizontal diffusivity (d_trop) ──────────
        # Above the PBL the Hanna/Thomson schemes supply only σ_w (with
        # σ_w = 0.1 m/s, τ_Lw = 10 s ⇒ K_v = 0.1 m²/s = FLEXPART d_strat);
        # their horizontal turbulent component is ZERO by construction
        # (uout/vout above).  FLEXPART instead applies a constant horizontal
        # diffusivity as a memoryless random walk (Stohl et al. 2005,
        # d_trop = 50 m²/s): δ = √(2·K_h·Δt)·N(0,1), once per sync step —
        # a random walk has no velocity memory, so nothing is stored in
        # uturb/vturb.  Applied above the PBL top only; the BL schemes own
        # turbulence below.  TRAJECTORY route: build_run refuses the config
        # if [physics.mixing] is also armed (one K either moves parcels or
        # exchanges mass, never both — double counting).
        if rc_n.kh_ft > 0.0f0
            ix_s = clamp(trajsin.indx[i], 1.0f0, Float32(size(zlev_arr, 1) - 1))
            iy_s = clamp(trajsin.indy[i], 1.0f0, Float32(size(zlev_arr, 2) - 1))
            il_s = clamp(trajsin.indl[i], 1.0f0, indl_top)
            zpar = interpolation4Dlevsingle(zlev_arr, ix_s, iy_s, il_s, indt_n)
            pblh_here2 = interpolate2Dsingle(metfields2Din.pblh, ix_s, iy_s, indt_n)
            if isfinite(zpar) && zpar > max(pblh_here2, 50.0f0)
                sk = sqrt(2.0f0 * rc_n.kh_ft * float_timestep)
                dxsave += sk * randn(Float32)
                dysave += sk * randn(Float32)
            end
        end

        if hemi == 0.0f0
            trajsin.indx[i] += dxtodindx(gridconfig, dxsave, trajsin.indx[i], trajsin.indy[i])
            trajsin.indy[i] += dytodindy(gridconfig, dysave, trajsin.indx[i], trajsin.indy[i])
        else
            # Stereo-plane update — exact across the pole; |lat| ≤ 90 by
            # construction so indy lands in [1, Nj].  Nudge off the exact pole
            # row so the interpolation stencil (rows j, j+1) stays in bounds.
            iy_new, ix_new = polar_stereo_update(gridconfig, hemi,
                                                 trajsin.indy[i], trajsin.indx[i],
                                                 dxsave, dysave)
            trajsin.indy[i] = clamp(iy_new, 1.0f0, Float32(Nj) - 1.0f-3)
            trajsin.indx[i] = ix_new
        end

        trajsin.uturb[i] = uout
        trajsin.vturb[i] = vout
        trajsin.wturb[i] = wout
        trajsin.dt[i]    = ldt
        trajsin.miss2[i] = tke_z
        trajsin.miss3[i] = yl
        trajsin.miss4[i] = totdelz

        # Wind dimensionality: freeze suppressed indices (2D/1D). Restore BEFORE the
        # boundary checks so a frozen dimension can't trigger a spurious OOB kill.
        @inbounds begin
            mx, my, ml = rc_n.wind_mask
            mx || (trajsin.indx[i] = ix0)
            my || (trajsin.indy[i] = iy0)
            ml || (trajsin.indl[i] = il0)
        end

        # A parcel still on the lat/lon path (equatorward of the switch row) can
        # overshoot the pole row in one step when the switch sits close to the
        # pole; mirror it across the pole instead of killing it (FLEXPART's
        # HSO/AL pole fix): lat → 180 − lat, lon → lon + 180.  Cap parcels never
        # trigger this (clamped inside [1, Nj) above).
        if rc_n.polar_n && trajsin.indy[i] > Float32(Nj)
            trajsin.indy[i] = 2.0f0 * Float32(Nj) - trajsin.indy[i]
            trajsin.indx[i] += 180.0f0 / gridconfig.dlon
        elseif rc_n.polar_s && trajsin.indy[i] < 1.0f0
            trajsin.indy[i] = 2.0f0 - trajsin.indy[i]
            trajsin.indx[i] += 180.0f0 / gridconfig.dlon
        end

        # Boundary conditions
        if globalx
            # Periodic wrap into [1, Ni+1).  Upper threshold MUST be Ni+1, not Ni:
            # the halo column at Ni+1 makes [Ni, Ni+1) a valid seam, and wrapping
            # at >= Ni instead dumped those parcels into [0,1) → killed next step
            # (the 180° void).  The lower wrap (+Ni) lands in [1, Ni+1) and is not
            # re-wrapped because indx < Ni+1 there.
            trajsin.indx[i] < 1.0f0                && (trajsin.indx[i] += 360.0f0 / gridconfig.dlon)
            trajsin.indx[i] >= Float32(Ni) + 1.0f0 && (trajsin.indx[i] -= 360.0f0 / gridconfig.dlon)
        else
            (trajsin.indx[i] < 1 || trajsin.indx[i] >= Ni) && (trajsin.Dtime[i] = DateTime(1))
        end
        (trajsin.indy[i] < 1 || trajsin.indy[i] >= Nj) && (trajsin.Dtime[i] = DateTime(2))

        trajsin.Dtime[i] += timestep
    end
    return nothing
end
