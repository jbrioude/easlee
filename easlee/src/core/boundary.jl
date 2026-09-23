
# ── Boundary-flux parcel injection ───────────────────────────────────────────
#
# Implements the FLEXPART DOMAINFILL strategy: at every timestep inject new
# parcels at domain boundaries proportional to the incoming mass flux.
# Fractional parcel counts are resolved by stochastic rounding so the
# expectation is exact and no mass bias accumulates.
#
# Applicable to:
#   • Well-mixed criterion (3D wind)  — keeps the parcel ensemble stationary
#   • Chemistry                       — boundary inflow of tracers
#   • Any limited-area domain-filling run
#
# Usage:
#   1. Allocate extra state slots at startup via extra_slots in build_runtime.
#   2. Compute mass_per_parcel = compute_domain_mass(nestctx, indt) / NPARTS.
#   3. In the time loop, BEFORE step_runtime!(rt, currentdate):
#        indt_array = MetAdapter.met_time_index_array(currentdate, ctx.nests)
#        inject_boundary_parcels!(state, ctx, currentdate, indt_array, mass_per_parcel)
#
# The injected parcels have Dtime = currentdate so advance_active_parcels!
# picks them up in the same step.

# ── Face-metric helpers ───────────────────────────────────────────────────────
# Physical width [m] of a West/East boundary face cell (its north–south extent).
@inline _face_dy_m(gc::Rectilinear, ::Float32) = 1.0f0 / gc.dyconst   # = dlat·111320 m
@inline _face_dy_m(gc::Lambert,     ::Float32) = gc.dy                 # constant [m]

# Physical width [m] of a South/North boundary face cell (its west–east extent).
@inline function _face_dx_m(gc::Rectilinear, indy::Float32)
    lat = gc.latmin + (indy - 1.0f0) * gc.dlat
    return cosd(lat) / gc.dxconst       # = dlon·111320·cos(lat) m
end
@inline _face_dx_m(gc::Lambert, ::Float32) = gc.dx                     # constant [m]

# Layer thickness Δz [m] between level indices k and k+1.
@inline function _dz_lev(zarr, indx::Float32, indy::Float32, k::Int, indt::Float32)
    z_lo = interpolation4Dlevsingle(zarr, indx, indy, Float32(k),   indt)
    z_hi = interpolation4Dlevsingle(zarr, indx, indy, Float32(k+1), indt)
    return max(z_hi - z_lo, 1.0f0)
end

# ── Parcel-pool helpers ───────────────────────────────────────────────────────
# Linear scan for the first unused slot.  Boundary injection is rare (a few
# hundred parcels per step) so the O(n) cost is negligible in practice.
@inline function _find_unused_boundary(state::State.ParcelState)
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x00 && return i
    end
    return 0   # pool exhausted
end

# Write all fields for a boundary-injected parcel.
# is_new = false: initialize_new_parcels! will skip it (position already set).
# relbox = 0:    sentinel for "injected from boundary, no release box".
@inline function _write_boundary_parcel!(state, slot::Int, nest_id::Int,
                                          indx::Float32, indy::Float32, indl::Float32,
                                          mass::Float32, currentdate::DateTime)
    F = eltype(state.indx)
    I = eltype(state.nest)
    @inbounds begin
        state.status[slot]    = 0x01
        state.nest[slot]      = I(nest_id)
        state.indx[slot]      = F(indx)
        state.indy[slot]      = F(indy)
        state.indl[slot]      = F(indl)
        state.mass[slot]      = F(mass)
        state.Dtime[slot]     = currentdate
        state.Dtimeinit[slot] = currentdate
        state.Dtime_ns[slot]  = Int64(Dates.value(currentdate))
        state.relbox[slot]    = I(0)
        state.wturb[slot]     = F(0)
        state.uturb[slot]     = F(0)
        state.vturb[slot]     = F(0)
        state.ugrad[slot]     = F(0)
        state.cxx[slot]       = F(0)   # covariance kernel: re-init on first touch
        state.cxy[slot]       = F(0)
        state.cyy[slot]       = F(0)
        state.vgrad[slot]     = F(0)
        state.wgrad[slot]     = F(0)
        state.dt[slot]        = F(1)
        state.miss[slot]      = F(indl)
        state.miss2[slot]     = F(0)
        state.miss3[slot]     = F(0)
        state.miss4[slot]     = F(0)
        state.airm[slot]      = F(0)   # ice closure: unassigned on slot reuse
        state.is_new[slot]    = false
        state.rng_state[slot] = UInt64(0)
    end
    return nothing
end

# ── Per-face injection kernels ────────────────────────────────────────────────
# West (indx_face=1, u_sign=+1) and East (indx_face=Ni, u_sign=-1) faces.
# The incoming normal velocity is u * u_sign; positive → air flows in.
function _inject_x_face!(state, gc, met, Nk::Int, Nj::Int, nest_id::Int,
                           indx_face::Float32, u_sign::Float32,
                           currentdate::DateTime, indt::Float32,
                           dt::Float32, mass_per_parcel::Float32, rng)
    n_injected = 0
    uarr   = met.Ux
    rhoarr = met.rho
    zarr   = met.zlev

    @inbounds for j in 1:(Nj - 1)
        jf = Float32(j) + 0.5f0        # column centre in y
        dy = _face_dy_m(gc, jf)

        for k in 1:(Nk - 1)
            kf  = Float32(k) + 0.5f0
            u   = interpolation4Dlevsingle(uarr,   indx_face, jf, kf, indt)
            u_in = u * u_sign
            # NaN <= 0 is false in IEEE 754, so isfinite guard is required —
            # WRF terrain-masked cells can produce NaN winds.
            (isfinite(u_in) && u_in > 0.0f0) || continue

            rho  = interpolation4Dlevsingle(rhoarr, indx_face, jf, kf, indt)
            (isfinite(rho) && rho > 0.0f0) || continue
            dz   = _dz_lev(zarr, indx_face, jf, k, indt)
            flux = rho * u_in * dy * dz * dt   # incoming mass [kg]

            n_real = flux / mass_per_parcel
            isfinite(n_real) || continue
            n_inj  = floor(Int, n_real)
            # Stochastic rounding: unbiased handling of the fractional part.
            rand(rng, Float32) < (n_real - Float32(n_inj)) && (n_inj += 1)
            n_inj <= 0 && continue

            # Place the parcel(s) at a random position on the face cell, nudged
            # slightly inside the domain (0.01 cells) to avoid immediate exit.
            indx_inj = indx_face + u_sign * 0.01f0

            for _ in 1:n_inj
                slot = _find_unused_boundary(state)
                slot == 0 && return n_injected   # state pool exhausted

                indy_inj = jf + (rand(rng, Float32) - 0.5f0)
                indl_inj = Float32(k) + rand(rng, Float32)
                _write_boundary_parcel!(state, slot, nest_id,
                                         indx_inj, indy_inj, indl_inj,
                                         mass_per_parcel, currentdate)
                n_injected += 1
            end
        end
    end
    return n_injected
end

# South (indy_face=1, v_sign=+1) and North (indy_face=Nj, v_sign=-1) faces.
function _inject_y_face!(state, gc, met, Nk::Int, Ni::Int, nest_id::Int,
                           indy_face::Float32, v_sign::Float32,
                           currentdate::DateTime, indt::Float32,
                           dt::Float32, mass_per_parcel::Float32, rng)
    n_injected = 0
    varr   = met.Vy
    rhoarr = met.rho
    zarr   = met.zlev

    @inbounds for i in 1:(Ni - 1)
        if_ = Float32(i) + 0.5f0       # column centre in x
        dx  = _face_dx_m(gc, indy_face)

        for k in 1:(Nk - 1)
            kf   = Float32(k) + 0.5f0
            v    = interpolation4Dlevsingle(varr,   if_, indy_face, kf, indt)
            v_in = v * v_sign
            (isfinite(v_in) && v_in > 0.0f0) || continue

            rho  = interpolation4Dlevsingle(rhoarr, if_, indy_face, kf, indt)
            (isfinite(rho) && rho > 0.0f0) || continue
            dz   = _dz_lev(zarr, if_, indy_face, k, indt)
            flux = rho * v_in * dx * dz * dt

            n_real = flux / mass_per_parcel
            isfinite(n_real) || continue
            n_inj  = floor(Int, n_real)
            rand(rng, Float32) < (n_real - Float32(n_inj)) && (n_inj += 1)
            n_inj <= 0 && continue

            indy_inj = indy_face + v_sign * 0.01f0

            for _ in 1:n_inj
                slot = _find_unused_boundary(state)
                slot == 0 && return n_injected

                indx_inj = if_ + (rand(rng, Float32) - 0.5f0)
                indl_inj = Float32(k) + rand(rng, Float32)
                _write_boundary_parcel!(state, slot, nest_id,
                                         indx_inj, indy_inj, indl_inj,
                                         mass_per_parcel, currentdate)
                n_injected += 1
            end
        end
    end
    return n_injected
end

# ── Public API ────────────────────────────────────────────────────────────────

"""
    compute_domain_mass(nestctx, indt) → Float32

Estimate the total air mass [kg] in the domain by integrating ρ·dV over all
grid cells.  Call once at startup to set:

    mass_per_parcel = compute_domain_mass(nestctx, indt) / NPARTS
"""
function compute_domain_mass(nestctx, indt::Float32)
    gc     = nestctx.gridconfig
    rhoarr = nestctx.metfields3D.rho
    zarr   = nestctx.metfields3D.zlev
    Ni     = Int(gc.Ni)
    Nj     = Int(gc.Nj)
    Nk     = Int(nestctx.levelconfig.Nk)

    total = 0.0f0
    @inbounds for j in 1:(Nj - 1), i in 1:(Ni - 1), k in 1:(Nk - 1)
        jf  = Float32(j) + 0.5f0
        if_ = Float32(i) + 0.5f0
        kf  = Float32(k) + 0.5f0
        rho = interpolation4Dlevsingle(rhoarr, if_, jf, kf, indt)
        (isfinite(rho) && rho > 0.0f0) || continue   # skip terrain-masked NaN cells
        dz  = _dz_lev(zarr, if_, jf, k, indt)
        isfinite(dz) || continue
        dx  = _face_dx_m(gc, jf)
        dy  = _face_dy_m(gc, jf)
        total += rho * dx * dy * dz
    end
    return total
end

"""
    inject_boundary_parcels!(state, ctx, currentdate, indt_array, mass_per_parcel; rng)

Inject new parcels at all 4 lateral domain boundaries proportional to the
incoming mass flux, implementing the FLEXPART DOMAINFILL strategy.

Call BEFORE step_runtime!(rt, currentdate).  Injected parcels have
Dtime = currentdate and are advanced in the same timestep.

Returns the total number of injected parcels across all nests.

Pre-requisites:
  • state was built with enough extra_slots (build_runtime(...; extra_slots=N))
  • mass_per_parcel = compute_domain_mass(nestctx, indt) / NPARTS  [kg/parcel]
"""
function inject_boundary_parcels!(state::State.ParcelState, ctx,
                                    currentdate::DateTime,
                                    indt_array,
                                    mass_per_parcel::Float32;
                                    rng = Random.default_rng())
    timestep = ctx.runconfig.timestep
    dt       = Float32(abs(Dates.value(timestep)) / 1000.0f0)
    n_total  = 0

    oc    = ctx.output
    Niout = Int(oc.Niout)
    Njout = Int(oc.Njout)

    for n in 1:length(ctx.nests)
        nestctx = ctx.nests[n]
        gc      = nestctx.gridconfig
        met     = nestctx.metfields3D
        indt    = indt_array[n]

        Ni = Int(gc.Ni)
        Nj = Int(gc.Nj)
        Nk = Int(nestctx.levelconfig.Nk)

        # ── OUTGRID-aligned face positions ────────────────────────────────────
        # The coordinate mapping is: xout = indx * coeffx + offsetx
        #                             yout = indy * coeffy + offsety
        # so: indx = (xout - offsetx) / coeffx
        #
        # We want injection faces to coincide with the OUTGRID domain boundary
        # (as FLEXPART does via nx_we/ny_sn) rather than the raw MET boundary.
        # This ensures parcels injected from, e.g., the south face actually land
        # in OUTGRID cell 1 (jj=1), not outside the grid (jj=0).
        #
        # Each face position is clamped to stay within the MET domain so that
        # @inbounds interpolation in the inner loops remains safe.
        #
        # West  OUTGRID face at xout = 1.0 → indx = (1 - offsetx)/coeffx
        indx_w = clamp((1.0f0 - gc.offsetx) / gc.coeffx, 1.0f0, Float32(Ni - 1))
        # East  OUTGRID face at xout = Niout + 0.5 (centre of last column, inflow from east)
        indx_e = clamp((Float32(Niout) + 0.5f0 - gc.offsetx) / gc.coeffx,
                        1.0f0, Float32(Ni) - 0.5f0)
        # South OUTGRID face at yout = 1.0 → indy = (1 - offsety)/coeffy
        indy_s = clamp((1.0f0 - gc.offsety) / gc.coeffy, 1.0f0, Float32(Nj - 1))
        # North OUTGRID face at yout = Njout + 0.5 (centre of last row, inflow from north)
        indy_n = clamp((Float32(Njout) + 0.5f0 - gc.offsety) / gc.coeffy,
                        1.0f0, Float32(Nj) - 0.5f0)

        # West  (inflow when U > 0)
        n_total += _inject_x_face!(state, gc, met, Nk, Nj, n,
                                    indx_w, +1.0f0,
                                    currentdate, indt, dt, mass_per_parcel, rng)
        # East  (inflow when U < 0)
        n_total += _inject_x_face!(state, gc, met, Nk, Nj, n,
                                    indx_e, -1.0f0,
                                    currentdate, indt, dt, mass_per_parcel, rng)
        # South (inflow when V > 0)
        n_total += _inject_y_face!(state, gc, met, Nk, Ni, n,
                                    indy_s, +1.0f0,
                                    currentdate, indt, dt, mass_per_parcel, rng)
        # North (inflow when V < 0)
        n_total += _inject_y_face!(state, gc, met, Nk, Ni, n,
                                    indy_n, -1.0f0,
                                    currentdate, indt, dt, mass_per_parcel, rng)
    end
    return n_total
end
