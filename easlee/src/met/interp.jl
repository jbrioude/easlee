
# ── Met-field interpolation kernels ───────────────────────────────────────────
# Linear/trilinear/quadrilinear interpolation of the archived met fields at a
# parcel's fractional grid position (indx, indy, indl, indt), plus the
# staggered-wind wrappers (interpolateU/V/W — FLEXPART Arakawa-C handling), the
# irregular-z level⇄height inversions, and the in-place vertical-profile
# extraction (full column and windowed) used by the turbulence step.
# All kernels are allocation-free and GPU-safe.

# ── Interpolation weights ─────────────────────────────────────────────────────
# 1-D
@inline c0(di) = (1.0f0 - di)
@inline c1(di) = di

# 3-D trilinear
@inline c000(di,dj,dk) = (1.0f0-di)*(1.0f0-dj)*(1.0f0-dk)
@inline c001(di,dj,dk) = (1.0f0-di)*(1.0f0-dj)*dk
@inline c010(di,dj,dk) = (1.0f0-di)*dj*(1.0f0-dk)
@inline c011(di,dj,dk) = (1.0f0-di)*dj*dk
@inline c100(di,dj,dk) = di*(1.0f0-dj)*(1.0f0-dk)
@inline c101(di,dj,dk) = di*(1.0f0-dj)*dk
@inline c110(di,dj,dk) = di*dj*(1.0f0-dk)
@inline c111(di,dj,dk) = di*dj*dk

# 4-D quadrilinear (time dimension last).
# Orientation: `indt` is SLOT-SPACE (indt = 1 ⇒ at slot l's valid time, see
# core/metadapter.jl), so the trailing-0 coefficients (applied to field[…, l])
# carry (1 − dl) and the trailing-1 coefficients (field[…, l+1]) carry dl.
# Until 2026-08-11 the two groups were swapped — the winds were time-MIRRORED
# inside every met bracket (exact at mid-bracket, worst at analysis times)
# while the profile family was correct.  Pinned by
# test/test_time_interp_orientation.jl.
@inline c0000(di,dj,dk,dl) = (1.0f0-di)*(1.0f0-dj)*(1.0f0-dk)*(1.0f0-dl)
@inline c0010(di,dj,dk,dl) = (1.0f0-di)*(1.0f0-dj)*dk*(1.0f0-dl)
@inline c0100(di,dj,dk,dl) = (1.0f0-di)*dj*(1.0f0-dk)*(1.0f0-dl)
@inline c0110(di,dj,dk,dl) = (1.0f0-di)*dj*dk*(1.0f0-dl)
@inline c1000(di,dj,dk,dl) = di*(1.0f0-dj)*(1.0f0-dk)*(1.0f0-dl)
@inline c1010(di,dj,dk,dl) = di*(1.0f0-dj)*dk*(1.0f0-dl)
@inline c1100(di,dj,dk,dl) = di*dj*(1.0f0-dk)*(1.0f0-dl)
@inline c1110(di,dj,dk,dl) = di*dj*dk*(1.0f0-dl)
@inline c0001(di,dj,dk,dl) = (1.0f0-di)*(1.0f0-dj)*(1.0f0-dk)*dl
@inline c0011(di,dj,dk,dl) = (1.0f0-di)*(1.0f0-dj)*dk*dl
@inline c0101(di,dj,dk,dl) = (1.0f0-di)*dj*(1.0f0-dk)*dl
@inline c0111(di,dj,dk,dl) = (1.0f0-di)*dj*dk*dl
@inline c1001(di,dj,dk,dl) = di*(1.0f0-dj)*(1.0f0-dk)*dl
@inline c1011(di,dj,dk,dl) = di*(1.0f0-dj)*dk*dl
@inline c1101(di,dj,dk,dl) = di*dj*(1.0f0-dk)*dl
@inline c1111(di,dj,dk,dl) = di*dj*dk*dl

# ── Primitive kernels ─────────────────────────────────────────────────────────

@inline interpolate1D(field, di, i) =
    @inbounds c0(di)*field[i] + c1(di)*field[i+1]

@inline interpolate3D(field, di, dj, dk, i, j, k) =
    @inbounds (  c000(di,dj,dk)*field[i,  j,  k  ]
               + c001(di,dj,dk)*field[i,  j,  k+1]
               + c010(di,dj,dk)*field[i,  j+1,k  ]
               + c011(di,dj,dk)*field[i,  j+1,k+1]
               + c100(di,dj,dk)*field[i+1,j,  k  ]
               + c101(di,dj,dk)*field[i+1,j,  k+1]
               + c110(di,dj,dk)*field[i+1,j+1,k  ]
               + c111(di,dj,dk)*field[i+1,j+1,k+1])

# 3D interpolation on an irregular z field stored as (x, y, z, t).
# Arguments: l = time index (integer), k = level index (integer),
#            dk = time fraction for interpolation between t=l and t=l+1.
@inline @views interpolate3Dirreg(field, di, dj, dk, i, j, l, k) =
    @inbounds (  c000(di,dj,dk)*field[i,  j,  k,  l  ]
               + c001(di,dj,dk)*field[i,  j,  k,  l+1]
               + c010(di,dj,dk)*field[i,  j+1,k,  l  ]
               + c011(di,dj,dk)*field[i,  j+1,k,  l+1]
               + c100(di,dj,dk)*field[i+1,j,  k,  l  ]
               + c101(di,dj,dk)*field[i+1,j,  k,  l+1]
               + c110(di,dj,dk)*field[i+1,j+1,k,  l  ]
               + c111(di,dj,dk)*field[i+1,j+1,k,  l+1])

@inline interpolate4D(field, di, dj, dk, dl, i, j, k, l) =
    @inbounds (  c0000(di,dj,dk,dl)*field[i,  j,  k,  l  ]
               + c0010(di,dj,dk,dl)*field[i,  j,  k+1,l  ]
               + c0100(di,dj,dk,dl)*field[i,  j+1,k,  l  ]
               + c0110(di,dj,dk,dl)*field[i,  j+1,k+1,l  ]
               + c1000(di,dj,dk,dl)*field[i+1,j,  k,  l  ]
               + c1010(di,dj,dk,dl)*field[i+1,j,  k+1,l  ]
               + c1100(di,dj,dk,dl)*field[i+1,j+1,k,  l  ]
               + c1110(di,dj,dk,dl)*field[i+1,j+1,k+1,l  ]
               + c0001(di,dj,dk,dl)*field[i,  j,  k,  l+1]
               + c0011(di,dj,dk,dl)*field[i,  j,  k+1,l+1]
               + c0101(di,dj,dk,dl)*field[i,  j+1,k,  l+1]
               + c0111(di,dj,dk,dl)*field[i,  j+1,k+1,l+1]
               + c1001(di,dj,dk,dl)*field[i+1,j,  k,  l+1]
               + c1011(di,dj,dk,dl)*field[i+1,j,  k+1,l+1]
               + c1101(di,dj,dk,dl)*field[i+1,j+1,k,  l+1]
               + c1111(di,dj,dk,dl)*field[i+1,j+1,k+1,l+1])

# ── Public single-parcel wrappers ─────────────────────────────────────────────
# These match the legacy API used in timeloop.jl and operators.jl.

@inline function _split_t(t::Real)
    t >= 2f0 && return (1.0f0, 1)
    return (Float32(mod(t, 1)), Base.unsafe_trunc(Int, t))
end

# review 2026-09-04 A1: staggered read at the model top.  For indl in
# [nz-0.5, nz) _split_z_stagger gives kk = nz, so the k+1 stencil point is
# field[…, nz+1, …] — under @inbounds that aliases the next time slot's
# surface level (weight up to 0.5) or runs past the array end (l = 2).
# Shift the base down one level with fraction 1 (the value is then exactly
# the top mass level), as interpolation4Dlevsingle already does.
@inline function _clamp_top_k(dk::Float32, kk::Int, nk::Int)
    kk >= nk && return (1.0f0, nk - 1)
    return (dk, kk)
end

@inline function _split_z_stagger(z::Real)
    if z >= 1.5f0
        # General case: stagger level k maps to mass level k + 0.5.
        # e.g. indl=2 → between mass levels 2 and 3 at weight 0.5.
        return (Float32(mod(z + 0.5f0, 1)), Base.unsafe_trunc(Int, z + 0.5f0))
    else
        # Near-surface: particle is between the surface stagger level (indl=1, z=0)
        # and the midpoint (indl=1.5 = first full mass level).
        # Correct linear map: mass_idx = 2*indl - 1
        #   indl=1.0 → mass_idx=1 (surface extra level)
        #   indl=1.5 → mass_idx=2 (first full level, matches main branch exactly)
        # Old formula (2z+0.5) was wrong: discontinuous at z=1.5 and accessed
        # mass levels 2–3 instead of 1–2 for near-surface particles.
        val = max(1.0f0, 2f0 * Float32(z) - 1.0f0)
        return (Float32(mod(val, 1)), Base.unsafe_trunc(Int, val))
    end
end

function interpolation1Dsingle(field, i)
    dt = Float32(mod(i, 1))
    l  = Base.unsafe_trunc(Int, i)
    return interpolate1D(field, dt, l)
end

function interpolation1Dsingleirreg(field, z)
    Nz = Int32(length(field))
    k1 = Int32(2)
    @inbounds for k in 2:Nz
        z <= field[k] && (k1 = k; break)
    end
    k0 = k1 - Int32(1)
    z0 = Float32(field[k0])
    z1 = Float32(field[k1])
    α  = (z - z0) / max(z1 - z0, 1.0f-6)
    return Float32(k0) + α
end

function interpolation3Dsingle(field, i, j, k)
    dk, kk = _split_z_stagger(k)
    dk, kk = _clamp_top_k(dk, kk, size(field, 3))   # review 2026-09-04 A1
    return interpolate3D(field, Float32(mod(i,1)), Float32(mod(j,1)), dk,
                         Base.unsafe_trunc(Int,i), Base.unsafe_trunc(Int,j), kk)
end

function interpolation4Dsingle(field, x, y, z, t)
    dl, l = _split_t(t)
    dz, kk = _split_z_stagger(z)
    dz, kk = _clamp_top_k(dz, kk, size(field, 3))   # review 2026-09-04 A1
    return interpolate4D(field,
                         Float32(mod(x,1)), Float32(mod(y,1)), dz, dl,
                         Base.unsafe_trunc(Int,x), Base.unsafe_trunc(Int,y), kk, l)
end

function interpolation4Dlevsingle(field, i, j, k, t)
    dl, l = _split_t(t)
    dk = Float32(mod(k, 1))
    kk = Base.unsafe_trunc(Int32, k)
    # Exactly-integral top-level query (k = Nk): the k+1 stencil point has
    # weight 0 but still DEREFERENCES field[…, Nk+1, …] under @inbounds — an
    # out-of-bounds touch that segfaults when the address is unmapped (hit by
    # build_domainfill's Δz at the top level).  Shift the base down one level
    # with fraction 1: the value is identically field[…, Nk, …], in bounds.
    if kk >= Int32(size(field, 3))
        kk = Int32(size(field, 3)) - Int32(1)
        dk = 1.0f0
    end
    return interpolate4D(field,
                         Float32(mod(i,1)), Float32(mod(j,1)), dk, dl,
                         Base.unsafe_trunc(Int32,i), Base.unsafe_trunc(Int32,j),
                         kk, l)
end

# ── Optimized multi-field helpers ─────────────────────────────────────────────

# Decompose parcel position once; call interpolate4D for each field.
# Avoids repeating mod/unsafe_trunc for every field at the same parcel location.
# `nk` (optional) = size(field, 3) of the fields the indices will be applied
# to; when given, the top-level stencil is clamped (review 2026-09-04 A1).
@inline function _parcel_indices(x, y, z, t, nk::Int = typemax(Int))
    dl, l   = _split_t(t)
    dz, kk  = _split_z_stagger(z)
    dz, kk  = _clamp_top_k(dz, kk, nk)
    di = Float32(mod(x, 1))
    dj = Float32(mod(y, 1))
    ii = Base.unsafe_trunc(Int, x)
    jj = Base.unsafe_trunc(Int, y)
    return di, dj, dz, dl, ii, jj, kk, l
end

@inline function interp4(field, di, dj, dk, dl, i, j, k, l)
    interpolate4D(field, di, dj, dk, dl, i, j, k, l)
end

# Extract a vertical profile in-place — no allocation.
# `buf` must have length >= nz; caller owns the buffer.
@inline function get_vertical_profile!(buf::AbstractVector,
                                       field::AbstractArray{<:AbstractFloat,4},
                                       indx, indy, indt, nz::Int)
    i1 = Base.unsafe_trunc(Int32, indx)
    i2 = i1 + Int32(1)
    j1 = Base.unsafe_trunc(Int32, indy)
    j2 = j1 + Int32(1)
    t1 = Base.unsafe_trunc(Int32, indt)
    t2 = min(t1 + Int32(1), Int32(size(field, 4)))
    di = Float32(mod(indx, 1.0f0))
    dj = Float32(mod(indy, 1.0f0))
    dt = Float32(mod(indt, 1.0f0))
    w00 = (1.0f0 - di) * (1.0f0 - dj)
    w01 = (1.0f0 - di) * dj
    w10 = di * (1.0f0 - dj)
    w11 = di * dj
    @inbounds for k in 1:nz
        v1 = w00*field[i1,j1,k,t1] + w01*field[i1,j2,k,t1] +
             w10*field[i2,j1,k,t1] + w11*field[i2,j2,k,t1]
        v2 = w00*field[i1,j1,k,t2] + w01*field[i1,j2,k,t2] +
             w10*field[i2,j1,k,t2] + w11*field[i2,j2,k,t2]
        buf[k] = (1.0f0 - dt) * v1 + dt * v2
    end
    return buf
end

# Nearest-column in-place vertical profile — NO horizontal bilinear blend: the
# column of the nearest grid point, time-interpolated.  FLEXPART samples its
# turbulence coefficient profiles the same way (per grid column); for the
# stochastic Langevin coefficients (theta/tke/rho/drhodz) the horizontal blend
# is sub-noise, and this is 4× fewer array reads than the 4-neighbour version —
# used by the once-per-step hoisted extraction in advance.jl (height keeps the
# full bilinear column: its ladder is searched end-to-end by the vertical index
# inversion, where a seam between neighbouring columns would bias z→indl).
@inline function get_vertical_profile_nearest!(buf::AbstractVector,
                                               field::AbstractArray{<:AbstractFloat,4},
                                               indx, indy, indt, nz::Int)
    i  = clamp(Base.unsafe_trunc(Int32, indx + 0.5f0), Int32(1), Int32(size(field, 1)))
    j  = clamp(Base.unsafe_trunc(Int32, indy + 0.5f0), Int32(1), Int32(size(field, 2)))
    t1 = clamp(Base.unsafe_trunc(Int32, indt), Int32(1), Int32(size(field, 4)))
    t2 = min(t1 + Int32(1), Int32(size(field, 4)))
    w  = Float32(mod(indt, 1.0f0))
    wm = 1.0f0 - w
    @inbounds for k in 1:nz
        buf[k] = wm * field[i, j, k, t1] + w * field[i, j, k, t2]
    end
    return buf
end

# Windowed in-place vertical profile — only levels [kc-W, kc+W] are bilinearly
# interpolated; the rest of the buffer is constant-extended from the window edges
# (a cheap store, not an interpolation) so any out-of-window read stays valid and
# bounded.  The turbulence scheme reads these coefficient profiles (theta/tke/rho/
# drhodz) only at the parcel's own level ±1, so for any W covering the per-substep
# vertical excursion the result is identical to the full extraction — at a fraction
# of the bilinear cost (the dominant term).  `height` must still use the full
# extraction (its column is searched end-to-end by interpolation1Dsingleirreg).
@inline function get_vertical_profile_window!(buf::AbstractVector,
                                              field::AbstractArray{<:AbstractFloat,4},
                                              indx, indy, indt, kc::Int, W::Int, nz::Int)
    i1 = Base.unsafe_trunc(Int32, indx); i2 = i1 + Int32(1)
    j1 = Base.unsafe_trunc(Int32, indy); j2 = j1 + Int32(1)
    t1 = Base.unsafe_trunc(Int32, indt)
    t2 = min(t1 + Int32(1), Int32(size(field, 4)))
    di = Float32(mod(indx, 1.0f0)); dj = Float32(mod(indy, 1.0f0)); dt = Float32(mod(indt, 1.0f0))
    w00 = (1.0f0 - di)*(1.0f0 - dj); w01 = (1.0f0 - di)*dj
    w10 = di*(1.0f0 - dj);           w11 = di*dj
    klo = max(1, kc - W); khi = min(nz, kc + W)
    @inbounds for k in klo:khi
        v1 = w00*field[i1,j1,k,t1] + w01*field[i1,j2,k,t1] +
             w10*field[i2,j1,k,t1] + w11*field[i2,j2,k,t1]
        v2 = w00*field[i1,j1,k,t2] + w01*field[i1,j2,k,t2] +
             w10*field[i2,j1,k,t2] + w11*field[i2,j2,k,t2]
        buf[k] = (1.0f0 - dt) * v1 + dt * v2
    end
    @inbounds begin
        lo = buf[klo]; for k in 1:(klo-1); buf[k] = lo; end
        hi = buf[khi]; for k in (khi+1):nz; buf[k] = hi; end
    end
    return buf
end

# Allocating version kept for compatibility with call sites not yet updated.
function get_vertical_profile(field::Array{Float32,4}, indx, indy, indt, nz::Int32)
    buf = Vector{Float32}(undef, nz)
    return get_vertical_profile!(buf, field, indx, indy, indt, Int(nz))
end

# ── Irregular-z functions ─────────────────────────────────────────────────────

function interpolation4D_irregzsingle(field, indx, indy, z, indt, Nkmin, Nkmax)
    dl, l = _split_t(indt)
    di = Float32(mod(indx, 1))
    dj = Float32(mod(indy, 1))
    ii = Base.unsafe_trunc(Int32, indx)
    jj = Base.unsafe_trunc(Int32, indy)
    z1 = zero(Float32)
    # 2026-09-22: the loop stopped at Nkmax-1, so the TOP layer [z(Nk-1), z(Nk)] was
    # never bracketed: every height inside it fell through to the fallback and was
    # placed on the layer's LOWER boundary.  The top model layer was unreachable on
    # the KA path (well-mixed test: a -20 % void in the top bin, zero parcels above
    # z(Nk) - 50 m, while the CPU's interpolation1Dsingleirreg, which loops 2:Nz,
    # fills it).  Loop to Nkmax; the fallback now only serves z above z(Nk).
    @inbounds for k in Nkmin:Nkmax
        val = interpolate3Dirreg(field, di, dj, dl, ii, jj, l, k)
        if val < z
            z1 = val
        else
            return k - 1 + (z - z1) / max(val - z1, 1.0f-6)
        end
    end
    return Float32(Nkmax) - 1.0f-3
end

# ── Staggered-grid wind interpolation ─────────────────────────────────────────
# U is staggered half a cell in x; V in y.
# W (vertical wind / levdot) MUST use interpolation4Dlevsingle — never
# interpolation4Dsingle — because levdot shares the same k=1…Nk index space as
# zlev.  Applying _split_z_stagger to levdot would shift the level index by 0.5
# and read from the wrong vertical layer.
# Wedge is always set to `true` in initgridconfig for both Lambert and
# Rectilinear grids, so this path (interpolation4Dlevsingle) is always taken.

@inline function interpolateU(edge::Bool, array, x, y, z, t)
    edge ? interpolation4Dsingle(array, x + 0.5f0, y, z, t) :
           interpolation4Dsingle(array, x, y, z, t)
end

@inline function interpolateV(edge::Bool, array, x, y, z, t)
    edge ? interpolation4Dsingle(array, x, y + 0.5f0, z, t) :
           interpolation4Dsingle(array, x, y, z, t)
end

@inline function interpolateW(edge::Bool, array, x, y, z, t)
    # edge=true  → interpolation4Dlevsingle (no _split_z_stagger) — always used
    # edge=false → interpolation4Dsingle    (applies _split_z_stagger) — never used
    edge ? interpolation4Dlevsingle(array, x, y, z, t) :
           interpolation4Dsingle(array, x, y, z, t)
end

