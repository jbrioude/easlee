# ── Monin-Obukhov surface-layer wind ([physics].surface_layer_wind) ──────────
#
# The archived wind is constant-extended to the ground: the reader copies the
# lowest full level into the surface slot (layout contract 2026-10-22: U, V
# slot 1 = slot 2), so the interpolated profile below the lowest full level is a
# SLAB.  That level is read at indl = 1.5 (MO_ANCHOR_INDL), height
# zf(2) = ½(zlev[1] + zlev[2]) = ½·zlev[2] — on ERA5 L137 about 10 m, half the
# ~19.6 m of zlev[2] that the ETEX atlas (docs/INTERP_ERROR_ATLAS.md) and the
# pre-contract code took as the lowest wind level.  For indl in [1.5, 2) the
# wind already blends the two lowest full levels and is left alone.
#
# Similarity theory instead gives, in the constant-flux layer,
#
#     u(z) = (u*/κ) [ ln(z/z0) − ψ_m(z/L) ]
#
# Anchoring on the model's own lowest-level wind U1 at z1 removes z0 entirely:
#
#     u(z) = U1 − (u*/κ) [ ln(z1/z) − ψ_m(z1/L) + ψ_m(z/L) ]
#
# so the profile passes exactly through the archived wind at z1 and needs only
# u* and L — both already carried in MetFields2D.  This also keeps the wind
# consistent with whatever u* the run has (archived stress, or the reader's own
# neutral-log-law inversion), instead of introducing a second roughness length.
#
# Direction is preserved: the surface layer's Ekman turning is small compared
# with the speed shear below the lowest full level, and the archived wind
# carries no information to resolve it.
#
# ψ_m: Businger-Dyer / Paulson (1970), as in FLEXPART.

const _MO_KARMAN = 0.4f0

# Level index of the anchor (the lowest full level, class F slot 2): the rescaling
# applies for indl < MO_ANCHOR_INDL, with z1 = zlev interpolated there.
const MO_ANCHOR_INDL = 1.5f0

@inline function _psi_m(zeta::Float32)
    if zeta < 0.0f0                       # unstable
        x = (1.0f0 - 16.0f0 * zeta)^0.25f0
        return 2.0f0 * log(0.5f0 * (1.0f0 + x)) + log(0.5f0 * (1.0f0 + x * x)) -
               2.0f0 * atan(x) + 1.5707963f0
    elseif zeta > 0.0f0                   # stable
        return -5.0f0 * min(zeta, 1.0f0)  # cap: the linear form breaks down for z/L > 1
    else
        return 0.0f0
    end
end

"""
    mo_wind_factor(z, z1, ustar, ol, U1) → Float32

Ratio u(z)/U1 on the Monin-Obukhov profile anchored at (z1, U1).  Returns 1
above z1 (no correction) and for degenerate inputs — the correction may only
ever slow a parcel toward the ground, never invent wind.
"""
@inline function mo_wind_factor(z::Float32, z1::Float32, ustar::Float32,
                                ol::Float32, U1::Float32)::Float32
    (isfinite(z) && isfinite(z1) && isfinite(ustar) && isfinite(U1)) || return 1.0f0
    (z1 > 0.5f0 && U1 > 0.01f0 && ustar > 0.0f0) || return 1.0f0
    z >= z1 && return 1.0f0
    # Floor: below the roughness length the profile is meaningless.  0.5 m is
    # well inside the deposition layer and above any sensible z0 here.
    zc = max(z, 0.5f0)
    invL = (isfinite(ol) && abs(ol) > 1.0f0) ? 1.0f0 / ol : 0.0f0
    du = (ustar / _MO_KARMAN) *
         (log(z1 / zc) - _psi_m(z1 * invL) + _psi_m(zc * invL))
    return clamp((U1 - du) / U1, 0.0f0, 1.0f0)
end
