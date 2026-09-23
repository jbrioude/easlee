# output/cov_kernel.jl — flow-evolved per-parcel covariance for the output
# estimator ([output] kernel = "covariance").
#
# Instead of the kNN STATISTICAL bandwidth (spacing among trajectories), each
# parcel carries the PHYSICAL 2×2 horizontal covariance of the air it
# represents, evolved every sync step by
#
#     P ← F·P·Fᵀ + 2·D·dt          F = I + ∇u·dt   (resolved deformation:
#                                                    wind shear stretches and
#                                                    rotates the footprint)
#     D = diag(σ_v²·T_L, σ_v²·T_L)                  (unresolved/turbulent
#                                                    relative-dispersion
#                                                    diffusivity, Taylor limit)
#
# and used as the Gaussian width only WHEN a concentration is evaluated:
# the grid deposit takes the axis-aligned widths (σx=√Pxx, σy=√Pyy — the
# separable renormalised deposit machinery, conservation exact), the gridless
# receptor evaluates the FULL bivariate form exp(−½·mahal)/2π√det.  R(x)
# reports the area-equivalent width (det P)^(1/4).
#
# ∇u is estimated by centred finite differences of the archived wind at the
# parcel (±1 grid index at its own level/time).  The velocity-gradient tensor
# entries have units 1/s directly from Δ(m/s)/Δ(m).
#
# Initialisation: cxx = 0 marks a fresh parcel — its first evolution step sets
# the release footprint σ0² I (config `sigma0`).  Positivity is preserved by
# construction (F P Fᵀ is SPD for det F ≠ 0; 2D dt ≥ 0); a growth cap keeps a
# pathological gradient from exploding a single parcel.

"""
    update_parcel_cov!(state, met, gconfig, indt, nn, dt_s;
                       sigma0=50.0, sigma_v=0.5, tl=300.0, sigma_max=5.0e5)
        → n_updated

One covariance evolution step for every active parcel of nest `nn` (threaded).
`met` is the nest's MetFields3D (Ux/Vy in m/s), `indt` the fractional met time
index.
"""
function update_parcel_cov!(state::State.ParcelState, met, gconfig, indt,
                            nn::Int, dt_s::Float64;
                            sigma0::Float64 = 50.0, sigma_v::Float64 = 0.5,
                            tl::Float64 = 300.0, sigma_max::Float64 = 5.0e5)
    Ni = Int(gconfig.Ni); Nj = Int(gconfig.Nj)
    Ddt  = Float32(2.0 * sigma_v^2 * tl * dt_s)     # 2D·dt per axis [m²]
    s02  = Float32(sigma0^2)
    smax = Float32(sigma_max^2)
    dtf  = Float32(dt_s)
    indtf = Float32(indt)
    nn8 = Int8(nn)
    # metres per MET index (equator value for Rectilinear — the per-parcel lon
    # convergence comes from _xcos; metres directly for Lambert, _xcos = 1)
    mx0, my0 = _out_spacing_m(gconfig)
    n = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in eachindex(state.status)
        (state.status[i] == 0x01 && state.nest[i] == nn8) || continue
        @inbounds begin
            if state.cxx[i] <= 0.0f0                 # fresh parcel / reused slot
                state.cxx[i] = s02; state.cyy[i] = s02; state.cxy[i] = 0.0f0
            end
            x = state.indx[i]; y = state.indy[i]; l = state.indl[i]
            # metre-per-index at the parcel latitude
            mx = Float32(mx0 * _xcos(gconfig, y))
            my = Float32(my0)
            # centred differences of the archived wind (±1 index, clamped)
            xm = max(x - 1.0f0, 1.0f0); xp = min(x + 1.0f0, Float32(Ni) - 1.0f-3)
            ym = max(y - 1.0f0, 1.0f0); yp = min(y + 1.0f0, Float32(Nj) - 1.0f-3)
            u_xp = interpolation4Dsingle(met.Ux, xp, y, l, indtf)
            u_xm = interpolation4Dsingle(met.Ux, xm, y, l, indtf)
            u_yp = interpolation4Dsingle(met.Ux, x, yp, l, indtf)
            u_ym = interpolation4Dsingle(met.Ux, x, ym, l, indtf)
            v_xp = interpolation4Dsingle(met.Vy, xp, y, l, indtf)
            v_xm = interpolation4Dsingle(met.Vy, xm, y, l, indtf)
            v_yp = interpolation4Dsingle(met.Vy, x, yp, l, indtf)
            v_ym = interpolation4Dsingle(met.Vy, x, ym, l, indtf)
            (isfinite(u_xp) && isfinite(u_xm) && isfinite(u_yp) && isfinite(u_ym) &&
             isfinite(v_xp) && isfinite(v_xm) && isfinite(v_yp) && isfinite(v_ym)) || continue
            ux = (u_xp - u_xm) / ((xp - xm) * mx)    # ∂u/∂x [1/s]
            uy = (u_yp - u_ym) / ((yp - ym) * my)
            vx = (v_xp - v_xm) / ((xp - xm) * mx)
            vy = (v_yp - v_ym) / ((yp - ym) * my)
            # F = I + A·dt
            f11 = 1.0f0 + ux * dtf; f12 = uy * dtf
            f21 = vx * dtf;         f22 = 1.0f0 + vy * dtf
            a = state.cxx[i]; b = state.cxy[i]; d = state.cyy[i]
            # F·P·Fᵀ (symmetric result)
            na = f11 * (f11 * a + f12 * b) + f12 * (f11 * b + f12 * d) + Ddt
            nb = f21 * (f11 * a + f12 * b) + f22 * (f11 * b + f12 * d)
            nd = f21 * (f21 * a + f22 * b) + f22 * (f21 * b + f22 * d) + Ddt
            # growth cap (per-axis) against pathological gradients
            na = min(na, smax); nd = min(nd, smax)
            # keep positive-definite after the cap: |cxy| ≤ 0.99·√(cxx·cyy)
            blim = 0.99f0 * sqrt(na * nd)
            state.cxx[i] = na
            state.cyy[i] = nd
            state.cxy[i] = clamp(nb, -blim, blim)
        end
        Threads.atomic_add!(n, 1)
    end
    return n[]
end

"Area-equivalent Gaussian width of the parcel's covariance: (det P)^(1/4) [m]."
@inline function cov_sigma_eff(state, i::Int)
    det = Float64(state.cxx[i]) * Float64(state.cyy[i]) -
          Float64(state.cxy[i])^2
    det > 0.0 ? det^0.25 : max(sqrt(max(Float64(state.cxx[i]), 0.0)), 1.0)
end
