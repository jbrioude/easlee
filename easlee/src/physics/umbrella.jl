
# ── Umbrella-intrusion spreading (VOLCANIC_PLUME.md Level 2, transport phase) ──
#
# After detrainment the umbrella spreads as a GRAVITY CURRENT — the neutral-
# buoyancy intrusion driven by its own density anomaly and continuously fed by
# the column's volume flux.  The archived met (≥0.25°) does not contain this
# outflow, so without it the modelled cloud is far too narrow in the first
# hours.  A parcel model has no pressure forces; we impose the SELF-SIMILAR
# radial flow of the classic umbrella law (Woods & Kienle 1994; Costa et al.
# 2013):
#
#     R(t)   = (3 λ N Q_V t² / 2π)^(1/3)        (front radius)
#     Ṙ(t)   = (2/3) R/t                          (front speed)
#     u_r(r) = max( R²/(2·t·max(r, r_src)),       (interior: SOURCE-FED flow —
#                   Ṙ · min(r/R, 1) )              continuity u·2πr·h̄ = Q_V with
#                                                  the mean thickness h̄ =
#                                                  Q_V t/(π R²) ⇒ u = R²/(2 r t);
#                                                  saturated inside the collapse
#                                                  dome r < r_src — parcels there
#                                                  ride the dome outflow; the Ṙ
#                                                  term keeps the head feed and
#                                                  the beyond-front carry)
#
# The 2026-07-07 0.25° verdict run showed why the interior term matters: the
# original u ∝ r profile conserves each parcel's front-fraction f = r/R
# EXACTLY (df/dt = 0), so column material released near the axis could never
# migrate to the head — the parcel footprint tracked only ~½ the analytic
# front and missed the first-hour impulsive expansion entirely (43 km vs the
# observed 199 km at +30 min: the overshoot-dome collapse ejects material at
# the dome radius with tens of m/s).  The source-fed term restores the
# vent→head transport: u(r_src, 10 min) ≈ 40–60 m/s for HT (r_src ≈ 40 km,
# the observed collapse dome) — the observed front-speed class, untuned.
#
# Q_V = MER·amp/ρ(NBH) comes from the SAME 1-D column already validated on the
# heights (amp = mass-flux amplification from the entrainment bookkeeping) and
# N from the met sounding at the NBH — no free parameters beyond the standard
# intrusion constant λ ≈ 0.2.  For Hunga Tonga this predicts R(1 h) ≈ 230 km,
# matching the observed 200–250 km umbrella with nothing tuned.
#
# Handoff.  The PHYSICAL criterion is the intrusion THICKNESS, not a speed:
# closing the current with h = V/(πR²) (V = Q_V·tfeed, constant-volume coasting)
# and u_r = λ N h/2 reproduces the law above exactly, so u_r ∝ h — the current
# spreads for as long as it remains a distinct layer, and stops when it thins to
# the ambient mixing scale.  `hmin` encodes that (h < hmin ⇒ hand off).  Shear
# does NOT stop it: u_r/(|du/dz|·h) = λN/(2|du/dz|) is independent of t and R
# (1.9 for HT with ERA5's 1.36e-3 s⁻¹) — the current always outruns the tearing.
#
# `umin`/`tmax_s` are the older, cruder handoffs kept for compatibility.  Beware:
# umin conflates "slower than the wind" with "stopped spreading" — the wind
# ADVECTS the cloud, it does not halt the radial expansion.  For HT the 5 m/s
# default fires at 8 h (R = 432 km, σ_lat = 1.94°) and freezes the cloud at 55 %
# of the MLS width; the thickness closure run to h = 350 m tracks MLS to ~10 %
# on every day (2026-07-27, see results/HT_EXPERIMENT_LOG.md).  Momentum bookkeeping: the outflow momentum is
# supplied by the intrusion's pressure gradient (continuously fed), NOT by the
# parcels' detrainment w — which is ≈ 0 by construction (see plume_phase.jl).

"""
Per-release umbrella-intrusion parameters (registered at first release of a
`plume_rise.mode="parcels"` box with `umbrella_spread = true`).
"""
struct UmbrellaIntrusion
    lam::Float64          # intrusion constant λ (≈ 0.2, Costa et al. 2013)
    qv_flux::Float64      # volume flux into the intrusion at the NBH [m³/s]
    nbv::Float64          # Brunt–Väisälä frequency at the NBH [1/s]
    t0::DateTime          # eruption (phase) start — intrusion age reference
    tmax_s::Float64       # stop pushing after this age [s]
    umin::Float64         # stop pushing once Ṙ < umin [m/s] (met takes over)
    vlon::Float64         # vent [°E, °N]
    vlat::Float64
    rsrc::Float64         # source/collapse-dome radius [m]: the source-fed
                          # interior flow saturates inside it (u(r<rsrc) =
                          # u(rsrc)) — HT-class ≈ 40 km, default 10 km
    wind_shield::Bool     # shield in-cloud parcels from the RESOLVED wind by
                          # s = min(1, u_r/|u_amb|): the intrusion is a density
                          # current that EXPELS the ambient flow — a parcel
                          # model superimposing met advection double-counts it
                          # (HT: the observed cloud spread as a vent-centred
                          # circle while the outflow ≫ ambient wind; the drift
                          # begins as the current weakens — self-regulating,
                          # no time constant)
    hmin::Float64         # PHYSICAL handoff [m]: stop once the intrusion has
                          # thinned to h = V/(πR²) < hmin (it is no longer a
                          # distinct layer).  0 ⇒ disabled (use umin/tmax only).
    # Altitude-dependent spread amplification.  The single-level law uses N and
    # Q_V at the NBH and pushes every parcel with the same u_r(r,t) regardless
    # of z.  The observed HT cloud is WIDER aloft: the 2026-07-29 backward
    # footprints put the strongest sonde filaments at ±5-6° of latitude at
    # 27-31 km where the modelled cloud ends at ±3.5-4° (miss rate 79 % beyond
    # 4°, 62-66 % above 26 km) — physically, the upper levels are fed by the
    # dome-collapse outflow, more energetic than the NBH-fed intrusion the law
    # describes.  u_r is scaled by 1 + zgain·ramp(z; z0→z1); zgain = 0 ⇒ off,
    # bit-exact legacy.  The target is observed GEOMETRY; mass is untouched.
    zgain::Float64        # amplification at/above zg_z1 (0 = off)
    zg_z0::Float64        # ramp start [m ASL]
    zg_z1::Float64        # ramp end   [m ASL]
    # Azimuth-dependent spread amplification.  The axisymmetric law cannot make
    # the MERIDIONAL wings the HT observations demand: the 2026-07 inversion
    # (x̂ = f·x_b + a at +32 h) put 53 % of its additive water on the northern
    # flank (−12/−16°) and 9 % on the southern (−26°) at 24-30 km — water the
    # disc solution never delivers.  Rather than re-initializing the state at
    # t_ref (which skips the 0-32 h microphysics for the added wings), speed up
    # the radial push in the meridional quadrants: u_r × (1 + g·cos²θ·ramp(z)),
    # θ the azimuth from north, g = azgain_n (ry ≥ 0) or azgain_s (ry < 0).
    # cos²θ = (ry/r)² — full gain on the N-S axis, none on the E-W axis, so the
    # zonal spread (already right vs MLS) is untouched.  Mass is untouched; the
    # wings now RIDE the full microphysics (freeze/fall/sublimate) from t = 0.
    azgain_n::Float64     # meridional amplification, northern half (0 = off)
    azgain_s::Float64     # meridional amplification, southern half (0 = off)
    az_z0::Float64        # ramp start [m ASL]; az_z1 ≤ az_z0 ⇒ no ramp (all z)
    az_z1::Float64        # ramp end   [m ASL]
    tfeed_s::Float64      # FEEDING duration [s] (mass-weighted release-end age,
                          # set at registration).  The constant-flux law only
                          # holds while the column feeds; afterwards the
                          # intrusion COASTS as a constant-volume current
                          # R ∝ t^(1/3) and the source-fed interior flow dies.
                          # Inf → constant feeding (legacy behaviour).
end

UmbrellaIntrusion(; lam=0.2, qv_flux, nbv, t0, tmax_s=12*3600.0, umin=5.0,
                    vlon, vlat, rsrc=10.0e3, wind_shield=false, hmin=0.0,
                    zgain=0.0, zg_z0=26.0e3, zg_z1=30.0e3,
                    azgain_n=0.0, azgain_s=0.0, az_z0=0.0, az_z1=0.0,
                    tfeed_s=Inf) =
    UmbrellaIntrusion(lam, qv_flux, nbv, t0, tmax_s, umin, vlon, vlat, rsrc,
                      wind_shield, hmin, zgain, zg_z0, zg_z1,
                      azgain_n, azgain_s, az_z0, az_z1, tfeed_s)

"""
Intrusion thickness h(t) = V/(πR²) [m], V = Q_V·tfeed the volume delivered at the
neutral level.  Returns Inf while the column still feeds (V is not yet fixed).
"""
umbrella_thickness(u::UmbrellaIntrusion, t::Real) =
    (t <= u.tfeed_s || !isfinite(u.tfeed_s)) ? Inf :
    u.qv_flux * u.tfeed_s / (π * umbrella_radius(u, t)^2)

"""
Front radius R(t) [m]: constant-flux (3λNQ_V t²/2π)^(1/3) while the column
feeds (t ≤ tfeed), then constant-volume coasting R(tfeed)·(t/tfeed)^(1/3) —
the HT stereo track follows exactly this transition (272 km at +1.5 h ×
(4.5/1.5)^(1/3) = 392 km vs the observed 389 km plateau).
"""
umbrella_radius(u::UmbrellaIntrusion, t::Real) =
    t <= 0.0 ? 0.0 :
    t <= u.tfeed_s ? cbrt(3.0 * u.lam * u.nbv * u.qv_flux * t^2 / (2π)) :
    cbrt(3.0 * u.lam * u.nbv * u.qv_flux * u.tfeed_s^2 / (2π)) * cbrt(t / u.tfeed_s)

"Front speed Ṙ(t) [m/s]: (2/3)R/t while fed, R/(3t) coasting."
umbrella_edge_speed(u::UmbrellaIntrusion, t::Real) =
    t <= 0.0 ? 0.0 :
    t <= u.tfeed_s ? (2.0 / 3.0) * umbrella_radius(u, t) / t :
                     umbrella_radius(u, t) / (3.0 * t)

"""
    apply_umbrella_spread!(state, gridconfig, globalx, Ni, umb, boxset,
                           currentdate, dt_s) → n_pushed

Add the intrusion's radial displacement to every active parcel of the release
boxes in `boxset` (the expanded boxes of one volcanic release).  Positions are
grid indices; the index offset from the vent is wrapped across the dateline on
global grids, converted to metres, pushed radially by u_r·dt, and converted
back.  No vertical component (the intrusion spreads at its neutral level).
"""
function apply_umbrella_spread!(state::State.ParcelState, gridconfig, globalx::Bool,
                                Ni::Integer, umb::UmbrellaIntrusion,
                                boxset::AbstractSet{<:Integer},
                                currentdate::DateTime, dt_s::Real;
                                met = nothing, indt::Real = 1.0)
    age = Dates.value(currentdate - umb.t0) / 1000.0
    (age <= 0.0 || age > umb.tmax_s) && return 0
    R  = umbrella_radius(umb, age)
    Rd = umbrella_edge_speed(umb, age)
    (R <= 0.0 || Rd < umb.umin) && return 0
    (umb.hmin > 0.0 && umbrella_thickness(umb, age) < umb.hmin) && return 0
    # vent in grid indices
    vy, vx = latlontoindyindx(gridconfig, umb.vlat, umb.vlon)
    npush = 0
    @inbounds for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        Int(state.relbox[i]) in boxset || continue
        dxi = Float64(state.indx[i]) - vx
        if globalx                              # shortest way around the seam
            per = 360.0 / Float64(gridconfig.dlon)
            dxi > per / 2 && (dxi -= per)
            dxi < -per / 2 && (dxi += per)
        end
        dyi = Float64(state.indy[i]) - vy
        # index offsets → metres (local flat-earth at the parcel latitude)
        la, _ = indyindxtolatlon(gridconfig, Float64(state.indy[i]),
                                 Float64(state.indx[i]))
        mx = 111_320.0 * Float64(gridconfig.dlon) * cosd(la)   # m per x-index
        my = 111_320.0 * Float64(gridconfig.dlat)              # m per y-index
        rx = dxi * mx; ry = dyi * my
        r = sqrt(rx * rx + ry * ry)
        r < 1.0 && continue                    # at the axis: no direction
        # source-fed interior (u = R²/2rt, saturated inside the dome; only
        # WHILE the column feeds) with the front-speed floor near/beyond the
        # head — see the header derivation
        ur = Rd * min(r / R, 1.0)
        if age <= umb.tfeed_s
            ur = max(R * R / (2.0 * age * max(r, umb.rsrc)), ur)
        end
        # altitude-dependent amplification (see header): upper levels spread more
        azg = ry >= 0.0 ? umb.azgain_n : umb.azgain_s
        zp = NaN
        if met !== nothing &&
           (umb.zgain > 0.0 || (azg > 0.0 && umb.az_z1 > umb.az_z0))
            zp = Float64(interpolation4Dsingle(met.zlev, state.indx[i], state.indy[i],
                                               state.indl[i], Float32(indt)))
        end
        if umb.zgain > 0.0 && isfinite(zp) && zp > umb.zg_z0
            ur *= 1.0 + umb.zgain *
                  clamp((zp - umb.zg_z0) / max(umb.zg_z1 - umb.zg_z0, 1.0), 0.0, 1.0)
        end
        # azimuth-dependent amplification (see struct): meridional wings
        if azg > 0.0
            rz = umb.az_z1 <= umb.az_z0 ? 1.0 :
                 isfinite(zp) ? clamp((zp - umb.az_z0) /
                                      max(umb.az_z1 - umb.az_z0, 1.0), 0.0, 1.0) : 0.0
            ur *= 1.0 + azg * rz * (ry * ry) / (r * r)
        end
        f  = ur * dt_s / r                     # fractional radial stretch
        state.indx[i] += Float32(f * rx / mx)
        state.indy[i] += Float32(f * ry / my)
        # wind shielding: cancel the fraction of THIS step's resolved advection
        # that the density current expels (s = min(1, u_r/|u_amb|)); only
        # in-cloud (r ≲ front).  Turbulent displacement is left untouched.
        if umb.wind_shield && met !== nothing && r <= 1.2 * R
            ua = Float64(interpolation4Dsingle(met.Ux, state.indx[i], state.indy[i],
                                               state.indl[i], Float32(indt)))
            va = Float64(interpolation4Dsingle(met.Vy, state.indx[i], state.indy[i],
                                               state.indl[i], Float32(indt)))
            if isfinite(ua) && isfinite(va)
                ss = min(1.0, ur / max(hypot(ua, va), 1.0e-6))
                state.indx[i] -= Float32(ss * ua * dt_s / mx)
                state.indy[i] -= Float32(ss * va * dt_s / my)
            end
        end
        if globalx                              # re-wrap into [1, Ni+1)
            state.indx[i] < 1.0f0 && (state.indx[i] += Float32(360.0 / gridconfig.dlon))
            state.indx[i] >= Float32(Ni) + 1.0f0 && (state.indx[i] -= Float32(360.0 / gridconfig.dlon))
        end
        npush += 1
    end
    return npush
end
