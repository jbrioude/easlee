# physics/seaspray.jl
# ─────────────────────────────────────────────────────────────────────────────
# OSSA sea-spray aerosol source function — Ovadnevaite et al. (2014), ACP 14,
# 1837–1852: five lognormal modes driven by the wave Reynolds number
#
#     Re_Hw = u* · H_s / ν_w        (their Eq. 2)
#
# instead of wind speed alone — encapsulating wave state (rising/waning wind,
# fetch), friction velocity, and the SST/salinity dependence through the
# kinematic viscosity of seawater.  Mode parameters verbatim from their
# Table 1 (dry count-median diameters; number fluxes in m⁻² s⁻¹):
#
#   i   σ_i    CMD_i [µm]   F_i(Re_Hw)
#   1   1.37   0.018        104.5 ·(Re_Hw − 1×10⁵)^0.556
#   2   1.5    0.041        0.0442·(Re_Hw − 1×10⁵)^1.08
#   3   1.42   0.09         149.6 ·(Re_Hw − 1×10⁵)^0.545
#   4   1.53   0.23         2.96  ·(Re_Hw − 1×10⁵)^0.79
#   5   1.85   0.83         0.51  ·(Re_Hw − 2×10⁵)^0.87
#
# Thresholds: Re_Hw = 1×10⁵ for the submicron modes (≈ onset of wave breaking,
# U10 ≈ 3.7 m s⁻¹ at H_s = 1.23 m, SST = 9 °C), 2×10⁵ for the supermicron
# mode 5.  Validity: derived from 15 nm–6 µm dry diameter, U10 up to ~18 m/s.
#
# MEAN vs TRANSIENT decomposition (this model's extension): the flux is
# strongly nonlinear in the wind, so sub-grid/gust variability contributes
# beyond the mean-field flux — most where the mean wind is low (< 10–15 m/s)
# and the threshold truncates the PDF's lower tail.  `ossa_flux_gusty`
# integrates the mode fluxes over a Gaussian wind distribution
# N(U10_mean, σ_gust) by Gauss–Hermite quadrature and reports BOTH the
# mean-field flux and the gust-integrated flux, so the transient enhancement
# ⟨F⟩/F(Ū) is a first-class diagnostic, not a hidden correction.
#
# Auxiliary closures (each replaceable by archived met fields when available):
#   ν_w(SST, S)  — Sharqawy, Lienhard & Zubair (2010) seawater correlation
#                  (reproduces the paper's reference ν = 1.34e-6 at 9 °C, 35 ‰);
#   H_s          — archived significant wave height, or fully-developed
#                  Pierson–Moskowitz H_s = 0.0246·U10² as fallback;
#   C_d          — archived, constant, or the bulk C_d = (0.8 + 0.065·U10)e-3
#                  (Wu 1980-class); u* = √C_d · U10.
# ─────────────────────────────────────────────────────────────────────────────

# OSSA Table 1 (Ovadnevaite et al. 2014): σ_g, dry CMD [m], amplitude,
# exponent, Re threshold
const OSSA_MODES = (
    (sigma = 1.37, cmd = 0.018e-6, a = 104.5,  b = 0.556, re0 = 1.0e5),
    (sigma = 1.50, cmd = 0.041e-6, a = 0.0442, b = 1.08,  re0 = 1.0e5),
    (sigma = 1.42, cmd = 0.09e-6,  a = 149.6,  b = 0.545, re0 = 1.0e5),
    (sigma = 1.53, cmd = 0.23e-6,  a = 2.96,   b = 0.79,  re0 = 1.0e5),
    (sigma = 1.85, cmd = 0.83e-6,  a = 0.51,   b = 0.87,  re0 = 2.0e5),
)
const OSSA_NMODES = length(OSSA_MODES)

"""
    seawater_kinematic_viscosity(sst_K; salinity_permil=35.0) → ν [m²/s]

Kinematic viscosity of seawater, Sharqawy, Lienhard & Zubair (2010,
Desalin. Water Treat. 16) dynamic-viscosity correlation over ρ_sw.
At 9 °C / 35 ‰ this gives 1.34e-6 m²/s — the OSSA reference value.
"""
function seawater_kinematic_viscosity(sst_K::Real; salinity_permil::Real = 35.0)
    t = clamp(Float64(sst_K) - 273.15, -2.0, 40.0)     # °C
    S = salinity_permil / 1000.0                        # kg/kg
    μw = 4.2844e-5 + 1.0 / (0.157 * (t + 64.993)^2 - 91.296)
    A  = 1.541 + 1.998e-2 * t - 9.52e-5 * t^2
    B  = 7.974 - 7.561e-2 * t + 4.724e-4 * t^2
    μ  = μw * (1.0 + A * S + B * S^2)
    ρ  = 1000.0 * (1.0 + 0.8 * S) * (1.0 - 0.00025 * (t - 10.0))  # light ρ_sw closure
    return μ / ρ
end

"""
    ossa_reynolds(u10, sst_K; hs=nothing, cd=nothing, salinity_permil=35) → Re_Hw

Wave Reynolds number Re_Hw = u*·H_s/ν_w (OSSA Eq. 2).  `hs`/`cd` use the
archived wave fields when given; otherwise H_s falls back to the
fully-developed Pierson–Moskowitz sea H_s = 0.0246·U10² and C_d to the bulk
(0.8 + 0.065·U10)×10⁻³.
"""
function ossa_reynolds(u10::Real, sst_K::Real;
                       hs = nothing, cd = nothing, salinity_permil::Real = 35.0)
    u  = max(Float64(u10), 0.0)
    Cd = cd === nothing ? (0.8 + 0.065 * u) * 1e-3 : Float64(cd)
    Hs = hs === nothing ? 0.0246 * u^2 : Float64(hs)
    ν  = seawater_kinematic_viscosity(sst_K; salinity_permil)
    return sqrt(Cd) * u * Hs / ν
end

"""
    ossa_mode_fluxes(re_hw) → NTuple{5,Float64}  [m⁻² s⁻¹]

Per-mode number fluxes F_i(Re_Hw), zero below each mode's threshold.
"""
@inline function ossa_mode_fluxes(re_hw::Real)
    ntuple(OSSA_NMODES) do i
        m = OSSA_MODES[i]
        x = Float64(re_hw) - m.re0
        x > 0.0 ? m.a * x^m.b : 0.0
    end
end

"""
    ossa_number_flux(re_hw) → Σ F_i  [m⁻² s⁻¹]
    ossa_mass_flux(re_hw; rho_dry=2017.0) → Σ M_i  [kg m⁻² s⁻¹]

Total number flux, and total dry-mass flux using the lognormal third moment
about the count-median diameter: M_i = F_i · (π/6)·ρ·CMD_i³·exp(4.5 ln²σ_i)
(dry sea-salt density 2017 kg/m³ by default).
"""
ossa_number_flux(re_hw::Real) = sum(ossa_mode_fluxes(re_hw))

function ossa_mass_flux(re_hw::Real; rho_dry::Real = 2017.0)
    F = ossa_mode_fluxes(re_hw)
    M = 0.0
    for i in 1:OSSA_NMODES
        m = OSSA_MODES[i]
        M += F[i] * (π / 6) * rho_dry * m.cmd^3 * exp(4.5 * log(m.sigma)^2)
    end
    return M
end

"""
    ossa_dfdlogd(D_m, re_hw) → dF/dlog₁₀D  [m⁻² s⁻¹]

The five-modal size distribution at dry diameter `D_m` [m].
"""
function ossa_dfdlogd(D_m::Real, re_hw::Real)
    F = ossa_mode_fluxes(re_hw)
    s = 0.0
    for i in 1:OSSA_NMODES
        m  = OSSA_MODES[i]
        lσ = log10(m.sigma)
        s += F[i] / (sqrt(2π) * lσ) *
             exp(-(log10(Float64(D_m) / m.cmd))^2 / (2 * lσ^2))
    end
    return s
end

# ── Mean-field vs transient decomposition ────────────────────────────────────

# 8-point Gauss–Hermite nodes/weights for ∫ f(x)·N(0,1) dx
const _GH8_X = (-2.9306374202572440, -1.9816567566958429, -1.1571937124467802,
                -0.3811869902073221,  0.3811869902073221,  1.1571937124467802,
                 1.9816567566958429,  2.9306374202572440)
const _GH8_W = (1.9960407221136762e-4, 1.7077983007413475e-2,
                2.0780232581489188e-1, 6.6114701255824129e-1,
                6.6114701255824129e-1, 2.0780232581489188e-1,
                1.7077983007413475e-2, 1.9960407221136762e-4)
const _GH8_NORM = sqrt(π)

"""
    ossa_flux_gusty(u10_mean, sigma_u, sst_K; hs=nothing, cd=nothing,
                    salinity_permil=35) → (mean = F̄, gusty = ⟨F⟩, modes = ...)

The OSSA mode fluxes evaluated (a) at the MEAN wind (the resolved / mean-field
contribution: wind + wave state) and (b) integrated over a Gaussian wind
distribution N(u10_mean, sigma_u) by 8-point Gauss–Hermite quadrature (the
TRANSIENT contribution: sub-grid variability, gustiness).  The nonlinearity
F ∝ (Re − Re₀)^b with Re ∝ u²·√Cd(u) makes ⟨F⟩ > F(ū), most strongly at low
mean winds where part of the distribution sits below the production
threshold while the mean is above it — or vice versa.

When `hs` is supplied (archived wave field) it is held at its MEAN-state
value inside the quadrature: waves respond to the wind history, not to the
gusts — exactly the mean/transient separation this decomposition is for.
With the Pierson–Moskowitz fallback H_s(u) responds fully (worst case).

Returns per-mode and total number fluxes for both evaluations, so the
enhancement ⟨F⟩/F̄ is available as a diagnostic per mode.
"""
function ossa_flux_gusty(u10_mean::Real, sigma_u::Real, sst_K::Real;
                         hs = nothing, cd = nothing, salinity_permil::Real = 35.0)
    ū  = max(Float64(u10_mean), 0.0)
    σu = max(Float64(sigma_u), 0.0)
    Fm = ossa_mode_fluxes(ossa_reynolds(ū, sst_K; hs, cd, salinity_permil))
    if σu == 0.0
        return (mean = sum(Fm), gusty = sum(Fm), modes_mean = Fm, modes_gusty = Fm)
    end
    hs_fixed = hs === nothing ? nothing : hs   # wave state frozen at the mean
    acc = ntuple(_ -> 0.0, OSSA_NMODES)
    for k in eachindex(_GH8_X)
        u  = max(ū + sqrt(2.0) * σu * _GH8_X[k], 0.0)
        Fk = ossa_mode_fluxes(ossa_reynolds(u, sst_K;
                                            hs = hs_fixed, cd, salinity_permil))
        w  = _GH8_W[k] / _GH8_NORM
        acc = ntuple(i -> acc[i] + w * Fk[i], OSSA_NMODES)
    end
    return (mean = sum(Fm), gusty = sum(acc), modes_mean = Fm, modes_gusty = acc)
end

# ─────────────────────────────────────────────────────────────────────────────
# Engine source term ([[release]] kind = "seaspray")
#
# The release box reserves a relbox for species attribution (0 scheduled
# particles/mass); build_seaspray precomputes the ocean cells (met grid ∩ box,
# lsm < 0.5) once, and apply_seaspray! runs each transport step: OSSA mass
# flux per cell — mean-field AND gust-integrated when gust_sigma > 0 — is
# accumulated, and a parcel of one mass quantum is spawned at the surface of
# a cell whenever its accumulator fills.  The mean-field and gust-integrated
# emissions are tracked SEPARATELY (the transient-enhancement diagnostic).
# ─────────────────────────────────────────────────────────────────────────────

struct SeaSprayBox
    relbox::Int32
    lon0::Float64; lon1::Float64             # release box (ocean-cell scan bounds)
    lat0::Float64; lat1::Float64
    initialized::Base.RefValue{Bool}         # cells scanned on the FIRST apply —
                                             # at build time the met (incl. lsm)
                                             # has not been read yet
    ii::Vector{Int32}; jj::Vector{Int32}     # ocean met cells within the box
    area::Vector{Float64}                    # cell areas [m²]
    acc::Vector{Float64}                     # unspawned accumulated mass [kg]
    quantum::Float64                         # mass per spawned parcel [kg]
    t1::DateTime; t2::DateTime
    sst::Float64                             # K; < 0 → met 2t per cell
    salinity::Float64                        # ‰
    hs::Float64                              # m; < 0 → Pierson–Moskowitz
    cd::Float64                              # < 0 → bulk (0.8 + 0.065·U10)e-3
    gust_sigma::Float64                      # m/s; 0 = mean-field only
    rho_dry::Float64                         # dry sea-salt density [kg/m³]
    emitted::Base.RefValue{Float64}          # actually spawned mass [kg]
    emitted_mean::Base.RefValue{Float64}     # mean-field-only counterfactual [kg]
    emitted_gusty::Base.RefValue{Float64}    # gust-integrated production [kg]
    nspawn::Base.RefValue{Int}
end

# dry-mass flux [kg m⁻² s⁻¹] from per-mode number fluxes (lognormal 3rd moment)
@inline function _ossa_massflux(F::NTuple{5,Float64}, rho_dry::Float64)
    M = 0.0
    for i in 1:OSSA_NMODES
        m = OSSA_MODES[i]
        M += F[i] * (π / 6) * rho_dry * m.cmd^3 * exp(4.5 * log(m.sigma)^2)
    end
    return M
end

"""
    build_seaspray(ctx, m2raw, relbox, r; nest=1) → SeaSprayBox or nothing

`m2raw` is the RAW archived 2-D field dict (rt.metfields2D[nest]) — the typed
nestctx.metfields2D struct only carries the PBL fields.

Precompute the ocean-cell list for a `kind="seaspray"` release `r`
(ReleaseBox from the RunSpec).  Returns nothing (with a warning) when the met
has no land-sea mask or the box contains no ocean cells.
"""
function build_seaspray(ctx, m2::AbstractDict, relbox::Integer, r; nest::Int = 1)
    haskey(m2, "lsm") ||
        (@warn "seaspray: met has no land-sea mask (lsm) — the archive needs surface fields (e.g. an ERA5 namelist with _SFC_ companions); source disabled"; return nothing)
    q = r.mass_per_parcel_kg
    q > 0 || error("seaspray release \"$(r.species)\": mass_per_parcel must be set (the spawn quantum)")
    return SeaSprayBox(Int32(relbox),
                       minimum(r.lon), maximum(r.lon), minimum(r.lat), maximum(r.lat),
                       Ref(false), Int32[], Int32[], Float64[], Float64[], q,
                       r.time[1], r.time[2],
                       r.spray_sst, r.spray_salinity, r.spray_hs, r.spray_cd,
                       r.spray_gust_sigma, r.spray_rho,
                       Ref(0.0), Ref(0.0), Ref(0.0), Ref(0))
end

# Ocean-cell scan, deferred to the first apply (the met is read by then).
function _ss_init_cells!(ssb::SeaSprayBox, gc, m2)
    lsm = m2["lsm"]
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    for j in 1:Int(gc.Nj), i in 1:Int(gc.Ni)
        la, lo = indyindxtolatlon(gc, Float64(j), Float64(i))
        (ssb.lat0 <= la <= ssb.lat1) || continue
        lom = mod(lo - ssb.lon0, 360.0) + ssb.lon0     # dateline-safe box test
        (ssb.lon0 <= lom <= ssb.lon1) || continue
        Float64(lsm[i, j, 1]) < 0.5 || continue         # ocean only
        push!(ssb.ii, Int32(i)); push!(ssb.jj, Int32(j))
        push!(ssb.area, dlon * _DF_MDEG * cosd(la) * dlat * _DF_MDEG)
        push!(ssb.acc, 0.0)
    end
    ssb.initialized[] = true
    if isempty(ssb.ii)
        @warn "seaspray: no ocean cells inside the release box — source inactive"
    else
        @info "sea-spray cells initialised (first step)" ocean_cells=length(ssb.ii)
    end
    return nothing
end

"""
    _u10_from_ustar(ustar) → U10 [m/s]

Invert the neutral bulk drag law u* = √(C_d)·U10 with the OSSA C_d =
(0.8 + 0.065·U10)·10⁻³ (Large & Pond form, as used in ossa_flux_gusty).
Monotone in U10; a few fixed-point iterations converge to <1e-4 rel.
"""
@inline function _u10_from_ustar(ustar::Float64)
    ustar > 0 || return 0.0
    U = ustar / 0.035                       # C_d ≈ 1.2e-3 seed
    for _ in 1:8
        cd = (0.8 + 0.065 * U) * 1e-3
        U  = ustar / sqrt(cd)
    end
    return U
end

"""
    apply_seaspray!(state, ctx, m2raw, ssb, d, indt, dt_s; nest=1, rng) → nspawned

One source step: per ocean cell, the OSSA dry-mass flux × area × dt is added
to the cell accumulator (gust-integrated when gust_sigma > 0; the mean-field
counterfactual is tracked in parallel), and parcels of one quantum spawn at
the cell surface while the accumulator holds one.
"""
function apply_seaspray!(state::State.ParcelState, ctx, m2::AbstractDict,
                         ssb::SeaSprayBox, d::DateTime, indt::Real, dt_s::Real;
                         nest::Int = 1, rng = Random.default_rng())
    (ssb.t1 <= d <= ssb.t2) || return 0
    nestctx = ctx.nests[nest]
    gc  = nestctx.gridconfig
    ssb.initialized[] || _ss_init_cells!(ssb, gc, m2)
    Nk  = Int(nestctx.levelconfig.Nk)
    # 10 m wind: prefer archived 10u/10v; else reconstruct U10 from the friction
    # velocity u* ("ustar", e.g. ERA5 FULL_SFC zust) via the neutral drag law
    # (the surface quantity OSSA is built around) — some surface archives ship
    # u* but no 10 m winds.
    has_w10 = haskey(m2, "10u") && haskey(m2, "10v")
    ustar   = has_w10 ? nothing : (haskey(m2, "ustar") ? m2["ustar"] : nothing)
    (has_w10 || ustar !== nothing) ||
        error("seaspray: 2-D met fields carry neither 10u/10v nor ustar — cannot form U10")
    u2  = has_w10 ? m2["10u"] : ustar
    v2  = has_w10 ? m2["10v"] : nothing
    t2m = ssb.sst < 0 && haskey(m2, "2t") ? m2["2t"] : nothing
    it  = clamp(round(Int, indt), 1, size(u2, 3))
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    nsp = 0; slot = 0
    @inbounds for c in eachindex(ssb.ii)
        i = Int(ssb.ii[c]); j = Int(ssb.jj[c])
        u10 = has_w10 ? hypot(Float64(u2[i, j, it]), Float64(v2[i, j, it])) :
                        _u10_from_ustar(Float64(u2[i, j, it]))
        isfinite(u10) || continue
        sst = t2m === nothing ? abs(ssb.sst) : Float64(t2m[i, j, it])
        hs  = ssb.hs < 0 ? nothing : ssb.hs
        cd  = ssb.cd < 0 ? nothing : ssb.cd
        g   = ossa_flux_gusty(u10, ssb.gust_sigma, sst;
                              hs = hs, cd = cd, salinity_permil = ssb.salinity)
        Mm  = _ossa_massflux(g.modes_mean,  ssb.rho_dry)
        Mg  = _ossa_massflux(g.modes_gusty, ssb.rho_dry)
        dm_mean  = Mm * ssb.area[c] * dt_s
        dm_gusty = Mg * ssb.area[c] * dt_s
        ssb.emitted_mean[]  += dm_mean
        ssb.emitted_gusty[] += dm_gusty
        ssb.acc[c] += ssb.gust_sigma > 0 ? dm_gusty : dm_mean
        while ssb.acc[c] >= ssb.quantum
            slot = _df_next_free(state, slot)
            slot == 0 && return nsp                   # parcel pool exhausted
            la, lo = indyindxtolatlon(gc, Float64(j), Float64(i))
            jlo = lo + (rand(rng) - 0.5) * dlon
            jla = la + (rand(rng) - 0.5) * dlat
            indy, indx = latlontoindyindx(gc, jla, jlo)
            _write_boundary_parcel!(state, slot, nest, Float32(indx), Float32(indy),
                                    Float32(Nk) - 0.2f0,   # lowest model level (~10 m)
                                    Float32(ssb.quantum), d)
            state.relbox[slot] = eltype(state.relbox)(ssb.relbox)
            ssb.acc[c] -= ssb.quantum
            ssb.emitted[] += ssb.quantum
            ssb.nspawn[] += 1
            nsp += 1
        end
    end
    return nsp
end
