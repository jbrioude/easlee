
# ── Aerosol transport: size bins, settling, dry/wet deposition, AOD ───────────
#
# The aerosol population is discretised into size bins.  Each parcel is
# assigned ONE bin (sampled ∝ the bin mass fractions at release) — the
# standard Lagrangian treatment: the parcel's `mass` is bulk aerosol mass and
# the ensemble reproduces the size distribution in expectation.
#
# Processes (operator-split, applied after advection):
#   • gravitational settling  — Stokes–Cunningham with Sutherland viscosity
#     and an iterative Clift–Gauvin transition-regime correction (needed for
#     d ≳ 20 µm);
#   • dry deposition          — resistance model vd = vs + 1/(ra + rb + ra·rb·vs)
#     (Slinn/Seinfeld-Pandis), applied to parcels in the lowest `href` metres;
#     ra from Monin–Obukhov similarity (ustar/L derived for every dataset),
#     rb from Brownian diffusion (Sc^-2/3) + impaction (10^(-3/St));
#   • wet deposition          — size-dependent: below-cloud Slinn collection
#     efficiency E(d) (Greenfield gap), in-cloud activation fraction rising
#     with d; driven by the 3-D precipitation flux + cloud fraction
#     (ConvFields), like the bulk scheme in wetdep.jl;
#   • optical depth           — per-bin Mie extinction efficiency Q_ext(x, m)
#     (Bohren & Huffman series, computed once at setup) → column AOD on the
#     output grid:  τ = Σ N_p·Q_ext·(π/4)d² / A_cell,
#     N_p = mass / (ρ_p·π/6·d³).

# ── Refractive indices at 550 nm by aerosol type ─────────────────────────────
# κ-Köhler hygroscopicity parameter by aerosol type (Petters & Kreidenweis
# 2007): growth factor GF(RH) = (1 + κ·RH/(1−RH))^(1/3) on the dry diameter.
const AEROSOL_KAPPA = Dict{Symbol,Float32}(
    :dust    => 0.03f0,
    :sulfate => 0.61f0,
    :bc      => 0.00f0,
    :oc      => 0.10f0,
    :seasalt => 1.28f0,
    :ash     => 0.00f0,
    :water   => 0.00f0,
)

# RH grid for the per-bin wet-Mie lookup table.
const AEROSOL_RH_GRID = Float32[0.0, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98]

const AEROSOL_REFRACTIVE = Dict{Symbol,ComplexF64}(
    :dust    => 1.53 - 0.0030im,   # mineral dust (OPAC)
    :sulfate => 1.43 - 1.0e-8im,   # ammonium sulfate / H2SO4 solution
    :bc      => 1.95 - 0.79im,     # black carbon (Bond & Bergstrom 2006)
    :oc      => 1.53 - 0.006im,    # organic carbon
    :seasalt => 1.50 - 1.0e-8im,
    :ash     => 1.52 - 0.0015im,   # volcanic ash
    :water   => 1.33 - 1.0e-9im,
)

# ── Mie extinction efficiency (Bohren & Huffman, Qext only) ──────────────────
"""
    mie_qext(x, m) → Float64

Mie extinction efficiency Q_ext for size parameter `x = πd/λ` and complex
refractive index `m` (Bohren & Huffman series; either sign convention for
the imaginary part is accepted).
"""
function mie_qext(x::Float64, m::ComplexF64)
    x < 1.0e-6 && return 0.0
    # The BH series assumes the n + ik convention (absorbing ⇒ Im > 0); the
    # indices are tabulated as n − ik, so flip when needed.
    imag(m) < 0 && (m = conj(m))
    nstop = ceil(Int, x + 4.0 * x^(1/3) + 2.0)
    y    = m * x
    nmx  = ceil(Int, max(nstop, abs(y))) + 16
    # Downward recurrence for the logarithmic derivative D_n(y):
    # D_n = (n+1)/y − 1/(D_{n+1} + (n+1)/y)
    dn = zeros(ComplexF64, nmx + 1)
    for n in nmx:-1:1
        a = (n + 1) / y
        dn[n] = a - 1.0 / (dn[n+1] + a)
    end
    psi0 = cos(x);  psi1 = sin(x)
    chi0 = -sin(x); chi1 = cos(x)
    xi1  = complex(psi1, -chi1)
    qext = 0.0
    for n in 1:nstop
        fn  = (2.0n - 1.0) / x
        psi = fn * psi1 - psi0
        chi = fn * chi1 - chi0
        xi  = complex(psi, -chi)
        ta  = dn[n] / m + n / x
        an  = (ta * psi - psi1) / (ta * xi - xi1)
        tb  = dn[n] * m + n / x
        bn  = (tb * psi - psi1) / (tb * xi - xi1)
        qext += (2.0n + 1.0) * real(an + bn)
        psi0 = psi1; psi1 = psi
        chi0 = chi1; chi1 = chi
        xi1  = xi
    end
    return 2.0 / x^2 * qext
end

# ── Size-bin definition ───────────────────────────────────────────────────────
struct AerosolBins
    d::Vector{Float32}          # DRY bin mass-mean diameters [m]
    mass_frac::Vector{Float32}  # mass fraction per bin (sums to 1)
    rho_p::Float32              # dry particle density [kg/m³]
    m_ri::ComplexF64            # dry refractive index at `lambda`
    lambda::Float32             # wavelength for optics [m]
    kappa::Float32              # κ-Köhler hygroscopicity (0 = no growth)
    qext::Matrix{Float32}       # Mie Qext, (nbins × length(AEROSOL_RH_GRID)):
                                # wet diameter + water-mixed refractive index
    species::Symbol
end

# Hygroscopic growth factor on the dry diameter (κ-Köhler, a_w ≈ RH).
"""
    growth_factor(kappa, rh) → Float32

κ-Köhler hygroscopic diameter growth factor d_wet/d_dry at relative humidity
`rh` (0–1, capped at 0.98); `kappa ≤ 0` → 1 (no growth).
"""
@inline function growth_factor(kappa::Float32, rh::Float32)
    kappa <= 0.0f0 && return 1.0f0
    rh = clamp(rh, 0.0f0, 0.98f0)
    return cbrt(1.0f0 + kappa * rh / max(1.0f0 - rh, 0.02f0))
end

# Wet particle density from volume mixing with water.
@inline function _wet_density(rho_dry::Float32, gf::Float32)
    gf3 = gf^3
    return (rho_dry + 1000.0f0 * (gf3 - 1.0f0)) / gf3
end

function _finish_bins(d, mf, rho_p, species, m_ri, lambda, kappa)
    mf  = Float32.(mf ./ sum(mf))
    ri  = m_ri === nothing ? AEROSOL_REFRACTIVE[species] : ComplexF64(m_ri)
    kap = kappa === nothing ? get(AEROSOL_KAPPA, species, 0.0f0) : Float32(kappa)
    # Wet-Mie lookup: for each bin and RH-grid point, grow the diameter and
    # volume-mix the refractive index with water (1.33).
    nb  = length(d)
    qe  = zeros(Float32, nb, length(AEROSOL_RH_GRID))
    m_w = 1.33 - 1.0e-9im
    for (jr, rh) in enumerate(AEROSOL_RH_GRID)
        gf = growth_factor(kap, rh)
        fw = 1.0 - 1.0 / Float64(gf)^3          # water volume fraction
        mi = ri * (1.0 - fw) + m_w * fw
        for b in 1:nb
            x = Float64(π) * Float64(d[b]) * Float64(gf) / Float64(lambda)
            qe[b, jr] = Float32(mie_qext(x, mi))
        end
    end
    return AerosolBins(Float32.(d), mf, Float32(rho_p), ri, Float32(lambda), kap, qe, species)
end

"""
    AerosolBins(; dg, sigma_g, nbins=8, rho_p=2000.0, species=:dust,
                  m_ri=nothing, lambda=550e-9)

Lognormal MASS distribution: median diameter `dg` [m], geometric std `sigma_g`.
Bins are log-spaced over ±3 ln σg and weighted by the analytic lognormal mass
fractions.
"""
function AerosolBins(; dg::Real, sigma_g::Real, nbins::Int=8,
                       rho_p::Real=2000.0, species::Symbol=:dust,
                       m_ri=nothing, lambda::Real=550.0e-9, kappa=nothing)
    lnsg  = log(sigma_g)
    edges = exp.(range(log(dg) - 3lnsg, log(dg) + 3lnsg, length=nbins + 1))
    d  = sqrt.(edges[1:end-1] .* edges[2:end])                 # geometric mid
    # mass fraction between edges from the lognormal CDF
    cdf(x) = 0.5 * (1.0 + erf((log(x) - log(dg)) / (sqrt(2.0) * lnsg)))
    mf = [cdf(edges[k+1]) - cdf(edges[k]) for k in 1:nbins]
    return _finish_bins(d, mf, rho_p, species, m_ri, lambda, kappa)
end

"""
    AerosolBins(edges; mass_frac=nothing, rho_p=2000.0, species=:dust, ...)

Explicit size bins from `edges` [m] (length nbins+1).  `mass_frac` defaults
to equal mass per bin.
"""
function AerosolBins(edges::AbstractVector; mass_frac=nothing,
                     rho_p::Real=2000.0, species::Symbol=:dust,
                     m_ri=nothing, lambda::Real=550.0e-9, kappa=nothing)
    nb = length(edges) - 1
    d  = sqrt.(edges[1:end-1] .* edges[2:end])
    mf = mass_frac === nothing ? fill(1.0 / nb, nb) : collect(mass_frac)
    return _finish_bins(d, mf, rho_p, species, m_ri, lambda, kappa)
end

# erf without pulling in SpecialFunctions: Abramowitz & Stegun 7.1.26 (|ε|<1.5e-7)
function erf(x::Float64)
    s = sign(x); x = abs(x)
    t = 1.0 / (1.0 + 0.3275911x)
    y = 1.0 - (((((1.061405429t - 1.453152027)t) + 1.421413741)t - 0.284496736)t + 0.254829592)t * exp(-x^2)
    return s * y
end

# ── Per-run aerosol state ─────────────────────────────────────────────────────
"""
Per-run aerosol state: the size-bin definition, each parcel's assigned bin,
and the per-nest dry/wet deposition accumulation grids [kg per met cell].
`href` is the dry-deposition application depth [m].
"""
mutable struct AerosolModel
    bins::AerosolBins
    bin::Vector{Int16}                 # per-parcel bin index (0 = unassigned)
    # Per-parcel aerosol mass [kg] for EMISSION-DRIVEN aerosol (domainfill air
    # parcels whose state.mass is the AIR mass).  EMPTY ⇒ use state.mass (the
    # release-box convention, where the parcel IS the aerosol).
    pmass::Vector{Float64}
    drydep::Vector{Matrix{Float64}}    # per-nest deposited mass (met grid)
    wetdep::Vector{Matrix{Float64}}
    href::Float32                      # dry-deposition layer depth [m]
    # Release boxes the aerosol physics applies to (mixed gas+aerosol runs);
    # EMPTY = every box (legacy single-species behaviour).
    boxon::Vector{Bool}
    # Boxes whose species has [species].settling = true (EMPTY = all settle).
    # Was silently ignored before 2026-07-18 — settling ran unconditionally,
    # which the sulfate-descent no-settling twin exposed.
    settleon::Vector{Bool}
end

@inline _aero_on(aero::AerosolModel, b::Integer) =
    isempty(aero.boxon) || (1 <= b <= length(aero.boxon) && aero.boxon[b])

@inline _aero_settles(aero::AerosolModel, b::Integer) =
    isempty(aero.settleon) || (1 <= b <= length(aero.settleon) && aero.settleon[b])

# Aerosol mass of parcel i: the emission-driven pmass when present, else state.mass.
@inline _aero_mass(aero::AerosolModel, state, i::Integer) =
    isempty(aero.pmass) ? Float64(state.mass[i]) : aero.pmass[i]

# Remove `dm` kg from parcel i's aerosol mass (the pmass array or state.mass).
@inline function _aero_mass_dec!(aero::AerosolModel, state, i::Integer, dm::Real)
    if isempty(aero.pmass)
        state.mass[i] -= eltype(state.mass)(dm)
    else
        aero.pmass[i] = max(aero.pmass[i] - Float64(dm), 0.0)
    end
end

function AerosolModel(n_max::Integer, ctx, bins::AerosolBins; href::Real=30.0,
                      boxon::Vector{Bool}=Bool[], settleon::Vector{Bool}=Bool[],
                      emit::Bool=false)
    AerosolModel(bins, zeros(Int16, Int(n_max)),
        emit ? zeros(Float64, Int(n_max)) : Float64[],
        [zeros(Float64, Int(n.gridconfig.Ni), Int(n.gridconfig.Nj)) for n in ctx.nests],
        [zeros(Float64, Int(n.gridconfig.Ni), Int(n.gridconfig.Nj)) for n in ctx.nests],
        Float32(href), boxon, settleon)
end

# Sample a size bin for every active parcel that doesn't have one yet
# (newly released).  Probability ∝ bin mass fraction.
function assign_aerosol_bins!(aero::AerosolModel, state)
    mf = aero.bins.mass_frac
    nb = length(mf)
    @inbounds for i in eachindex(state.status)
        (state.status[i] == 0x01 && aero.bin[i] == 0) || continue
        _aero_on(aero, Int(state.relbox[i])) || continue
        u = _conv_rand(state, i)
        acc = 0.0f0
        b = nb
        for k in 1:nb
            acc += mf[k]
            if u <= acc
                b = k
                break
            end
        end
        aero.bin[i] = Int16(b)
    end
    return nothing
end

# Relative humidity at the parcel, reconstructed from the typed fields
# (q from Tv = p/(R_d·ρ), saturation from Magnus over water).
@inline function _parcel_rh(m3, indx::Float32, indy::Float32, indl::Float32, indt)
    T   = interpolation4Dsingle(m3.T,   indx, indy, indl, indt)
    p   = interpolation4Dsingle(m3.puv, indx, indy, indl, indt)
    rho = interpolation4Dsingle(m3.rho, indx, indy, indl, indt)
    (isfinite(T) && isfinite(p) && isfinite(rho) && T > 150.0f0 && rho > 1.0f-6) ||
        return 0.0f0
    Tv = p / (287.05f0 * rho)
    q  = clamp((Tv / T - 1.0f0) / 0.608f0, 0.0f0, 0.05f0)
    tc = T - 273.15f0
    es = 611.2f0 * exp(17.67f0 * tc / (tc + 243.5f0))
    e  = q * p / (0.622f0 + 0.378f0 * q)
    return clamp(e / es, 0.0f0, 0.98f0)
end

# ── Terminal settling velocity ────────────────────────────────────────────────
"""
    settling_velocity(d, rho_p, T, p) → Float32 [m/s]

Terminal fall speed of a particle of diameter `d` [m] and density `rho_p`
[kg/m³] at temperature `T` [K] and pressure `p` [Pa]: Stokes–Cunningham with
Sutherland viscosity, iterated with Clift–Gauvin drag in the transition
regime (Re > 0.01, relevant for d ≳ 20 µm).
"""
@inline function settling_velocity(d::Float32, rho_p::Float32, T::Float32, p::Float32)
    mu   = 1.458f-6 * T * sqrt(T) / (T + 110.4f0)            # Sutherland [Pa s]
    lam  = 6.6f-8 * (101325.0f0 / max(p, 100.0f0)) * (T / 293.15f0)
    kn   = 2.0f0 * lam / d
    cc   = 1.0f0 + kn * (1.257f0 + 0.4f0 * exp(-1.1f0 / kn)) # Cunningham slip
    rhoa = p / (287.05f0 * T)
    vs   = rho_p * d * d * 9.81f0 * cc / (18.0f0 * mu)       # Stokes
    # Transition regime (Re > 0.01): Clift–Gauvin drag, fixed-point iteration
    @inbounds for _ in 1:6
        re = vs * d * rhoa / mu
        re < 0.01f0 && break
        cd = 24.0f0 / re * (1.0f0 + 0.15f0 * re^0.687f0)
        vs_new = sqrt(4.0f0 * 9.81f0 * d * rho_p * cc / (3.0f0 * cd * rhoa))
        conv = abs(vs_new - vs) < 1.0f-3 * vs
        vs = vs_new
        conv && break
    end
    return vs
end

"""
    apply_settling!(aero, state, ctx, dt_s, t_array)

Gravitational settling for all active parcels (z ← z − v_t·Δt), with v_t from
the parcel's size bin and local T/p.  Assigns bins to new parcels first.
"""
function apply_settling!(aero::AerosolModel, state, ctx, dt_s::Real, t_array)
    assign_aerosol_bins!(aero, state)
    dt = Float32(dt_s)
    for n in 1:length(ctx.nests)
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8 &&
                           _aero_on(aero, Int(state.relbox[i])) &&
                           _aero_settles(aero, Int(state.relbox[i])),
                      eachindex(state.status))
        isempty(idx) && continue
        m3   = nestctx.metfields3D
        nz   = size(m3.zlev, 3)
        Ni   = Int(nestctx.gridconfig.Ni)
        Nj   = Int(nestctx.gridconfig.Nj)
        indt = t_array[n]
        nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
        zbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
        Threads.@threads :static for i in idx
            tid  = Threads.threadid()
            indx = state.indx[i]
            indy = state.indy[i]
            (indx < 1.0f0 || indx >= Float32(Ni) ||
             indy < 1.0f0 || indy >= Float32(Nj)) && continue
            b = aero.bin[i]
            b == 0 && continue
            dbin = aero.bins.d[b]
            indl = clamp(state.indl[i], 1.0f0, Float32(nz - 1))
            T = interpolation4Dsingle(m3.T,   indx, indy, indl, indt)
            p = interpolation4Dsingle(m3.puv, indx, indy, indl, indt)
            (isfinite(T) && isfinite(p) && T > 150.0f0 && p > 10.0f0) || continue
            # Hygroscopic growth: settle with the WET diameter and density.
            gf = growth_factor(aero.bins.kappa,
                               _parcel_rh(m3, indx, indy, indl, indt))
            vt = settling_velocity(dbin * gf, _wet_density(aero.bins.rho_p, gf), T, p)
            z  = zbufs[tid]
            get_vertical_profile!(z, m3.zlev, indx, indy, indt, nz)
            z_new = _z_at(z, indl, nz) - vt * dt
            isfinite(z_new) || continue
            state.indl[i] = _z_to_indl(z, max(z_new, 0.1f0), nz)
        end
    end
    return nothing
end

# ── Dry deposition (resistance model) ─────────────────────────────────────────
@inline function _psi_h(zeta::Float32)
    if zeta < 0.0f0
        a = sqrt(max(1.0f0 - 16.0f0 * zeta, 1.0f0))
        return 2.0f0 * log((1.0f0 + a) / 2.0f0)
    end
    return -5.0f0 * zeta
end

"""
    deposition_velocity(d, rho_p, T, p, ustar, ol, href; z0=0.03) → Float32 [m/s]

Aerosol dry-deposition velocity vd = vs + 1/(ra + rb + ra·rb·vs)
(Slinn / Seinfeld–Pandis resistance form): ra from Monin–Obukhov similarity
between `z0` and `href` [m] using `ustar` [m/s] and Obukhov length `ol` [m];
rb from Brownian diffusion (Sc^-2/3) + impaction (10^(-3/St)); vs = settling.
"""
@inline function deposition_velocity(d::Float32, rho_p::Float32, T::Float32, p::Float32,
                                     ustar::Float32, ol::Float32, href::Float32;
                                     z0::Float32 = 0.03f0)
    mu   = 1.458f-6 * T * sqrt(T) / (T + 110.4f0)
    rhoa = p / (287.05f0 * T)
    nu   = mu / rhoa
    lam  = 6.6f-8 * (101325.0f0 / max(p, 100.0f0)) * (T / 293.15f0)
    cc   = 1.0f0 + 2.0f0 * lam / d * (1.257f0 + 0.4f0 * exp(-1.1f0 * d / (2.0f0 * lam)))
    vs   = settling_velocity(d, rho_p, T, p)
    us   = max(ustar, 0.05f0)
    # Aerodynamic resistance, MO similarity between z0 and href.
    # |L| floored at 1 m with the SIGN preserved (sign(0)=0 would give
    # ψ_h(±Inf) and ra = Inf − Inf = NaN).
    aol  = abs(ol) < 1.0f0 ? (ol < 0.0f0 ? -1.0f0 : 1.0f0) : ol
    ra   = (log(href / z0) - _psi_h(href / aol) + _psi_h(z0 / aol)) / (0.4f0 * us)
    isfinite(ra) || (ra = 100.0f0)
    ra   = clamp(ra, 1.0f0, 1.0f4)
    # Quasi-laminar resistance: Brownian diffusion + impaction (Slinn 1980)
    Db   = 1.380649f-23 * T * cc / (3.0f0 * Float32(π) * mu * d)
    sc   = nu / Db
    st   = vs * us * us / (9.81f0 * nu)
    rb   = 1.0f0 / (us * (sc^(-2.0f0/3.0f0) + 10.0f0^(-3.0f0 / max(st, 0.05f0))))
    rb   = clamp(rb, 1.0f0, 1.0f5)
    return vs + 1.0f0 / (ra + rb + ra * rb * vs)
end

"""
    dry_deposition_aerosol!(aero, state, ctx, dt_s, t_array)

Resistance-model dry deposition for parcels in the lowest `aero.href` metres:
mass decays as exp(−vd·Δt/href); removed mass accumulates on `aero.drydep[n]`.
Returns total mass removed this step.
"""
function dry_deposition_aerosol!(aero::AerosolModel, state, ctx, dt_s::Real, t_array)
    Int(ctx.runconfig.direction) == -1 && return 0.0
    dt = Float32(dt_s)
    total = Threads.Atomic{Float64}(0.0)
    for n in 1:length(ctx.nests)
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8 &&
                           _aero_on(aero, Int(state.relbox[i])),
                      eachindex(state.status))
        isempty(idx) && continue
        m3   = nestctx.metfields3D
        m2   = nestctx.metfields2D
        nz   = size(m3.zlev, 3)
        Ni   = Int(nestctx.gridconfig.Ni)
        Nj   = Int(nestctx.gridconfig.Nj)
        indt = t_array[n]
        dep  = aero.drydep[n]
        nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
        zbufs   = [Vector{Float32}(undef, nz) for _ in 1:nthr]
        removed = zeros(Float64, nthr)
        Threads.@threads :static for i in idx
            tid  = Threads.threadid()
            indx = state.indx[i]
            indy = state.indy[i]
            (indx < 1.0f0 || indx >= Float32(Ni) ||
             indy < 1.0f0 || indy >= Float32(Nj)) && continue
            b = aero.bin[i]
            b == 0 && continue
            z = zbufs[tid]
            get_vertical_profile!(z, m3.zlev, indx, indy, indt, nz)
            zp = _z_at(z, clamp(state.indl[i], 1.0f0, Float32(nz - 1)), nz)
            zp > aero.href && continue
            # lowest full level = slot 2, read at indl 1.5 (layout contract)
            T  = interpolation4Dsingle(m3.T,   indx, indy, 1.5f0, indt)
            p  = interpolation4Dsingle(m3.puv, indx, indy, 1.5f0, indt)
            us = interpolate2Dsingle(m2.ustar, indx, indy, indt)
            ol = interpolate2Dsingle(m2.ol,    indx, indy, indt)
            (isfinite(T) && isfinite(p) && T > 150.0f0) || continue
            isfinite(us) || (us = 0.2f0)
            isfinite(ol) || (ol = 1.0f6)
            gf = growth_factor(aero.bins.kappa, _parcel_rh(m3, indx, indy, 1.5f0, indt))
            vd = deposition_velocity(aero.bins.d[b] * gf, _wet_density(aero.bins.rho_p, gf),
                                     T, p, us, ol, aero.href)
            frem = -expm1(-vd * dt / aero.href)
            dm = _aero_mass(aero, state, i) * frem
            _aero_mass_dec!(aero, state, i, dm)
            removed[tid] += Float64(dm)
            ii = clamp(round(Int, indx), 1, Ni)
            jj = clamp(round(Int, indy), 1, Nj)
            Atomix.@atomic dep[ii, jj] += Float64(dm)
        end
        Threads.atomic_add!(total, sum(removed))
    end
    return total[]
end

# ── Size-dependent wet deposition ─────────────────────────────────────────────
# Below-cloud collection efficiency E(d): Brownian capture below ~0.1 µm,
# Greenfield gap minimum near 0.5–1 µm, impaction above (Slinn-shaped fit).
@inline function _collection_eff(d::Float32)
    dum = d * 1.0f6     # µm
    return clamp(9.0f-4 * dum^(-1.2f0) + 9.0f-4 * dum^2.5f0, 1.0f-3, 1.0f0)
end

# In-cloud activated (nucleation-scavenged) fraction: rises with size, ~0.45
# at the 0.1 µm activation diameter, → 0.9 for accumulation/coarse modes.
@inline function _activated_frac(d::Float32)
    d3 = d^3
    return 0.9f0 * d3 / (d3 + (1.0f-7)^3)
end

"""
    wet_deposition_aerosol!(aero, state, ctx, convfields, dt_s, t_array;
                            cc_threshold=0.05, I_min=0.01)

Size-dependent wet scavenging from the 3-D precipitation flux (tpf) and cloud
fraction (cc) in `convfields`:
  in-cloud   : Λ = f_act(d) · 3.36e-4 · I^0.79          (nucleation)
  below-cloud: Λ = 1.5 · E(d) · I_SI / d_drop            (impaction, d_drop=1 mm)
Removed mass accumulates on `aero.wetdep[n]`.  Returns total removed.
"""
function wet_deposition_aerosol!(aero::AerosolModel, state, ctx, convfields,
                                 dt_s::Real, t_array;
                                 cc_threshold::Float32 = 0.05f0,
                                 I_min::Float32 = 0.01f0)
    Int(ctx.runconfig.direction) == -1 && return 0.0
    dt = Float32(dt_s)
    total = Threads.Atomic{Float64}(0.0)
    for n in 1:length(ctx.nests)
        cf = convfields[n]
        cf === nothing && continue
        nestctx = ctx.nests[n]
        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8 &&
                           _aero_on(aero, Int(state.relbox[i])),
                      eachindex(state.status))
        isempty(idx) && continue
        m3   = nestctx.metfields3D
        Ni   = Int(nestctx.gridconfig.Ni)
        Nj   = Int(nestctx.gridconfig.Nj)
        nz   = size(cf.cc, 3)
        indt = t_array[n]
        dep  = aero.wetdep[n]
        nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
        ccbufs  = [Vector{Float32}(undef, nz) for _ in 1:nthr]
        removed = zeros(Float64, nthr)
        Threads.@threads :static for i in idx
            tid  = Threads.threadid()
            indx = state.indx[i]
            indy = state.indy[i]
            (indx < 1.0f0 || indx >= Float32(Ni) ||
             indy < 1.0f0 || indy >= Float32(Nj)) && continue
            b = aero.bin[i]
            b == 0 && continue
            tpf = interpolation4Dsingle(cf.tpf, indx, indy,
                                        clamp(state.indl[i], 1.0f0, Float32(nz - 1)), indt)
            (isfinite(tpf) && tpf * 3600.0f0 > I_min) || continue
            I_mmh = tpf * 3600.0f0
            ccp = ccbufs[tid]
            get_vertical_profile!(ccp, cf.cc, indx, indy, indt, nz)
            # cloud fraction of the layer containing the parcel: full level k
            # spans zlev[k-1]..zlev[k] (layout contract) → slot floor(indl) + 1
            k0 = clamp(Base.unsafe_trunc(Int, state.indl[i]) + 1, 2, nz)
            cc_here = ccp[k0]
            cc_max  = 0.0f0
            @inbounds for k in 1:nz
                ccp[k] > cc_max && (cc_max = ccp[k])
            end
            gf   = growth_factor(aero.bins.kappa,
                                 _parcel_rh(m3, indx, indy,
                                            clamp(state.indl[i], 1.0f0, Float32(nz - 1)), indt))
            dbin = aero.bins.d[b] * gf
            lam = if isfinite(cc_here) && cc_here > cc_threshold
                _activated_frac(dbin) * 3.36f-4 * I_mmh^0.79f0          # in-cloud
            else
                # Λ = (3/2)·E·I/d_drop, I = tpf·1e-3 [m/s], d_drop = 1 mm
                # → Λ = 1.5·E·tpf  [1/s]
                1.5f0 * _collection_eff(dbin) * tpf                     # below-cloud
            end
            fp   = clamp(isfinite(cc_max) ? cc_max : 0.05f0, 0.05f0, 1.0f0)
            frem = fp * -expm1(-lam * dt / fp)     # sub-grid precipitating fraction
            dm = _aero_mass(aero, state, i) * frem
            _aero_mass_dec!(aero, state, i, dm)
            removed[tid] += Float64(dm)
            ii = clamp(round(Int, indx), 1, Ni)
            jj = clamp(round(Int, indy), 1, Nj)
            Atomix.@atomic dep[ii, jj] += Float64(dm)
        end
        Threads.atomic_add!(total, sum(removed))
    end
    return total[]
end

# ── Aerosol optical depth on the output grid ──────────────────────────────────
"""
    compute_aod(aero, state, ctx, t_array) → Matrix{Float64} (Niout × Njout)

Column aerosol optical depth at `aero.bins.lambda`, INCLUDING hygroscopic
growth: per parcel, the dry-mass particle number N = m/(ρ_dry·(π/6)·d_dry³)
extinguishes with the WET cross-section Q_ext(d_wet, m_wet)·(π/4)·d_wet²
(Q_ext from the per-bin RH lookup built at setup).
"""
function compute_aod(aero::AerosolModel, state, ctx, t_array)
    oc    = ctx.output
    Niout = Int(oc.Niout)
    Njout = Int(oc.Njout)
    aod   = zeros(Float64, Niout, Njout)
    bins  = aero.bins
    nrh   = length(AEROSOL_RH_GRID)
    for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        _aero_on(aero, Int(state.relbox[i])) || continue
        b = aero.bin[i]
        b == 0 && continue
        n  = Int(state.nest[i])
        gc = ctx.nests[n].gridconfig
        m3 = ctx.nests[n].metfields3D
        nz = size(m3.zlev, 3)
        rh = _parcel_rh(m3, state.indx[i], state.indy[i],
                        clamp(state.indl[i], 1.0f0, Float32(nz - 1)), t_array[n])
        gf = growth_factor(bins.kappa, rh)
        # Qext interpolated on the RH grid
        jr = 1
        while jr < nrh - 1 && AEROSOL_RH_GRID[jr+1] < rh
            jr += 1
        end
        fr = clamp((rh - AEROSOL_RH_GRID[jr]) /
                   max(AEROSOL_RH_GRID[jr+1] - AEROSOL_RH_GRID[jr], 1.0f-3), 0.0f0, 1.0f0)
        qe = (1 - fr) * bins.qext[b, jr] + fr * bins.qext[b, jr+1]
        # extinction cross-section per unit DRY mass [m²/kg]
        sigma_m = Float64(qe) * (π/4) * (Float64(bins.d[b]) * gf)^2 /
                  (Float64(bins.rho_p) * (π/6) * Float64(bins.d[b])^3)
        xout = indxtoxout(gc, state.indx[i], state.indy[i])
        yout = indytoyout(gc, state.indx[i], state.indy[i])
        ii = Base.unsafe_trunc(Int, xout)
        jj = Base.unsafe_trunc(Int, yout)
        (1 <= ii <= Niout && 1 <= jj <= Njout) || continue
        aod[ii, jj] += _aero_mass(aero, state, i) * sigma_m / Float64(oc.area[ii, jj])
    end
    return aod
end

# Backward-compatible method (dry/last-bracket assumption: RH from slot 1).
compute_aod(aero::AerosolModel, state, ctx) =
    compute_aod(aero, state, ctx, fill(1.0f0, length(ctx.nests)))
