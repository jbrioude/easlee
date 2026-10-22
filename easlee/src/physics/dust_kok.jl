# physics/dust_kok.jl
# ──────────────────────────────────────────────────────────────────────────────
# Mineral-dust emission after Kok et al. (2014), "An improved dust emission model –
# Part 1: Model description and comparison against measurements", ACP 14, 13023–13041,
# with the brittle-fragmentation emitted size distribution of Kok (2011), PNAS 108, 1016.
#
# Vertical dust flux (Kok et al. 2014, Eq. 18; u* > u*_t, else 0):
#
#     F_d = C_d · f_bare · f_clay · ρ_a · (u*² − u*_t²) / u*_st · (u* / u*_t)^(C_α · (u*_st − u*_st0) / u*_st0)
#     C_d = C_d0 · exp(−C_e · (u*_st − u*_st0) / u*_st0)
#
# with u*_t the fluid threshold friction velocity, u*_st = u*_t · sqrt(ρ_a/ρ_a0) the
# "standardised" threshold (ρ_a0 = 1.225 kg m⁻³), and the constants of their Table 1:
# C_d0 = 4.4e-5 (±0.5e-5), C_e = 2.0 (±0.3), C_α = 2.7 (±1.0), u*_st0 = 0.16 m s⁻¹.
# The soil-state dependence is carried by u*_st (soil erodibility increases as the
# threshold decreases), the clay fraction f_clay (the dust-producing fraction of the
# saltating soil), and the bare-soil fraction f_bare (vegetation / roughness masking).
#
# Fluid threshold: Shao & Lu (2000), JGR 105, 22437, for the most erodible particle size
# (D ≈ 127 µm, Kok et al. 2014):
#     u*_t0 = sqrt( A_N · (ρ_p g D / ρ_a + γ / (ρ_a D)) ),   A_N = 0.0123, γ = 1.65e-4 kg s⁻²
# and the soil-moisture increase of Fécan, Marticorena & Bergametti (1999), Ann. Geophys.
# 17, 149:  u*_t = u*_t0 · sqrt(1 + 1.21 (w − w′)^0.68)  for w > w′,  w′ = 0.0014 c² + 0.17 c
# (w = gravimetric soil moisture in %, c = clay content in %).
#
# Emitted size distribution (Kok 2011, Eq. 6), volume per ln D, 0.2 ≤ D ≤ 20 µm:
#     dV/dlnD = (D / c_V) · [1 + erf( ln(D/D_s) / (√2 ln σ_s) )] · exp(−(D/λ)³)
#     D_s = 3.4 µm, σ_s = 3.0, λ = 12 µm, c_V = 12.62 µm (normalisation over 0.2–20 µm)
# — independent of wind speed and soil, which is the point of the brittle-fragmentation
# theory; used here to set the mass fractions of the species' size bins (species
# `size.distribution = "kok"`).
#
# Engine source term ([[release]] kind = "dust"): like the sea-spray source, the release
# box only reserves a relbox; `build_dust` records the box and `apply_dust!` runs each
# transport step over the LAND cells of the box (lsm ≥ 0.5), accumulating F_d · area · Δt
# per cell and spawning a parcel of one mass quantum at the surface when the accumulator
# fills.  Inputs per cell: u* (2-D met "ustar"), air density from surface pressure and 2-m
# temperature when present (else 1.2 kg m⁻³), soil moisture from the 2-D met "swvl1"
# (volumetric, converted with a bulk density of 1.3 g cm⁻³) when present, else the
# configured constant; clay and bare fractions constant or from lon/lat maps.
# ──────────────────────────────────────────────────────────────────────────────

const KOK_CD0    = 4.4e-5        # dimensionless
const KOK_CE     = 2.0
const KOK_CALPHA = 2.7
const KOK_UST0   = 0.16          # m/s
const KOK_RHOA0  = 1.225         # kg/m³
const SHAO_AN    = 0.0123
const SHAO_GAMMA = 1.65e-4       # kg/s²
const DUST_RHO_P = 2650.0        # kg/m³ (quartz-like soil grains, threshold only)
const KOK11_DS   = 3.4e-6        # m
const KOK11_SIGS = 3.0
const KOK11_LAM  = 12.0e-6       # m
const KOK11_CV   = 12.62e-6      # m
const KOK11_DMIN = 0.2e-6
const KOK11_DMAX = 20.0e-6

"""
    dust_threshold_dry(rho_a; D=127e-6, rho_p=2650) → u*_t0 [m/s]

Shao & Lu (2000) fluid threshold for grains of diameter `D` in air of density `rho_a`.
≈ 0.21 m/s at 1.225 kg/m³ for D = 127 µm.
"""
@inline function dust_threshold_dry(rho_a::Real; D::Real = 127.0e-6, rho_p::Real = DUST_RHO_P)
    ra = max(Float64(rho_a), 0.1)
    return sqrt(SHAO_AN * (Float64(rho_p) * 9.80665 * Float64(D) / ra + SHAO_GAMMA / (ra * Float64(D))))
end

"""
    fecan_moisture_factor(w_pct, clay_frac) → f_w ≥ 1

Fécan et al. (1999): threshold multiplier for gravimetric soil moisture `w_pct` [%] and
clay fraction `clay_frac` (0–1).  1 below the residual moisture w′.
"""
@inline function fecan_moisture_factor(w_pct::Real, clay_frac::Real)
    c  = 100.0 * clamp(Float64(clay_frac), 0.0, 1.0)
    wp = 0.0014 * c^2 + 0.17 * c
    w  = max(Float64(w_pct), 0.0)
    return w > wp ? sqrt(1.0 + 1.21 * (w - wp)^0.68) : 1.0
end

"""
    dust_flux_kok(ustar, rho_a; clay=0.2, bare=1.0, w_pct=0.0, D=127e-6,
                  cd0=KOK_CD0, ce=KOK_CE, calpha=KOK_CALPHA) → F_d [kg m⁻² s⁻¹]

Kok et al. (2014) vertical dust flux (Eq. 18) for friction velocity `ustar` and air
density `rho_a`; 0 when u* ≤ u*_t.
"""
function dust_flux_kok(ustar::Real, rho_a::Real; clay::Real = 0.2, bare::Real = 1.0,
                       w_pct::Real = 0.0, D::Real = 127.0e-6,
                       cd0::Real = KOK_CD0, ce::Real = KOK_CE, calpha::Real = KOK_CALPHA,
                       ut::Union{Nothing,Real} = nothing)
    us = Float64(ustar); ra = Float64(rho_a)
    (isfinite(us) && isfinite(ra) && us > 0 && ra > 0) || return 0.0
    # `ut` = an explicit fluid threshold [m/s] (verification against the published curves at
    # fixed u*_st); default = Shao & Lu (2000) at D with the Fécan et al. (1999) moisture increase
    ut  = ut === nothing ? dust_threshold_dry(ra; D = D) * fecan_moisture_factor(w_pct, clay) : Float64(ut)
    us <= ut && return 0.0
    ust = ut * sqrt(ra / KOK_RHOA0)
    x   = (ust - KOK_UST0) / KOK_UST0
    cd  = Float64(cd0) * exp(-Float64(ce) * x)
    return cd * clamp(Float64(bare), 0.0, 1.0) * clamp(Float64(clay), 0.0, 1.0) * ra *
           (us^2 - ut^2) / ust * (us / ut)^(Float64(calpha) * x)
end

"""
    kok_dvdlnd(D) → dV/dlnD  [m³ per m³ of emitted dust, per unit ln D]

Kok (2011) emitted volume size distribution (Eq. 6); unnormalised units cancel in
`kok_mass_fractions`.
"""
@inline function kok_dvdlnd(D::Real)
    d = Float64(D)
    return (d / KOK11_CV) * (1.0 + erf(log(d / KOK11_DS) / (sqrt(2.0) * log(KOK11_SIGS)))) *
           exp(-(d / KOK11_LAM)^3)
end

"""
    kok_mass_fractions(edges) → Vector (sums to 1)

Mass fraction of each bin [edges[k], edges[k+1]] under the Kok (2011) distribution,
integrated in ln D (401-point trapezoid per bin).  Bins outside 0.2–20 µm receive the
distribution's own (small) tail values; the vector is renormalised over the given bins.
"""
function kok_mass_fractions(edges::AbstractVector)
    nb = length(edges) - 1
    nb >= 1 || error("kok_mass_fractions: need ≥ 2 edges")
    mf = zeros(Float64, nb)
    for k in 1:nb
        a = log(Float64(edges[k])); b = log(Float64(edges[k+1]))
        xs = range(a, b; length = 401)
        ys = [kok_dvdlnd(exp(x)) for x in xs]
        mf[k] = (b - a) / 400 * (0.5 * (ys[1] + ys[end]) + sum(ys[2:end-1]))
    end
    s = sum(mf); s > 0 || error("kok_mass_fractions: zero mass in the given bins")
    return mf ./ s
end

# ── Engine source term ([[release]] kind = "dust") ─────────────────────────────
struct DustBox
    relbox::Int32
    lon0::Float64; lon1::Float64
    lat0::Float64; lat1::Float64
    initialized::Base.RefValue{Bool}
    ii::Vector{Int32}; jj::Vector{Int32}     # land met cells inside the box
    area::Vector{Float64}                    # [m²]
    acc::Vector{Float64}                     # unspawned accumulated mass [kg]
    clay::Vector{Float64}                    # per cell (filled at init)
    bare::Vector{Float64}
    quantum::Float64                         # kg per spawned parcel
    t1::DateTime; t2::DateTime
    clay0::Float64                           # constants when no map
    bare0::Float64
    w_pct::Float64                           # gravimetric soil moisture [%] when no met swvl1
    D::Float64                               # threshold grain diameter [m]
    cd0::Float64; ce::Float64; calpha::Float64
    clay_map::Any                            # nothing | (lons, lats, vals)
    bare_map::Any
    emitted::Base.RefValue{Float64}
    nspawn::Base.RefValue{Int}
    flux_max::Base.RefValue{Float64}         # diagnostic: max cell flux seen [kg m⁻² s⁻¹]
end

function build_dust(ctx, m2::AbstractDict, relbox::Integer, r; clay_map = nothing, bare_map = nothing)
    haskey(m2, "lsm") ||
        (@warn "dust: met has no land-sea mask (lsm) — source disabled"; return nothing)
    haskey(m2, "ustar") ||
        (@warn "dust: met has no friction velocity (ustar) — source disabled"; return nothing)
    q = r.mass_per_parcel_kg
    q > 0 || error("dust release \"$(r.species)\": mass_per_parcel must be set (the spawn quantum)")
    return DustBox(Int32(relbox), minimum(r.lon), maximum(r.lon), minimum(r.lat), maximum(r.lat),
                   Ref(false), Int32[], Int32[], Float64[], Float64[], Float64[], Float64[], q,
                   r.time[1], r.time[2], r.dust_clay, r.dust_bare, r.dust_w_pct, r.dust_d_m,
                   r.dust_cd0, r.dust_ce, r.dust_calpha, clay_map, bare_map,
                   Ref(0.0), Ref(0), Ref(0.0))
end

# nearest-neighbour sample of a (lons, lats, vals) map; lons may be 0–360 or −180–180
@inline function _map_sample(mp, lon::Float64, lat::Float64, default::Float64)
    mp === nothing && return default
    lons, lats, vals = mp
    i = argmin(abs.(mod.(lons .- lon .+ 180.0, 360.0) .- 180.0))
    j = argmin(abs.(lats .- lat))
    v = Float64(vals[i, j])
    return isfinite(v) ? clamp(v, 0.0, 1.0) : default
end

function _dust_init_cells!(db::DustBox, gc, m2)
    lsm = m2["lsm"]
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    for j in 1:Int(gc.Nj), i in 1:Int(gc.Ni)
        la, lo = indyindxtolatlon(gc, Float64(j), Float64(i))
        (db.lat0 <= la <= db.lat1) || continue
        lom = mod(lo - db.lon0, 360.0) + db.lon0
        (db.lon0 <= lom <= db.lon1) || continue
        Float64(lsm[i, j, 1]) >= 0.5 || continue          # land only
        push!(db.ii, Int32(i)); push!(db.jj, Int32(j))
        push!(db.area, dlon * _DF_MDEG * cosd(la) * dlat * _DF_MDEG)
        push!(db.acc, 0.0)
        push!(db.clay, _map_sample(db.clay_map, lo, la, db.clay0))
        push!(db.bare, _map_sample(db.bare_map, lo, la, db.bare0))
    end
    db.initialized[] = true
    if isempty(db.ii)
        @warn "dust: no land cells inside the release box — source inactive"
    else
        @info "dust (Kok 2014) cells initialised (first step)" land_cells=length(db.ii) clay_mean=round(sum(db.clay)/length(db.clay); digits=3) bare_mean=round(sum(db.bare)/length(db.bare); digits=3) soil_moisture_source=(haskey(m2, "swvl1") ? "met swvl1" : "constant $(db.w_pct) %")
    end
    return nothing
end

"""
    apply_dust!(state, ctx, m2raw, db, d, indt, dt_s; nest=1, rng) → nspawned

One source step of the Kok et al. (2014) dust flux over the land cells of the box.
"""
function apply_dust!(state::State.ParcelState, ctx, m2::AbstractDict, db::DustBox,
                     d::DateTime, indt::Real, dt_s::Real; nest::Int = 1, rng = Random.default_rng())
    (db.t1 <= d <= db.t2) || return 0
    nestctx = ctx.nests[nest]
    gc  = nestctx.gridconfig
    db.initialized[] || _dust_init_cells!(db, gc, m2)
    us2 = m2["ustar"]
    it  = clamp(round(Int, indt), 1, size(us2, 3))
    sp  = haskey(m2, "sp") ? m2["sp"] : (haskey(m2, "lnsp") ? m2["lnsp"] : nothing)
    t2m = haskey(m2, "2t") ? m2["2t"] : nothing
    swv = haskey(m2, "swvl1") ? m2["swvl1"] : nothing
    dlon = _df_dlon(gc); dlat = _df_dlat(gc)
    nsp = 0; slot = 0
    @inbounds for c in eachindex(db.ii)
        i = Int(db.ii[c]); j = Int(db.jj[c])
        us = Float64(us2[i, j, it])
        isfinite(us) || continue
        ra = 1.2
        if sp !== nothing && t2m !== nothing
            p = Float64(sp[i, j, it]); p < 2000.0 && (p = exp(p))      # lnsp → Pa
            T = Float64(t2m[i, j, it])
            (isfinite(p) && isfinite(T) && T > 150.0) && (ra = p / (287.05 * T))
        end
        w = swv === nothing ? db.w_pct : 100.0 * max(Float64(swv[i, j, it]), 0.0) / 1.3
        F = dust_flux_kok(us, ra; clay = db.clay[c], bare = db.bare[c], w_pct = w, D = db.D,
                          cd0 = db.cd0, ce = db.ce, calpha = db.calpha)
        F > db.flux_max[] && (db.flux_max[] = F)
        db.acc[c] += F * db.area[c] * dt_s
        while db.acc[c] >= db.quantum
            slot = _df_next_free(state, slot)
            slot == 0 && return nsp
            la, lo = indyindxtolatlon(gc, Float64(j), Float64(i))
            jlo = lo + (rand(rng) - 0.5) * dlon
            jla = la + (rand(rng) - 0.5) * dlat
            indy, indx = latlontoindyindx(gc, jla, jlo)
            # keep the jittered point inside the (cropped) met grid; level 1 is the
            # surface after the reader's level shift.  indl 1.2 is inside the lowest
            # model layer, at 0.2·zlev[2] above the ground (layout contract) — about
            # 4 m on ERA5 L137, whose lowest full level (~10 m) sits at indl 1.5
            indx = clamp(Float64(indx), 1.0, Float64(gc.Ni) - 1.0e-3)
            indy = clamp(Float64(indy), 1.0, Float64(gc.Nj) - 1.0e-3)
            _write_boundary_parcel!(state, slot, nest, Float32(indx), Float32(indy),
                                    1.2f0, Float32(db.quantum), d)
            state.relbox[slot] = eltype(state.relbox)(db.relbox)
            db.acc[c] -= db.quantum
            db.emitted[] += db.quantum
            db.nspawn[] += 1
            nsp += 1
        end
    end
    return nsp
end
