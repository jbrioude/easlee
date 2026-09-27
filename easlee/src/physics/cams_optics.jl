# physics/cams_optics.jl
# ──────────────────────────────────────────────────────────────────────────────
# Optics of the CAMS aerosol classes (aermr01–20): the per-class table, the single-mode
# and lognormal mass-extinction coefficients and their κ-Köhler / Mie wet enhancement.
# TRANSPORT side (no chemistry dependency) since 2026-09-26, so the AOD optics can be
# validated against a CAMS analysis without the chemistry module
# (validation/check_aod_optics_cams.jl).  Consumers: chemistry/cams_background.jl
# (`cams_aod550_above`, the per-parcel overhead column) and chemistry/hetchem.jl
# (wet surface area).
# ──────────────────────────────────────────────────────────────────────────────

# ── CAMS aerosol fields (aermr01–20) → surface area + extinction ─────────────
#
# Per-type single-mode effective properties (IFS aerosol climatology class):
# dry effective radius r_eff [m], particle density ρ_p [kg/m³], κ-Köhler
# hygroscopicity, and dry 550-nm mass extinction coefficient MEC [m²/g].
# Surface-to-mass of a mode is 3/(ρ_p·r_eff); hygroscopic growth scales the
# surface (and, to leading order for r ≳ λ, the extinction) by gf².
# NOTE (paper): IFS sea-salt MMR is defined at 80 % RH — treated as dry here
# (few-% high bias on SS surface); refine with the IFS RH-banded optics later.
const CAMS_AERO = (
    # shortName    r_eff     ρ_p     κ     MEC550
    ("aermr01", 0.20e-6, 2160.0, 1.10,  3.0),   # sea salt 0.03–0.5 µm
    ("aermr02", 1.50e-6, 2160.0, 1.10,  0.4),   # sea salt 0.5–5 µm
    ("aermr03", 7.50e-6, 2160.0, 1.10,  0.10),  # sea salt 5–20 µm
    ("aermr04", 0.20e-6, 2610.0, 0.03,  2.5),   # dust 0.03–0.55 µm
    ("aermr05", 0.70e-6, 2610.0, 0.03,  1.0),   # dust 0.55–0.9 µm
    ("aermr06", 4.00e-6, 2610.0, 0.03,  0.25),  # dust 0.9–20 µm
    ("aermr07", 0.15e-6, 1400.0, 0.20,  4.5),   # OM hydrophilic
    ("aermr08", 0.15e-6, 1400.0, 0.00,  4.0),   # OM hydrophobic
    ("aermr09", 0.05e-6, 1000.0, 0.10,  8.0),   # BC hydrophilic
    ("aermr10", 0.05e-6, 1000.0, 0.00,  8.0),   # BC hydrophobic
    ("aermr11", 0.15e-6, 1770.0, 0.60,  4.5),   # sulfate
    ("aermr16", 0.15e-6, 1720.0, 0.67,  4.5),   # nitrate fine
    ("aermr17", 1.50e-6, 1720.0, 0.67,  0.6),   # nitrate coarse
    ("aermr18", 0.15e-6, 1770.0, 0.60,  4.5),   # ammonium
    ("aermr19", 0.15e-6, 1400.0, 0.20,  4.5),   # SOA biogenic
    ("aermr20", 0.15e-6, 1400.0, 0.20,  4.5),   # SOA anthropogenic
)
const CAMS_AERO_NAMES = String[a[1] for a in CAMS_AERO]

# ── RH-grown mass-extinction enhancement, per CAMS aerosol class ──────────────
# f_class(RH) = [Q_ext(x_wet, m_wet) / Q_ext(x_dry, m_dry)] · (d_wet/d_dry)²
#             = [Q_ext(x_wet, m_wet) / Q_ext(x_dry, m_dry)] · gf²
# the standard wet/dry MEC ratio: extinction cross-section per particle scales
# as Q_ext·d², so per-unit-DRY-mass extinction (MEC) scales by that ratio.
# Precomputed once on AEROSOL_RH_GRID (same table the size-bin optics use) so the
# column loop is a cheap RH-grid interpolation — no Mie call in the inner loop.
#
# Per-class refractive index for the wet volume-mixing: CAMS classes → the
# AEROSOL_REFRACTIVE symbols (sea salt / dust / OM / BC / sulfate /
# nitrate ≈ sulfate optically / ammonium ≈ sulfate).  κ and dry effective radius
# come straight from CAMS_AERO (cols 4 and 2), so the mapping only
# needs the refractive-index symbol.
const _CAMS_AERO_RI_SYM = (
    :seasalt, :seasalt, :seasalt,          # aermr01–03 sea salt
    :dust, :dust, :dust,                    # aermr04–06 dust
    :oc, :oc,                               # aermr07–08 organic matter
    :bc, :bc,                               # aermr09–10 black carbon
    :sulfate,                              # aermr11 sulfate
    :sulfate, :sulfate,                    # aermr16–17 nitrate (sulfate-like RI)
    :sulfate,                              # aermr18 ammonium
    :oc, :oc,                              # aermr19–20 SOA (organic RI)
)

# _CAMS_AERO_FRH[t] :: Vector{Float64} over AEROSOL_RH_GRID for CAMS class t.
const _CAMS_AERO_FRH = let
    lambda = 550.0e-9
    m_w    = 1.33 - 1.0e-9im
    tbl = Vector{Vector{Float64}}(undef, length(CAMS_AERO))
    for (t, a) in enumerate(CAMS_AERO)
        d_dry = 2.0 * Float64(a[2])                      # dry diameter = 2·r_eff [m]
        kap   = Float32(a[4])                            # κ-Köhler for this class
        sym   = _CAMS_AERO_RI_SYM[t]
        ri    = AEROSOL_REFRACTIVE[sym]
        x_dry = π * d_dry / lambda
        q_dry = mie_qext(x_dry, ri)
        f = Vector{Float64}(undef, length(AEROSOL_RH_GRID))
        for (jr, rh) in enumerate(AEROSOL_RH_GRID)
            gf = Float64(growth_factor(kap, rh))   # RH already capped at 0.98
            fw = 1.0 - 1.0 / gf^3                             # water volume fraction
            mi = ri * (1.0 - fw) + m_w * fw                  # volume-mixed wet RI
            x_wet = π * d_dry * gf / lambda
            q_wet = mie_qext(x_wet, mi)
            f[jr] = q_dry > 0.0 ? (q_wet / q_dry) * gf^2 : 1.0
        end
        tbl[t] = f
    end
    tbl
end

# Enhancement factor for CAMS class `t` at relative humidity `rh` (0–1), linearly
# interpolated on AEROSOL_RH_GRID.  RH capped at 0.98 (deliquescence hysteresis
# not modelled — the growth-factor table itself only spans up to 0.98).
@inline function _cams_aero_frh(t::Int, rh::Float64) :: Float64
    grid = AEROSOL_RH_GRID
    f    = _CAMS_AERO_FRH[t]
    nrh  = length(grid)
    r    = clamp(rh, 0.0, 0.98)
    jr = 1
    @inbounds while jr < nrh - 1 && Float64(grid[jr+1]) < r
        jr += 1
    end
    @inbounds begin
        g0 = Float64(grid[jr]); g1 = Float64(grid[jr+1])
        w  = clamp((r - g0) / max(g1 - g0, 1.0e-3), 0.0, 1.0)
        return (1.0 - w) * f[jr] + w * f[jr+1]
    end
end

# ── Per-class optics: single mode (the v1 table) or a lognormal size distribution ──────
# `cams_aero_optics(size_dist, sigma_g, nodes)` → (mec [m²/g] per class, f(RH) table per class)
#   "single"    : the CAMS_AERO table MEC and the single-radius Mie enhancement (_CAMS_AERO_FRH);
#   "lognormal" : each class is a lognormal number distribution of geometric width sigma_g
#                 whose EFFECTIVE radius is the table r_eff (count median r_g = r_eff·exp(−2.5 ln²σ));
#                 MEC = ∫Q_ext πr² n dln r / ∫(4/3)πr³ρ n dln r, integrated on `nodes` radii —
#                 Gauss–Hermite in ln r for nodes ≤ 9 (3 = the "three radii" version, Jerome
#                 25 Sept 2026: EASLEE/CAMS column ratio 1.1 for σ = 2; the converged integral
#                 gives 1.00), trapezoid over ±4 ln σ otherwise.  Wet: r → r·gf, RI volume-mixed.
const _GH_NODES = Dict(
    3 => ([-1.2247448713915890, 0.0, 1.2247448713915890], [0.2954089751509193, 1.1816359006036774, 0.2954089751509193]),
    5 => ([-2.0201828704560856, -0.9585724646138185, 0.0, 0.9585724646138185, 2.0201828704560856],
          [0.0199532420590459, 0.3936193231522412, 0.9453087204829419, 0.3936193231522412, 0.0199532420590459]),
    7 => ([-2.6519613568352334, -1.6735516287674714, -0.8162878828589647, 0.0, 0.8162878828589647, 1.6735516287674714, 2.6519613568352334],
          [0.0009717812450995, 0.0545155828191271, 0.4256072526101278, 0.8102646175568073, 0.4256072526101278, 0.0545155828191271, 0.0009717812450995]),
)
function _lognormal_optics(t::Int, sigma_g::Float64, nodes::Int)
    lambda = 550.0e-9; m_w = 1.33 - 1.0e-9im
    a = CAMS_AERO[t]; reff = Float64(a[2]); rho = Float64(a[3]); kap = Float32(a[4])
    ri = AEROSOL_REFRACTIVE[_CAMS_AERO_RI_SYM[t]]
    lns = log(sigma_g); rg = reff * exp(-2.5 * lns^2)
    if nodes <= 9
        xs, ws = get(_GH_NODES, nodes, _GH_NODES[3])
        rs = rg .* exp.(sqrt(2.0) * lns .* xs); ws = ws ./ sqrt(pi)
    else
        xs = range(-4.0, 4.0; length=nodes); rs = rg .* exp.(lns .* xs)
        ws = exp.(-xs .^ 2 ./ 2) ./ sqrt(2pi) .* step(xs); ws = ws ./ sum(ws)
    end
    massint = sum(w * (4 / 3) * pi * r^3 * rho for (w, r) in zip(ws, rs))
    ext(gf, mi) = sum(w * mie_qext(2pi * r * gf / lambda, mi) * pi * (r * gf)^2 for (w, r) in zip(ws, rs))
    ed = ext(1.0, ri)
    mec = ed / massint / 1000.0                                  # m²/kg → m²/g
    f = Vector{Float64}(undef, length(AEROSOL_RH_GRID))
    for (jr, rh) in enumerate(AEROSOL_RH_GRID)
        gf = Float64(growth_factor(kap, rh)); fw = 1.0 - 1.0 / gf^3
        mi = ri * (1.0 - fw) + m_w * fw
        f[jr] = ed > 0.0 ? ext(gf, mi) / ed : 1.0
    end
    return mec, f
end
function cams_aero_optics(size_dist::AbstractString = "single", sigma_g::Real = 2.0, nodes::Integer = 3)
    n = length(CAMS_AERO)
    if lowercase(size_dist) == "single"
        return [Float64(a[5]) for a in CAMS_AERO], [copy(_CAMS_AERO_FRH[t]) for t in 1:n]
    end
    lowercase(size_dist) == "lognormal" || error("chemistry.aerosols.size_dist must be single | lognormal")
    lg = [_lognormal_optics(t, Float64(sigma_g), Int(nodes)) for t in 1:n]
    return first.(lg), last.(lg)
end

@inline function _frh_table(f::Vector{Float64}, rh::Float64) :: Float64
    grid = AEROSOL_RH_GRID; nrh = length(grid); r = clamp(rh, 0.0, 0.98); jr = 1
    @inbounds while jr < nrh - 1 && Float64(grid[jr+1]) < r
        jr += 1
    end
    @inbounds begin
        g0 = Float64(grid[jr]); g1 = Float64(grid[jr+1])
        w  = clamp((r - g0) / max(g1 - g0, 1.0e-3), 0.0, 1.0)
        return (1.0 - w) * f[jr] + w * f[jr+1]
    end
end

