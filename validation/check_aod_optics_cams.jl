# check_aod_optics_cams.jl — validate EASLEE's aerosol optics against a CAMS analysis,
# WITHOUT transport or chemistry (Jerome, 26 Sept 2026): take the CAMS aerosol mixing ratios
# (aermr01–20) of one analysis, compute the 550-nm column AOD with the model's per-class
# optics (Mie, κ-Köhler growth at the ERA5 relative humidity of the same hour), and compare
# with the 532-nm extinction column (aerext532) the CAMS system itself carries in the same file.
#
#   julia --project=easlee validation/check_aod_optics_cams.jl \
#        --cams /data/chemistry/2025/All/cams_20250920_1200_3.grib \
#        --era5 /data/ERA5/All/20250920/IFS_EA_MET_20250920_600_9.grb \
#        [--size lognormal|single] [--sigma 2.0] [--nodes 3] [--dry] [--window lon0,lon1,lat0,lat1]
#
# Both files must be the same valid time (here 15 UTC: CAMS 12Z + 3 h, ERA5 06Z + 9 h).  The
# ERA5 file provides t, q (137 model levels) and lnsp; the CAMS file the 25 pressure levels.
# PASS: mean EASLEE/CAMS ratio in [0.7, 1.4] and spatial correlation > 0.7 over the window
# (default: the paper's Africa / SW Indian Ocean window −13…73 °E, −38…13 °N; 26 Sept 2026,
# lognormal σ 2, 3 radii: ratio 1.11, r 0.75; dry single mode: 0.55, r 0.96).
global gpu = false; global chem = false
using Dates, Printf, Statistics, GRIB
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: CAMS_AERO, cams_aero_optics, _frh_table     # physics/cams_optics.jl (transport side)

function argval(flag, default)
    i = findfirst(==(flag), ARGS); i === nothing || i == length(ARGS) ? default : ARGS[i+1]
end
CAMS = argval("--cams", "/data/chemistry/2025/All/cams_20250920_1200_3.grib")
ERA5 = argval("--era5", "/data/ERA5/All/20250920/IFS_EA_MET_20250920_600_9.grb")
SIZE = argval("--size", "lognormal"); SIG = parse(Float64, argval("--sigma", "2.0")); NODES = parse(Int, argval("--nodes", "3"))
DRY  = "--dry" in ARGS
WIN  = parse.(Float64, split(argval("--window", "-13,73,-38,13"), ","))
const g = 9.80665

# ── CAMS: aermr fields on pressure levels + aerext532 ───────────────────────────────────
names = String[a[1] for a in CAMS_AERO]
fields = Dict{String,Dict{Int,Matrix{Float64}}}()   # shortName → level → values (lon × lat)
const CGRID = Ref{Any}(nothing)                                    # (lons, lats) — set inside the closure
GribFile(CAMS) do f
    for msg in f
        m = Message(msg); sn = String(m["shortName"])
        (sn in names || sn == "aerext532") || continue
        lv = Int(m["level"])
        CGRID[] === nothing && (CGRID[] = (Float64.(m["distinctLongitudes"]), Float64.(m["distinctLatitudes"])))
        get!(fields, sn, Dict{Int,Matrix{Float64}}())[lv] = Matrix{Float64}(m["values"])
    end
end
clon, clat = CGRID[]
haskey(fields, "aerext532") || error("no aerext532 in $CAMS")
plevs = sort(collect(keys(fields["aerext532"])))                 # hPa, ascending = top first
carried = [n for n in names if haskey(fields, n)]
Ni, Nj = size(fields["aerext532"][plevs[1]])
@printf("CAMS %s: %d×%d, %d levels (%d–%d hPa), classes %s\n", basename(CAMS), Ni, Nj, length(plevs), plevs[1], plevs[end], join(carried, " "))

# ── ERA5: T, q on model levels + lnsp → RH on the CAMS levels ──────────────────────────
T = Dict{Int,Matrix{Float64}}(); Q = Dict{Int,Matrix{Float64}}()
const EGRID = Ref{Any}(nothing); const LNSP = Ref{Any}(nothing)
GribFile(ERA5) do f
    for msg in f
        m = Message(msg); sn = String(m["shortName"])
        if sn == "t" || sn == "q"
            (sn == "t" ? T : Q)[Int(m["level"])] = Matrix{Float64}(m["values"])
            EGRID[] === nothing && (EGRID[] = (Float64.(m["pv"]), Float64.(m["distinctLongitudes"]), Float64.(m["distinctLatitudes"])))
        elseif sn == "lnsp"
            LNSP[] = Matrix{Float64}(m["values"])
        end
    end
end
(LNSP[] !== nothing && EGRID[] !== nothing) || error("ERA5 file needs t, q and lnsp")
pv, elon, elat = EGRID[]; lnsp = LNSP[]
nlev = length(T); A = pv[1:nlev+1]; B = pv[nlev+2:2nlev+2]
@printf("ERA5 %s: %d×%d, %d levels\n", basename(ERA5), size(lnsp)..., nlev)
# nearest ERA5 index for each CAMS point (both grids are regular; longitudes may differ in origin)
ei = [argmin(abs.(mod.(elon .- lo .+ 180, 360) .- 180)) for lo in clon]
ej = [argmin(abs.(elat .- la)) for la in clat]
lvl = sort(collect(keys(T)))
function rh_column!(out, i, j)
    ie = ei[i]; je = ej[j]; ps = exp(lnsp[ie, je])
    lp = Vector{Float64}(undef, nlev); rh = Vector{Float64}(undef, nlev)
    for (k, l) in enumerate(lvl)
        p = 0.5 * ((A[k] + B[k] * ps) + (A[k+1] + B[k+1] * ps)); lp[k] = log(p)
        t = T[l][ie, je]; q = max(Q[l][ie, je], 0.0)
        e = q * p / (0.622 + 0.378 * q); es = 611.2 * exp(17.67 * (t - 273.15) / (t - 29.65))
        rh[k] = clamp(e / max(es, 1e-3), 0.0, 0.98)
    end
    for (k, ph) in enumerate(plevs)
        x = log(ph * 100.0)
        kk = clamp(searchsortedlast(lp, x), 1, nlev - 1)
        w = clamp((x - lp[kk]) / (lp[kk+1] - lp[kk]), 0.0, 1.0)
        out[k] = (1 - w) * rh[kk] + w * rh[kk+1]
    end
    return out
end

# ── optics + column integrals ─────────────────────────────────────────────────────────
mec_all, ftab_all = cams_aero_optics(SIZE, SIG, NODES)
cls = [findfirst(==(n), names) for n in carried]
inwin = [(WIN[1] <= mod(lo + 180, 360) - 180 <= WIN[2]) for lo in clon], [(WIN[3] <= la <= WIN[4]) for la in clat]
E = fill(NaN, Ni, Nj); R = fill(NaN, Ni, Nj); rhp = zeros(length(plevs))
Rd = 287.05; T0 = 260.0
edges = vcat(0.5 * plevs[1], 0.5 .* (plevs[1:end-1] .+ plevs[2:end]), 1013.0) .* 100.0   # Pa, top → surface
dz = Rd * T0 / g .* log.(edges[2:end] ./ max.(edges[1:end-1], 1.0))
for j in 1:Nj, i in 1:Ni
    (inwin[1][i] && inwin[2][j]) || continue
    DRY || rh_column!(rhp, i, j)
    aod = 0.0; eprev = 0.0; pprev = plevs[1]
    for (k, ph) in enumerate(plevs)
        e = 0.0; rh = DRY ? 0.0 : rhp[k]
        for (m, n) in enumerate(carried)
            q = fields[n][ph][i, j]; q > 0 || continue
            e += q * 1e3 * mec_all[cls[m]] * (DRY ? 1.0 : _frh_table(ftab_all[cls[m]], rh))
        end
        k > 1 && (aod += 0.5 * (e + eprev) * (ph - pprev) * 100.0 / g)
        eprev = e; pprev = ph
    end
    E[i, j] = aod
    R[i, j] = sum(fields["aerext532"][ph][i, j] * dz[k] for (k, ph) in enumerate(plevs))
end
m = isfinite.(E) .& isfinite.(R)
ratio = mean(E[m]) / mean(R[m]); r = cor(vec(E[m]), vec(R[m]))
@printf("optics %s (σ %.1f, %d nodes)%s — window %s: EASLEE AOD550 mean %.3f vs CAMS AOD532 %.3f, ratio %.2f, corr %.2f, n %d\n",
        SIZE, SIG, NODES, DRY ? ", DRY" : ", RH from ERA5", join(WIN, ","), mean(E[m]), mean(R[m]), ratio, r, count(m))
ok = 0.7 <= ratio <= 1.4 && r > 0.7
println(ok ? "PASS" : "FAIL (need 0.7 ≤ ratio ≤ 1.4 and corr > 0.7)")
exit(ok ? 0 : 1)
