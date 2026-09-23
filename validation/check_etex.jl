#!/usr/bin/env julia
# check_etex.jl — self-check an EASLEE transport run against the ETEX-1
# observations, WITHOUT the validation task: reads the time-framed grid
# (grid_<name>.nc from configs/output_etex.toml), samples the surface layer
# bilinearly at every station for every 3-h ETEX window, writes a per-record
# CSV you can inspect/plot, and prints the ATMES-II / FLEXPART-11 statistics
# (same datem_stats as the validation task — scripts/datem.jl).
#
#   julia --project=easlee scripts/check_etex.jl [grid.nc] [obs_dir]
#
# defaults: results/etex/grid_output_etex.nc and the ETEX obs dir.
using Dates, Printf, Statistics, NetCDF

include(joinpath(@__DIR__, "datem.jl"))

const GRID = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "results", "etex", "grid_output_etex.nc")
const OBS  = length(ARGS) >= 2 ? ARGS[2] :
    "/media/jbrioude/952dd8b2-6e30-4f01-acdb-26dc78abd6c8/jbrioude/Desktop/julia/ETEX"
isfile(GRID) || error("grid file not found: $GRID — run configs/output_etex.toml first")

# ── model frames ──────────────────────────────────────────────────────────────
lon   = NetCDF.ncread(GRID, "lon")
lat   = NetCDF.ncread(GRID, "lat")
hgt   = NetCDF.ncread(GRID, "height")
tsec  = NetCDF.ncread(GRID, "time")
conc  = NetCDF.ncread(GRID, "conc")                     # (lon, lat, height, time) [kg/m³]
tunits = NetCDF.ncgetatt(GRID, "time", "units")         # "seconds since YYYY-MM-DDTHH:MM:SS"
t0 = DateTime(match(r"since (.+)$", tunits).captures[1])
avg_s = Float64(NetCDF.ncgetatt(GRID, "global", "average_s"))
frame_end = [t0 + Second(round(Int, s)) for s in tsec]
klay = size(conc, 3) >= 2 ? 2 : 1                       # first physical layer (0–100 m)
dlon = lon[2] - lon[1]; dlat = lat[2] - lat[1]
@printf("grid: %d×%d @ %.2f°, layer %.0f m, %d frames (%s → %s, avg %.0f s)\n",
        length(lon), length(lat), dlon, hgt[klay], length(tsec),
        frame_end[1], frame_end[end], avg_s)

# bilinear sample of frame `ti` at (la, lo) → pg/m³
function model_pg(lo::Float64, la::Float64, ti::Int)
    gx = (lo - lon[1]) / dlon + 1.0
    gy = (la - lat[1]) / dlat + 1.0
    (1.0 <= gx <= length(lon) && 1.0 <= gy <= length(lat)) || return NaN
    i0 = clamp(floor(Int, gx), 1, length(lon) - 1); wx = gx - i0
    j0 = clamp(floor(Int, gy), 1, length(lat) - 1); wy = gy - j0
    v = (1-wx)*(1-wy)*conc[i0,   j0,   klay, ti] + wx*(1-wy)*conc[i0+1, j0,   klay, ti] +
        (1-wx)*wy    *conc[i0,   j0+1, klay, ti] + wx*wy    *conc[i0+1, j0+1, klay, ti]
    return Float64(v) * 1.0e15                          # kg/m³ → pg/m³
end

# ── pair the observations with the matching frames ────────────────────────────
obs = read_meas(joinpath(OBS, "meas-t1.txt"))
frame_of = Dict(t => i for (i, t) in enumerate(frame_end))
O = Float64[]; M = Float64[]
csv = joinpath(dirname(GRID), "etex_station_check.csv")
open(csv, "w") do io
    println(io, "stn,lat,lon,window_start,window_end,obs_pg_m3,model_pg_m3")
    for ob in obs
        # the frame ending at the obs window's end holds the same 3-h mean
        ti = get(frame_of, ob.stop, 0)
        ti == 0 && continue                              # outside the run window
        Second(round(Int, avg_s)) == Second(ob.stop - ob.start) ||
            @warn "averaging window ≠ obs window — set output.average = obs duration" maxlog=1
        m = model_pg(ob.lon, ob.lat, ti)
        isfinite(m) || continue
        push!(O, ob.conc_pg); push!(M, m)
        @printf(io, "%d,%.3f,%.3f,%s,%s,%.2f,%.2f\n",
                ob.stn, ob.lat, ob.lon, ob.start, ob.stop, ob.conc_pg, m)
    end
end
isempty(O) && error("no obs window matched a model frame — check run window / interval")

# ── statistics ────────────────────────────────────────────────────────────────
s = datem_stats(O, M)
pos  = [(O[i] > 0 && M[i] > 0) for i in eachindex(O)]
fa5  = any(pos) ? count(i -> pos[i] && 0.2 <= M[i]/O[i] <= 5.0, eachindex(O)) / count(pos) : 0.0
println("\nETEX-1 vs EASLEE (pooled, $(s.n) station-windows, incl. zeros)")
@printf("  obs mean %.1f pg/m³   model mean %.1f pg/m³\n", s.obar, s.mbar)
@printf("  R=%.3f  FB=%+.3f  NMSE=%.2f  FAC2=%.2f  FA5=%.2f  FMS=%.1f%%  KSP=%.1f  RANK=%.2f\n",
        s.r, s.fb, s.nmse, s.fac2, fa5, s.fms, s.ksp, s.rank)
println("  reference — FLEXPART 11 (ETEX-1): PCC=0.84 FA2=0.74 FA5=0.87 FMS=0.68 NMSE=0.43")
println("  EASLEE validation task (ETEX_RESULTS.md): FMS 0.81, rank-corr 0.79")
println("\nper-record pairs → $csv")

# ── receptor estimator (gridless KDE at the station coordinates) ─────────────
RCSV = length(ARGS) >= 3 ? ARGS[3] :
    joinpath(dirname(GRID), "receptors_" * splitext(basename(GRID))[1][6:end] * ".csv")
if isfile(RCSV)
    rc = Dict{Tuple{Int,DateTime},Float64}()
    rres = Float64[]
    for ln in eachline(RCSV)
        (startswith(ln, "#") || startswith(ln, "time_s")) && continue
        f = split(ln, ",")
        length(f) >= 5 || continue
        t = t0 + Second(parse(Int, f[1]))
        stn = tryparse(Int, f[2]); stn === nothing && continue
        rc[(stn, t)] = parse(Float64, f[3]) * 1.0e15
        r = tryparse(Float64, f[4]); r === nothing || push!(rres, r)
    end
    Or = Float64[]; Mr = Float64[]
    for ob in obs
        haskey(rc, (ob.stn, ob.stop)) || continue
        push!(Or, ob.conc_pg); push!(Mr, rc[(ob.stn, ob.stop)])
    end
    if !isempty(Or)
        sr = datem_stats(Or, Mr)
        posr = [(Or[i] > 0 && Mr[i] > 0) for i in eachindex(Or)]
        fa5r = any(posr) ? count(i -> posr[i] && 0.2 <= Mr[i]/Or[i] <= 5.0, eachindex(Or)) / count(posr) : 0.0
        println("\nETEX-1 vs EASLEE — RECEPTORS (gridless, $(sr.n) station-windows)")
        @printf("  obs mean %.1f pg/m³   model mean %.1f pg/m³\n", sr.obar, sr.mbar)
        @printf("  R=%.3f  FB=%+.3f  NMSE=%.2f  FAC2=%.2f  FA5=%.2f  FMS=%.1f%%  KSP=%.1f  RANK=%.2f\n",
                sr.r, sr.fb, sr.nmse, sr.fac2, fa5r, sr.fms, sr.ksp, sr.rank)
        rrf = filter(isfinite, rres)
        isempty(rrf) || @printf("  native resolution at the stations: median %.0f km (p10 %.0f, p90 %.0f)\n",
                median(rrf)/1e3, quantile(rrf, 0.1)/1e3, quantile(rrf, 0.9)/1e3)
    end
else
    println("(no receptor CSV at $RCSV — run with [output.receptors])")
end
