#!/usr/bin/env julia
# check_captex.jl — score a CAPTEX A+B run against the DATEM observations,
# for BOTH output estimators side by side:
#   • gridded : bilinear sample of the 3-h grid frames at each station
#   • receptor: the gridless KDE time series dumped AT the station coordinates
#               (results/captex/receptors_output_captex_ab.csv)
# Prints ATMES-II / FLEXPART-11 statistics for each (shared scripts/datem.jl).
#
#   julia --project=easlee scripts/check_captex.jl [grid.nc] [obs_dir] [receptors.csv]
using Dates, Printf, Statistics, NetCDF

include(joinpath(@__DIR__, "datem.jl"))

const GRID = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(@__DIR__, "..", "results", "captex", "grid_output_captex_ab.nc")
const OBS  = length(ARGS) >= 2 ? ARGS[2] :
    get(ENV, "EASLEE_CAPTEX_DIR", joinpath(@__DIR__, "..", "data", "CAPTEX"))
const RCSV = length(ARGS) >= 3 ? ARGS[3] :
    joinpath(dirname(GRID), "receptors_" * splitext(basename(GRID))[1][6:end] * ".csv")
isfile(GRID) || error("grid file not found: $GRID — run configs/output_captex_ab.toml first")

obs = read_meas(joinpath(OBS, "meas-t1.txt"))

function report(tag, O, M)
    s = datem_stats(O, M)
    pos = [(O[i] > 0 && M[i] > 0) for i in eachindex(O)]
    fa5 = any(pos) ? count(i -> pos[i] && 0.2 <= M[i]/O[i] <= 5.0, eachindex(O)) / count(pos) : 0.0
    @printf("\n%s  (%d station-windows, incl. zeros)\n", tag, s.n)
    @printf("  obs mean %.1f pg/m³   model mean %.1f pg/m³\n", s.obar, s.mbar)
    @printf("  R=%.3f  FB=%+.3f  NMSE=%.2f  FAC2=%.2f  FA5=%.2f  FMS=%.1f%%  KSP=%.1f  RANK=%.2f\n",
            s.r, s.fb, s.nmse, s.fac2, fa5, s.fms, s.ksp, s.rank)
    return s
end

# ── gridded estimator ─────────────────────────────────────────────────────────
lon   = NetCDF.ncread(GRID, "lon");  lat = NetCDF.ncread(GRID, "lat")
hgt   = NetCDF.ncread(GRID, "height")
tsec  = NetCDF.ncread(GRID, "time")
conc  = NetCDF.ncread(GRID, "conc")                     # (lon, lat, height, time) [kg/m³]
t0 = DateTime(match(r"since (.+)$", NetCDF.ncgetatt(GRID, "time", "units")).captures[1])
avg_s = Float64(NetCDF.ncgetatt(GRID, "global", "average_s"))
frame_end = [t0 + Second(round(Int, s)) for s in tsec]
klay = size(conc, 3) >= 2 ? 2 : 1
dlon = lon[2] - lon[1]; dlat = lat[2] - lat[1]
@printf("grid: %d×%d @ %.2f°, layer %.0f m, %d frames (%s → %s)\n",
        length(lon), length(lat), dlon, hgt[klay], length(tsec),
        frame_end[1], frame_end[end])

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

frame_of = Dict(t => i for (i, t) in enumerate(frame_end))
# CAPTEX obs windows are 3 h or 6 h; a 6-h window = mean of its two 3-h frames
function grid_pg(ob)
    dur = Second(ob.stop - ob.start)
    if haskey(frame_of, ob.stop) && dur == Second(round(Int, avg_s))
        return model_pg(ob.lon, ob.lat, frame_of[ob.stop])
    elseif dur == Second(2 * round(Int, avg_s)) &&
           haskey(frame_of, ob.stop) && haskey(frame_of, ob.stop - Second(round(Int, avg_s)))
        a = model_pg(ob.lon, ob.lat, frame_of[ob.stop - Second(round(Int, avg_s))])
        b = model_pg(ob.lon, ob.lat, frame_of[ob.stop])
        return 0.5 * (a + b)
    end
    return NaN
end

Og = Float64[]; Mg = Float64[]
csv = joinpath(dirname(GRID), "captex_station_check.csv")
open(csv, "w") do io
    println(io, "stn,lat,lon,window_start,window_end,obs_pg_m3,grid_pg_m3")
    for ob in obs
        m = grid_pg(ob)
        isfinite(m) || continue
        push!(Og, ob.conc_pg); push!(Mg, m)
        @printf(io, "%d,%.3f,%.3f,%s,%s,%.2f,%.2f\n",
                ob.stn, ob.lat, ob.lon, ob.start, ob.stop, ob.conc_pg, m)
    end
end
isempty(Og) && error("no obs window matched a model frame — check run window / interval")
report("CAPTEX R1 — GRIDDED (grid frames, bilinear at stations)", Og, Mg)

# ── receptor estimator ────────────────────────────────────────────────────────
if isfile(RCSV)
    # rows: time_s,station,conc_kgm3,resolution_m,nsamp  (time = window END)
    rc = Dict{Tuple{Int,DateTime},Float64}()
    rres = Float64[]
    for ln in eachline(RCSV)
        (startswith(ln, "#") || startswith(ln, "time_s")) && continue
        f = split(ln, ",")
        length(f) >= 5 || continue
        t = t0 + Second(parse(Int, f[1]))
        stn = tryparse(Int, f[2]); stn === nothing && continue
        rc[(stn, t)] = parse(Float64, f[3]) * 1.0e15    # kg/m³ → pg/m³
        r = tryparse(Float64, f[4]); r === nothing || push!(rres, r)
    end
    Or = Float64[]; Mr = Float64[]
    for ob in obs
        dur = Second(ob.stop - ob.start)
        m = if dur == Second(round(Int, avg_s)) && haskey(rc, (ob.stn, ob.stop))
            rc[(ob.stn, ob.stop)]
        elseif dur == Second(2 * round(Int, avg_s)) &&
               haskey(rc, (ob.stn, ob.stop)) &&
               haskey(rc, (ob.stn, ob.stop - Second(round(Int, avg_s))))
            0.5 * (rc[(ob.stn, ob.stop)] + rc[(ob.stn, ob.stop - Second(round(Int, avg_s)))])
        else
            NaN
        end
        isfinite(m) || continue
        push!(Or, ob.conc_pg); push!(Mr, m)
    end
    if !isempty(Or)
        report("CAPTEX R1 — RECEPTORS (gridless KDE at the station coordinates)", Or, Mr)
        rrf = filter(isfinite, rres)
        isempty(rrf) || @printf("  native resolution at the stations: median %.0f km (p10 %.0f, p90 %.0f)\n",
                median(rrf)/1e3, quantile(rrf, 0.1)/1e3, quantile(rrf, 0.9)/1e3)
    else
        println("\n(receptor CSV present but no windows matched)")
    end
else
    println("\n(no receptor CSV at $RCSV — run with [output.receptors])")
end
println("\nreference — validation task (CAPTEX_RESULTS.md): RANK ≈ 2.50")
println("per-record pairs → $csv")
