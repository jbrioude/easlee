# check_wellmixed_met.jl — well-mixed criterion on REAL meteorological data.
#
#     julia -t 16 --project=easlee validation/check_wellmixed_met.jl configs/<run>.toml \
#           HOURS NPARCELS [lon0 lon1 lat0 lat1]
#
# Parcels of equal mass are seeded uniformly in AIR MASS (cell weight ρ·area·Δz) from the
# ground to 15 km over a box of the met domain (the whole domain when no box is given,
# minus 2 cells), then advanced HOURS with the configuration's transport — advection,
# vertical coordinate velocity, turbulence — convection, deposition and releases off.
# A transport consistent with the met fields keeps the mixing ratio χ = parcel mass /
# air mass uniform.  χ is measured at the end, per height band and 8×8-cell block, in
# an interior that air from outside the seeded volume cannot have reached (margin =
# 1.3 × the largest wind below 8 km × HOURS horizontally; bands below 6 km).
# Reported per band: mean χ (must be 1) and the block-to-block spread over the Poisson
# spread (≈ 1 if uniform).  Initial χ is reported too (checks the seeding).
# WM_TURB=none: advection only (winds, density, vertical coordinate velocity), no turbulence.
# Written 2026-10-13 (Jerome: run the well-mixed criterion on every met dataset).

using Dates, Printf, Random, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE, .EASLEE.State
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")
global gpu = false
const OP = EASLEE.Operators

cfg = ARGS[1]; HOURS = parse(Float64, ARGS[2]); NP = parse(Int, ARGS[3])
box = length(ARGS) >= 7 ? parse.(Float64, ARGS[4:7]) : nothing
spec = read_runspec(cfg)
spec.physics.convection.scheme = "off"
spec.physics.deposition.wet = false; spec.physics.deposition.dry = false
for s in spec.species; s.wet_deposition = false; s.dry_deposition = false; end
for r in spec.release; r.particles = 1; r.time = [spec.stop + Hour(1), spec.stop + Hour(1)]; end
spec.stop = max(spec.stop, spec.start + Minute(round(Int, 60HOURS)) + Hour(1))
rt, obj = build_run(spec; prefetch=false, extra_slots=NP + 10)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
read_met_step!(rt, rc.startdate)
if get(ENV, "WM_TURB", "") == "none"                 # advection only (winds, ρ, coordinate velocity)
    for n in 1:length(ctx.nests); ctx.nests[n].nestrunconfig.turbulence = EASLEE.Types.NoTurb(); end
end
nc_ = ctx.nests[1]; gc = nc_.gridconfig; m3 = nc_.metfields3D
Ni, Nj = Int(gc.Ni), Int(gc.Nj); Nk = size(m3.zlev, 3)
println("met ", nc_.metconfig.model, "  grid ", nameof(typeof(gc)), " $Ni×$Nj×$Nk  turbulence ", typeof(nc_.nestrunconfig.turbulence))

# column area [m²] from the lon/lat of the cell corners (any projection)
function colarea(i, j)
    ll(x, y) = OP.indyindxtolatlon(gc, Float32(y), Float32(x))
    a = ll(i - 0.5, j - 0.5); b = ll(i + 0.5, j - 0.5); c = ll(i - 0.5, j + 0.5)
    dx = 6.371e6 * deg2rad(mod(b[2] - a[2] + 180, 360) - 180) * cosd(a[1]); dy0 = 6.371e6 * deg2rad(b[1] - a[1])
    ex = 6.371e6 * deg2rad(mod(c[2] - a[2] + 180, 360) - 180) * cosd(a[1]); ey = 6.371e6 * deg2rad(c[1] - a[1])
    return abs(dx * ey - dy0 * ex)
end
# box in index space
if box === nothing
    i0, i1, j0, j1 = 3, Ni - 3, 3, Nj - 3
else
    jy0, ix0 = OP.latlontoindyindx(gc, box[3], box[1]); jy1, ix1 = OP.latlontoindyindx(gc, box[4], box[2])
    i0, i1 = round(Int, min(ix0, ix1)), round(Int, max(ix0, ix1)); j0, j1 = round(Int, min(jy0, jy1)), round(Int, max(jy0, jy1))
    i0 = max(i0, 3); j0 = max(j0, 3); i1 = min(i1, Ni - 3); j1 = min(j1, Nj - 3)
end
area = [colarea(i, j) for i in 1:Ni, j in 1:Nj]
ZFILL = 15000.0
# layer air masses (slot k → k+1) at time index tt, cut at ZFILL.  Stacked layout
# (the transport code's): zlev on half levels, the layer zlev[k]..zlev[k+1] holds
# the full level of slot k+1, whose density is rho[k+1].
zl(i, j, k, tt) = OP.interpolation4Dlevsingle(m3.zlev, Float32(i), Float32(j), Float32(k), tt)
rh(i, j, k, tt) = OP.interpolation4Dlevsingle(m3.rho,  Float32(i), Float32(j), Float32(k), tt)
function layer_masses(tt)
    W = Float64[]; C = NTuple{3,Int}[]
    for j in j0:j1, i in i0:i1
        for k in 1:Nk-1                      # separate loop: `break` must end this column only
            za = zl(i, j, k, tt); zb = zl(i, j, k + 1, tt); za >= ZFILL && break
            zb = min(zb, ZFILL)
            push!(W, rh(i, j, k + 1, tt) * (zb - za) * area[i, j]); push!(C, (i, j, k))
        end
    end
    return W, C
end
t0 = EASLEE.MetAdapter.met_time_index_array(rc.startdate, ctx.nests)[1]
W, C = layer_masses(t0); cw = cumsum(W); M0 = cw[end]; mp = M0 / NP
rng = MersenneTwister(11)
fill!(state.status, 0x00)
for p in 1:NP
    q = searchsortedfirst(cw, rand(rng) * M0); i, j, k = C[q]
    za = zl(i, j, k, t0); zb = zl(i, j, k + 1, t0); zt = min(zb, ZFILL)
    f = (rand(rng) * (zt - za)) / max(zb - za, 1e-3)              # uniform in z inside the layer
    state.status[p] = 0x01; state.nest[p] = 1; state.is_new[p] = false; state.relbox[p] = 1
    state.indx[p] = Float32(i - 0.5 + rand(rng)); state.indy[p] = Float32(j - 0.5 + rand(rng))
    state.indl[p] = Float32(k + f); state.mass[p] = 1.0f0; state.Dtime[p] = rc.startdate
end
State.sync_dtime_ns!(state)

# evaluation interior and bands
umax = 0.0
for j in j0:4:j1, i in i0:4:i1
    for k in 1:Nk
        zl(i, j, k, t0) > 8000 && break
        global umax = max(umax, hypot(m3.Ux[i, j, k, 1], m3.Vy[i, j, k, 1]))
    end
end
spacing = sqrt(area[(i0 + i1) ÷ 2, (j0 + j1) ÷ 2])
mg = ceil(Int, 1.3 * umax * HOURS * 3600 / spacing)
e0, e1, f0, f1 = i0 + mg, i1 - mg, j0 + mg, j1 - mg
(e1 - e0 >= 8 && f1 - f0 >= 8) || error("box too small for the margin ($mg cells at umax $(round(umax, digits=1)) m/s)")
bands = [0, 100, 300, 600, 1000, 2000, 3500, 6000.0]
println(@sprintf("seeded %d parcels (%.3g kg of air), box i %d-%d j %d-%d; interior i %d-%d j %d-%d (margin %d cells, umax %.1f m/s)",
                 NP, M0, i0, i1, j0, j1, e0, e1, f0, f1, mg, umax))

function chi_table(tt, label)
    bi(x) = (x - e0) ÷ 8; nbx = (e1 - e0) ÷ 8; nby = (f1 - f0) ÷ 8
    nb = length(bands) - 1
    cnt = zeros(nbx, nby, nb); air = zeros(nbx, nby, nb)
    for p in eachindex(state.status)
        state.status[p] == 0x01 || continue
        x = state.indx[p]; y = state.indy[p]
        # cell i spans x ∈ [i − ½, i + ½), as in the seeding (2026-10-22: the
        # blocks were taken from floor(x), half a cell off the air-mass columns)
        (e0 - 0.5 <= x < e0 + 8nbx - 0.5 && f0 - 0.5 <= y < f0 + 8nby - 0.5) || continue
        z = OP.interpolation4Dlevsingle(m3.zlev, x, y, state.indl[p], tt)
        b = searchsortedlast(bands, z); (1 <= b <= nb) || continue
        cnt[floor(Int, x - e0 + 0.5) ÷ 8 + 1, floor(Int, y - f0 + 0.5) ÷ 8 + 1, b] += 1
    end
    for bx in 1:nbx, by in 1:nby, ii in 0:7, jj in 0:7
        i = e0 + 8(bx - 1) + ii; j = f0 + 8(by - 1) + jj
        for k in 1:Nk-1
            za = zl(i, j, k, tt); zb = zl(i, j, k + 1, tt); za >= bands[end] && break
            r = rh(i, j, k + 1, tt)     # stacked layout: the layer zlev[k]..zlev[k+1] holds full level k+1
            for b in 1:nb
                lo = max(za, bands[b]); hi = min(zb, bands[b+1])
                hi > lo && (air[bx, by, b] += r * (hi - lo) * area[i, j])
            end
        end
    end
    println("  $label")
    println("    band [m AGL]      mean χ   spread/Poisson   parcels/block")
    for b in 1:nb
        chi = vec(cnt[:, :, b] .* mp ./ air[:, :, b]); n = vec(cnt[:, :, b])
        ok = n .> 0
        mchi = sum(n[ok]) * mp / sum(vec(air[:, :, b])[ok])
        pois = sqrt(mean(1 ./ n[ok]))
        @printf("    %5.0f-%-5.0f       %6.3f   %6.2f           %6.0f\n", bands[b], bands[b+1], mchi, std(chi[ok]) / mchi / pois, mean(n))
    end
    # the two lowest bands split by the relief of the block: range of the surface
    # pressure (slot 1) over its 8×8 cells, 6 hPa ≈ 50 m, 60 hPa ≈ 500 m
    ps = @view m3.puv[:, :, 1, clamp(round(Int, tt), 1, size(m3.puv, 4))]
    if minimum(ps[e0:e0+8nbx-1, f0:f0+8nby-1]) > 1.0e4
        rel = [let blk = ps[e0+8(bx-1):e0+8bx-1, f0+8(by-1):f0+8by-1]; maximum(blk) - minimum(blk); end
               for bx in 1:nbx, by in 1:nby]
        for (cl, lo, hi) in (("flat  (<50 m)", 0.0, 600.0), ("hilly (50-500 m)", 600.0, 6000.0), ("steep (>500 m)", 6000.0, Inf))
            sel = vec((rel .>= lo) .& (rel .< hi)); count(sel) == 0 && continue
            s = join([@sprintf("%5.0f-%-5.0f χ %6.3f", bands[b], bands[b+1],
                               sum(vec(cnt[:, :, b])[sel]) * mp / sum(vec(air[:, :, b])[sel])) for b in 1:2], "   ")
            @printf("    relief %-17s %4d blocks   %s\n", cl, count(sel), s)
        end
    end
end
if get(ENV, "WM_DEBUG", "0") == "1"
    println("DEBUG t0=$t0 active ", count(==(0x01), state.status))
    for p in (1, 2, 3, 1000)
        x = state.indx[p]; y = state.indy[p]
        println("DEBUG p$p x=$x y=$y l=$(state.indl[p]) z=", OP.interpolation4Dlevsingle(m3.zlev, x, y, state.indl[p], t0), " e0=$e0 e1=$e1 f0=$f0 f1=$f1")
    end
end
chi_table(t0, "t = 0 (seeding)")
d = rc.startdate; tend = rc.startdate + Minute(round(Int, 60HOURS))
while d < tend
    global d
    read_met_step!(rt, d); step_runtime!(rt, d; do_output=false); State.sync_dtime_ns!(state)
    d += rc.timestep
end
read_met_step!(rt, d)
chi_table(EASLEE.MetAdapter.met_time_index_array(d, ctx.nests)[1], "t = $(HOURS) h")
# WM_DEBUG=2: where the near-ground parcels of the steep blocks sit (level index
# fraction above the ground slot, and their level rate)
if get(ENV, "WM_DEBUG", "0") == "2"
    ps = @view m3.puv[:, :, 1, 1]
    nbx = (e1 - e0) ÷ 8; nby = (f1 - f0) ÷ 8
    steep = falses(nbx, nby)
    for bx in 1:nbx, by in 1:nby
        blk = ps[e0+8(bx-1):e0+8bx-1, f0+8(by-1):f0+8by-1]; steep[bx, by] = maximum(blk) - minimum(blk) > 6000
    end
    edges = (1.0, 1.02, 1.1, 1.3, 1.5, 2.0, 3.0, 5.0)
    h = zeros(Int, length(edges) - 1); hf = zeros(Int, length(edges) - 1)
    for p in eachindex(state.status)
        state.status[p] == 0x01 || continue
        x = state.indx[p]; y = state.indy[p]
        (e0 - 0.5 <= x < e0 + 8nbx - 0.5 && f0 - 0.5 <= y < f0 + 8nby - 0.5) || continue
        b = searchsortedlast(collect(edges), state.indl[p]); (1 <= b < length(edges)) || continue
        steep[floor(Int, x - e0 + 0.5) ÷ 8 + 1, floor(Int, y - f0 + 0.5) ÷ 8 + 1] ? (h[b] += 1) : (hf[b] += 1)
    end
    println("  level-index histogram (indl bins ", edges, ")")
    println("    steep blocks: ", h); println("    other blocks: ", hf)
end
