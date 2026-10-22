# check_output_mapping.jl — does a parcel's mass land in the output cell that
# contains it?
#
#     julia -t 4 --project=easlee validation/check_output_mapping.jl configs/<run>.toml [N]
#
# One parcel at a time at a random position of the met domain, gridded by the
# CPU (deposit_to_grid!) and the KernelAbstractions (deposit_to_grid_ka!) paths
# with the point estimator ([output].kernel = "none" is forced).  The cell that
# receives the mass must be the one containing the parcel's lon/lat, cell ii
# spanning [lon0 + (ii-1)·dx, lon0 + ii·dx) as the NetCDF writer labels it.
# Exact check — every parcel, no tolerance.  Written 2026-10-13 after the
# output offset (deposits shifted by (met/output spacing − 1) cells, i.e. 1.5
# cells on the ETEX grid) and the Lambert output grid (no deposit at all).

using Dates, Printf, Random
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE, .EASLEE.State
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")
global gpu = false

cfg = ARGS[1]; NP = length(ARGS) > 1 ? parse(Int, ARGS[2]) : 300
spec = read_runspec(cfg)
spec.output.kernel = "none"
rt, obj = build_run(spec; prefetch=false, extra_slots=10)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
read_met_step!(rt, rc.startdate)
rc.currentdate = rc.startdate
gc = ctx.nests[1].gridconfig; oc = ctx.output; g = spec.output.grid
Ni, Nj = Int(gc.Ni), Int(gc.Nj)
println("grid ", nameof(typeof(gc)), "  crs ", g.crs, "  output lon0 $(g.lon0) lat0 $(g.lat0) dx $(g.dx) dy $(g.dy) nx $(g.nx) ny $(g.ny)")

fill!(state.status, 0x00)
slot = 1
rng = MersenneTwister(7)
bad = Dict("cpu" => 0, "ka" => 0); n = 0; outside = 0
for q in 1:NP
    global n, outside
    x = Float32(2 + rand(rng) * (Ni - 4)); y = Float32(2 + rand(rng) * (Nj - 4))
    lat, lon = EASLEE.Operators.indyindxtolatlon(gc, y, x)
    if g.crs == "native"          # projection metres from the first mass point
        fx = ((x - 1) * gc.dx - g.lon0) / g.dx; fy = ((y - 1) * gc.dy - g.lat0) / g.dy
    else
        fx = mod(lon - g.lon0, 360.0) / g.dx; fy = (lat - g.lat0) / g.dy
    end
    ie = floor(Int, fx) + 1; je = floor(Int, fy) + 1
    (1 <= ie <= g.nx && 1 <= je <= g.ny) || (outside += 1; continue)
    # keep clear of cell borders by more than float rounding
    (min(fx - floor(fx), ceil(fx) - fx) > 1e-3 && min(fy - floor(fy), ceil(fy) - fy) > 1e-3) || continue
    n += 1
    state.status[slot] = 0x01; state.nest[slot] = 1; state.is_new[slot] = false; state.relbox[slot] = 1
    state.indx[slot] = x; state.indy[slot] = y; state.indl[slot] = 2.5f0; state.mass[slot] = 1.0f0
    state.Dtime[slot] = rc.startdate
    State.sync_dtime_ns!(state)                 # the KA deposit selects parcels by Dtime_ns
    for (path, f!) in (("cpu", EASLEE.Operators.deposit_to_grid!), ("ka", EASLEE.Operators.deposit_to_grid_ka!))
        fill!(oc.massarray, 0.0); fill!(oc.concarray, 0.0f0)
        f!(ctx, state)
        m2 = dropdims(sum(oc.massarray; dims=(3, 4, 5)); dims=(3, 4, 5))
        idx = findall(>(0), m2)
        ok = length(idx) == 1 && Tuple(idx[1]) == (ie, je)
        if !ok
            bad[path] += 1
            bad[path] <= 3 && println(@sprintf("  %s: parcel lon %.4f lat %.4f → expected cell (%d,%d), mass in %s", path, lon, lat, ie, je,
                                               isempty(idx) ? "no cell" : string(Tuple.(idx))))
        end
    end
end
@printf("%d parcels checked (%d outside the output grid skipped): wrong cell CPU %d, KA %d\n", n, outside, bad["cpu"], bad["ka"])
println(bad["cpu"] == 0 && bad["ka"] == 0 && n > 0 ? "PASS" : "FAIL")
exit(bad["cpu"] == 0 && bad["ka"] == 0 && n > 0 ? 0 : 1)
