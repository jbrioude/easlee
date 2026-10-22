# check_vertical_kinematics.jl — over one short step, does a parcel's height change
# equal the archived vertical velocity?
#
#     julia -t 8 --project=easlee validation/check_vertical_kinematics.jl configs/<run>.toml MODE
#
# MODE = mesonh | arome | merra2.  The configuration must have a short timestep (60 s)
# and start at a met time.  Parcels at several level indices over land are advanced
# ONE step without turbulence; whatever the vertical coordinate, the height above sea
# level of a parcel then changes by w·dt (MERRA-2: its pressure by ω·dt), w taken
# straight from the met file at the start position:
#   regression slope of Δz on w·dt ≈ 1  → the terrain-following motion is counted once
#                                ≈ 2  → counted twice (physical w used as levdot)
# The reference w: MesoNH WT (w-levels, averaged onto the mass levels), AROME NH_dW
# integrated from the kinematic ground value with the file's sign, MERRA-2 OMEGA.
# Written 2026-10-13 with the coordinate-velocity fix.

using Dates, Printf, Random, Statistics, NetCDF
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE, .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")
global gpu = false
const OP = EASLEE.Operators

cfg, mode = ARGS[1], ARGS[2]
spec = read_runspec(cfg)
rt, obj = build_run(spec; prefetch=false, extra_slots=20_000)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
read_met_step!(rt, rc.startdate)
for n in 1:length(ctx.nests); ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb(); end
nc_ = ctx.nests[1]; gc = nc_.gridconfig; Ni, Nj = Int(gc.Ni), Int(gc.Nj)
d3 = rt.metfields3D[1]; d2 = rt.metfields2D[1]; m3 = nc_.metfields3D
Nk = size(d3["zlev"], 3)
tix = nc_.nestrunconfig.t_index
slot = findfirst(s -> nc_.metconfig.filedate[tix[s]] == rc.startdate, 1:2)
slot === nothing && error("start the run at a met time (bracket $(nc_.metconfig.filedate[tix[1]]) – $(nc_.metconfig.filedate[tix[2]]))")
fname = nc_.metconfig.filename[tix[slot]]; vdate = rc.startdate; TS = Float32(slot)   # time index of that slot
dt = Dates.value(rc.timestep) / 1000
println("file ", basename(fname), "  dt = $dt s")

# reference w (or ω) in the EASLEE slot layout (i, j, k), STACKED: the vertical
# motion on the half levels (MesoNH w-level s+1, AROME the interface above layer s−1)
ref = zeros(Float32, Ni, Nj, Nk); terrain = zeros(Float32, Ni, Nj)
if mode == "mesonh"
    W = Float32.(ncread(fname, "WT")[:, :, :, 1])                 # (i, j, level_w)
    for s in 1:Nk-1; ref[:, :, s] .= W[:, :, s+1]; end       # stacked: half level s = w-level s+1
    terrain .= Float32.(ncread(fname, "ZS")[:, :])
elseif mode == "arome"
    D = Float32.(ncread(fname, "NH_dW")); D = D[:, :, end:-1:1]  # bottom-up layers
    terrain .= Float32.(ncread(fname, "Zg")) ./ 9.81f0
    U = Float32.(ncread(fname, "U_X")); V = Float32.(ncread(fname, "V_Y"))
    dx = Float64(rt.metconfig[1]["dx"])
    for j in 2:Nj-1, i in 2:Ni-1
        wb = U[i, j, end] * (terrain[i+1, j] - terrain[i-1, j]) / (2dx) + V[i, j, end] * (terrain[i, j+1] - terrain[i, j-1]) / (2dx)
        for k in 1:Nk-1
            wt = wb + D[i, j, k] / 9.81f0; ref[i, j, k+1] = wt; wb = wt   # stacked: w on the interfaces
        end
    end
elseif mode == "merra2"
    t = ncread(fname, "time"); it = argmin(abs.(t .- Dates.value(Minute(vdate - DateTime(Dates.year(vdate), Dates.month(vdate), Dates.day(vdate), 1, 30)))))
    O = Float32.(ncread(fname, "OMEGA")[:, :, end:-1:1, it])        # bottom-up mid levels
    for k in 2:Nk; ref[1:size(O, 1), :, k] .= O[:, :, k-1]; end      # full levels
    for k in 2:Nk-1; ref[:, :, k] .= 0.5f0 .* (ref[:, :, k] .+ ref[:, :, k+1]); end   # → half levels (stacked)
    ref[:, :, Nk] .= 0
    terrain .= Float32.(ncread(fname, "PHIS")[:, :, it]) ./ 9.81f0
else
    error("MODE = mesonh | arome | merra2")
end
ref4 = cat(ref, ref; dims=4)                                       # for the 4-D interpolator
ter3 = cat(terrain, terrain; dims=3)
vcoord(x, y, l) = mode == "merra2" ? OP.interpolation4Dsingle(m3.puv, x, y, l, TS) :
                  OP.interpolation4Dlevsingle(m3.zlev, x, y, l, TS) + OP.interpolate2Dsingle(ter3, x, y, TS)
# the reference rate is interpolated like the model's levdot (level interpolator)
refat(x, y, l) = OP.interpolation4Dlevsingle(ref4, x, y, l, TS)

rng = MersenneTwister(9); P0 = length(state.status) - 15000
fill!(state.status, 0x00)
rows = Tuple{Int,Float32,Float64,Float64}[]          # (parcel, level, start coordinate, reference rate·dt)
for lvl in (2.5f0, 4.5f0, 8.5f0, 15.5f0, 25.5f0)
    n = 0
    for q in 1:200000
        n >= 1500 && break
        x = Float32(4 + rand(rng) * (Ni - 8)); y = Float32(4 + rand(rng) * (Nj - 8))
        i, j = round(Int, x), round(Int, y)
        hypot(terrain[i+1, j] - terrain[i-1, j], terrain[i, j+1] - terrain[i, j-1]) > 100 || continue
        k = P0 + length(rows) + 1; n += 1
        state.status[k] = 0x01; state.nest[k] = 1; state.is_new[k] = false; state.relbox[k] = 1
        state.indx[k] = x; state.indy[k] = y; state.indl[k] = lvl; state.mass[k] = 1.0f0; state.Dtime[k] = rc.startdate
        push!(rows, (k, lvl, vcoord(x, y, lvl), refat(x, y, lvl) * dt))
    end
end
State.sync_dtime_ns!(state)
step_runtime!(rt, rc.startdate; do_output=false)
@printf("%s: Δ%s over one %.0f-s step vs the archived %s·dt, parcels over terrain (> 100 m across two cells)\n",
        mode, mode == "merra2" ? "p" : "z(ASL)", dt, mode == "merra2" ? "ω" : "w")
for lvl in unique(r[2] for r in rows)
    sel = [r for r in rows if r[2] == lvl && state.status[r[1]] == 0x01]
    dz = [vcoord(state.indx[r[1]], state.indy[r[1]], state.indl[r[1]]) - r[3] for r in sel]
    wd = [r[4] for r in sel]
    b = cov(wd, dz) / var(wd)
    @printf("  level index %5.1f: n=%4d  slope %+.3f  corr %+.3f  (1 = counted once, 2 = twice)\n", lvl, length(sel), b, cor(wd, dz))
end
