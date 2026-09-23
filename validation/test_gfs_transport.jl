# test_gfs_transport.jl — end-to-end transport on GFS ISOBARIC levels.
# Confirms the isobaric → index-space reduction carries through advection:
#   • horizontal parcel displacement ≈ local wind × t (1 h, pure advection);
#   • vertical positions stay finite/sane (no level-ladder artifacts);
#   • parcels at jet level move ≈ jet speed.
#
# Usage:  julia --threads=4 --project=../easlee test_gfs_transport.jl

using Dates, Printf, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")

global gpu = false
npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-52s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

spec = read_runspec(joinpath(@__DIR__,"..","configs","gfs.toml"))
rt, _obj = build_run(spec; prefetch=false, extra_slots=5_000)
ctx   = rt.ctx
state = rt.state
rc    = ctx.runconfig
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end
read_met_step!(rt, rc.startdate)

m3 = ctx.nests[1].metfields3D
gc = ctx.nests[1].gridconfig
nz = size(m3.zlev, 3)

# Manual parcels at three regimes: BL (~900 hPa), mid-trop (~500), jet (~250)
P0 = 2001
levels = [findfirst(k -> m3.puv[720, 360, k, 1] < p, 1:nz) for p in (90000.0f0, 50000.0f0, 25000.0f0)]
spots  = [(400.0f0, 500.0f0), (900.0f0, 300.0f0), (1200.0f0, 450.0f0)]
np = 0
for (ix, iy) in spots, kl in levels
    global np += 1
    i = P0 + np
    state.status[i] = 0x01; state.nest[i] = 1; state.is_new[i] = false
    state.indx[i] = ix; state.indy[i] = iy; state.indl[i] = Float32(kl)
    state.mass[i] = 1.0f0; state.Dtime[i] = rc.startdate
end
ids = collect(P0+1:P0+np)
# Reference: wind at the initial positions; expected displacement over 1 h
uref = [EASLEE.Operators.interpolateU(gc.Uedge, m3.Ux, state.indx[i], state.indy[i], state.indl[i], 1.0f0) for i in ids]
vref = [EASLEE.Operators.interpolateV(gc.Vedge, m3.Vy, state.indx[i], state.indy[i], state.indl[i], 1.0f0) for i in ids]
x0 = [state.indx[i] for i in ids]; y0 = [state.indy[i] for i in ids]
z0 = [EASLEE.Operators.interpolation4Dlevsingle(m3.zlev, state.indx[i], state.indy[i], state.indl[i], 1.0f0) for i in ids]

ts = rc.timestep
for itime in 0:3                                   # 4 × 15 min = 1 h
    d = rc.startdate + itime * ts
    read_met_step!(rt, d)
    step_runtime!(rt, d; do_output=false)
    State.sync_dtime_ns!(state)
end

errs = Float64[]
for (q, i) in enumerate(ids)
    state.status[i] == 0x01 || continue
    # expected index displacement from the initial wind (1 h; tolerance covers
    # wind shear along the trajectory and time evolution)
    dxe = EASLEE.Operators.dxtodindx(gc, uref[q] * 3600.0f0, x0[q], y0[q])
    dye = EASLEE.Operators.dytodindy(gc, vref[q] * 3600.0f0, x0[q], y0[q])
    de  = sqrt(dxe^2 + dye^2)
    da  = sqrt((state.indx[i] - x0[q])^2 + (state.indy[i] - y0[q])^2)
    de > 0.3 && push!(errs, abs(da - de) / de)
end
check("horizontal displacement ≈ u·t (1 h, 9 parcels)",
      !isempty(errs) && median(errs) < 0.3,
      @sprintf("median rel. error %.2f (n=%d)", median(errs), length(errs)))

zend = [EASLEE.Operators.interpolation4Dlevsingle(m3.zlev, state.indx[i], state.indy[i],
        clamp(state.indl[i], 1.0f0, Float32(nz-1)), 1.0f0) for i in ids if state.status[i] == 0x01]
check("vertical positions finite & sane", all(isfinite, zend) && all(0 .<= zend .<= 30000),
      @sprintf("z range %.0f–%.0f m (started %.0f–%.0f)", minimum(zend), maximum(zend),
               minimum(z0), maximum(z0)))
dz = [abs(zend[q] - z0[q]) for q in eachindex(zend)]
check("vertical drift bounded (ω transport, 1 h)", maximum(dz) < 2000,
      @sprintf("max |Δz| = %.0f m", maximum(dz)))

println()
@printf("=== %d passed, %d failed ===\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
