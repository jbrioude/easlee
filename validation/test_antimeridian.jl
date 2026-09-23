# test_antimeridian.jl ── regression test for the 180° (antimeridian) parcel void.
#
# Bug (fixed 2026-06-19): on a GLOBAL grid the met arrays carry a wrap (halo) column
# at Ni+1 (= copy of column 1), so a parcel at grid index indx ∈ [Ni, Ni+1) is a VALID
# seam position that interpolates across 180°.  advance.jl / advance_ka.jl used `indx >= Ni`
# for BOTH the OOB-skip guard and the periodic wrap, which (a) skipped advection for that
# band — freezing Dtime so reassign_escaped_parcels! freed the slot — and (b) wrapped
# legitimate crossings into [0,1) where the next interpolation reads column 0.  Net effect:
# a ~0.25° parcel void at 180° (parcels could not cross going east).  See domainfill /
# chem_cams_pacific.
#
# This test loads the REAL global ERA5 grid, places control + seam-band parcels straddling
# 180°, advances one step through the real advance! kernel, and asserts the seam parcels
# survive with a finite, in-range grid index.  Before the fix the seam parcels are freed
# (status → 0x00) and the test FAILS; after the fix all survive.
#
#   julia -t 1 --project=easlee test/test_antimeridian.jl
#   CFG=configs/era5_conv.toml julia -t 1 --project=easlee test/test_antimeridian.jl

using Dates, Printf
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl"))
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))
using .EASLEE, .EASLEE.State, .EASLEE.MetAdapter, .EASLEE.Operators
const idx2ll = Operators.indyindxtolatlon

const CFG = get(ENV, "CFG", joinpath(REPO, "configs", "era5_conv.toml"))

# ── 1. Pure-arithmetic check of the wrap invariant (no met needed) ────────────
# The fixed wrap must map any indx into [1, Ni+1); the OLD wrap (threshold Ni)
# left values in [Ni, Ni+1) → subtracted to [0,1) → dead next step.
function wrap_fixed(indx::Float32, Ni::Int)
    indx < 1.0f0          && (indx += Float32(Ni))
    indx >= Float32(Ni)+1 && (indx -= Float32(Ni))
    return indx
end
let Ni = 1440, fails = 0
    for x in Float32[0.3, 0.999, 1.0, 720.0, 1439.5, 1440.0, 1440.5, 1440.999, 1441.2, 1441.999]
        w = wrap_fixed(x, Ni)
        ok = 1.0f0 <= w < Float32(Ni)+1
        ok || (fails += 1; @printf("  WRAP FAIL: indx=%.3f → %.3f (outside [1,%d))\n", x, w, Ni+1))
    end
    @printf("[1] wrap arithmetic: %s\n", fails == 0 ? "PASS" : "FAIL ($fails)")
end

# ── 2. Integration: advance seam parcels on the real global ERA5 grid ─────────
println("\n-- build_run (loads first met frames) --"); flush(stdout)
spec = read_runspec(CFG)
for r in spec.release; r.particles = 1; end          # we place our own parcels
rt, _obj = build_run(spec; extra_slots=64, prefetch=false)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
nestctx = ctx.nests[1]; gc = nestctx.gridconfig
Ni = Int(gc.Ni); Nj = Int(gc.Nj); Nk = Int(nestctx.levelconfig.Nk)
loopstart = rc.startdate; ts = rc.timestep

@printf("  grid: Ni=%d Nj=%d Nk=%d globalx=%s lonmin=%.3f dlon=%.4f\n",
        Ni, Nj, Nk, gc.globalx, gc.lonmin, gc.dlon)
if !gc.globalx
    @error "config $(CFG) is not a global-in-longitude grid (globalx=false) — cannot test the antimeridian seam."
    exit(1)
end

read_met_step!(rt, loopstart)
indt_arr = MetAdapter.met_time_index_array(loopstart, ctx.nests)

# Reset every slot to free, sentinel Dtime, so ONLY our placed parcels advance.
fill!(state.status, 0x00)
fill!(state.Dtime, DateTime(1))

# Place parcels across 3 latitude rows × a set of x-positions spanning the seam.
# Positions are grid indices: lon = lonmin + (indx-1)*dlon.  The band [Ni, Ni+1)
# (lon ∈ [180-dlon, 180)) is the one the OLD code killed.
jrows  = Float32[div(Nj,3), div(Nj,2), 2*div(Nj,3)]
xpos   = Float32[Ni-1.5, Ni-0.5,           # controls, inside the data columns
                 Ni-0.001, Ni+0.001, Ni+0.5, Ni+0.999]  # seam band [Ni, Ni+1)
klev   = Float32(max(2, div(Nk,2)))
slots  = Int[]; s = 0
for jr in jrows, xp in xpos
    global s += 1
    @inbounds begin
        state.status[s]    = 0x01
        state.nest[s]      = Int8(1)
        state.indx[s]      = xp
        state.indy[s]      = jr
        state.indl[s]      = klev
        state.Dtime[s]     = loopstart
        state.Dtimeinit[s] = loopstart
        state.is_new[s]    = false
        state.mass[s]      = 1.0f0
        state.dt[s]        = 1.0f0
        state.uturb[s]=0f0; state.vturb[s]=0f0; state.wturb[s]=0f0
        state.miss2[s]=0f0; state.miss3[s]=0f0; state.miss4[s]=0f0
    end
    push!(slots, s)
end
@printf("  placed %d test parcels (%d lat rows × %d x-positions) at level %.0f\n",
        length(slots), length(jrows), length(xpos), klev)

# One real advance step + the escape/free pass that would cull a frozen parcel.
Operators.advance_active_parcels!(state, ctx, loopstart, indt_arr)
Operators.reassign_escaped_parcels!(state, ctx, loopstart)

# ── Verdict ───────────────────────────────────────────────────────────────────
next_date = loopstart + ts
nfail = 0; nseam = 0; nseam_ok = 0
println("\n  slot   indx0       lon0      indx1       lon1     alive  note")
for (k, sl) in enumerate(slots)
    global nfail, nseam, nseam_ok
    x0 = xpos[((k-1) % length(xpos)) + 1]
    seam = x0 >= Float32(Ni)             # in the band the old code killed
    seam && (nseam += 1)
    alive = state.status[sl] == 0x01 && state.Dtime[sl] == next_date
    x1 = state.indx[sl]
    lon0 = gc.lonmin + (x0 - 1) * gc.dlon
    lon1 = alive ? gc.lonmin + (x1 - 1) * gc.dlon : NaN
    inrange = alive && isfinite(x1) && 1.0f0 <= x1 < Float32(Ni)+1
    ok = alive && inrange
    seam && ok && (nseam_ok += 1)
    ok || (nfail += 1)
    note = !alive ? "DEAD (freed)" : !inrange ? "OUT OF RANGE" : (seam ? "seam ok" : "ok")
    @printf("  %4d  %9.3f  %8.3f  %9.3f  %8.3f   %-5s  %s\n",
            sl, x0, lon0, x1, lon1, alive ? "yes" : "NO", note)
end

@printf("\n[2] seam parcels: %d/%d survived & in-range\n", nseam_ok, nseam)
if nfail == 0
    println("\n✅ ANTIMERIDIAN TEST PASS — no 180° void (all seam parcels advected & survived).")
    exit(0)
else
    @printf("\n❌ ANTIMERIDIAN TEST FAIL — %d parcel(s) died or left valid range at the 180° seam.\n", nfail)
    exit(1)
end
