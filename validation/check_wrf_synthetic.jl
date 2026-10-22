# check_wrf_synthetic.jl — EASLEE's WRF reader and Lambert transport against known truth.
#
#     python3 validation/make_synthetic_wrfout.py GFS_f003 GFS_f006 /data/.../WRF_synth
#     julia -t 8 --project=easlee validation/check_wrf_synthetic.jl configs/<wrf>.toml /data/.../WRF_synth_truth
#
# The synthetic wrfout files carry WRF's conventions (staggered grid-relative winds,
# theta' = theta - 300, P + PB, PH + PHB on the interfaces, QVAPOR mixing ratio) and
# the truth on the mass points next to them.  Checked here:
#   1. reader: T, p, height above ground, q at slot k+1 against mass level k of the
#      truth; u, v (grid-relative, destaggered) against the truth rotated to the
#      grid; vertical speed implied by levdot against w;
#   2. transport: parcels advected one step without turbulence must move by the
#      local wind on the sphere — speed ratio 1 and no direction error.  Without
#      the map factor the ratio is 1/m (1.01-1.04 on this grid).
# Exit status 1 when a check fails.

using Dates, Printf, Statistics, Random, NetCDF
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

cfg, tdir = ARGS[1], ARGS[2]
nfail = 0
function check(label, ok, detail)
    global nfail
    ok || (nfail += 1)
    @printf("  %-60s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end
rel(a, b) = sqrt(mean((a .- b) .^ 2)) / max(sqrt(mean(b .^ 2)), 1e-12)

spec = read_runspec(cfg)
rt, obj = build_run(spec; prefetch=false, extra_slots=20_000)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig
read_met_step!(rt, rc.startdate)
nc_ = ctx.nests[1]; gc = nc_.gridconfig
d3 = rt.metfields3D[1]; d2 = rt.metfields2D[1]
Ni, Nj = Int(gc.Ni), Int(gc.Nj)
vdate = nc_.metconfig.filedate[nc_.nestrunconfig.t_index[1]]
tf = joinpath(tdir, "truth_" * Dates.format(vdate, "yyyy-mm-dd_HH:MM:SS") * ".nc")
println("truth: ", tf)
tr(name) = Float64.(ncread(tf, name))                    # (i, j, k) / (i, j)
Tt, pt, zt, qt, ut, vt, wt = tr("T"), tr("p"), tr("zm"), tr("q"), tr("u"), tr("v"), tr("w")
lon = tr("lon"); mf = tr("mapfac")
cone = Float64(ncgetatt(tf, "Global", "cone")); slon = Float64(ncgetatt(tf, "Global", "stand_lon"))
Km = size(Tt, 3)
I, J = 2:Ni-1, 2:Nj-1                                     # interior (destaggering is exact only inside)
s(name) = Float64.(d3[name][I, J, 2:Km+1, 1])            # EASLEE slots 2 … Km+1 ↔ mass levels 1 … Km
t(a) = a[I, J, :]

println("\n1. reader against the truth (mass levels 1-$Km, interior points)")
check("temperature T (not theta - 300)", rel(s("T"), t(Tt)) < 1e-3, @sprintf("rel RMS %.2e; EASLEE %.1f-%.1f K", rel(s("T"), t(Tt)), extrema(s("T"))...))
check("pressure puv = P + PB", rel(s("puv"), t(pt)) < 1e-4, @sprintf("rel RMS %.2e", rel(s("puv"), t(pt))))
# stacked layout (the transport code's): zlev on the HALF levels = the w-levels,
# T, p, u, v, q on the mass levels stacked between them (slot k+1 = mass level k)
zwt = tr("zw"); wit = tr("wi")                          # truth interfaces (1 = ground)
zE = Float64.(d3["zlev"][I, J, 1:Km+1, 1])
check("zlev = heights of the w-levels above ground (stacked)", rel(zE[:, :, 2:end], zwt[I, J, 2:Km+1]) < 1e-3,
      @sprintf("rel RMS %.2e; first w-level above ground: EASLEE %.1f m, truth %.1f m (median)",
               rel(zE[:, :, 2:end], zwt[I, J, 2:Km+1]), median(zE[:, :, 2]), median(zwt[I, J, 2])))
check("specific humidity (QVAPOR mixing ratio converted)", rel(s("qv"), t(qt)) < 1e-3, @sprintf("rel RMS %.2e", rel(s("qv"), t(qt))))
α = cone .* deg2rad.(((lon .- slon .+ 180) .% 360) .- 180)
ug = ut .* cos.(α) .- vt .* sin.(α); vg = ut .* sin.(α) .+ vt .* cos.(α)
check("u (grid-relative, destaggered)", rel(s("Ux"), t(ug)) < 0.05, @sprintf("rel RMS %.3f", rel(s("Ux"), t(ug))))
check("v (grid-relative, destaggered)", rel(s("Vy"), t(vg)) < 0.05, @sprintf("rel RMS %.3f", rel(s("Vy"), t(vg))))
z3 = Float64.(d3["zlev"][:, :, :, 1]); lv = Float64.(d3["levdot"][:, :, :, 1])
# levdot on the w-levels s = 2 … Km (half levels), as a vertical speed
wE = [lv[i, j, s] * (z3[i, j, s+1] - z3[i, j, s-1]) / 2 for i in I, j in J, s in 2:Km]
# it is the COORDINATE velocity on the terrain-following levels:
# w − u·∂z_s/∂x − v·∂z_s/∂y on the w-level, u, v the mean of the two mass levels
# around it, z_s its height ASL, slopes over the true ground distance DX/m
hgt = tr("hgt"); dxm = 30000.0
zwa = zwt .+ reshape(hgt, size(hgt)..., 1)                       # (i, j, s) ASL
uh(i, j, s) = 0.5 * (ug[i,j,s-1] + ug[i,j,s]); vh(i, j, s) = 0.5 * (vg[i,j,s-1] + vg[i,j,s])
wc = [wit[i,j,s] - uh(i,j,s) * (zwa[i+1,j,s] - zwa[i-1,j,s]) * mf[i,j] / (2dxm) -
                   vh(i,j,s) * (zwa[i,j+1,s] - zwa[i,j-1,s]) * mf[i,j] / (2dxm) for i in I, j in J, s in 2:Km]
check("vertical speed from levdot = coordinate velocity (w-levels)", rel(wE, wc) < 0.1,
      @sprintf("rel RMS %.3f (vs physical w: %.3f)", rel(wE, wc), rel(wE, wit[I, J, 2:Km])))
check("tke filled from UST/PBLH/HFX", maximum(d3["tke"][:, :, :, 1]) > 0, @sprintf("max %.3f", maximum(d3["tke"][:, :, :, 1])))
check("theta filled", minimum(d3["theta"][I, J, 2:Km+1, 1]) > 250, @sprintf("min %.1f K", minimum(d3["theta"][I, J, 2:Km+1, 1])))

println("\n2. one advection step on the Lambert grid (no turbulence)")
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end
rng = MersenneTwister(1); P0 = length(state.status) - 5000
ids = Int[]; lat0 = Float64[]; lon0 = Float64[]; ue = Float64[]; ve = Float64[]; m0 = Float64[]
for q in 1:4000
    i = P0 + q
    x = Float32(10 + rand(rng) * (Ni - 20)); y = Float32(10 + rand(rng) * (Nj - 20)); l = Float32(8 + rand(rng) * 22)
    state.status[i] = 0x01; state.nest[i] = 1; state.is_new[i] = false; state.relbox[i] = 1
    state.indx[i] = x; state.indy[i] = y; state.indl[i] = l; state.mass[i] = 1.0f0; state.Dtime[i] = rc.startdate
    la, lo = EASLEE.Operators.indyindxtolatlon(gc, y, x)
    u = EASLEE.Operators.interpolateU(gc.Uedge, nc_.metfields3D.Ux, x, y, l, 1.0f0)
    v = EASLEE.Operators.interpolateV(gc.Vedge, nc_.metfields3D.Vy, x, y, l, 1.0f0)
    a = cone * deg2rad(((lo - slon + 180) % 360) - 180)
    push!(ids, i); push!(lat0, la); push!(lon0, lo)
    push!(ue, u * cos(a) + v * sin(a)); push!(ve, -u * sin(a) + v * cos(a))
    # map factor from the truth file (independent of the model's own formula)
    push!(m0, mf[clamp(round(Int, x), 1, Ni), clamp(round(Int, y), 1, Nj)])
end
step_runtime!(rt, rc.startdate; do_output=false)
dt = Dates.value(rc.timestep) / 1000
ratio = Float64[]; dang = Float64[]
for (q, i) in enumerate(ids)
    state.status[i] == 0x01 || continue
    la, lo = EASLEE.Operators.indyindxtolatlon(gc, state.indy[i], state.indx[i])
    dN = deg2rad(la - lat0[q]) * 6.37e6; dE = deg2rad(((lo - lon0[q] + 180) % 360) - 180) * 6.37e6 * cosd((la + lat0[q]) / 2)
    sp = hypot(ue[q], ve[q]) * dt
    sp > 2000 || continue                                  # only well-defined displacements (> 2 km)
    push!(ratio, hypot(dE, dN) / sp)
    push!(dang, rad2deg(angle(complex(dE, dN) / complex(ue[q], ve[q]))))
end
@printf("  parcels %d; map factor at the start points %.4f-%.4f\n", length(ratio), extrema(m0)...)
check("displacement / (wind x dt) = 1 (map factor applied)", abs(median(ratio) - 1) < 0.005,
      @sprintf("median %.4f, p5-p95 %.4f-%.4f  (no map factor would give 1/m: %.4f)", median(ratio),
               quantile(ratio, 0.05), quantile(ratio, 0.95), median(1 ./ m0)))
check("direction error below 1 degree (grid-relative winds)", abs(median(dang)) < 1.0,
      @sprintf("median %.3f deg, p95 |err| %.2f deg", median(dang), quantile(abs.(dang), 0.95)))
println(nfail == 0 ? "\nALL CHECKS PASS" : "\n$nfail CHECK(S) FAIL")
exit(nfail == 0 ? 0 : 1)
