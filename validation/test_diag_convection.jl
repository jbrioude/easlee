# test_diag_convection.jl — CAPE-gated diagnostic convection (MPTRAC-style)
# on GFS data: find a high-CAPE column, fill it with parcels uniform in air
# mass, apply the scheme, and verify (a) redistribution happens below the
# equilibrium level, (b) levels above the EL are untouched, (c) repeated
# application keeps the distribution well-mixed (uniform in air mass).
#
# Usage:  julia --threads=8 --project=../easlee test_diag_convection.jl

using Dates, Printf, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.MetAdapter
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")

global gpu = false

spec = read_runspec(joinpath(@__DIR__,"..","configs","gfs.toml"))
rt, _obj = build_run(spec; prefetch=false, extra_slots=50_000)
ctx   = rt.ctx
state = rt.state
read_met_step!(rt, ctx.runconfig.startdate)

diagconv = [EASLEE.DiagConvFields(rt.metfields2D[n]) for n in 1:rt.nest]
@assert diagconv[1] !== nothing "cape not mapped in the GFS namelist?"

m2 = rt.metfields2D[1]; m3 = rt.metfields3D[1]
cape = m2["cape"][:,:,1]
Ni, Nj = size(cape)
nz = size(m3["zlev"], 3)

# Strongest-CAPE column (away from the grid edge)
cmax, idxm = findmax([isfinite(cape[i,j]) ? cape[i,j] : -1.0f0
                      for i in 2:Ni-1, j in 2:Nj-1])
ci, cj = Tuple(idxm) .+ 1
@printf("max CAPE = %.0f J/kg at i=%d j=%d  (cin = %.0f)\n", cmax, ci, cj, m2["cin"][ci,cj,1])

zprof   = [m3["zlev"][ci,cj,k,1] for k in 1:nz]
rhoprof = [m3["rho"][ci,cj,k,1]  for k in 1:nz]
layerm  = [0.5f0*(rhoprof[k]+rhoprof[k+1])*(zprof[k+1]-zprof[k]) for k in 1:nz-1]
cm      = cumsum(layerm); cmn = cm ./ cm[end]

N = 50_000
for p in 1:N
    u = rand(Float32)
    k = searchsortedfirst(cmn, u)
    f = k == 1 ? u/cmn[1] : (u-cmn[k-1])/(cmn[k]-cmn[k-1])
    state.status[p] = 0x01
    state.nest[p]   = 1
    state.indx[p]   = Float32(ci)
    state.indy[p]   = Float32(cj)
    state.indl[p]   = Float32(k) + f
    state.mass[p]   = 1.0f0
end

t_array = MetAdapter.met_time_index_array(ctx.runconfig.startdate, ctx.nests)
params  = EASLEE.DiagConvParams(cape_min=100.0f0, cin_max=1.0f6)  # ignore CIN gate here

z_of(p) = begin
    k = clamp(floor(Int, state.indl[p]), 1, nz-1)
    zprof[k] + (state.indl[p]-k)*(zprof[k+1]-zprof[k])
end
z0 = [z_of(p) for p in 1:N]

EASLEE.diagnostic_convection!(state, ctx, diagconv, 900.0, t_array; params=params)
z1 = [z_of(p) for p in 1:N]

moved = count(p -> abs(z1[p]-z0[p]) > 1, 1:N)
zel_est = maximum(z1[p] for p in 1:N if abs(z1[p]-z0[p]) > 1)
@printf("\nafter 1 application: %.1f%% of parcels redistributed; max destination z = %.0f m (≈ EL)\n",
        100moved/N, zel_est)
@assert moved > 0 "no parcels redistributed — EL/threshold logic broken?"
above_moved = count(p -> z0[p] > zel_est*1.05 && abs(z1[p]-z0[p]) > 1, 1:N)
println("parcels above the EL that moved: $above_moved (must be 0)")

# Well-mixed: 20 more applications, compare histogram below EL vs air mass
for s in 1:20
    EASLEE.diagnostic_convection!(state, ctx, diagconv, 900.0, t_array; params=params)
end
z2 = [z_of(p) for p in 1:N]
nb = 12
edges = range(0, zel_est, length=nb+1)
hsel  = [z for z in z2 if z < zel_est]
hcnt  = [count(z -> edges[b] <= z < edges[b+1], hsel) for b in 1:nb]
expm  = Float64[]
for b in 1:nb
    mb = 0.0
    for k in 1:nz-1
        lo = max(zprof[k], edges[b]); hi = min(zprof[k+1], edges[b+1])
        hi > lo && (mb += layerm[k] * (hi-lo)/(zprof[k+1]-zprof[k]))
    end
    push!(expm, mb)
end
expn = expm ./ sum(expm) .* length(hsel)
println("\n  z bin [m]        expected   final     dev")
maxdev = 0.0
for b in 1:nb
    dev = abs(hcnt[b]-expn[b])/expn[b]
    global maxdev = max(maxdev, dev)
    @printf("  %5.0f–%-6.0f  %8.0f  %7d  %5.1f%%\n", edges[b], edges[b+1], expn[b], hcnt[b], 100dev)
end
tol = 5/sqrt(minimum(expn)) + 0.03
@printf("\nmax deviation %.1f%%  (tolerance %.1f%%)  —  %s\n",
        100maxdev, 100tol, maxdev < tol ? "PASS" : "FAIL")
exit(maxdev < tol && above_moved == 0 ? 0 : 1)
