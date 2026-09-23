# plot_convection_wetdep.jl — rerun the 3 h convection+wetdep test and dump
# arrays for plotting (see plot_convection_wetdep.py).
using Dates, Printf, Statistics
include("../../easlee/src/io/initialize.jl")
include("../../easlee/src/io/tools.jl")
include("../../easlee/src/io/output.jl")
include("../../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../../easlee/runtime.jl")
include("../../easlee/build_run.jl")

global gpu  = false
global chem = false

outdir = "../juliaout/output_era5_conv"
mkpath(outdir)

spec = read_runspec(joinpath(@__DIR__,"..","..","configs","era5_conv.toml")); spec.physics.deposition.wet = true
rt, _obj = build_run(spec; prefetch=false)
ctx   = rt.ctx
state = rt.state
rc    = ctx.runconfig
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end

convfields = [EASLEE.ConvFields(rt.metfields3D[n]) for n in 1:rt.nest]
wdparams   = EASLEE.WetDepParams()
wdstate    = EASLEE.WetDepState(ctx)

# Parcel height above ground from its column's zlev profile.
function parcel_heights(state, ctx, indt)
    zarr = ctx.nests[1].metfields3D.zlev
    nz   = size(zarr, 3)
    hs = Float32[]
    for i in eachindex(state.status)
        state.status[i] == 0x01 || continue
        k  = clamp(floor(Int, state.indl[i]), 1, nz - 1)
        f  = state.indl[i] - k
        ii = clamp(round(Int, state.indx[i]), 1, size(zarr, 1))
        jj = clamp(round(Int, state.indy[i]), 1, size(zarr, 2))
        push!(hs, zarr[ii, jj, k, indt] * (1 - f) + zarr[ii, jj, k + 1, indt] * f)
    end
    hs
end

ts     = rc.timestep
nsteps = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(ts)))
h_init = Float32[]
for itime in 0:nsteps
    currentdate = rc.startdate + itime * ts
    read_met_step!(rt, currentdate)
    EASLEE.step_model!(state, ctx, currentdate;
                       convfields=convfields, wetdepparams=wdparams, wetdepstate=wdstate)
    State.sync_dtime_ns!(state)
    itime == 4 && (global h_init = parcel_heights(state, ctx, 1))  # all parcels released by 08:00
end
h_final = parcel_heights(state, ctx, 1)

# ── Dumps (raw Float32/Float64 binary; dims encoded in filename) ──────────────
m3  = rt.metfields3D[1]
dep = wdstate.depgrid[1]                       # (Ni+1, Nj) Float64
# Column-max precip rate [mm/h] over the run window (use last loaded slot 1)
pmax = dropdims(maximum(m3["tpf"][:,:,:,1]; dims=3); dims=3) .* 3600f0

write(joinpath(outdir, "depgrid_$(size(dep,1))x$(size(dep,2)).f64"), Float64.(dep))
write(joinpath(outdir, "precmax_$(size(pmax,1))x$(size(pmax,2)).f32"), Float32.(pmax))
write(joinpath(outdir, "h_init_$(length(h_init)).f32"),  h_init)
write(joinpath(outdir, "h_final_$(length(h_final)).f32"), h_final)

# Convective column profiles at the strongest column
umf = m3["umf"]; dmf = m3["dmf"]; zl = m3["zlev"]; rho = m3["rho"]; udr = m3["udr"]
colsum = dropdims(sum(umf[:,:,:,1]; dims=3); dims=3)
ci, cj = Tuple(argmax(colsum))
nz = size(zl, 3)
prof = zeros(Float32, nz, 5)
for k in 1:nz
    prof[k,1] = zl[ci,cj,k,1]
    prof[k,2] = max(umf[ci,cj,k,1], 0f0)
    prof[k,3] = max(-dmf[ci,cj,k,1], 0f0)
    prof[k,4] = udr[ci,cj,k,1]
    prof[k,5] = -(max(umf[ci,cj,k,1],0f0) - max(-dmf[ci,cj,k,1],0f0)) / max(rho[ci,cj,k,1],1f-3)
end
write(joinpath(outdir, "colprof_$(nz)x5.f32"), prof)

# Grid metadata for georeferencing the maps
lons = rt.metconfig[1]["longitudes"]; lats = rt.metconfig[1]["latitudes"]
open(joinpath(outdir, "gridmeta.txt"), "w") do io
    println(io, lons[1], " ", lons[2]-lons[1], " ", length(lons))
    println(io, lats[1], " ", lats[2]-lats[1], " ", length(lats))
    println(io, ci, " ", cj)
end
@printf("dumps written to %s  (deposited=%.4f)\n", outdir, sum(dep))
