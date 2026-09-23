# test_backward_srr.jl — backward source-receptor correctness on the OUTPUT grid.
#
# Short check (≈12 steps).  Replicates step_model!'s internal order so the
# residence time can be tallied INDEPENDENTLY of conccalc! (same Dtime match,
# own geometry test), then verifies:
#   1. Σ_cells massarray  ==  independent Σ_t Σ_{in-domain} m_p·Δt/M
#      (the grid faithfully records source-receptor residence);
#   2. Σ sens·V = mean in-domain residence time (the S&F 2004 identity);
#   3. sens_mr / sens ≈ 1/ρ  (mixing-ratio vs concentration receptor);
#   4. parcels move UPWIND (footprint opposite the mean flow).
#
# Usage:  julia --threads=8 --project=../easlee test_backward_srr.jl

using Dates, Printf, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/src/io/output.jl")
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.Operators
using .EASLEE.MetAdapter
using .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")

global gpu = false
npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-50s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

spec = read_runspec(joinpath(@__DIR__, "..", "configs", "era5_backward.toml"))
rt, _obj = build_run(spec; prefetch=false)
ctx   = rt.ctx
state = rt.state
rc    = ctx.runconfig
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end

@assert rc.direction == -1 "this test requires a backward (LDIRECT=-1) config"
oc  = ctx.output
zi  = oc.outinterfaces
M   = Float64(oc.relmass[1])
dt  = abs(rc.timestep.value) / 1000.0      # s

# Independent residence tally: replicate conccalc's accept test for one parcel.
function deposits_indomain(i)
    n  = Int(state.nest[i]); n in 1:length(ctx.nests) || return false
    nc = ctx.nests[n]
    gc = nc.gridconfig
    indt = MetAdapter.met_time_index_array(rc.currentdate, ctx.nests)[n]
    indx = state.indx[i]; indy = state.indy[i]; indl = state.indl[i]
    xout = Operators.indxtoxout(gc, indx, indy)
    yout = Operators.indytoyout(gc, indx, indy)
    ii = Base.unsafe_trunc(Int, xout); jj = Base.unsafe_trunc(Int, yout)
    (1 <= ii <= Int(oc.Niout) && 1 <= jj <= Int(oc.Njout)) || return false
    zout = Operators.indltozlev(nc.levelconfig, nc.metfields3D.zlev, indx, indy, indl, indt)
    (zout >= 0.0f0 && zi[1] <= zout <= zi[end]) || return false
    return true
end

loopstart = rc.enddate
nsteps    = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(rc.timestep)))
ts        = rc.timestep

ref_resid = 0.0        # independent Σ m_p·Δt/M over in-domain deposits
lon_hist  = Float64[]  # mean parcel lon per step (upwind check)
u_hist    = Float64[]

for itime in 0:nsteps
    d = loopstart + itime * ts
    read_met_step!(rt, d)
    rc.currentdate = d
    t_array = MetAdapter.met_time_index_array(d, ctx.nests)

    rel = Operators.release_parcels!(state, ctx, d)
    rel > 0 && Operators.initialize_new_parcels!(state, ctx, d, t_array)

    Operators.deposit_to_grid!(ctx, state)             # → massarray / concarray

    # Independent tally over the SAME parcels conccalc just matched.
    nin = 0; meanlon = 0.0; meanu = 0.0
    @inbounds for i in eachindex(state.Dtime)
        (state.status[i] == 0x01 && state.Dtime[i] == d) || continue
        if deposits_indomain(i)
            global ref_resid += Float64(state.mass[i]) * dt / M
        end
        # diagnostics (all alive parcels, any nest=1)
        if state.nest[i] == 1
            n = 1; gc = ctx.nests[n].gridconfig
            lat, lon = Operators.indyindxtolatlon(gc, Float64(state.indy[i]), Float64(state.indx[i]))
            meanlon += lon; nin += 1
            meanu += Operators.interpolateU(gc.Uedge, ctx.nests[n].metfields3D.Ux,
                        state.indx[i], state.indy[i], clamp(state.indl[i], 1.0f0,
                        Float32(size(ctx.nests[n].metfields3D.zlev,3)-1)), t_array[n])
        end
    end
    if nin > 0
        push!(lon_hist, meanlon / nin)
        push!(u_hist, meanu / nin)
    end

    Operators.advance_active_parcels!(state, ctx, d, t_array)
    Operators.reassign_escaped_parcels!(state, ctx, d)
    State.sync_dtime_ns!(state)
end

# Grid totals
grid_mass = sum(oc.massarray)               # Σ m_p·(Δt/M)  (= Σ sens·V)
grid_vol  = sum(oc.concarray)               # Σ (m_p/ρ)·(Δt/M)

println("── backward source-receptor on the output grid ───────────────────")
@printf("  released particles N        : %d\n", sum(Int, ctx.releases.nbpart))
@printf("  Σ massarray (grid)          : %.4f s\n", grid_mass)
@printf("  independent residence tally : %.4f s\n", ref_resid)
@printf("  mean in-domain residence    : %.1f s of %.0f s run window\n",
        grid_mass, nsteps * dt)

check("grid SRR == independent residence tally",
      isapprox(grid_mass, ref_resid; rtol=1e-4),
      @sprintf("%.5f vs %.5f", grid_mass, ref_resid))
check("Σ sens·V positive and within run window",
      0.0 < grid_mass < nsteps * dt,
      @sprintf("%.0f s", grid_mass))

# sens_mr / sens ≈ mass-weighted <1/ρ> — must be in a physical air-density band
ratio = grid_vol / max(grid_mass, 1e-30)
check("sens_mr/sens = <1/ρ> in physical band (0.5–2 m³/kg)",
      0.5 < ratio < 2.0, @sprintf("%.3f m³/kg → ρ≈%.2f kg/m³", ratio, 1/ratio))

# Upwind: with mean u (here ~westerly/easterly), the parcel centroid must
# drift OPPOSITE to the flow as time runs backward.
dlon = lon_hist[end] - lon_hist[1]
umean = mean(u_hist)
check("parcels advected upwind (Δlon opposes mean u)",
      sign(dlon) == -sign(umean) || abs(umean) < 0.5,
      @sprintf("Δlon=%+.3f°, mean u=%+.1f m/s over %d steps", dlon, umean, length(lon_hist)))

# Write the file and confirm the backward variables/units land on disk.
outfile = "../juliaout/output_era5_back/grid_test_srr.nc"
initnetcdf!(oc, outfile)
writenetcdf(oc, outfile, nsteps + 1)
using NetCDF
sens = ncread(outfile, "spec024_sens")
volc = Float64.(oc.volume)
goodc = isfinite.(volc) .& (volc .> 0)
sumSV = sum(Float64.(sens)[goodc] .* volc[goodc])
check("written sens·V reproduces grid residence",
      isapprox(sumSV, grid_mass; rtol=1e-3), @sprintf("%.3f vs %.3f s", sumSV, grid_mass))
check("backward units attribute on disk",
      NetCDF.open(outfile).vars["spec024_sens"].atts["units"] == "s m-3", "")

println()
@printf("=== %d passed, %d failed ===\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
