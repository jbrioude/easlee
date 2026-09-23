#!/usr/bin/env julia
# test_wind_dim.jl — validate [run].wind freezes the suppressed parcel indices.
# Builds a tiny run for each wind mode, steps a few times, and checks that the
# frozen dimensions (indx/indy/indl) do NOT move while the free ones do.
#   julia -t 2 --project=easlee test/test_wind_dim.jl
using Dates, Printf
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO,"easlee","src","config","EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO,"easlee","src","io","initialize.jl"))
include(joinpath(REPO,"easlee","src","io","tools.jl"))
include(joinpath(REPO,"easlee","src","io","output.jl"))
include(joinpath(REPO,"easlee","EASLEE.jl")); using .EASLEE
include(joinpath(REPO,"easlee","runtime.jl"))
include(joinpath(REPO,"easlee","build_run.jl"))

const CFG = joinpath(REPO, "configs", "test_wind.toml")

function probe(wind::String)
    spec = read_runspec(CFG); spec.wind = wind
    rt, obj = build_run(spec; prefetch=false)
    ctx = rt.ctx; st = rt.state; rc = ctx.runconfig; ts = rc.timestep
    kw = step_kwargs(obj)
    t0 = rc.startdate
    read_met_step!(rt, t0)
    EASLEE.step_model!(st, ctx, t0; kw...)                 # release + first advance
    EASLEE.State.sync_dtime_ns!(st)
    act = findall(i -> st.status[i] == 0x01, eachindex(st.status))
    ix = copy(st.indx[act]); iy = copy(st.indy[act]); il = copy(st.indl[act])
    for k in 1:6
        t = t0 + k*ts
        read_met_step!(rt, t)
        EASLEE.step_model!(st, ctx, t; kw...)
        EASLEE.State.sync_dtime_ns!(st)
    end
    keep = [a for a in eachindex(act) if st.status[act[a]] == 0x01]   # still alive
    didx = maximum(abs.(st.indx[act[keep]] .- ix[keep]); init=0.0f0)
    didy = maximum(abs.(st.indy[act[keep]] .- iy[keep]); init=0.0f0)
    didl = maximum(abs.(st.indl[act[keep]] .- il[keep]); init=0.0f0)
    @printf("wind=%-12s  N=%-5d  max|Δindx|=%.4f  max|Δindy|=%.4f  max|Δindl|=%.4f\n",
            wind, length(keep), didx, didy, didl)
    return (didx, didy, didl)
end

println("\n── wind dimensionality check (frozen index should give Δ≈0) ──")
d3 = probe("3D")
dz = probe("1D-vertical")    # freeze indx & indy
dh = probe("2D-lonlat")      # freeze indl
println("\nexpect: 3D moves all; 1D-vertical Δindx≈Δindy≈0 (Δindl>0); 2D-lonlat Δindl≈0")
ok = dz[1] < 1e-3 && dz[2] < 1e-3 && dz[3] > 1e-3 && dh[3] < 1e-3 && d3[1] > 1e-3
println(ok ? "PASS ✓" : "FAIL ✗")
