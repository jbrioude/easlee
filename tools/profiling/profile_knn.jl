#!/usr/bin/env julia
# profile_knn.jl — allocation/time profiling of the adaptive-kernel kNN pass
# (update_bandwidths!) and the adaptive deposit (_conc_nest! :adaptive), on a
# synthetic 400k ETEX-like filament.  Motivation: the etex_ab bench shows the
# :output group at 522.9 s with 1.35 TiB allocated — the intended buffers are
# ~16 MB/refresh, so ~everything else is a hot-loop allocation bug.
#
#   julia --project=easlee -t 4 scripts/profile_knn.jl
#
# Env knobs: KNN_N (default 400000), ALLOC_RATE (Profile.Allocs sample rate,
# default 5e-5), REPS (default 3).
using Printf, Random, Statistics, Dates
using Profile

global gpu = false; global chem = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: update_bandwidths!, set_output_popmap!, _conc_nest!

const N     = parse(Int, get(ENV, "KNN_N", "400000"))
const RATE  = parse(Float64, get(ENV, "ALLOC_RATE", "5e-5"))
const REPS  = parse(Int, get(ENV, "REPS", "3"))

# ── synthetic ETEX-like setup ─────────────────────────────────────────────────
# 0.1°-class output geometry: 400×300 met grid, 1° spacing constants as the
# unit tests (only the metre constants matter for the algorithm).
G = EASLEE.Types.Rectilinear(Int32(400), Int32(300), Int32(10), 1.0f0, 1.0f0,
        -10.0f0, 40.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
zlev = Array{Float32,4}(undef, 400, 300, 10, 2)
for k in 1:10; zlev[:, :, k, :] .= Float32((k - 1) * 500.0); end
rho = fill(1.0f0, 400, 300, 10, 2)
zi  = Float32[500.0f0 * k for k in 1:2]
t0  = DateTime(2025, 1, 1)

function filament_state(n)
    st = EASLEE.State.ParcelState(n)
    st.status .= 0x01; st.relbox .= Int32(1); st.nest .= Int32(1)
    st.indl .= 2.5f0; st.mass .= 1.0f0; st.Dtime .= t0
    Random.seed!(7)
    nd = n ÷ 10                      # 10 % dense near-source clump (stride path)
    ct, sn = cosd(30.0), sind(30.0)
    for i in 1:n
        if i <= nd
            st.indx[i] = 60.0f0 + Float32(0.05 * randn())
            st.indy[i] = 60.0f0 + Float32(0.05 * randn())
        else
            a = 120.0 * randn(); b = 8.0 * randn()      # filament: 120×8 idx units
            st.indx[i] = clamp(Float32(200.0 + a * ct - b * sn), 2.0f0, 398.0f0)
            st.indy[i] = clamp(Float32(150.0 + a * sn + b * ct), 2.0f0, 298.0f0)
        end
    end
    return st
end

println("── profile_knn: N=$N threads=$(Threads.nthreads()) rate=$RATE ──")
set_output_popmap!([1])              # mirror the etex_ab run (one population)

# warmup / compile on a small state, then the real one
stw = filament_state(20_000)
update_bandwidths!(stw, G, 1)
st = filament_state(N)
update_bandwidths!(st, G, 1)         # warm caches, compile the big path

println("\n[1] kNN pass (update_bandwidths!) — time + allocations")
for r in 1:REPS
    tm = @timed update_bandwidths!(st, G, 1)
    @printf("  refresh %d: %6.2f s   %8.3f GiB alloc   gc %4.1f%%\n",
            r, tm.time, tm.bytes / 2^30, 100 * tm.gctime / max(tm.time, 1e-9))
end
@printf("  → ETEX projection: ×24 refreshes = the run's kNN share\n")

println("\n[2] adaptive deposit (_conc_nest! :adaptive) — time + allocations")
grid    = zeros(Float32, 400, 300, 2, 1, 1)
massarr = zeros(Float64, 400, 300, 2, 1, 1)
resnum  = zeros(Float64, 400, 300); resden = zeros(Float64, 400, 300)
_conc_nest!(grid, massarr, st, t0, 1.0f0, 400, 300, 1, 1, zi, 1,
            zlev, rho, LC, G, false, 600.0, Float32[1.0f0], :adaptive, resnum, resden)
for r in 1:REPS
    fill!(grid, 0.0f0); fill!(massarr, 0.0)
    tm = @timed _conc_nest!(grid, massarr, st, t0, 1.0f0, 400, 300, 1, 1, zi, 1,
                            zlev, rho, LC, G, false, 600.0, Float32[1.0f0],
                            :adaptive, resnum, resden)
    @printf("  deposit %d: %6.2f s   %8.3f GiB alloc   (144 calls in the run)\n",
            r, tm.time, tm.bytes / 2^30)
end

println("\n[3] uniform deposit (reference)")
for r in 1:REPS
    fill!(grid, 0.0f0); fill!(massarr, 0.0)
    tm = @timed _conc_nest!(grid, massarr, st, t0, 1.0f0, 400, 300, 1, 1, zi, 1,
                            zlev, rho, LC, G, false, 600.0, Float32[1.0f0], :uniform)
    @printf("  deposit %d: %6.2f s   %8.3f GiB alloc\n", r, tm.time, tm.bytes / 2^30)
end

println("\n[4] allocation sites (Profile.Allocs, kNN pass)")
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=RATE update_bandwidths!(st, G, 1)
res = Profile.Allocs.fetch()
agg = Dict{String,Tuple{Int,Int}}()          # site → (count, bytes)
typ = Dict{String,Tuple{Int,Int}}()          # alloc type → (count, bytes)
for a in res.allocs
    tk = string(a.type)
    c, b = get(typ, tk, (0, 0)); typ[tk] = (c + 1, b + a.size)
    site = "runtime/unattributed"
    for f in a.stacktrace
        s = string(f.file)
        if occursin("easlee", s) || occursin("conc.jl", s)
            site = basename(s) * ":" * string(f.line) * "  (" * string(f.func) * ")"
            break
        end
    end
    c, b = get(agg, site, (0, 0)); agg[site] = (c + 1, b + a.size)
end
println("  sampled $(length(res.allocs)) allocs (rate $RATE → ≈×$(round(Int,1/RATE)) real)")
if isempty(agg)
    println("  (no allocations sampled — the pass is allocation-free at this rate)")
    set_output_popmap!(Int[])
    println("\ndone")
    exit(0)
end
tot = sum(b for (_, b) in values(agg))
println("  top sites by sampled bytes:")
for (s, (c, b)) in sort(collect(agg); by=x -> -x[2][2])[1:min(12, length(agg))]
    @printf("   %6.1f%%  %9.2f MiB  n=%-8d %s\n", 100b / tot, b / 2^20, c, s)
end
println("  top allocation TYPES:")
for (s, (c, b)) in sort(collect(typ); by=x -> -x[2][2])[1:min(8, length(typ))]
    @printf("   %6.1f%%  %9.2f MiB  n=%-8d %s\n", 100b / tot, b / 2^20, c, s)
end
set_output_popmap!(Int[])
println("\ndone")
