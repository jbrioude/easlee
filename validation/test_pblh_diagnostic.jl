#!/usr/bin/env julia
# test_pblh_diagnostic.jl — validate the PBL-height diagnostics:
#   1. _pblh_from_tke! on a synthetic TKE profile (threshold = max(5 % sfc, 0.01))
#   2. _pblh_from_rib! (bulk Richardson, Ri_c = 0.25) against ERA5's own archived
#      blh over the Pacific domain — ERA5's blh is itself a Rib(0.25) diagnostic,
#      so the two should correlate strongly (r > 0.7) with a median ratio ≈ 1.
#
#   julia --project=easlee test/test_pblh_diagnostic.jl
using Dates, Printf, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

println("="^78)
println("[1] TKE-threshold PBL height (synthetic profile)")
println("="^78)
# 1-column dicts: TKE 0.5 m²/s² in the lowest 800 m decaying to ~0 above 1 km
Nk = 20
z  = Float32.(0:100:1900)
tke = Float32[zz <= 800 ? 0.5f0 : 0.5f0 * exp(-(zz - 800) / 80) for zz in z]
m3 = Dict{String,Any}("tke"  => reshape(tke, 1, 1, Nk, 1),
                      "zlev" => reshape(z,   1, 1, Nk, 1))
m2 = Dict{String,Any}("pblh" => zeros(Float32, 1, 1, 1))
Main._pblh_from_tke!(m3, m2, 1, 1, 1, Nk)
h = m2["pblh"][1, 1, 1]
# threshold 0.05·0.5 = 0.025 → crossed at 800 + 80·ln(0.5/0.025) ≈ 1040 m; the
# first LEVEL below threshold is the one at or above that → expect 1000–1200 m
check("TKE PBL top at the decay level", 900 <= h <= 1300; info=@sprintf("(h = %.0f m)", h))
tke0 = zeros(Float32, 1, 1, Nk, 1); m3["tke"] = tke0
Main._pblh_from_tke!(m3, m2, 1, 1, 1, Nk)
check("zero TKE → 50 m floor", m2["pblh"][1,1,1] == 50.0f0)

println()
println("="^78)
println("[2] Bulk-Richardson PBL vs the archived ERA5 blh (Pacific, 1 met step)")
println("="^78)
spec = read_runspec(joinpath(REPO, "configs", "era5_pacific.toml"))
for r in spec.release; r.particles = 10; end          # tiny — met fields only
rt, _ = build_run(spec; prefetch=false)
read_met_step!(rt, spec.start)
indt = EASLEE.MetAdapter.met_time_index_array(spec.start, rt.ctx.nests)[1]
it   = clamp(round(Int, indt), 1, 2)
m3d  = rt.metfields3D[1]; m2d = rt.metfields2D[1]
gc   = rt.ctx.nests[1].gridconfig
Ni, Nj, Nk2 = Int(gc.Ni), Int(gc.Nj), size(m3d["zlev"], 3)
blh_arch = copy(m2d["pblh"][:, :, it])                # ERA5's own blh
fill!(view(m2d["pblh"], :, :, it), 0.0f0)
Main._pblh_from_rib!(m3d, m2d, it, Ni, Nj, Nk2)
blh_rib = m2d["pblh"][:, :, it]
m2d["pblh"][:, :, it] .= blh_arch                     # restore

# compare over the Pacific box (150–230°E, 0–50°N)
pairs = Tuple{Float64,Float64}[]
for j in 1:Nj, i in 1:Ni
    lat, lon = EASLEE.Operators.indyindxtolatlon(gc, Float64(j), Float64(i))
    lonw = mod(lon - 150.0, 360.0) + 150.0
    (150.0 <= lonw <= 230.0 && 0.0 <= lat <= 50.0) || continue
    a, r = Float64(blh_arch[i, j]), Float64(blh_rib[i, j])
    (isfinite(a) && isfinite(r) && a > 10) || continue
    push!(pairs, (a, r))
end
A = [p[1] for p in pairs]; R = [p[2] for p in pairs]
r  = cor(A, R)
mr = median(R ./ A)
@printf("  n=%d cells   archived blh median %.0f m   Rib median %.0f m\n",
        length(A), median(A), median(R))
check("correlation vs ERA5 blh > 0.7", r > 0.7; info=@sprintf("(r = %.3f)", r))
check("median ratio Rib/ERA5 in [0.5, 2]", 0.5 < mr < 2.0; info=@sprintf("(ratio %.2f)", mr))

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
