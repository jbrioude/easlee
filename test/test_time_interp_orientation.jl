# test_time_interp_orientation.jl — pin the TIME orientation of every met sampler.
#
# Contract (slot space): `indt` from MetAdapter.met_time_index_array is defined
# so that indt = 1 means "at the valid time of the field in MEMORY SLOT 1" and
# indt = 2 "at slot 2's valid time", REGARDLESS of which slot holds the earlier
# file (readmetfiles! alternates the slot that receives new data; the fractional
# index formula (currentdate - t1)/(t2 - t1) + 1 with t1 = filedate[t_index[1]]
# handles both orientations through a negative dt).  Every sampler must therefore
# weight slot l by (1 - frac(indt)) and slot l+1 by frac(indt).
#
# History: the 4-D c-table (interpolate4D) had the time axis inverted — slot l
# weighted by frac instead of (1 - frac) — so all winds (interpolateU/V/W) were
# time-MIRRORED inside each met bracket: at an analysis time the parcel got 100%
# of the OTHER bracket end.  Zero error at mid-bracket, maximal at the ends.
# The profile family, interpolate2Dsingle and interpolate3Dirreg were correct,
# i.e. winds and coefficient profiles disagreed in time.  This test pins all of
# them to the same contract.

const REPO = abspath(joinpath(@__DIR__, ".."))

using Dates

include(joinpath(REPO, "easlee", "src", "met", "interp.jl"))
# turbulence_hanna.jl needs the Hanna dispatch tag (types/grid.jl defines it in
# the full model); a structurally identical local definition is enough here.
struct Hanna
    turbswitch::Bool
    sigw_ft::Float32
end
include(joinpath(REPO, "easlee", "src", "physics", "turbulence_hanna.jl"))  # interpolate2Dsingle
include(joinpath(REPO, "easlee", "src", "core", "metadapter.jl"))
using .MetAdapter: met_time_index_array

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name)
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

# ── Synthetic two-slot fields: value == 100 in slot 1, 200 in slot 2 ─────────
# Any temporally linear sampler must return 100 + 100*(indt-1).
N = 6
f4 = Array{Float32,4}(undef, N, N, N, 2); f4[:,:,:,1] .= 100f0; f4[:,:,:,2] .= 200f0
f3 = Array{Float32,3}(undef, N, N, 2);    f3[:,:,1]   .= 100f0; f3[:,:,2]   .= 200f0

expected(t) = 100f0 + 100f0 * (t - 1f0)
ts = Float32[1.0, 1.25, 1.5, 1.75, 2.0]

x = 2.5f0; y = 2.5f0; z = 2.5f0   # interior point, away from stagger edges

for t in ts
    e = expected(t)
    check("interpolation4Dsingle    t=$t",
          isapprox(interpolation4Dsingle(f4, x, y, z, t), e; atol = 1e-3),
          info = "got $(interpolation4Dsingle(f4, x, y, z, t)) want $e")
    check("interpolation4Dlevsingle t=$t",
          isapprox(interpolation4Dlevsingle(f4, x, y, z, t), e; atol = 1e-3),
          info = "got $(interpolation4Dlevsingle(f4, x, y, z, t)) want $e")
    check("interpolateU (edge)      t=$t",
          isapprox(interpolateU(true, f4, x, y, z, t), e; atol = 1e-3),
          info = "got $(interpolateU(true, f4, x, y, z, t)) want $e")
    check("interpolateV (edge)      t=$t",
          isapprox(interpolateV(true, f4, x, y, z, t), e; atol = 1e-3),
          info = "got $(interpolateV(true, f4, x, y, z, t)) want $e")
    check("interpolateW (edge)      t=$t",
          isapprox(interpolateW(true, f4, x, y, z, t), e; atol = 1e-3),
          info = "got $(interpolateW(true, f4, x, y, z, t)) want $e")
    check("interpolate2Dsingle      t=$t",
          isapprox(interpolate2Dsingle(f3, x, y, t), e; atol = 1e-3),
          info = "got $(interpolate2Dsingle(f3, x, y, t)) want $e")
    # trilinear table used with time as its 3rd axis (interpolate3Dirreg route)
    dl, l = _split_t(t)
    check("interpolate3Dirreg       t=$t",
          isapprox(interpolate3Dirreg(f4, 0.5f0, 0.5f0, dl, 2, 2, l, 3), e; atol = 1e-3),
          info = "got $(interpolate3Dirreg(f4, 0.5f0, 0.5f0, dl, 2, 2, l, 3)) want $e")
    # profile family
    buf = zeros(Float32, N)
    get_vertical_profile!(buf, f4, x, y, t, N)
    check("get_vertical_profile!    t=$t", isapprox(buf[3], e; atol = 1e-3),
          info = "got $(buf[3]) want $e")
    get_vertical_profile_nearest!(buf, f4, x, y, t, N)
    check("profile_nearest!         t=$t", isapprox(buf[3], e; atol = 1e-3),
          info = "got $(buf[3]) want $e")
    get_vertical_profile_window!(buf, f4, x, y, t, 3, 1, N)
    check("profile_window!          t=$t", isapprox(buf[3], e; atol = 1e-3),
          info = "got $(buf[3]) want $e")
end

# ── Level-dependent slot values: catch axis cross-talk (z vs t) ──────────────
# slot 1: value = k;   slot 2: value = k + 10.  At indt=1.5, level k=3.0 exactly
# → 3 + 5 = 8.  A z/t axis swap would break this even if the constant-field
# tests pass.
g4 = Array{Float32,4}(undef, N, N, N, 2)
for k in 1:N
    g4[:,:,k,1] .= Float32(k)
    g4[:,:,k,2] .= Float32(k) + 10f0
end
check("levsingle z/t independence",
      isapprox(interpolation4Dlevsingle(g4, x, y, 3.0f0, 1.5f0), 8.0f0; atol = 1e-3),
      info = "got $(interpolation4Dlevsingle(g4, x, y, 3.0f0, 1.5f0)) want 8.0")

# ── MetAdapter slot-space contract, both slot orientations ───────────────────
struct FakeNRC; t_index::NTuple{2,Int}; end
struct FakeMC;  filedate::Vector{DateTime}; end
struct FakeNest; nestrunconfig::FakeNRC; metconfig::FakeMC; end

t0 = DateTime(1994, 10, 23, 12)
fd = [t0, t0 + Hour(6), t0 + Hour(12)]

# fresh bracket: slot 1 = file 1 (earlier), slot 2 = file 2 (later)
fresh = [FakeNest(FakeNRC((1, 2)), FakeMC(fd))]
check("metadapter fresh @t1",   met_time_index_array(t0,           fresh)[1] == 1.0f0)
check("metadapter fresh @mid",  met_time_index_array(t0 + Hour(3), fresh)[1] == 1.5f0)
check("metadapter fresh @t2",   met_time_index_array(t0 + Hour(6), fresh)[1] == 2.0f0)

# rotated bracket: slot 1 = file 3 (LATER), slot 2 = file 2 (earlier)
rot = [FakeNest(FakeNRC((3, 2)), FakeMC(fd))]
check("metadapter rotated @slot1-time", met_time_index_array(t0 + Hour(12), rot)[1] == 1.0f0)
check("metadapter rotated @mid",        met_time_index_array(t0 + Hour(9),  rot)[1] == 1.5f0)
check("metadapter rotated @slot2-time", met_time_index_array(t0 + Hour(6),  rot)[1] == 2.0f0)

println()
println(npass, " passed, ", nfail, " failed")
nfail == 0 || println("FAIL")
