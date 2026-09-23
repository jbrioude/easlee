# test_substeps_clamp.jl — pins the CPU hot-path fixes of the 2026-09-04 review
# (physics/advance.jl, met/interp.jl, output/conc.jl):
#   A1  staggered interpolation at the model top: indl in [nz-0.5, nz) must
#       return the top mass level and never touch field[…, nz+1, …] (which
#       aliased the next time slot's surface under @inbounds).
#   A2  the CTL sub-step loop runs EXACTLY ctl sub-steps (the Float32
#       `while inc_time < Δt` accumulator ran ctl+1 for ~20 % of (Δt, ctl)).
#   A3  the Monin-Obukhov guard `indl < 2` is equivalent to `z_sl < z1_sl`.
#   A5  the grid deposit skips culled parcels (status 0x00, Dtime == now).
# Synthetic data only — no met archive needed.

const REPO = abspath(joinpath(@__DIR__, ".."))

using Dates

include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: _nsubsteps, _conc_nest!
# NOTE: io/tools.jl defines a legacy `interpolation4Dsingle` in Main — always
# call the engine's kernels through the module below.
const OP = EASLEE.Operators

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name)
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

# ═════════════════════════════════════════════════════════════════════════════
println("[A2] exact integer sub-step count")
# ═════════════════════════════════════════════════════════════════════════════
# Replica of the former loop, for the record.
function old_loop_count(Δt_s::Real, ctl::Real)
    float_timestep::Float32 = Float32(abs(Δt_s))
    dt_inc::Float32 = float_timestep / Float32(ctl)
    inc_time::Float32 = 0.0f0
    n = 0
    while inc_time < float_timestep
        inc_time += dt_inc
        n += 1
    end
    return n
end
function new_loop_count(Δt_s::Real, ctl::Real)
    nsub, dt = _nsubsteps(Float32(abs(Δt_s)), Float32(ctl))
    n = 0
    for _ in 1:nsub; n += 1; end
    return n, dt
end
# The worst previous cases (verified ctl+1 under the old loop) + a few sane ones.
cases = [(3, 9), (1, 12), (900, 7), (3600, 7), (600, 5), (60, 1), (1800, 10)]
for (Δt, ctl) in cases
    nold = old_loop_count(Δt, ctl)
    nnew, dt = new_loop_count(Δt, ctl)
    println("      Δt=$(rpad(Δt, 5)) ctl=$(rpad(ctl, 3)) old=$(nold)  new=$(nnew)  dt=$(dt)")
    check("Δt=$(Δt)s ctl=$(ctl): new loop runs exactly ctl sub-steps", nnew == ctl,
          info = "got $nnew")
    check("Δt=$(Δt)s ctl=$(ctl): nsub·dt == Δt (Float32)",
          isapprox(nnew * dt, Float32(Δt); rtol = 2 * eps(Float32)),
          info = "nsub·dt = $(nnew * dt)")
end
for (Δt, ctl) in ((3, 9), (1, 12), (900, 7), (3600, 7))
    check("Δt=$(Δt)s ctl=$(ctl): old loop DID run ctl+1 (defect reproduced)",
          old_loop_count(Δt, ctl) == ctl + 1, info = "old = $(old_loop_count(Δt, ctl))")
end
# Exhaustive sweep: the new count is ctl for every (Δt, ctl) the old one missed.
nbad_old = 0; nbad_new = 0
for Δt in 1:3600, ctl in 1:20
    old_loop_count(Δt, ctl) == ctl || (global nbad_old += 1)
    new_loop_count(Δt, ctl)[1] == ctl || (global nbad_new += 1)
end
println("      sweep Δt=1..3600 s × ctl=1..20: old loop wrong in $nbad_old / 72000, new in $nbad_new")
check("sweep: new loop never deviates from ctl", nbad_new == 0)
check("sweep: old loop was wrong for a sizeable fraction (>10 %)", nbad_old > 7200)
check("ctl given as Float32 integer round-trips", _nsubsteps(600.0f0, 5.0f0) == (5, 120.0f0))
check("ctl < 1 is floored at one sub-step", _nsubsteps(600.0f0, 0.0f0) == (1, 600.0f0))

# ═════════════════════════════════════════════════════════════════════════════
println("[A1] staggered interpolation at the model top")
# ═════════════════════════════════════════════════════════════════════════════
# A checked wrapper: getindex ALWAYS bounds-checks (not under @boundscheck), so
# an @inbounds read past size(field, 3) throws instead of aliasing.
struct Checked{T,N} <: AbstractArray{T,N}
    a::Array{T,N}
end
Base.size(c::Checked) = size(c.a)
Base.IndexStyle(::Type{<:Checked}) = IndexCartesian()
Base.@propagate_inbounds function Base.getindex(c::Checked{T,N}, I::Vararg{Int,N}) where {T,N}
    checkbounds(c.a, I...)            # unconditional — @inbounds cannot elide it
    return c.a[I...]
end

const NX = 6; const NZT = 7
# slot 1: value = mass-level index k; slot 2: k + 100 (so time weighting shows).
f4 = Array{Float32,4}(undef, NX, NX, NZT, 2)
for k in 1:NZT; f4[:, :, k, 1] .= Float32(k); f4[:, :, k, 2] .= Float32(k + 100); end
# Sentinel variant: the level that field[…, NZT+1, 1] aliases under linear
# indexing is field[…, 1, 2] (slot 2's surface).  Poison the WHOLE slot-2
# surface so any leak is unmistakable.
f4s = copy(f4); f4s[:, :, 1, 2] .= 1.0f30

x = 2.5f0; y = 3.25f0
indl_top = Float32(NZT) - 1.0f-3          # the parcel clamp used in advance.jl

# (i) no leak of the aliased slot: value == top mass level at t = 1 (dl = 0)
v = OP.interpolation4Dsingle(f4s, x, y, indl_top, 1.0f0)
check("4D stagger @ indl = nz-1e-3, t=1: finite", isfinite(v), info = "got $v")
check("4D stagger @ indl = nz-1e-3, t=1: == top mass level (no slot-2 surface leak)",
      isapprox(v, Float32(NZT); atol = 2 * eps(Float32) * NZT), info = "got $v want $NZT")
# (ii) time weighting still applies at the clamped level
for t in (1.0f0, 1.25f0, 1.5f0, 2.0f0)
    want = Float32(NZT) + 100f0 * (t - 1f0)
    got  = OP.interpolation4Dsingle(f4, x, y, indl_top, t)
    check("4D stagger @ indl = nz-1e-3, t=$t: == top level, time-weighted",
          isapprox(got, want; atol = 1e-3), info = "got $got want $want")
end
# (iii) the whole band [nz-0.5, nz) is the top-level plateau, continuous with
#       the regular branch at nz-0.5-ε
for z in (Float32(NZT) - 0.5f0, Float32(NZT) - 0.25f0, indl_top)
    got = OP.interpolation4Dsingle(f4s, x, y, z, 1.0f0)
    check("4D stagger @ indl=$z: top-level plateau, no leak",
          isapprox(got, Float32(NZT); atol = 1e-4), info = "got $got")
end
zm = Float32(NZT) - 0.5f0 - 1.0f-3
got = OP.interpolation4Dsingle(f4s, x, y, zm, 1.0f0)
check("4D stagger just below nz-0.5: continuous with the plateau",
      isapprox(got, Float32(NZT); atol = 5e-3), info = "got $got")
# (iv) interior points unchanged by the clamp
for (z, want) in ((2.0f0, 2.5f0), (3.5f0, 4.0f0), (1.0f0, 1.0f0), (1.25f0, 1.5f0))
    g = OP.interpolation4Dsingle(f4, x, y, z, 1.0f0)
    check("4D stagger interior indl=$z unchanged (== $want)",
          isapprox(g, want; atol = 1e-4), info = "got $g")
end
# (vi) the defect, for the record: the raw stencil with the UNclamped stagger
#      index reads field[…, nz+1, 1] == field[…, 1, 2] (the poisoned surface)
#      (a third time slot keeps the l+1 = 2, k+1 = nz+1 read inside the
#      allocation — it aliases slot 3's surface — so the replica is not UB).
let (dzr, kkr) = OP._split_z_stagger(indl_top)
    f4s3 = cat(f4s, f4s[:, :, :, 2:2]; dims = 4); f4s3[:, :, 1, 3] .= 0.0f0
    raw = OP.interpolate4D(f4s3, Float32(mod(x, 1)), Float32(mod(y, 1)), dzr, 0.0f0,
                           Base.unsafe_trunc(Int, x), Base.unsafe_trunc(Int, y), kkr, 1)
    println("      unclamped stencil (kk=$kkr, dk=$dzr) on the poisoned array → $raw")
    check("defect reproduced: unclamped stencil leaks the slot-2 surface sentinel",
          kkr == NZT && raw > 1.0f29, info = "got $raw")
end
# (v) never indexes beyond size(field, 3) — checked wrapper throws on OOB
function no_oob(fn)
    try
        fn(); return true
    catch e
        e isa BoundsError && return false
        rethrow()
    end
end
c4 = Checked(f4)
check("4D stagger @ indl = nz-1e-3, t=1: no read beyond size(field,3)",
      no_oob(() -> OP.interpolation4Dsingle(c4, x, y, indl_top, 1.0f0)))
check("4D stagger @ indl = nz-1e-3, t=2 (l=1, dl=1): no read beyond the array",
      no_oob(() -> OP.interpolation4Dsingle(c4, x, y, indl_top, 2.0f0)))
check("4D stagger @ indl = nz-0.5: no read beyond size(field,3)",
      no_oob(() -> OP.interpolation4Dsingle(c4, x, y, Float32(NZT) - 0.5f0, 1.0f0)))
# 3-D variant (x, y, z) and the shared index decomposition
f3 = Array{Float32,3}(undef, NX, NX, NZT)
for k in 1:NZT; f3[:, :, k] .= Float32(k); end
check("3D stagger @ indl = nz-1e-3: == top mass level",
      isapprox(OP.interpolation3Dsingle(f3, x, y, indl_top), Float32(NZT); atol = 1e-4))
check("3D stagger @ indl = nz-1e-3: no read beyond size(field,3)",
      no_oob(() -> OP.interpolation3Dsingle(Checked(f3), x, y, indl_top)))
di, dj, dz, dl, ii, jj, kk, l = OP._parcel_indices(x, y, indl_top, 1.0f0, NZT)
check("_parcel_indices(…, nk): kk clamped to nk-1 with dz = 1",
      kk == NZT - 1 && dz == 1.0f0, info = "kk=$kk dz=$dz")
check("_parcel_indices(…, nk) + interp4: == top mass level, no leak",
      isapprox(OP.interp4(f4s, di, dj, dz, dl, ii, jj, kk, l), Float32(NZT); atol = 1e-4))
# interpolation4Dlevsingle (unstaggered, already guarded) is untouched
check("4Dlevsingle @ k = nz exactly: still == top level",
      isapprox(OP.interpolation4Dlevsingle(f4s, x, y, Float32(NZT), 1.0f0), Float32(NZT); atol = 1e-4))

# ═════════════════════════════════════════════════════════════════════════════
println("[A3] Monin-Obukhov guard: indl < 2  ⇔  z_sl < z1_sl")
# ═════════════════════════════════════════════════════════════════════════════
zlev = Array{Float32,4}(undef, NX, NX, NZT, 2)
for k in 1:NZT
    zlev[:, :, k, 1] .= Float32(20.0 * (k - 1)^1.3)          # monotone, irregular
    zlev[:, :, k, 2] .= Float32(24.0 * (k - 1)^1.3)          # different column at slot 2
end
# add a gentle horizontal gradient (terrain) to make the bilinear part non-trivial
for i in 1:NX, j in 1:NX; zlev[i, j, :, :] .+= Float32(0.7 * i - 0.3 * j); end
nmismatch = 0; ntested = 0
for t in (1.0f0, 1.4f0, 2.0f0), xi in (1.0f0, 2.3f0, 4.9f0), yj in (1.5f0, 3.8f0),
    indl in Float32[1.0, 1.001, 1.5, 1.999, 2.0, 2.001, 2.5, 3.0, 4.7, NZT - 0.5, NZT - 1e-3]
    z_sl  = OP.interpolation4Dlevsingle(zlev, xi, yj, indl, t)
    z1_sl = OP.interpolation4Dlevsingle(zlev, xi, yj, 2.0f0, t)
    global ntested += 1
    (isfinite(z_sl) && z_sl < z1_sl) == (indl < 2.0f0) || (global nmismatch += 1)
end
println("      $ntested (t, x, y, indl) combinations, $nmismatch mismatches")
check("guard `indl < 2` ⇔ former `z_sl < z1_sl` on every tested position", nmismatch == 0)

# ═════════════════════════════════════════════════════════════════════════════
println("[A5] grid deposit skips culled parcels (status 0x00, Dtime == now)")
# ═════════════════════════════════════════════════════════════════════════════
const NZ = 10
dims = (12, 10, NZ, 2)
zl = zeros(Float32, dims); rho = fill(1.0f0, dims)
for k in 1:NZ; zl[:, :, k, :] .= Float32((k - 1) * 30.0); end
GEO = EASLEE.Types.Rectilinear(Int32(12), Int32(10), Int32(NZ), 1.0f0, 1.0f0,
        0.0f0, -50.0f0, Float32(1/111320.0), Float32(1/111320.0),
        1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, false)
LC  = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
t0  = DateTime(2025, 1, 1)
st = EASLEE.State.ParcelState(4)
for i in 1:3
    st.status[i] = 0x01; st.nest[i] = Int32(1); st.relbox[i] = Int32(1)
    st.indx[i] = 5.0f0; st.indy[i] = 5.0f0; st.indl[i] = 3.0f0
    st.mass[i] = 2.0f0; st.Dtime[i] = t0
end
st.status[2] = 0x00      # culled this step (domainfill_cull!): Dtime still == t0
st.status[3] = 0x02      # deposited/inactive flag — also not active
zi = Float32[0.0, 100.0, 200.0]
g5 = zeros(Float32, 12, 10, 3, 1, 1); m5 = zeros(Float64, size(g5))
d = _conc_nest!(g5, m5, st, t0, 1.0f0, 12, 10, 1, 1, zi, 1,
                zl, rho, LC, GEO, false, 60.0, Float32[1.0f0], :none)
println("      deposited parcels = $(d.n_deposited), mass = $(sum(m5)) kg (active parcel carries 2 kg)")
check("only the ACTIVE parcel is deposited (1 parcel, 2 kg)",
      d.n_deposited == 1 && isapprox(sum(m5), 2.0; rtol = 1e-12),
      info = "n=$(d.n_deposited) mass=$(sum(m5))")

println("\n$npass passed, $nfail failed")
nfail == 0 || exit(1)
