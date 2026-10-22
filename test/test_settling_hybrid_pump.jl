#!/usr/bin/env julia
# test_settling_hybrid_pump.jl — ISOLATION TEST for the ASCENT PUMP.
#
# Diagnosed 2026-08-04 at run scale: the population that cycles
# freeze -> settle -> sublimate MOVES UP (~300 m/h; water 20-26 km -> 32-40 km
# in ~2 h) under [physics.vertical] mode = "hybrid".  Removing EITHER settling
# or the converts removes the ascent (runs ht_asc4/ht_asc5).  Dynamics are
# innocent: SO2 (never freezes) stays pinned; ERA5 etadot is 0 at the vent.
#
# This test reproduces the cycle for ONE parcel on a realistic stratospheric
# theta ladder (measured ERA5 vent column, 2022-01-15 18:00: theta 390 K at
# 17.6 km -> 1250 K at 38.5 km, dtheta/dz 20-60 K/km) and prints z AND theta_p
# every step.  Settling is DOWNWARD only, so:
#
#   THE ASSERTION: z must be monotone NON-INCREASING across the cycle.
#
# A step where theta_p jumps by tens-to-hundreds of K localises the writer.
# Suspects, in order: the settling rebase (species_physics.jl), the hybrid
# blend-back (isentropic.jl:341), and the two SILENT fallbacks _z_to_indl /
# _theta_to_indl, which both return nz-1 (the MODEL TOP) when their argument
# leaves the ladder.
#
#   julia --project=easlee test/test_settling_hybrid_pump.jl
using Printf, Dates

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: apply_convert_thermo!, apply_species_settling!,
                         esat_ice, interpolation4Dsingle

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic met: REALISTIC stratospheric column (ERA5 at the HT vent) ──────
# 21 levels, 1 km apart, 16 -> 36 km, on the stacked layout contract: zlev[k]
# = zk(k) are the half levels (class H); the full-level fields theta, rho, p
# (class F) hold the column base in slot 1 and full level k at zk(k) - DZ/2 in
# slot k >= 2.  T and p interpolated from the measured profile; theta derived,
# so the ladder is the real one the model inverts.
const NZ  = 21
const Z0  = 16_000.0
const DZ  = 1_000.0
zk(k) = Z0 + (k - 1) * DZ
# measured (z_km, T_K, p_hPa) at the vent, 2022-01-15 18:00
const ZM = [15.3, 16.4, 17.6, 18.9, 20.3, 21.9, 23.7, 25.9, 28.3, 31.1, 34.4, 38.5]
const TM = [200.0, 194.0, 192.2, 194.5, 202.1, 212.7, 215.2, 224.0, 226.8, 228.8, 240.6, 247.0]
const PM = [127.33, 103.71, 83.88, 67.19, 52.86, 40.53, 30.18, 21.71, 14.98, 9.82, 6.04, 3.43]
function interp1(xs, ys, x)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs, x)
    f = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + f * (ys[k+1] - ys[k])
end
dims  = (8, 8, NZ, 2)
zlev  = zeros(Float32, dims); theta = zeros(Float32, dims)
rho   = zeros(Float32, dims); puv   = zeros(Float32, dims)
Tcol  = zeros(NZ); pcol = zeros(NZ); thcol = zeros(NZ)
for k in 1:NZ
    z = k == 1 ? zk(1) : zk(k) - DZ / 2          # full level k
    T = interp1(ZM, TM, z/1000); p = interp1(ZM, PM, z/1000) * 100.0
    Tcol[k] = T; pcol[k] = p; thcol[k] = T * (1.0e5 / p)^0.2857
    zlev[:, :, k, :]  .= Float32(zk(k))
    theta[:, :, k, :] .= Float32(thcol[k])
    rho[:, :, k, :]   .= Float32(p / (287.05 * T))
    puv[:, :, k, :]   .= Float32(p)
end
m3 = (theta = theta, rho = rho, zlev = zlev, puv = puv, T = Float32.(theta) .* 0 .+ 220f0)
LC = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
GC = nothing
fctx = (nests = [(metfields3D = m3, levelconfig = LC,)],)

@printf("theta ladder: %.0f K at %.0f km -> %.0f K at %.0f km  (dtheta/dz %.0f-%.0f K/km)\n",
        thcol[1], zk(1)/1000, thcol[NZ], zk(NZ)/1000,
        (thcol[2]-thcol[1])/1, (thcol[NZ]-thcol[NZ-1])/1)

# ── the cycle: ICE box 1 (settles, sublimates -> box 2 = vapour) ─────────────
const NB = 2
const R_ICE = 4.24e-6                       # bin b3 of the HT ladder
function mkcv(; nslots = 4)
    EASLEE.SpeciesConvert(Float32[1.0, 0.0], Int32[2, 0], Float32[1.0, 1.0],
                          zeros(Float64, NB), Int8[2, 0], Float32[1.5, 1.5],
                          zeros(Float32, NB), Float32[R_ICE, 0.0],
                          zeros(Float32, nslots), Int32[0, 0],
                          zeros(Float32, NB))
end
function mksp(; settle = true)
    z2 = zeros(Float32, NB)
    EASLEE.SpeciesPhysics(z2, copy(z2), copy(z2), zeros(Bool, NB), copy(z2),
                          copy(z2), copy(z2), copy(z2), 0f0,
                          [zeros(Float64, 8, 8)],
                          settle ? Float32[R_ICE, 0.0] : copy(z2),
                          Float32[917.0, 917.0])
end

# start at 22 km (level index for 22 km on this ladder), saturated so the ice
# does NOT vanish instantly (S just under 1 keeps the r² law gentle)
const K_START = (22_000.0 - Z0) / DZ + 1
esi0 = esat_ice(interp1(ZM, TM, 22.0))
qsi0 = 0.622 * esi0 / (interp1(ZM, PM, 22.0)*100 - esi0)
function mkstate(; nslots = 4, S = 0.9999)
    st = EASLEE.State.ParcelState(nslots)
    st.status[1] = 0x01; st.nest[1] = Int32(1); st.relbox[1] = Int32(1)
    st.indx[1] = 4.0f0; st.indy[1] = 4.0f0; st.indl[1] = Float32(K_START)
    st.mass[1] = 1.0f6; st.airm[1] = Float32(1.0e6 / (S * qsi0))
    return st
end

zof(st) = begin
    k = clamp(Int(floor(st.indl[1])), 1, NZ-1); f = st.indl[1] - k
    (zk(k) + f * DZ) / 1000.0
end

println("="^78)
println("[1] settling ALONE (no hybrid step): z must DECREASE")
println("="^78)
st = mkstate(); sp = mksp(); cv = mkcv()
θ_start = EASLEE.Operators.interpolation4Dsingle(theta, 4.0f0, 4.0f0, Float32(K_START), 1.0f0)
θp = zeros(Float32, 4); θp[1] = θ_start
z0 = zof(st)
for s in 1:5
    apply_species_settling!(st, fctx, sp, 300.0, [1.0]; thetap = θp, rice = cv.rice)
end
@printf("  z: %.3f -> %.3f km   theta_p: %.1f -> %.1f K\n", z0, zof(st), θ_start, θp[1])
check("settling alone moves the parcel DOWN", zof(st) <= z0 + 1e-6,
      info=@sprintf("dz = %+.3f km", zof(st) - z0))
const FALL_STEP = (z0 - zof(st)) / 5           # km per 300 s step, settling alone

println("="^78)
println("[2] THE CYCLE: settle + convert + hybrid, 24 steps of 300 s")
println("="^78)
# σ_r = 0.9, σ_t = 0.8: σ = p/p(16 km) ≈ 0.36 at 22 km, so w = 1 — the
# hybrid re-imposes the parcel's isentrope after EVERY settling step, which is
# where a rebase/ladder mismatch turns into motion (with the operational
# 0.3/0.1 this column sits at w = 0 and the hybrid never relocates anything).
hyb = EASLEE.HybridVert(Union{EASLEE.Operators.IsentropicFields,Nothing}[nothing],
                        nothing, zeros(Float32, 4), 0.9f0, 0.8f0, false)
# S = 3 (unambiguously supersaturated): the sublimate gate's eligible fraction
# is 0, so the parcel stays ICE and settles on EVERY step — precisely the
# pump's duty cycle.  S is set with margin because qsi0 is computed at 22 km
# from the measured profile while the gate reads the stacked grid T (a nominal
# S = 1.0 is not S = 1.0 as the gate sees it).  With S < 1 the crystal
# vanishes in a single step at these sizes and the test degenerates to one
# settling call.
st = mkstate(S = 3.0); sp = mksp(); cv = mkcv()
hyb.theta_p[1] = 0.0f0                       # unset -> initialised from the grid
zs = Float64[]; ths = Float64[]
println("  step   z_km    theta_p_K   relbox   rice_um    d(theta_p)")
global thprev = NaN
for s in 1:24
    apply_species_settling!(st, fctx, sp, 300.0, [1.0];
                            thetap = hyb.theta_p, rice = cv.rice)
    apply_convert_thermo!(st, fctx, cv, 300.0, [1.0])
    EASLEE.Operators.hybrid_vertical!(hyb, st,
        (nests = [(metfields3D = m3, levelconfig = LC,
                   gridconfig = (Ni = Int32(8), Nj = Int32(8)))],
         runconfig = (direction = 1, currentdate = DateTime(2022,1,15,18))),
        300.0, [1.0])
    z = zof(st); th = Float64(hyb.theta_p[1])
    push!(zs, z); push!(ths, th)
    dth = isnan(thprev) ? 0.0 : th - thprev
    (s <= 8 || abs(dth) > 20.0 || s == 24) &&
        @printf("  %4d  %6.3f   %8.1f   %6d   %7.3f   %+9.2f\n",
                s, z, th, Int(st.relbox[1]), cv.rice[1]*1e6, dth)
    global thprev = th
end
@printf("\n  NET over 24 steps (2 h): z %.3f -> %.3f km  (%+.3f km)\n", zs[1], zs[end], zs[end]-zs[1])
check("z is monotone NON-INCREASING through the cycle",
      all(diff(zs) .<= 1e-6), info=@sprintf("max rise in one step = %+.4f km", maximum(diff(zs))))
check("net displacement is DOWNWARD", zs[end] <= zs[1] + 1e-6,
      info=@sprintf("%+.3f km in 2 h", zs[end]-zs[1]))
check("the parcel stayed ICE for all 24 steps (settling was exercised)",
      Int(st.relbox[1]) == 1,
      info="relbox = $(Int(st.relbox[1])) (1 = ice, 2 = vapour)")
check("the fall is the physical one (v_t·dt·24, not zero)",
      (zs[1] - zs[end]) > 0.01,
      info=@sprintf("fell %.3f km in 2 h", zs[1]-zs[end]))
# Under the stacked contract a rebase read on the PLAIN convention sits half a
# level LOW, so the mismatch now pumps DOWN (measured with that read restored:
# -4.9 km in 2 h) and the monotone checks above cannot see it.  Bound the fall
# by the settling-alone rate of [1] (23 intervals between zs[1] and zs[end]).
check("the fall is not amplified by the hybrid (≤ 2× settling alone)",
      (zs[1] - zs[end]) < 2 * 23 * FALL_STEP,
      info=@sprintf("fell %.3f km vs %.3f km from settling alone",
                    zs[1]-zs[end], 23 * FALL_STEP))

println("="^78)
println("[3] THE INVARIANT that was violated: theta_p must be what the hybrid inverts")
println("="^78)
# The 2026-08-05 pump was a half-level mismatch between the theta written by the
# settling rebase and the theta ladder the hybrid inverts.  theta is a
# full-level field: under the stacked layout contract (2026-10-22) both read it
# with the stacked map (slot k at indl k - 1/2), as interpolation4Dsingle does.
# (Updated 2026-10-22: this block used to assert the co-located "plain" read,
# theta[k] at indl k, and flag the stacked read as the bug.)  The guard is the
# round trip: a w = 1 hybrid step right after a settling rebase, with no
# diabatic term, must leave the parcel exactly where settling put it.
st3 = mkstate(); sp3 = mksp(); cv3 = mkcv()
θ3 = zeros(Float32, 4); θ3[1] = 1.0f0
apply_species_settling!(st3, fctx, sp3, 300.0, [1.0]; thetap = θ3, rice = cv3.rice)
indl3 = st3.indl[1]
th_stack = Float64(EASLEE.Operators.interpolation4Dsingle(theta, 4.0f0, 4.0f0, indl3, 1.0f0))
k0 = clamp(Int(Base.unsafe_trunc(Int, indl3)), 1, NZ-1)
th_plain = thcol[k0] + (indl3 - k0) * (thcol[k0+1] - thcol[k0])   # pre-contract read
θ3_set = θ3[1]
hyb3 = EASLEE.HybridVert(Union{EASLEE.Operators.IsentropicFields,Nothing}[nothing],
                         nothing, θ3, 0.9f0, 0.8f0, false)
EASLEE.Operators.hybrid_vertical!(hyb3, st3,
    (nests = [(metfields3D = m3, levelconfig = LC,
               gridconfig = (Ni = Int32(8), Nj = Int32(8)))],
     runconfig = (direction = 1, currentdate = DateTime(2022,1,15,18))),
    300.0, [1.0])
@printf("  indl after settling = %.4f, after the w = 1 hybrid step = %.4f\n", indl3, st3.indl[1])
@printf("  theta_p written by settling = %.2f K\n", θ3_set)
@printf("  stacked grid theta (CONTRACT) = %.2f K   <- must match\n", th_stack)
@printf("  plain read (pre-contract)     = %.2f K   (offset %+.2f K)\n",
        th_plain, th_plain - th_stack)
check("settling's theta_p == the stacked grid theta",
      isapprox(Float64(θ3_set), th_stack; rtol=1e-4),
      info=@sprintf("%.2f vs %.2f K", θ3_set, th_stack))
check("a w = 1 hybrid step does not move the settled parcel",
      abs(st3.indl[1] - indl3) < 1e-4,
      info=@sprintf("Δindl = %+.2e", st3.indl[1] - indl3))
# What matters for the guard is that a read on the other convention is wrong by
# MANY KELVIN — that is what the hybrid would turn into a vertical displacement.
check("the plain read is wrong by many K (the mismatch's signature)",
      abs(th_plain - th_stack) > 5.0,
      info=@sprintf("%+.1f K = %.0f m per settling step",
                    th_plain - th_stack,
                    1000 * (th_plain - th_stack) / ((thcol[k0+1]-thcol[k0]) / (DZ/1000))))

@printf("\n%d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
