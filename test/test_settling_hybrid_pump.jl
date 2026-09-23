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
# 21 levels, 1 km apart, 16 -> 36 km.  T and p interpolated from the measured
# profile; theta derived, so the ladder is the real one the model inverts.
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
    z = zk(k); T = interp1(ZM, TM, z/1000); p = interp1(ZM, PM, z/1000) * 100.0
    Tcol[k] = T; pcol[k] = p; thcol[k] = T * (1.0e5 / p)^0.2857
    zlev[:, :, k, :]  .= Float32(z)
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
θp = zeros(Float32, 4); θp[1] = Float32(thcol[Int(K_START)])
z0 = zof(st)
for s in 1:5
    apply_species_settling!(st, fctx, sp, 300.0, [1.0]; thetap = θp, rice = cv.rice)
end
@printf("  z: %.3f -> %.3f km   theta_p: %.1f -> %.1f K\n", z0, zof(st), thcol[Int(K_START)], θp[1])
check("settling alone moves the parcel DOWN", zof(st) <= z0 + 1e-6,
      info=@sprintf("dz = %+.3f km", zof(st) - z0))

println("="^78)
println("[2] THE CYCLE: settle + convert + hybrid, 24 steps of 300 s")
println("="^78)
hyb = EASLEE.HybridVert(Union{EASLEE.Operators.IsentropicFields,Nothing}[nothing],
                        nothing, zeros(Float32, 4), 0.3f0, 0.1f0, false)
# S = 3 (unambiguously supersaturated): the sublimate gate's eligible fraction
# is 0, so the parcel stays ICE and settles on EVERY step — precisely the
# pump's duty cycle.  S must be set with margin because the GATE reads T via
# the staggered map while this test's qsi is computed on the plain one, so a
# nominal S = 1.0 is not S = 1.0 as the gate sees it (see journal 08-05, the
# residual convention split).  With S < 1 the crystal vanishes in a single
# step at these sizes and the test degenerates to one settling call.
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

println("="^78)
println("[3] THE INVARIANT that was violated: theta_p must use the CONSUMER's convention")
println("="^78)
# The hybrid computes its grid theta as  th[k0] + (indl-k0)*(th[k0+1]-th[k0])
# (isentropic.jl).  Whatever settling writes into theta_p MUST agree with that,
# else the hybrid relocates the parcel.  `interpolation4Dsingle` applies the
# STAGGER map (indl -> mass level indl+0.5) and does NOT agree — that mismatch
# was the pump.  This assertion is the regression guard.
st3 = mkstate(); sp3 = mksp(); cv3 = mkcv()
θ3 = zeros(Float32, 4); θ3[1] = 1.0f0
apply_species_settling!(st3, fctx, sp3, 300.0, [1.0]; thetap = θ3, rice = cv3.rice)
indl3 = st3.indl[1]
k0 = clamp(Int(Base.unsafe_trunc(Int, indl3)), 1, NZ-1)
th_grid = thcol[k0] + (indl3 - k0) * (thcol[k0+1] - thcol[k0])   # hybrid's formula
th_stag = Float64(EASLEE.Operators.interpolation4Dsingle(theta, 4.0f0, 4.0f0, indl3, 1.0f0))
@printf("  indl = %.4f\n  theta_p written by settling = %.2f K\n", indl3, θ3[1])
@printf("  hybrid's th_grid (CONSUMER)  = %.2f K   <- must match\n", th_grid)
@printf("  staggered read (the old bug) = %.2f K   (offset %+.2f K = %.2f level)\n",
        th_stag, th_stag - th_grid, (th_stag - th_grid)/(thcol[k0+1]-thcol[k0]))
check("settling's theta_p == the hybrid's own grid theta",
      isapprox(Float64(θ3[1]), th_grid; rtol=1e-4),
      info=@sprintf("%.2f vs %.2f K", θ3[1], th_grid))
# The offset is ~0.5 level by construction, but theta is CURVED in k, so in
# local-gradient units it lands near 0.35 here.  What matters for the guard is
# that the wrong read is wrong by MANY KELVIN — that is what the hybrid turns
# into a vertical displacement.
check("the staggered read is wrong by many K (the bug's signature)",
      (th_stag - th_grid) > 5.0,
      info=@sprintf("+%.1f K = %.0f m of ascent per settling step",
                    th_stag - th_grid,
                    1000 * (th_stag - th_grid) / ((thcol[k0+1]-thcol[k0]) / (DZ/1000))))

@printf("\n%d passed, %d failed\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
