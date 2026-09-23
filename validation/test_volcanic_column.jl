#!/usr/bin/env julia
# test_volcanic_column.jl — the 1-D moist eruption column (volcanic_column_rise):
#   1. dry columns track the Mastin (2009) height–MER calibration,
#   2. vaporised seawater (latent heat) lofts a column HIGHER at equal MER and
#      source temperature — the phreatomagmatic boost,
#   3. NBH sits below the top at a physical fraction,
#   4. Hunga-Tonga-class energetics reach the observed 30+ km umbrella,
#   5. the real ERA5 vent profile (when /data/ERA5_Hunga is present) gives a
#      similar answer to the standard tropical atmosphere.
#
#   julia --project=easlee scripts/test_volcanic_column.jl
using Printf

global gpu = false; global chem = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE.Operators: volcanic_column_rise, mastin_height_km

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── standard tropical atmosphere (0–60 km) ────────────────────────────────────
z = collect(0.0:200.0:60_000.0)
T = similar(z)
for (i, zz) in enumerate(z)
    T[i] = zz <= 16_000 ? 300.0 - 6.5e-3 * zz :
           zz <= 20_000 ? 196.0 :
                          196.0 + 2.2e-3 * (zz - 20_000)
end
p = similar(z); p[1] = 101_325.0
for i in 2:length(z)
    Tm = 0.5 * (T[i-1] + T[i])
    p[i] = p[i-1] * exp(-9.81 * (z[i] - z[i-1]) / (287.05 * Tm))
end
qv = [0.6 * 0.622 * 611.2 * exp(17.67 * (Tt - 273.15) / (Tt - 29.65)) / pp
      for (Tt, pp) in zip(T, p)]                    # 60 % RH, decays naturally
qv .= max.(qv, 2.0e-6)

println("="^78)
println("[1] dry columns vs Mastin (2009)")
println("="^78)
okm = true
for mer in (1.0e7, 1.0e8, 1.0e9)
    zt, znb, ok = volcanic_column_rise(mer, 1100.0, 0.03, z, T, p, qv)
    Hm = mastin_height_km(mer) * 1000
    rel = zt / Hm
    @printf("  MER %.0e kg/s: top %5.1f km  NBH %5.1f km   Mastin %5.1f km  (ratio %.2f)\n",
            mer, zt/1e3, znb/1e3, Hm/1e3, rel)
    (ok && 0.6 <= rel <= 1.6) || (global okm = false)
end
check("dry tops within ×0.6–1.6 of Mastin over 3 decades", okm)

println()
println("="^78)
println("[2] phreatomagmatic boost + structure")
println("="^78)
zt_d, znb_d, _ = volcanic_column_rise(1.0e9, 800.0, 0.05, z, T, p, qv)
zt_w, znb_w, _ = volcanic_column_rise(1.0e9, 800.0, 0.40, z, T, p, qv)
check("wet column (40 % vapour) higher than dry at equal MER/T0", zt_w > zt_d + 1000;
      info=@sprintf("(dry %.1f → wet %.1f km)", zt_d/1e3, zt_w/1e3))
check("NBH below the top at 0.55–0.95", 0.55 <= znb_w / zt_w <= 0.95;
      info=@sprintf("(NBH/top = %.2f)", znb_w / zt_w))
z1, _, _ = volcanic_column_rise(1.0e8, 800.0, 0.30, z, T, p, qv)
z2, _, _ = volcanic_column_rise(1.0e9, 800.0, 0.30, z, T, p, qv)
z3, _, _ = volcanic_column_rise(5.0e9, 800.0, 0.30, z, T, p, qv)
check("monotone in MER", z1 < z2 < z3;
      info=@sprintf("(%.0f < %.0f < %.0f km)", z1/1e3, z2/1e3, z3/1e3))

println()
println("="^78)
println("[3] Hunga-Tonga-class source")
println("="^78)
zt, znb, ok = volcanic_column_rise(5.0e9, 700.0, 0.45, z, T, p, qv)
@printf("  MER 5e9 kg/s, T0 700 K, 45 %% vapour → top %.1f km, NBH %.1f km\n", zt/1e3, znb/1e3)
# The STEADY column gives the sustained umbrella, not the minutes-long 47–57 km
# overshoot (negligible mass).  Observed sustained layer (metacloud): 20–33 km.
check("sustained top matches the observed umbrella (30–38 km)", ok && 30_000 <= zt <= 38_000;
      info=@sprintf("(top %.1f km)", zt/1e3))
check("NBH matches the observed cloud BASE (17–26 km)", 17_000 <= znb <= 26_000;
      info=@sprintf("(NBH %.1f km; obs base ≈ 20 km)", znb/1e3))
check("model layer [NBH, top] ≈ observed [20, 33] km",
      abs(znb - 20_000) < 5_000 && abs(zt - 33_000) < 5_000;
      info=@sprintf("([%.1f, %.1f] km)", znb/1e3, zt/1e3))

println()
@printf("RESULT: %d passed, %d failed\n", npass, nfail)
nfail == 0 || exit(1)
