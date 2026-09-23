# test_surface_layer.jl — Monin-Obukhov surface-layer wind factor
# (physics/surface_layer.jl).  The archived wind is a SLAB below the lowest
# model level; this puts a parcel back on the similarity profile anchored at
# that level.  Pinned properties: exact log law in neutral conditions, correct
# stability ordering, no correction above the anchor, and hard safety limits.

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "physics", "surface_layer.jl"))

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name)
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

κ  = 0.4f0
z1 = 19.6f0        # ETEX archive median lowest-level height
U1 = 8.0f0
us = 0.35f0        # u* consistent with that wind over land

# ── neutral: u(z) = U1 − (u*/κ)·ln(z1/z) exactly ─────────────────────────────
for z in (2.0f0, 5.0f0, 10.0f0, 15.0f0)
    want = (U1 - (us / κ) * log(z1 / z)) / U1
    got  = mo_wind_factor(z, z1, us, 0.0f0, U1)
    check("neutral log law at z=$(z) m", abs(got - want) < 1e-6,
          info = "got $got want $want")
end

# equivalent statement: the implied roughness is a single consistent z0
z0 = z1 * exp(-κ * U1 / us)
check("implied z0 is physical (1e-4 … 1 m)", 1e-4 < z0 < 1.0, info = "z0=$z0")
for z in (2.0f0, 10.0f0)
    want = Float32(log(z / z0) / log(z1 / z0))
    got  = mo_wind_factor(z, z1, us, 0.0f0, U1)
    check("neutral ≡ ln(z/z0)/ln(z1/z0) at z=$(z) m", abs(got - want) < 1e-5,
          info = "got $got want $want")
end

# ── no correction at or above the anchor ─────────────────────────────────────
check("factor == 1 exactly at z1", mo_wind_factor(z1, z1, us, 0.0f0, U1) == 1.0f0)
check("factor == 1 above z1", mo_wind_factor(50.0f0, z1, us, 0.0f0, U1) == 1.0f0)

# ── stability ordering: unstable (well mixed) > neutral > stable ─────────────
z = 5.0f0
f_unst = mo_wind_factor(z, z1, us, -50.0f0, U1)    # L < 0
f_neut = mo_wind_factor(z, z1, us,   0.0f0, U1)
f_stab = mo_wind_factor(z, z1, us,  50.0f0, U1)    # L > 0
check("unstable > neutral > stable", f_unst > f_neut > f_stab,
      info = "unstable=$f_unst neutral=$f_neut stable=$f_stab")
check("unstable flattens the profile vs neutral", f_unst > f_neut + 0.03f0,
      info = "unstable=$f_unst neutral=$f_neut")
check("stable profile has strong shear", f_stab < f_neut - 0.03f0,
      info = "stable=$f_stab neutral=$f_neut")

# ── magnitude reproduces the atlas measurement for the LAND case ─────────────
# The atlas quotes 0.568 at 2 m over land, which corresponds to z0 = 0.1 m —
# i.e. u* = κ·U1/ln(z1/z0).  The factor is set by u*/U1, not by u* alone, so
# the check must use a land-consistent pair.
us_land = Float32(0.4 * U1 / log(z1 / 0.1f0))
f2 = mo_wind_factor(2.0f0, z1, us_land, 0.0f0, U1)
check("2 m land factor reproduces the atlas 0.568", abs(f2 - 0.568f0) < 0.01f0,
      info = "$f2 with u*=$(round(us_land, digits=3))")
f10 = mo_wind_factor(10.0f0, z1, us_land, 0.0f0, U1)
check("10 m land factor reproduces the atlas 0.873", abs(f10 - 0.873f0) < 0.01f0,
      info = "$f10")

# ── safety: never amplifies, never negative, degenerate inputs are no-ops ────
check("never > 1 (very large u*)", mo_wind_factor(2.0f0, z1, 3.0f0, 0.0f0, U1) <= 1.0f0)
check("never < 0 (very large u*)", mo_wind_factor(2.0f0, z1, 3.0f0, 0.0f0, U1) >= 0.0f0)
check("calm wind → no-op", mo_wind_factor(2.0f0, z1, us, 0.0f0, 0.0f0) == 1.0f0)
check("NaN height → no-op", mo_wind_factor(NaN32, z1, us, 0.0f0, U1) == 1.0f0)
check("NaN ustar → no-op", mo_wind_factor(2.0f0, z1, NaN32, 0.0f0, U1) == 1.0f0)
check("zero ustar → no-op", mo_wind_factor(2.0f0, z1, 0.0f0, 0.0f0, U1) == 1.0f0)
check("shallow anchor → no-op", mo_wind_factor(0.2f0, 0.3f0, us, 0.0f0, U1) == 1.0f0)
check("below the 0.5 m floor is clamped, not extrapolated",
      mo_wind_factor(0.01f0, z1, us, 0.0f0, U1) ==
      mo_wind_factor(0.5f0,  z1, us, 0.0f0, U1))

# ── monotone in height, finite everywhere on a fine sweep ───────────────────
let prev = -1.0f0, mono = true, finite = true
    for z in 0.5f0:0.25f0:19.5f0
        f = mo_wind_factor(z, z1, us, 20.0f0, U1)
        isfinite(f) || (finite = false)
        f < prev - 1f-6 && (mono = false)
        prev = f
    end
    check("monotone increasing with height (stable)", mono)
    check("finite over the whole layer", finite)
end

println()
println(npass, " passed, ", nfail, " failed")
nfail == 0 || println("FAIL")
