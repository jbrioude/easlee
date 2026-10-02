#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# test_levdot_spacing.jl — the metre-per-second to level-per-second conversion.
#
# Every reader that is handed a vertical velocity in m/s has to divide it by the
# LOCAL level spacing to get a level-index rate.  For an interior level that
# spacing is the centred half span, (z[k+1] - z[k-1]) / 2.  Dividing by the full
# two-level span instead makes the rate four times too small, which is what the
# WRF reader (initialize.jl, the `W` branch) and the AROME reader (the `NH_dW`
# branch) both did until 2026-10-02: `w / (z[k+1] - z[k-1]) / 2` parses as
# `w / (2·span)`, because `/` is left-associative.
#
# The checks below are the ones that would have caught it:
#   [1] uniform column: the conversion returns w/dz exactly,
#   [2] uniform column: an interior level agrees with the one-sided boundary
#       formula the same readers use at k = 1 and k = Nk,
#   [3] stretched column: the rate follows the local spacing, not the span,
#   [4] regression guard: the result is NOT the old four-times-too-small value,
#   [5] a centred DERIVATIVE on the same stencil still divides by the FULL span
#       (drho/dz and friends must not pick up the factor).
#
#   julia --project=easlee test/test_levdot_spacing.jl
# ─────────────────────────────────────────────────────────────────────────────
using Printf

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    if ok
        npass += 1; println("  PASS  ", name)
    else
        nfail += 1; println("  FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

≈ₜ(a, b) = isapprox(a, b; rtol=1e-6)

println("levdot spacing — m/s to level/s")

# ── [1] uniform column ──────────────────────────────────────────────────────
dz = 10.0
z  = collect(0.0:dz:100.0)          # z[k] = (k-1)·dz
w  = 2.5
k  = 5
rate = levrate_from_w(w, z[k-1], z[k+1])
check("uniform column: w/dz exactly", rate ≈ₜ w/dz;
      info=@sprintf("got %.6f, expected %.6f", rate, w/dz))

# ── [2] interior agrees with the one-sided boundary formula ─────────────────
# The readers use w/(z[2]-z[1]) at the bottom and w/(z[Nk]-z[Nk-1]) at the top.
# On a uniform grid an interior level must give the same number: a jump there
# would mean the discretisation is inconsistent with its own boundaries.
bottom = w / (z[2] - z[1])
top    = w / (z[end] - z[end-1])
check("interior matches the one-sided bottom", rate ≈ₜ bottom;
      info=@sprintf("interior %.6f vs bottom %.6f", rate, bottom))
check("interior matches the one-sided top", rate ≈ₜ top;
      info=@sprintf("interior %.6f vs top %.6f", rate, top))

# ── [3] stretched column: local spacing, not the span ───────────────────────
zs = [0.0, 10.0, 25.0, 50.0, 90.0, 140.0]      # growing spacing
for kk in 2:length(zs)-1
    local_dz = (zs[kk+1] - zs[kk-1]) / 2
    got      = levrate_from_w(w, zs[kk-1], zs[kk+1])
    check(@sprintf("stretched column, level %d", kk), got ≈ₜ w/local_dz;
          info=@sprintf("got %.6f, expected %.6f", got, w/local_dz))
end

# ── [4] regression guard on the exact historical bug ────────────────────────
wrong = w / (z[k+1] - z[k-1]) / 2          # what the two readers used to compute
check("not the four-times-too-small form", !(rate ≈ₜ wrong);
      info=@sprintf("rate %.6f, old buggy value %.6f", rate, wrong))
check("the old form was wrong by exactly 4", rate / wrong ≈ₜ 4.0;
      info=@sprintf("ratio %.6f", rate / wrong))

# ── [5] a centred derivative keeps the FULL span ────────────────────────────
# Same stencil, different quantity: d(rho)/dz is (rho[k+1]-rho[k-1])/(z[k+1]-z[k-1]).
# With rho = a·z the derivative must come back as a, with no factor of two.
a    = 0.3
rho  = a .* z
drho = (rho[k+1] - rho[k-1]) / (z[k+1] - z[k-1])
check("centred derivative divides by the full span", drho ≈ₜ a;
      info=@sprintf("got %.6f, expected %.6f", drho, a))

# ── [6] broadcast form, as the readers call it ─────────────────────────────
W3   = fill(w, 2, 2, length(z)-2)
zb   = reshape(repeat(z[1:end-2], inner=4), 2, 2, :)
za   = reshape(repeat(z[3:end],   inner=4), 2, 2, :)
rate3 = levrate_from_w.(W3, zb, za)
check("broadcast over a 3-D slab", all(r -> r ≈ₜ w/dz, rate3);
      info=@sprintf("min %.6f max %.6f", minimum(rate3), maximum(rate3)))

@printf("\n%d passed, %d failed\n", npass, nfail)
nfail == 0 || println("FAIL")
exit(nfail == 0 ? 0 : 1)
