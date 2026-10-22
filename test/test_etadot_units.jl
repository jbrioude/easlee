# test_etadot_units.jl — the archived vertical velocity comes in TWO conventions
# for the same ECMWF parameter 77:
#   raw MARS / ERA5 dissemination : eta-dot            [s^-1]     ~1e-6
#   flex_extract (calc_etadot)    : eta-dot * dp/deta  [Pa s^-1]   ~1e-1
# and on different levels (full for raw MARS, half for flex_extract: pinned at the end).
# They differ by dp/deta ~ 1e5, so the conversion to a level-index rate must
# pick the right divisor.  This pins the detector and both branches.

const REPO = abspath(joinpath(@__DIR__, ".."))

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name)
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

# Pull the two functions out of initialize.jl without its GRIB/NetCDF stack.
src = read(joinpath(REPO, "easlee", "src", "io", "initialize.jl"), String)
function extract(fname)
    i = findfirst("function $fname(", src)
    i === nothing && error("$fname not found in initialize.jl")
    j = first(i); depth = 0; k = j
    while true
        m = findnext(r"\bfunction\b|\bend\b|\bfor\b|\bif\b|\blet\b|\bdo\b|\bwhile\b", src, k)
        m === nothing && error("unbalanced")
        tok = src[m]
        depth += (tok == "end") ? -1 : 1
        k = last(m) + 1
        depth == 0 && return src[j:last(m)]
    end
end
using Logging
Logging.disable_logging(Logging.Info)
eval(Meta.parse(extract("_etadot_is_per_second")))
eval(Meta.parse(extract("_etadot_use_deta!")))
eval(Meta.parse(extract("_levdrho_chunk!")))

Ni, Nj, Nk = 6, 5, 30
mkfield(v) = (a = zeros(Float32, Ni, Nj, Nk, 2); a[:, :, :, 1] .= v; a)

# ── detector ─────────────────────────────────────────────────────────────────
check("detects raw eta-dot [s^-1]",  _etadot_is_per_second(mkfield(2.0f-6), Nk, 1))
check("detects flex_extract [Pa/s]", !_etadot_is_per_second(mkfield(1.0f-1), Nk, 1))
# two orders of margin either side of the 1e-2 threshold
check("margin: 1e-4 s^-1 still per-second", _etadot_is_per_second(mkfield(1.0f-4), Nk, 1))
check("margin: 1.0 Pa/s still Pa/s",       !_etadot_is_per_second(mkfield(1.0f0),  Nk, 1))
# a quiescent Pa/s domain must not be misread: the detector uses the MAX, and a
# genuinely tiny Pa/s field is indistinguishable — that is why the namelist
# override exists.  Pin the documented behaviour rather than pretend otherwise.
check("all-zero field falls to per-second (documented)",
      _etadot_is_per_second(mkfield(0.0f0), Nk, 1))
check("NaNs are skipped, not propagated",
      (f = mkfield(2.0f-6); f[1, 1, 5, 1] = NaN32; _etadot_is_per_second(f, Nk, 1)))

# ── the decision is made ONCE per met source, then cached ───────────────────
let mc = Dict{String,Any}("etadot_per_second" => nothing)
    first  = _etadot_use_deta!(mc, mkfield(2.0f-6), Nk, 1)     # s^-1 file
    cachedv = mc["etadot_per_second"]
    # a later file whose magnitude says Pa/s must NOT flip the cached answer
    second = _etadot_use_deta!(mc, mkfield(1.0f-1), Nk, 1)
    check("first file decides", first === true && cachedv === true)
    check("later files reuse the cached decision (no rescan)", second === true)
end
let mc = Dict{String,Any}("etadot_per_second" => false)         # explicit pa/s
    check("explicit namelist setting is never overridden",
          _etadot_use_deta!(mc, mkfield(2.0f-6), Nk, 1) === false)
end

# ── both conversion branches ─────────────────────────────────────────────────
# puv decreasing with k (k=1 surface), 1000 Pa per level; deta 0.01 per level.
puv = zeros(Float32, Ni, Nj, Nk, 2)
for k in 1:Nk; puv[:, :, k, 1] .= 100000.0f0 - 1000.0f0 * (k - 1); end
rho = ones(Float32, Ni, Nj, Nk, 2); drhodz = zeros(Float32, Ni, Nj, Nk, 2)
zlev = zeros(Float32, Ni, Nj, Nk, 2)
for k in 1:Nk; zlev[:, :, k, 1] .= 100.0f0 * (k - 1); end
deta = fill(-0.01, Nk)

f1 = mkfield(2.0f-6)
_levdrho_chunk!(2:Nk-1, f1, puv, rho, drhodz, zlev, deta, 1)
check("s^-1 branch divides by deta", isapprox(f1[1, 1, 10, 1], 2.0f-6 / -0.01f0; rtol = 1e-5),
      info = "got $(f1[1,1,10,1]) want $(2.0f-6 / -0.01f0)")

f2 = mkfield(1.0f-1)
_levdrho_chunk!(2:Nk-1, f2, puv, rho, drhodz, zlev, nothing, 1)
check("Pa/s branch divides by dp", isapprox(f2[1, 1, 10, 1], 1.0f-1 / -1000.0f0; rtol = 1e-5),
      info = "got $(f2[1,1,10,1]) want $(1.0f-1 / -1000.0f0)")

# the two conventions describing the SAME motion must give the SAME level rate:
# eta-dot = 2e-6 s^-1 with dp/deta = 1000/0.01 = 1e5 Pa  <=>  0.2 Pa/s
f3 = mkfield(2.0f-6);  _levdrho_chunk!(2:Nk-1, f3, puv, rho, drhodz, zlev, deta, 1)
f4 = mkfield(2.0f-1);  _levdrho_chunk!(2:Nk-1, f4, puv, rho, drhodz, zlev, nothing, 1)
check("both conventions agree on the same physical motion",
      isapprox(f3[1, 1, 10, 1], f4[1, 1, 10, 1]; rtol = 1e-4),
      info = "s^-1 → $(f3[1,1,10,1]) vs Pa/s → $(f4[1,1,10,1])")

# drho/dz is computed identically in both branches
check("drhodz unaffected by the branch", drhodz[1, 1, 10, 1] == 0.0f0)

# ── Post-conversion sanity check on the resulting vertical velocity ─────────
# _etadot_is_per_second judges the INPUT magnitude; _levdot_sanity! judges the
# OUTPUT, which is what catches an archive in neither convention, or a
# convention forced wrongly through `#info etadot_units`.
for c in ("_LEVDOT_W_LO", "_LEVDOT_W_HI")
    m = match(Regex("const $c = ([0-9.eEf+-]+)"), src)
    m === nothing && error("$c not found in initialize.jl")
    eval(Meta.parse("const $c = $(m.captures[1])"))
end
eval(Meta.parse(extract("_levdot_sanity!")))

# Collect log records so the test can assert on warn-vs-info.
struct Collector <: Logging.AbstractLogger; recs::Vector{Any}; end
Logging.min_enabled_level(::Collector) = Logging.Debug
Logging.shouldlog(::Collector, args...) = true
Logging.catch_exceptions(::Collector) = false
Logging.handle_message(c::Collector, lvl, msg, _mod, _grp, _id, _file, _line; kw...) =
    push!(c.recs, (lvl, string(msg)))
function run_check(levdot_val, dz_m; per_second = true)
    Ni = Nj = 4; Nk = 30
    ld = fill(Float32(levdot_val), Ni, Nj, Nk, 1)
    zl = Array{Float32}(undef, Ni, Nj, Nk, 1)
    for k in 1:Nk; zl[:, :, k, 1] .= Float32((k - 1) * dz_m); end
    mc = Dict{String,Any}("etadot_per_second" => per_second)
    col = Collector(Any[])
    # the file disables Info globally; re-enable it around the collector so the
    # "passed" branch is observable, then put the filter back.
    Logging.disable_logging(Logging.Debug)
    try
        Logging.with_logger(col) do
            _levdot_sanity!(mc, ld, zl, Nk, 1)
        end
    finally
        Logging.disable_logging(Logging.Info)
    end
    return mc, col.recs
end

# A physical field: 1e-4 level/s over 300 m layers -> |w| = 0.03 m/s.
mc1, r1 = run_check(1.0f-4, 300.0)
check("plausible w passes without a warning",
      !any(x -> x[1] == Logging.Warn, r1) && any(x -> x[1] == Logging.Info, r1),
      info = string(r1))
check("the check records that it ran", get(mc1, "levdot_checked", false) === true)

# The same motion with the wrong divisor: 1e5 times too small.
mc2, r2 = run_check(1.0f-9, 300.0)
check("w that is 1e5 too small raises a warning",
      any(x -> x[1] == Logging.Warn, r2), info = string(r2))
check("the warning names the namelist override",
      any(x -> x[1] == Logging.Warn && occursin("etadot_units", x[2]), r2))

# And 1e5 too large (the mirror error).
_, r3 = run_check(1.0f0, 300.0)
check("w that is far too large raises a warning",
      any(x -> x[1] == Logging.Warn, r3), info = string(r3))

# Runs once per met source: a second call on the same metconfig is a no-op.
begin
    Ni = Nj = 4; Nk = 30
    ld = fill(1.0f-9, Ni, Nj, Nk, 1)
    zl = Array{Float32}(undef, Ni, Nj, Nk, 1)
    for k in 1:Nk; zl[:, :, k, 1] .= Float32((k - 1) * 300.0); end
    mc = Dict{String,Any}("etadot_per_second" => true)
    Logging.disable_logging(Logging.Debug)
    c1 = Collector(Any[]); Logging.with_logger(c1) do; _levdot_sanity!(mc, ld, zl, Nk, 1); end
    c2 = Collector(Any[]); Logging.with_logger(c2) do; _levdot_sanity!(mc, ld, zl, Nk, 1); end
    Logging.disable_logging(Logging.Info)
    check("the check runs once per met source", !isempty(c1.recs) && isempty(c2.recs))
end

# The bounds bracket the documented ERA5 range (mean 1.9 cm/s, peaks 0.36 m/s).
check("bounds bracket the documented ERA5 peak",
      _LEVDOT_W_LO < 0.36 < _LEVDOT_W_HI)

# ── Levels the archived eta-dot sits on (2026-10-22) ────────────────────────
# flex_extract: half levels, level 137 = the ground (rms ratio 137/136 = 0.05);
# raw MARS: full levels (0.32) — validation/check_etadot_levels.py.
let m = match(r"const _ETADOT_HALF_RATIO = ([0-9.]+)", src)
    m === nothing && error("_ETADOT_HALF_RATIO not found in initialize.jl")
    eval(Meta.parse("const _ETADOT_HALF_RATIO = $(m.captures[1])"))
end
eval(Meta.parse(extract("_etadot_levels!")))
function etacol(r)          # slot 1 = r × slot 2, random sign pattern
    a = Float32.(randn(6, 5, 4, 1)); a[:, :, 1, 1] .= r .* a[:, :, 2, 1]; return a
end
mc = Dict{String,Any}("etadot_levels" => "auto")
check("rms ratio 0.05 (flex_extract) → half levels", _etadot_levels!(mc, etacol(0.05f0), 1) == "half")
check("…cached once per met source", get(mc, "_etadot_levels_auto", "") == "half" &&
      _etadot_levels!(mc, etacol(0.32f0), 1) == "half")
mc = Dict{String,Any}("etadot_levels" => "auto")
check("rms ratio 0.32 (raw MARS) → full levels", _etadot_levels!(mc, etacol(0.32f0), 1) == "full")
mc = Dict{String,Any}("etadot_levels" => "full")
check("#info etadot_levels = full wins over a half-level signature", _etadot_levels!(mc, etacol(0.05f0), 1) == "full")
mc = Dict{String,Any}("etadot_levels" => "half")
check("#info etadot_levels = half wins over a full-level signature", _etadot_levels!(mc, etacol(0.32f0), 1) == "half")
mc = Dict{String,Any}("etadot_levels" => "auto")
check("all-zero eta-dot: undecided, not cached",
      _etadot_levels!(mc, zeros(Float32, 6, 5, 4, 1), 1) == "half" && !haskey(mc, "_etadot_levels_auto"))

println()
println(npass, " passed, ", nfail, " failed")
nfail == 0 || println("FAIL")
