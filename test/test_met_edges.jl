#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# test_met_edges.jl — met-reader edge cases fixed in the 2026-09-04 review (A3):
#   [1] surface-flux unit detector (_sfcflux_is_accumulated / _sfcflux_divisor /
#       _sfcflux_convention!): per-second vs accumulated values give the same W/m²,
#       header period honoured, 3600 s fallback, #info override, once-per-run cache;
#       the windowed-accumulation decision (_flux_window_mode!) and no second
#       division of window-converted fields (2026-10-22);
#   [2] WRF time-slice matcher (_wrf_time_slice): loops over the TIME dimension,
#       errors naming the file/date/span when no slice matches;
#   [3] first-call met bracket (_first_met_bracket): clamped at both archive
#       ends, explicit error outside the span;
#   [4] one-pass GRIB store (_store_grib_slice!): bit-identical to the old
#       map+reverse expression, flip / no-flip / crop / missing values;
#   [5] prefetch end-of-archive guard (readmetfiles_prefetch! holds the bracket).
# Synthetic only — no met archive needed.
# ─────────────────────────────────────────────────────────────────────────────
using Dates, Printf, Logging

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end
errs(f) = try; f(); false; catch; true; end
errmsg(f) = try; f(); ""; catch e; sprint(showerror, e); end

println("="^78); println("[1] surface-flux unit detector"); println("="^78)
Ni, Nj = 12, 7
sshf_w = Float32.(200.0 .* (rand(Ni, Nj) .- 0.6))            # W/m², max |.| ≤ 200
sshf_j = sshf_w .* 3600.0f0                                   # J/m² accumulated over 1 h
tau_n  = Float32.(0.8 .* rand(Ni, Nj))                        # N/m²
tau_a  = tau_n .* 3600.0f0                                    # N s/m²
check("W/m² field (max 200) → not accumulated", !_sfcflux_is_accumulated(maximum(abs, sshf_w), _SFCFLUX_SSHF_THRESHOLD))
check("J/m² field (max ~7e5) → accumulated",     _sfcflux_is_accumulated(maximum(abs, sshf_j), _SFCFLUX_SSHF_THRESHOLD))
check("N/m² stress (max 0.8) → not accumulated", !_sfcflux_is_accumulated(maximum(abs, tau_n), _SFCFLUX_TAU_THRESHOLD))
check("N s/m² stress (max ~3e3) → accumulated",  _sfcflux_is_accumulated(maximum(abs, tau_a), _SFCFLUX_TAU_THRESHOLD))
d_w = _sfcflux_divisor(maximum(abs, sshf_w), _SFCFLUX_SSHF_THRESHOLD, 3600)
d_j = _sfcflux_divisor(maximum(abs, sshf_j), _SFCFLUX_SSHF_THRESHOLD, 3600)
check("divisor = 1 for W/m², 3600 for J/m²", d_w == 1.0f0 && d_j == 3600.0f0; info="$d_w $d_j")
check("both representations give the same W/m² field",
      maximum(abs.(sshf_j ./ d_j .- sshf_w ./ d_w)) < 1e-4; info=@sprintf("max diff %.2e", maximum(abs.(sshf_j ./ d_j .- sshf_w ./ d_w))))

# full helper: decision + cache + header period
mk2D(sshf, tau) = Dict{String,Any}("sshf" => reshape(copy(sshf), Ni, Nj, 1),
                                   "ewss" => reshape(copy(tau), Ni, Nj, 1),
                                   "nsss" => reshape(copy(tau) .* 0.5f0, Ni, Nj, 1))
Logging.with_logger(Logging.NullLogger()) do
    mc = Dict{Any,Any}("sfcflux_accum" => nothing, "sfcflux_step_s" => 0)
    div = _sfcflux_convention!(mc, mk2D(sshf_w, tau_n), 1, true, true)
    check("per-second archive (zion EAg): divisors (1, 1)", div == (1.0f0, 1.0f0); info="$div")
    check("decision cached in metconfig", haskey(mc, "sfcflux_div"))
    div2 = _sfcflux_convention!(mc, mk2D(sshf_j, tau_a), 1, true, true)
    check("cache wins on later files (per run, not per file)", div2 == div)

    mc = Dict{Any,Any}("sfcflux_accum" => nothing, "sfcflux_step_s" => 3 * 3600)
    div = _sfcflux_convention!(mc, mk2D(sshf_j .* 3, tau_a .* 3), 1, true, true)
    check("accumulated archive: period from GRIB header (3 h)", div == (10800.0f0, 10800.0f0); info="$div")

    mc = Dict{Any,Any}("sfcflux_accum" => nothing, "sfcflux_step_s" => 0)
    div = _sfcflux_convention!(mc, mk2D(sshf_j, tau_a), 1, true, true)
    check("accumulated archive, header step 0 → 3600 s assumed", div == (3600.0f0, 3600.0f0); info="$div")

    mc = Dict{Any,Any}("sfcflux_accum" => nothing)          # no step key at all
    div = _sfcflux_convention!(mc, mk2D(sshf_j, tau_n), 1, true, true)
    check("mixed: sshf accumulated, stress per-second", div == (3600.0f0, 1.0f0); info="$div")

    mc = Dict{Any,Any}("sfcflux_accum" => true, "sfcflux_step_s" => 0)
    div = _sfcflux_convention!(mc, mk2D(sshf_w, tau_n), 1, true, true)
    check("#info sfcflux_units = accum overrides the magnitude test", div == (3600.0f0, 3600.0f0); info="$div")
    mc = Dict{Any,Any}("sfcflux_accum" => false, "sfcflux_step_s" => 3600)
    div = _sfcflux_convention!(mc, mk2D(sshf_j, tau_a), 1, true, true)
    check("#info sfcflux_units = per-s overrides the magnitude test", div == (1.0f0, 1.0f0); info="$div")

    # WRF HFX / AROME: instantaneous W/m² with a strong convective peak (~600 W/m²)
    # — must stay undivided (the old unconditional ÷3600 broke these runs too)
    hfx = Float32.(600.0 .* rand(Ni, Nj)); hfx[3, 3] = 600.0f0
    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    div = _sfcflux_convention!(mc, mk2D(hfx, tau_n), 1, true, true)
    check("WRF-magnitude sshf (max 600 W/m²) → no division, model-agnostic", div == (1.0f0, 1.0f0); info="$div")
    check("…same decision through the pure helper", _sfcflux_divisor(600.0, _SFCFLUX_SSHF_THRESHOLD, 3600) == 1.0f0)

    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    div = _sfcflux_convention!(mc, mk2D(zeros(Float32, Ni, Nj), tau_n), 1, true, true)
    check("all-zero flux slot: neutral divisors, NOT cached", div == (1.0f0, 1.0f0) && !haskey(mc, "sfcflux_div"))
    div = _sfcflux_convention!(mc, mk2D(sshf_j, tau_a), 1, true, true)
    check("…and the next file with data decides", div == (3600.0f0, 3600.0f0); info="$div")

    # 2026-10-22: fields already converted from their GRIB window (raw MARS/CDS
    # stepRange s0-e) must not be divided again, even under sfcflux_units = accum
    mc = Dict{Any,Any}("sfcflux_accum" => true, "sfcflux_step_s" => 3600,
                       "flux_window_converted" => ["sshf", "ewss", "nsss"])
    div = _sfcflux_convention!(mc, mk2D(sshf_w, tau_n), 1, true, true)
    check("accum declared, fields window-converted → divisors (1, 1)", div == (1.0f0, 1.0f0); info="$div")

    # one decision per met source for the windowed stepType=accum fields
    cn = Dict("sshf" => "sshf", "ewss" => "ewss", "nsss" => "nsss", "avg_ishf" => "shf", "prate" => "prate")
    function wscan(kinds)
        sc = GribScan(String[])
        for (n, k) in kinds; sc.acckind[n] = k; sc.acc2d[n] = (8, 9, 3600, 202201141800); end
        return sc
    end
    todo3 = [(n, (8, 9, 3600, 202201141800)) for n in ("sshf", "ewss", "nsss")]
    sc3 = wscan(["sshf" => "accum", "ewss" => "accum", "nsss" => "accum"])
    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(sshf_j, tau_a), 1, Ni)
    check("windowed accum, J/m² magnitudes → follow the header, cached", m == "header" && get(mc, "flux_window", "") == "header")
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(sshf_w, tau_n), 1, Ni)
    check("…the cached decision holds for later files", m == "header")
    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(sshf_w, tau_n), 1, Ni)
    check("windowed accum header over W/m² values → rates, cached", m == "rate" && get(mc, "flux_window", "") == "rate")
    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(zeros(Float32, Ni, Nj), zeros(Float32, Ni, Nj)), 1, Ni)
    check("all-zero windowed fields → header, NOT cached", m == "header" && !haskey(mc, "flux_window"))
    mc = Dict{Any,Any}("sfcflux_accum" => false)
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(sshf_j, tau_a), 1, Ni)
    check("sfcflux_units = per-s wins over J/m² magnitudes", m == "rate")
    mc = Dict{Any,Any}("sfcflux_accum" => true)
    m = _flux_window_mode!(mc, todo3, sc3, cn, mk2D(sshf_w, tau_n), 1, Ni)
    check("sfcflux_units = accum wins over W/m² magnitudes", m == "header")
    # GFS: only stepType=avg windows (W/m² bucket means) — no magnitude evidence
    sca = wscan(["avg_ishf" => "avg", "prate" => "avg"])
    f2 = Dict{String,Any}("shf" => reshape(copy(sshf_w), Ni, Nj, 1), "prate" => zeros(Float32, Ni, Nj, 1))
    mc = Dict{Any,Any}("sfcflux_accum" => nothing)
    m = _flux_window_mode!(mc, [("avg_ishf", (0, 6, 3600, 0)), ("prate", (0, 6, 3600, 0))], sca, cn, f2, 1, Ni)
    check("avg-only windows (GFS) → header, NOT cached", m == "header" && !haskey(mc, "flux_window"))
end

println("="^78); println("[2] WRF time-slice matcher"); println("="^78)
dates = ["2020-03-01_00:00:00", "2020-03-01_06:00:00", "2020-03-01_12:00:00"]
Times = Matrix{Char}(undef, 19, 3)
for (k, d) in enumerate(dates); Times[:, k] .= collect(d); end
check("third slice found (index 3)", _wrf_time_slice(Times, DateTime(2020, 3, 1, 12), "wrfout_d01") == 3)
check("first slice found (index 1)", _wrf_time_slice(Times, DateTime(2020, 3, 1, 0), "wrfout_d01") == 1)
m = errmsg(() -> _wrf_time_slice(Times, DateTime(2020, 3, 1, 18), "/x/wrfout_d01"))
check("missing date → error", !isempty(m))
check("…naming file, requested date and the file's first/last dates",
      occursin("wrfout_d01", m) && occursin("2020-03-01T18:00:00", m) &&
      occursin("2020-03-01T00:00:00", m) && occursin("2020-03-01T12:00:00", m); info=m)
# the old loop iterated the 19 characters of Times[:,1]: with 3 slices it would
# have walked past the last column — make sure a 1-slice file also works
check("single-slice file", _wrf_time_slice(Times[:, 1:1], DateTime(2020, 3, 1, 0), "f") == 1)

println("="^78); println("[3] first-call met bracket"); println("="^78)
fd = [DateTime(2013, 9, 17, h) for h in 0:5]           # 6 hourly files
nf = length(fd)
check("forward, start exactly at first file → (1,2)",  _first_met_bracket(fd, fd[1], 1) == 1)
check("forward, start inside → i-1",                    _first_met_bracket(fd, DateTime(2013, 9, 17, 2, 30), 1) == 3)
check("forward, start ON an analysis → that file first", _first_met_bracket(fd, fd[3], 1) == 2)
check("forward, start exactly at LAST file → (nf-1, nf)", _first_met_bracket(fd, fd[end], 1) == nf - 1)
check("backward, start exactly at LAST file → (nf, nf-1)", _first_met_bracket(fd, fd[end], -1) == nf)
check("backward, start inside → i+1",                   _first_met_bracket(fd, DateTime(2013, 9, 17, 2, 30), -1) == 4)
check("backward, start exactly at first file → (2,1)", _first_met_bracket(fd, fd[1], -1) == 2)
m = errmsg(() -> _first_met_bracket(fd, DateTime(2013, 9, 18, 3), 1))
check("forward start AFTER the archive → error (was MethodError nothing+1)", !isempty(m))
check("…naming start date and archive span",
      occursin("2013-09-18T03:00:00", m) && occursin("2013-09-17T00:00:00", m) && occursin("2013-09-17T05:00:00", m); info=m)
m = errmsg(() -> _first_met_bracket(fd, DateTime(2013, 9, 16, 23), -1))
check("backward start BEFORE the archive → error", occursin("outside the met archive span", m))
m = errmsg(() -> _first_met_bracket(fd, DateTime(2013, 9, 16, 23), 1))
check("forward start BEFORE the archive → error (was silent clamp to 1)", occursin("outside", m))
check("single-file archive → error", errs(() -> _first_met_bracket(fd[1:1], fd[1], 1)))

println("="^78); println("[4] one-pass GRIB store vs old map+reverse expression"); println("="^78)
Nm, Mm = 9, 5
vals = 1000.0 .* rand(Nm, Mm) .- 500.0
mv = 9999.0
vals[2, 3] = mv; vals[7, 1] = mv
old_store(v, has_mv, flip) = begin
    w = has_mv ? map(x -> x == mv ? NaN : x, v) : v
    Float32.(flip ? reverse(w, dims=2) : w)
end
same(a, b) = size(a) == size(b) && all(reinterpret(UInt32, vec(a)) .== reinterpret(UInt32, vec(b)))
for (has_mv, flip) in ((true, true), (false, true), (true, false), (false, false))
    dst = zeros(Float32, Nm, Mm)
    _store_grib_slice!(dst, vals, has_mv, mv, flip, nothing)
    check("bit-identical store (missing=$has_mv, flip=$flip)", same(dst, old_store(vals, has_mv, flip)))
end
iidx = Int32[8, 9, 1, 2]; jr = 2:4                     # seam-crossing lon window, lat sub-range
dst = zeros(Float32, length(iidx), length(jr))
_store_grib_slice!(dst, vals, true, mv, true, (iidx, jr))
check("cropped store == old `@view reverse(v,dims=2)[iidx, jr]`",
      same(dst, Float32.(reverse(map(x -> x == mv ? NaN : x, vals), dims=2)[iidx, jr])))
check("wrong row count → error", errs(() -> _store_grib_slice!(zeros(Float32, Nm, Mm + 1), vals, false, mv, true, nothing)))

println("="^78); println("[5] prefetch end-of-archive guard"); println("="^78)
# minimal duck-typed nest context: the guard runs before any file access
mutable struct _NRC; t_index::NTuple{2,Int}; end
struct _MC; filedate::Vector{DateTime}; filename::Vector{String}; metfileformat::String; model::String; end
struct _GC; globalx::Bool; end
struct _Nest; metconfig::_MC; nestrunconfig::_NRC; gridconfig::_GC; end
files = ["f$h" for h in 0:nf-1]
pf = MetPrefetcher(Dict{Any,Any}("T" => zeros(Float32, 2, 2, 2, 2)), Dict{Any,Any}("sp" => zeros(Float32, 2, 2, 2)))
pf.enabled = false
nest = _Nest(_MC(fd, files, "GRIB", "ERA5"), _NRC((nf - 1, nf)), _GC(false))
m = errmsg(() -> readmetfiles_prefetch!(pf, nest, fd[end] + Hour(3), 1, Dict{Any,Any}(), Dict{Any,Any}(), Dict{Any,Any}(), nothing))
check("forward past the last file: no BoundsError, bracket held", isempty(m) && nest.nestrunconfig.t_index == (nf - 1, nf); info=m)
nest = _Nest(_MC(fd, files, "GRIB", "ERA5"), _NRC((2, 1)), _GC(false))
m = errmsg(() -> readmetfiles_prefetch!(pf, nest, fd[1] - Hour(3), -1, Dict{Any,Any}(), Dict{Any,Any}(), Dict{Any,Any}(), nothing))
check("backward before the first file: no BoundsError, bracket held", isempty(m) && nest.nestrunconfig.t_index == (2, 1); info=m)

println("\n", npass, " passed, ", nfail, " failed")
nfail == 0 || (println("FAIL"); exit(1))
