# datem.jl — DATEM/ATMES-II observation parsers + statistics, shared by
# validate_captex.jl (CAPTEX/ETEX validation task) and check_etex.jl (the
# self-check of a plain transport run's time-framed grid output).
# Requires: using Dates, Statistics in the includer.

# ── DATEM parsers ─────────────────────────────────────────────────────────────
_hhmm(s) = (h = parse(Int, s[1:2]); m = parse(Int, s[3:4]); (h, m))

"Releases: (start, stop, lat, lon, mass_g)."
function read_emit(path)
    rels = NamedTuple[]
    for ln in eachline(path)
        occursin(r"^\s*\d{4}\s", ln) || continue
        f = split(ln)
        yr, mo, dy = parse.(Int, f[1:3])
        (sh, sm) = _hhmm(f[4]); (dh, dm) = _hhmm(f[5])
        start = DateTime(yr, mo, dy, sh, sm)
        push!(rels, (start = start, stop = start + Hour(dh) + Minute(dm),
                     lat = parse(Float64, f[6]), lon = parse(Float64, f[7]),
                     mass_g = parse(Float64, f[8])))
    end
    rels
end

"Stations: Dict(id => (lat, lon))."
function read_stations(path)
    st = Dict{Int,Tuple{Float64,Float64}}()
    for ln in eachline(path)
        occursin(r"^\s*\d", ln) || continue
        f = split(ln); length(f) >= 3 || continue
        st[parse(Int, f[1])] = (parse(Float64, f[2]), parse(Float64, f[3]))
    end
    st
end

"Measurements: (start, stop, lat, lon, conc_pg, stn)."
function read_meas(path)
    obs = NamedTuple[]
    for ln in eachline(path)
        occursin(r"^\s*\d{4}\s", ln) || continue
        f = split(ln); length(f) >= 9 || continue
        yr, mo, dy = parse.(Int, f[1:3])
        (sh, sm) = _hhmm(f[4]); (dh, dm) = _hhmm(f[5])
        start = DateTime(yr, mo, dy, sh, sm)
        push!(obs, (start = start, stop = start + Hour(dh) + Minute(dm),
                    lat = parse(Float64, f[6]), lon = parse(Float64, f[7]),
                    conc_pg = parse(Float64, f[8]), stn = parse(Int, f[9])))
    end
    obs
end

# ── ATMES-II / DATEM statistics (Mosca et al. 1998; Draxler 2001) ─────────────
# Pairs (o, m) include zeros.  thr = threshold for the spatial overlap (FMS).
function datem_stats(o::Vector{<:Real}, m::Vector{<:Real}; thr = 0.0)
    n = length(o); @assert n == length(m) && n > 0
    ō = mean(o); m̄ = mean(m)
    fb   = (m̄ + ō) > 0 ? 2(m̄ - ō) / (m̄ + ō) : 0.0
    nmse = (m̄ * ō) > 0 ? mean((m .- o).^2) / (m̄ * ō) : NaN
    so = std(o); sm = std(m)
    r  = (so > 0 && sm > 0) ? cor(o, m) : 0.0
    # FAC2: fraction of pairs within a factor of 2 (both positive)
    pos = [(o[i] > 0 && m[i] > 0) for i in 1:n]
    fac2 = any(pos) ? count(i -> pos[i] && 0.5 <= m[i]/o[i] <= 2.0, 1:n) / count(pos) : 0.0
    # FMS (station overlap): both-above / either-above, in %
    both = count(i -> o[i] > thr && m[i] > thr, 1:n)
    either = count(i -> o[i] > thr || m[i] > thr, 1:n)
    fms = either > 0 ? 100 * both / either : 0.0
    # KSP: Kolmogorov–Smirnov parameter (max CDF difference, %)
    ksp = _ksp(o, m)
    # DATEM rank (0–4): R² + (1−|FB|/2) + FMS/100 + (1−KSP/100)
    rank = r^2 + (1 - abs(fb)/2) + fms/100 + (1 - ksp/100)
    return (; n, obar = ō, mbar = m̄, fb, nmse, r, fac2, fms, ksp, rank)
end

function _ksp(o, m)
    xs = sort(unique(vcat(o, m)))
    no = length(o); nm = length(m)
    d = 0.0
    for x in xs
        co = count(<=(x), o) / no
        cm = count(<=(x), m) / nm
        d = max(d, abs(co - cm))
    end
    return 100 * d
end

