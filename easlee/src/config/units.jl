
# ── Human-readable unit parsing for the run config ────────────────────────────
# Config values may be written with a unit suffix ("15min", "2.0um", "1kg") so
# the file is self-documenting.  Each parser returns the SI base value as the
# concrete type EASLEE uses internally, and accepts a bare number (already SI)
# unchanged.

module Units

export parse_duration_s, parse_length_m, parse_mass_kg, format_duration

const _DUR = Dict("ms"=>1e-3, "s"=>1.0, "sec"=>1.0, "min"=>60.0, "m"=>60.0,
                  "h"=>3600.0, "hr"=>3600.0, "hour"=>3600.0, "d"=>86400.0, "day"=>86400.0)
const _LEN = Dict("nm"=>1e-9, "um"=>1e-6, "µm"=>1e-6, "mm"=>1e-3, "cm"=>1e-2,
                  "m"=>1.0, "km"=>1e3)
const _MASS = Dict("ng"=>1e-12, "ug"=>1e-9, "µg"=>1e-9, "mg"=>1e-6, "g"=>1e-3,
                   "kg"=>1.0, "t"=>1e3, "ton"=>1e3,
                   # large source terms (volcanic/fire): kilo/mega/giga/teragram
                   "kt"=>1e6, "gg"=>1e6, "mt"=>1e9, "tg"=>1e9)

# Split a string like "2.0um" / "15 min" into (number, unit).  Returns
# (value, "") when there is no unit suffix.
function _split_unit(s::AbstractString)
    m = match(r"^\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([a-zA-Zµ]*)\s*$", s)
    m === nothing && error("cannot parse value with unit: \"$s\"")
    return parse(Float64, m.captures[1]), lowercase(String(m.captures[2]))
end

_parse(x::Real, ::Dict) = Float64(x)
function _parse(x::AbstractString, table::Dict)
    val, u = _split_unit(x)
    isempty(u) && return val
    haskey(table, u) || error("unknown unit \"$u\" in \"$x\"")
    return val * table[u]
end

"""Duration → seconds.  Accepts "15min", "1h", "900", 900."""
parse_duration_s(x) = _parse(x, _DUR)

"""Length → metres.  Accepts "2.0um", "500m", 0.3e-6."""
parse_length_m(x) = _parse(x, _LEN)

"""Mass → kilograms.  Accepts "1kg", "5g", 1.0."""
parse_mass_kg(x) = _parse(x, _MASS)

"""Seconds → compact human string for writing configs back out."""
function format_duration(sec::Real)
    sec == 0 && return "0s"
    for (u, f) in (("d", 86400.0), ("h", 3600.0), ("min", 60.0), ("s", 1.0))
        if sec >= f && isapprox(sec / f, round(sec / f); atol=1e-9)
            n = round(Int, sec / f)
            return "$(n)$(u)"
        end
    end
    return "$(sec)s"
end

end # module Units
