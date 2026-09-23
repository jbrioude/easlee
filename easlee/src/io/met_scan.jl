
# ── Automatic AVAILABLE generation ────────────────────────────────────────────
# Scans a meteorology directory and builds the (date → file) mapping the model
# needs, so a user never has to hand-write an AVAILABLE file.
#
# Strategy per file:
#   1. read the AUTHORITATIVE valid time from the header
#      (GRIB validityDate/validityTime — already accounts for forecast step;
#       NetCDF time/DTCUR/Times variable);
#   2. fall back to parsing the file name when the header read fails;
#   multi-time files (one file, several valid times — e.g. AROME-HDF, or any
#   NetCDF with a time dimension) expand to one entry per time.
#
# Header reads decode only keys (not field values), so scanning is cheap even
# for files with many messages.  Loaded into memory by build_run; the optional
# write_available! emits a legacy AVAILABLE file for inspection.

using Dates, GRIB, NetCDF

# (date, filepath) entries sorted by date, de-duplicated on (date, file).
function scan_met_directory(directory::AbstractString; model::AbstractString="",
                            format::AbstractString="", multistep::Bool=false,
                            recursive::Bool=true,
                            include_only::Vector{Regex}=Regex[],
                            exclude::Vector{Regex}=Regex[])
    files = String[]
    if recursive
        for (root, _, fs) in walkdir(directory), f in fs
            push!(files, joinpath(root, f))
        end
    else
        append!(files, readdir(directory; join=true))
    end
    files = filter(files) do f
        b = basename(f)
        isfile(f) || return false
        occursin(r"AVAILABLE|\.vars$|\.txt$|README|\.pdf$"i, b) && return false
        # When companions are declared, the MAIN files are those matching the
        # companion source token (_MET_, __HP1__, _asm_Nv.); companion files
        # (_CLD_, __SP1__, _mst_Ne.) and auxiliaries (PT) are read via the
        # main-file mechanism and must not appear as standalone entries.
        isempty(include_only) || any(p -> occursin(p, b), include_only) || return false
        any(p -> occursin(p, b), exclude) && return false
        return true
    end
    fmt = isempty(format) ? _guess_format(files, model) : lowercase(format)
    # (validtime, inittime, file) — inittime = forecast initialisation; for
    # reanalysis/analysis it equals validtime.
    raw = Tuple{DateTime,DateTime,String}[]
    for f in sort(files)
        try
            vts = fmt == "grib" ? _grib_validtimes(f; multistep=multistep) :
                                  [(d, d) for d in _netcdf_validtimes(f; model=model)]
            isempty(vts) && (d = _date_from_filename(basename(f)); d !== nothing && (vts = [(d, d)]))
            for (vt, it) in vts
                push!(raw, (vt, it, f))
            end
        catch err
            d = _date_from_filename(basename(f))
            d === nothing ? (@warn "skipping unreadable met file" file=basename(f) err) :
                            push!(raw, (d, d, f))
        end
    end
    # Dedupe on valid time, keeping the FRESHEST forecast (latest init time →
    # shortest lead).  This resolves overlapping cycles, e.g. 19:00 valid from
    # both the 06z (+13h) and the 18z (+1h) IFS run → keep the 18z file.
    best = Dict{DateTime,Tuple{DateTime,String}}()
    for (vt, it, f) in raw
        (!haskey(best, vt) || it > best[vt][1]) && (best[vt] = (it, f))
    end
    entries = [(vt, f) for (vt, (it, f)) in best]
    sort!(entries; by = first)
    return entries
end

# WRF/MESONH met is NetCDF but often extension-less (wrfout_d01_…), so extension
# sniffing would mis-guess GRIB — key off the model for those.
_guess_format(files, model="") =
    uppercase(model) in ("WRF", "MESONH") ? "netcdf" :
    any(f -> occursin(r"\.gri?b2?$"i, f), files) ? "grib" :
    any(f -> occursin(r"\.nc4?$"i, f), files)    ? "netcdf" : "grib"

# ── GRIB: (validity, initialisation) times from the message headers ───────────
function _grib_validtimes(path::AbstractString; multistep::Bool=false)
    out = Tuple{DateTime,DateTime}[]
    GribFile(path) do f
        for msg in f
            m = Message(msg)
            vd = Int(m["validityDate"]); vt = Int(m["validityTime"])
            valid = DateTime(string(vd), dateformat"yyyymmdd") +
                    Hour(div(vt, 100)) + Minute(mod(vt, 100))
            dd = Int(m["dataDate"]); dt = Int(m["dataTime"])
            init = DateTime(string(dd), dateformat"yyyymmdd") +
                   Hour(div(dt, 100)) + Minute(mod(dt, 100))
            any(x -> x[1] == valid, out) || push!(out, (valid, init))
            multistep || break          # one valid time per file → first msg is enough
        end
    end
    return out
end

# ── NetCDF: model-aware time extraction ───────────────────────────────────────
function _netcdf_validtimes(path::AbstractString; model::AbstractString="")
    v = NetCDF.open(path)
    M = uppercase(model)
    if M == "WRF" && haskey(v.vars, "Times")
        T = v.vars["Times"]
        return [DateTime(join(T[:, k]), dateformat"yyyy-mm-dd_HH:MM:SS") for k in 1:size(T, 2)]
    elseif M == "MERRA2" && haskey(v.vars, "time")
        tv = v.vars["time"]
        m = match(r"since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})", get(tv.atts, "units", ""))
        m === nothing && return DateTime[]
        base = DateTime(m.captures[1] * "T" * m.captures[2])
        return [base + Minute(round(Int, t)) for t in tv[:]]
    elseif M == "MESONH" && haskey(v.vars, "DTCUR")
        # seconds since the segment start time (DTCUR units carry the origin)
        u = get(v.vars["DTCUR"].atts, "units", "")
        m = match(r"since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})", u)
        base = m === nothing ? DateTime(2000) : DateTime(m.captures[1] * "T" * m.captures[2])
        return [base + Second(round(Int, v.vars["DTCUR"][1]))]
    elseif haskey(v.vars, "time")
        tv = v.vars["time"]
        m = match(r"(seconds|minutes|hours|days) since (\d{4}-\d{2}-\d{2})[ T]?(\d{2}:\d{2}:\d{2})?", get(tv.atts, "units", ""))
        if m !== nothing
            base = DateTime(m.captures[2] * "T" * (m.captures[3] === nothing ? "00:00:00" : m.captures[3]))
            step = Dict("seconds"=>Second, "minutes"=>Minute, "hours"=>Hour, "days"=>Day)[m.captures[1]]
            return [base + step(round(Int, t)) for t in tv[:]]
        end
    end
    return DateTime[]
end

# ── Filename fallback (operational naming conventions) ────────────────────────
function _date_from_filename(name::AbstractString)
    try
        # ERA5 IFS_*_YYYYMMDD_<cycleHHMM>_<step>.grb  → base + step hours
        m = match(r"_(\d{8})_(\d{3,4})_(\d{1,2})\.gr", name)
        if m !== nothing
            base = DateTime(m.captures[1], dateformat"yyyymmdd")
            cyc  = parse(Int, m.captures[2]); step = parse(Int, m.captures[3])
            return base + Hour(div(cyc, 100)) + Minute(mod(cyc, 100)) + Hour(step)
        end
        # AROME …YYYY-MM-DDTHH-MM-SSZ…  (base time; steps via the header path)
        m = match(r"(\d{4}-\d{2}-\d{2})T(\d{2})-(\d{2})-(\d{2})Z", name)
        m !== nothing && return DateTime(m.captures[1] * "T" * m.captures[2] * ":" * m.captures[3] * ":" * m.captures[4])
        # Plain YYYYMMDD[HH] (require a non-digit or end on both sides so we
        # don't grab a version token like MERRA2_400).
        m = match(r"(?<!\d)(\d{4})(\d{2})(\d{2})(\d{2})?(?!\d)", name)
        if m !== nothing
            y,mo,da = parse.(Int, (m.captures[1], m.captures[2], m.captures[3]))
            h = m.captures[4] === nothing ? 0 : parse(Int, m.captures[4])
            return DateTime(y, mo, da, h)
        end
    catch
        return nothing
    end
    return nothing
end

"""
    write_available!(path, entries, variablesfile)

Write a FLEXPART-style AVAILABLE file (line 1 = variables namelist) from the
scanned entries.  Optional: build_run keeps the mapping in memory instead.
"""
function write_available!(path::AbstractString, entries, variablesfile::AbstractString,
                          metdir::AbstractString)
    # readavailable! resolves each entry as  metdir * entry, so write the path
    # RELATIVE to the met directory (handles nested layouts like GFS date/cycle).
    open(path, "w") do io
        println(io, variablesfile)
        println(io, "YYYYMMDD HHMISS")
        println(io, "________ ______      __________      __________")
        for (d, f) in entries
            rel = relpath(f, metdir)
            println(io, Dates.format(d, "yyyymmdd"), " ", Dates.format(d, "HHMMSS"),
                    "      ", rel, "  ON DISC")
        end
    end
    return path
end
