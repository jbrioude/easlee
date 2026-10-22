
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
                  uppercase(model) == "AROME" ? _arome_validinit(f) :
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

# Météo-France AROME NetCDF (AR_T<HHMM>P_<yyyymmdd>_<hh>.nc): `date` (yyyymmdd)
# is the RUN date, T<HHMM>P the run hour, `forecast` (= <hh>) the LEAD hour:
# valid = date + run hour + lead.  The archive holds leads 06-17 of the 00 and
# 12 UTC runs (which tile the day) and leads 06-11 of the 06 and 18 UTC runs;
# overlapping valid times resolve to the freshest run like any other cycle.
# Checked on La Réunion 2024-06-29/30: under this reading land T2m is
# continuous across runs (steps ≤ 2.7 K, the morning rise), while reading
# `forecast` as the valid hour puts four runs' "10 UTC" 4.8 K apart.
# History: until 2026-10-13 every file was read as <date> 00:00 (a whole day on
# one valid time); a first fix the same day read `forecast` as the valid hour,
# which was wrong for every run but 00 UTC.
function _arome_validinit(path::AbstractString)
    v = NetCDF.open(path)
    g = v.gatts
    (haskey(g, "date") && haskey(g, "forecast")) || return Tuple{DateTime,DateTime}[]
    day  = DateTime(strip(string(g["date"])), dateformat"yyyymmdd")
    m    = match(r"T(\d{2})(\d{2})P", basename(path))
    init = m === nothing ? day : day + Hour(parse(Int, m.captures[1])) + Minute(parse(Int, m.captures[2]))
    return [(init + Hour(parse(Int, strip(string(g["forecast"])))), init)]
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
        # AROME AR_T<HHMM>P_<yyyymmdd>_<hh>.nc: run date + run hour + lead hour
        m = match(r"T(\d{2})(\d{2})P_(\d{8})_(\d{2})\.\w+$", name)
        m !== nothing && return DateTime(m.captures[3], dateformat"yyyymmdd") + Hour(parse(Int, m.captures[1])) +
                                Minute(parse(Int, m.captures[2])) + Hour(parse(Int, m.captures[4]))
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

# ── Pre-generated AVAILABLE in the met directory (2026-09-24) ─────────────────────
# `[meteorology].available = "auto"` first looks for <path>/AVAILABLE — written by the
# download tools (download/*.jl) or by a previous run — and uses it when it still
# describes what is on disk; otherwise the directory is scanned and the index written
# there for the next run.  "scan" forces a rescan.  Line 1 (the variables namelist)
# is rewritten with the run's own namelist when the index is copied into the work dir.
const AVAILABLE_INDEX_NAME = "AVAILABLE"

# The MAIN met files of a directory (the same filter as scan_met_directory), as paths
# relative to `directory`.
function list_main_met_files(directory::AbstractString; include_only::Vector{Regex}=Regex[],
                             exclude::Vector{Regex}=Regex[], recursive::Bool=true)
    files = String[]
    if recursive
        for (root, _, fs) in walkdir(directory), f in fs
            push!(files, joinpath(root, f))
        end
    else
        append!(files, readdir(directory; join=true))
    end
    out = String[]
    for f in files
        b = basename(f)
        isfile(f) || continue
        occursin(r"AVAILABLE|\.vars$|\.txt$|README|\.pdf$|\.part$|\.toml$"i, b) && continue
        isempty(include_only) || any(p -> occursin(p, b), include_only) || continue
        any(p -> occursin(p, b), exclude) && continue
        push!(out, relpath(f, directory))
    end
    return sort!(out)
end

# (date, relative path) entries of an AVAILABLE file (line 1 = namelist, 2-3 = header).
function read_available_entries(path::AbstractString)
    entries = Tuple{DateTime,String}[]
    for (k, line) in enumerate(eachline(path))
        k <= 3 && continue
        s = split(line)
        length(s) >= 3 || continue
        m = match(r"^\d{8}$", s[1]); m === nothing && continue
        push!(entries, (DateTime(s[1] * s[2], dateformat"yyyymmddHHMMSS"), String(s[3])))
    end
    return entries
end

# Sidecar <index>.files: the MAIN files (relative path, size) the index was built from.
# Written next to the index by the model and by download/easlee_available.jl; it is what
# makes the validity check exact for archives whose overlapping cycles leave some files
# out of the index (the scanner keeps the freshest cycle per valid time).
sidecar_path(index::AbstractString) = index * ".files"
# Version of the scanner that wrote an index, stamped as the sidecar's first line
# ("# easlee-scan N", no tab: readers of the file list skip it).  3 = AROME valid
# times = run date + run hour + lead (2026-10-13; version 2, the same day, read
# the lead as the valid hour).  An older AROME index (no stamp, or 2) would
# otherwise be reused forever, since its files and sizes still match — it is
# rescanned instead.
const SCAN_VERSION = 3
const _SCAN_VERSION_REQUIRED = Dict("AROME" => 3)
function _sidecar_scan_version(path::AbstractString)
    isfile(path) || return 0
    for l in eachline(path)
        m = match(r"^# easlee-scan (\d+)", l)
        m === nothing || return parse(Int, m.captures[1])
    end
    return 1
end
function write_available_sidecar!(index::AbstractString, metdir::AbstractString;
                                  include_only::Vector{Regex}=Regex[], exclude::Vector{Regex}=Regex[])
    open(sidecar_path(index), "w") do io
        println(io, "# easlee-scan ", SCAN_VERSION)
        for f in list_main_met_files(metdir; include_only=include_only, exclude=exclude)
            println(io, f, "\t", filesize(joinpath(metdir, f)))
        end
    end
end
_read_sidecar(path) = [(String(split(l, '\t')[1]), parse(Int, split(l, '\t')[2])) for l in eachline(path) if occursin('\t', l)]

# Does the pre-generated index still describe the directory?  Every listed file must
# exist and be non-empty; then, with a sidecar, the MAIN files on disk (names and sizes)
# must be exactly those recorded; without one, no main file may be newer than the index.
# Returns (ok, reason).
function cached_available_ok(index::AbstractString, metdir::AbstractString;
                             include_only::Vector{Regex}=Regex[], exclude::Vector{Regex}=Regex[],
                             model::AbstractString="")
    isfile(index) || return false, "no $(basename(index)) in $metdir"
    need = get(_SCAN_VERSION_REQUIRED, uppercase(model), 0)
    if need > 0
        v = _sidecar_scan_version(sidecar_path(index))
        v >= need || return false, "index written by an older scanner (scan version $(v == 0 ? "unknown, no sidecar" : v); $(uppercase(model)) needs $need)"
    end
    entries = read_available_entries(index)
    isempty(entries) && return false, "$(basename(index)) lists no files"
    listed = unique!(sort!([normpath(f) for (_, f) in entries]))
    for f in listed
        p = joinpath(metdir, f)
        isfile(p) || return false, "listed file missing: $f"
        filesize(p) > 0 || return false, "listed file is empty: $f"
    end
    ondisk = list_main_met_files(metdir; include_only=include_only, exclude=exclude)
    sc = sidecar_path(index)
    if isfile(sc)
        rec = _read_sidecar(sc)
        now = [(f, filesize(joinpath(metdir, f))) for f in ondisk]
        if now != rec
            added = setdiff(first.(now), first.(rec)); gone = setdiff(first.(rec), first.(now))
            why = !isempty(added) ? "$(length(added)) file(s) added (e.g. $(first(added)))" :
                  !isempty(gone)  ? "$(length(gone)) file(s) removed (e.g. $(first(gone)))" :
                                    "file sizes changed"
            return false, why
        end
    else
        tidx = mtime(index)
        newer = filter(f -> mtime(joinpath(metdir, f)) > tidx, ondisk)
        isempty(newer) || return false, "$(length(newer)) file(s) newer than the index (e.g. $(first(newer))); no $(basename(sc)) sidecar"
    end
    return true, "$(length(entries)) time steps, $(first(entries)[1]) → $(last(entries)[1])"
end

# Copy an index into the work dir with line 1 replaced by the run's variables namelist.
function copy_available_with_vars(src::AbstractString, dst::AbstractString, varsfile::AbstractString)
    lines = readlines(src)
    isempty(lines) && error("empty AVAILABLE: $src")
    lines[1] = String(varsfile)
    open(dst, "w") do io
        foreach(l -> println(io, l), lines)
    end
    return dst
end
