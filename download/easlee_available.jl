#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# easlee_available.jl — write the EASLEE AVAILABLE index of a meteorology
# directory from the file names.
#
# `[meteorology].available = "auto"` otherwise scans the directory and opens
# every file to read its valid time, which costs minutes per hundred GRIB
# files.  The download scripts here call this after each download so the index
# is ready before the first run, and the model reuses it as long as it still
# matches what is on disk.
#
#     julia --project=easlee download/easlee_available.jl --model GFS    --dir ./data/GFS
#     julia --project=easlee download/easlee_available.jl --model MERRA2 --dir ./data/MERRA2
#     julia --project=easlee download/easlee_available.jl --model ERA5   --dir ./data/ERA5_ETEX --vars ERA5_CAPTEX
#
# Format (FLEXPART lineage, read by initialize.jl `readavailable!`):
#     line 1   absolute path of the variables namelist (listvariables/<MODEL>.vars)
#     line 2-3 header
#     then     YYYYMMDD HHMMSS      <path relative to the met dir>  ON DISC
# One line per VALID TIME; a file holding several times (MERRA-2) appears once
# per time.  When two files give the same valid time the shorter forecast lead
# wins, as the model's own scanner does.
# ─────────────────────────────────────────────────────────────────────────────

using Dates
using Printf

include(joinpath(@__DIR__, "cli.jl"))

const HERE       = @__DIR__
const REPO       = abspath(joinpath(HERE, ".."))
const INDEX_NAME = "AVAILABLE"

"listvariables/<name>.vars of this EASLEE tree, or an explicit path."
function vars_path(name::AbstractString)
    isfile(name) && return abspath(name)
    fname = endswith(name, ".vars") ? name : name * ".vars"
    cand  = joinpath(REPO, "listvariables", fname)
    isfile(cand) || error("variables namelist not found: $name (looked for $cand)")
    return cand
end

"""
    main_tokens(varsfile) -> Vector{String}

`companions = A>B, C>D` in the namelist gives the MAIN-file tokens [A, C], the
same way the model's scanner reads it.  Empty when no companions are declared.
"""
function main_tokens(varsfile::AbstractString)
    toks = String[]
    isfile(varsfile) || return toks
    for line in eachline(varsfile)
        m = match(r"^\s*companions\s*=\s*(.+?)\s*(?:#.*)?$"i, line)
        m === nothing && continue
        for pat in split(m.captures[1], ",")
            parts = split(pat, ">")
            length(parts) == 2 && push!(toks, strip(parts[1]))
        end
    end
    return toks
end

# ── per-model (valid time, lead hours, relative path) ────────────────────────

const GFS_PAT = r"^gfs\.t(\d{2})z\.pgrb2\.(0p25|0p50|1p00)\.f(\d{3})$"

"<dir>/YYYYMMDD/HH/gfs.tHHz.pgrb2.0p25.fSSS  ->  valid = date + HH + SSS h."
function gfs_entries(metdir::AbstractString)
    out = Tuple{DateTime,Int,String}[]
    for (root, _, files) in walkdir(metdir), f in files
        m = match(GFS_PAT, f)
        m === nothing && continue
        rel   = relpath(joinpath(root, f), metdir)
        parts = splitpath(rel)
        idx   = findfirst(p -> occursin(r"^\d{8}$", p), parts)
        idx === nothing && continue
        base = DateTime(parts[idx], dateformat"yyyymmdd") + Hour(parse(Int, m.captures[1]))
        step = parse(Int, m.captures[3])
        push!(out, (base + Hour(step), step, rel))
    end
    return out
end

const MERRA_PAT = r"^MERRA2_\d{3}\.(inst\d|tavg\d)_(\d)d_([a-z]+)_N([a-z])\.(\d{8})\.nc4$"

"""
    merra2_entries(metdir, main_token)

One entry per valid time of each MAIN file; the companion collections are read
through the main file (`listvariables/MERRA2.vars` `companions`).  Times come
from the file itself when NetCDF can read it, else from the collection
convention (tavgN carries an N/2 h offset).
"""
function merra2_entries(metdir::AbstractString, main_token::AbstractString="_asm_Nv.")
    out = Tuple{DateTime,Int,String}[]
    for (root, _, files) in walkdir(metdir), f in sort(files)
        m = match(MERRA_PAT, f)
        (m === nothing || !occursin(main_token, f)) && continue
        path = joinpath(root, f)
        rel  = relpath(path, metdir)
        day  = DateTime(m.captures[5], dateformat"yyyymmdd")
        times = merra2_times_from_file(path)
        if times === nothing
            kind = m.captures[1]
            n    = parse(Int, kind[end:end])
            off  = startswith(kind, "tavg") ? Minute(30n) : Minute(0)
            times = [day + off + Hour(n * k) for k in 0:(24 ÷ n - 1)]
        end
        append!(out, ((t, 0, rel) for t in times))
    end
    return out
end

"Valid times from a MERRA-2 file's `time` axis, or `nothing` if it cannot be read."
function merra2_times_from_file(path::AbstractString)
    try
        @eval Main using NetCDF
        units = Base.invokelatest(Main.NetCDF.ncgetatt, path, "time", "units")
        u = match(r"since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})", String(units))
        u === nothing && return nothing
        base = DateTime(u.captures[1] * "T" * u.captures[2], dateformat"yyyy-mm-ddTHH:MM:SS")
        tv   = Base.invokelatest(Main.NetCDF.ncread, path, "time")
        return [base + Minute(round(Int, t)) for t in tv]
    catch
        return nothing
    end
end

const ERA5_PAT     = r"^IFS_[A-Z0-9]+(?:_[A-Z]+)?_(\d{8})_(\d{3,4})_(\d{1,2})\.grb$"
const ERA5_EXCLUDE = (r"IFS_PT_MET", r"FULL_PT", r"^IFS_SFC_")

"""
    era5_entries(metdir, tokens)

`IFS_EA_MET_YYYYMMDD_HHMM_S.grb` or `IFS_EA_YYYYMMDD_HHMM_S.grb` -> base + S h,
the base hour written without a leading zero.  Only MAIN files are listed; the
`_CLD_` / `_SFC_` companions come with them.
"""
function era5_entries(metdir::AbstractString, tokens=String[])
    out = Tuple{DateTime,Int,String}[]
    for (root, _, files) in walkdir(metdir), f in files
        any(x -> occursin(x, f), ERA5_EXCLUDE) && continue
        (!isempty(tokens) && !any(t -> occursin(t, f), tokens)) && continue
        m = match(ERA5_PAT, f)
        m === nothing && continue
        hhmm = m.captures[2]
        hh   = length(hhmm) > 2 ? parse(Int, hhmm[1:end-2]) : 0
        mi   = parse(Int, hhmm[end-1:end])
        base = DateTime(m.captures[1], dateformat"yyyymmdd") + Hour(hh) + Minute(mi)
        step = parse(Int, m.captures[3])
        push!(out, (base + Hour(step), step, relpath(joinpath(root, f), metdir)))
    end
    return out
end

const MODELS       = ("GFS", "MERRA2", "ERA5")
const DEFAULT_VARS = Dict("GFS" => "GFS", "MERRA2" => "MERRA2", "ERA5" => "ERA5_CONV")

"""
    write_available(metdir, model; varsfile=nothing, out=nothing, quiet=false)

Build and write `<metdir>/AVAILABLE` plus its `.files` sidecar.
Returns `(path, n_times, first, last)`.
"""
function write_available(metdir::AbstractString, model::AbstractString;
                         varsfile=nothing, out=nothing, quiet::Bool=false)
    model = uppercase(model)
    model in MODELS || error("unknown model $model; known: $(join(MODELS, ", "))")
    vp   = vars_path(something(varsfile, DEFAULT_VARS[model]))
    toks = main_tokens(vp)

    entries = if model == "ERA5"
        era5_entries(metdir, toks)
    elseif model == "MERRA2"
        merra2_entries(metdir, isempty(toks) ? "_asm_Nv." : toks[1])
    else
        gfs_entries(metdir)
    end
    if isempty(entries)
        hint = isempty(toks) ? "" :
               "  The namelist $(basename(vp)) selects main files by $(join(toks, ", ")); " *
               "if these files use a different token, pass the matching --vars."
        error("no $model files found under $metdir$hint")
    end

    # shorter forecast lead wins on a tie, as the model's scanner does
    best = Dict{DateTime,Tuple{Int,String}}()
    for (valid, step, rel) in entries
        if !haskey(best, valid) || step < best[valid][1]
            best[valid] = (step, rel)
        end
    end

    outpath = something(out, joinpath(metdir, INDEX_NAME))
    ts = sort!(collect(keys(best)))
    open(outpath, "w") do f
        println(f, vp)
        println(f, "YYYYMMDD HHMISS")
        println(f, "________ ______      __________      __________")
        for valid in ts
            @printf(f, "%s %s      %s  ON DISC\n",
                    Dates.format(valid, "yyyymmdd"), Dates.format(valid, "HHMMSS"),
                    best[valid][2])
        end
    end

    # sidecar <index>.files: the MAIN files (relative path, size) — what the
    # model checks the index against before trusting it
    mains = sort!(unique(rel for (_, _, rel) in entries))
    open(outpath * ".files", "w") do f
        for rel in mains
            println(f, rel, "\t", filesize(joinpath(metdir, rel)))
        end
    end

    if !quiet
        println("$outpath: $(length(ts)) valid times, $(ts[1]) -> $(ts[end])  (namelist $(basename(vp)))")
        if length(ts) > 2
            base = ts[2] - ts[1]
            gaps = [(a, b) for (a, b) in zip(ts, ts[2:end]) if (b - a) > base]
            isempty(gaps) ||
                println("  note: $(length(gaps)) gap(s) longer than the first interval, e.g. $(gaps[1][1]) -> $(gaps[1][2])")
        end
    end
    return outpath, length(ts), ts[1], ts[end]
end

# ── CLI ──────────────────────────────────────────────────────────────────────

function available_main(args=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        usage: julia --project=easlee download/easlee_available.jl --model MODEL --dir DIR [--vars NAME] [--out PATH]

          --model   $(join(MODELS, " | "))
          --dir     meteorology directory
          --vars    variables namelist name or path (default per model)
          --out     index path (default <dir>/AVAILABLE)
        """)
        return 0
    end
    _, o = parse_flags(args, Dict{String,Any}("model" => nothing, "dir" => nothing,
                                              "vars" => nothing, "out" => nothing))
    o["model"] === nothing && error("--model is required")
    o["dir"]   === nothing && error("--dir is required")
    write_available(o["dir"], o["model"]; varsfile=o["vars"], out=o["out"])
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(available_main())
end
