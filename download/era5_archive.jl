#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# era5_archive.jl — fetch and unpack a prepared ERA5 archive of one of the
# validation cases, then write its AVAILABLE index.
#
#     julia --project=easlee download/era5_archive.jl --dataset etex
#     julia --project=easlee download/era5_archive.jl --dataset captex
#     julia --project=easlee download/era5_archive.jl --url https://.../ERA5_ETEX.tar
#
# Download links live in links.toml next to this script; fill them in once the
# archives are hosted.  General ERA5 retrievals are not scripted here: EASLEE
# reads flex_extract-style GRIB (IFS_EA_MET_* + IFS_EA_CLD_*), see
# https://www.flexpart.eu/flex_extract.  This script only unpacks a prepared
# archive.
#
# `--dataset etex` unpacks into data/ERA5_ETEX and indexes it with the
# ERA5_CAPTEX namelist, which is what the configs/etex_* and configs/validate_*
# files expect, so a fresh clone needs no path editing.
# ─────────────────────────────────────────────────────────────────────────────

using Downloads
using Tar
using p7zip_jll

include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "easlee_available.jl"))

const LINKS = joinpath(@__DIR__, "links.toml")

# The prepared archives are flex_extract IFS_EA_* files with IFS_SFC_*
# companions, which is what ERA5_CAPTEX.vars declares; ERA5_CONV.vars expects
# _MET_/_CLD_ and would write an index that matches nothing.
const ARCHIVE_VARS = Dict("etex" => "ERA5_CAPTEX", "captex" => "ERA5_CAPTEX")

"Read links.toml without pulling in a TOML parser for three lines of `key = \"value\"`."
function read_links()
    links = Dict{String,String}()
    isfile(LINKS) || return links
    for line in eachline(LINKS)
        line = strip(split(line, "#")[1])
        occursin("=", line) || continue
        k, v = split(line, "=", limit=2)
        links[strip(k)] = strip(strip(v), ['"', '\''])
    end
    return links
end

"""
    unpack(archive, out)

Unpack `.tar`, `.tar.gz`/`.tgz` or `.zip` into `out`.  Plain tar goes through
the `Tar` standard library; the compressed forms go through p7zip, which ships
with EASLEE, so nothing has to be installed on the system.
"""
function unpack(archive::AbstractString, out::AbstractString)
    low = lowercase(archive)
    # Tar.extract insists on an empty destination, and `out` already holds the
    # downloaded archive, so everything is staged next door and moved in.  The
    # staging directory is a sibling of `out` to keep the move on one filesystem.
    stage = mktempdir(dirname(abspath(out)))
    try
        if endswith(low, ".tar")
            rm(stage; recursive=true)           # Tar.extract creates it itself
            Tar.extract(archive, stage)
        elseif endswith(low, ".zip")
            run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$stage $archive`, stdout=devnull))
        elseif endswith(low, ".tar.gz") || endswith(low, ".tgz")
            mktempdir() do tmp                  # decompress, then untar
                run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$tmp $archive`, stdout=devnull))
                inner = only(readdir(tmp; join=true))
                rm(stage; recursive=true)
                Tar.extract(inner, stage)
            end
        else
            error("unknown archive type: $archive (expected .tar, .tar.gz, .tgz or .zip)")
        end
        for name in readdir(stage)
            mv(joinpath(stage, name), joinpath(out, name); force=true)
        end
    finally
        rm(stage; recursive=true, force=true)
    end
end

"""
    strip_top_level(out, keep) -> String or nothing

If the archive wrapped everything in one directory, lift its contents up.
`tar -cf a.tar DIR` stores `DIR/...` while `tar -cf a.tar -C DIR .` does not,
and the meteorology should sit directly under `out` either way.  `keep` is the
downloaded archive, which must not be mistaken for the wrapper.
"""
function strip_top_level(out::AbstractString, keep::AbstractString)
    entries = filter(e -> joinpath(out, e) != keep, readdir(out))
    length(entries) == 1 || return nothing
    inner = joinpath(out, entries[1])
    isdir(inner) || return nothing
    for name in readdir(inner)
        mv(joinpath(inner, name), joinpath(out, name); force=true)
    end
    rm(inner)
    return entries[1]
end

function main(args=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        usage: julia --project=easlee download/era5_archive.jl (--dataset etex|captex | --url URL) [options]

          --dataset D     named archive from links.toml (etex or captex)
          --url U         archive URL (.tar / .tar.gz / .tgz / .zip), overrides --dataset
          --out DIR       target directory (default ./data/ERA5_<DATASET>)
          --vars NAME     variables namelist for the index (default ERA5_CAPTEX for these archives)
          --keep-archive  do not delete the downloaded archive after unpacking
        """)
        return 0
    end

    _, o = parse_flags(args, Dict{String,Any}(
        "dataset" => nothing, "url" => nothing, "out" => nothing,
        "vars" => nothing, "keep-archive" => false))

    ds = o["dataset"]
    ds === nothing || ds in ("etex", "captex") || error("--dataset must be etex or captex")
    url = something(o["url"], get(read_links(), "$(ds)_era5", ""), "")
    isempty(url) &&
        error("no URL: give --url, or --dataset etex|captex with the link filled in download/links.toml")

    vars_ = something(o["vars"], get(ARCHIVE_VARS, something(ds, ""), "ERA5_CONV"))
    out   = something(o["out"], "./data/ERA5_$(uppercase(something(ds, "archive")))")
    mkpath(out)

    base = basename(split(url, "?")[1])
    arc  = joinpath(out, isempty(base) ? "era5_archive.tar" : base)
    if !isfile(arc)
        println("get $url")
        Downloads.download(url, arc; progress=progress_printer())
        println()
    else
        println("exists: $arc")
    end

    println("unpacking into $out")
    unpack(arc, out)
    wrapper = strip_top_level(out, arc)
    wrapper === nothing || println("lifted the contents of $wrapper/ up into $out")
    o["keep-archive"] || rm(arc; force=true)

    write_available(out, "ERA5"; varsfile=vars_)
    println("[[meteorology]]  model = \"ERA5\"  path = \"$out\"  variables = \"$vars_\"  available = \"auto\"")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
