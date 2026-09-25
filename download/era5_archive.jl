#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# era5_archive.jl — fetch and unpack a prepared ERA5 archive of one of the
# validation cases, then write its AVAILABLE index.
#
#     julia --project=easlee download/era5_archive.jl --dataset etex
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
#
# Only the ETEX-1 meteorology is distributed.  The CAPTEX-83 case is supported by
# the model and its configs, but its ERA5 archive is three times the size and is
# not hosted; retrieve it yourself with flex_extract and point
# `[[meteorology]] path` at it, or pass `--url` if you have it somewhere.
# ─────────────────────────────────────────────────────────────────────────────

using Downloads
using Tar
using p7zip_jll

include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "easlee_available.jl"))

const LINKS = joinpath(@__DIR__, "links.toml")

# The prepared archive holds flex_extract IFS_EA_* files with IFS_SFC_*
# companions, which is what ERA5_CAPTEX.vars declares (the namelist is named for
# the campaign it was first written for, and serves the ETEX archive too);
# ERA5_CONV.vars expects _MET_/_CLD_ and would write an index matching nothing.
const ARCHIVE_VARS = Dict("etex" => "ERA5_CAPTEX")

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
    resolve_archive_url(url) -> String

A DOI, or a Zenodo record page, points at a landing page rather than a file.
Given either, ask Zenodo which files the record holds and return the direct
download link (preferring a .tar).  Anything else is passed through unchanged,
so a plain file URL still works.

The reply is scanned with a regex rather than parsed as JSON: it saves EASLEE a
dependency for the two fields we need, and a change in Zenodo's schema shows up
as "could not find a file", not as wrong data.
"""
function resolve_archive_url(url::AbstractString)
    m = match(r"zenodo\.org/records?/(\d+)", url)
    m === nothing && (m = match(r"10\.5281/zenodo\.(\d+)"i, url))
    m === nothing && return String(url)
    id = m.captures[1]

    io = IOBuffer()
    try
        Downloads.download("https://zenodo.org/api/records/$id", io)
    catch e
        error("Zenodo record $id could not be read ($(nameof(typeof(e)))).  If the " *
              "deposition is still a draft its DOI is not registered yet; publish it, " *
              "or give the direct file URL with --url.")
    end
    body = String(take!(io))

    links = [x.captures[1] for x in eachmatch(r"\"(https://zenodo\.org/api/records/\d+/files/[^\"]+/content)\"", body)]
    isempty(links) &&
        error("Zenodo record $id lists no downloadable file.  Give the direct file URL with --url.")
    tar = findfirst(l -> occursin(".tar", l), links)
    chosen = links[something(tar, 1)]
    println("resolved $url -> $chosen")
    return chosen
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
        usage: julia --project=easlee download/era5_archive.jl (--dataset etex | --url URL) [options]

          --dataset D     named archive from links.toml (etex; CAPTEX is not distributed)
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
    ds === nothing || ds == "etex" ||
        error("--dataset only knows etex: the CAPTEX meteorology is not distributed " *
              "(retrieve it with flex_extract, or give --url)")
    url = something(o["url"], get(read_links(), "$(ds)_era5", ""), "")
    isempty(url) &&
        error("no URL: give --url, or --dataset etex with the link filled in download/links.toml")

    vars_ = something(o["vars"], get(ARCHIVE_VARS, something(ds, ""), "ERA5_CONV"))
    out   = something(o["out"], "./data/ERA5_$(uppercase(something(ds, "archive")))")
    mkpath(out)

    url  = resolve_archive_url(url)
    base = basename(split(url, "?")[1])
    endswith(base, "/content") && (base = basename(dirname(split(url, "?")[1])))
    arc  = joinpath(out, isempty(base) ? "era5_archive.tar" : base)
    if !isfile(arc)
        println("get $url")
        Downloads.download(url, arc; progress=progress_printer())
        println()
        # A landing page or an error page saved under an archive name fails later
        # with a confusing "unknown archive type"; say what actually happened.
        magic = open(arc, "r") do io; read(io, min(filesize(arc), 512)); end
        if occursin("<html", lowercase(String(copy(magic))))
            rm(arc; force=true)
            error("that URL returned a web page, not an archive.  A DOI or record page " *
                  "is a landing page: give the direct file link (Zenodo: the link behind " *
                  "the Download button), or a DOI this script can resolve.")
        end
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
