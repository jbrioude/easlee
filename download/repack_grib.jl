#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# repack_grib.jl — shrink a GRIB meteorological archive with 12-bit JPEG packing.
#
# ERA5 archives as retrieved are packed with `grid_simple` at 16 bits per value.
# Repacking them to `grid_jpeg` at 12 bits makes them roughly 3x smaller on
# disk, which is what turns a 33 GB archive into something one can publish.  The
# cost is a quantisation of each field to 4096 levels across its own range, i.e.
# a maximum absolute error of range/4096: about 0.004 m/s on a wind field and
# 0.004 K on a temperature field.  EASLEE scores are unchanged by it
# (docs/BENCHMARKS.md §0).
#
#     julia --project=easlee download/repack_grib.jl /data/ETEX --out ./data/ERA5_ETEX
#     julia --project=easlee download/repack_grib.jl /data/ETEX --out ./out --jobs 8
#     julia --project=easlee download/repack_grib.jl /data/ETEX --dry-run
#
# The directory tree is mirrored: every GRIB file is repacked, everything else
# (AVAILABLE, station lists, notes) is copied unchanged.  The source is never
# modified.
#
# Messages that cannot be JPEG-packed — constant and near-constant fields, which
# the JPEG coder rejects — keep their original packing.  That is expected and is
# why the saving is ~3x rather than the 16/12 ratio.  Every output file is
# checked to hold exactly as many GRIB messages as its input; if one does not,
# the original is copied instead and the file is reported.
#
# The ecCodes tools come from `eccodes_jll`, the same artifact GRIB.jl reads
# with, so nothing has to be installed on the system.
# ─────────────────────────────────────────────────────────────────────────────

using Printf
using eccodes_jll

include(joinpath(@__DIR__, "cli.jl"))

"True if the file starts with the GRIB magic, whatever it is called."
function is_grib(path::AbstractString)
    try
        open(path, "r") do io
            return read(io, 4) == b"GRIB"
        end
    catch
        return false
    end
end

"Number of GRIB messages in a file, or -1 if it cannot be counted."
function count_messages(path::AbstractString)
    try
        out = readchomp(`$(grib_count()) $path`)
        return parse(Int, strip(out))
    catch
        return -1
    end
end

"""
    repack_one(src, dst, bits) -> (src, ok, note, size_in, size_out)

Repack one GRIB file, falling back to a plain copy whenever the result would
not be a faithful, smaller replacement.
"""
function repack_one(src::AbstractString, dst::AbstractString, bits::Int)
    mkpath(dirname(dst))
    size_in = filesize(src)
    err = IOBuffer()
    # -f keeps going past messages the JPEG coder rejects, copying them as they
    # are; without it grib_set stops and silently writes a truncated file.
    cmd = `$(grib_set()) -f -r -s packingType=grid_jpeg,bitsPerValue=$bits $src $dst`
    rc = try
        run(pipeline(cmd, stdout=devnull, stderr=err))
        0
    catch e
        e isa ProcessFailedException ? 1 : rethrow()
    end

    if !isfile(dst)
        cp(src, dst; force=true)
        return (src, false, "grib_set produced nothing ($rc)", size_in, size_in)
    end

    n_in, n_out = count_messages(src), count_messages(dst)
    if n_in < 0 || n_out < 0 || n_in != n_out
        cp(src, dst; force=true)
        return (src, false, "message count $n_in -> $n_out, original kept", size_in, size_in)
    end

    size_out = filesize(dst)
    if size_out >= size_in            # already packed, or packing made it worse
        cp(src, dst; force=true)
        return (src, true, "no gain, original kept", size_in, size_in)
    end

    skipped = count(x -> occursin("ECCODES ERROR", x), split(String(take!(err)), '\n'))
    note = skipped > 0 ? "$skipped message(s) left as-is" : ""
    return (src, true, note, size_in, size_out)
end

function main(args=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        usage: julia --project=easlee download/repack_grib.jl SRC [--out DIR] [--bits N] [--jobs N] [--dry-run]

          SRC         meteorological archive directory to repack
          --out DIR   target directory (required unless --dry-run); the source is never modified
          --bits N    bits per value (default 12)
          --jobs N    files to repack in parallel (default: half the cores)
          --dry-run   report what would be done, change nothing
        """)
        return 0
    end

    positional, o = parse_flags(args, Dict{String,Any}(
        "out" => nothing, "bits" => "12",
        "jobs" => string(max(1, Sys.CPU_THREADS ÷ 2)), "dry-run" => false))

    length(positional) == 1 || error("give exactly one source directory")
    src  = positional[1]
    bits = parse(Int, o["bits"])
    jobs = parse(Int, o["jobs"])
    isdir(src) || error("not a directory: $src")
    (o["out"] === nothing && !o["dry-run"]) &&
        error("give --out <dir> (the source is never modified), or --dry-run")

    gribs, others = String[], String[]
    for (root, _, files) in walkdir(src), f in sort(files)
        p = joinpath(root, f)
        push!(is_grib(p) ? gribs : others, p)
    end
    isempty(gribs) && error("no GRIB files found under $src")

    total_in = sum(filesize, gribs)
    @printf("%d GRIB files (%.2f GB), %d other files\n", length(gribs), total_in/1e9, length(others))
    if o["dry-run"]
        @printf("would repack to grid_jpeg / %d bits, expect roughly %.2f GB (~3x smaller)\n",
                bits, total_in/3e9)
        return 0
    end

    out = o["out"]
    rel(p) = relpath(p, src)

    done_in = Threads.Atomic{Int}(0)
    done_out = Threads.Atomic{Int}(0)
    ndone = Threads.Atomic{Int}(0)
    failures = Tuple{String,String}[]
    flock = ReentrantLock()

    # One worker task per job, each pulling from a shared channel: grib_set is
    # an external process, so this parallelises even on a single Julia thread.
    work = Channel{String}(length(gribs))
    foreach(p -> put!(work, p), gribs)
    close(work)

    @sync for _ in 1:jobs
        Threads.@spawn for p in work
            s, ok, note, s_in, s_out = repack_one(p, joinpath(out, rel(p)), bits)
            Threads.atomic_add!(done_in, s_in)
            Threads.atomic_add!(done_out, s_out)
            n = Threads.atomic_add!(ndone, 1) + 1
            ok || lock(flock) do
                push!(failures, (rel(s), note))
            end
            @printf("\r  %d/%d  %.2f -> %.2f GB", n, length(gribs),
                    done_in[]/1e9, done_out[]/1e9)
            flush(stdout)
        end
    end
    println()

    for p in others
        dst = joinpath(out, rel(p))
        mkpath(dirname(dst))
        cp(p, dst; force=true)
    end

    ratio = done_out[] > 0 ? done_in[]/done_out[] : 1.0
    println("\n$src -> $out")
    @printf("  %.2f GB -> %.2f GB  (%.2fx smaller)\n", done_in[]/1e9, done_out[]/1e9, ratio)
    println("  $(length(others)) non-GRIB file(s) copied unchanged")
    if !isempty(failures)
        println("  $(length(failures)) file(s) kept at their original packing:")
        for (name, note) in first(failures, 10)
            println("    $name: $note")
        end
        length(failures) > 10 && println("    ... and $(length(failures)-10) more")
    end
    println("\nCheck it before publishing:")
    println("  julia --project=easlee scripts/easlee.jl <a config pointed at $out>")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
