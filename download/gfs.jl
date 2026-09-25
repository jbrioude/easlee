#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# gfs.jl — download NCEP GFS 0.25/0.5/1.0° GRIB2 forecasts for EASLEE and write
# the AVAILABLE index.
#
# Two ways to say what you want:
#
#   one analysis and its forecast (what "test it now" needs):
#       julia --project=easlee download/gfs.jl --analysis 2026-06-03T00 --hours 48 --out ./data/GFS
#   a date range of analyses, cycles and steps (a reanalysis-like archive):
#       julia --project=easlee download/gfs.jl --start 2026-06-03 --end 2026-06-08 \
#             --cycles 00,12 --steps 3,6,9,12,15 --out ./data/GFS
#
# Files land in  <out>/YYYYMMDD/HH/gfs.tHHz.pgrb2.<res>.fSSS, the layout
# listvariables/GFS.vars expects.  Steps must be >= 3 h: the surface fluxes
# EASLEE reads are 3-h averages, which the analysis file (f000) does not carry.
#
# Sources: NOMADS (last ~10 days), then the public NOAA archive on AWS (2021
# onward).  No account is needed for either.  Point EASLEE at the result with
#     [[meteorology]]  model = "GFS"  path = "<out>"  variables = "GFS"  available = "auto"
#
# Options: --res 0p25|0p50|1p00 (default 0p25), --source auto|nomads|aws,
#          --dry-run, --no-available.  Complete files already present are skipped.
# ─────────────────────────────────────────────────────────────────────────────

using Dates
using Downloads
using Printf

include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "easlee_available.jl"))

const NOMADS = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
const AWS    = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"

"File name and the candidate URLs for one (day, cycle, step)."
function urls(date::Date, cycle::AbstractString, step::Int, res::AbstractString, source::AbstractString)
    fn   = @sprintf("gfs.t%sz.pgrb2.%s.f%03d", cycle, res, step)
    day  = Dates.format(date, "yyyymmdd")
    cand = String[]
    source in ("auto", "nomads") && push!(cand, "$NOMADS/gfs.$day/$cycle/atmos/$fn")
    source in ("auto", "aws")    && push!(cand, "$AWS/gfs.$day/$cycle/atmos/$fn")
    return fn, cand
end

"Content-Length of a URL, or -1 when the server will not say."
function remote_size(url::AbstractString)
    try
        sz = -1
        Downloads.request(url; method="HEAD", timeout=60,
                          output=devnull) do resp
            for (k, v) in resp.headers
                lowercase(k) == "content-length" && (sz = something(tryparse(Int, v), -1))
            end
        end
        return sz
    catch
        return -1
    end
end

"""
    fetch_file(url, dest; tries=3) -> Bool

Download to `dest` through a `.part` file so an interrupted transfer never
leaves a half file that later looks complete.  `false` on a 404, after
retrying anything else.
"""
function fetch_file(url::AbstractString, dest::AbstractString; tries::Int=3)
    tmp = dest * ".part"
    for k in 1:tries
        try
            Downloads.download(url, tmp; timeout=120, progress=progress_printer())
            println()
            mv(tmp, dest; force=true)
            return true
        catch e
            isfile(tmp) && rm(tmp; force=true)
            if e isa Downloads.RequestError && e.response.status == 404
                return false        # simply not published, no point retrying
            end
            println("  $(nameof(typeof(e))): retry $k/$tries")
            k < tries && sleep(5k)
        end
    end
    return false
end

function main(args=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        usage: julia --project=easlee download/gfs.jl (--analysis T | --start D [--end D]) [options]

          --analysis T    analysis time YYYY-MM-DDTHH (one cycle and its forecast)
          --hours N       forecast length with --analysis (default 48)
          --start D       first analysis day YYYY-MM-DD
          --end D         last analysis day YYYY-MM-DD (default: --start)
          --cycles L      cycles with --start (default 00,12)
          --steps L       forecast steps [h] with --start (default 3,6,9,12,15)
          --res R         0p25 | 0p50 | 1p00 (default 0p25)
          --out DIR       target directory (default ./data/GFS)
          --source S      auto | nomads | aws (default auto)
          --dry-run       print the URLs, download nothing
          --no-available  do not write the AVAILABLE index
        """)
        return 0
    end

    _, o = parse_flags(args, Dict{String,Any}(
        "analysis" => nothing, "hours" => "48", "start" => nothing, "end" => nothing,
        "cycles" => "00,12", "steps" => "3,6,9,12,15", "res" => "0p25",
        "out" => "./data/GFS", "source" => "auto",
        "dry-run" => false, "no-available" => false))

    o["res"] in ("0p25", "0p50", "1p00") || error("--res must be 0p25, 0p50 or 1p00")
    o["source"] in ("auto", "nomads", "aws") || error("--source must be auto, nomads or aws")

    jobs = Tuple{Date,String,Int}[]
    if o["analysis"] !== nothing
        t0  = DateTime(o["analysis"], dateformat"yyyy-mm-ddTHH")
        cyc = @sprintf("%02d", (Dates.hour(t0) ÷ 6) * 6)
        for s in 3:3:parse(Int, o["hours"])
            push!(jobs, (Date(t0), cyc, s))
        end
    elseif o["start"] !== nothing
        d0 = Date(o["start"], dateformat"yyyy-mm-dd")
        d1 = o["end"] === nothing ? d0 : Date(o["end"], dateformat"yyyy-mm-dd")
        steps  = sort!(unique(parse(Int, strip(s)) for s in split(o["steps"], ",") if !isempty(strip(s))))
        cycles = [lpad(strip(c), 2, '0') for c in split(o["cycles"], ",") if !isempty(strip(c))]
        for d in d0:Day(1):d1, c in cycles, s in steps
            push!(jobs, (d, c, s))
        end
    else
        error("give --analysis (one forecast) or --start [--end] (a date range)")
    end
    any(j -> j[3] < 3, jobs) &&
        @warn "steps below 3 h carry no 3-h averaged surface fluxes; EASLEE needs steps >= 3"

    n_ok = n_skip = n_miss = 0
    for (date, cyc, step) in jobs
        fn, cands = urls(date, cyc, step, o["res"], o["source"])
        dest = joinpath(o["out"], Dates.format(date, "yyyymmdd"), cyc, fn)
        if o["dry-run"]
            println(cands[1]); continue
        end
        mkpath(dirname(dest))
        got = false
        for url in cands
            if isfile(dest) && filesize(dest) > 0
                rs = remote_size(url)
                if rs <= 0 || rs == filesize(dest)
                    println("exists: $dest"); n_skip += 1; got = true; break
                end
            end
            println("get $url")
            if fetch_file(url, dest)
                n_ok += 1; got = true; break
            end
        end
        if !got
            println("  missing: $fn ($date $(cyc)z)")
            n_miss += 1
        end
    end
    println("done: $n_ok downloaded, $n_skip already present, $n_miss missing")
    if !o["dry-run"] && !o["no-available"] && (n_ok + n_skip) > 0
        write_available(o["out"], "GFS")
    end
    return n_miss == 0 ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
