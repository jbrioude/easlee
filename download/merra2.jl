#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# merra2.jl — download MERRA-2 3-hourly model-level files from NASA GES DISC
# for EASLEE and write the AVAILABLE index.
#
#     julia --project=easlee download/merra2.jl --start 2022-01-18 --end 2022-01-29 --out ./data/MERRA2
#
# Collections (the default pair is what listvariables/MERRA2.vars reads):
# tavg3_3d_asm_Nv (M2T3NVASM — winds, T, q, omega, pressure, heights, cloud) and
# its companion tavg3_3d_mst_Ne (M2T3NEMST — convective mass flux and
# precipitation fluxes).  `--collections` fetches others (tavg3_3d_trb_Ne,
# tavg3_3d_cld_Nv, inst3_3d_asm_Nv, tavg1_2d_flx_Nx, tavg1_2d_slv_Nx).
#
# Authentication needs a free NASA Earthdata account
# (https://urs.earthdata.nasa.gov), and the account must be authorised for
# "NASA GESDISC DATA ARCHIVE" under Earthdata profile > Applications.  Give it
# either way:
#   * a bearer token, through `--token` or the EARTHDATA_TOKEN environment
#     variable (never written to disk or shown in the log), or
#   * a ~/.netrc line:  machine urs.earthdata.nasa.gov login <user> password <pw>
#     (on Windows, _netrc in your home directory).
#
# Files land in <out>/MERRA2_<stream>.<collection>.YYYYMMDD.nc4.  Point EASLEE
# at the result with
#     [[meteorology]]  model = "MERRA2"  path = "<out>"  variables = "MERRA2"  available = "auto"
# ─────────────────────────────────────────────────────────────────────────────

using Dates
using Downloads

include(joinpath(@__DIR__, "cli.jl"))
include(joinpath(@__DIR__, "easlee_available.jl"))

const BASE = "https://data.gesdisc.earthdata.nasa.gov/data/MERRA2"

const PRODUCT = Dict(
    "tavg3_3d_asm_Nv" => "M2T3NVASM", "tavg3_3d_mst_Ne" => "M2T3NEMST",
    "tavg3_3d_trb_Ne" => "M2T3NETRB", "tavg3_3d_cld_Nv" => "M2T3NVCLD",
    "inst3_3d_asm_Nv" => "M2I3NVASM", "tavg1_2d_flx_Nx" => "M2T1NXFLX",
    "tavg1_2d_slv_Nx" => "M2T1NXSLV", "inst1_2d_asm_Nx" => "M2I1NXASM")

const DEFAULT_COLL = ["tavg3_3d_asm_Nv", "tavg3_3d_mst_Ne"]

"MERRA-2 stream number for a year (the leading 100/200/300/400 of the file name)."
stream(year::Int) = year <= 1991 ? 100 : year <= 2000 ? 200 : year <= 2010 ? 300 : 400

"Path of the user's netrc file, whatever the platform calls it."
netrc_path() = Sys.iswindows() ? joinpath(homedir(), "_netrc") : joinpath(homedir(), ".netrc")

"True when the netrc file holds an entry for the Earthdata login host."
function has_netrc_entry()
    p = netrc_path()
    isfile(p) || return false
    try
        return occursin("urs.earthdata.nasa.gov", read(p, String))
    catch
        return false
    end
end

"""
    make_downloader() -> Downloads.Downloader

GES DISC redirects every data URL through the Earthdata login host, so the
transfer needs a cookie engine and credentials that survive the redirect.
`CURLOPT_NETRC` supplies the login from ~/.netrc, `CURLOPT_UNRESTRICTED_AUTH`
lets it follow to the redirected host, and an empty cookie file turns on
libcurl's in-memory cookie engine.
"""
function make_downloader()
    d = Downloads.Downloader()
    d.easy_hook = (easy, _) -> begin
        C = Downloads.Curl
        C.setopt(easy, C.CURLOPT_NETRC, C.CURL_NETRC_OPTIONAL)
        C.setopt(easy, C.CURLOPT_UNRESTRICTED_AUTH, Clong(1))
        C.setopt(easy, C.CURLOPT_COOKIEFILE, "")
    end
    return d
end

"""
    fetch_file(url, dest, headers, downloader; tries=3) -> Bool

`false` when the file simply is not there (404).  A 401/403 is fatal and says
what to fix, since retrying a permissions problem only wastes time.  An HTML
body means the login page came back instead of data, which is the usual symptom
of an account not yet authorised for the GES DISC archive.
"""
function fetch_file(url, dest, headers, downloader; tries::Int=3)
    tmp = dest * ".part"
    for k in 1:tries
        try
            resp = Downloads.request(url; output=tmp, headers=headers,
                                     downloader=downloader, timeout=120,
                                     progress=progress_printer(), throw=true)
            println()
            ctype = ""
            for (h, v) in resp.headers
                lowercase(h) == "content-type" && (ctype = lowercase(v))
            end
            if occursin("text/html", ctype)
                rm(tmp; force=true)
                error("got an HTML page instead of data: check the Earthdata login, " *
                      "and that the account is authorised for 'NASA GESDISC DATA ARCHIVE'")
            end
            mv(tmp, dest; force=true)
            return true
        catch e
            isfile(tmp) && rm(tmp; force=true)
            if e isa Downloads.RequestError
                st = e.response.status
                st == 404 && return false
                if st in (401, 403)
                    error("HTTP $st from GES DISC: the Earthdata account must be authorised " *
                          "for 'NASA GESDISC DATA ARCHIVE' (Earthdata profile > Applications)")
                end
                println("  HTTP $st, retry $k/$tries")
            else
                k == tries && rethrow()
                println("  $(nameof(typeof(e))): retry $k/$tries")
            end
            k < tries && sleep(5k)
        end
    end
    return false
end

function main(args=ARGS)
    if isempty(args) || "--help" in args || "-h" in args
        println("""
        usage: julia --project=easlee download/merra2.jl --start YYYY-MM-DD [--end YYYY-MM-DD] [options]

          --start D        first day (required)
          --end D          last day (default: --start)
          --collections L  comma list (default $(join(DEFAULT_COLL, ",")))
          --out DIR        target directory (default ./data/MERRA2)
          --token T        Earthdata bearer token (or set EARTHDATA_TOKEN)
          --dry-run        print the URLs, download nothing
          --no-available   do not write the AVAILABLE index
        """)
        return 0
    end

    _, o = parse_flags(args, Dict{String,Any}(
        "start" => nothing, "end" => nothing,
        "collections" => join(DEFAULT_COLL, ","), "out" => "./data/MERRA2",
        "token" => get(ENV, "EARTHDATA_TOKEN", ""),
        "dry-run" => false, "no-available" => false))

    o["start"] === nothing && error("--start YYYY-MM-DD is required")
    d0 = Date(o["start"], dateformat"yyyy-mm-dd")
    d1 = o["end"] === nothing ? d0 : Date(o["end"], dateformat"yyyy-mm-dd")
    colls = [strip(c) for c in split(o["collections"], ",") if !isempty(strip(c))]
    for c in colls
        haskey(PRODUCT, c) || error("unknown collection $c; known: $(join(sort(collect(keys(PRODUCT))), ", "))")
    end

    token   = String(o["token"])
    headers = isempty(token) ? Pair{String,String}[] :
              ["Authorization" => "Bearer $token"]
    if !o["dry-run"] && isempty(token) && !has_netrc_entry()
        error("no Earthdata credentials: give --token, set EARTHDATA_TOKEN, " *
              "or add a urs.earthdata.nasa.gov line to $(netrc_path())")
    end
    dl = o["dry-run"] ? nothing : make_downloader()
    o["dry-run"] || mkpath(o["out"])

    n_ok = n_skip = n_miss = 0
    for d in d0:Day(1):d1, c in colls
        got = false
        for st in (stream(year(d)), 401)      # 401 = the reprocessed 2020-2021 months
            fn   = "MERRA2_$st.$c.$(Dates.format(d, "yyyymmdd")).nc4"
            url  = "$BASE/$(PRODUCT[c]).5.12.4/$(Dates.format(d, "yyyy"))/$(Dates.format(d, "mm"))/$fn"
            dest = joinpath(o["out"], fn)
            if o["dry-run"]
                println(url); got = true; break
            end
            if isfile(dest) && filesize(dest) > 0
                println("exists: $dest"); n_skip += 1; got = true; break
            end
            println("get $url")
            if fetch_file(url, dest, headers, dl)
                n_ok += 1; got = true; break
            end
        end
        if !got && !o["dry-run"]
            println("  missing: $c $d"); n_miss += 1
        end
    end
    println("done: $n_ok downloaded, $n_skip already present, $n_miss missing")
    if !o["dry-run"] && !o["no-available"] && (n_ok + n_skip) > 0
        write_available(o["out"], "MERRA2")
    end
    return n_miss == 0 ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
