# ─────────────────────────────────────────────────────────────────────────────
# cli.jl — the small argument parser the download scripts share.
#
# Deliberately not an argument-parsing package: these scripts run from a plain
# `julia --project=easlee download/<script>.jl`, and EASLEE's dependency list
# should not grow for the sake of `--flag value`.
#
#     positional, opts = parse_flags(ARGS, Dict{String,Any}(
#         "out" => nothing,        # takes a value, default nothing
#         "jobs" => "4",           # takes a value, default "4"
#         "dry-run" => false))     # a Bool default makes it a switch
#
# Accepts `--flag value` and `--flag=value`; a flag whose default is `Bool`
# takes no value and becomes `true` when present.  Anything not starting with
# `--` is returned, in order, as a positional argument.
# ─────────────────────────────────────────────────────────────────────────────

function parse_flags(args, spec::Dict{String,Any})
    opts = copy(spec)
    positional = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if !startswith(a, "--")
            push!(positional, a)
            i += 1
            continue
        end
        key = a[3:end]
        if occursin('=', key)
            key, val = split(key, '=', limit=2)
            haskey(opts, key) || error("unknown option --$key")
            opts[key] = val
        elseif haskey(opts, key) && spec[key] isa Bool
            opts[key] = true
        else
            haskey(opts, key) || error("unknown option --$key")
            i += 1
            i <= length(args) || error("--$key needs a value")
            opts[key] = args[i]
        end
        i += 1
    end
    return positional, opts
end

"""
    progress_printer() -> f(total, now)

A `Downloads` progress callback that redraws at most twice a second, so a long
transfer leaves a handful of lines in a log instead of thousands.
"""
function progress_printer()
    last = Ref(0.0)
    return function (total, now)
        total <= 0 && return
        t = time()
        (t - last[] < 0.5 && now < total) && return
        last[] = t
        print("\r  ", lpad(round(now/1e6, digits=1), 8), " / ",
              round(total/1e6, digits=1), " MB")
        flush(stdout)
    end
end

"""
    set_stall_abort!(easy)

Abort a transfer only when it STALLS (below 1 kB/s for 120 s), never because it
is long.  `Downloads`' `timeout` keyword is a limit on the WHOLE request: the
120 s the scripts used to pass killed every MERRA-2 file (2.2 GB) at ~1.8 GB,
three times, after which the file was reported as absent.  Call it from a
`Downloader`'s `easy_hook`.
"""
function set_stall_abort!(easy)
    C = Downloads.Curl
    C.setopt(easy, C.CURLOPT_LOW_SPEED_LIMIT, Clong(1024))
    C.setopt(easy, C.CURLOPT_LOW_SPEED_TIME, Clong(120))
    return nothing
end

"A `Downloader` with no total time limit that still aborts a stalled transfer."
function stall_safe_downloader()
    d = Downloads.Downloader()
    d.easy_hook = (easy, _) -> set_stall_abort!(easy)
    return d
end
