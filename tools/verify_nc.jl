#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# verify_nc.jl — verify and sanitise EASLEE NetCDF output (any layout: the
# engine's grid_<name>.nc, or the older single-species files).
#
# For each file it reports, per DATA variable, how many values are finite, NaN
# and Inf, and the range of the finite data.  With --fix it rewrites the file in
# place with every non-finite value replaced by NaN, so a file only ever holds
# finite numbers or the NaN fill, never Inf.
#
#     julia --project=easlee tools/verify_nc.jl juliaout/<run>/grid_<name>.nc [more.nc] [--fix]
# ─────────────────────────────────────────────────────────────────────────────

using NetCDF
using Printf

const COORDS = Set(["lon", "lat", "longitude", "latitude", "plev", "height", "time", "numpoint"])

function main(args=ARGS)
    paths = filter(a -> !startswith(a, "--"), args)
    fix   = "--fix" in args
    if isempty(paths) || "--help" in args
        println("usage: julia --project=easlee tools/verify_nc.jl <file.nc> [more.nc ...] [--fix]")
        return isempty(paths) ? 1 : 0
    end

    rc = 0
    for path in paths
        println("\n=== $path ===")
        isfile(path) || (println("  missing"); rc = 1; continue)
        bad = false
        info = ncinfo(path)
        for (name, _) in info.vars
            name in COORDS && continue
            a = Array{Float64}(ncread(path, name))
            ninf = count(isinf, a)
            nnan = count(isnan, a)
            nfin = count(isfinite, a)
            fin  = filter(isfinite, a)
            lo, hi = isempty(fin) ? ("nothing", "nothing") : (string(minimum(fin)), string(maximum(fin)))
            flag = ninf > 0 ? "  <-- HAS INF!" : ""
            @printf("  %-12s finite=%9d NaN=%9d Inf=%6d  range=[%s,%s]%s\n",
                    name, nfin, nnan, ninf, lo, hi, flag)
            ninf > 0 && (bad = true)
            if fix && ninf > 0
                replace!(x -> isfinite(x) ? x : NaN, a)
                ncwrite(a, path, name)
            end
        end
        if fix
            println(bad ? "  -> sanitized (non-finite → NaN), written in place" : "  -> already clean")
        else
            println(bad ? "  (Inf present — rerun with --fix)" : "  (clean)")
            bad && (rc = 2)
        end
    end
    return rc
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
