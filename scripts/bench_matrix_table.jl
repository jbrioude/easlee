#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# bench_matrix_table.jl — paper tables from the benchmark matrix
# (scripts/bench_matrix.sh).
#
#     julia --project=easlee scripts/bench_matrix_table.jl <matrix output dir>
#
# Reads the per-run logs in that directory and prints two markdown tables:
# (A) every row with its cost-centre split, (B) per-process time and speedup at
# 1 / 16 / 64 threads and on the GPU, plus a log-log scaling exponent per
# process fitted on the no-prefetch CPU rows.
#
# Groups (timeloop.jl): metread = GRIB read + field processing; advection =
# advance_dispatch! (advection + turbulence, on the GPU including the
# host<->device copies); output = deposit_to_grid! (the output KERNEL only);
# release = release + initialize_new_parcels!.  "untimed" = TOTAL(run) minus the
# sum of the groups: the NetCDF frame dump, nest reassignment, sync and loop
# overhead.  Convection and deposition are OFF in the ETEX protocol and do not
# appear.
# ─────────────────────────────────────────────────────────────────────────────

using Printf

const GROUPS = ["metread", "advection", "output", "release"]

const ORDER = ["cpu_t1", "cpu_t2", "cpu_t4", "cpu_t8", "cpu_t16", "cpu_t32", "cpu_t64",
               "cpu_t16_pf", "cpu_t32_pf", "cpu_t64_pf",
               "gpu_t4", "gpu_t16i", "gpu_t4_pf", "gpu_t16i_pf"]

const LABELS = Dict(
    "metread"   => "met fields read + process",
    "advection" => "advection + turbulence",
    "output"    => "output kernel (deposit_to_grid!)",
    "release"   => "release",
    "untimed"   => "untimed (frame dump, reassign, sync)",
    "total"     => "TOTAL(run)")

"Parse one run log into its mode, total and per-group seconds."
function read_log(path::AbstractString)
    s = read(path, String)
    m   = match(r"mode = :([a-z_]+)", s)
    tot = match(r"(?m)^\s*TOTAL\(run\)\s+([\d.]+)", s)
    (m === nothing || tot === nothing) && return nothing

    r = Dict{String,Any}("mode" => m.captures[1], "total" => parse(Float64, tot.captures[1]))
    for g in GROUPS
        mm = match(Regex("(?m)^\\s*$g\\s+([\\d.]+)"), s)
        r[g] = mm === nothing ? NaN : parse(Float64, mm.captures[1])
    end
    r["untimed"] = r["total"] - sum(r[g] for g in GROUPS)

    tag = basename(path)[1:end-4]
    mt = match(r"^(cpu|gpu)_t(\d+)(i?)(_pf)?$", tag)
    mt === nothing && return nothing
    r["thr"] = parse(Int, mt.captures[2])
    r["pf"]  = mt.captures[4] !== nothing
    r["gpu"] = mt.captures[1] == "gpu"
    return tag, r
end

"Speed-up of `b` relative to `a`, or an en dash when it cannot be formed."
speedup(a, b) = (isfinite(a) && isfinite(b) && a > 0) ? @sprintf("%.2f×", b/a) : "–"

"Least-squares slope of log(time) against log(threads): t ∝ N^b."
function exponent(cpu, g)
    pts = [(log(n), log(r[g])) for (n, r) in cpu if isfinite(r[g]) && r[g] > 0]
    length(pts) < 3 && return "–"
    mx = sum(first, pts) / length(pts)
    my = sum(last, pts) / length(pts)
    b  = sum((x-mx)*(y-my) for (x, y) in pts) / sum((x-mx)^2 for (x, _) in pts)
    ns = [n for (n, _) in cpu]
    return @sprintf("t ∝ N^%+.2f (%d pts, %d–%d t)", b, length(pts), minimum(ns), maximum(ns))
end

function main(args=ARGS)
    if isempty(args) || "--help" in args
        println("usage: julia --project=easlee scripts/bench_matrix_table.jl <matrix output dir>")
        return isempty(args) ? 1 : 0
    end
    dir = args[1]
    isdir(dir) || error("not a directory: $dir")

    rows = Dict{String,Dict{String,Any}}()
    for f in sort(filter(f -> endswith(f, ".log"), readdir(dir; join=true)))
        got = read_log(f)
        got === nothing || (rows[got[1]] = got[2])
    end
    isempty(rows) && error("no usable .log files in $dir")

    have = filter(t -> haskey(rows, t), ORDER)
    ref1  = get(rows, "cpu_t1", nothing)
    ref16 = get(rows, "cpu_t16", nothing)

    println("### Table A — every row (run-only seconds; ETEX 1 M parcels, 2 days; `$(basename(dir))`)\n")
    println("| row | backend | threads | prefetch | TOTAL | metread | advection | output kernel | release | untimed | vs 1 t | vs 16 t |")
    println("|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for t in have
        r = rows[t]
        @printf("| %s | %s | %d | %s | %.0f | %.1f | %.1f | %.1f | %.1f | %.1f | %s | %s |\n",
                t, r["mode"], r["thr"], r["pf"] ? "on" : "off", r["total"],
                r["metread"], r["advection"], r["output"], r["release"], r["untimed"],
                ref1  === nothing ? "–" : speedup(r["total"], ref1["total"]),
                ref16 === nothing ? "–" : speedup(r["total"], ref16["total"]))
    end

    println("\n### Table B — per process: time [s] and speedup (vs 1 thread; in brackets vs 16 threads)\n")
    cols = [(t, l) for (t, l) in (("cpu_t1", "1 t"), ("cpu_t16", "16 t"), ("cpu_t64", "64 t"),
                                  ("gpu_t4", "GPU (4 t)"), ("gpu_t4_pf", "GPU + prefetch"))
            if haskey(rows, t)]
    println("| process | " * join((l for (_, l) in cols), " | ") * " | scaling exponent (CPU, no prefetch) |")
    println("|---|" * repeat("---:|", length(cols) + 1))

    cpu = [(rows[t]["thr"], rows[t]) for t in have if !rows[t]["gpu"] && !rows[t]["pf"]]
    for g in ("metread", "advection", "output", "release", "untimed", "total")
        cells = String[]
        for (t, _) in cols
            v = rows[t][g]
            s = @sprintf("%.1f", v)
            ref1 === nothing || t == "cpu_t1" || (s *= " (" * speedup(v, ref1[g]))
            ref16 === nothing || t in ("cpu_t1", "cpu_t16") || (s *= "; [" * speedup(v, ref16[g]) * "]")
            ref1 === nothing || t == "cpu_t1" || (s *= ")")
            push!(cells, s)
        end
        println("| $(LABELS[g]) | " * join(cells, " | ") * " | " * exponent(cpu, g) * " |")
    end
    println("| convection | off in the ETEX protocol — see the full-physics table (2026-07-09): 4.7× at 6 threads |" * repeat(" |", length(cols)))
    println("| dry / wet deposition | off in the ETEX protocol — full-physics table: wet dep 1.5× at 6 threads |" * repeat(" |", length(cols)))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
