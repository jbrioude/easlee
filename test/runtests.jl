#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# runtests.jl — run the self-contained EASLEE test suite.
#
#     julia --project=easlee test/runtests.jl            # all tests
#     julia --project=easlee test/runtests.jl plume      # only tests matching "plume"
#
# Every test_*.jl in this directory is a standalone script using synthetic
# data only — NO meteorological archive is required, so this suite runs
# anywhere, right after `Pkg.instantiate()`.  Each test runs in its own Julia
# process; it fails if the process exits non-zero or prints "FAIL".
#
# Checks that DO need a met archive on disk live in ../validation (see
# ../validation/README.md) and are not run from here.
# ─────────────────────────────────────────────────────────────────────────────

const HERE   = @__DIR__
const REPO   = abspath(joinpath(HERE, ".."))
const FILTER = isempty(ARGS) ? "" : ARGS[1]

tests = sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"),
                    readdir(HERE)))
isempty(FILTER) || (tests = filter(t -> occursin(FILTER, t), tests))

# review 2026-09-04 A6: Base.julia_cmd() carries no -t, so every test process
# was single-threaded (test_wellmixed_pbl: 47 min).  Forward the thread count.
jl = `$(Base.julia_cmd()) -t $(Threads.nthreads()) --project=$(joinpath(REPO, "easlee"))`

passed, failed = String[], String[]
for t in tests
    print(rpad(t, 38), " ... "); flush(stdout)
    out = IOBuffer()
    t0  = time()
    ok  = try
        p = run(pipeline(`$jl $(joinpath(HERE, t))`; stdout = out, stderr = out);
                wait = true)
        success(p)
    catch
        false
    end
    txt = String(take!(out))
    ok  = ok && !occursin("FAIL", txt)
    dt  = round(time() - t0; digits = 1)
    if ok
        println("PASS  ($(dt)s)")
        push!(passed, t)
    else
        println("FAIL  ($(dt)s)")
        push!(failed, t)
        # replay the tail of the failing test's output for diagnosis
        lines = split(txt, '\n')
        for l in lines[max(1, end - 20):end]
            println("    │ ", l)
        end
    end
end

println("\n", length(passed), "/", length(tests), " tests passed")
if !isempty(failed)
    println("failed: ", join(failed, ", "))
    exit(1)
end
