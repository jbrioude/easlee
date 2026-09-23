
# ── Per-group run benchmarking ────────────────────────────────────────────────
# Accumulates wall time, allocations, and call counts for named phases
# (met read, advection, convection, mixing, chemistry, deposition, output …)
# across the whole run, then prints an @time-style summary table.
#
# Enabled via the config ([run] benchmark = true): build_run creates a
# BenchTimer and threads it into the time loop and step_model!.  When the timer
# is `nothing` the wrappers call straight through — zero overhead.

module Bench

using Printf
export BenchTimer, timeit!, bench_report

mutable struct BenchTimer
    time_ns::Dict{Symbol,UInt64}    # accumulated wall time per group
    bytes::Dict{Symbol,Int64}       # accumulated allocations per group
    calls::Dict{Symbol,Int}         # call count per group
    order::Vector{Symbol}           # first-seen order for the report
    t_start::UInt64
end
BenchTimer() = BenchTimer(Dict{Symbol,UInt64}(), Dict{Symbol,Int64}(),
                          Dict{Symbol,Int}(), Symbol[], time_ns())

@inline function _record!(b::BenchTimer, name::Symbol, dt::UInt64, db::Int64)
    if haskey(b.time_ns, name)
        b.time_ns[name] += dt; b.bytes[name] += db; b.calls[name] += 1
    else
        b.time_ns[name] = dt; b.bytes[name] = db; b.calls[name] = 1
        push!(b.order, name)
    end
    return nothing
end

# Time a phase.  `timer === nothing` → call straight through (no overhead).
@inline function timeit!(f, timer, name::Symbol)
    timer === nothing && return f()
    b0 = Base.gc_bytes()
    t0 = time_ns()
    v  = f()
    _record!(timer::BenchTimer, name, time_ns() - t0, Base.gc_bytes() - b0)
    return v
end

_humanbytes(n) = n < 1024 ? "$(n) B" :
                 n < 1024^2 ? string(round(n/1024, digits=1), " KiB") :
                 n < 1024^3 ? string(round(n/1024^2, digits=1), " MiB") :
                              string(round(n/1024^3, digits=2), " GiB")

# Pretty per-group summary, sorted by total time descending.
function bench_report(b::BenchTimer; io::IO = stdout)
    total = time_ns() - b.t_start
    println(io, "\n── benchmark (per-group totals) ──────────────────────────────────")
    @printf(io, "  %-14s %10s %8s %12s %9s\n", "group", "time[s]", "%", "alloc", "calls")
    for s in sort(collect(b.order); by = s -> -Int(b.time_ns[s]))
        t = b.time_ns[s] / 1e9
        @printf(io, "  %-14s %10.3f %7.1f%% %12s %9d\n",
            String(s), t, 100t / (total / 1e9), _humanbytes(b.bytes[s]), b.calls[s])
    end
    @printf(io, "  %-14s %10.3f\n", "TOTAL(run)", total / 1e9)
    return nothing
end

end # module Bench
