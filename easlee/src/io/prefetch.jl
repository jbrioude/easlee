# io/prefetch.jl
# Asynchronous met-file prefetching — third time slice loaded in the background.
#
# The model holds 2 met time slices in memory (slot dimension of the metfields
# arrays).  Legacy behaviour: when currentdate passes the newest loaded date,
# the model BLOCKS while the next file is read and decoded (GRIB/NetCDF) and
# the derived fields are computed — typically seconds per rotation.
#
# This module adds a third, staging copy of the met arrays and a background
# task (Threads.@spawn) that loads the *next* required file into staging while
# the model integrates.  At rotation time the main thread:
#   1. waits for the task (normally already finished → no stall),
#   2. memcpys staging → the stale slot (~10–100 ms, vs seconds for decode),
#   3. updates t_index with the SAME bracket rules as legacy readmetfiles!,
#   4. immediately spawns the prefetch of the following file.
#
# Forward (direction=+1) and backward (direction=−1) runs are both handled:
# the next file is filename[max(t_idx)+1] forward, filename[min(t_idx)−1]
# backward — the same rule the legacy rotation uses.
#
# Robustness:
#   • If the prefetched file is NOT the one needed (direction reversal,
#     restart, big time jump), it falls back to the legacy synchronous load.
#   • The fallback first waits for any in-flight task, so the GRIB library is
#     never called from two threads at once.
#   • With Threads.nthreads() == 1 prefetching is disabled and every load is
#     synchronous (identical to legacy).
#
# Cost: one extra in-memory copy of all met fields per nest (the staging
# arrays mirror metfields3D/2D shapes).
#
# Requires: include after initialize.jl (uses loadgribmetfiles!,
# loadnetcdfmetfiles!, metgribfieldcalculation!, metaromefieldcalculation!,
# checkmetfields!, readmetfiles!).

mutable struct MetPrefetcher
    staging3D :: Dict{Any,Any}
    staging2D :: Dict{Any,Any}
    task      :: Union{Task, Nothing}
    file_idx  :: Int     # file index being (or already) prefetched; 0 = none
    enabled   :: Bool
end

"""
    MetPrefetcher(metfields3D, metfields2D) → MetPrefetcher

Allocate staging buffers mirroring the shapes of the live met dicts.
Prefetching is enabled only when Julia runs with ≥ 2 threads.
"""
function MetPrefetcher(metfields3D::Dict, metfields2D::Dict)
    mk(d) = Dict{Any,Any}(k => (v isa AbstractArray ? zeros(eltype(v), size(v)) : deepcopy(v))
                          for (k, v) in d)
    MetPrefetcher(mk(metfields3D), mk(metfields2D), nothing, 0,
                  Threads.nthreads() > 1)
end

# ── Background load: file → staging slot 1 + per-slot field calculations ─────

function _prefetch_load!(pf::MetPrefetcher, nestctx, metconfig, fileidx::Int)
    mc      = nestctx.metconfig
    fname   = mc.filename[fileidx]
    globalx = nestctx.gridconfig.globalx
    _reset_derived_surface!(pf.staging2D, 1)
    _load_met_file!(fname, 1, mc.filedate[fileidx], mc.metfileformat, mc.model,
                    globalx, metconfig, pf.staging3D, pf.staging2D)
    return nothing
end

function _start_prefetch!(pf::MetPrefetcher, nestctx, metconfig, fileidx::Int)
    nfiles = length(nestctx.metconfig.filename)
    if !pf.enabled || fileidx < 1 || fileidx > nfiles
        pf.task = nothing
        pf.file_idx = 0
        return nothing
    end
    pf.file_idx = fileidx
    # Run the (heavy, ~minutes) met load on the :interactive threadpool so it does
    # NOT occupy one of the :default compute threads.  Otherwise advance_active_parcels!
    # (@threads :static over the :default pool) waits for the prefetch-busy thread's
    # static partition, stalling advection for the whole file-load time each rotation.
    # Requires an interactive thread: launch julia with -t N,1 (e.g. -t 6,1).
    pf.task = if isdefined(Threads, :maxthreadid) && Threads.nthreads(:interactive) >= 1
        Threads.@spawn :interactive _prefetch_load!(pf, nestctx, metconfig, fileidx)
    else
        Threads.@spawn _prefetch_load!(pf, nestctx, metconfig, fileidx)
    end
    return nothing
end

# Wait for any in-flight prefetch (used before synchronous fallback loads so
# the GRIB/NetCDF C libraries are never entered from two threads at once).
function _drain_prefetch!(pf::MetPrefetcher)
    if pf.task !== nothing
        try
            wait(pf.task)
        catch err
            @warn "met prefetch task failed; falling back to synchronous reads" err
        end
        pf.task = nothing
        pf.file_idx = 0
    end
    return nothing
end

# ── Commit: staging slot 1 → live slot indt ──────────────────────────────────

function _commit_prefetch!(pf::MetPrefetcher, metfields3D, metfields2D, indt::Int)
    wait(pf.task)   # normally already done — zero stall
    for (k, src) in pf.staging3D
        (src isa AbstractArray && haskey(metfields3D, k)) || continue
        dst = metfields3D[k]
        copyto!(selectdim(dst, ndims(dst), indt), selectdim(src, ndims(src), 1))
    end
    for (k, src) in pf.staging2D
        (src isa AbstractArray && haskey(metfields2D, k)) || continue
        dst = metfields2D[k]
        copyto!(selectdim(dst, ndims(dst), indt), selectdim(src, ndims(src), 1))
    end
    pf.task = nothing
    pf.file_idx = 0
    return nothing
end

# ── Main entry: drop-in replacement for readmetfiles! ─────────────────────────

"""
    readmetfiles_prefetch!(pf, nestctx, currentdate, direction, metconfig,
                           metfields3D, metfields2D, metfields3Dgpu)

Same contract and same t_index bracket rules as `readmetfiles!`, but rotations
consume the asynchronously prefetched file when available.  Handles forward
(direction = 1) and backward (direction = −1) runs.
"""
function readmetfiles_prefetch!(pf::MetPrefetcher, nestctx, currentdate, direction,
                                 metconfig, metfields3D, metfields2D, metfields3Dgpu)
    mc       = nestctx.metconfig
    filedate = mc.filedate
    t_idx    = Int[nestctx.nestrunconfig.t_index[1], nestctx.nestrunconfig.t_index[2]]

    # ── First call: synchronous double load (legacy), then start prefetching ──
    if t_idx[1] == 0 && t_idx[2] == 0
        readmetfiles!(nestctx, currentdate, direction, metconfig,
                      metfields3D, metfields2D, metfields3Dgpu)
        t1, t2 = nestctx.nestrunconfig.t_index
        nxt = direction == 1 ? max(t1, t2) + 1 : min(t1, t2) - 1
        _start_prefetch!(pf, nestctx, metconfig, nxt)
        return nothing
    end

    # ── Rotation needed?  (same test as legacy) ───────────────────────────────
    if direction == 1
        t_newest = maximum(t_idx)
        testdate = currentdate > filedate[t_newest]
    else
        t_newest = minimum(t_idx)
        testdate = currentdate < filedate[t_newest]
    end
    testdate || return nothing

    # review 2026-09-04 A3 (item 2): end-of-archive guard, same as legacy
    # readmetfiles! — forward needs max+1, backward min−1 to exist; otherwise
    # hold the boundary bracket instead of indexing filename[nfiles+1].
    nfiles  = length(filedate)
    next_ok = direction == 1 ? maximum(t_idx) + 1 <= nfiles :
                               minimum(t_idx) - 1 >= 1
    if !next_ok
        println("nothing to read (end of met archive — holding the last frame)")
        return nothing
    end

    # ── Which slot becomes stale, which file is needed (same rules as legacy) ─
    local indt::Int
    if t_idx[1] < t_idx[2]
        if direction == 1
            t_idx[1] = t_idx[2] + 1;  indt = 1
        else
            t_idx[2] = t_idx[1] - 1;  indt = 2
        end
    else
        if direction == 1
            t_idx[2] = t_idx[1] + 1;  indt = 2
        else
            t_idx[1] = t_idx[2] - 1;  indt = 1
        end
    end
    needed = t_idx[indt]

    if pf.task !== nothing && pf.file_idx == needed
        # Hot path: the needed file is already (being) loaded in staging.
        _commit_prefetch!(pf, metfields3D, metfields2D, indt)
        println("index=", indt, " prefetched file committed: ", mc.filename[needed])
    else
        # Fallback: synchronous load into the stale slot (legacy code path).
        _drain_prefetch!(pf)
        fname   = mc.filename[needed]
        globalx = nestctx.gridconfig.globalx
        _reset_derived_surface!(metfields2D, indt)
        println("index=", indt, " file to read (sync)=", fname)
        _load_met_file!(fname, indt, filedate[needed], mc.metfileformat, mc.model,
                        globalx, metconfig, metfields3D, metfields2D)
    end

    nestctx.nestrunconfig.t_index = (t_idx[1], t_idx[2])

    # ── Queue the next file ────────────────────────────────────────────────────
    nxt = direction == 1 ? needed + 1 : needed - 1
    _start_prefetch!(pf, nestctx, metconfig, nxt)
    return nothing
end
