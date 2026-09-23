# ──────────────────────────────────────────────────────────────────────────────
# core/metadapter.jl — met time interpolation index (MetAdapter module)
#
# Maps the current model date onto the fractional index between the two loaded
# met time slots, per nest (consumed by the interpolation kernels).
# ──────────────────────────────────────────────────────────────────────────────

module MetAdapter

using Dates

export met_time_index_array

"""
    met_time_index_array(currentdate, nests) → Vector{Float32}

Compute the fractional time-interpolation index for each nest.
`nests` is the `Vector{NestContext}` from `ModelContext`; each element carries
a typed `nestrunconfig.t_index` (synced from the legacy Dict after every
`readmetfiles!` call) and a `metconfig.filedate` vector.
"""
function met_time_index_array(currentdate, nests)
    nnest = length(nests)
    t_array = zeros(Float32, nnest)
    for n in 1:nnest
        t_idx = nests[n].nestrunconfig.t_index
        # (0, 0) is the initial sentinel before the first readmetfiles! call.
        # Return 1.0 so callers get a valid (but unused) interpolation weight.
        if t_idx[1] == 0 || t_idx[2] == 0
            t_array[n] = 1.0f0
            continue
        end
        t1 = nests[n].metconfig.filedate[t_idx[1]]
        t2 = nests[n].metconfig.filedate[t_idx[2]]
        dt = t2 - t1
        # Clamp to the bracket [1,2].  For a normal advancing run the weight is
        # already in range; it exceeds it only when the run has passed the last
        # met file (readmetfiles! then holds the boundary bracket) — clamping
        # here HOLDS that last frame instead of extrapolating the winds.
        t_array[n] = clamp(Float32((currentdate - t1) / dt + 1.0), 1.0f0, 2.0f0)
    end
    return t_array
end

end
