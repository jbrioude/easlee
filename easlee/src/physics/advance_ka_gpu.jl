# physics/advance_ka_gpu.jl
# CUDA side of the KA advance (EASLEE_ADVANCE=ka_gpu): resident device buffers
# + upload/launch/download for advance_parcel_kernel! (advance_ka.jl).
# Included by core/operators.jl ONLY when the GPU stack is loaded
# (EASLEE_ADVANCE=ka_gpu / EASLEE_GPU=1 — see EASLEE.jl); the CPU ka_cpu path
# never references these symbols.  Code moved verbatim from advance_ka.jl.

# ── GPU launch: resident device buffers, upload → kernel → download ───────────
# The kernel writes the 15 parcel-state arrays (indx..mesow) and reads the 11
# met arrays.  We keep one set of device buffers alive across steps (reallocated
# only when the parcel count or met dims change) and copyto! the STATE into them
# each step, because the CPU side mutates state (release/deposit/reassign)
# between advances.
#
# review 2026-09-04 A5: the 8 met arrays were re-uploaded every step although
# they only change when readmetfiles! refills a slot.  Change detection (this
# file cannot see initialize.jl/prefetch.jl): (a) nrc.t_index — the met bracket
# written back by readmetfiles!/prefetch at the exact moment a slot is refilled
# (initialize.jl:980, prefetch.jl:192) → bracket change = full re-upload; (b) a
# sampled fingerprint of each host array (~1024 strided elements, both slots,
# bit-hashed so NaN is stable) catches in-place edits that leave t_index alone
# (hanna_fill_tke!, isentropic θ̇ refresh, ...) → only the changed arrays are
# re-uploaded.  A full checksum would cost as much as the upload it avoids.
# One buffer set PER NEST (Dict keyed by nest index): the single global slot
# thrashed (P/metdims differ between nests) and reallocated every step.
mutable struct _KAGpuBuffers
    P       ::Int
    metdims ::NTuple{4,Int}
    s       ::Vector{Any}   # 15 device state arrays (12 + the mesoscale OU triple)
    m       ::Vector{Any}   # 11 device met arrays (8 3-D + 3 surface, the last three empty unless Hanna)
    idx     ::Any           # device Int index vector (resized as needed)
    t_index ::NTuple{2,Int}      # met bracket the resident met was uploaded from
    mfp     ::NTuple{11,UInt64}  # sampled host fingerprints of the met arrays
end
const _KAGPU = Dict{Int,_KAGpuBuffers}()

# Sampled fingerprint of a host met array: ~1024 strided elements hashed by
# their bit pattern (splitmix64, rng.jl).  Cost: microseconds; a changed slot
# alters ~half the samples.
function _met_fingerprint(a::AbstractArray{Float32})
    n  = length(a)
    st = max(1, n ÷ 1024)
    h  = UInt64(n)
    @inbounds for k in 1:st:n
        h = splitmix64(h ⊻ UInt64(reinterpret(UInt32, a[k])))
    end
    return h
end

function _ka_gpu_advance!(kernel!, backend, state, idx_arr,
                          uarray, varray, levdotarray,
                          zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
                          ustar_arr, pblh_arr, ol_arr, scal,
                          nest::Int = 1, t_index::NTuple{2,Int} = (0, 0))
    P  = length(state.indx)
    md = size(uarray)
    # 2026-09-21: the mesoscale OU state joins the resident state set (13-15);
    # the kernel writes it, so it comes back with the rest.
    sc = (state.indx, state.indy, state.indl, state.uturb, state.vturb,
          state.wturb, state.dt, state.miss2, state.miss3, state.miss4,
          state.Dtime_ns, state.rng_state,
          state.mesou, state.mesov, state.mesow)
    # 2026-09-21: the three surface arrays join the resident met set (positions
    # 9-11).  They are length-0 unless on-the-fly Hanna is active, so the
    # fingerprint and the upload are free in every other configuration.
    mc = (uarray, varray, levdotarray, zlev_arr, theta_arr, tke_arr, rho_arr, drhodz_arr,
          ustar_arr, pblh_arr, ol_arr)
    mfp = ntuple(k -> _met_fingerprint(mc[k]), 11)

    buf = get(_KAGPU, nest, nothing)
    if buf === nothing || buf.P != P || buf.metdims != md
        buf = _KAGpuBuffers(P, md,
                            Any[CuArray(a) for a in sc],
                            Any[CuArray(a) for a in mc],
                            CuArray(idx_arr), t_index, mfp)
        _KAGPU[nest] = buf
    else
        @inbounds for k in 1:15; copyto!(buf.s[k], sc[k]); end
        if buf.t_index != t_index
            @inbounds for k in 1:11; copyto!(buf.m[k], mc[k]); end
        else
            @inbounds for k in 1:11
                buf.mfp[k] == mfp[k] || copyto!(buf.m[k], mc[k])
            end
        end
        buf.t_index = t_index
        buf.mfp     = mfp
        if length(buf.idx) != length(idx_arr)
            buf.idx = CuArray(idx_arr)
        else
            copyto!(buf.idx, idx_arr)
        end
    end

    kernel!(buf.idx, buf.s[1], buf.s[2], buf.s[3], buf.s[4], buf.s[5], buf.s[6],
            buf.s[7], buf.s[8], buf.s[9], buf.s[10], buf.s[11], buf.s[12],
            buf.m[1], buf.m[2], buf.m[3], buf.m[4], buf.m[5], buf.m[6], buf.m[7], buf.m[8],
            buf.m[9], buf.m[10], buf.m[11],
            buf.s[13], buf.s[14], buf.s[15],
            scal...; ndrange = length(idx_arr))
    KernelAbstractions.synchronize(backend)

    # Download the 12 kernel-written state arrays back to the host.
    @inbounds for (k, dst) in enumerate(sc); copyto!(dst, buf.s[k]); end
    return nothing
end
