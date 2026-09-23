
# ── Portable per-thread PRNG ──────────────────────────────────────────────────
# xorshift64 + Box-Muller.  Works identically on CPU and GPU (no stdlib calls).
# Each parcel carries one UInt64 RNG state in `trajsin.rng_state[i]`.
#
# Usage inside a kernel:
#   state = rng_states[i]
#   state, z = xorshift_randn_clamped(state)
#   rng_states[i] = state

@inline function xorshift64(s::UInt64)
    # NB: Julia's ^ is exponentiation, not XOR — the C idiom `s ^= s << 13`
    # silently computed integer powers here and collapsed the state to 1.
    s = xor(s, s << 13)
    s = xor(s, s >> 7)
    s = xor(s, s << 17)
    return s
end

# Returns (new_state, u) where u ∈ (0,1).
@inline function xorshift_uniform01(s::UInt64)
    s  = xorshift64(s)
    u  = Float32(s) * Float32(5.421010862427522e-20) + Float32(1e-7)  # avoid 0
    return s, min(u, 1.0f0 - Float32(1e-7))
end

# Returns (new_state, z) where z ~ N(0,1) clamped to [-3,3].
@inline function xorshift_randn_clamped(s::UInt64)
    s, u1 = xorshift_uniform01(s)
    s, u2 = xorshift_uniform01(s)
    r  = sqrt(-2.0f0 * log(u1))
    θ  = 2.0f0 * Float32(π) * u2
    return s, clamp(r * cos(θ), -3.0f0, 3.0f0)
end

# review 2026-09-04 A7: per-parcel seed = f(slot, first-advance time, run seed).
# The old seed depended on the slot index only; with the slot zeroed on reuse
# (boundary.jl:88) every later occupant replayed the same noise, and there was
# no way to change the ensemble.  splitmix64 is the standard 64-bit finaliser
# (Steele et al. 2014) — a bijection, so distinct inputs give distinct seeds.
@inline function splitmix64(z::UInt64)
    z += 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# Run seed: no seed reaches advance_dispatch! (step_model! threads its `rng`
# only into release/initialize), so it is a module-level Ref read once from
# EASLEE_SEED (any integer), defaulting to a fixed constant (reproducible runs).
const RNG_RUN_SEED = Ref{UInt64}(
    let v = get(ENV, "EASLEE_SEED", "")
        isempty(v) ? UInt64(0x45a51ee5eed00001) :
                     UInt64(reinterpret(UInt64, parse(Int64, v)))
    end)

# Seed for parcel in slot i first advanced at dtime_ms (Dates.value(Dtime),
# ms since epoch) under run seed `run_seed`.  Never returns 0 (xorshift64's
# absorbing state).
@inline function rng_seed(i::Int, dtime_ms::Int64, run_seed::UInt64)
    s = splitmix64(UInt64(i) * UInt64(6364136223846793005) ⊻ run_seed)
    s = splitmix64(s ⊻ reinterpret(UInt64, dtime_ms))
    s == UInt64(0) && (s = UInt64(0x9e3779b97f4a7c15))
    return xorshift64(xorshift64(s))   # warm up
end

# Seed for parcel i (legacy 1-arg form, still used by convection.jl's lazy
# fallback): slot + run seed, no injection time.
@inline rng_seed(i::Int) = rng_seed(i, Int64(0), RNG_RUN_SEED[])
