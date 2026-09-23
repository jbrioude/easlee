# ──────────────────────────────────────────────────────────────────────────────
# core/state.jl — parcel state containers (State module)
#
# ParcelState: typed struct-of-arrays for all per-parcel variables (position in
# met-grid indices, turbulent/mesoscale velocities, mass, timing, status flags,
# output-kernel bandwidth/covariance, GPU RNG state).  ReleaseSchedule holds the
# per-box release definition; OutputState the output accumulators.
# ──────────────────────────────────────────────────────────────────────────────

module State

export ParcelState, ReleaseSchedule, OutputState, sync_dtime_ns!, isunused, isactive

using Dates
using StructArrays

"""
Typed struct-of-arrays container for Lagrangian parcels, parameterized over
floating-point precision F and integer precision I.

Default: F = Float32, I = Int32.
Switch to F = Float64 for a high-precision run:
    state = ParcelState(n; F = Float64)
"""
mutable struct ParcelState{F<:AbstractFloat, I<:Integer}
    indx::Vector{F}
    indy::Vector{F}
    indl::Vector{F}
    uturb::Vector{F}
    vturb::Vector{F}
    wturb::Vector{F}
    ugrad::Vector{F}
    vgrad::Vector{F}
    wgrad::Vector{F}
    mass::Vector{F}
    relbox::Vector{I}
    nest::Vector{I}
    id::Vector{I}
    dt::Vector{F}
    Dtimeinit::Vector{DateTime}
    Dtime::Vector{DateTime}
    status::Vector{UInt8}      # 0 unused, 1 active, 2 deposited, 3 outside
    is_new::BitVector          # true only on the time step of first activation
    miss::Vector{F}
    miss2::Vector{F}
    miss3::Vector{F}
    miss4::Vector{F}
    # Mesoscale velocity fluctuations (FLEXPART mesovel): per-parcel OU
    # velocities [m/s, m/s, level-index/s]; 0 = uninitialised (fresh draw).
    mesou::Vector{F}
    mesov::Vector{F}
    mesow::Vector{F}
    # Adaptive-output-kernel bandwidth h_i [m]: RMS distance of the k nearest
    # neighbours (statistical spacing, NOT a physical puff size); updated at
    # the output sampling cadence when output.kernel = "adaptive".
    hband::Vector{F}
    # Per-parcel horizontal covariance [m²] (output.kernel = "covariance"):
    # P = [cxx cxy; cxy cyy], evolved every sync step by the resolved wind
    # deformation P ← (I+∇u·dt)·P·(I+∇u·dt)ᵀ plus turbulent growth 2D·dt —
    # the PHYSICAL footprint each parcel represents (puff-style), used as the
    # Gaussian width at deposit/receptor evaluation.  cxx = 0 marks
    # uninitialised (set to the release footprint on first touch).
    cxx::Vector{F}
    cxy::Vector{F}
    cyy::Vector{F}
    # Air mass the parcel represents [kg] for the ice-microphysics closure
    # (mixing ratio q = mass/airm): set at release from the plume member's
    # detrainment water fraction, or lazily from [species.ice] q0.
    # 0 = unassigned.
    airm::Vector{F}
    # FLOOR on airm [kg] (0 = none): the seeded umbrella allotment
    # (mass/q_det).  update_airm! writes max(kNN estimate, airm0) — the kNN
    # of a dense fresh release cluster otherwise OVERWRITES the physical
    # seed value by orders of magnitude and manufactures a flash-freeze
    # (41–53 Tg in the first hours of every HT run; diagnosed 2026-08-07).
    # A parcel's represented air can grow by mixing; it cannot shrink below
    # its initial umbrella allotment.
    airm0::Vector{F}
    nest_buffers::Vector{Vector{Int32}}
    # GPU / KernelAbstractions support
    rng_state::Vector{UInt64}  # per-parcel xorshift64 state (0 = uninitialised)
    Dtime_ns::Vector{Int64}    # Dtime as Int64 milliseconds since epoch (GPU-safe)
end

function ParcelState(n::Integer; F::Type{<:AbstractFloat}=Float32, I::Type{<:Integer}=Int32)
    nI = Int(n)
    zero_time = DateTime(0)
    z = zero(F)
    return ParcelState{F, I}(
        fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI),
        fill(I(0), nI), fill(I(1), nI), I.(collect(1:nI)),
        fill(z, nI),
        fill(zero_time, nI), fill(zero_time, nI),
        fill(UInt8(0), nI), falses(nI),
        fill(z, nI), fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI),
        fill(z, nI), fill(z, nI), fill(z, nI),
        fill(z, nI),
        fill(z, nI),
        Vector{Vector{Int32}}(),
        fill(UInt64(0), nI),
        fill(Int64(0), nI),
    )
end

# Keep Dtime_ns in sync with Dtime after each modification of Dtime
function sync_dtime_ns!(state::ParcelState)
    @inbounds for i in eachindex(state.Dtime)
        state.Dtime_ns[i] = Int64(Dates.value(state.Dtime[i]))
    end
end

@inline isunused(state::ParcelState, i::Int) = state.status[i] == 0x00
@inline isactive(state::ParcelState, i::Int) = state.status[i] == 0x01

"""
Release schedule parameterized over float/int precision, matching ParcelState.
The x/y release bounds are stored as geographic lon/lat here; conversion to
grid indices happens in Operators.release_parcels!.
"""
struct ReleaseSchedule{F<:AbstractFloat, I<:Integer}
    xrelease1::Vector{F}
    xrelease2::Vector{F}
    yrelease1::Vector{F}
    yrelease2::Vector{F}
    zrelease1::Vector{F}
    zrelease2::Vector{F}
    timerelease1::Vector{DateTime}
    timerelease2::Vector{DateTime}
    nbpart::Vector{I}
    mass_per_box::Vector{F}
    zkind::Vector{String}
    boxid::Vector{I}
    species::Vector{I}      # FLEXPART species number per release box (SPECNUM_REL)
end

mutable struct OutputState
    concarray
    massarray
    outinterfaces::Vector{Float32}
end

end
