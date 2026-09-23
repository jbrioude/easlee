
# ── Wet deposition (in-cloud + below-cloud scavenging) ────────────────────────
#
# Driven by the 3-D ERA5/IFS fields from the CLD companion files:
#   tpf  total precipitation flux through each model level [kg m⁻² s⁻¹ = mm/s]
#   cc   cloud fraction (0–1)
#
# Because tpf is 3-D, the vertical structure of scavenging comes for free:
# a parcel is only scavenged where precipitation actually falls through its
# level — no cloud-base/cloud-top bookkeeping from 2-D precip is needed.
#
# Scavenging coefficient (Hertel et al. 1995, as in FLEXPART):
#     Λ = A · I^B        [s⁻¹],  I = precip rate at the parcel level [mm h⁻¹]
# with separate (A, B) for in-cloud (parcel inside cloud, cc > cc_threshold)
# and below-cloud (precipitation present but no cloud at parcel level).
#
# Sub-grid correction: the archived I is a grid-cell mean.  Precipitation
# falls in a fraction fp of the cell (estimated from the column-maximum cloud
# fraction); the local rate there is I/fp and only the fraction fp of parcels
# is exposed, giving the effective grid-mean coefficient
#     Λ_eff = fp · A · (I/fp)^B = fp^(1−B) · A · I^B
#
# Mass update per step:  m ← m · exp(−Λ_eff Δt).  Removed mass is accumulated
# on a per-nest 2-D deposition grid (met-grid resolution, kg per cell).
#
# Wet deposition is skipped in backward runs (rc.direction = −1): scavenging
# is not self-adjoint and backward-mode source-receptor runs handle it as a
# sensitivity, not by removing mass.

"""
Tunable scavenging parameters.  Defaults are the standard Hertel et al. (1995)
rain coefficients used by FLEXPART for a soluble aerosol-like tracer; set
`A_incloud`/`A_below` to 0 to switch one pathway off, or scale both for less
soluble gases.
"""
Base.@kwdef struct WetDepParams
    A_below::Float32      = 1.0f-4    # below-cloud  Λ = A·I^B  [s⁻¹ (mm/h)⁻ᴮ]
    B_below::Float32      = 0.62f0
    A_incloud::Float32    = 3.36f-4   # in-cloud (rainout)
    B_incloud::Float32    = 0.79f0
    cc_threshold::Float32 = 0.05f0    # cloud fraction above which parcel is "in cloud"
    I_min::Float32        = 0.01f0    # [mm/h] no scavenging below this rate
    # review 2026-09-04 A6: per-release-box mask (empty = every box).  Lets the bulk gas
    # scheme coexist with AerosolModel (which scavenges its own boxes) in mixed runs.
    boxon::Vector{Bool}   = Bool[]
end
@inline _wetdep_box_on(p::WetDepParams, b::Integer) =
    isempty(p.boxon) || (1 <= b <= length(p.boxon) && @inbounds p.boxon[b])

"""
Per-nest accumulation of deposited mass [same unit as parcel mass, per met
grid cell].  Query/plot `wd.depgrid[n]` after the run.
"""
struct WetDepState
    depgrid::Vector{Matrix{Float64}}
end

function WetDepState(ctx)
    WetDepState([zeros(Float64, Int(n.gridconfig.Ni), Int(n.gridconfig.Nj))
                 for n in ctx.nests])
end

"""
    wet_deposition!(state, ctx, convfields, dt_s, t_array;
                    params=WetDepParams(), wdstate=nothing)

Apply one wet-deposition step to all active parcels.  `convfields` is the same
per-nest Vector of `ConvFields`/`nothing` used by `convective_redistribution!`
(tpf and cc live there).  Returns the total mass removed this step.
"""
function wet_deposition!(state, ctx, convfields, dt_s::Real, t_array;
                         params::WetDepParams = WetDepParams(),
                         wdstate::Union{WetDepState,Nothing} = nothing,
                         sphys = nothing)
    Int(ctx.runconfig.direction) == -1 && return 0.0
    total_removed = Threads.Atomic{Float64}(0.0)

    for n in 1:length(ctx.nests)
        cf = convfields[n]
        cf === nothing && continue
        nestctx = ctx.nests[n]

        n8 = eltype(state.nest)(n)
        idx = findall(i -> state.status[i] == 0x01 && state.nest[i] == n8,
                      eachindex(state.status))
        isempty(idx) && continue

        dep = wdstate === nothing ? nothing : wdstate.depgrid[n]
        _wetdep_nest!(state, nestctx, cf, Float32(dt_s), t_array[n], idx,
                      params, dep, total_removed, sphys)
    end
    return total_removed[]
end

# `sphys` (a SpeciesPhysics, or nothing) overrides the BELOW-CLOUD (A, B) with
# the release-box species' own coefficients where the species DB defines them;
# the in-cloud (nucleation) pathway keeps the bulk coefficients until per-species
# CCN/IN scavenging is wired.
function _wetdep_nest!(state, nestctx, cf, dt::Float32, indt, idx,
                       p::WetDepParams, dep, total_removed, sphys = nothing)
    nab = sphys === nothing ? 0 : length(sphys.weta)
    Ni = Int(nestctx.gridconfig.Ni)
    Nj = Int(nestctx.gridconfig.Nj)
    nz = size(cf.cc, 3)

    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    ccbufs = [Vector{Float32}(undef, nz) for _ in 1:nthr]
    removed = zeros(Float64, nthr)

    Threads.@threads :static for i in idx
        tid = Threads.threadid()

        indx = state.indx[i]
        indy = state.indy[i]
        (indx < 1.0f0 || indx >= Float32(Ni) ||
         indy < 1.0f0 || indy >= Float32(Nj)) && continue
        _wetdep_box_on(p, Int(state.relbox[i])) || continue   # review 2026-09-04 A6

        # Precipitation rate at the parcel's level [mm/h].  tpf is collocated
        # with the other mass-level fields (T, q, rho) → same stagger handling.
        I_here = interpolation4Dsingle(cf.tpf, indx, indy,
                                       clamp(state.indl[i], 1.0f0, Float32(nz - 1)),
                                       indt) * 3600.0f0
        (isfinite(I_here) && I_here > p.I_min) || continue

        # Cloud context: fraction at parcel level + column maximum (for the
        # precipitating-area estimate fp).
        cc_prof = ccbufs[tid]
        get_vertical_profile!(cc_prof, cf.cc, indx, indy, indt, nz)
        k0 = clamp(Base.unsafe_trunc(Int, state.indl[i]), 1, nz - 1)
        cc_here = cc_prof[k0]
        cc_max  = 0.0f0
        @inbounds for k in 1:nz
            cc_prof[k] > cc_max && (cc_max = cc_prof[k])
        end

        in_cloud = cc_here > p.cc_threshold
        A = in_cloud ? p.A_incloud : p.A_below
        B = in_cloud ? p.B_incloud : p.B_below
        if !in_cloud && nab > 0
            b = Int(state.relbox[i])
            if 1 <= b <= nab && sphys.weta[b] >= 0.0f0
                A = sphys.weta[b]; B = sphys.wetb[b]   # species' own below-cloud A/B
            end
        end
        A <= 0.0f0 && continue

        # Sub-grid precipitating fraction and effective grid-mean coefficient.
        fp   = clamp(cc_max, 0.05f0, 1.0f0)
        lam  = fp^(1.0f0 - B) * A * I_here^B
        frem = -expm1(-lam * dt)                  # fraction removed, ∈ (0,1)

        dm = state.mass[i] * frem
        state.mass[i] -= dm
        removed[tid] += Float64(dm)

        if dep !== nothing
            ii = clamp(round(Int, indx), 1, Ni)
            jj = clamp(round(Int, indy), 1, Nj)
            Atomix.@atomic dep[ii, jj] += Float64(dm)
        end
    end

    Threads.atomic_add!(total_removed, sum(removed))
    return nothing
end
