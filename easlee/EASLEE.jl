
# ── EASLEE — the Lagrangian model engine module ───────────────────────────────
# Core engine: typed grid/met/state structures, met adapters, transport +
# physics operators, the inter-trajectory mixing schemes, and the time loop.
# NOTE: the run-configuration layer (src/config/EasleeConfig.jl) is a SEPARATE
# standalone module loaded on its own; the Main-scope bridge (runtime.jl +
# build_run.jl) ties the two together — the runtime path is
# read_runspec(TOML) → build_run / run_spec!.
#
# This is the TRANSPORT distribution of EASLEE: the reactive-chemistry stack
# (mechanisms, solvers, heterogeneous chemistry, emissions, adjoint) is not
# part of it.  Inter-trajectory mixing ([physics.mixing]) and the kNN
# adaptive-kernel output estimator ARE included.

# Scripts include this file into Main, so nothing is precompiled anyway; the
# directive keeps that explicit and stops Julia caching a stale image.
__precompile__(false)

module EASLEE

using Dates
using Random
using StructArrays

# ── Optional GPU stack (CUDA) — loaded ON DEMAND ─────────────────────────────
# Requested via EASLEE_ADVANCE=ka_gpu (GPU advance kernel).  CUDA is NOT a
# dependency of the easlee project: a CPU-only user installs and loads nothing.
# A GPU user installs it ONCE in the DEFAULT Julia environment
#     julia -e 'using Pkg; Pkg.add("CUDA")'
# and environment stacking makes it loadable from any project.  When CUDA loads
# but is non-functional, the advance dispatch in core/operators.jl routes every
# call back to the CPU path automatically.
const GPU_REQUESTED = lowercase(get(ENV, "EASLEE_ADVANCE", "")) == "ka_gpu" ||
                      get(ENV, "EASLEE_GPU", "0") == "1"
const GPU_ENABLED = if GPU_REQUESTED
    try
        @eval using CUDA
        true
    catch err
        @warn """EASLEE: GPU requested (EASLEE_ADVANCE=$(get(ENV, "EASLEE_ADVANCE", "")) EASLEE_GPU=$(get(ENV, "EASLEE_GPU", "0"))) but CUDA could not be loaded — CPU paths used instead.
                 Install it once in your DEFAULT Julia environment (it is then loadable from any project):
                     julia -e 'using Pkg; Pkg.add("CUDA")'""" exception=err
        false
    end
else
    false
end

# Types submodule — submodules reference these via `using ..Types: Lambert, ...`
module Types
    using Dates   # needed for MetConfig.filedate::Vector{DateTime}
    include("src/types/grid.jl")
    include("src/types/level.jl")
    include("src/types/met.jl")
    include("src/types/config.jl")   # RunConfig, NestRunConfig
    include("src/types/output.jl")   # OutputConfig
end
using .Types   # re-export into EASLEE namespace

include("src/core/state.jl")
include("src/core/bench.jl")
include("src/core/config.jl")
include("src/core/metadapter.jl")

include("src/core/operators.jl")
include("src/core/timeloop.jl")

using .State
using .Bench
import .Config: NestContext, ModelContext, build_context, MetConfig, RunConfig, NestRunConfig, OutputConfig

using .MetAdapter
using .Operators
# Grid-cell mapping helpers used by the gridded chem output (chemconc.jl). They
# live in Operators (met/coords.jl, output/conc.jl) and are not exported, so the
# blanket `using .Operators` above does not bring them into EASLEE scope.
using .Operators: indxtoxout, indytoyout, indyindxtolatlon, indltozlev, find_out_layer,
                  interpolateW, interpolate2Dsingle
using .TimeLoop

# ── Inter-trajectory mixing ([physics.mixing]) ───────────────────────────────
# Tracer mass exchange BETWEEN Lagrangian parcels — the alternative to a
# horizontal random walk (one K either moves parcels or exchanges mass between
# them, never both).  See the [physics.mixing] section of docs/CONFIG.md.
include("src/physics/mixing.jl")            # MixingConfig, apply_mixing!/_xyz!
include("src/physics/puff_mixing.jl")       # geometric puff-covariance mixing (v2)
include("src/output/trajectory.jl")         # per-trajectory CSV dump

import .Types: ConvFields, DiagConvFields

export ParcelState, ReleaseSchedule, OutputState, ModelContext,
       NestContext, MetConfig, RunConfig, NestRunConfig, OutputConfig,
       build_context,
       release_parcels!, initialize_new_parcels!, advance_active_parcels!,
       step_model!

# ── Fix C: driver-consistent vertical mixing closure (σ_w² / K_z semantics) ──
export VMixClosure, set_vertical_mixing!, vertical_mixing

# ── Domain fill (release kind = "domainfill") ─────────────────────────────────
export DomainFillBox, build_domainfill, domainfill_init!, domainfill_inject!,
       domainfill_cull!, inject_boundary_parcels!, compute_domain_mass

# ── Convection & wet deposition ───────────────────────────────────────────────
export ConvFields, convective_redistribution!, wet_deposition!,
       WetDepParams, WetDepState
export DiagConvFields, diagnostic_convection!, DiagConvParams

# ── Aerosol transport ─────────────────────────────────────────────────────────
export AerosolBins, AerosolModel, AEROSOL_REFRACTIVE
export apply_settling!, dry_deposition_aerosol!, wet_deposition_aerosol!
export compute_aod, settling_velocity, deposition_velocity, mie_qext
export growth_factor, AEROSOL_KAPPA, AEROSOL_RH_GRID

# ── Hybrid kinematic/isentropic vertical transport ────────────────────────────
export IsentropicFields, HybridVert, hybrid_vertical!

# ── Benchmarking ──────────────────────────────────────────────────────────────
export BenchTimer, bench_report, timeit!

# ── Inter-trajectory mixing exports ([physics.mixing]) ───────────────────────
export MixingConfig, apply_mixing_xyz!, apply_mixing_xyz_substeps!
export MixingPairlist, build_pairlist_xyz, needs_rebuild
# geometric puff-covariance mixing (puff_mixing.jl)
export Sym22, PuffCovarianceState, PuffMixingConfig
export apply_puff_mixing!, puff_shell_volume, overlap3d, delta_overlap
export run_puff_mixing_tests

end
