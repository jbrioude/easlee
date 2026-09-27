module Operators

using Dates
using Random
using KernelAbstractions
# Optional GPU stack: CUDA (CUDABackend / CuArray / CUDA.functional()) is
# loaded ON DEMAND — only when EASLEE_ADVANCE=ka_gpu or EASLEE_GPU=1 (see
# EASLEE.jl for the install hint).  Submodules cannot see the parent module's
# const at include time, so the gate is re-derived locally from ENV (cheap).
const GPU_REQUESTED = lowercase(get(ENV, "EASLEE_ADVANCE", "")) == "ka_gpu" ||
                      get(ENV, "EASLEE_GPU", "0") == "1"
const GPU_ENABLED = GPU_REQUESTED &&
    (try; @eval using CUDA; true
     catch; false end)          # load failure already warned about in EASLEE.jl
# CUDA.functional() without touching the CUDA symbol when the stack is absent
# (GPU_ENABLED === false short-circuits before the CUDA global is looked up).
_gpu_functional() = GPU_ENABLED && CUDA.functional()

# ── Startup GPU status line ──────────────────────────────────────────────────
# One unambiguous line at run start saying whether THIS process is using the
# GPU — requested by a first-time HPC user (2026-09-04) who had the GPU
# allocated to the job (nvidia-smi showed it) but no Julia process on it,
# because the run script did not export EASLEE_ADVANCE.  Being inside a GPU
# allocation and using GPU acceleration are different things; this line tells
# them apart.  On a functional GPU it also touches the device, so the CUDA
# context exists (and shows in nvidia-smi) from the first minute of the run.
const _GPU_INSTALL_HINT = "install CUDA once in your default Julia environment: " *
    "julia -e 'using Pkg; Pkg.add(\"CUDA\")'"

if GPU_ENABLED
    function _gpu_device_summary()
        dev = CUDA.device()
        mem = try; round(Int, CUDA.totalmem(dev) / 2^30); catch; missing; end
        rt  = try; string(CUDA.runtime_version()); catch; "?"; end
        drv = try; string(CUDA.driver_version()); catch; "?"; end
        name = CUDA.name(dev)
        mem === missing ? "$name (CUDA runtime $rt, driver $drv)" :
                          "$name, $(mem) GiB (CUDA runtime $rt, driver $drv)"
    end
else
    _gpu_device_summary() = ""
end

function gpu_status_banner()
    nt  = Threads.nthreads()
    adv = lowercase(get(ENV, "EASLEE_ADVANCE", ""))
    if !GPU_REQUESTED
        @info "GPU: not requested — CPU run on $nt threads " *
              "(set EASLEE_ADVANCE=ka_gpu for the GPU particle advance)"
    elseif !GPU_ENABLED
        @warn "GPU requested (EASLEE_ADVANCE=$adv) but the CUDA package could not be loaded " *
              "— CPU run on $nt threads; $_GPU_INSTALL_HINT"
    elseif !_gpu_functional()
        @warn "GPU requested (EASLEE_ADVANCE=$adv), CUDA loaded, but no usable NVIDIA device/driver " *
              "— CPU run on $nt threads (run julia -e 'using CUDA; CUDA.functional(true)' for the reason)"
    else
        what = adv == "ka_gpu" ? "particle advance on GPU" :
               "advance stays on CPU (EASLEE_ADVANCE not ka_gpu)"
        @info "GPU backend enabled: $(_gpu_device_summary()); CUDA initialized — " *
              "$what, $nt CPU threads"
    end
    return nothing
end
using Atomix
using GRIB   # PT (isentropic-level wind) reader in physics/isentropic.jl
using ..State
using ..Config
using ..MetAdapter

# Grid/level types live in EASLEE.Types — bring them in for dispatch.
using ..Types: Lambert, Rectilinear, Mercator, AbstractGrid
using ..Types: SigmaPressureLevel, AbstractLevel
using ..Types: Thomson, Hanna, NoTurb
using ..Types: MetFields3D, MetFields2D, ConvFields, DiagConvFields
using ..Types: RunConfig, NestRunConfig, OutputConfig

# Include ported src/ kernels directly into this submodule's namespace.
# Order matters: types → interp → coords → turbulence → rng → turbulence_4d → advance → boundary → conc
include("../met/interp.jl")
include("../physics/surface_layer.jl")
include("../met/coords.jl")
include("../physics/turbulence.jl")
include("../physics/turbulence_hanna.jl")
include("../gpu/rng.jl")
include("../physics/turbulence_4d.jl")
include("../physics/advance.jl")
include("../physics/advance_ka.jl")
# CUDA side of the KA advance (_KAGpuBuffers + _ka_gpu_advance!) — GPU runs only.
GPU_ENABLED && include("../physics/advance_ka_gpu.jl")
include("../physics/plume_rise.jl")
include("../physics/plume_phase.jl")
include("../physics/umbrella.jl")
include("../physics/track_nudge.jl")   # observed-track nudging + wind-bias drift (HT)
include("../physics/convection.jl")
include("../physics/wetdep.jl")
include("../physics/aerosol.jl")
include("../physics/species_physics.jl")
include("../physics/ice_physics.jl")   # transport-phase freezing/sedimentation ([species.ice])
include("../physics/dust_kok.jl")      # mineral-dust source (Kok et al. 2014) + Kok (2011) size distribution
include("../physics/cams_optics.jl")   # CAMS aerosol class optics (table, MEC, f(RH); transport side)
include("../physics/seaspray.jl")      # OSSA sea-spray source (Ovadnevaite et al. 2014)
include("../physics/plume_radiation.jl")  # in-plume LW cooling + radiative descent
include("../physics/transport_tags.jl")   # passive air-mass tags ([tracers]): STE, PBL, …
include("../physics/isentropic.jl")
include("boundary.jl")
include("domainfill.jl")
include("../output/conc.jl")
include("../output/cov_kernel.jl")   # flow-evolved per-parcel covariance (kernel="covariance")
include("../output/receptors.jl")
include("../output/conc_ka.jl")

export release_parcels!, initialize_new_parcels!
export VMixClosure, set_vertical_mixing!, vertical_mixing   # Fix C closure
export advance_active_parcels!, advance_active_parcels_ka!, advance_dispatch!, gpu_status_banner
export initialize_new_parcels_gpu!
export deposit_to_grid!, deposit_to_grid_ka!, deposit_conc_gpu!
export set_output_popmap!, update_bandwidths!, update_airm!
export reassign_escaped_parcels!
export set_plume_config!
export mastin_height_km, mastin_mer_kg_s, volcanic_column_rise, volcanic_column_provenance
export volcanic_column_profiles, PlumeParcels, PlumeProfiles, PlumeWashout,
       build_plume_profiles, plume_parcel_init!, plume_parcels_substep!,
       plume_parcels_substep_v2!, plume_umbrella_equilibrate!,
       plume_phase_stats, plume_release_ensemble,
       VentMix, vent_mix_state, vent_mix_means
export UmbrellaIntrusion, umbrella_radius, umbrella_edge_speed,
       umbrella_thickness, apply_umbrella_spread!

# Per-release-box plume-rise config (box index → (scheme,frp,area,profile,bottom)),
# set by build_run from the spec; consumed at release time in initialize_new_parcels!.
const PLUME_CFG = Ref{Dict{Int,NamedTuple}}(Dict{Int,NamedTuple}())
set_plume_config!(d::AbstractDict) = (PLUME_CFG[] = Dict{Int,NamedTuple}(d);
                                      empty!(PLUME_ENS_CACHE);
                                      empty!(UMBRELLA_REG))
# mode="parcels" ensemble cache: (spec release, met-time slot) → distributions.
# The met column only changes per met window, not per sync step.
const PLUME_ENS_CACHE = Dict{Tuple{Int,Int,Int},NamedTuple}()   # (specidx, phase, met-slot)
# umbrella_spread registry: spec release → (intrusion params, expanded box ids).
# Registered at the first release of a parcels-mode box; consumed each sync
# step by apply_umbrella_spread! in the engine loop.
const UMBRELLA_REG = Dict{Tuple{Int,Int},Tuple{UmbrellaIntrusion,Set{Int}}}()  # (specidx, phase)

# Ambient column (z ASL, T, p, q) at a release cell, for the volcanic schemes.
function _volcanic_column_arrays(nestctx, indx, indy, indt)
    g  = nestctx.gridconfig
    ix = clamp(round(Int, indx), 1, Int(g.Ni)); iy = clamp(round(Int, indy), 1, Int(g.Nj))
    m3 = nestctx.metfields3D
    z3 = m3.zlev; th = m3.theta
    nz = size(z3, 3); it = clamp(round(Int, indt), 1, size(z3, 4))
    zcol = Float64[z3[ix, iy, k, it] for k in 1:nz]
    θcol = Float64[th[ix, iy, k, it] for k in 1:nz]
    p0 = 100000.0; κ = 287.05 / 1005.0
    pcol = Float64[m3.puv[ix, iy, k, it] for k in 1:nz]
    all(>(0), pcol) || (pcol = p0 .* (Float64[m3.rho[ix, iy, k, it] for k in 1:nz] .*
                                      287.05 .* θcol ./ p0) .^ (1 / (1 - κ)))
    Tcol = Float64[m3.T[ix, iy, k, it] for k in 1:nz]
    all(>(100), Tcol) || (Tcol = θcol .* (pcol ./ p0) .^ κ)      # T zeroed (AROME) → from θ
    # ambient humidity is second-order for a Tg-scale column (the vent water
    # dominates); MetFields3D carries no qv → dry-ish ambient
    qcol = fill(1.0e-5, nz)
    return zcol, Tcol, pcol, qcol
end

# Injection layer [z_bottom, z_top] (m AGL) for a plume box, from the met column
# at the release cell.  Called once per box at first release (then cached).
function _compute_plume_layer(cfg, nestctx, indx, indy, indt)
    g  = nestctx.gridconfig
    ix = clamp(round(Int, indx), 1, Int(g.Ni)); iy = clamp(round(Int, indy), 1, Int(g.Nj))
    z3 = nestctx.metfields3D.zlev; th = nestctx.metfields3D.theta
    nz = size(z3, 3); it = clamp(round(Int, indt), 1, size(z3, 4))
    zcol = Float64[z3[ix, iy, k, it] for k in 1:nz]
    zagl = zcol .- zcol[1]
    θcol = Float64[th[ix, iy, k, it] for k in 1:nz]
    if cfg.scheme == "volcanic"
        # 1-D moist eruption column through the ambient met profile: injection
        # layer = [NBH, top] (the umbrella spreads at/above neutral buoyancy).
        zc, Tcol, pcol, qcol = _volcanic_column_arrays(nestctx, indx, indy, indt)
        ztop, znbh, ok = volcanic_column_rise(cfg.mer, cfg.temp, cfg.waterfrac,
                                              zc, Tcol, pcol, qcol; z_vent=zc[1])
        ok || @warn "volcanic column collapsed (MER=$(cfg.mer), T0=$(cfg.temp), n_w=$(cfg.waterfrac)) — injecting the buoyant part anyway" maxlog=2
        # return AGL (release machinery treats plume layers as AGL heights)
        return (max(znbh - zc[1], 500.0), max(ztop - zc[1], 1500.0))
    end
    pblh = Float64(nestctx.metfields2D.pblh[ix, iy, it])
    # Many readers (e.g. ERA5_CONV) don't archive a usable blh → derive it from the
    # θ profile so plume rise doesn't depend on a missing field.
    pblh > 50.0 || (pblh = estimate_pbl_height(zagl, θcol))
    return plume_injection_layer(cfg.scheme, cfg.frp, cfg.area, pblh, zagl, θcol; bottom=cfg.bottom)
end

# mode="parcels": Lagrangian plume-phase ensemble through the met column at the
# vent — released trajectories sample their detrainment height + washout
# survival from it.  Computed once per SPEC release (all its expanded boxes
# share the vent and source) per initialize_new_parcels! call (met-time fresh).
function _compute_plume_ensemble(cfg, nestctx, indx, indy, indt)
    zc, Tcol, pcol, qcol = _volcanic_column_arrays(nestctx, indx, indy, indt)
    ens = plume_release_ensemble(cfg.mer, cfg.temp, cfg.waterfrac,
                                 zc, Tcol, pcol, qcol;
                                 washout = cfg.washout ?
                                     PlumeWashout(qc_crit  = Float64(get(cfg, :qc_crit, 5.0e-3)),
                                                  umb_r_med = Float64(get(cfg, :umb_r_med, 0.0)),
                                                  umb_r_sg  = Float64(get(cfg, :umb_r_sg, 2.0)),
                                                  umb_h     = Float64(get(cfg, :umb_h, 2000.0)),
                                                  umb_t     = Float64(get(cfg, :umb_t, 3600.0)),
                                                  solid_r   = Float64(get(cfg, :solid_r, 0.0))) : nothing,
                                 entrainment = Symbol(get(cfg, :entrainment, "meanfield")),
                                 ventmix = get(cfg, :ventmix, nothing),
                                 vapor_cap = Float64(get(cfg, :vapor_cap, Inf)),
                                 # seeded per (release, phase, met-slot): the draw is
                                 # REPRODUCIBLE run-to-run, so single-factor pairs are
                                 # not contaminated by release noise (found 2026-08-04:
                                 # ~1.5-point spread in the >32 km fraction at +1 h)
                                 rng = Random.Xoshiro(hash((0x5EED, cfg.species,
                                                            get(cfg, :phase, 0)))))
    # survival factor per member for THIS box's species (by id; default 1).
    # With a vent mix members are not equal: gases ride the MELT (gasw),
    # water rides the member's lofted source-water share (w0src).
    sid = uppercase(cfg.species)
    w̄0 = sum(ens.w0src) / max(length(ens.w0src), 1)
    fmass = occursin("SO2", sid) ? ens.fso2 .* ens.gasw :
            occursin("HCL", sid) ? ens.fhcl .* ens.gasw :
            # coarse-mode box (ICE_C...): carries the vapour-cap excess — the
            # condensed load that leaves as large crystals (vapor_cap key)
            occursin("ICE_C", sid) ?
                (ens.fumb .+ ens.fcoarse) .* (ens.w0src ./ max(w̄0, 1e-12)) :
            (occursin("H2O", sid) || occursin("WATER", sid)) ?
                ens.fh2o .* (ens.w0src ./ max(w̄0, 1e-12)) :
            ones(length(ens.z_det))
    # early check (trap #1): the release-altitude distribution, at point of use.
    # A release-physics change that leaves these numbers untouched did nothing.
    let zs = sort(ens.z_det ./ 1e3), nm = length(zs)
        sr = Float64(get(cfg, :solid_r, 0.0))
        fcm = sum(ens.fcoarse .+ ens.fumb) / max(length(ens.fcoarse), 1)
        @info "plume ensemble sampled" species=sid tephra_r_mm=round(sr*1e3, digits=2) z_det_p50_km=round(zs[max(1, nm ÷ 2)], digits=2) z_det_p90_km=round(zs[max(1, 9nm ÷ 10)], digits=2) fmass_mean=round(sum(fmass)/nm, sigdigits=3) fcoarse_mean=round(fcm, sigdigits=3)
    end
    # umbrella-intrusion feed: Q_V = MER·amp/ρ(NBH) (amp = mass-flux
    # amplification from the 1-D entrainment bookkeeping) and N at the NBH.
    # With a vent mix the column runs on the AXIS (core) state (same as the
    # ensemble above — the ascending jet is the magmatic core).
    tk0, nw0 = cfg.temp, cfg.waterfrac
    if get(cfg, :ventmix, nothing) !== nothing
        s0 = vent_mix_state(cfg.ventmix, 0.0)
        tk0 = s0.T0; nw0 = s0.qv0
    end
    _, _, _, amp, _, _, _, _ = volcanic_column_provenance(
        cfg.mer, tk0, nw0, zc, Tcol, pcol, qcol)
    k = clamp(searchsortedlast(zc, ens.prof.z_nbh), 2, length(zc) - 1)
    ρnbh = pcol[k] / (287.05 * Tcol[k])
    θk(kk) = Tcol[kk] * (100000.0 / pcol[kk])^0.2857
    N2 = 9.81 / θk(k) * (θk(k + 1) - θk(k - 1)) / max(zc[k + 1] - zc[k - 1], 1.0)
    return (z_agl = ens.z_det .- zc[1], fmass = fmass, qdet = ens.qdet,
            qv_flux = cfg.mer * amp / max(ρnbh, 1e-6),
            nbv = sqrt(max(N2, 1.0e-6)))
end
export inject_boundary_parcels!, compute_domain_mass
export DomainFillBox, build_domainfill, domainfill_init!, domainfill_inject!, domainfill_cull!
export convective_redistribution!, wet_deposition!, WetDepParams, WetDepState
export SpeciesPhysics, apply_decay!, dry_deposition_gas!, any_decay, any_drydep
export IcePhysics, apply_ice!, any_ice, esat_ice, esat_mk, esat_liq_mk, ice_vterm
export SpeciesConvert, apply_convert!, any_convert, apply_convert_thermo!, any_thermo_convert, apply_species_settling!, any_settling
export PlumeRadiation, apply_plume_radiation!, any_prad
export TransportTags, apply_transport_tags!, strat_fraction, pbl_fraction,
       cont_fraction, reset_tags!, any_tags
export OSSA_MODES, ossa_reynolds, ossa_mode_fluxes, ossa_number_flux,
       ossa_mass_flux, ossa_dfdlogd, seawater_kinematic_viscosity
export SeaSprayBox, build_seaspray, apply_seaspray!, _ossa_massflux
export CAMS_AERO, CAMS_AERO_NAMES, cams_aero_optics, _frh_table, _CAMS_AERO_FRH, _cams_aero_frh
export DustBox, build_dust, apply_dust!, dust_flux_kok, dust_threshold_dry, fecan_moisture_factor, kok_dvdlnd, kok_mass_fractions
export diagnostic_convection!, DiagConvParams
export AerosolBins, AerosolModel, AEROSOL_REFRACTIVE
export apply_settling!, dry_deposition_aerosol!, wet_deposition_aerosol!
export compute_aod, settling_velocity, deposition_velocity, mie_qext
export growth_factor, AEROSOL_KAPPA, AEROSOL_RH_GRID
export IsentropicFields, HybridVert, hybrid_vertical!

# v1.1 (share review B2c): the free-slot scan used to restart from 1 for every
# release box every step — O(boxes × used slots).  A monotone cursor remembers
# where the last scan stopped; the scan runs cursor→end first and wraps to
# 1→cursor-1 only when the tail is exhausted (slots freed BELOW the cursor by the
# cull/boundary paths are found by that wrap, so no hook is needed at the
# freeing sites).  The cursor is a process-wide Ref, clamped to the state length
# so several ParcelStates (tests) can share it safely.
const _FREE_CURSOR = Ref(1)
@inline function _first_unused(state::State.ParcelState, nrelease::Int)
    idx = Int[]
    n = length(state.status)
    n == 0 && return idx
    start = clamp(_FREE_CURSOR[], 1, n)
    @inbounds for i in start:n
        if state.status[i] == 0x00
            push!(idx, i)
            length(idx) == nrelease && (_FREE_CURSOR[] = i + 1; return idx)
        end
    end
    # tail exhausted: wrap to the head (slots freed below the cursor)
    @inbounds for i in 1:start-1
        if state.status[i] == 0x00
            push!(idx, i)
            length(idx) == nrelease && (_FREE_CURSOR[] = i + 1; return idx)
        end
    end
    _FREE_CURSOR[] = 1       # fewer free slots than asked: rescan from 1 next time
    return idx
end

# v1.1: number of release steps inside [t1, t2] without materialising
# `t1 .<= tloops .<= t2` (a BitVector over the whole run per box per step).
@inline function _release_nsteps(tloops::AbstractRange{DateTime}, t1::DateTime, t2::DateTime)
    a = searchsortedfirst(tloops, t1); b = searchsortedlast(tloops, t2)
    return max(b - a + 1, 0)
end

function release_parcels!(state::State.ParcelState, ctx::Config.ModelContext, currentdate::DateTime; rng=Random.default_rng())
    rc = ctx.runconfig
    timestep = rc.timestep
    tloops = rc.startdate:abs(timestep):rc.enddate
    rs = ctx.releases
    F = eltype(state.indx)
    I = eltype(state.nest)

    released = 0
    for ibox in eachindex(rs.boxid)
        t1 = rs.timerelease1[ibox]
        t2 = rs.timerelease2[ibox]
        if !(t1 <= currentdate <= t2)
            continue
        end

        nbsteps = _release_nsteps(tloops, t1, t2)     # v1.1: O(log n), no allocation
        nbase = floor(Int, rs.nbpart[ibox] / max(nbsteps, 1))
        nextra = (currentdate - timestep < t1) ? (rs.nbpart[ibox] - nbase * nbsteps) : 0
        nrelease = nbase + nextra
        nrelease <= 0 && continue

        idx = _first_unused(state, nrelease)
        isempty(idx) && continue

        dx = rs.xrelease2[ibox] - rs.xrelease1[ibox]
        dy = rs.yrelease2[ibox] - rs.yrelease1[ibox]

        for i in idx
            # Sample a random position in geographic coordinates (lon/lat)
            xpos = rs.xrelease1[ibox] + rand(rng, F) * dx
            ypos = rs.yrelease1[ibox] + rand(rng, F) * dy

            # Convert to grid indices; assign to the finest nest that contains (xpos,ypos).
            # In the pathnames file, nest 1 is the inner/finest domain (e.g. AROME) and
            # higher-numbered nests are progressively coarser (e.g. ERA5).  We search from
            # finest (nest 1) to coarsest so that parcels inside the fine domain are
            # assigned to it, and only fall through to a coarser nest when the point lies
            # outside the fine domain.
            # Default: use nest 1 coordinates (coarsest fallback is last nest checked below)
            assigned_nest = I(length(ctx.nests))   # fallback = outermost nest
            gridcfg_last = ctx.nests[end].gridconfig
            indy_last, indx_last = latlontoindyindx(gridcfg_last, ypos, xpos)
            assigned_indx = F(indx_last)
            assigned_indy = F(indy_last)
            for n in 1:length(ctx.nests)
                gridcfg = ctx.nests[n].gridconfig
                indy_n, indx_n = latlontoindyindx(gridcfg, ypos, xpos)
                Ni = gridcfg.Ni
                Nj = gridcfg.Nj
                if 1 <= indx_n <= Ni && 1 <= indy_n <= Nj
                    assigned_nest = I(n)
                    assigned_indx = F(indx_n)
                    assigned_indy = F(indy_n)
                    break
                end
            end

            state.relbox[i]    = rs.boxid[ibox]
            state.Dtimeinit[i] = currentdate
            state.Dtime[i]     = currentdate
            state.status[i]    = 0x01
            state.is_new[i]    = true
            state.nest[i]      = assigned_nest
            state.indx[i]      = assigned_indx
            state.indy[i]      = assigned_indy
            state.mass[i]      = rs.mass_per_box[ibox] / max(F(rs.nbpart[ibox]), F(1))
            state.cxx[i]       = F(0)   # covariance kernel: re-init on first touch
            state.cxy[i]       = F(0)
            state.cyy[i]       = F(0)
            state.airm[i]      = F(0)   # ice closure: set by the plume sample / q0
            released += 1
        end
    end
    return released
end

#function initialize_new_parcels!(state::State.ParcelState, ctx::Config.ModelContext, currentdate::DateTime; rng=Random.default_rng())
function initialize_new_parcels!(state::State.ParcelState, ctx::Config.ModelContext, currentdate::DateTime, t_array; rng=Random.default_rng())
    F = eltype(state.indx)
    ctl = F(ctx.runconfig.ctl)
    # Plume-rise injection layer per box, computed once at release (cached here).
    plume_zr = Dict{Int,Tuple{F,F}}()
    have_plume = !isempty(PLUME_CFG[])

    for nestctx in ctx.nests
        n = nestctx.id
        rhoarray = nestctx.metfields3D.rho
        tkearray = nestctx.metfields3D.tke
        zarray   = nestctx.metfields3D.zlev
        theta    = nestctx.metfields3D.theta
        nz = Int32(size(zarray, 3))
        indt = t_array[n]

        for i in eachindex(state.indx)
            if !(state.status[i] == 0x01 && state.is_new[i])
                continue
            end
            state.nest[i] != n && continue
            ibox = Int(state.relbox[i])

            # Plume rise: override the release z-range with the computed injection
            # layer (parcels then fill it via the existing uniform draw below) —
            # or, in mode="parcels", sample the Lagrangian plume-phase ensemble:
            # emergent detrainment height + mass scaled by washout survival.
            if have_plume && haskey(PLUME_CFG[], ibox) &&
               get(PLUME_CFG[][ibox], :mode, "layer") == "parcels"
                cfg = PLUME_CFG[][ibox]
                sidx = get(cfg, :specidx, ibox)
                ph   = get(cfg, :phase, 0)
                key = (sidx, ph, floor(Int, indt))
                ens = get!(PLUME_ENS_CACHE, key) do
                    _compute_plume_ensemble(cfg, nestctx,
                                            state.indx[i], state.indy[i], indt)
                end
                # umbrella-intrusion registration (once per spec release): the
                # gravity-current spreading applied during transport
                if get(cfg, :umbrella, false) && !haskey(UMBRELLA_REG, (sidx, ph))
                    boxset = Set{Int}(b for (b, c) in PLUME_CFG[]
                                      if get(c, :specidx, b) == sidx &&
                                         get(c, :phase, 0) == ph)
                    la, lo = indyindxtolatlon(nestctx.gridconfig,
                                              Float64(state.indy[i]),
                                              Float64(state.indx[i]))
                    t0 = minimum(ctx.releases.timerelease1[b] for b in boxset)
                    # feeding duration = MASS-WEIGHTED release-end age over the
                    # phase boxes: the constant-flux umbrella law holds only
                    # while the column feeds; afterwards the intrusion coasts
                    # as a constant-volume current (R ∝ t^(1/3)) and the
                    # source-fed interior flow dies (see umbrella.jl header)
                    msum = 0.0; tsum = 0.0
                    for b in boxset
                        mb = Float64(ctx.releases.mass_per_box[b])
                        msum += mb
                        tsum += mb * (Dates.value(ctx.releases.timerelease2[b] - t0) / 1000.0)
                    end
                    tfeed = msum > 0.0 ? max(tsum / msum, 60.0) : Inf
                    UMBRELLA_REG[(sidx, ph)] = (UmbrellaIntrusion(
                        lam = Float64(get(cfg, :spread_lambda, 0.2)),
                        qv_flux = ens.qv_flux, nbv = ens.nbv, t0 = t0,
                        vlon = lo, vlat = la,
                        rsrc = Float64(get(cfg, :umbrella_rsource, 10.0e3)),
                        wind_shield = Bool(get(cfg, :umbrella_windshield, false)),
                        umin = Float64(get(cfg, :spread_umin, 5.0)),
                        tmax_s = Float64(get(cfg, :spread_tmax, 43200.0)),
                        hmin = Float64(get(cfg, :spread_hmin, 0.0)),
                        zgain = Float64(get(cfg, :spread_zgain, 0.0)),
                        azgain_n = Float64(get(cfg, :spread_azgain_n, 0.0)),
                        azgain_s = Float64(get(cfg, :spread_azgain_s, 0.0)),
                        az_z0 = Float64(get(cfg, :spread_az_z0, 0.0)),
                        az_z1 = Float64(get(cfg, :spread_az_z1, 0.0)),
                        zg_z0 = Float64(get(cfg, :spread_zg_z0, 26.0e3)),
                        zg_z1 = Float64(get(cfg, :spread_zg_z1, 30.0e3)),
                        tfeed_s = tfeed), boxset)
                    @info "umbrella intrusion registered" release=sidx phase=ph Qv_m3s=round(ens.qv_flux, sigdigits=3) N=round(ens.nbv, sigdigits=3) R_1h_km=round(umbrella_radius(UMBRELLA_REG[(sidx, ph)][1], 3600.0)/1e3, digits=0) tfeed_h=round(tfeed/3600, digits=2) umin=UMBRELLA_REG[(sidx, ph)][1].umin hmin_m=UMBRELLA_REG[(sidx, ph)][1].hmin zgain=UMBRELLA_REG[(sidx, ph)][1].zgain azgain_n=UMBRELLA_REG[(sidx, ph)][1].azgain_n azgain_s=UMBRELLA_REG[(sidx, ph)][1].azgain_s tmax_h=round(UMBRELLA_REG[(sidx, ph)][1].tmax_s/3600, digits=1)
                end
                j = rand(rng, 1:length(ens.z_agl))
                z_release = F(ens.z_agl[j])
                state.mass[i] *= F(ens.fmass[j])
                # ice-microphysics closure: the air mass this parcel represents,
                # from its member's detrainment water fraction (q = mass/airm)
                state.airm[i] = F(Float64(state.mass[i]) /
                                  max(Float64(ens.qdet[j]), 1e-9))
                # floor: the kNN closure may enlarge but never shrink this
                state.airm0[i] = state.airm[i]
                # dome-wide horizontal placement: emergent parcels are spread
                # AREA-UNIFORMLY across the collapse dome (radius = the
                # umbrella r_source) — the real umbrella starts dome-wide, not
                # as a point (HT: the imaged dome is ~40 km radius)
                if get(cfg, :umbrella, false)
                    rdome = Float64(get(cfg, :umbrella_rsource, 10.0e3))
                    if rdome > 0.0
                        rr = rdome * sqrt(rand(rng)); ph = 2π * rand(rng)
                        gcn = nestctx.gridconfig
                        lav, _ = indyindxtolatlon(gcn, Float64(state.indy[i]),
                                                  Float64(state.indx[i]))
                        mxr = 111_320.0 * Float64(gcn.dlon) * max(cosd(lav), 0.05)
                        myr = 111_320.0 * Float64(gcn.dlat)
                        state.indx[i] += F(rr * cos(ph) / mxr)
                        state.indy[i] += F(rr * sin(ph) / myr)
                    end
                end
            else
                if have_plume && haskey(PLUME_CFG[], ibox)
                    zr = get!(plume_zr, ibox) do
                        Tuple{F,F}(F.(_compute_plume_layer(PLUME_CFG[][ibox], nestctx,
                                                           state.indx[i], state.indy[i], indt)))
                    end
                    z1, z2 = zr
                else
                    z1 = ctx.releases.zrelease1[ibox]
                    z2 = ctx.releases.zrelease2[ibox]
                end
                z_release = z1 + rand(rng, F) * (z2 - z1)
            end

            # Clamp to valid grid range before calling interpolation functions.
            # Parcels released outside the domain (indx > Ni) are possible when the
            # release point is at the grid edge; updatepos! will kill them next step.
            safe_indx = clamp(state.indx[i], F(1), F(nestctx.gridconfig.Ni - 1))
            safe_indy = clamp(state.indy[i], F(1), F(nestctx.gridconfig.Nj - 1))

            state.indl[i] = ztoindl(nestctx.levelconfig, safe_indx, safe_indy, z_release, indt, nestctx.metfields3D)
            # ztoindl can return NaN where the zlev column is masked/missing
            # (clamp does NOT sanitise NaN → floor(Int, NaN) would crash).  A
            # release over such a cell defaults to the surface level.
            isfinite(state.indl[i]) || (state.indl[i] = F(1))
            state.indl[i] = clamp(state.indl[i], F(1), F(nz - 1))

            indl1 = clamp(floor(Int, state.indl[i]), 1, Int(nz - 1))
            indl2 = indl1 + 1

            rho = interpolation4Dsingle(rhoarray, safe_indx, safe_indy, state.indl[i], indt)
            z2m = interpolation4Dsingle(zarray, safe_indx, safe_indy, Float32(indl2), indt)
            z1m = interpolation4Dsingle(zarray, safe_indx, safe_indy, Float32(indl1), indt)

            if state.mass[i] <= zero(F)
                state.mass[i] = rho * max(z2m - z1m, F(1e-3))
            end

            tkeprof = get_vertical_profile(tkearray, safe_indx, safe_indy, indt, nz)
            zprof = get_vertical_profile(zarray, safe_indx, safe_indy, indt, nz)
            thetaprof = get_vertical_profile(theta, safe_indx, safe_indy, indt, nz)
            tke_z = interpolation1Dsingle(tkeprof, state.indl[i])
            z = interpolation1Dsingle(zprof, state.indl[i])
            indz = clamp(Base.unsafe_trunc(Int, state.indl[i]), 1, length(zprof) - 1)
            indzp = indz + 1
            yl = turbulent_length_scale_v2(indz, indzp, tke_z, z, zprof, thetaprof)
            tlw, sigw = tke_partition(tke_z, yl)

            state.wturb[i] = randn(rng, F)
            state.dt[i] = min(tlw / ctl / F(5), F(10))
            state.miss[i] = state.indl[i]
            state.miss2[i] = tke_z
            state.miss3[i] = yl
            state.is_new[i] = false
        end
    end
    return nothing
end
function bin_active_parcels_by_nest!(buffers::Vector{Vector{Int32}}, state)
    nbuff = length(buffers)

    for b in buffers
        empty!(b)
    end

    npar = length(state.nest)

    @inbounds for i in 1:npar
        # adapt these conditions to your real activity flags
        # here I assume:
        #   state.nest[i] > 0     => parcel assigned to a nest
        #   state.active[i] == 1  => active parcel
        #
        # If your state uses another flag name, replace it.
        nestid = state.nest[i]

        if nestid >= 1 && nestid <= nbuff
            # replace this with your true activity test
            if state.active[i] == 1
                push!(buffers[nestid], Int32(i))
            end
        end
    end

    return nothing
end

function advance_active_parcels!(state, ctx, currentdate, indt_array)
    timestep = ctx.runconfig.timestep
    nnest = length(ctx.nests)

    for n in 1:nnest
        # Single findall — one allocation, returns Vector{Int} directly.
        # Avoids 3 temporary BitVectors from the (.== currentdate) .&& (.== n) dance.
        n8 = Int8(n)
        idx = findall(i -> state.Dtime[i] == currentdate && state.nest[i] == n8,
                      eachindex(state.Dtime))

        isempty(idx) && continue

        nestctx = ctx.nests[n]

        advance!(
            state,
            timestep,
            nestctx.nestrunconfig,      # typed NestRunConfig — no Dict lookup
            nestctx.metfields3D,
            nestctx.metfields2D,
            nestctx.gridconfig,
            nestctx.levelconfig,
            indt_array[n],
            idx,
        )
    end

    return nothing
end

# ── Advance-backend dispatch ──────────────────────────────────────────────────
# EASLEE_ADVANCE selects how the parcel advance runs:
#   "threads" (default) — multithreaded CPU path (advance!/updatepos!), full
#                         physics incl. Hanna; convection handled by the caller.
#   "ka_cpu"            — KernelAbstractions kernel on the CPU() backend.
#   "ka_gpu"            — same kernel on CUDABackend() (NVIDIA GPU).
# The KA kernel supports Thomson / no-turbulence only (Hanna is CPU-only) and
# does NOT include convection — so ka_* are for turbulence/advection studies and
# GPU benchmarking, not for convective runs.
const ADVANCE_MODE = Ref(:unset)
const _WIND_WARNED = Ref(false)   # warn once if 2D/1D wind forces the threads path
const _KA_PHYS_WARNED = Ref(false) # warn once if a threads-only switch forces the threads path

# review 2026-09-04 A3: physics the KA kernel does not implement.  Returns the
# name of the first armed switch (any nest), or nothing.  Only on-the-fly Hanna
# is turbulence-related: Hanna(false) runs the Thomson body on both paths.
function _ka_unsupported_switch(ctx)
    for nc in ctx.nests
        rc = nc.nestrunconfig
        # 2026-09-21: on-the-fly Hanna now runs on the KA path (turbulence_hanna_4d
        # + the three surface arrays), so it is no longer a threads-only switch.
        # 2026-09-21: mesoscale, kh_ft, horizontal_tke, horizontal_hanna and the
        # Monin-Obukhov surface wind all run on the KA path now — nothing in this
        # list forces the threads advance any more.  Keep the function: the next
        # physics switch added to the CPU advance belongs here until it is ported.
    end
    return nothing
end

# true when the user set EASLEE_ADVANCE explicitly — then an unsupported switch is
# an ERROR, not a silent downgrade (2026-09-21).
const _ADVANCE_REQUESTED = Ref(false)

function _advance_mode()
    if ADVANCE_MODE[] === :unset
        v = lowercase(get(ENV, "EASLEE_ADVANCE", "threads"))
        _ADVANCE_REQUESTED[] = haskey(ENV, "EASLEE_ADVANCE") && v in ("ka_gpu", "ka_cpu")
        m = v == "ka_gpu" ? :ka_gpu : v == "ka_cpu" ? :ka_cpu : :threads
        if m === :ka_gpu && !_gpu_functional()
            hint = GPU_ENABLED ? "" : " ($_GPU_INSTALL_HINT)"
            error("EASLEE_ADVANCE=ka_gpu was requested but CUDA is not functional$hint. " *
                  "Refusing to run on the CPU under a GPU label — set EASLEE_ADVANCE=threads " *
                  "if that is what you want.")
        end
        ADVANCE_MODE[] = m
        @info "advance backend selected" mode=m
    end
    ADVANCE_MODE[]
end

function advance_dispatch!(state, ctx, currentdate, indt_array)
    m = _advance_mode()
    # Wind dimensionality (2D/1D) is honoured only by the threads advance
    # (updatepos! reads nestrunconfig.wind_mask). The KA kernels don't yet, so
    # force the threads path when any dimension is frozen — correctness first.
    if m !== :threads && !all(ctx.nests[1].nestrunconfig.wind_mask)
        _ADVANCE_REQUESTED[] && error("EASLEE_ADVANCE=$m was requested but wind " *
            "dimensionality (2D/1D) requires the threads advance.  Refusing to downgrade silently.")
        if !_WIND_WARNED[]
            @warn "wind dimensionality (2D/1D) requires the threads advance — overriding $m with threads"
            _WIND_WARNED[] = true
        end
        m = :threads
    end
    # review 2026-09-04 A3: the KA kernel silently DROPPED these switches (no
    # Hanna on-the-fly, no mesoscale OU term, no K_h random walk, no MO surface
    # layer) — a run asking for them ran different physics without notice.
    #
    # 2026-09-21 (Jerome): a silent downgrade is a bug factory.  `surface_layer_wind`
    # DEFAULTS TO TRUE, so for two weeks every "GPU" run of a default config was
    # actually the threads advance, evidenced only by one warning line among
    # hundreds — it wasted an afternoon here on 21 Sept and it would waste a
    # reader's day in the paper's benchmark table.  When the user ASKED for a
    # backend, we now refuse rather than quietly run something else.  The auto
    # default (no EASLEE_ADVANCE set) still falls back, because nobody asked.
    if m !== :threads
        sw = _ka_unsupported_switch(ctx)
        if sw !== nothing
            if _ADVANCE_REQUESTED[]
                error("EASLEE_ADVANCE=$m was requested but $sw is implemented by the " *
                      "threads advance only.  Either switch it off in the config, or run " *
                      "with EASLEE_ADVANCE=threads.  (Refusing to downgrade silently: a " *
                      "run that says GPU must BE GPU.)")
            end
            if !_KA_PHYS_WARNED[]
                @warn "$sw is implemented by the threads advance only — overriding $m with threads"
                _KA_PHYS_WARNED[] = true
            end
            m = :threads
        end
    end
    if m === :threads
        advance_active_parcels!(state, ctx, currentdate, indt_array)
    elseif m === :ka_cpu
        advance_active_parcels_ka!(state, ctx, currentdate, indt_array;
                                   backend = KernelAbstractions.CPU(; static = true))
    else  # :ka_gpu
        advance_active_parcels_ka!(state, ctx, currentdate, indt_array;
                                   backend = CUDABackend())
    end
    return nothing
end
"""
    reassign_escaped_parcels!(state, ctx, currentdate)

After `advance_active_parcels!`, parcels that exited their current (inner) nest
have `Dtime != currentdate+|ts|`.  For each such parcel, convert its grid position to
geographic coordinates and check whether a coarser outer nest contains the point.
If yes, update the parcel's nest assignment and grid indices, then reset Dtime to
`currentdate + |timestep|` so it is picked up by the outer nest at the next step.

Parcels that exit the outermost nest remain with Dtime == DateTime(1/2) (truly dead).
"""
function reassign_escaped_parcels!(state::State.ParcelState, ctx::Config.ModelContext, currentdate::DateTime)
    nnest = length(ctx.nests)
    F = eltype(state.indx)
    I = eltype(state.nest)
    timestep    = ctx.runconfig.timestep
    # SIGNED step: a successful advance sets Dtime = currentdate + timestep.
    # In backward runs timestep < 0, so using abs(timestep) here made every
    # survivor (Dtime = currentdate − |ts|) look escaped and freed it every
    # step — backward transport never advanced past the release step.
    next_date   = currentdate + timestep

    for i in eachindex(state.status)
        state.status[i] == 0x01 || continue      # only active parcels
        state.Dtime[i]  == next_date && continue  # still alive in its nest

        current_nest = Int(state.nest[i])

        if current_nest >= nnest
            # Already in outermost nest and escaped it — or single-nest run.
            # Free the slot so it can be reused (boundary injection, new releases).
            state.status[i] = 0x00
            continue
        end

        # ── Convert current (out-of-bounds) grid position to lat/lon ─────────
        src_gc = ctx.nests[current_nest].gridconfig
        lat, lon = indyindxtolatlon(src_gc, Float64(state.indy[i]), Float64(state.indx[i]))

        # ── Search coarser nests (inner → outer) for one that contains (lat,lon) ─
        new_nest  = 0
        new_indx  = F(0)
        new_indy  = F(0)
        for n in (current_nest + 1):nnest
            tgt_gc  = ctx.nests[n].gridconfig
            indy_n, indx_n = latlontoindyindx(tgt_gc, lat, lon)
            Ni = tgt_gc.Ni
            Nj = tgt_gc.Nj
            # Wrap longitude for global grids (Rectilinear only; Lambert.globalx is always false)
            if tgt_gc.globalx && tgt_gc isa Rectilinear
                indx_n < 1   && (indx_n += 360.0 / tgt_gc.dlon)
                indx_n >= Ni && (indx_n -= 360.0 / tgt_gc.dlon)
            end
            if 1 <= indx_n <= Ni && 1 <= indy_n <= Nj
                new_nest = n
                new_indx = F(indx_n)
                new_indy = F(indy_n)
                break
            end
        end

        if new_nest == 0
            # Exited all nests — free the slot.
            state.status[i] = 0x00
            continue
        end

        # ── Clamp vertical level to the new nest's range ─────────────────────
        Nk_new  = Int32(size(ctx.nests[new_nest].metfields3D.zlev, 3))
        new_indl = clamp(state.indl[i], F(1), F(Nk_new - 1))

        state.nest[i]  = I(new_nest)
        state.indx[i]  = new_indx
        state.indy[i]  = new_indy
        state.indl[i]  = new_indl
        # The stored turbulent velocities were drawn under the OLD nest's
        # closure and its sigma; on the new nest they are meaningless numbers of
        # the wrong scale (AROME's native-e sigma_u vs ERA5's Hanna-fill sigma_w).
        # Zero them: the OU process re-establishes the correct variance within
        # one tau_L, which is short against the time a parcel spends in the
        # outer nest, and zero is the same state a freshly released parcel has.
        state.uturb[i] = F(0)
        state.vturb[i] = F(0)
        state.wturb[i] = F(0)
        # Set Dtime to next_date so the parcel is picked up at the next step
        # (matching what a successful advance!() would have produced).
        state.Dtime[i] = next_date
    end
    return nothing
end

function deposit_to_grid!(ctx::Config.ModelContext, state::State.ParcelState)
    rc      = ctx.runconfig
    oc      = ctx.output
    t_array = MetAdapter.met_time_index_array(rc.currentdate, ctx.nests)
    met3    = [n.metfields3D for n in ctx.nests]
    gc      = [n.gridconfig  for n in ctx.nests]
    lc      = [n.levelconfig for n in ctx.nests]
    conccalc!(oc, state, rc, t_array, met3, gc, lc, false)
    return nothing
end

end
