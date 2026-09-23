
# ── Typed run-configuration containers ───────────────────────────────────────
#
# RunConfig     — global simulation parameters (direction, dates, timestep …)
#                 plus mutable state updated every step (currentdate, dateoutput).
# NestRunConfig — per-nest mutable state: time-bracket indices into metfile
#                 array, turbulence scheme, and the small set of scalars that
#                 advance! needs (ctl, direction).
#
# Both are `mutable struct` because a handful of fields are updated in-place
# each step (currentdate, t_index, turbulence).  All other fields are set once
# at startup and then read-only in practice.
#
# Thomson, Hanna, NoTurb are defined in src/types/grid.jl (included first).

using Dates

# ── Global run parameters ─────────────────────────────────────────────────────

mutable struct RunConfig
    # Static (set once at startup)
    direction::Int8           # +1 forward, -1 backward
    startdate::DateTime
    enddate::DateTime
    timestep::Millisecond
    outstep::Millisecond
    ctl::Float32              # Lagrangian time-scale factor (global default)
    nest::Int                 # number of nests
    nbreleasebox::Int
    Nc::Int                   # number of output classes (= release boxes)
    iout::Int8                # output kind: 1 concentration, 2 mixing ratio, 3 both
    # Mutable (updated each step)
    currentdate::DateTime
    dateoutput::DateTime
end

# ── Per-nest run state ────────────────────────────────────────────────────────

mutable struct NestRunConfig
    # Mutable
    t_index::NTuple{2,Int}           # bracket indices into metconfig.filedate
    turbulence::Union{Thomson,Hanna,NoTurb} # turbulence scheme for this nest
    # Static (per-nest scalars needed by advance!)
    ctl::Float32
    direction::Int8
    # Wind dimensionality: which parcel indices the advance may update —
    # (update_indx, update_indy, update_indl).  3D = (true,true,true); freezing an
    # index keeps each parcel's own value → 2D/1D runs (mass-conservation tests).
    wind_mask::NTuple{3,Bool}
    # Mesoscale velocity fluctuations (FLEXPART usig/vsig/wsig): per-parcel OU
    # process with memory r = exp(−2Δt/windint_s) and amplitude
    # fturbmeso · σ(surrounding grid-point winds) — the sub-grid wind
    # uncertainty the archived met cannot resolve (see turbulence_mesoscale).
    mesoscale::Bool
    fturbmeso::Float32
    windint_s::Float32               # met archive interval [s] (OU decorrelation)
    # FLEXPART free-atmosphere horizontal diffusivity D_h [m²/s]
    # ([physics.turbulence].kh_ft): memoryless random walk √(2·K·Δt)·N(0,1)
    # applied ABOVE the PBL only (the Hanna/Thomson BL owns turbulence below;
    # above it they supply σ_w only — K_v = σ_w²·τ_Lw = 0.1 m²/s = FLEXPART
    # d_strat — and no horizontal term at all).  0 = off (legacy).
    # TRAJECTORY route: [physics.mixing] must stay OFF when this is on —
    # one K either moves parcels or exchanges mass, never both.
    kh_ft::Float32
    # TKE-driven HORIZONTAL turbulence ([physics.turbulence].horizontal_tke).
    # Native-e drivers (AROME native e, MesoNH) archive the full
    # e = ½(σ_u²+σ_v²+σ_w²), so the horizontal variance is already in the field
    # the vertical Langevin step reads: σ_u² = σ_v² = 2·e·(1−r_w)/2, τ_Lh =
    # K_h/σ_u² (tke_partition_h).  Unlike kh_ft this is an OU velocity WITH
    # memory, applied at every level (not only above the PBL), and it is the
    # only horizontal turbulent term the trajectories have.  false = legacy.
    # TRAJECTORY route: [physics.mixing] must stay OFF when this is on.
    horiz_tke::Bool
    # Hanna horizontal Langevin ([physics.turbulence].horizontal_hanna): σ_u, σ_v,
    # T_Lu, T_Lv from Hanna (1982) in the PBL, on Hanna-fill drivers (ERA5, GFS,
    # MERRA-2, AROME-HDF) — FLEXPART's in-PBL horizontal turbulence.  The
    # native-e counterpart is horiz_tke; build_run arms each per nest.
    horiz_hanna::Bool
    # Per-nest vertical-mixing closure (Fix C), taken from THIS nest's variables
    # namelist (#info tke_kind / vmix_rw / vmix_ck).  Nests may legitimately
    # differ — AROME (native e) inside ERA5 (Hanna fill, tke = 0.5*sigma_w^2) is
    # the TROPIVOC 48-h setup — because nothing crosses a nest boundary except
    # the parcel itself: each nest's turbulence is derived from its own met, so
    # a run-global closure was never physically required.  Stored as scalars
    # rather than a VMixClosure because types/config.jl is included before
    # physics/turbulence.jl defines that struct; advance rebuilds it (isbits).
    vmix_native::Bool
    vmix_rw::Float32
    vmix_ck::Float32
    # Polar-stereographic caps ([physics].polar_stereo_lat — FLEXPART
    # switchnorth/switchsouth): poleward of the switch row the horizontal
    # position update runs in a pole-tangent stereographic plane (see
    # polar_stereo_update in met/coords.jl).  Only ever true on a global
    # Rectilinear grid whose latitude rows reach the pole (build_run guards).
    polar_n::Bool
    polar_s::Bool
    switchnorth_j::Float32           # indy above which the north cap applies
    switchsouth_j::Float32           # indy below which the south cap applies
    # Monin-Obukhov surface-layer wind ([physics].surface_layer_wind).  Below
    # the lowest model level the archived wind is constant-extended to the
    # ground (metcalc_models.jl copies level 2 into the surface slot), so a
    # parcel at 2 m carries the ~20 m wind.  With this on, the wind is put back
    # on the similarity profile anchored at that same level.
    mo_surface::Bool
end

# ── Constructors from legacy Dicts ────────────────────────────────────────────

function RunConfig(rc::Dict, Nc::Int)
    dir = Int8(rc["direction"])
    sd  = rc["startdate"]
    RunConfig(
        dir,
        sd,
        rc["enddate"],
        Millisecond(rc["timestep"]),
        Millisecond(abs(rc["outstep"])),
        Float32(rc["ctl"]),
        rc["nest"],
        rc["nbreleasebox"],
        Nc,
        Int8(get(rc, "iout", 1)),
        get(rc, "currentdate", sd),   # may already be set
        get(rc, "dateoutput",  sd),
    )
end
