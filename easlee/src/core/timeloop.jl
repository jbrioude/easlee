# ──────────────────────────────────────────────────────────────────────────────
# core/timeloop.jl — one model time step (TimeLoop module)
#
# step_model! sequences the operator-split step: release → output deposit →
# advection → nest reassignment → optional hybrid vertical / convection /
# wet deposition / aerosol / per-species physics.
# FLEXPART lineage: timemanager.
# ──────────────────────────────────────────────────────────────────────────────

module TimeLoop

using Dates
using Random
using ..State
using ..Config
using ..Operators
using ..MetAdapter
using ..Bench

export step_model!

"""
    step_model!(state, ctx, currentdate; do_output, rng,
                convfields, wetdepparams, wetdepstate)

Advance the model by one time step
(release → deposit → advect → convection → wet deposition → species physics).

Optional convection / wet-deposition arguments:
- `convfields`   : Vector (one per nest) of `ConvFields` or `nothing`.  When
                   provided, parcels in convective columns are redistributed
                   by the archived mass-flux scheme after advection, and wet
                   scavenging is applied from the 3-D precipitation flux.
- `wetdepparams` : a `WetDepParams`; pass `nothing` to run convection without
                   wet deposition.
- `wetdepstate`  : optional `WetDepState` accumulating deposited mass on the
                   met grid.

`t_abs` (absolute time [s since midnight UTC]) is accepted for interface
compatibility and is unused by the transport distribution.
"""
function step_model!(state::State.ParcelState, ctx::Config.ModelContext, currentdate::DateTime;
                     t_abs::Float64 = 0.0,
                     do_output::Bool = true,
                     rng = Random.default_rng(),
                     convfields = nothing,
                     wetfields = nothing,
                     diagconv = nothing,
                     diagconvparams = nothing,
                     wetdepparams = nothing,
                     wetdepstate = nothing,
                     aerosol = nothing,
                     aerosol_drydep::Bool = true,   # review 2026-09-04 A6: [species].dry_deposition / [physics.deposition].dry
                     aerosol_wetdep::Bool = true,   # review 2026-09-04 A6: [species].wet_deposition / [physics.deposition].wet
                     hybridvert = nothing,
                     species_physics = nothing,
                     ice_physics = nothing,
                     convert_physics = nothing,
                     plume_radiation = nothing,
                     transport_tags = nothing,
                     do_airm::Bool = false,          # [physics].prognostic_airm
                     airm_boxes = nothing,           # release boxes concerned (Set)
                     bench = nothing)
    # Update typed RunConfig.
    ctx.runconfig.currentdate = currentdate

    t_array = MetAdapter.met_time_index_array(currentdate, ctx.nests)

    released = Bench.timeit!(bench, :release) do
        r = Operators.release_parcels!(state, ctx, currentdate; rng=rng)
        r > 0 && Operators.initialize_new_parcels!(state, ctx, currentdate, t_array; rng=rng)
        r
    end

    # deposit BEFORE advance so conccalc! sees Dtime == currentdate (matching old main.jl order)
    do_output && Bench.timeit!(() -> Operators.deposit_to_grid!(ctx, state), bench, :output)

    Bench.timeit!(() -> Operators.advance_dispatch!(state, ctx, currentdate, t_array),
                  bench, :advection)

    # Nest-switching: parcels that exited their inner nest are reassigned to the
    # next coarser nest (if one exists) that contains their current position.
    Bench.timeit!(() -> Operators.reassign_escaped_parcels!(state, ctx, currentdate),
                  bench, :reassign)

    # ── Optional hybrid kinematic/isentropic vertical transport ───────────────
    if hybridvert !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.hybrid_vertical!(hybridvert, state, ctx, dt_s, t_array),
                      bench, :isentropic)
    end

    # ── Optional convection + wet deposition (operator-split after advection) ─
    if convfields !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)  # ms → s
        Bench.timeit!(() -> Operators.convective_redistribution!(state, ctx, convfields, dt_s, t_array;
                      conv_count = transport_tags === nothing ? nothing : transport_tags.conv_count),
                      bench, :convection)
    elseif diagconv !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.diagnostic_convection!(state, ctx, diagconv, dt_s, t_array;
            params = diagconvparams === nothing ? Operators.DiagConvParams() : diagconvparams),
            bench, :convection)
    end
    # Wet deposition needs tpf/cc only, never umf (2026-10-13: it ran inside the
    # convection branch and so vanished on datasets without a mass flux).
    # Callers that pass only convfields keep the old behaviour.
    wf = wetfields === nothing ? convfields : wetfields
    if wf !== nothing && wetdepparams !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.wet_deposition!(state, ctx, wf, dt_s, t_array;
                      params=wetdepparams, wdstate=wetdepstate, sphys=species_physics),
                      bench, :wetdep)
    end

    # ── Optional aerosol processes (settling → dry dep → size-dependent wet dep)
    if aerosol !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.apply_settling!(aerosol, state, ctx, dt_s, t_array), bench, :settling)
        # review 2026-09-04 A6: deposition ran regardless of the config switches; now gated.
        aerosol_drydep &&
            Bench.timeit!(() -> Operators.dry_deposition_aerosol!(aerosol, state, ctx, dt_s, t_array), bench, :drydep)
        if wf !== nothing && aerosol_wetdep
            Bench.timeit!(() -> Operators.wet_deposition_aerosol!(aerosol, state, ctx, wf, dt_s, t_array),
                          bench, :wetdep)
        end
    end

    # ── Optional per-species processes (decay + gas dry deposition) ──────────
    if species_physics !== nothing
        dt_s = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Operators.any_decay(species_physics) &&
            Bench.timeit!(() -> Operators.apply_decay!(state, species_physics, dt_s),
                          bench, :decay)
        Operators.any_drydep(species_physics) &&
            Bench.timeit!(() -> Operators.dry_deposition_gas!(state, ctx, species_physics,
                          dt_s, t_array), bench, :drydep)
        # [species.settling] : chute gravitationnelle d'une espèce condensée.
        # thetap : la sédimentation est CROSS-ISENTROPIQUE — sans rebasage du θ
        # prognostique le pas hybride annule exactement la chute (voir docstring)
        # rice : rayon PROPRE de la parcelle ([species.convert].sectional) —
        # la sédimentation devient résolue en taille (un cristal sublimé tombe
        # moins vite) au lieu d'être figée sur le r_eff nominal de la boîte
        ricev = (convert_physics !== nothing && !isempty(convert_physics.rice)) ?
                convert_physics.rice : nothing
        Operators.any_settling(species_physics) &&
            Bench.timeit!(() -> Operators.apply_species_settling!(state, ctx,
                          species_physics, dt_s, t_array;
                          thetap = hybridvert === nothing ? nothing : hybridvert.theta_p,
                          rice = ricev),
                          bench, :settling)
    end
    # [species.convert]: stochastic in-flight species conversion (before the
    # per-species processes so the daughter's physics applies this step)
    # T′ radiatif ([physics.plume_radiation].tprime_tau) : les gates de
    # saturation évaluent e_si à T_ambiant + T′ (anomalie froide du panache) —
    # convert tourne AVANT apply_plume_radiation!, donc T′ du pas précédent
    tanomv = (plume_radiation !== nothing && !isempty(plume_radiation.tanom)) ?
             plume_radiation.tanom : nothing
    if convert_physics !== nothing && Operators.any_convert(convert_physics)
        dt_c = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        if Operators.any_thermo_convert(convert_physics)
            # gate thermodynamique (freeze / sublimate) : a besoin de T,p locales
            Bench.timeit!(() -> Operators.apply_convert_thermo!(state, ctx,
                          convert_physics, dt_c, t_array; tanom = tanomv),
                          bench, :convert)
        else
            Bench.timeit!(() -> Operators.apply_convert!(state, convert_physics, dt_c),
                          bench, :convert)
        end
    end
    # [physics].prognostic_airm: refresh the air mass each parcel represents
    # from the CURRENT parcel spacing (kNN) — the dilution is read from the
    # trajectories, not prescribed.  MUST run before the two consumers of
    # q = mass/airm below (plume_radiation, then ice).
    if do_airm
        # fill the per-parcel layer thickness for whichever consumer holds it:
        # [species.ice] when armed, else plume_radiation (they SHARE the same
        # vector — build_run allocates it once from prognostic_airm).  Without
        # this fallback a convert/fall ice run leaves hv empty and the LW
        # cooling silently uses the constant h_ref (handoff §2.3).
        hvv = (ice_physics !== nothing && !isempty(ice_physics.hv)) ?
              ice_physics.hv :
              (plume_radiation !== nothing && !isempty(plume_radiation.hv)) ?
              plume_radiation.hv : nothing
        Bench.timeit!(() -> Operators.update_airm!(state, ctx.nests[1].gridconfig,
                      ctx.nests[1].metfields3D, t_array[1], 1; boxes = airm_boxes,
                      hv_out = hvv), bench, :airm)
    end
    # [physics.plume_radiation]: LW cooling from excess plume water (+descent)
    if plume_radiation !== nothing && Operators.any_prad(plume_radiation)
        dt_r = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        # hybrid/isentropic vertical: the radiative descent must be applied to
        # the prognostic parcel θ (else the next hybrid step undoes the move)
        θp = hybridvert === nothing ? nothing : hybridvert.theta_p
        Bench.timeit!(() -> Operators.apply_plume_radiation!(state, ctx, plume_radiation,
                      dt_r, t_array; thetap = θp), bench, :radiation)
    end
    # [tracers]: passive transport-history tags (strat residence, …) — diagnostic
    if transport_tags !== nothing
        dt_t = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.apply_transport_tags!(state, ctx, transport_tags,
                      dt_t, t_array), bench, :tags)
    end
    # [species.ice]: freezing + sedimentation of tagged water, after transport
    if ice_physics !== nothing && Operators.any_ice(ice_physics)
        dt_i = Float64(abs(ctx.runconfig.timestep.value) / 1000.0)
        Bench.timeit!(() -> Operators.apply_ice!(state, ctx, ice_physics, dt_i,
                      t_array; tanom = tanomv), bench, :ice)
    end

    return released
end

end

