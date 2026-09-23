
# ── runtime — Main-scope bridge: legacy met ingestion → typed engine ──────────
# Included by the scripts/CLI (alongside the legacy readers) before
# build_run.jl.  build_runtime wires the legacy met/grid initialisation to the
# EASLEE typed state/context layer; read_met_step! / step_runtime! drive one
# model step.  Internal machinery of the read_runspec(TOML) → build_run /
# run_spec! path — build_run feeds it the input dict it emitted itself, so the
# user never sees this layer.

using Dates

"""
    build_runtime(input; prefetch=true, extra_slots=0, F=Float32, I=Int32, …)
        → NamedTuple (rt)

Build the runtime from a readflexpartinput-style `input` dict: initialise the
legacy met/grid structures, then wire the typed state/context layer on top.
Returns the runtime tuple (ctx, state, met arrays, prefetchers, configs).

Meteorological reading stays in the legacy layer (readmetfiles!).
Physics and orchestration run through EASLEE typed structs.
`extra_slots` reserves parcel slots beyond the release schedule (boundary
inflow, domainfill headroom); `prefetch` enables the async met prefetcher.
"""
function build_runtime(input; gpu::Bool=false,
                                  F::Type{<:AbstractFloat}=Float32,
                                  I::Type{<:Integer}=Int32,
                                  extra_slots::Int=0,
                                  prefetch::Bool=true,
                                  crop=nothing)
    nest = input["nest"]

    metconfig = Vector{Any}(undef, nest)
    for n in 1:nest
        metconfig[n] = gridcheck_onefile!(input[n]["pathmetfields"], input[n])
        # Domain crop (domain-fill runs): shrink the met to box + halo BEFORE
        # the field arrays are allocated.  crop = (lon0, lat0, lon1, lat1).
        crop === nothing ||
            _crop_metconfig!(metconfig[n], crop[1], crop[2], crop[3], crop[4])
    end

    # ── Fix C: publish the per-driver vertical-mixing closure ────────────────
    # gridcheck_onefile! derived (native_e, r_w, c_k) from each nest's variables
    # namelist (metconfig["vmix"], see vmix_namelist.jl); the closure is
    # run-global, so nests must agree.  Hanna-fill drivers get
    # (false, NaN, NaN) — numerically identical to the historic behaviour by
    # construction (see tke_partition in physics/turbulence.jl).
    let vm1 = get(metconfig[1], "vmix", (false, NaN32, NaN32))
        for n in 2:nest
            v = get(metconfig[n], "vmix", (false, NaN32, NaN32))
            same = v[1] == vm1[1] && (!v[1] || (v[2] == vm1[2] && v[3] == vm1[3]))
            same || @info("vertical mixing: nests declare different closures — " *
                          "each nest uses its OWN (NestRunConfig.vmix_*); the " *
                          "global default below is nest 1's",
                          nest_1 = vm1, nest_n = v, n = n)
        end
        EASLEE.set_vertical_mixing!(vm1[1]; r_w=vm1[2], c_k=vm1[3])
    end

    runconfig = initializerun(input, metconfig)
    outputconfig = initoutgrid(runconfig, input, metconfig)

    gridconfig = Vector{Any}(undef, nest)
    metfields3D = Vector{Any}(undef, nest)
    metfields2D = Vector{Any}(undef, nest)
    metfields3Dgpu = Vector{Any}(undef, nest)
    levelconfig = Vector{Any}(undef, nest)

    # outinterfaces == outheight: OUTHEIGHTS values ARE the layer edge heights.
    outputconfig["outinterfaces"] = outputconfig["outheight"]
    outputconfig["massarray"] = zeros(Float64, size(outputconfig["concarray"]))

    for n in 1:nest
        gridconfig[n] = initgridconfig(metconfig[n], outputconfig, n)
        levelconfig[n] = initlevelconfig(metconfig[n], outputconfig, n)
        metfields3D[n], metfields2D[n], metfields3Dgpu[n] = initarrays(metconfig[n], gridconfig[n], levelconfig[n])
    end

    ctx = EASLEE.build_context(input, runconfig, outputconfig,
                                   metconfig, gridconfig, levelconfig,
                                   metfields3D, metfields2D; F=F, I=I)

    # Allocate ParcelState directly from the typed release schedule.
    # All parcels start inactive (Dtime = DateTime(0)); release_parcels! activates them.
    # extra_slots reserves additional unused slots for boundary injection or chemistry inflow.
    maxtraj = sum(Int, ctx.releases.nbpart) + extra_slots
    state = EASLEE.ParcelState(maxtraj; F=F, I=I)
    state.nest_buffers = [Int32[] for _ in 1:length(ctx.nests)]

    # Async met prefetching: one staging copy of the met arrays per nest.
    # Effective only with ≥ 2 Julia threads; nothing to disable otherwise
    # (the prefetcher falls back to synchronous loads on its own).
    prefetchers = prefetch ?
        [MetPrefetcher(metfields3D[n], metfields2D[n]) for n in 1:nest] :
        nothing

    return (
        input = input,
        nest = nest,
        outputconfig = outputconfig,
        metconfig = metconfig,
        gridconfig = gridconfig,
        levelconfig = levelconfig,
        metfields3D = metfields3D,
        metfields2D = metfields2D,
        metfields3Dgpu = metfields3Dgpu,
        prefetchers = prefetchers,
        ctx = ctx,
        state = state,
    )
end

"""
    read_met_step!(rt, currentdate) → nothing

Read (or prefetch-swap) the meteorology of every nest for one model step at
`currentdate`.  Calls the legacy readmetfiles! with typed NestContext; t_index
is written back to the typed NestRunConfig inside readmetfiles! — no separate
sync needed.
"""
function read_met_step!(rt, currentdate::DateTime)
    rc = rt.ctx.runconfig
    rc.currentdate = currentdate
    direction = Int(rc.direction)
    for n in 1:rt.nest
        if rt.prefetchers !== nothing
            readmetfiles_prefetch!(rt.prefetchers[n],
                          rt.ctx.nests[n], currentdate, direction,
                          rt.metconfig[n],
                          rt.metfields3D[n], rt.metfields2D[n], rt.metfields3Dgpu[n])
        else
            readmetfiles!(rt.ctx.nests[n], currentdate, direction,
                          rt.metconfig[n],
                          rt.metfields3D[n], rt.metfields2D[n], rt.metfields3Dgpu[n])
        end
    end
    return nothing
end

"""
    step_runtime!(rt, currentdate; do_output=true)

Advance the model one full synchronisation step through the typed
orchestration (release → advance → deposit; see step_model!).
"""
function step_runtime!(rt, currentdate::DateTime; do_output::Bool=true)
    return EASLEE.step_model!(rt.state, rt.ctx, currentdate; do_output=do_output)
end

"""
Initialize output NetCDF using the same legacy writer.
"""
function initialize_output_file!(rt)
    rc = rt.ctx.runconfig
    rc.dateoutput = rc.startdate
    date = Dates.format(rc.dateoutput + rc.outstep, "yyyymmddHHMMSS")
    filename = rt.input["pathoutput"] * "grid_conc_" * date * ".nc"
    initnetcdf!(rt.ctx.output, filename)
    return filename
end

"""
Reset accumulation arrays between output writes.
"""
function reset_output_accumulators!(rt)
    EASLEE.Types.reset_accumulators!(rt.ctx.output)
    return nothing
end
