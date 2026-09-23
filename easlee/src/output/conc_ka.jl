
# ── KernelAbstractions concentration kernel ───────────────────────────────────
# conc_parcel_kernel! deposits one parcel's mass per work-item.
#
# Race-condition safety: multiple parcels can fall in the same grid cell.
# The += on grid/massarray uses Atomix.@atomic for correctness on GPU.
# On CPU this compiles to a plain non-atomic add (no lock overhead) because
# the CPU backend runs each work-item sequentially (or via @threads, where
# we pass a unique-cell guarantee via pre-sorting — for now always atomic).
#
# Grid arrays are Float32 (concarray) and Float64 (massarray).
# GPU Float64 atomics require sm_60+ — RTX A1000 is sm_86 → fine.

using KernelAbstractions
using Atomix

@kernel function conc_parcel_kernel!(
    # output grids (atomic write)
    grid      ::AbstractArray{Float32,5},
    massarray ::AbstractArray{Float64,5},
    # parcel state arrays (read-only)
    @Const(Dtime_ns_arr),   # Int64 ms since epoch
    @Const(nest_arr),        # Int8
    @Const(indx_arr),
    @Const(indy_arr),
    @Const(indl_arr),
    @Const(mass_arr),
    @Const(relbox_arr),
    @Const(outclass_arr),   # release box → output column ([output].group);
                            # length-0 = identity (per-box columns)
    # met arrays (read-only)
    @Const(zlev_arr),
    @Const(rho_arr),
    # output grid geometry
    zi        ::AbstractVector{Float32},  # output level interfaces
    Niout     ::Int32,
    Njout     ::Int32,
    Nc        ::Int32,
    # met time index
    indt_n    ::Float32,
    nn        ::Int32,      # nest index this kernel handles
    currentdate_ns::Int64,  # current date in ms since epoch
    # grid config (isbitstype)
    gridconfig,
)
    j = @index(Global, Linear)
    # KA kernels cannot have bare return — nest all guards inside if blocks
    if j <= length(indx_arr) &&
       Dtime_ns_arr[j] == currentdate_ns &&
       Int32(nest_arr[j]) == nn

        indx = indx_arr[j]
        indy = indy_arr[j]
        indl = indl_arr[j]
        mass = mass_arr[j]

        zout ::Float32 = interpolation4Dlevsingle(zlev_arr, indx, indy, indl, indt_n)
        xout ::Float32 = indxtoxout(gridconfig, indx, indy)
        yout ::Float32 = indytoyout(gridconfig, indx, indy)

        if 1.0f0 <= xout <= Float32(Niout) && 1.0f0 <= yout <= Float32(Njout) &&
           zout >= 0.0f0

            l = find_out_layer(zout, zi)
            if l != 0
                rho ::Float32 = interpolation4Dsingle(rho_arr, indx, indy, indl, indt_n)
                if isfinite(rho) && rho > 0.0f0
                    ii  = Base.unsafe_trunc(Int32, xout)
                    jj  = Base.unsafe_trunc(Int32, yout)
                    box = Nc == Int32(1) ? Int32(1) : begin
                        rb = max(Int32(1), Int32(relbox_arr[j]))
                        length(outclass_arr) == 0 ? rb :
                            Int32(outclass_arr[min(rb, Int32(length(outclass_arr)))])
                    end
                    age = Int32(1)

                    addvol  = mass / rho
                    addmass = Float64(mass)

                    # Atomic accumulation — safe for concurrent writes.
                    Atomix.@atomic grid[ii, jj, l, box, age]      += addvol
                    Atomix.@atomic massarray[ii, jj, l, box, age] += addmass
                end
            end
        end
    end
end

# ── CPU/GPU dispatch wrapper ───────────────────────────────────────────────────

function deposit_to_grid_ka!(ctx, state; backend=nothing)
    rc  = ctx.runconfig
    oc  = ctx.output

    grid      = oc.concarray
    massarray = oc.massarray
    Niout     = oc.Niout
    Njout     = oc.Njout
    Nc        = Int32(rc.Nc)
    zi        = oc.outinterfaces

    currentdate    = rc.currentdate
    t_array        = MetAdapter.met_time_index_array(currentdate, ctx.nests)
    currentdate_ns = Int64(Dates.value(currentdate))

    ka_backend = isnothing(backend) ? CPU(; static=true) : backend
    block_size = ka_backend isa CPU ? 256 : 256

    for n in 1:length(ctx.nests)
        nestctx   = ctx.nests[n]
        met       = nestctx.metfields3D
        zlev_arr  = met.zlev
        rho_arr   = met.rho
        gridconfig = nestctx.gridconfig
        indt_n    = t_array[n]

        npar = length(state.indx)

        kernel! = conc_parcel_kernel!(ka_backend, block_size)
        kernel!(
            grid, massarray,
            state.Dtime_ns, state.nest,
            state.indx, state.indy, state.indl, state.mass, state.relbox,
            oc.outclass,
            zlev_arr, rho_arr,
            zi, Niout, Njout, Nc,
            indt_n, Int32(n), currentdate_ns,
            gridconfig;
            ndrange = npar,
        )
        KernelAbstractions.synchronize(ka_backend)
    end

    if _DEBUG_CONC
        n_dep = sum(state.Dtime_ns .== currentdate_ns)
        maxval, maxidx = findmax(view(grid, :, :, :, 1, 1))
        println("DEBUG conc_ka (Niout=", Niout, " Njout=", Njout, ")")
        println("  parcels processed = ", n_dep)
        println("  concarray max     = ", maxval, " at ", maxidx)
    end
    return nothing
end

# ── Explicit-arrays GPU variant ────────────────────────────────────────────────
# Accepts pre-allocated device arrays (CuArray or plain Array) for all inputs.
# Call this when the parcel state and output grids already live on the GPU to
# avoid any CPU→GPU copies per step.
#
# Arguments:
#   backend       — KA backend (CUDABackend(), CPU(), …)
#   grid          — AbstractArray{Float32,5} on device  (output accumulator)
#   massarray     — AbstractArray{Float64,5} on device
#   gs            — NamedTuple of device parcel-state arrays
#   gpu_met_nest  — MetFields3D with device arrays (CuArray on GPU, Array on CPU)
#   zi_dev        — AbstractVector{Float32} on device  (output level interfaces)
#   Niout, Njout, Nc — Int32 output-grid dimensions
#   indt_n        — Float32 met-time interpolation index
#   nn            — Int32 nest index
#   currentdate_ns — Int64 ms since epoch
#   gridconfig    — isbitstype grid config struct
function deposit_conc_gpu!(
    backend,
    grid, massarray,
    gs,
    gpu_met_nest,
    zi_dev,
    Niout     ::Int32,
    Njout     ::Int32,
    Nc        ::Int32,
    indt_n    ::Float32,
    nn        ::Int32,
    currentdate_ns::Int64,
    gridconfig,
)
    zlev_arr  = gpu_met_nest.zlev
    rho_arr   = gpu_met_nest.rho
    npar      = length(gs.indx)

    kernel! = conc_parcel_kernel!(backend, 256)
    # device-resident state: outclass must live on the same device (haskey via
    # hasproperty — legacy gs tuples without it fall back to identity mapping)
    ocl_dev = hasproperty(gs, :outclass) ? gs.outclass : similar(gs.relbox, Int32, 0)
    kernel!(
        grid, massarray,
        gs.Dtime_ns, gs.nest,
        gs.indx, gs.indy, gs.indl, gs.mass, gs.relbox,
        ocl_dev,
        zlev_arr, rho_arr,
        zi_dev, Niout, Njout, Nc,
        indt_n, nn, currentdate_ns,
        gridconfig;
        ndrange = npar,
    )
    KernelAbstractions.synchronize(backend)
    return nothing
end
