#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# ka_compare_grids.jl — compare two EASLEE grid outputs: total mass and Pearson
# correlation.  Used for the backend validations (threads vs ka_gpu) recorded in
# docs/BENCHMARKS.md.
#
#     julia --project=easlee scripts/ka_compare_grids.jl <dir A> <dir B>
#
# Each directory is searched for grid_*.nc and the last one by name is read.
# The comparison uses the last index along the file's final dimension of `conc`,
# which is the last release box for a multi-box run and the top output level for
# a single-box one.  That is what the Python version this replaces did, so the
# numbers in BENCHMARKS.md stay reproducible; pass --all to compare the whole
# array instead.
# ─────────────────────────────────────────────────────────────────────────────

using NetCDF
using Printf
using Statistics

ENV["HDF5_USE_FILE_LOCKING"] = "FALSE"

"Newest grid_*.nc of a directory, and its `conc` with the fill value zeroed."
function load(dir::AbstractString; all::Bool=false)
    isdir(dir) || error("not a directory: $dir")
    files = sort(filter(f -> startswith(f, "grid_") && endswith(f, ".nc"), readdir(dir)))
    isempty(files) && error("no grid_*.nc in $dir")
    path = joinpath(dir, files[end])

    a = Array{Float64}(ncread(path, "conc"))
    fill = try
        Float64(ncgetatt(path, "conc", "_FillValue"))
    catch
        nothing
    end
    fill === nothing || replace!(x -> x == fill ? 0.0 : x, a)
    replace!(x -> isfinite(x) ? x : 0.0, a)

    slice = all ? a : selectdim(a, ndims(a), size(a, ndims(a)))
    return path, vec(collect(slice))
end

function main(args=ARGS)
    opts = filter(a -> startswith(a, "--"), args)
    dirs = filter(a -> !startswith(a, "--"), args)
    if length(dirs) != 2 || "--help" in opts
        println("usage: julia --project=easlee scripts/ka_compare_grids.jl [--all] <dir A> <dir B>")
        return length(dirs) == 2 ? 0 : 1
    end
    all = "--all" in opts

    fa, a = load(dirs[1]; all=all)
    fb, b = load(dirs[2]; all=all)
    length(a) == length(b) ||
        error("grids differ in size: $(length(a)) vs $(length(b)) values")

    ma, mb = sum(a), sum(b)
    r = cor(a, b)
    @printf("%s vs %s: mass %.6e vs %.6e (%.2e apart), correlation %.6f\n",
            basename(fa), basename(fb), ma, mb, abs(mb/ma - 1), r)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
