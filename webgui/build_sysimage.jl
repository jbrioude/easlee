#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# build_sysimage.jl — bake a fast-start sysimage for the EASLEE web GUI.
#
# PackageCompiler cannot cross-compile (a Linux Julia builds a Linux image, a
# Windows Julia builds a Windows .dll), so run this with the SAME Julia you use
# to launch the GUI.  On Windows just double-click  build_sysimage.bat .
#
#   julia --project=webgui webgui/build_sysimage.jl
#
# Produces  easlee_sysimage.<so|dll|dylib>  in the project root.  run_gui.bat /
# run_gui.sh pick it up automatically and the GUI then starts in a few seconds
# instead of a minute, because HTTP, CairoMakie, NetCDF, GRIB and the EASLEE
# kernels are already compiled into the image.
# ─────────────────────────────────────────────────────────────────────────────

import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "webgui"))
try
    @eval import PackageCompiler
catch
    @info "Installing PackageCompiler…"
    Pkg.add("PackageCompiler")
    @eval import PackageCompiler
end

const ROOT   = normpath(joinpath(@__DIR__, ".."))
const EXT    = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
const OUTIMG = joinpath(ROOT, "easlee_sysimage.$EXT")

# Packages to bake into the image (everything the GUI loads at run time).
const PKGS = [:HTTP, :JSON3, :TOML, :CairoMakie, :FileIO, :NetCDF, :StructArrays]

# A tiny workload so the tracer captures the hot methods (config parse, schema,
# figure rendering).  It must not need real met data — it only touches code.
const PRECOMPILE = joinpath(@__DIR__, "precompile_workload.jl")

@info "Building EASLEE sysimage → $OUTIMG  (this takes 10–25 min, once)"
PackageCompiler.create_sysimage(
    PKGS;
    sysimage_path        = OUTIMG,
    precompile_execution_file = PRECOMPILE,
    cpu_target           = PackageCompiler.default_app_cpu_target(),
)
@info "Done.  Launch with:  julia --sysimage \"$OUTIMG\" --project=webgui webgui/easlee_app.jl"
