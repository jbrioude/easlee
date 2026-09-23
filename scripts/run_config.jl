#!/usr/bin/env julia
# run_config.jl — run any EASLEE TOML config to NetCDF (legacy-free).
#
#   julia -t 1 --project=easlee scripts/run_config.jl configs/arome_hdf_france.toml
#   CFG=configs/era5_conv.toml julia -t 1 --project=easlee scripts/run_config.jl
#
# Wraps read_runspec → run_spec! (build_run + step loop + grid_<name>.nc).
using Dates, Printf
global gpu  = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))

cfg = get(ENV, "CFG", isempty(ARGS) ? error("usage: run_config.jl <config.toml>  (or CFG=...)") : ARGS[1])
isfile(cfg) || error("config not found: $cfg")
spec = read_runspec(cfg)
@info "EASLEE run \"$(spec.name)\": $(spec.direction)  $(spec.start) → $(spec.stop)  Δt=$(spec.timestep_s)s"
EASLEE.Operators.gpu_status_banner()   # one line: GPU used or not (and why not)
prefetch = lowercase(get(ENV, "PREFETCH", "0")) in ("1", "true")
run_spec!(spec; prefetch=prefetch, verbose=true)
@info "done → $(joinpath(spec.output.path, "grid_$(isempty(spec.output.name) ? spec.name : spec.output.name).nc"))"
