#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# EASLEE Web GUI  —  easlee_app.jl
#
# Launch (from the repository root):
#   julia --threads=auto --project=webgui webgui/easlee_app.jl
# ─────────────────────────────────────────────────────────────────────────────

const WEBGUI_DIR = @__DIR__
const EASLEE_DIR = normpath(joinpath(WEBGUI_DIR, "..", "easlee"))

# Legacy globals referenced by initoutgrid and other initialize.jl routines
gpu  = false

# ── EASLEE source ──────────────────────────────────────────────────────────────
include(joinpath(EASLEE_DIR, "src/io/initialize.jl"))
include(joinpath(EASLEE_DIR, "src/io/tools.jl"))
include(joinpath(EASLEE_DIR, "src/io/output.jl"))
include(joinpath(EASLEE_DIR, "EASLEE.jl"))
include(joinpath(EASLEE_DIR, "runtime.jl"))
using .EASLEE, .EASLEE.State, .EASLEE.MetAdapter, .EASLEE.Operators

# ── New TOML/RunSpec configuration layer + data-driven build_run ───────────────
include(joinpath(EASLEE_DIR, "src/config/EasleeConfig.jl"))
using .EasleeConfig
include(joinpath(EASLEE_DIR, "build_run.jl"))

# ── WebGUI modules ─────────────────────────────────────────────────────────────
include(joinpath(WEBGUI_DIR, "model_runner.jl"))
include(joinpath(WEBGUI_DIR, "figure_manager.jl"))
include(joinpath(WEBGUI_DIR, "repl_session.jl"))
include(joinpath(WEBGUI_DIR, "server.jl"))

using Printf, Dates

const PORT = 8080

function open_browser(url::String)
    sleep(1.0)   # give HTTP.serve a moment to bind
    if Sys.islinux()
        run(`xdg-open $url`, wait=false)
    elseif Sys.isapple()
        run(`open $url`, wait=false)
    elseif Sys.iswindows()
        run(`cmd /c start "" "$url"`, wait=false)
    end
end

# ── Auto-load a TOML config from a CLI arg (legacy-free) ───────────────────────
if !isempty(ARGS) && isfile(ARGS[1])
    try
        SPEC_DICT[] = runspec_to_dict(read_runspec(ARGS[1]))
        @info "Auto-loaded config: $(ARGS[1])"
    catch e
        @warn "Auto-load failed: $e"
    end
end

@info "EASLEE Web GUI → http://localhost:$PORT"
Threads.@spawn open_browser("http://localhost:$PORT")
start_server(PORT)   # blocks
