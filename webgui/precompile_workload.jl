# Precompile workload for build_sysimage.jl — touches the GUI's hot code paths
# (config parse, schema JSON, JSON round-trip, a CairoMakie render) WITHOUT
# needing any met data, so the sysimage bakes in the methods that otherwise
# compile on first use.  Best-effort: every step is guarded so a missing optional
# dependency never aborts the trace.
try
    const _ROOT = normpath(joinpath(@__DIR__, ".."))

    # ── config layer: schema + TOML round-trip ────────────────────────────────
    try
        include(joinpath(_ROOT, "easlee", "src", "config", "EasleeConfig.jl"))
        @eval using .EasleeConfig
        sj = Base.invokelatest(schema_json)
        for cfg in ("era5_conv.toml", "example_dust.toml")
            p = joinpath(_ROOT, "configs", cfg)
            isfile(p) || continue
            s  = Base.invokelatest(read_runspec, p)
            d  = Base.invokelatest(runspec_to_dict, s)
            Base.invokelatest(dict_to_runspec, d)
        end
    catch e
        @info "precompile: config step skipped ($e)"
    end

    # ── HTTP + JSON3 + TOML serialisation ─────────────────────────────────────
    try
        @eval import JSON3, TOML, HTTP
        io = IOBuffer(); Base.invokelatest(TOML.print, io, Dict("a"=>1, "b"=>[1,2,3]))
        Base.invokelatest(JSON3.write, Dict("x"=>1, "y"=>[1.0,2.0]))
    catch e
        @info "precompile: http/json step skipped ($e)"
    end

    # ── CairoMakie render → PNG bytes (the figure pipeline) ────────────────────
    try
        @eval import CairoMakie
        CairoMakie.activate!()
        fig = CairoMakie.Figure(size=(400, 300))
        ax  = CairoMakie.Axis(fig[1, 1])
        CairoMakie.heatmap!(ax, rand(Float32, 20, 20))
        CairoMakie.scatter!(ax, rand(10), rand(10))
        mktempdir() do dir
            path = joinpath(dir, "p.png")
            CairoMakie.save(path, fig)
            read(path)
        end
    catch e
        @info "precompile: render step skipped ($e)"
    end
catch e
    @info "precompile workload: top-level skip ($e)"
end
