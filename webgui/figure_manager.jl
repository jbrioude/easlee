using CairoMakie, NetCDF, Statistics
import JSON3

CairoMakie.activate!()

# ── Coastlines (same polylines the Setup map uses) ────────────────────────────
const _COAST = Ref{Vector{Tuple{Vector{Float32},Vector{Float32}}}}(Tuple{Vector{Float32},Vector{Float32}}[])
function _coastlines()
    if isempty(_COAST[])
        try
            d = JSON3.read(read(joinpath(@__DIR__, "frontend", "coastline.json"), String))
            for ln in d["lines"]
                lons = Float32[p[1] for p in ln]; lats = Float32[p[2] for p in ln]
                push!(_COAST[], (lons, lats))
            end
        catch e
            @warn "coastline load failed: $(typeof(e))"
        end
    end
    return _COAST[]
end

# Overlay coastlines clipped to the axis lon/lat window.  Preferred source: the
# 0.5 contour of the met LAND-SEA MASK of the live runtime (exact for the run's
# own domain, no external data); fallback: the bundled coastline.json segments.
function _draw_coast!(ax, lonmin, lonmax, latmin, latmax)
    drew = false
    try
        if LIVE.rt !== nothing && haskey(LIVE.rt.metfields2D[1], "lsm")
            lsm = LIVE.rt.metfields2D[1]["lsm"]
            gc  = LIVE.rt.ctx.nests[1].gridconfig
            Ni  = min(Int(gc.Ni), size(lsm, 1)); Nj = min(Int(gc.Nj), size(lsm, 2))
            lons = Float32[EASLEE.Operators.indyindxtolatlon(gc, 1.0, Float64(i))[2] for i in 1:Ni]
            lats = Float32[EASLEE.Operators.indyindxtolatlon(gc, Float64(j), 1.0)[1] for j in 1:Nj]
            if issorted(lons) && issorted(lats)
                contour!(ax, lons, lats, Float32.(lsm[1:Ni, 1:Nj, 1]);
                         levels=Float32[0.5], color=(:gray70, 0.9), linewidth=0.8)
                drew = true
            end
        end
    catch e
        @warn "coast contour failed: $(typeof(e))" maxlog=3
    end
    if !drew
        for (lons, lats) in _coastlines()
            (maximum(lats) < latmin || minimum(lats) > latmax) && continue
            (maximum(lons) < lonmin || minimum(lons) > lonmax) && continue
            lines!(ax, lons, lats; color=(:gray60, 0.5), linewidth=0.7)
        end
    end
    xlims!(ax, lonmin, lonmax); ylims!(ax, latmin, latmax)
end

# Log-scale color mapping for plume fields: returns (plotdata, colorrange, label)
# — log10 of the positive values, NaN elsewhere (rendered transparent).
function _log10field(a::AbstractMatrix, unitlabel::String)
    pos = filter(x -> isfinite(x) && x > 0, vec(a))
    isempty(pos) && return (Float32.(a), (0.0f0, 1.0f0), unitlabel)
    hi = log10(maximum(pos))
    lo = max(log10(minimum(pos)), hi - 6.0)     # ≤ 6 decades of dynamic range
    la = map(x -> (isfinite(x) && x > 0) ? Float32(log10(x)) : NaN32, a)
    return (la, (Float32(lo), Float32(hi)), "log₁₀ $unitlabel")
end

# ── Per-view data cache (updated after each NC write) ─────────────────────────
const FIG_CACHE = Dict{Symbol,Any}(
    :column    => nothing,
    :cross     => nothing,
    :particles => nothing,
    :met       => nothing,
)
const FIG_LOCK = ReentrantLock()

# ── Data loaders (called from model_runner after each output step) ─────────────
# `conc` is (lon,lat,height) for a single release box, or (lon,lat,height,relbox)
# when per_release keeps several boxes.  Read the whole array and sum over any
# trailing release dimension so the column/cross-section figures work for both.
# EASLEE writes two output conventions: transport grids use dims longitude/latitude/
# height + a "conc" var; other grids may use lon/lat/plev(+time) + named vars
# (O3, …).  These helpers read EITHER so the figures work for any output.
const _COORDNAMES = Set(["lon","lat","longitude","latitude","height","plev","level","time"])
_findvar(v, cands) = (for c in cands; haskey(v.vars, c) && return c; end; nothing)
function _datavar(v)
    haskey(v.vars, "spec001_conc") && return "spec001_conc"  # physical kg/m³
    haskey(v.vars, "conc") && return "conc"        # transport grid (framed: kg/m³)
    for k in sort(collect(keys(v.vars)))            # else first real ≥3-D species field
        (k in _COORDNAMES) && continue
        (endswith(k, "_pass") || endswith(k, "_cams") || endswith(k, "_ref") || endswith(k, "_mr")) && continue
        ndims(v.vars[k]) >= 3 && return k
    end
    return nothing
end
# read the data var as (lon, lat, [vertical], [time]) — NetCDF.jl reverses dim order
_readvar(v, name) = Float32.(v.vars[name][ntuple(_->Colon(), ndims(v.vars[name]))...])
# NaN → 0 (NaN fills would poison sums and the Colorbar)
_nan0(a) = map(x -> isnan(x) ? zero(x) : x, a)
# collapse to (lon,lat): take the last time frame, then sum over the vertical axis
function _column(raw)
    a = ndims(raw) >= 4 ? raw[:, :, :, end] : raw
    a = _nan0(a)
    return ndims(a) >= 3 ? dropdims(sum(a, dims=3), dims=3) : a
end

# Plottable fields in an output file (≥3-D non-coordinate vars), for the GUI dropdown.
function list_species_in_nc(nc_path::String)
    isfile(nc_path) || return String[]
    try
        v = NetCDF.open(nc_path)
        ks = [k for k in keys(v.vars) if !(k in _COORDNAMES) && ndims(v.vars[k]) >= 3]
        NetCDF.close(v)
        pref = ["spec001_conc","conc","O3"]        # show the usual defaults first
        return vcat([k for k in pref if k in ks], sort([k for k in ks if !(k in pref)]))
    catch
        return String[]
    end
end

# units attribute of a variable ("a.u." when absent)
_var_units(v, name) = String(get(v.vars[name].atts, "units", "a.u."))

function load_column_from_nc!(nc_path::String; species=nothing)
    isfile(nc_path) || return
    try
        v = NetCDF.open(nc_path)
        lonn = _findvar(v, ("longitude","lon")); latn = _findvar(v, ("latitude","lat"))
        dv   = (species !== nothing && haskey(v.vars, String(species))) ? String(species) : _datavar(v)
        (lonn === nothing || latn === nothing || dv === nothing) &&
            (NetCDF.close(v); @warn "column: no lon/lat/data var in $(basename(nc_path))"; return)
        lon = Float32.(v.vars[lonn][:]); lat = Float32.(v.vars[latn][:])
        units = _var_units(v, dv)
        col = _column(_readvar(v, dv))
        NetCDF.close(v)
        lock(FIG_LOCK) do; FIG_CACHE[:column] = (lon, lat, col, dv, units); end
    catch e
        @warn "column NC load error: $(typeof(e))"
    end
end

function load_crosssection_from_nc!(nc_path::String, j_frac::Float64=0.5; species=nothing)
    isfile(nc_path) || return
    try
        v = NetCDF.open(nc_path)
        lonn = _findvar(v, ("longitude","lon")); altn = _findvar(v, ("height","plev","level"))
        dv   = (species !== nothing && haskey(v.vars, String(species))) ? String(species) : _datavar(v)
        (lonn === nothing || altn === nothing || dv === nothing) && (NetCDF.close(v); return)
        lon = Float32.(v.vars[lonn][:]); alt = Float32.(v.vars[altn][:])
        units = _var_units(v, dv)
        latn = _findvar(v, ("latitude","lat"))
        latv = latn === nothing ? Float32[] : Float32.(v.vars[latn][:])
        raw = _readvar(v, dv)                          # (lon,lat,vert[,time])
        raw = ndims(raw) >= 4 ? raw[:, :, :, end] : raw
        Nj  = size(raw, 2); j0 = clamp(round(Int, j_frac * Nj), 1, Nj)
        lat0 = isempty(latv) ? NaN32 : latv[clamp(j0, 1, length(latv))]
        xsec = _nan0(Float32.(raw[:, j0, :]))          # (lon, vert)
        NetCDF.close(v)
        lock(FIG_LOCK) do; FIG_CACHE[:cross] = (lon, alt, xsec, dv, units, lat0); end
    catch e
        @warn "cross-section NC load error: $(typeof(e))"
    end
end

function load_particles_from_live!()
    LIVE.rt === nothing && return
    try
        st   = LIVE.rt.state
        gc   = LIVE.rt.ctx.nests[1].gridconfig
        idx  = findall(i -> st.status[i] == 0x01 && st.nest[i] == Int8(1),
                       eachindex(st.status))
        n    = length(idx)
        lons = Vector{Float32}(undef, n); lats = Vector{Float32}(undef, n)
        levs = Vector{Float32}(undef, n)
        @inbounds for (k, i) in enumerate(idx)
            la, lo = EASLEE.Operators.indyindxtolatlon(gc, Float64(st.indy[i]), Float64(st.indx[i]))
            lons[k] = Float32(lo); lats[k] = Float32(la); levs[k] = st.indl[i]
        end
        lock(FIG_LOCK) do; FIG_CACHE[:particles] = (lons, lats, levs); end
    catch e
        @warn "particles load error: $(typeof(e))"
    end
end

function load_met_from_live!(field::String="Ux", lev::Int=1)
    LIVE.rt === nothing && return
    try
        ctx     = LIVE.rt.ctx
        nestctx = ctx.nests[1]
        hasproperty(nestctx.metfields3D, Symbol(field)) ||
            (@warn "met view: no field \"$field\""; return)
        arr = getfield(nestctx.metfields3D, Symbol(field))
        k   = clamp(lev, 1, size(arr, 3))
        # 2-slot rotating time axis: the CURRENT slot comes from the fractional
        # met time index, NOT nestrunconfig.t_index (those are met-FILE indices).
        it  = try
            clamp(round(Int, EASLEE.MetAdapter.met_time_index_array(
                  ctx.runconfig.currentdate, ctx.nests)[1]), 1, size(arr, 4))
        catch
            1
        end
        gc   = nestctx.gridconfig
        Ni   = min(Int(gc.Ni), size(arr, 1)); Nj = min(Int(gc.Nj), size(arr, 2))
        lons = Vector{Float32}(undef, Ni); lats = Vector{Float32}(undef, Nj)
        for i in 1:Ni
            la, lo = EASLEE.Operators.indyindxtolatlon(gc, 1.0, Float64(i)); lons[i] = Float32(lo)
        end
        for j in 1:Nj
            la, lo = EASLEE.Operators.indyindxtolatlon(gc, Float64(j), 1.0); lats[j] = Float32(la)
        end
        sl = Float32.(arr[1:Ni, 1:Nj, k, it])
        lock(FIG_LOCK) do; FIG_CACHE[:met] = (lons, lats, sl, field, k); end
    catch e
        @warn "met field load error: $(typeof(e))"
    end
end

# ── Render figure → PNG bytes ─────────────────────────────────────────────────
const BG  = RGBf(0.07, 0.07, 0.10)
const AXC = RGBf(0.10, 0.10, 0.14)

function _styled_axis(fig, pos; xlabel="", ylabel="", title="")
    Axis(fig[pos...];
         backgroundcolor = AXC,
         xticklabelcolor = :gray75, yticklabelcolor = :gray75,
         xlabelcolor = :gray75,     ylabelcolor = :gray75,
         titlecolor  = :gray90,
         xlabel = xlabel, ylabel = ylabel, title = title)
end

function render_figure(view::Symbol;
                       met_field::String = "Ux",
                       level::Int        = 1,
                       lat_frac::Float64 = 0.5) :: Vector{UInt8}
    data = lock(FIG_LOCK) do; FIG_CACHE[view]; end

    fig = Figure(size=(820, 520), backgroundcolor=BG)

    if data === nothing
        ax = _styled_axis(fig, [1,1]; title="$view — run the model first")
        text!(ax, 0.5, 0.5; text="No data yet", color=:gray55,
              align=(:center,:center), fontsize=16)
        ax.limits[] = (0,1,0,1)
    elseif view == :column
        lon, lat, col, dv, units = data
        ax = _styled_axis(fig, [1,1];
            xlabel="Longitude [°]", ylabel="Latitude [°]",
            title="$dv — column (Σ layers), log scale")
        la, crange, clabel = _log10field(col, units)
        hm = heatmap!(ax, lon, lat, la; colormap=:inferno, interpolate=false,
                      colorrange=crange, nan_color=:transparent)
        _draw_coast!(ax, minimum(lon), maximum(lon), minimum(lat), maximum(lat))
        Colorbar(fig[1,2], hm; label=clabel,
                 labelcolor=:gray75, ticklabelcolor=:gray75, width=14)
    elseif view == :cross
        lon, alt, xsec, dv, units, lat0 = data
        ax = _styled_axis(fig, [1,1];
            xlabel="Longitude [°]", ylabel="Altitude",
            title=isfinite(lat0) ? "$dv — cross-section at $(round(lat0, digits=1))°N, log scale" :
                                   "$dv — vertical cross-section, log scale")
        la, crange, clabel = _log10field(xsec, units)
        hm = heatmap!(ax, lon, alt, la; colormap=:inferno, interpolate=false,
                      colorrange=crange, nan_color=:transparent)
        Colorbar(fig[1,2], hm; label=clabel,
                 labelcolor=:gray75, ticklabelcolor=:gray75, width=14)
    elseif view == :particles
        x, y, cv = data
        ax = _styled_axis(fig, [1,1];
            xlabel="Longitude [°]", ylabel="Latitude [°]",
            title="Particle positions  (n=$(length(x)))")
        if !isempty(x)
            # pad the view around the cloud so the plume has context
            lo1, lo2 = extrema(x); la1, la2 = extrema(y)
            padx = max(0.5f0, 0.15f0*(lo2-lo1)); pady = max(0.5f0, 0.15f0*(la2-la1))
            sc = scatter!(ax, x, y; color=cv, colormap=:viridis, markersize=3)
            _draw_coast!(ax, lo1-padx, lo2+padx, la1-pady, la2+pady)
            Colorbar(fig[1,2], sc; label="model level",
                     labelcolor=:gray75, ticklabelcolor=:gray75, width=14)
        end
    elseif view == :met
        lon, lat, mdata, fname, flev = data
        ax = _styled_axis(fig, [1,1];
            xlabel="Longitude [°]", ylabel="Latitude [°]",
            title="Met field: $fname  level $flev")
        hm = heatmap!(ax, lon, lat, mdata; colormap=:RdBu, interpolate=false)
        _draw_coast!(ax, minimum(lon), maximum(lon), minimum(lat), maximum(lat))
        Colorbar(fig[1,2], hm;
                 labelcolor=:gray75, ticklabelcolor=:gray75, width=14)
    end

    mktempdir() do dir
        path = joinpath(dir, "fig.png")
        save(path, fig)
        read(path)
    end
end
