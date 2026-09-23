# validate_captex_run.jl — heavy model path for validate_captex.jl (MODE=run).
# Multi-release driver, forward AND backward.  Included by validate_captex.jl, so
# read_emit/read_meas/read_stations/datem_stats and OBS_DIR/MET_ROOT are in scope.
#
#   DIRECTION=forward  RELEASE=all|N|"1,3,5"  RUNDAYS=4 NP=200000
#       run each release from its source; sample modelled surface conc at every
#       station over each obs window; pool all (obs,model) pairs → ATMES-II stats.
#   DIRECTION=backward MAXBACK=8
#       adjoint check: release backward from the strongest receptors, read the
#       source-receptor sensitivity, convolve with the 7 emissions → modelled conc.

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))
import NetCDF

const DIRECTION = lowercase(get(ENV, "DIRECTION", "forward"))
const RUNDAYS   = parse(Float64, get(ENV, "RUNDAYS", "4"))
const NP        = parse(Int, get(ENV, "NP", "200000"))
const MAXBACK   = parse(Int, get(ENV, "MAXBACK", "8"))
const RELZ      = (parse(Float64, get(ENV, "RELZ0", "0.0")),
                   parse(Float64, get(ENV, "RELZ1", "20.0")))   # release layer [m AGL]
const CTL       = parse(Float64, get(ENV, "CTL", "5.0"))        # adaptive substep: dt ≤ tL/CTL
const IFINE     = parse(Int,    get(ENV, "IFINE", "4"))         # turbulence sub-steps per dt
# Output domain + longitude convention are auto-detected from the met grid
# (CAPTEX: 0–360°, 260–300°E; ETEX: −180…180 across the prime meridian).  obs are
# always −180…180, so convert them to the grid's convention with conv().
const DOMAIN = Ref{NamedTuple}((lon0=260.0, lat0=30.0, nx=160, ny=100, conv360=true))
conv(x) = DOMAIN[].conv360 ? (x < 0 ? x + 360.0 : Float64(x)) : Float64(x)

function _met_year()
    yrs = Int[]
    for d in readdir(MET_ROOT)
        m = match(r"^(\d{4})\d{4}$", d); m === nothing || push!(yrs, parse(Int, m.captures[1]))
    end
    isempty(yrs) ? error("no yyyymmdd dirs in $MET_ROOT") : maximum(yrs)
end

# Symlink every IFS_*.grb from the day-dirs spanning [t0,t1] into one flat dir.
# Start one day EARLY: a day-dir's earliest forecast valid time is 07:00, so an
# early-morning start (e.g. 04:00) needs the previous day's 18:00-base steps.
function _link_met_window(t0, t1)
    dest = mktempdir(prefix = "captex_met_")
    d = Date(min(t0, t1)) - Day(1)
    while d <= Date(max(t0, t1))
        src = joinpath(MET_ROOT, Dates.format(d, "yyyymmdd"))
        isdir(src) && for f in readdir(src)
            endswith(f, ".grb") && symlink(joinpath(src, f), joinpath(dest, f))
        end
        d += Day(1)
    end
    return dest
end

# Tuning knobs (convergence): NP particles, DXOUT output-grid spacing, TS time step.
const TS    = get(ENV, "TS", "15min")
const DXOUT = parse(Float64, get(ENV, "DXOUT", "0.25"))

# Read the met grid's lon/lat extent + longitude convention from one IFS_EA file,
# so the output grid auto-fits the region (CAPTEX or ETEX) with no hard-coding.
function _met_extent()
    f = ""
    for d in sort(readdir(MET_ROOT; join=true))
        isdir(d) || continue
        for g in readdir(d; join=true)
            (occursin("IFS_EA", basename(g)) && endswith(g, ".grb")) && (f = g; break)
        end
        f != "" && break
    end
    f == "" && error("no IFS_EA*.grb under $MET_ROOT")
    # Use the DECODED lon/lat axes (distinctLongitudes/Latitudes), not the header
    # longitudeOfFirstGridPoint key: ecCodes normalises a meridian-crossing grid
    # to -180..180 (ETEX header says 350->30 but decodes to -10..30), and that
    # decoded convention is the one EASLEE's gridconfig (and the parcels) use.
    lonmin = lonmax = latmin = latmax = 0.0
    GribFile(f) do gf
        for msg in gf
            m = Message(msg)
            lons = Float64.(m["distinctLongitudes"]); lats = Float64.(m["distinctLatitudes"])
            lonmin = minimum(lons); lonmax = maximum(lons)
            latmin = minimum(lats); latmax = maximum(lats)
            break
        end
    end
    conv360 = lonmin >= 0.0 && lonmax > 180.0
    nx = max(2, round(Int, (lonmax - lonmin) / DXOUT)); ny = max(2, round(Int, (latmax - latmin) / DXOUT))
    return (lon0=lonmin, lat0=latmin, nx=nx, ny=ny, conv360=conv360)
end

# Common spec; `kind` = :forward (release at source) or :backward (at receptor).
function build_spec(name, lon, lat, rstart, rstop, t0, t1, metdir, outdir, mass_kg, dir)
    Dict{String,Any}(
        "name" => name,
        "run"  => Dict("direction"=>dir, "start"=>min(t0,t1), "stop"=>max(t0,t1), "timestep"=>TS),
        "meteorology" => [Dict("model"=>"ERA5", "path"=>metdir, "variables"=>"ERA5_CAPTEX",
                               "available"=>"auto", "max_gap"=>get(ENV, "MAXGAP", "6h"))],
        "species" => [Dict("id"=>"pmch", "kind"=>"gas", "specnum"=>24,
                           "dry_deposition"=>false, "wet_deposition"=>false)],
        "release" => [Dict("species"=>"pmch", "geometry"=>"box",
                           "lon"=>[lon, lon], "lat"=>[lat, lat],
                           "z"=>[RELZ[1], RELZ[2]], "zkind"=>"agl",
                           "time"=>[rstart, rstop], "particles"=>NP, "mass"=>"$(mass_kg)kg")],
        "physics" => Dict("ctl"=>CTL, "ifine"=>IFINE,
                          # MO=false turns the Monin-Obukhov surface-layer wind off (threads-only
                          # physics; needed for a same-physics threads-vs-ka_gpu benchmark)
                          "surface_layer_wind"=>parse(Bool, get(ENV, "MO", "true")),
                          "turbulence"=>Dict("scheme"=>get(ENV,"TURB","hanna")),
                          "convection"=>Dict("scheme"=>get(ENV,"CONV","massflux")),
                          "vertical"=>Dict("mode"=>"kinematic")),
        "output" => Dict("quantity"=>"concentration", "per_species"=>true, "per_release"=>true,
                         "kernel"=>get(ENV, "KERNEL", "uniform"),
                         "path"=>outdir, "interval"=>"3h", "average"=>"3h", "sample"=>TS,
                         "grid"=>Dict("crs"=>"lonlat", "lon0"=>DOMAIN[].lon0, "lat0"=>DOMAIN[].lat0,
                                      "dx"=>DXOUT, "dy"=>DXOUT, "nx"=>DOMAIN[].nx, "ny"=>DOMAIN[].ny,
                                      "heights"=>[0.0, 100.0])),   # auto-fit to the met grid
    )
end

# Surface field (lon×lat) of the first conc/sens variable in a frame NetCDF.
function _read_surface(nc, want)   # want = "_conc" or "_sens"
    v = NetCDF.open(nc); name = ""
    for k in keys(v.vars)
        (occursin(want, k) && !occursin("_mr", k)) && (name = k; break)
    end
    name == "" && (NetCDF.close(v); error("no $want var in $nc"))
    units = String(get(v.vars[name].atts, "units", "kg m-3"))
    raw = v.vars[name][ntuple(_ -> Colon(), ndims(v.vars[name]))...]
    NetCDF.close(v)
    a = Float32.(raw); nd = ndims(a)
    # height index 1 is the surface-edge MARKER (always empty); the real layers
    # (e.g. 0–100 m) start at index 2 — sum those for the surface field.
    hi = size(a, 3); lays = hi >= 2 ? (2:hi) : (1:hi)
    surf = nd == 3 ? dropdims(sum(a[:, :, lays], dims=3), dims=3) :
                     dropdims(sum(a[:, :, lays, :], dims=(3,4)), dims=(3,4))
    return surf, units
end

# Step a built runtime, writing per-interval frames; return (frames, ftimes, lons, lats, units).
function _run_frames(spec, want)
    rt, obj = build_run(spec; prefetch=false)
    ctx = rt.ctx; st = rt.state; rc = ctx.runconfig; ts = rc.timestep; kw = step_kwargs(obj)
    dir = Int(rc.direction)
    t0  = dir == -1 ? rc.enddate : rc.startdate
    nsteps = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(ts)))
    lstep  = max(1, round(Int, abs(Dates.value(rc.outstep)) / abs(Dates.value(ts))))
    read_met_step!(rt, t0)
    frames = Matrix{Float32}[]; ftimes = DateTime[]; units = "kg m-3"
    lons = Float64[]; lats = Float64[]
    mkpath(spec.output.path)
    acc = 0                                   # samples accumulated in the current interval
    for itime in 1:nsteps
        cur = t0 + itime * ts
        try                                       # stop gracefully when the met window runs out
            read_met_step!(rt, cur)
            EASLEE.step_model!(st, ctx, cur; do_output=true, kw...)   # sample conc EVERY step
        catch e
            @warn "met/step stopped at $cur ($(typeof(e))) — using $(length(frames)) frames so far"
            break
        end
        EASLEE.State.sync_dtime_ns!(st)
        acc += 1
        if itime % lstep == 0                 # interval end → write the average, then reset
            # a fresh file per frame: reusing one name trips NetCDF.jl's handle cache
            fnc = joinpath(spec.output.path, "_frame_$(itime).nc")
            isfile(fnc) && rm(fnc)
            initnetcdf!(ctx.output, fnc); writenetcdf(ctx.output, fnc, acc)
            sf, units = _read_surface(fnc, want); push!(frames, sf); push!(ftimes, cur)
            if isempty(lons)                  # capture grid axes once
                vv = NetCDF.open(fnc); lons = Float64.(vv.vars["longitude"][:]); lats = Float64.(vv.vars["latitude"][:]); NetCDF.close(vv)
            end
            rm(fnc; force=true)
            EASLEE.Types.reset_accumulators!(ctx.output)
            acc = 0
        end
    end
    isempty(lons) && (lons = collect(260.0:0.25:300.0); lats = collect(30.0:0.25:55.0))
    return frames, ftimes, lons, lats, units
end

_scale_pg(units) = occursin("kg", units) ? 1.0f15 : (occursin("ng", units) ? 1.0f3 : 1.0f0)

function _bilin(fr, lons, lats, lon, lat)
    (lon < lons[1] || lon > lons[end] || lat < lats[1] || lat > lats[end]) && return 0.0
    i = clamp(searchsortedlast(lons, lon), 1, length(lons)-1); j = clamp(searchsortedlast(lats, lat), 1, length(lats)-1)
    tx = (lon-lons[i])/(lons[i+1]-lons[i]); ty = (lat-lats[j])/(lats[j+1]-lats[j])
    return (1-tx)*(1-ty)*fr[i,j] + tx*(1-ty)*fr[i+1,j] + (1-tx)*ty*fr[i,j+1] + tx*ty*fr[i+1,j+1]
end

# ── parse RELEASE selection ───────────────────────────────────────────────────
function _release_indices(nrel)
    s = lowercase(get(ENV, "RELEASE", "all"))
    s == "all" && return collect(1:nrel)
    occursin(",", s) && return parse.(Int, split(s, ","))
    return [parse(Int, s)]
end

# ── FORWARD: each release → stations ──────────────────────────────────────────
function run_forward(rels, stns, obs, off)
    O = Float64[]; M = Float64[]; tags = Int[]; S = Int[]
    per = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
    for ri in _release_indices(length(rels))
        rel = rels[ri]
        rstart = rel.start + Year(off); rstop = rel.stop + Year(off)
        runend = rstart + Day(round(Int, RUNDAYS))
        outdir = joinpath(REPO, "results", "captex_fwd_rel$(ri)"); mkpath(outdir)
        metdir = _link_met_window(rstart, runend)
        @printf("[fwd] release %d  %s @ (%.2f,%.2f) %.0f g → %s … %s\n",
                ri, rstart, rel.lat, rel.lon, rel.mass_g, rstart, runend)
        spec = validate_runspec(dict_to_runspec(build_spec(
            "captex_fwd_rel$(ri)", conv(rel.lon), rel.lat, rstart, rstop, rstart, runend,
            metdir, outdir, rel.mass_g/1000, "forward")))
        frames, ftimes, lons, lats, units = _run_frames(spec, "_conc")
        sc = _scale_pg(units)
        oo = Float64[]; mm = Float64[]
        for o in obs
            os = o.start + Year(off); oe = o.stop + Year(off)
            (os >= rstart && oe <= runend) || continue
            idx = findall(t -> os < t <= oe + Minute(1), ftimes); isempty(idx) && continue
            mval = mean(_bilin(frames[k], lons, lats, conv(o.lon), o.lat) for k in idx) * sc
            push!(O, o.conc_pg); push!(M, mval); push!(tags, ri); push!(S, o.stn)
            push!(oo, o.conc_pg); push!(mm, mval)
        end
        per[ri] = (oo, mm)
    end
    return O, M, tags, S, per
end

# Per-station time-integrated dosage (Σ conc over each station's windows) and the
# FLEXPART-11 metric set on it (Bakels et al. 2024, Table 1: FMS, PCC, FA2, FA5,
# NMSE).  This matches FLEXPART's n = #stations pairing, unlike the per-window
# datem_stats (which pairs every station × 3-h window incl. zeros).
function dosage_report(O, M, S)
    od = Dict{Int,Float64}(); md = Dict{Int,Float64}()
    for i in eachindex(S)
        od[S[i]] = get(od, S[i], 0.0) + O[i]
        md[S[i]] = get(md, S[i], 0.0) + M[i]
    end
    ks = sort(collect(keys(od)))
    o = [od[k] for k in ks]; m = [md[k] for k in ks]; n = length(o)
    ob = sum(o)/n; mb = sum(m)/n
    fb = (mb+ob) > 0 ? 2(mb-ob)/(mb+ob) : 0.0
    nmse = (mb*ob) > 0 ? sum((m.-o).^2)/n/(mb*ob) : NaN
    pcc  = (std(o) > 0 && std(m) > 0) ? cor(o, m) : 0.0
    fac(k) = (p=[(o[i]>0&&m[i]>0) for i in 1:n]; any(p) ? count(i->p[i]&&1/k<=m[i]/o[i]<=k,1:n)/count(p) : 0.0)
    both = count(i->o[i]>0&&m[i]>0,1:n); either = count(i->o[i]>0||m[i]>0,1:n)
    fms = either>0 ? both/either : 0.0
    return (; n, fms, pcc, fa2=fac(2.0), fa5=fac(5.0), nmse, fb)
end

# ── BACKWARD: receptor → source-receptor sensitivity → conc ───────────────────
# C(receptor) ≈ Σ_r  SRR(source_r) [s m⁻³] · Q_r/Δt_r [kg s⁻¹]  → pg m⁻³.
# (Seibert & Frank 2004 convention — VERIFY the absolute factor against forward
# on the first real run; the adjoint should reproduce the forward concentrations.)
function run_backward(rels, stns, obs, off)
    cand = sort([o for o in obs if o.conc_pg > 0]; by = o -> -o.conc_pg)[1:min(MAXBACK, count(o->o.conc_pg>0, obs))]
    O = Float64[]; M = Float64[]
    for (n, o) in enumerate(cand)
        tm1 = o.stop + Year(off); tm0 = o.start + Year(off)
        runstart = tm1 - Day(round(Int, RUNDAYS))
        outdir = joinpath(REPO, "results", "captex_bwd_$(o.stn)_$(n)"); mkpath(outdir)
        metdir = _link_met_window(runstart, tm1)
        @printf("[bwd] %d/%d  stn %d @ (%.2f,%.2f)  obs=%.1f  back %s … %s\n",
                n, length(cand), o.stn, o.lat, o.lon, o.conc_pg, runstart, tm1)
        spec = validate_runspec(dict_to_runspec(build_spec(
            "captex_bwd_$(o.stn)_$(n)", conv(o.lon), o.lat, tm0, tm1, runstart, tm1,
            metdir, outdir, 1.0, "backward")))
        frames, ftimes, lons, lats, units = _run_frames(spec, "_sens")   # s m⁻³ (time-integrated)
        srr = isempty(frames) ? zeros(Float32, length(lons), length(lats)) : sum(frames)  # ∫ over run
        mod_pg = 0.0
        for r in rels
            rs = r.start + Year(off); re = r.stop + Year(off)
            (rs >= runstart && re <= tm1) || continue
            qrate = (r.mass_g/1000) / max(Dates.value(Second(re - rs)), 1)   # kg/s
            mod_pg += _bilin(srr, lons, lats, conv(r.lon), r.lat) * qrate * 1.0f15
        end
        push!(O, o.conc_pg); push!(M, mod_pg)
    end
    return O, M
end

# ── driver ────────────────────────────────────────────────────────────────────
function drive()
    rels = read_emit(joinpath(OBS_DIR, "emit-t1.txt"))
    stns = read_stations(joinpath(OBS_DIR, "stations.txt"))
    obs  = read_meas(joinpath(OBS_DIR, "meas-t1.txt"))
    off  = _met_year() - year(rels[1].start)
    off != 0 && @warn "YEAR MISMATCH: obs $(year(rels[1].start)) vs met $(_met_year()) — shifting +$(off) yr to run; NOT a valid validation until $(year(rels[1].start)) ERA5 is supplied."
    DOMAIN[] = _met_extent()
    @info "met domain" lon0=DOMAIN[].lon0 lat0=DOMAIN[].lat0 nx=DOMAIN[].nx ny=DOMAIN[].ny conv360=DOMAIN[].conv360 DXOUT=DXOUT
    @info "run knobs" NP TS CTL IFINE RELZ TURB=get(ENV,"TURB","hanna") CONV=get(ENV,"CONV","massflux") MO=get(ENV,"MO","true") ADVANCE=get(ENV,"EASLEE_ADVANCE","threads")

    if DIRECTION == "backward"
        O, M = run_backward(rels, stns, obs, off)
        isempty(O) && (println("no backward pairs"); return)
        s = datem_stats(O, M; thr = 0.0)
        @printf("\n══ CAPTEX BACKWARD (adjoint) — %d receptors ══\n", s.n)
        @printf("  R=%.3f FB=%+.3f FAC2=%.3f FMS=%.1f%% KSP=%.1f RANK=%.3f\n", s.r, s.fb, s.fac2, s.fms, s.ksp, s.rank)
        @printf("  (verify the absolute backward scaling against the forward run)\n")
        return
    end

    O, M, tags, S, per = run_forward(rels, stns, obs, off)
    isempty(O) && (println("no forward pairs"); return)
    outdir = joinpath(REPO, "results"); mkpath(outdir)
    open(joinpath(outdir, "captex_paired.txt"), "w") do io
        println(io, "# release  stn  obs_pg/m3  model_pg/m3")
        for i in eachindex(O); @printf(io, "%d  %d  %.4f  %.4f\n", tags[i], S[i], O[i], M[i]); end
    end
    @printf("\n══ CAPTEX FORWARD ══   releases: %s\n", sort(unique(tags)))
    @printf("%-8s %6s %7s %7s %6s %6s %6s %6s\n", "release", "N", "R", "FB", "FAC2", "FMS", "KSP", "RANK")
    for ri in sort(collect(keys(per)))
        oo, mm = per[ri]; isempty(oo) && continue
        s = datem_stats(oo, mm)
        @printf("%-8d %6d %7.3f %+7.3f %6.2f %6.1f %6.1f %6.2f\n", ri, s.n, s.r, s.fb, s.fac2, s.fms, s.ksp, s.rank)
    end
    s = datem_stats(O, M)
    @printf("%-8s %6d %7.3f %+7.3f %6.2f %6.1f %6.1f %6.2f\n", "POOLED", s.n, s.r, s.fb, s.fac2, s.fms, s.ksp, s.rank)

    # FLEXPART-11-style per-station integrated dosage (Bakels et al. 2024, Table 1)
    d = dosage_report(O, M, S)
    @printf("\n══ per-station integrated dosage (FLEXPART-11 metric set) ══\n")
    @printf("  n=%d  FMS=%.2f  PCC=%.2f  FA2=%.2f  FA5=%.2f  NMSE=%.2f  FB=%+.2f\n",
            d.n, d.fms, d.pcc, d.fa2, d.fa5, d.nmse, d.fb)
    @printf("  FLEXPART 11 (ETEX-1):   FMS=0.68  PCC=0.84  FA2=0.74  FA5=0.87  NMSE=0.43  (n=168)\n")
    @printf("→ %s\n", joinpath(outdir, "captex_paired.txt"))
end

drive()
