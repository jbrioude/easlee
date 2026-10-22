# output/trajectory.jl
# ─────────────────────────────────────────────────────────────────────────────
# Per-trajectory CSV dump: one row per sampled parcel per output interval, with
# position + time, mass + age and met at the parcel.  Driven by
# [output.trajectory] in the TOML.  A fixed random subsample of parcel SLOTS is
# chosen once (so each row series follows the SAME physical parcel), capped at
# max_parcels to keep the file manageable at 1e6 parcels.
#
# Two row layouts ([output.trajectory].format):
#   "full"    — headered CSV: time,id,lon,lat,height_m,nest,mass_kg,age_s,rho,pblh
#   "demovor" — the legacy DEMOVOR/AROME column set, no header, %.5g reals:
#               yyyy-mm-dd,HH:MM:SS,id,lon,lat,height_m,tke,pblh_m
#               (pblh_m appended 2026-09-09 for the PBL / free-troposphere
#               partition; the seven DEMOVOR fields ahead of it are unchanged)
# ─────────────────────────────────────────────────────────────────────────────
using Random: MersenneTwister, randperm
using Printf: @sprintf

mutable struct TrajWriter
    io       :: IO
    sel      :: Vector{Int}      # sampled parcel slots (fixed across the run)
    every    :: Int              # write every `every` model steps
    lon0     :: Float64          # lon wrap origin: rows report lon ∈ [lon0, lon0+360)
                                 # (the output grid's convention, so a dateline-
                                 # crossing domain stays continuous)
    demovor  :: Bool             # legacy DEMOVOR row layout (see header)
end

"""
    traj_open(spec, state, dt_s) -> TrajWriter | nothing

Open `<output.path>/traj_<name>.csv`, write the header (none for the demovor
layout), and pick the subsample.
"""
function traj_open(spec, state, dt_s::Float64)
    tj = spec.output.trajectory
    tj.enabled || return nothing
    outdir = spec.output.path; mkpath(outdir)
    nm   = isempty(spec.output.name) ? spec.name : spec.output.name
    path = joinpath(outdir, "traj_$(nm).csv")
    io   = open(path, "w")

    intv  = tj.interval_s > 0 ? tj.interval_s : spec.output.sample_s
    every = max(1, round(Int, intv / dt_s))

    # sample from the parcels that are active now (domain-fill places them before the
    # loop); fall back to all slots for a release run where none are active yet.
    rng    = MersenneTwister(20240520)
    active = findall(==(0x01), state.status)
    pool   = isempty(active) ? collect(1:length(state.status)) : active
    sel    = (tj.max_parcels > 0 && tj.max_parcels < length(pool)) ?
             sort!(pool[randperm(rng, length(pool))[1:tj.max_parcels]]) : copy(pool)

    demovor = tj.format == "demovor"
    if !demovor
        hdr = String["time","id","lon","lat","height_m","nest","mass_kg","age_s",
                     "rho_kgm3","pblh_m"]
        println(io, join(hdr, ","))
    end
    @info "trajectory CSV → $path" parcels=length(sel) every_steps=every format=tj.format
    return TrajWriter(io, sel, every, spec.output.grid.lon0, demovor)
end

# Write one block of rows (all sampled active parcels) at `currentdate`.
# Call every step; it self-gates on `every`.  `indt` may be a scalar (single
# met source) or the per-nest array from `met_time_index_array`.
function traj_write!(tw, state, ctx, currentdate, indt, istep::Int)
    tw === nothing && return
    (istep % tw.every == 0) || return
    # Per-nest views, selected below by the parcel's OWN nest.  This used to read
    # ctx.nests[1] for every parcel: in a nested run (e.g. AROME inside ERA5) a
    # parcel that has left the inner domain carries indices on the OUTER grid,
    # and putting those through the inner grid's projection reports it somewhere
    # else entirely — silently, since the bogus indices usually still fall inside
    # the inner grid's bounds.
    nn    = length(ctx.nests)
    gcs   = [ctx.nests[n].gridconfig       for n in 1:nn]
    lcs   = [ctx.nests[n].levelconfig      for n in 1:nn]
    zlevs = [ctx.nests[n].metfields3D.zlev for n in 1:nn]
    rhoFs = [ctx.nests[n].metfields3D.rho  for n in 1:nn]
    pblFs = [ctx.nests[n].metfields2D.pblh for n in 1:nn]
    tkeFs = [ctx.nests[n].metfields3D.tke  for n in 1:nn]
    its   = indt isa AbstractVector ? Float32.(collect(indt)) :
                                      fill(Float32(indt), nn)
    tstr = tw.demovor ? Dates.format(currentdate, "yyyy-mm-dd,HH:MM:SS") :
                        Dates.format(currentdate, "yyyy-mm-ddTHH:MM:SS")
    io   = tw.io
    @inbounds for i in tw.sel
        state.status[i] == 0x01 || continue   # active parcels only
        nst = Int(state.nest[i])
        (1 <= nst <= nn) || continue
        gc = gcs[nst]; lc = lcs[nst]; zlev = zlevs[nst]
        rhoF = rhoFs[nst]; pblhF = pblFs[nst]; tkeF = tkeFs[nst]
        it = its[nst]
        Ni = Float32(gc.Ni); Nj = Float32(gc.Nj)
        ix = state.indx[i]; iy = state.indy[i]; il = state.indl[i]
        (ix < 1f0 || ix >= Ni || iy < 1f0 || iy >= Nj) && continue   # OOB → skip
        # Met index → geographic, via the grid's exact inverse.  (The old
        # OUTGRID-cell route `lat0 + (indytoyout−1)·dy_out` carried the deposit
        # convention's +one-MET-spacing offset, biasing positions by
        # (d_met − d_out) whenever the output spacing differed from the met's —
        # e.g. +0.25° lat on 0.5° met with a 0.25° grid, reporting lat > 90.)
        # Wrap lon into [lon0, lon0+360) so a dateline-crossing domain (e.g.
        # Pacific 150–230) stays continuous instead of splitting to ±180.
        lat, lon_geo = indyindxtolatlon(gc, iy, ix)
        lon = tw.lon0 + mod(lon_geo - tw.lon0, 360.0)
        h   = indltozlev(lc, zlev, ix, iy, il, it)
        if tw.demovor
            # height_m is AGL, so height_m < pblh_m is the boundary-layer test
            # the DEMOVOR PBL / free-troposphere fractions use.
            # tke is a half-level (class H) field: plain read at indl
            tke_p  = Operators.interpolation4Dlevsingle(tkeF, ix, iy, il, it)
            pblh_p = interpolate2Dsingle(pblhF, ix, iy, it)
            println(io, tstr, ',', state.id[i], ',',
                    @sprintf("%.5g,%.5g,%.5g,%.5g,%.5g", lon, lat, h, tke_p, pblh_p))
            continue
        end
        age = Dates.value(state.Dtime[i] - state.Dtimeinit[i]) / 1000.0   # s
        print(io, tstr, ',', state.id[i], ',', lon, ',', lat, ',', h, ',',
              state.nest[i], ',', state.mass[i], ',', age)
        rho_p  = Operators.interpolation4Dsingle(rhoF, ix, iy, il, it)   # class F: stacked
        pblh_p = interpolate2Dsingle(pblhF, ix, iy, it)
        println(io, ',', rho_p, ',', pblh_p)
    end
    flush(io)
    return
end

traj_close(tw) = tw === nothing ? nothing : close(tw.io)
