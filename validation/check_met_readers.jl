# check_met_readers.jl — reader conformance: does a meteorological reader hand the
# physics PHYSICAL fields in the common EASLEE layout?
#
#     julia -t 8 --project=easlee validation/check_met_readers.jl configs/<run>.toml [more.toml ...]
#
# Builds the runtime from each configuration (the first two met files go through
# the real reader and its derived-field calculation) and checks, on both loaded
# time slots, what every reader must deliver whatever the dataset:
#
#   layout      zlev(k=1) = 0, zlev strictly increasing, pressure decreasing with k
#   ranges      T, p, rho, qv, winds, vertical speed inside physical bounds
#   thermo      theta = T (1e5/p)^0.2857;  rho = p / (Rd T (1 + 0.608 q))
#   hydrostatic dp/dz = -rho g between levels (catches heights and pressures that
#               belong to different levels, e.g. interface heights on mass fields)
#   turbulence  tke filled (the Hanna fill ran), ustar and pblh in range
#   stability   the Obukhov length is negative (unstable) around local noon and
#               positive (stable) around local midnight over land — the sign test
#               that an inverted heat-flux convention fails
#   wet dep     tpf >= 0 and cc in [0, 1] when the reader provides them
#
# Written 2026-10-13 after three reader defects reached users (WRF thermodynamics
# never written back, wet deposition tied to the convective mass flux, ERA5 heat
# flux sign).  Each would have failed one of these checks on the first file.
# Exit status 1 when any check fails.

using Dates, Printf, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")
global gpu = false

const RD = 287.05; const G0 = 9.81
nfail_total = 0

pct(x, q) = isempty(x) ? NaN : quantile(x, q)
fin(x) = filter(isfinite, vec(x))

function check!(res, label, ok, detail)
    push!(res, (label, ok, detail))
    @printf("    %-58s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

function check_slot(d3, d2, gc, Ni, Nj, Nk, t, vdate)
    res = Any[]
    T  = Float64.(d3["T"][1:Ni, 1:Nj, :, t])
    p  = Float64.(d3["puv"][1:Ni, 1:Nj, :, t])
    z  = Float64.(d3["zlev"][1:Ni, 1:Nj, :, t])
    th = Float64.(d3["theta"][1:Ni, 1:Nj, :, t])
    rho= Float64.(d3["rho"][1:Ni, 1:Nj, :, t])
    q  = haskey(d3, "qv") ? Float64.(d3["qv"][1:Ni, 1:Nj, :, t]) : zeros(size(T))
    u  = Float64.(d3["Ux"][1:Ni, 1:Nj, :, t]); v = Float64.(d3["Vy"][1:Ni, 1:Nj, :, t])
    lev= Float64.(d3["levdot"][1:Ni, 1:Nj, :, t])

    # ── missing values: bitmap-missing points (AROME-HDF corners outside the
    # native domain) are NaN by design; every check below uses finite values only
    nnan = count(isnan, T) + count(isnan, u)
    nnan > 0 && @printf("    %-58s %s   %s\n", "NaN (bitmap-missing) points", "INFO",
                        @sprintf("T %.1f %%, u %.1f %% of the values", 100count(isnan, T) / length(T), 100count(isnan, u) / length(u)))
    # ── layout
    check!(res, "zlev(k=1) = 0", maximum(abs, z[:,:,1]) < 1.0, @sprintf("max |z1| = %.2f m", maximum(abs, z[:,:,1])))
    dz = diff(z; dims=3)
    check!(res, "zlev strictly increasing", count(<=(0.0), dz) == 0,
           @sprintf("%d non-increasing steps of %d; top median %.1f km", count(<=(0.0), dz), length(dz), median(fin(z[:,:,Nk])) / 1e3))
    dp = diff(p; dims=3)
    check!(res, "pressure decreasing with k", count(>=(0.0), dp) <= 0.001 * length(dp),
           @sprintf("%d non-decreasing steps of %d", count(>=(0.0), dp), length(dp)))
    if haskey(d2, "sp")
        sp = Float64.(d2["sp"][1:Ni, 1:Nj, t])
        r = fin(p[:,:,1] ./ sp)
        check!(res, "puv(k=1) = surface pressure", pct(abs.(r .- 1), 0.99) < 0.01, @sprintf("p1/sp p1-p99 %.4f-%.4f", pct(r, 0.01), pct(r, 0.99)))
    end
    # ── ranges
    fT, fp, fr, fq = fin(T), fin(p), fin(rho), fin(q)
    check!(res, "T in [150, 340] K", minimum(fT) >= 150 && maximum(fT) <= 340, @sprintf("[%.1f, %.1f]", minimum(fT), maximum(fT)))
    check!(res, "p in (0, 110000] Pa", minimum(fp) > 0 && maximum(fp) <= 1.1e5, @sprintf("[%.3g, %.4g]", minimum(fp), maximum(fp)))
    check!(res, "rho in (0, 1.6] kg/m3", minimum(fr) > 0 && maximum(fr) <= 1.6, @sprintf("[%.3g, %.3f]", minimum(fr), maximum(fr)))
    check!(res, "qv in [0, 0.04]", minimum(fq) >= -1e-5 && maximum(fq) <= 0.04, @sprintf("[%.2g, %.4f]", minimum(fq), maximum(fq)))
    # the winter polar-night jet reaches 150-190 m/s near the stratopause
    # (ERA5 June 2020: 183 m/s at 0.6 hPa, 42 S), so the tight bound is below 30 km
    spd = hypot.(u, v); lo = copy(spd); lo[z .> 30e3] .= 0
    check!(res, "horizontal wind < 150 m/s below 30 km, < 250 m/s above", maximum(fin(lo)) < 150 && maximum(fin(spd)) < 250,
           @sprintf("max %.1f m/s below 30 km, %.1f overall", maximum(fin(lo)), maximum(fin(spd))))
    # vertical speed implied by levdot [level/s] and the local level spacing
    w = Float64[]
    for k in 2:Nk-1
        append!(w, vec(lev[:,:,k] .* (z[:,:,k+1] .- z[:,:,k-1]) ./ 2))
    end
    w = fin(w)
    check!(res, "implied vertical speed: 99.9 % below 5 m/s", pct(abs.(w), 0.999) < 5.0,
           @sprintf("p99.9 |w| = %.3f m/s, max %.2f", pct(abs.(w), 0.999), maximum(abs, w)))
    # ── thermodynamics
    thc = T .* (1e5 ./ p) .^ 0.2857
    e = fin(abs.(th ./ thc .- 1))
    check!(res, "theta = T (1e5/p)^0.2857 (0.3 %)", pct(e, 0.999) < 3e-3, @sprintf("p99.9 rel err %.2e", pct(e, 0.999)))
    rc = p ./ (RD .* T .* (1 .+ 0.608 .* max.(q, 0)))
    rd = p ./ (RD .* T)
    e1 = fin(abs.(rho ./ rc .- 1)); e2 = fin(abs.(rho ./ rd .- 1))
    check!(res, "rho = p / (Rd T [1 + 0.608 q]) (2 %)", min(pct(e1, 0.99), pct(e2, 0.99)) < 0.02,
           @sprintf("p99 rel err moist %.2e, dry %.2e", pct(e1, 0.99), pct(e2, 0.99)))
    # hydrostatic balance between levels k-1 and k+1, above 1 km and below 30 km
    hy = Float64[]
    for k in 2:Nk-1, j in 1:Nj, i in 1:Ni
        zc = z[i,j,k]; (1000 < zc < 30000) || continue
        dzz = z[i,j,k+1] - z[i,j,k-1]; dzz > 0 || continue
        push!(hy, ((p[i,j,k+1] - p[i,j,k-1]) / dzz) / (-rho[i,j,k] * G0) - 1)
    end
    hy = fin(hy)
    # The transport code reads zlev and levdot on HALF levels and u, v, T, p, rho on
    # the FULL levels stacked between them (interpolation4Dsingle reads index
    # z + 0.5; the ERA5 reader's layout): the full level of slot s lies between
    # zlev[s-1] and zlev[s].  Hydrostatic balance with those full-level heights is
    # the layout contract; the same-slot version is printed for comparison (a
    # reader that puts zlev on the full levels passes that one and fails this one).
    hs = Float64[]
    for k in 3:Nk-1, j in 1:Nj, i in 1:Ni
        zf(s) = 0.5 * (z[i,j,s-1] + z[i,j,s])
        zc = zf(k); (1000 < zc < 30000) || continue
        dzz = zf(k+1) - zf(k-1); dzz > 0 || continue
        push!(hs, ((p[i,j,k+1] - p[i,j,k-1]) / dzz) / (-rho[i,j,k] * G0) - 1)
    end
    hs = fin(hs)
    check!(res, "hydrostatic, stacked layout (full level between zlev[s-1], zlev[s])", abs(median(hs)) < 0.01,
           @sprintf("median %+.4f, p5-p95 %+.3f..%+.3f (n=%d)", median(hs), pct(hs, 0.05), pct(hs, 0.95), length(hs)))
    @printf("    %-58s %s   %s\n", "hydrostatic, same-slot layout (for comparison)", "INFO",
            @sprintf("median %+.4f, p5-p95 %+.3f..%+.3f", median(hy), pct(hy, 0.05), pct(hy, 0.95)))
    # ── turbulence fields
    tke = Float64.(d3["tke"][1:Ni, 1:Nj, :, t])
    ftke = fin(tke)
    check!(res, "tke filled (Hanna fill or native)", !isempty(ftke) && maximum(ftke) > 0, @sprintf("max %.3g m2/s2", isempty(ftke) ? NaN : maximum(ftke)))
    us = fin(Float64.(d2["ustar"][1:Ni, 1:Nj, t])); bl = fin(Float64.(d2["pblh"][1:Ni, 1:Nj, t]))
    check!(res, "ustar in (0, 3] m/s", !isempty(us) && minimum(us) >= 0 && maximum(us) <= 3 && median(us) > 0.01,
           @sprintf("median %.3f, max %.2f", median(us), maximum(us)))
    check!(res, "pblh in [10, 8000] m", !isempty(bl) && minimum(bl) >= 9.5 && maximum(bl) <= 8000,
           @sprintf("[%.0f, %.0f], median %.0f", minimum(bl), maximum(bl), median(bl)))
    # ── stability sign by local solar time (over land when a land mask exists)
    ol = Float64.(d2["ol"][1:Ni, 1:Nj, t])
    lsm = haskey(d2, "lsm") ? Float64.(d2["lsm"][1:Ni, 1:Nj, t]) : nothing
    if (lsm === nothing || maximum(lsm) <= 0.5)
        # no land mask (MERRA-2): surface elevation above 100 m as the land proxy
        for (key, scale) in (("phis", 1 / G0), ("orog", 1.0), ("z", 1 / G0))
            haskey(d2, key) || continue
            zs = Float64.(d2[key][1:Ni, 1:Nj, t]) .* scale
            maximum(zs) > 100 || continue
            lsm = Float64.(zs .> 100); break
        end
    end
    island = lsm !== nothing && maximum(lsm) > 0.5
    day = Float64[]; night = Float64[]; dayblh = Float64[]; nightblh = Float64[]
    hutc = Dates.hour(vdate) + Dates.minute(vdate) / 60
    for j in 1:Nj, i in 1:Ni
        island && lsm[i,j] < 0.5 && continue
        lat, lon = EASLEE.Operators.indyindxtolatlon(gc, Float32(j), Float32(i))
        abs(lat) < 60 || continue
        lst = mod(hutc + lon / 15, 24)
        if 11 <= lst <= 14
            push!(day, ol[i,j]); push!(dayblh, d2["pblh"][i,j,t])
        elseif lst >= 23 || lst <= 2
            push!(night, ol[i,j]); push!(nightblh, d2["pblh"][i,j,t])
        end
    end
    # noon and midnight tested separately: a regional domain (WRF, AROME) has
    # only one of them in a given slot (2026-10-21: both were required, so the
    # sign was never checked on regional data)
    day = fin(day); night = fin(night)
    if length(day) > 20
        check!(res, "Obukhov L < 0 (unstable) at local noon over land", median(day) < 0,
               @sprintf("median L %.0f m (n=%d), unstable share %.2f", median(day), length(day), count(<(0), day) / length(day)))
    end
    if length(night) > 20
        check!(res, "Obukhov L > 0 (stable) at local midnight over land", median(night) > 0,
               @sprintf("median L %.0f m (n=%d), stable share %.2f", median(night), length(night), count(>(0), night) / length(night)))
    end
    dayblh = fin(dayblh); nightblh = fin(nightblh)
    if length(dayblh) > 20 && length(nightblh) > 20
        check!(res, "pblh deeper at local noon than at midnight (land)", median(dayblh) > median(nightblh),
               @sprintf("median %.0f m vs %.0f m", median(dayblh), median(nightblh)))
    end
    (length(day) > 20 || length(night) > 20) ||
        @printf("    %-58s %s   %s\n", "Obukhov sign by local time", "SKIP", "no local-noon or local-midnight land cells in this slot")
    # ── wet deposition inputs
    if haskey(d3, "tpf")
        tpf = fin(Float64.(d3["tpf"][1:Ni, 1:Nj, :, t]))
        check!(res, "tpf >= 0 and not all zero", minimum(tpf) >= 0 && maximum(tpf) > 0,
               @sprintf("max %.3g kg/m2/s = %.1f mm/h", maximum(tpf), 3600 * maximum(tpf)))
    end
    if haskey(d3, "cc")
        cc = fin(Float64.(d3["cc"][1:Ni, 1:Nj, :, t]))
        check!(res, "cc in [0, 1] and not all zero", minimum(cc) >= 0 && maximum(cc) <= 1 && maximum(cc) > 0,
               @sprintf("[%.3f, %.3f], mean %.3f", minimum(cc), maximum(cc), mean(cc)))
    end
    return res
end

for cfg in ARGS
    global nfail_total
    println("\n══ ", cfg)
    spec = read_runspec(cfg)
    rt, obj = build_run(spec; prefetch=false, extra_slots=10)
    ctx = rt.ctx
    read_met_step!(rt, ctx.runconfig.startdate)        # loads the two bracketing met files
    for n in 1:length(ctx.nests)
        nc = ctx.nests[n]; gc = nc.gridconfig
        Ni, Nj = Int(gc.Ni), Int(gc.Nj); Nk = size(rt.metfields3D[n]["zlev"], 3)
        tix = nc.nestrunconfig.t_index; fd = nc.metconfig.filedate
        println(@sprintf("  nest %d: %s  Ni=%d Nj=%d Nk=%d", n, nc.metconfig.model, Ni, Nj, Nk))
        wf = hasproperty(obj, :wetfields) ? obj.wetfields : obj.convfields   # (runs on older trees too)
        @printf("  wet deposition fields: %s; convection fields: %s\n",
                wf === nothing ? "none" : "present", obj.convfields === nothing ? "none" : "present")
        for t in 1:2
            vdate = fd[tix[t]]
            println("  slot $t, valid $vdate")
            res = check_slot(rt.metfields3D[n], rt.metfields2D[n], gc, Ni, Nj, Nk, t, vdate)
            nfail_total += count(r -> !r[2], res)
        end
    end
end
println(nfail_total == 0 ? "\nALL CHECKS PASS" : "\n$nfail_total CHECK(S) FAIL")
exit(nfail_total == 0 ? 0 : 1)
