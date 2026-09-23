
# ── Concentration grid accumulation ───────────────────────────────────────────
# Deposits parcel mass into the 5-D concentration array
# [Niout, Njout, Nkout, Nc, nage].
#
# Type discipline:
# - grid/massarray are type-asserted to Array{Float32,5}/Array{Float64,5}.
# - Niout/Njout/Nc/age are asserted to Int so comparisons don't box Float32.
# - The outer loop is over nests (not parcels) so all per-nest Dict/Vector
#   lookups happen once outside the hot parcel loop.  Inside the inner loop
#   every variable is concretely typed → zero boxing → near-zero allocations.

# ── Vertical output-layer lookup ─────────────────────────────────────────────

function find_out_layer(z::Float32, outheight::AbstractVector{Float32})
    # outheight is a list of N interface heights (lower-to-upper edges).
    # This creates N-1 physical layers, each spanning [outheight[k], outheight[k+1]].
    # Starting the search from k=2 means z in [outheight[1], outheight[2]) goes to
    # layer k=2 — which has the correct dz = outheight[2]-outheight[1].
    # The degenerate k=1 slot (a zero-thickness "point" at the surface) is never
    # populated, so the sentinel dz=1m in initoutgrid never inflates any cell.
    #
    # Returns 0 if z is strictly below the bottom edge or above the top edge.
    (z < outheight[1] || z > outheight[end]) && return 0
    @inbounds for k in 2:length(outheight)
        outheight[k] >= z && return k
    end
    return 0
end

# ── Adaptive output kernel (ADAPTIVE_KERNEL design, work item A) ──────────────
#
# Each parcel carries a STATISTICAL bandwidth h_i [m] — the RMS distance of
# its k nearest neighbours (adaptive KDE spacing, k = 12), NOT a physical puff
# size.  The gridded field is the sum of per-parcel truncated Gaussians
# renormalised over the in-domain support, so mass is conserved exactly
# (including at domain edges, unlike the uniform kernel) and the FLEXPART 3-h
# age rule disappears: at a fresh point release the spacing → 0, the kernel
# degenerates continuously to single-cell deposition.  R(x) ∝ N^(−1/2):
# resolution is bought with trajectories, not with the output grid.

# erf via Abramowitz–Stegun 7.1.26 (max abs error 1.5e-7 — far below the
# Float32 deposit precision; avoids a SpecialFunctions dependency).
@inline function _erf_as(x::Float64)
    s = x < 0.0 ? -1.0 : 1.0
    ax = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    y = 1.0 - (((((1.061405429t - 1.453152027)t) + 1.421413741)t - 0.284496736)t +
               0.254829592) * t * exp(-ax * ax)
    return s * y
end
@inline _Φ(t::Float64) = 0.5 * (1.0 + _erf_as(t / 1.4142135623730951))

# metres per index unit of the met grid (mirrors mixing.jl's _grid_spacing_m;
# duplicated because conc.jl lives in Operators and mixing.jl in EASLEE).
_out_spacing_m(g::Lambert)  = (Float64(g.dx), Float64(g.dy))
_out_spacing_m(g::Mercator) = (Float64(g.dx), Float64(g.dy))
function _out_spacing_m(g::Rectilinear)
    # EQUATOR spacing — the lon convergence is applied PER PARCEL via _xcos
    # (a single grid-mid cos compressed the kNN metric anisotropically for any
    # population away from the grid middle; the old mid-lat formula also mixed
    # output-index units with degrees).
    return 1.0 / Float64(g.dxconst), 1.0 / Float64(g.dyconst)
end

# per-parcel lon-convergence factor: metres per index unit of x at the
# parcel's own latitude (local flat approximation — kNN and the deposit
# stencil are local)
@inline _xcos(g::Rectilinear, indy) =
    max(cosd(Float64(g.latmin) + (Float64(indy) - 1.0) * Float64(g.dlat)), 0.05)
@inline _xcos(g, indy) = 1.0

# OUTPUT-grid cell size in metres: one output index = 1/coeff met indices
# (indxtoxout: xout = indx·coeffx + offsetx).  Mercator maps through lon/lat
# (no coeff fields) → fall back to met spacing (small-domain MesoNH; the
# adaptive kernel is validated on Rectilinear/Lambert).
_out_cell_m(g::Rectilinear) = ((a, b) = _out_spacing_m(g);
                                     (a / Float64(g.coeffx), b / Float64(g.coeffy)))
_out_cell_m(g::Lambert)     = ((a, b) = _out_spacing_m(g);
                                     (a / Float64(g.coeffx), b / Float64(g.coeffy)))
_out_cell_m(g)                    = _out_spacing_m(g)

# Physical bandwidth floor (relative dispersion): the kNN hband measures the
# SAMPLING spacing of the parcel cloud, which collapses to metres inside a
# young dense plume — far below the width of atmosphere one parcel represents
# after turbulent transport.  A point evaluation there reports single-parcel
# spikes (ETEX station 8021: h ≈ 60 m ⇒ one 0.85 g parcel ⇒ ~7e-10 kg/m³),
# and the grid deposit paints a one-cell filament.  When enabled
# ([output] bandwidth_floor = "dispersion"), h is floored by the Taylor (1921)
# absolute-dispersion growth
#     σ²(t) = 2 σv² T_L² (t/T_L − 1 + e^(−t/T_L))
# (ballistic σ ≈ σv·t, diffusive σ ≈ σv·√(2 T_L t)), σv an EFFECTIVE
# unresolved-relative-dispersion velocity.  Defaults σv = 0.5 m/s, T_L = 300 s
# → ≈ 0.8 km at 1 h, 3.6 km at 24 h, 6.2 km at 72 h.  R(x) then reports
# max(statistical, physical) — the true resolution either way.
@inline function _hfloor_taylor(age_s::Float64, sv::Float64, tl::Float64)
    (age_s > 0.0 && sv > 0.0 && tl > 0.0) || return 0.0
    x = age_s / tl
    return sv * tl * sqrt(max(2.0 * (x - 1.0 + exp(-min(x, 60.0))), 0.0))
end
@inline _parcel_age_s(state, i) =
    (Dates.value(state.Dtime[i]) - Dates.value(state.Dtimeinit[i])) * 1.0e-3

"""
    update_bandwidths!(state, gconfig, nn; k=12, floor_sv=0, floor_tl=300) → n_updated

kNN bandwidth pass (output sampling cadence): for every active parcel of nest
`nn`, h_i² = mean squared horizontal distance of its k nearest neighbours.
Cell-list build at L₀ = √(A_bbox/N)·√k (mean-density guess); parcels with
fewer than k neighbours in their 3×3 block take the RMS of what was found —
that isolated a parcel genuinely has that spacing (no second pass, no special
case).  Floor h ≥ 1 m (0/0 guard only).  Do NOT reuse the mixing pairlist:
its neighbour cap keeps encounter order, not the nearest (would bias h high
exactly where it must collapse).
"""
# Population map for per-population bandwidths (task #35): release box index →
# population id (usually the SPEC release index, set by build_run).  Empty →
# every parcel is one population.  With ≥2 co-located populations (domainfill
# background + concentrated release plume) a pooled kNN is wrong both ways:
# background parcels inside the plume would collapse their h to the plume
# spacing (spurious background structure), and the plume h would inflate near
# its edge.  Each population is a separate KDE — its bandwidth must come from
# its OWN sample.  Parcels with relbox outside the map (0 = boundary/domainfill
# injections) form population 0.
const OUT_POPMAP = Ref(Int[])
set_output_popmap!(v::AbstractVector{<:Integer}) = (OUT_POPMAP[] = collect(Int, v))

function update_bandwidths!(state::State.ParcelState, gconfig, nn::Int; k::Int = 12,
                            floor_sv::Float64 = 0.0, floor_tl::Float64 = 300.0)
    popmap = OUT_POPMAP[]
    if isempty(popmap)
        return _update_bandwidths_pop!(state, gconfig, nn, nothing; k=k,
                                       floor_sv=floor_sv, floor_tl=floor_tl)
    end
    pops = Set{Int}()
    @inbounds for i in eachindex(state.relbox)
        state.status[i] == 0x01 || continue
        rb = Int(state.relbox[i])
        push!(pops, 1 <= rb <= length(popmap) ? popmap[rb] : 0)
    end
    n = 0
    for p in pops
        n += _update_bandwidths_pop!(state, gconfig, nn, p; k=k,
                                     floor_sv=floor_sv, floor_tl=floor_tl)
    end
    return n
end

# One kNN pass restricted to population `pop` (nothing → all active parcels).
function _update_bandwidths_pop!(state::State.ParcelState, gconfig, nn::Int,
                                 pop::Union{Nothing,Int}; k::Int = 12,
                                 floor_sv::Float64 = 0.0, floor_tl::Float64 = 300.0)
    popmap = OUT_POPMAP[]
    @inline _pop(i) = begin
        rb = Int(state.relbox[i])
        1 <= rb <= length(popmap) ? popmap[rb] : 0
    end
    dxm, dym = _out_spacing_m(gconfig)
    nn8 = Int8(nn)
    N = length(state.indx)
    x = Vector{Float64}(undef, N); y = Vector{Float64}(undef, N)
    act = falses(N)
    # Accumulate the bbox in loop-local names and bind xmin/ymin ONCE below:
    # the @threads body is a closure, and a captured variable that is
    # REASSIGNED anywhere in the function is captured as Core.Box — the cell
    # index computed from it goes dynamic and poisons the whole inner distance
    # loop (measured 167 GiB boxed Float64/Int64 allocs and ~80 s per 400k
    # refresh; single-assignment bindings: same algorithm, allocation-free).
    xmn = Inf; xmx = -Inf; ymn = Inf; ymx = -Inf
    nact = 0
    @inbounds for i in 1:N
        (state.status[i] == 0x01 && state.nest[i] == nn8) || continue
        pop === nothing || _pop(i) == pop || continue
        act[i] = true
        x[i] = Float64(state.indx[i]) * dxm
        y[i] = Float64(state.indy[i]) * dym
        xmn = min(xmn, x[i]); xmx = max(xmx, x[i])
        ymn = min(ymn, y[i]); ymx = max(ymx, y[i])
        nact += 1
    end
    nact <= 1 && return 0
    xmin = xmn; ymin = ymn        # single-assignment (typed closure capture)
    xmax = xmx; ymax = ymx
    L0 = max(sqrt(max(xmx - xmn, 1.0) * max(ymx - ymn, 1.0) / nact) * sqrt(k), 1.0)

    # CSR cell list at cell size L0 (same layout as mixing's _build_cell_list_2d)
    ncx = max(1, ceil(Int, (xmax - xmin) / L0) + 1)
    ncy = max(1, ceil(Int, (ymax - ymin) / L0) + 1)
    cid = zeros(Int, N)
    starts = zeros(Int, ncx * ncy + 1)
    @inbounds for i in 1:N
        act[i] || continue
        ix = clamp(floor(Int, (x[i] - xmin) / L0) + 1, 1, ncx)
        iy = clamp(floor(Int, (y[i] - ymin) / L0) + 1, 1, ncy)
        c = ix + (iy - 1) * ncx
        cid[i] = c; starts[c + 1] += 1
    end
    starts[1] = 1
    @inbounds for c in 1:(ncx * ncy); starts[c + 1] += starts[c]; end
    items = Vector{Int}(undef, nact)
    cursor = copy(starts)
    @inbounds for i in 1:N
        cid[i] == 0 && continue
        items[cursor[cid[i]]] = i; cursor[cid[i]] += 1
    end

    # Per-parcel k-smallest squared distances over the 3×3 block, with a HARD
    # WORK BOUND for clustered populations.  L0 is a global mean-density
    # scale: a point-release plume (HT) puts ~all N parcels in ONE cell and
    # the plain block scan goes O(N²) — the 2026-07-06 HT run sat here for
    # hours.  Over-occupied blocks are scanned with a systematic STRIDE
    # subsample instead: the k-th neighbour among a 1/s subsample estimates
    # the (k·s)-th of the full set, and for a locally uniform 2-D field
    # r²(m) ∝ m, so h²_true = h²_sub / s — unbiased, no extra memory, and
    # ≤ 9·CAPC distance evaluations per parcel whatever the crowding.
    CAPC = 6 * k                             # per-cell scan budget (tight: a
                                             # dense cell is stride-subsampled
                                             # to ~6k candidates, the rank
                                             # correction recovers the density)
    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    bufs = [Vector{Float64}(undef, k) for _ in 1:nthr]
    Threads.@threads :static for i in 1:N
        act[i] || continue
        buf = bufs[Threads.threadid()]
        nf = 0                                 # entries in buf (ascending)
        cx = clamp(floor(Int, (x[i] - xmin) / L0) + 1, 1, ncx)
        cy = clamp(floor(Int, (y[i] - ymin) / L0) + 1, 1, ncy)
        # lon convergence at the parcel's OWN latitude: the cell list lives in
        # equator coordinates (a per-parcel-scaled embedding would fake
        # distances between parcels at different lats); distances apply ci at
        # comparison time, and the x-scan reach widens by 1/ci so the block
        # still covers 3·L0 PHYSICAL metres.
        ci = _xcos(gconfig, state.indy[i])
        ex = ci >= 0.999 ? 1 : min(ceil(Int, 1.0 / ci), 8)
        # one stride for the whole block (max cell occupancy) so a single
        # rank correction applies to the mixed buffer
        occmax = 0
        @inbounds for oy in -1:1, ox in -ex:ex
            gx = cx + ox; gy = cy + oy
            (1 <= gx <= ncx && 1 <= gy <= ncy) || continue
            c = gx + (gy - 1) * ncx
            o = starts[c + 1] - starts[c]
            o > occmax && (occmax = o)
        end
        stride = occmax > CAPC ? cld(occmax, CAPC) : 1
        off = stride > 1 ? Int(mod(hash(i), stride)) : 0   # decorrelate parcels
        @inbounds for oy in -1:1, ox in -ex:ex
            gx = cx + ox; gy = cy + oy
            (1 <= gx <= ncx && 1 <= gy <= ncy) || continue
            c = gx + (gy - 1) * ncx
            for m in (starts[c] + off):stride:(starts[c + 1] - 1)
                j = items[m]
                j == i && continue
                dxp = (x[j] - x[i]) * ci; dyp = y[j] - y[i]
                r2 = dxp * dxp + dyp * dyp
                if nf < k
                    p = nf
                    while p > 0 && buf[p] > r2; buf[p + 1] = buf[p]; p -= 1; end
                    buf[p + 1] = r2; nf += 1
                elseif r2 < buf[k]
                    p = k - 1
                    while p > 0 && buf[p] > r2; buf[p + 1] = buf[p]; p -= 1; end
                    buf[p + 1] = r2
                end
            end
        end
        h2 = 0.0
        if nf > 0
            @inbounds for m in 1:nf; h2 += buf[m]; end
            h2 /= (nf * stride)                # rank correction (stride=1 → exact)
        else
            h2 = L0 * L0                       # isolated beyond the block: spacing ≥ L0
        end
        hf = floor_sv > 0.0 ?
             _hfloor_taylor(_parcel_age_s(state, i), floor_sv, floor_tl) : 0.0
        state.hband[i] = Float32(max(sqrt(h2), 1.0, hf))
    end
    return nact
end

# ── Per-nest inner kernel ──────────────────────────────────────────────────────
# All arguments are concretely typed; the compiler can specialise fully and
# emit purely stack-allocated code with no dynamic dispatch.

function _conc_nest!(
    grid      ::Array{Float32,5},
    massarray ::Array{Float64,5},
    trajsin,
    currentdate,
    indt      ::Float32,
    Niout     ::Int,
    Njout     ::Int,
    Nc        ::Int,
    age       ::Int,
    zi        ::Vector{Float32},
    nn        ::Int,
    zlev_mat  ::Array{Float32,4},
    rho_mat   ::Array{Float32,4},
    lconfig,          # SigmaPressureLevel (concrete at call site via dispatch)
    gconfig,          # Rectilinear/Lambert (concrete at call site via dispatch)
    backward  ::Bool,            # backward run → sensitivity accumulation
    dtsync    ::Float64,         # |timestep| [s] (backward weight)
    relmass   ::Vector{Float32}, # released mass per output column [kg] (backward norm)
    kernelmode::Symbol = :uniform,   # :none | :uniform | :adaptive | :covariance
    resnum    ::Union{Nothing,Matrix{Float64}} = nothing,   # R(x) accumulators
    resden    ::Union{Nothing,Matrix{Float64}} = nothing,   # (adaptive kernel)
    bwfloor_sv::Float64 = 0.0,       # >0 → Taylor dispersion bandwidth floor σv [m/s]
    bwfloor_tl::Float64 = 300.0,     #      (also covers parcels born after the
                                     #      last kNN refresh, whose hband is 0)
    # diagnostic accumulators (returned as a NamedTuple to avoid mutation of
    # scalars across a potential inlining barrier)
    ;
    # release box → output column ([output].group); empty = identity (per-box)
    outclass  ::Vector{Int32} = Int32[],
)
    deposited_mass    = 0.0
    deposited_raw     = 0.0
    n_outside_xy      = 0
    n_outside_z       = 0
    n_below_ground    = 0
    n_below_first     = 0
    n_deposited       = 0
    mass_below_ground = 0.0
    mass_below_first  = 0.0

    nn8 = Int8(nn)

    # adaptive kernel: OUTPUT-cell size in metres (met spacing / cells-per-met-
    # index coeff) + reusable per-axis weight buffers (single-threaded loop)
    dxo, dyo = (kernelmode === :adaptive || kernelmode === :covariance ||
                kernelmode === :covshape) ?
               _out_cell_m(gconfig) : (1.0, 1.0)
    Wxb = Vector{Float64}(undef, 256)
    Wyb = Vector{Float64}(undef, 256)

    @inbounds for i in eachindex(trajsin.Dtime)
        trajsin.Dtime[i] == currentdate || continue
        # review 2026-09-04 A5: culled domain-fill parcels (status 0x00, Dtime
        # left at currentdate) were deposited one extra step — gate on status
        # like receptors.jl does.
        trajsin.status[i] == 0x01       || continue
        trajsin.nest[i]  == nn8         || continue

        indx = trajsin.indx[i]
        indy = trajsin.indy[i]
        indl = trajsin.indl[i]
        mass = trajsin.mass[i]

        zout::Float32 = indltozlev(lconfig, zlev_mat, indx, indy, indl, indt)
        xout::Float32 = indxtoxout(gconfig, indx, indy)
        yout::Float32 = indytoyout(gconfig, indx, indy)

        # Convert fractional OUTGRID coordinates to integer cell indices using floor
        # (column-major convention: cell ii covers xout ∈ [ii, ii+1)).
        # Important: use the integer index for the bounds check, NOT xout directly.
        # Reason: with a small positive offsetx, an injected parcel at the inflow
        # boundary (e.g. west MET face, indx≈0.5) maps to xout≈0.54 — which fails
        # the naive "xout >= 1" check even though it physically belongs to cell 1.
        # Similarly, parcels near the east outflow face map to xout slightly above
        # Niout but still belong in cell Niout.
        # Using ii = trunc(xout) and checking 1 ≤ ii ≤ Niout correctly handles both
        # boundary cases: xout ∈ [0.5, 1) → ii=0 rejected; xout ∈ [Niout, Niout+1) → ii=Niout accepted.
        ii  = Base.unsafe_trunc(Int, xout)
        jj  = Base.unsafe_trunc(Int, yout)

        if !(1 <= ii <= Niout && 1 <= jj <= Njout)
            n_outside_xy += 1
            continue
        end

        if zout < 0.0f0
            mass_below_ground += mass
            n_below_ground += 1
            continue
        end

        l = find_out_layer(zout, zi)
        if l == 0
            if zout < zi[1]
                mass_below_first += mass
                n_below_first += 1
            else
                n_outside_z += 1
            end
            continue
        end

        rho::Float32 = interpolation4Dsingle(rho_mat, indx, indy, indl, indt)
        (isfinite(rho) && rho > 0.0f0) || continue
        # Boundary-injected parcels have relbox=0; clamp to 1 so they contribute
        # to the total count. Python sums all boxes, so the box assignment is moot.
        # outclass maps the box to its accumulator column ([output].group —
        # identity under per-release grouping, species index under "species").
        box = Nc == 1 ? 1 : begin
            rb = max(1, Int(trajsin.relbox[i]))
            isempty(outclass) ? rb : Int(outclass[min(rb, length(outclass))])
        end

        # Forward: plain mass deposition (→ conc = Σm/V, mr = Σ(m/ρ)/V).
        # Backward: source-receptor sensitivity — each sample contributes the
        # residence weight Δt·(m_p/M_released) (Seibert & Frank 2004), so the
        # written grids are  sens = ΣΔt·(m/M)/V [s m⁻³]  and
        # sens_mr = ΣΔt·(m/(ρM))/V [s kg⁻¹].
        w = 1.0
        if backward
            M = box <= length(relmass) ? relmass[box] : 1.0f0
            w = dtsync / max(Float64(M), 1.0e-30)
        end
        addvol  = Float64(mass) / rho * w
        addmass = Float64(mass) * w

        if kernelmode === :adaptive || kernelmode === :covariance ||
           kernelmode === :covshape
            # Per-parcel truncated Gaussian (±3σ), renormalised over the
            # IN-DOMAIN touched cells → Σ weights = 1 exactly, incl. edges.
            # σ → floor collapses the support to one cell (single-cell limit,
            # no age branch).
            local sxm::Float64, sym::Float64, heff::Float64
            if kernelmode === :adaptive
                # bandwidth = kNN spacing (statistical), isotropic
                h  = Float64(trajsin.hband[i])
                h < 1.0 && (h = 1.0)
                if bwfloor_sv > 0.0
                    hf = _hfloor_taylor(_parcel_age_s(trajsin, i), bwfloor_sv, bwfloor_tl)
                    h < hf && (h = hf)
                end
                sxm = h; sym = h; heff = h
            elseif kernelmode === :covshape
                # SHAPE mode (no double counting): the kernel AREA is the kNN
                # spacing h² — the statistical smoothing scale, vanishes with
                # N — and only the ANISOTROPY comes from the flow-evolved P
                # (shear-aligned smoothing: along the filament, not across).
                # Aspect σx/σy = √(Pxx/Pyy), capped 8:1; heff = h (R(x) keeps
                # its resolution meaning).
                h = Float64(trajsin.hband[i])
                h < 1.0 && (h = 1.0)
                if bwfloor_sv > 0.0
                    hf = _hfloor_taylor(_parcel_age_s(trajsin, i), bwfloor_sv, bwfloor_tl)
                    h < hf && (h = hf)
                end
                ca = max(Float64(trajsin.cxx[i]), 1.0)
                cd = max(Float64(trajsin.cyy[i]), 1.0)
                asp = clamp(sqrt(ca / cd), 0.125, 8.0)   # σx/σy
                sxm = h * sqrt(asp); sym = h / sqrt(asp)
                heff = h
            else
                # flow-evolved covariance (physical footprint): per-axis widths
                # √Pxx / √Pyy — the deposit stays separable → conservation is
                # exact; the cxy tilt is honoured at the receptors (bivariate).
                # R(x) reports the area-equivalent width (det P)^(1/4).
                ca = Float64(trajsin.cxx[i]); cd = Float64(trajsin.cyy[i])
                cb = Float64(trajsin.cxy[i])
                ca < 1.0 && (ca = 1.0); cd < 1.0 && (cd = 1.0)
                sxm = sqrt(ca); sym = sqrt(cd)
                cdet = ca * cd - cb * cb
                heff = cdet > 0.0 ? cdet^0.25 : min(sxm, sym)
            end
            # metres → output-cell units, with the lon convergence at the
            # parcel's own latitude (dxo is the EQUATOR cell size)
            sx = sxm / (dxo * _xcos(gconfig, indy));  sy = sym / dyo
            xf = Float64(xout); yf = Float64(yout)
            # Support = ±3σ, but HARD-CAPPED at ±32 cells per axis.  Without
            # the cap, sparse-tail parcels (h ≫ cell) touch up to 256×256
            # cells EACH — the 2025-07-05 ETEX-A+B / HT runs spent hours in
            # this loop.  The weights are renormalised over the touched
            # window anyway (mass exact); for an ultra-diffuse parcel the
            # capped window only trims a near-uniform tail.
            i_lo = max(1,     Base.unsafe_trunc(Int, floor(max(xf - 3sx, xf - 32.0))))
            i_hi = min(Niout, Base.unsafe_trunc(Int, floor(min(xf + 3sx, xf + 32.0))))
            j_lo = max(1,     Base.unsafe_trunc(Int, floor(max(yf - 3sy, yf - 32.0))))
            j_hi = min(Njout, Base.unsafe_trunc(Int, floor(min(yf + 3sy, yf + 32.0))))
            if i_hi < i_lo || j_hi < j_lo
                n_outside_xy += 1
                continue
            end
            nx = min(i_hi - i_lo + 1, length(Wxb))
            ny = min(j_hi - j_lo + 1, length(Wyb))
            sWx = 0.0
            @inbounds for m in 1:nx
                w = _Φ((i_lo + m - xf) / sx) - _Φ((i_lo + m - 1.0 - xf) / sx)
                Wxb[m] = w; sWx += w
            end
            sWy = 0.0
            @inbounds for m in 1:ny
                w = _Φ((j_lo + m - yf) / sy) - _Φ((j_lo + m - 1.0 - yf) / sy)
                Wyb[m] = w; sWy += w
            end
            Wtot = sWx * sWy
            if !(Wtot > 0.0)
                # numerically empty support (h ≪ cell, xf at a cell centre far
                # from edges cannot happen — but guard): nearest-cell deposit
                grid[ii, jj, l, box, age]      += Float32(addvol)
                massarray[ii, jj, l, box, age] += addmass
            else
                @inbounds for mj in 1:ny, mi in 1:nx
                    w = Wxb[mi] * Wyb[mj] / Wtot
                    w > 0.0 || continue
                    ci = i_lo + mi - 1; cj = j_lo + mj - 1
                    grid[ci, cj, l, box, age]      += Float32(addvol * w)
                    massarray[ci, cj, l, box, age] += addmass * w
                    if resnum !== nothing
                        resnum[ci, cj] += addmass * w * heff   # R(x) accumulation
                        resden[ci, cj] += addmass * w
                    end
                end
            end
        elseif kernelmode === :uniform
            # FLEXPART uniform kernel: spread the mass bilinearly over the 4 nearest
            # cell centres (centre of cell k is at k+0.5, so shift by −0.5 to bracket).
            gx = xout - 0.5f0;  gy = yout - 0.5f0
            i0 = Base.unsafe_trunc(Int, floor(gx));  j0 = Base.unsafe_trunc(Int, floor(gy))
            dx = Float64(gx - i0);  dy = Float64(gy - j0)
            @inbounds for (ci, cj, wgt) in ((i0, j0, (1-dx)*(1-dy)), (i0+1, j0, dx*(1-dy)),
                                            (i0, j0+1, (1-dx)*dy),   (i0+1, j0+1, dx*dy))
                (wgt > 0 && 1 <= ci <= Niout && 1 <= cj <= Njout) || continue
                grid[ci, cj, l, box, age]      += Float32(addvol * wgt)
                massarray[ci, cj, l, box, age] += addmass * wgt
            end
        else
            grid[ii, jj, l, box, age]      += Float32(addvol)
            massarray[ii, jj, l, box, age] += addmass
        end

        deposited_mass += addmass
        deposited_raw  += addvol
        n_deposited    += 1
    end

    return (deposited_mass=deposited_mass, deposited_raw=deposited_raw,
            n_outside_xy=n_outside_xy, n_outside_z=n_outside_z,
            n_below_ground=n_below_ground, n_below_first=n_below_first,
            n_deposited=n_deposited,
            mass_below_ground=mass_below_ground, mass_below_first=mass_below_first)
end

# ── Main accumulation function ────────────────────────────────────────────────

const _DEBUG_CONC = false   # set true to print per-step accumulation diagnostics

function conccalc!(oc::OutputConfig, trajsin, rc::RunConfig, t_array, metdata3Din, gridconfig, levelconfig, gpu)
    currentdate = rc.currentdate

    grid      = oc.concarray
    massarray = oc.massarray
    Niout::Int = Int(oc.Niout)
    Njout::Int = Int(oc.Njout)
    Nc::Int    = rc.Nc
    age::Int   = 1
    zi::Vector{Float32} = oc.outinterfaces

    nnest = length(metdata3Din)
    bw_refresh = mod(oc.bw_count, max(oc.bw_stride, 1)) == 0   # this deposit's cadence
    oc.bw_count += 1

    # Accumulate diagnostics across nests.
    total_mass    = 0.0; total_raw     = 0.0
    tot_dep       = 0;   tot_oxy       = 0; tot_oxz = 0
    tot_bg        = 0;   tot_bf        = 0
    mass_bg       = 0.0; mass_bf       = 0.0

    for nn in 1:nnest
        # Extract all per-nest objects ONCE, outside the parcel loop.
        # metdata3Din[nn] is now a MetFields3D struct — field access is a direct
        # memory read with no dict-lookup overhead.
        zlev_mat = metdata3Din[nn].zlev::Array{Float32,4}
        rho_mat  = metdata3Din[nn].rho::Array{Float32,4}
        lconfig  = levelconfig[nn]   # SigmaPressureLevel in practice
        gconfig  = gridconfig[nn]    # Rectilinear/Lambert in practice
        indt     = t_array[nn]

        # adaptive kernel: refresh per-parcel bandwidths only every bw_stride
        # deposit calls (once per output frame by default — the kNN pass is the
        # dominant output cost; hband persists on the ParcelState between).
        if (oc.kernel === :adaptive || oc.kernel === :covshape) && bw_refresh
            update_bandwidths!(trajsin, gconfig, nn;
                               floor_sv=oc.bwfloor_sv, floor_tl=oc.bwfloor_tl)
        end

        # _conc_nest! is specialised for the concrete types of lconfig/gconfig
        # at the call site (Julia method dispatch): this is the "dispatch barrier"
        # that gives the inner loop fully concrete types.
        d = _conc_nest!(grid, massarray, trajsin, currentdate,
                        indt, Niout, Njout, Nc, age, zi, nn,
                        zlev_mat, rho_mat, lconfig, gconfig,
                        rc.direction == -1, abs(rc.timestep.value) / 1000.0,
                        oc.relmass, oc.kernel,
                        oc.resolution ? oc.resnum : nothing,
                        oc.resolution ? oc.resden : nothing,
                        oc.bwfloor_sv, oc.bwfloor_tl;
                        outclass = oc.outclass)

        total_mass += d.deposited_mass;  total_raw += d.deposited_raw
        tot_dep    += d.n_deposited;     tot_oxy   += d.n_outside_xy
        tot_oxz    += d.n_outside_z
        tot_bg     += d.n_below_ground;  mass_bg   += d.mass_below_ground
        tot_bf     += d.n_below_first;   mass_bf   += d.mass_below_first
    end

    if _DEBUG_CONC
        println("DEBUG conccalc (Niout=", Niout, " Njout=", Njout, ")")
        println("  deposited mass        = ", total_mass)
        println("  deposited raw (m/rho) = ", total_raw)
        println("  parcels deposited     = ", tot_dep)
        println("  parcels outside xy    = ", tot_oxy)
        println("  parcels outside z     = ", tot_oxz)
        println("  below ground          = ", tot_bg, " ; mass = ", mass_bg)
        println("  below first layer     = ", tot_bf, " ; mass = ", mass_bf)
        if tot_dep > 0
            maxval, maxidx = findmax(view(grid, :, :, :, 1, 1))
            println("  concarray max = ", maxval, " at ", maxidx)
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# PROGNOSTIC parcel air mass — the dilution READ from the LPDM's own dispersion
# ─────────────────────────────────────────────────────────────────────────────
"""
    update_airm!(state, gconfig, m3, indt, nn; k=12, boxes=nothing) → n_updated

Set each parcel's represented air mass `state.airm` from the LOCAL parcel
spacing, so the mixing ratio `q = mass/airm` used by the ice closure and by
`plume_radiation` DILUTES as the parcel cloud spreads.

Why this and not an entrainment rate: in an LPDM the dispersion is ALREADY
carried by the stochastic trajectories (turbulence + shear + mesoscale).
Prescribing an entrainment `dV` on top would count the same process twice.
Here nothing is prescribed — the spreading is only *measured*, with the same
k-nearest-neighbour estimator the adaptive output kernel uses, so the
microphysics sees exactly the concentration the model reports.

    V_par = 4π√3 · h_h² · h_v / k          airm = ρ_met · V_par

`h_h²` = mean squared HORIZONTAL distance and `h_v²` = mean squared VERTICAL
offset of the same k neighbours.  The vertical scale is essential: the plume
is a pancake (100s of km wide, ~1 km thick) and an isotropic `h³` would
inflate `V` by orders of magnitude.  Derivation (locally uniform slab): the k
neighbours fill a cylinder of radius R, half-height H with ⟨r_h²⟩ = R²/2 and
⟨Δz²⟩ = H²/3, hence k = n·πR²·2H = n·4π√3·h_h²·h_v.

No tracer mass is destroyed — only the air it is diluted into grows.
`boxes` (a Set of relbox indices) restricts the update to the plume
populations; the domain-fill parcels keep their release air mass.
"""
function update_airm!(state::State.ParcelState, gconfig, m3, indt::Real, nn::Int;
                      k::Int = 12, boxes = nothing,
                      hv_out::Union{Nothing,Vector{Float32}} = nothing)
    dxm, dym = _out_spacing_m(gconfig)
    nn8 = Int8(nn)
    indtf = Float32(indt)
    N = length(state.indx)
    x = Vector{Float64}(undef, N); y = Vector{Float64}(undef, N)
    z = Vector{Float64}(undef, N); rho = Vector{Float64}(undef, N)
    act = falses(N)
    xmn = Inf; xmx = -Inf; ymn = Inf; ymx = -Inf
    nact = 0
    @inbounds for i in 1:N
        (state.status[i] == 0x01 && state.nest[i] == nn8) || continue
        boxes === nothing || Int(state.relbox[i]) in boxes || continue
        act[i] = true
        x[i] = Float64(state.indx[i]) * dxm
        y[i] = Float64(state.indy[i]) * dym
        z[i] = Float64(interpolation4Dlevsingle(m3.zlev, state.indx[i], state.indy[i],
                                                state.indl[i], indtf))
        rho[i] = max(Float64(interpolation4Dsingle(m3.rho, state.indx[i], state.indy[i],
                                                   state.indl[i], indtf)), 1e-6)
        xmn = min(xmn, x[i]); xmx = max(xmx, x[i])
        ymn = min(ymn, y[i]); ymx = max(ymx, y[i])
        nact += 1
    end
    nact <= k && return 0
    xmin = xmn; ymin = ymn; xmax = xmx; ymax = ymx
    L0 = max(sqrt(max(xmx - xmn, 1.0) * max(ymx - ymn, 1.0) / nact) * sqrt(k), 1.0)

    ncx = max(1, ceil(Int, (xmax - xmin) / L0) + 1)
    ncy = max(1, ceil(Int, (ymax - ymin) / L0) + 1)
    cid = zeros(Int, N)
    starts = zeros(Int, ncx * ncy + 1)
    @inbounds for i in 1:N
        act[i] || continue
        ix = clamp(floor(Int, (x[i] - xmin) / L0) + 1, 1, ncx)
        iy = clamp(floor(Int, (y[i] - ymin) / L0) + 1, 1, ncy)
        c = ix + (iy - 1) * ncx
        cid[i] = c; starts[c + 1] += 1
    end
    starts[1] = 1
    @inbounds for c in 1:(ncx * ncy); starts[c + 1] += starts[c]; end
    items = Vector{Int}(undef, nact)
    cursor = copy(starts)
    @inbounds for i in 1:N
        cid[i] == 0 && continue
        items[cursor[cid[i]]] = i; cursor[cid[i]] += 1
    end

    # k smallest horizontal distances, carrying the matching vertical offset.
    # Same stride subsample as the bandwidth pass: a point release puts ~all
    # parcels in one cell and a plain block scan would go O(N²).
    CAPC = 6 * k
    nthr = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads() + 4
    bufr = [Vector{Float64}(undef, k) for _ in 1:nthr]
    bufz = [Vector{Float64}(undef, k) for _ in 1:nthr]
    nup = Threads.Atomic{Int}(0)
    Threads.@threads :static for i in 1:N
        act[i] || continue
        br = bufr[Threads.threadid()]; bz = bufz[Threads.threadid()]
        nf = 0
        cx = clamp(floor(Int, (x[i] - xmin) / L0) + 1, 1, ncx)
        cy = clamp(floor(Int, (y[i] - ymin) / L0) + 1, 1, ncy)
        ci = _xcos(gconfig, state.indy[i])
        ex = ci >= 0.999 ? 1 : min(ceil(Int, 1.0 / ci), 8)
        occmax = 0
        @inbounds for oy in -1:1, ox in -ex:ex
            gx = cx + ox; gy = cy + oy
            (1 <= gx <= ncx && 1 <= gy <= ncy) || continue
            c = gx + (gy - 1) * ncx
            o = starts[c + 1] - starts[c]
            o > occmax && (occmax = o)
        end
        stride = occmax > CAPC ? cld(occmax, CAPC) : 1
        off = stride > 1 ? Int(mod(hash(i), stride)) : 0
        @inbounds for oy in -1:1, ox in -ex:ex
            gx = cx + ox; gy = cy + oy
            (1 <= gx <= ncx && 1 <= gy <= ncy) || continue
            c = gx + (gy - 1) * ncx
            for m in (starts[c] + off):stride:(starts[c + 1] - 1)
                j = items[m]
                j == i && continue
                dxp = (x[j] - x[i]) * ci; dyp = y[j] - y[i]
                r2 = dxp * dxp + dyp * dyp
                dz2 = (z[j] - z[i])^2
                if nf < k
                    p = nf
                    while p > 0 && br[p] > r2
                        br[p + 1] = br[p]; bz[p + 1] = bz[p]; p -= 1
                    end
                    br[p + 1] = r2; bz[p + 1] = dz2; nf += 1
                elseif r2 < br[k]
                    p = k - 1
                    while p > 0 && br[p] > r2
                        br[p + 1] = br[p]; bz[p + 1] = bz[p]; p -= 1
                    end
                    br[p + 1] = r2; bz[p + 1] = dz2
                end
            end
        end
        nf > 0 || continue
        h2 = 0.0; v2 = 0.0
        @inbounds for m in 1:nf; h2 += br[m]; v2 += bz[m]; end
        h2 /= (nf * stride)          # rank correction (stride=1 → exact)
        v2 /= nf                     # vertical spread is a mean, not a rank
        hh = sqrt(max(h2, 1.0)); hv = sqrt(max(v2, 1.0))
        Vpar = 4.0 * pi * sqrt(3.0) * hh * hh * hv / k
        # never below the seeded umbrella allotment (state.airm0; 0 = no floor)
        state.airm[i] = eltype(state.airm)(max(rho[i] * Vpar,
                                               Float64(state.airm0[i])))
        # local vertical spread → settling escape scale for the ice closure
        hv_out === nothing || (hv_out[i] = Float32(hv))
        Threads.atomic_add!(nup, 1)
    end
    return nup[]
end
