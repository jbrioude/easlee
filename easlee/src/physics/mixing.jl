# physics/mixing.jl
# Pairwise mass exchange between overlapping parcels — ported from easlee_2D.
#
# Method:
#   1. Build a 3-D cell-list spatial hash of active parcel positions (x, y, z) [m].
#   2. For each pair (i,j) within nsig·σ of each other, compute
#        λ_ij = α·K_h / (π·σ_h²) · exp(-r_h²/2σ_h² - r_v²/2σ_v²)
#      where σ_v is the *per-parcel* Lagrangian turbulent length scale yl
#      stored in state.miss3.  σ_h is a fixed horizontal scale (≈ grid spacing).
#   3. Renormalise λ so that Σ_j λ_ij·r²_h,ij = 4·K_h  (2-D horizontal second-moment).
#   4. Flux: Fij = λ_ij·(c_j − c_i)  where c [ppb] is already a mixing ratio;
#      apply θ-limiter for positivity.
#   5. Neighbor cap (max_neighbors per parcel) prevents O(N²) memory when
#      parcels are spatially concentrated near a release point.
#
# Public API:
#   MixingConfig                  — parameters
#   apply_mixing_xyz!(q, x, y, z, ...) — one mixing sub-step for all parcels

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

"""
    MixingConfig

Parameters for the pairwise mass-exchange step.

Fields
──────
- `sigma0_h`      horizontal Gaussian kernel half-width [m] (typ. ≈ grid spacing)
- `K_h`           horizontal diffusivity [m²/s]
- `alpha`         overall rate prefactor (default 1.0)
- `nsig`          pair-search radius in units of σ₀_h (default 3.0)
- `max_neighbors` max pairs per parcel — prevents O(N²) for concentrated releases
                  (default 30; increase for sparser, more uniform distributions)

Vertical mixing length scale is taken from `state.miss3[i]` (the Lagrangian
turbulent length scale `yl` [m] written there by the advance! step).
Pairs with no turbulence history (yl = 0) fall back to a 10 m floor.
"""
struct MixingConfig
    sigma0_h      :: Float64
    K_h           :: Float64
    alpha         :: Float64
    nsig          :: Float64
    max_neighbors :: Int
    # exchange = :equal (current behaviour, equal air-mass quanta) or
    # :massweighted — the reduced-mass DIRECTED exchange
    #   dq_i = Δt·λ_ij·[2m_j/(m_i+m_j)]·(q_j−q_i)   (mirrored for j)
    # required for two-population runs (plume + background): species mass is
    # conserved exactly for ANY mass ratio, uniform q stays an exact fixed
    # point, and the equal-mass limit is identical to :equal (regression
    # guarantee).  The asymmetry IS the entrainment physics: small plume
    # parcels dilute fast, large background parcels barely move.
    exchange      :: Symbol
end

function MixingConfig(;
    sigma0_h      :: Float64 = 2500.0,
    K_h           :: Float64 = 10.0,
    alpha         :: Float64 = 1.0,
    nsig          :: Float64 = 3.0,
    max_neighbors :: Int     = 30,
    exchange      :: Symbol  = :equal,
)
    exchange in (:equal, :massweighted) ||
        error("MixingConfig: exchange must be :equal or :massweighted")
    MixingConfig(sigma0_h, K_h, alpha, nsig, max_neighbors, exchange)
end

# ─────────────────────────────────────────────────────────────────────────────
# Grid-spacing helpers (metres per index unit)
# ─────────────────────────────────────────────────────────────────────────────

function _grid_spacing_m(g::Types.Lambert, ::State.ParcelState)
    return Float64(g.dx), Float64(g.dy)
end
function _grid_spacing_m(g::Types.Mercator, ::State.ParcelState)
    return Float64(g.dx), Float64(g.dy)   # projection metres per index (Meso-NH etc.)
end
function _grid_spacing_m(g::Types.Rectilinear, ::State.ParcelState)
    dx_eq   = 1.0 / Float64(g.dxconst)
    dy_m    = 1.0 / Float64(g.dyconst)
    mid_j   = 0.5 * Float64(g.Nj)
    lat_mid = Float64(g.offsety) + Float64(g.coeffy) * mid_j
    return dx_eq * cosd(lat_mid), dy_m
end

# ─────────────────────────────────────────────────────────────────────────────
# Cell-list hash (horizontal cells only — vertical handled by yl filter)
# ─────────────────────────────────────────────────────────────────────────────
#
# CSR (compressed-sparse-row) layout: `items` holds all active parcel indices
# grouped by cell; cell c spans items[starts[c] : starts[c+1]-1].  Parcels
# within a cell are in ascending index order (stable counting sort), matching
# the insertion order of the previous Vector-of-Vectors implementation.
#
# A concrete struct (no closures) keeps everything type-stable: an earlier
# version returned a `cell_id` closure capturing reassigned variables, which
# boxed them and made every downstream arithmetic op dynamically dispatched
# (~10× slowdown and ~400 MB of allocations per pairlist build at N = 50k).

struct CellList2D
    ncx    :: Int
    ncy    :: Int
    x_min  :: Float64
    y_min  :: Float64
    starts :: Vector{Int}   # length ncx*ncy + 1 — CSR offsets
    items  :: Vector{Int}   # active parcel indices grouped by cell
end

function _build_cell_list_2d(x::AbstractVector{Float64},
                              y::AbstractVector{Float64},
                              active::BitVector,
                              cellsize::Float64)
    N = length(x)
    x_min = Inf;  x_max = -Inf;  y_min = Inf;  y_max = -Inf
    @inbounds for i in 1:N
        active[i] || continue
        x_min = min(x_min, x[i]);  x_max = max(x_max, x[i])
        y_min = min(y_min, y[i]);  y_max = max(y_max, y[i])
    end
    ncx = max(1, ceil(Int, (x_max - x_min + 2*cellsize) / cellsize))
    ncy = max(1, ceil(Int, (y_max - y_min + 2*cellsize) / cellsize))
    ncell = ncx * ncy

    # Pass 1: cell id per parcel + per-cell counts
    cid    = Vector{Int}(undef, N)
    starts = zeros(Int, ncell + 1)
    n_act  = 0
    @inbounds for i in 1:N
        if active[i]
            ix = clamp(floor(Int, (x[i] - x_min) / cellsize) + 1, 1, ncx)
            iy = clamp(floor(Int, (y[i] - y_min) / cellsize) + 1, 1, ncy)
            c  = ix + (iy-1)*ncx
            cid[i] = c
            starts[c+1] += 1
            n_act += 1
        else
            cid[i] = 0
        end
    end

    # Prefix sum → CSR offsets
    starts[1] = 1
    @inbounds for c in 1:ncell
        starts[c+1] += starts[c]
    end

    # Pass 2: stable fill (ascending parcel index within each cell)
    items  = Vector{Int}(undef, n_act)
    cursor = copy(starts)
    @inbounds for i in 1:N
        c = cid[i]
        c == 0 && continue
        items[cursor[c]] = i
        cursor[c] += 1
    end

    return CellList2D(ncx, ncy, x_min, y_min, starts, items)
end

# ─────────────────────────────────────────────────────────────────────────────
# Packed-cluster pre-conditioning
# ─────────────────────────────────────────────────────────────────────────────
#
# When many trajectories are collocated (distance < σ₀_h/pack_factor), the
# pairwise kernel gives them weight ≈ 1 and mixing is effectively instantaneous.
# Rather than let them explode the pair-list, we:
#   1. Detect clusters via a coarse cell-list (cell size = σ₀_h/pack_factor).
#   2. Average the mixing ratios across each cluster (intensive → mean).
#   3. Redistribute the k cluster members on a regular ceil(√k)×ceil(√k) grid
#      centred on the cluster centroid, with spacing σ₀_h.
#      → every parcel then has 4–8 well-defined grid neighbours; no cap issues.
#
# This is called automatically by apply_mixing_xyz! and
# modifies x, y, chem in-place.  z is left unchanged (vertical handled by
# the Langevin scheme).  Returns the number of parcels that were redistributed.

const _DEFAULT_PACK_FACTOR = 4   # cluster cell = σ₀_h / 4

# AbstractMatrix{<:Real} so this works whether chem is Float32 or Float64.
function _merge_packed_clusters!(
    x      :: Vector{Float64},
    y      :: Vector{Float64},
    chem   :: AbstractMatrix{<:Real},   # (N × S) — modified in place
    active :: BitVector,
    mcfg   :: MixingConfig,
    pack_factor :: Int = _DEFAULT_PACK_FACTOR,
)
    cell_size = mcfg.sigma0_h / pack_factor

    cl = _build_cell_list_2d(x, y, active, cell_size)

    n_species   = size(chem, 2)
    n_relocated = 0

    @inbounds for cid in 1:(cl.ncx * cl.ncy)
        lo = cl.starts[cid];  hi = cl.starts[cid+1] - 1
        k  = hi - lo + 1
        k <= 1 && continue              # singleton — nothing to do
        members = view(cl.items, lo:hi)

        # ── 1. Average chemistry (mixing ratio = intensive → plain mean) ──────
        #
        # c is a mixing ratio (ppb): intensive variable → mean is the correct
        # operation (no mass weighting needed, all parcels carry the same
        # air-mass element by construction).
        c_mean = zeros(Float64, n_species)
        for i in members
            for s in 1:n_species
                c_mean[s] += Float64(chem[i, s])
            end
        end
        inv_k = 1.0 / k
        for s in 1:n_species
            c_mean[s] *= inv_k
        end
        for i in members
            for s in 1:n_species
                chem[i, s] = c_mean[s]
            end
        end

        # ── 2. Cluster centroid ───────────────────────────────────────────────
        cx = 0.0;  cy = 0.0
        for i in members
            cx += x[i];  cy += y[i]
        end
        cx *= inv_k;  cy *= inv_k

        # ── 3. Redistribute on ceil(√k)×ceil(√k) regular grid ────────────────
        #
        # Spacing = σ₀_h so that the subsequent pairwise step sees each parcel
        # with 4–8 well-separated, well-defined neighbours.  This turns the
        # O(N²) near-collocated case into O(k) pairs and eliminates loop-order
        # bias from the cap.
        nx  = ceil(Int, sqrt(k))
        sp  = mcfg.sigma0_h                    # grid spacing = kernel width
        x0  = cx - 0.5 * (nx - 1) * sp        # grid origin (left column)
        y0  = cy - 0.5 * (nx - 1) * sp        # grid origin (bottom row)

        for (idx, i) in enumerate(members)
            row = (idx - 1) ÷ nx
            col = (idx - 1) % nx
            x[i] = x0 + col * sp
            y[i] = y0 + row * sp
        end

        n_relocated += k
    end

    return n_relocated
end

# ─────────────────────────────────────────────────────────────────────────────
# Coverage-based simple background relaxation  (q_bg vector, no CAMS needed)
# ─────────────────────────────────────────────────────────────────────────────
#
# For each active parcel i compute the coverage fraction from the pre-renorm
# second moment:
#   cov_i = min(mom2[i] × σ₀² / (4 α K_h), 1)
# An interior parcel fully embedded in the plume has cov ≈ 1 → no relaxation.
# An isolated / edge parcel has cov < 1 → relaxes toward q_bg at rate β_i:
#   β_i = (1 - cov_i) × K_h / σ₀²                           [1/s]
#   q_i += β_i · dt · (q_bg - q_i)
#
# Called by apply_mixing_xyz! / apply_mixing_xyz_substeps! when q_bg ≠ nothing.

function _apply_simple_background!(
    chem  :: Matrix{Float64},
    active:: BitVector,
    mom2  :: Vector{Float64},
    q_bg  :: Vector{Float64},
    mcfg  :: MixingConfig,
    dt    :: Float64,
)
    ideal_mom2  = 4.0 * mcfg.alpha * mcfg.K_h / mcfg.sigma0_h^2
    K_over_sig2 = mcfg.K_h / mcfg.sigma0_h^2
    n_bg        = length(q_bg)

    @inbounds for i in eachindex(active)
        active[i] || continue
        cov_i = clamp(mom2[i] / ideal_mom2, 0.0, 1.0)
        β     = (1.0 - cov_i) * K_over_sig2
        β * dt < 1e-12 && continue
        Δ = β * dt
        for s in 1:n_bg
            ci         = chem[i, s]
            chem[i, s] = max(ci + Δ * (q_bg[s] - ci), 0.0)
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Pairlist (cell-list + Gaussian filter + neighbor cap + 2nd-moment renorm)
# ─────────────────────────────────────────────────────────────────────────────

function _build_pairlist(
    x      :: Vector{Float64},
    y      :: Vector{Float64},
    z      :: Vector{Float64},
    yl     :: Vector{Float64},   # per-parcel vertical length scale [m]
    active :: BitVector,
    mcfg   :: MixingConfig,
)
    N         = length(x)
    rmax_h    = mcfg.nsig * mcfg.sigma0_h
    sig2_h    = mcfg.sigma0_h^2

    cl = _build_cell_list_2d(x, y, active, rmax_h)
    ncx    = cl.ncx
    ncy    = cl.ncy
    starts = cl.starts
    items  = cl.items
    n_act  = length(items)

    # Hoisted config fields (avoid per-iteration struct loads in the hot loop)
    max_nb   = mcfg.max_neighbors
    nsig2    = mcfg.nsig^2
    lam_pref = mcfg.alpha * mcfg.K_h / (π * sig2_h^2)

    pairs_i  = Int[];     sizehint!(pairs_i,  16 * n_act)
    pairs_j  = Int[];     sizehint!(pairs_j,  16 * n_act)
    lambdas  = Float64[]; sizehint!(lambdas,  16 * n_act)
    n_nbrs   = zeros(Int, N)   # neighbor count per parcel (for the cap)

    @inbounds for iy in 1:ncy, ix in 1:ncx
        cid0 = ix + (iy-1)*ncx
        p_lo = starts[cid0]
        p_hi = starts[cid0+1] - 1
        p_hi < p_lo && continue

        for jy in max(1,iy-1):min(ncy,iy+1),
            jx in max(1,ix-1):min(ncx,ix+1)

            cid1  = jx + (jy-1)*ncx
            # Skip lower-triangle of cell pairs: each (cid0,cid1) pair with
            # cid1 < cid0 would duplicate a (cid1,cid0) pair already processed.
            # Same-cell pairs (cid1==cid0) are deduplicated by bstart=ap+1.
            cid1 < cid0 && continue
            q_lo = starts[cid1]
            q_hi = starts[cid1+1] - 1

            for ap in p_lo:p_hi
                i = items[ap]
                # Fast-path: if i is already at cap (filled from being a j-target
                # in an earlier ap iteration), skip this outer iteration entirely.
                n_nbrs[i] >= max_nb && continue

                bstart = (cid0 == cid1) ? ap + 1 : q_lo
                xi = x[i];  yi = y[i];  zi = z[i];  yli = yl[i]

                for bp in bstart:q_hi
                    # i-side cap: break (not continue) so we stop scanning
                    # further j's entirely.  Without this the inner loop runs
                    # O(N) even after i is full — O(N²) for a point source.
                    n_nbrs[i] >= max_nb && break

                    j = items[bp]

                    # j-side cap check
                    n_nbrs[j] >= max_nb && continue

                    dxij = x[j] - xi
                    dyij = y[j] - yi

                    # horizontal Mahalanobis check first (cheap)
                    mah2_h = (dxij^2 + dyij^2) / sig2_h
                    mah2_h > nsig2 && continue

                    # vertical: use arithmetic mean of per-parcel yl
                    dzij   = z[j] - zi
                    sig_v  = 0.5 * (yli + yl[j])
                    mah2_v = sig_v > 1e-3 ? dzij^2 / sig_v^2 : 0.0

                    wij = exp(-0.5 * (mah2_h + mah2_v))

                    push!(pairs_i, i);  push!(pairs_j, j)
                    push!(lambdas, lam_pref * wij)
                    n_nbrs[i] += 1;     n_nbrs[j] += 1
                end
            end
        end
    end

    # Second-moment renormalisation: Σ_j λ_ij·r²_h,ij = 4·K_h  (2-D horizontal)
    #
    # c [ppb] is an intensive mixing ratio.  The 2-D Taylor expansion gives
    #   dc_i/dt = Σ_j λ_ij (c_j−c_i) ≈ ½ Σ_j λ_ij (Δx²+Δy²) ∇²_h c = K_h ∇²_h c
    # when Σ_j λ_ij r²_h,ij = 4·K_h (factor 4 = 2 dimensions × 2 for Taylor prefactor).
    # No parcel-mass factor appears because c is already intensive.
    npairs = length(lambdas)
    mom2   = zeros(Float64, N)
    @inbounds for p in 1:npairs
        i = pairs_i[p];  j = pairs_j[p]
        r2_h = (x[j]-x[i])^2 + (y[j]-y[i])^2   # horizontal distance only
        mom2[i] += lambdas[p] * r2_h
        mom2[j] += lambdas[p] * r2_h
    end
    @inbounds for p in 1:npairs
        i  = pairs_i[p];  j = pairs_j[p]
        si = 4.0 * mcfg.K_h / max(mom2[i], 1e-30)
        sj = 4.0 * mcfg.K_h / max(mom2[j], 1e-30)
        lambdas[p] *= 0.5 * (si + sj)
    end

    # Return mom2 (pre-renorm second moment, units 1/s) alongside pairs.
    # Callers use it to compute per-parcel coverage fraction:
    #   cov_i = min(mom2[i] × σ₀² / (4 α K_h), 1)
    # An interior parcel fully surrounded by plume neighbours has mom2[i] ≈ 4αK_h/σ₀².
    # An edge / isolated parcel has mom2[i] < 4αK_h/σ₀² → cov_i < 1 → background mixing.
    return pairs_i, pairs_j, lambdas, mom2
end

# ─────────────────────────────────────────────────────────────────────────────
# CSR adjacency + gather exchange
# ─────────────────────────────────────────────────────────────────────────────
#
# The exchange used to be a scatter over pairs (dc[i]+=F, dc[j]-=F) with one
# N×S accumulation buffer per thread, summed at the end.  The gather form
# stores each pair twice in a CSR adjacency (i→j and j→i, same λ) and computes
#   dc[i,s] = dt · Σ_k λ_k (q[nbr_k, s] − q[i, s])
# with one thread per parcel — no write conflicts, no per-thread buffers.
# Mathematically identical to the scatter (each pair contributes λ(q_j−q_i)
# to i and the negative to j); only the fp summation order differs.

"""
    _pairs_to_adjacency(N, pairs_i, pairs_j, lambdas)
        → (adj_start, adj_nbr, adj_lam)

Build a symmetric CSR adjacency from a pair list.  Parcel i's neighbours are
`adj_nbr[adj_start[i] : adj_start[i+1]-1]` with rates in `adj_lam`.
"""
function _pairs_to_adjacency(N       :: Int,
                             pairs_i :: Vector{Int},
                             pairs_j :: Vector{Int},
                             lambdas :: Vector{Float64};
                             masses  :: Union{Nothing,AbstractVector} = nothing)
    npairs = length(lambdas)
    adj_start = zeros(Int, N + 1)
    @inbounds for p in 1:npairs
        adj_start[pairs_i[p] + 1] += 1
        adj_start[pairs_j[p] + 1] += 1
    end
    adj_start[1] = 1
    @inbounds for i in 1:N
        adj_start[i+1] += adj_start[i]
    end

    adj_nbr = Vector{Int}(undef, 2*npairs)
    adj_lam = Vector{Float64}(undef, 2*npairs)
    cursor  = copy(adj_start)
    @inbounds for p in 1:npairs
        i = pairs_i[p];  j = pairs_j[p];  λ = lambdas[p]
        # mass-weighted DIRECTED exchange (exchange = :massweighted): each side
        # stores its own reduced-mass factor 2m_other/(m_i+m_j) ∈ (0,2) —
        # conserves Σ m·q exactly for any mass ratio; equal masses → both
        # factors 1 → identical to the unweighted path.
        λi = λ; λj = λ
        if masses !== nothing
            mi = Float64(masses[i]); mj = Float64(masses[j])
            s = mi + mj
            if s > 0.0
                λi = λ * (2.0 * mj / s)
                λj = λ * (2.0 * mi / s)
            end
        end
        ci = cursor[i];  adj_nbr[ci] = j;  adj_lam[ci] = λi;  cursor[i] = ci + 1
        cj = cursor[j];  adj_nbr[cj] = i;  adj_lam[cj] = λj;  cursor[j] = cj + 1
    end
    return adj_start, adj_nbr, adj_lam
end

"""
    _gather_exchange_dc!(dc, chem, adj_start, adj_nbr, adj_lam, dt)

Fill `dc[i,s] = dt · Σ_k λ_k (chem[nbr_k, s] − chem[i, s])` — one thread per
parcel, race-free.  `chem` may be Float32 or Float64; accumulation is Float64.
"""
function _gather_exchange_dc!(dc        :: Matrix{Float64},
                              chem      :: AbstractMatrix{<:Real},
                              adj_start :: Vector{Int},
                              adj_nbr   :: Vector{Int},
                              adj_lam   :: Vector{Float64},
                              dt        :: Float64)
    N, n_species = size(chem)
    Threads.@threads :static for i in 1:N
        @inbounds begin
            lo = adj_start[i];  hi = adj_start[i+1] - 1
            for s in 1:n_species
                ci  = Float64(chem[i, s])
                acc = 0.0
                for ap in lo:hi
                    acc += adj_lam[ap] * (Float64(chem[adj_nbr[ap], s]) - ci)
                end
                dc[i, s] = dt * acc
            end
        end
    end
    return dc
end

# ─────────────────────────────────────────────────────────────────────────────
# apply_mixing_xyz!  — position-direct entry point (testing / offline use)
# ─────────────────────────────────────────────────────────────────────────────

"""
    apply_mixing_xyz!(chem, x, y, z, yl, active, dt, mcfg;
                      q_bg=nothing, merge_packed=true, masses=nothing) → npairs

Pairwise exchange on pre-computed physical positions
[m] directly, bypassing the met-field interpolation step.  Useful for offline
validation tests that run without AROME data.

Arguments
──────────
- `chem`        (N × S) matrix [ppb or any consistent mixing-ratio unit]
- `x,y,z`       parcel positions [m] (inactive entries may hold garbage)
- `yl`          per-parcel vertical length scale [m] (10 m floor applied)
- `active`      BitVector marking the parcels that take part
- `dt`          sub-step [s]
- `mcfg`        MixingConfig

Optional keyword arguments
──────────────────────────
- `q_bg`         constant ambient tracer vector (length S) toward which
                 non-fully-surrounded parcels relax (coverage-based edge
                 entrainment); `nothing` = no background exchange.
- `merge_packed` merge sub-σ₀/4 clusters before pairing (default true).
- `masses`       per-parcel air masses for the `:massweighted` exchange.
"""
function apply_mixing_xyz!(chem        :: Matrix{Float64},
                            x           :: Vector{Float64},
                            y           :: Vector{Float64},
                            z           :: Vector{Float64},
                            yl          :: Vector{Float64},
                            active      :: BitVector,
                            dt          :: Float64,
                            mcfg        :: MixingConfig;
                            q_bg        :: Union{Vector{Float64}, Nothing}  = nothing,
                            merge_packed :: Bool                            = true,
                            masses      :: Union{Vector{Float64}, Nothing}  = nothing)
    count(active) < 2 && return 0
    N = length(x)

    # Merge packed clusters: averages chemistry in-place and redistributes
    # x, y onto a σ₀_h-spaced grid for any cluster (distance < σ₀_h/4).
    # Disable with merge_packed=false for fine-scale regular-grid simulations
    # where the cluster threshold σ₀_h/4 is comparable to the parcel spacing —
    # the merger would otherwise relocate parcels off the output grid.
    merge_packed && _merge_packed_clusters!(x, y, chem, active, mcfg)

    yl_fl = [active[i] ? max(yl[i], 10.0) : 10.0 for i in 1:N]
    pairs_i, pairs_j, lambdas, mom2 = _build_pairlist(x, y, z, yl_fl, active, mcfg)

    npairs = length(lambdas)
    npairs == 0 && return 0

    # Gather exchange (see _gather_exchange_dc! — race-free, one thread/parcel)
    n_species = size(chem, 2)
    adj_start, adj_nbr, adj_lam = _pairs_to_adjacency(N, pairs_i, pairs_j, lambdas;
        masses = mcfg.exchange === :massweighted ? masses : nothing)
    dc = Matrix{Float64}(undef, N, n_species)
    _gather_exchange_dc!(dc, chem, adj_start, adj_nbr, adj_lam, dt)

    # Per-species θ-limiter.
    theta = ones(Float64, n_species)
    @inbounds for i in 1:N
        active[i] || continue
        for s in 1:n_species
            dci = dc[i, s]
            dci >= 0.0 && continue
            θs = chem[i, s] / max(-dci, 1e-100)
            if θs < theta[s]
                theta[s] = θs
            end
        end
    end

    @inbounds for i in 1:N
        active[i] || continue
        for s in 1:n_species
            chem[i, s] = max(chem[i, s] + theta[s] * dc[i, s], 0.0)
        end
    end

    # ── Simple background relaxation (constant q_bg, no CAMS) ────────────────
    if !isnothing(q_bg)
        _apply_simple_background!(chem, active, mom2, q_bg, mcfg, dt)
    end

    return npairs
end

# ─────────────────────────────────────────────────────────────────────────────
# MixingPairlist — cached pairlist for sub-step reuse (Verlet-list approach)
# ─────────────────────────────────────────────────────────────────────────────

"""
    MixingPairlist

Pre-built pairlist for reuse across multiple mixing sub-steps.

The list is built once, then reused as long as no parcel has moved more than
`skin_h` from its position at build time.  The search radius at build time is
extended by `2·skin_h` beyond the normal `nsig·σ₀_h`, so that pairs that drift
into range between rebuilds are never missed.

Fields
──────
- `pairs_i`, `pairs_j`  parcel index pairs
- `lambdas`             renormalised exchange rates [1/s]
- `mom2`                pre-renorm 2nd-moment per parcel [m²/s] (background mixing)
- `adj_start`, `adj_nbr`, `adj_lam`
                        symmetric CSR adjacency (each pair stored both ways) —
                        used by the gather exchange, built once per pairlist
- `x0`, `y0`, `z0`     positions at build time (for staleness check)
- `skin_h`              max safe horizontal displacement before rebuild [m]
                        (= 0 → staleness check disabled, list always valid)
"""
struct MixingPairlist
    pairs_i   :: Vector{Int}
    pairs_j   :: Vector{Int}
    lambdas   :: Vector{Float64}
    mom2      :: Vector{Float64}
    adj_start :: Vector{Int}
    adj_nbr   :: Vector{Int}
    adj_lam   :: Vector{Float64}
    x0        :: Vector{Float64}
    y0        :: Vector{Float64}
    z0        :: Vector{Float64}
    skin_h    :: Float64
end

"""
    build_pairlist_xyz(x, y, z, yl, active, mcfg; skin_frac=0.0) → MixingPairlist

Build a pairlist suitable for reuse over multiple sub-steps.

`skin_frac ≥ 0`: extend the search radius by `skin_frac · σ₀_h` beyond the
normal `nsig · σ₀_h`.  The resulting `skin_h = skin_frac · σ₀_h / 2` is the
maximum displacement any parcel may undergo before the list must be rebuilt.

- `skin_frac = 0.0` — no extension; use with static positions (test/benchmark).
- `skin_frac = 0.5` — 25 % radius extension; `skin_h = σ₀_h/4`, safe for
  ~10 exchange sub-steps at typical wind speeds.
"""
function build_pairlist_xyz(
    x         :: Vector{Float64},
    y         :: Vector{Float64},
    z         :: Vector{Float64},
    yl        :: Vector{Float64},
    active    :: BitVector,
    mcfg      :: MixingConfig;
    skin_frac :: Float64 = 0.0,
)
    mcfg_ext = skin_frac > 0.0 ?
        MixingConfig(sigma0_h      = mcfg.sigma0_h,
                     K_h           = mcfg.K_h,
                     alpha         = mcfg.alpha,
                     nsig          = mcfg.nsig + skin_frac,
                     max_neighbors = mcfg.max_neighbors) :
        mcfg

    pairs_i, pairs_j, lambdas, mom2 =
        _build_pairlist(x, y, z, yl, active, mcfg_ext)

    adj_start, adj_nbr, adj_lam =
        _pairs_to_adjacency(length(x), pairs_i, pairs_j, lambdas)

    skin_h = skin_frac * mcfg.sigma0_h / 2.0

    return MixingPairlist(pairs_i, pairs_j, lambdas, mom2,
                          adj_start, adj_nbr, adj_lam,
                          copy(x), copy(y), copy(z), skin_h)
end

"""
    needs_rebuild(pl, x, y, z, active) → Bool

Return `true` if any active parcel has moved more than `pl.skin_h` in the
horizontal since the list was built.  Always `false` when `pl.skin_h == 0`.
"""
function needs_rebuild(pl     :: MixingPairlist,
                       x      :: Vector{Float64},
                       y      :: Vector{Float64},
                       z      :: Vector{Float64},
                       active :: BitVector)
    pl.skin_h == 0.0 && return false
    skin2 = pl.skin_h^2
    @inbounds for i in eachindex(x)
        active[i] || continue
        dx = x[i] - pl.x0[i];  dy = y[i] - pl.y0[i]
        dx*dx + dy*dy > skin2 && return true
    end
    return false
end

# ── Inner exchange + θ-limiter (shared by cached and non-cached paths) ────────

function _exchange_step!(
    chem   :: Matrix{Float64},
    active :: BitVector,
    pl     :: MixingPairlist,
    dt     :: Float64,
)
    N         = size(chem, 1)
    n_species = size(chem, 2)
    npairs    = length(pl.lambdas)
    npairs == 0 && return

    # Gather over the cached CSR adjacency (built once in build_pairlist_xyz)
    dc = Matrix{Float64}(undef, N, n_species)
    _gather_exchange_dc!(dc, chem, pl.adj_start, pl.adj_nbr, pl.adj_lam, dt)

    # Per-species θ-limiter (positivity without cross-species throttling)
    theta = ones(Float64, n_species)
    @inbounds for i in 1:N
        active[i] || continue
        for s in 1:n_species
            dci = dc[i, s]
            dci >= 0.0 && continue
            θs = chem[i, s] / max(-dci, 1e-100)
            θs < theta[s] && (theta[s] = θs)
        end
    end
    @inbounds for i in 1:N
        active[i] || continue
        for s in 1:n_species
            chem[i, s] = max(chem[i, s] + theta[s] * dc[i, s], 0.0)
        end
    end
end

"""
    apply_mixing_xyz_substeps!(chem, x, y, z, yl, active, dt, n_substeps, mcfg;
                                skin_frac=0.0, merge_packed=true, bg=nothing, ...)
        → (npairs, n_builds)

Build the pairlist **once**, then run `n_substeps` exchange sub-steps of size
`dt`.  Total simulated mixing time = `n_substeps × dt`.

Equivalent to calling `apply_mixing_xyz!` `n_substeps` times, but the pairlist
build cost (O(N · k) cell-list scan) is paid only when the list becomes stale.

With `skin_frac = 0` (default) the list is assumed permanently valid — correct
when parcel positions are not updated between sub-steps (pure-mixing test, or
when advection and mixing are operator-split at a coarser time scale than
`dt`).  Set `skin_frac > 0` to enable staleness detection via `needs_rebuild`.

Returns `(npairs, n_builds)`:
- `npairs`   number of active pairs at the first build
- `n_builds` actual number of pairlist rebuilds (1 if positions never change)
"""
function apply_mixing_xyz_substeps!(
    chem         :: Matrix{Float64},
    x            :: Vector{Float64},
    y            :: Vector{Float64},
    z            :: Vector{Float64},
    yl           :: Vector{Float64},
    active       :: BitVector,
    dt           :: Float64,
    n_substeps   :: Int,
    mcfg         :: MixingConfig;
    skin_frac    :: Float64                          = 0.0,
    use_gpu      :: Bool                             = false,
    merge_packed :: Bool                             = true,
    q_bg         :: Union{Vector{Float64}, Nothing}  = nothing,
)
    count(active) < 2 && return 0, 0
    N = length(x)

    merge_packed && _merge_packed_clusters!(x, y, chem, active, mcfg)

    yl_fl = [active[i] ? max(yl[i], 10.0) : 10.0 for i in 1:N]
    pl = build_pairlist_xyz(x, y, z, yl_fl, active, mcfg; skin_frac=skin_frac)

    npairs_0 = length(pl.lambdas)
    npairs_0 == 0 && return 0, 0

    n_builds = 1

    if false
    else
        # ── CPU path (threaded per-thread accumulation buffers) ────────────────
        for sub in 1:n_substeps
            # Staleness check: rebuild if positions have drifted outside skin.
            # (With skin_frac=0 this is a no-op — needs_rebuild always returns false.)
            if sub > 1 && needs_rebuild(pl, x, y, z, active)
                yl_fl2 = [active[i] ? max(yl[i], 10.0) : 10.0 for i in 1:N]
                pl = build_pairlist_xyz(x, y, z, yl_fl2, active, mcfg;
                                        skin_frac=skin_frac)
                n_builds += 1
            end
            _exchange_step!(chem, active, pl, dt)
        end
    end

    # Background mixing: relax over the full window dt*n_substeps.
    if !isnothing(q_bg)
        _apply_simple_background!(chem, active, pl.mom2, q_bg, mcfg,
                                  dt * n_substeps)
    end

    return npairs_0, n_builds
end

