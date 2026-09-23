
# ── Typed output-grid container ───────────────────────────────────────────────
#
# OutputConfig replaces the legacy outputconfig Dict{String,Any}.
# Array fields (concarray, massarray, area, volume) share the SAME underlying
# memory as the Dict entries created by initoutgrid — no copying, no extra
# allocation.  The struct is `mutable` so ncid and the accumulators can be
# updated in-place each output cycle.

mutable struct OutputConfig
    # Accumulators — written every step by conccalc!
    concarray::Array{Float32,5}    # (Niout, Njout, Nkout, Nc, Nage)
    massarray::Array{Float64,5}

    # Grid geometry — read-only after construction
    area::Matrix{Float32}          # horizontal cell area [m²]   (Niout × Njout)
    volume::Array{Float32,3}       # cell volume [m³]             (Niout × Njout × Nkout)
    outheight::Vector{Float32}     # output layer edge heights [m]
    outinterfaces::Vector{Float32} # layer boundaries used for z-binning (= outheight)
    dlon::Float32
    dlat::Float32
    lon0::Float32
    lat0::Float32
    Niout::Int32
    Njout::Int32
    Nkout::Int32
    Nc::Int32                      # number of OUTPUT COLUMNS (dim 4): boxes
                                   # ("release" grouping), species ("species"),
                                   # or 1 ("all" collapse)

    # Output semantics (set from COMMAND / RELEASES at build time)
    iout::Int8                     # 1 concentration, 2 mixing ratio, 3 both
    direction::Int8                # +1 forward, -1 backward (sensitivity output)
    boxspecies::Vector{Int32}      # FLEXPART species number per release box
    relmass::Vector{Float32}       # released mass per OUTPUT COLUMN [kg]
                                   # (backward norm; per box under "release")
    # ── Output grouping ([output].group / COMMAND OUTPUTGROUP) ───────────────
    # outclass maps a parcel's relbox to its accumulator column; classspecies
    # is the species number of each column (what the per-species writers group
    # by).  Under "release" grouping outclass is the identity and classspecies
    # == boxspecies — bit-identical legacy.  Grouping by species makes the
    # grid memory independent of release-box count (~18 MB/column on the HT
    # grid; a 3400-box file release OOM'd at 70 GB with per-box columns).
    outgroup::Symbol               # :species | :release | :all
    outclass::Vector{Int32}        # release box → output column (length = nb)
    classspecies::Vector{Int32}    # output column → species number (length = Nc)

    # ── Horizontal CRS of the output grid ────────────────────────────────────
    # proj_m = true when lon0/lat0/dlon/dlat are PROJECTION METRES, not degrees
    # (output.grid.crs = "native" on a Lambert nest — the same discriminator
    # initoutgrid uses to pick the cell-area formula).  In that case the axes
    # are written as CF x/y in metres and the true geographic coordinates go
    # out as the 2-D auxiliary variables lon2d/lat2d, sampled from the nest's
    # own LON/LAT.  Writing metres under a "degrees_east" label put the grid
    # at 0N 0E with 1300-degree spacing for every external reader.
    proj_m::Bool
    lon2d::Matrix{Float32}         # (Niout × Njout) true longitude, native CRS only
    lat2d::Matrix{Float32}         # (Niout × Njout) true latitude,  native CRS only

    # NetCDF file handle — set by initnetcdf!, read by writenetcdf
    ncid::Any                      # NcFile or nothing

    # Gridding kernel:
    #   :none     — nearest-cell deposition
    #   :uniform  — FLEXPART uniform kernel (bilinear over the 4 nearest centres)
    #   :adaptive — per-parcel truncated Gaussian with statistical bandwidth
    #               hband (kNN spacing, ADAPTIVE_KERNEL design): resolution is
    #               bought with trajectories, not with the output grid, and the
    #               3-h age rule disappears (h→0 at a fresh point release).
    kernel::Symbol
    # R(x) true-resolution field (output.resolution = true, adaptive kernel):
    # mass-weighted mean kernel width of everything depositing into the cell.
    resolution::Bool
    resnum::Matrix{Float64}        # Σ m·w·h   (Niout × Njout)
    resden::Matrix{Float64}        # Σ m·w
    # Adaptive-kernel bandwidth refresh cadence: the kNN pass is O(N) with a
    # large constant, so recompute hband only every `bw_stride` deposit calls
    # (build_run sets this to samples-per-output-frame — bandwidths are a
    # slowly-varying smoothing scale, refreshing once per frame is ample).
    bw_stride::Int                 # ≤1 → every deposit (legacy)
    bw_count::Int                  # deposit-call counter
    # Physical bandwidth floor (output.bandwidth_floor = "dispersion"):
    # per-parcel h ≥ Taylor σ(age; σv, T_L) — see _hfloor_taylor in conc.jl.
    bwfloor_sv::Float64            # σv [m/s]; 0 → floor off
    bwfloor_tl::Float64            # T_L [s]
end

# ── Constructor from legacy Dict — shares array references, O(1) ──────────────
function OutputConfig(oc::Dict)
    OutputConfig(
        oc["concarray"]::Array{Float32,5},
        oc["massarray"]::Array{Float64,5},
        oc["area"]   isa Matrix{Float32}    ? oc["area"]   : Matrix{Float32}(oc["area"]),
        oc["volume"] isa Array{Float32,3}   ? oc["volume"] : Array{Float32,3}(oc["volume"]),
        Vector{Float32}(oc["outheight"]),
        Vector{Float32}(oc["outinterfaces"]),
        Float32(oc["dlon"]),
        Float32(oc["dlat"]),
        Float32(oc["lon0"]),
        Float32(oc["lat0"]),
        Int32(oc["Niout"]),
        Int32(oc["Njout"]),
        Int32(length(oc["outheight"])),
        Int32(get(oc, "relbox", size(oc["concarray"], 4))),
        Int8(get(oc, "iout", 1)),
        Int8(get(oc, "direction", 1)),
        Vector{Int32}(get(oc, "boxspecies", Int32[1])),
        Vector{Float32}(get(oc, "relmass", Float32[1.0])),
        Symbol(get(oc, "outgroup", "release")),
        # legacy dicts (tests, old checkpoints): identity map over the columns,
        # classspecies == boxspecies — bit-identical per-box behaviour
        Vector{Int32}(get(oc, "outclass",
                          collect(Int32, 1:size(oc["concarray"]::Array{Float32,5}, 4)))),
        Vector{Int32}(get(oc, "classspecies",
                          Vector{Int32}(get(oc, "boxspecies", Int32[1])))),
        Bool(get(oc, "proj_m", false)),
        Matrix{Float32}(get(oc, "lon2d", Matrix{Float32}(undef, 0, 0))),
        Matrix{Float32}(get(oc, "lat2d", Matrix{Float32}(undef, 0, 0))),
        nothing,
        let k = get(oc, "kernel", true)   # legacy Bool or new String/Symbol
            k isa Bool ? (k ? :uniform : :none) : Symbol(k)
        end,
        Bool(get(oc, "resolution", false)),
        zeros(Float64, Int(oc["Niout"]), Int(oc["Njout"])),
        zeros(Float64, Int(oc["Niout"]), Int(oc["Njout"])),
        Int(get(oc, "bw_stride", 1)),
        0,
        Float64(get(oc, "bwfloor_sv", 0.0)),
        Float64(get(oc, "bwfloor_tl", 300.0)),
    )
end

function reset_accumulators!(oc::OutputConfig)
    fill!(oc.concarray,  0f0)
    fill!(oc.massarray, 0.0)
    fill!(oc.resnum, 0.0)
    fill!(oc.resden, 0.0)
end
