#!/usr/bin/env julia
# test_ka_parity.jl — threads vs KernelAbstractions (ka_cpu) advance parity on a
# SYNTHETIC met (no archive needed).  Pins the KA review fixes of 2026-09-04:
#   [1] A1  constant vertical wind, NoTurb, ctl=5: the KA kernel applied the
#           accumulated vertical displacement on every sub-step → (ctl+1)/2 = 3x
#           the CPU descent.  Now: per-parcel positions identical to threads.
#   [1b] A5 sub-step count: Δt=3 s, ctl=9 — the Float32 accumulating while-loop
#           ran 10 sub-steps (10/9 the displacement); now an integer nsub.
#   [2] A4  entry guard: OOB parcels (indx < 1 / >= Ni, indy >= Nj) are skipped
#           and their Dtime NOT advanced (CPU `continue` semantics) on both paths.
#   [3] A4  NaN in the met (masked cell) → zero displacement, finite positions.
#   [4] A2  Hanna(false) (pre-computed σ_w/T_L, the ERA5/GFS default) ran as
#           NoTurb in the KA kernel.  Now: the Thomson body on both paths →
#           comparable vertical spread inside a synthetic PBL.
#   [5] A7  per-parcel RNG seed depends on (slot, first-advance time, run seed):
#           the same slot re-injected later gets a different stream.
#   [6] A5  optional (EASLEE_GPU=1 + functional CUDA): same checks on CUDABackend
#           (exercises the per-nest resident buffers + met-change detection).
#
#   julia --project=easlee test/test_ka_parity.jl
#   EASLEE_GPU=1 julia --project=easlee test/test_ka_parity.jl   # + GPU case
using Printf, Dates, Random, Statistics

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
# tools.jl (legacy) defines Rectilinear/Thomson/Hanna in Main too — qualify the
# typed ones the kernels dispatch on.
const Ty = EASLEE.Types
using .EASLEE.Types: SigmaPressureLevel, MetFields3D, MetFields2D, NestRunConfig
using .EASLEE.Operators
using KernelAbstractions

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-70s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# ── synthetic met: flat terrain, 250 m levels, neutral column, PBL to 1 km ────
# Layout contract 2026-10-22: zlev and tke on the half levels (zlev[1] = ground),
# theta/rho class F (slot 1 = surface, slot k at ½(zlev[k-1] + zlev[k])).  rho was
# co-located with zlev until then; drhodz is no longer read by the turbulence.
const NI = 16; const NJ = 16; const NZ = 40; const DZ = 250.0
const PBLH = 1000.0                       # σ_w = 1 m/s below, floor above
const W0   = 1.0f-3                       # levdot [level/s] → 0.9 level per 900 s step
const DT   = Millisecond(900_000); const CTL = 5.0f0
dims = (NI, NJ, NZ, 2)
zlev = zeros(Float32, dims); theta = zeros(Float32, dims); rho = zeros(Float32, dims)
drhodz = zeros(Float32, dims); tke = zeros(Float32, dims)
for k in 1:NZ
    z  = (k - 1) * DZ
    zf = k == 1 ? 0.0 : (k - 1.5) * DZ                  # full-level height of slot k
    zlev[:, :, k, :]   .= Float32(z)
    theta[:, :, k, :]  .= 300.0f0                       # neutral: N² = 0
    rho[:, :, k, :]    .= Float32(1.2 * exp(-zf / 8000.0))
    drhodz[:, :, k, :] .= Float32(-1.2 / 8000.0 * exp(-zf / 8000.0))
    # Hanna-fill convention: field = 0.5·σ_w²  (σ_w = 1 m/s in the PBL)
    tke[:, :, k, :]    .= z < PBLH ? 0.5f0 : 0.0f0
end
Ux = fill(0.5f0, dims); Vy = fill(-0.3f0, dims)
levdot_w = fill(W0, dims)            # constant vertical wind (cases 1-3)
levdot_0 = zeros(Float32, dims)      # no mean vertical wind   (case 4)
# NaN "terrain" block for case 3: columns 13-14, all levels, both slots
levdot_nan = copy(levdot_w); levdot_nan[13:14, :, :, :] .= NaN32
zero4() = zeros(Float32, dims)
mk_m3(lev) = MetFields3D{Array{Float32,4}}(Ux, Vy, lev, rho, theta, tke, zlev, drhodz,
                                            zero4(), zero4(), zero4(), zero4())
m2 = MetFields2D{Array{Float32,3}}(fill(Float32(PBLH), NI, NJ, 2), fill(0.3f0, NI, NJ, 2),
                                   fill(1.0f0, NI, NJ, 2), fill(-100.0f0, NI, NJ, 2))
const DLON = 0.25f0; const DLAT = 0.25f0; const R_EARTH = 6.371e6
GC = Ty.Rectilinear(Int32(NI), Int32(NJ), Int32(NZ), DLON, DLAT, 0.0f0, 40.0f0,
                 Float32(180.0 / (DLON * R_EARTH * pi)), Float32(180.0 / (DLAT * R_EARTH * pi)),
                 1.0f0, 0.0f0, 1.0f0, 0.0f0, true, true, true, DLON, DLAT, false)
LC = SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
mk_nrc(turb; ctl = CTL, mesoscale = false) =
    NestRunConfig((1, 2), turb, Float32(ctl), Int8(1), (true, true, true),
                  mesoscale, 0.16f0, 3600.0f0, 0.0f0,
                  false, false,                    # horiz_tke, horiz_hanna
                  false, NaN32, NaN32,             # vmix_native / rw / ck (this distribution has no vmix_kz field)
                  false, false, Inf32, -Inf32, false)
mk_ctx(turb, lev; ctl = CTL, dt = DT, mesoscale = false) =
    (runconfig = (timestep = dt,),
     nests = [(nestrunconfig = mk_nrc(turb; ctl, mesoscale), metfields3D = mk_m3(lev),
               metfields2D = m2, gridconfig = GC, levelconfig = LC)])
const T0 = DateTime(2022, 1, 15, 4)
const INDT = [1.5f0]

# ── parcel population (deterministic) ────────────────────────────────────────
function mk_state(N; indl_lo, indl_hi, seed = 7)
    rng = MersenneTwister(seed)
    st = EASLEE.State.ParcelState(N + 8)
    for i in 1:N
        st.status[i] = 0x01; st.nest[i] = Int32(1)
        st.indx[i] = 4.0f0 + 8.0f0 * rand(rng, Float32)
        st.indy[i] = 4.0f0 + 8.0f0 * rand(rng, Float32)
        st.indl[i] = Float32(indl_lo + (indl_hi - indl_lo) * rand(rng, Float32))
        st.Dtime[i] = T0; st.Dtimeinit[i] = T0; st.mass[i] = 1.0f0; st.dt[i] = 1.0f0
    end
    return st
end
mk_state(N, indl) = mk_state(N; indl_lo = indl, indl_hi = indl)

"Advance `st` for `nsteps` on `backend` (:threads, :ka_cpu, or a KA backend object)."
function run!(st, ctx, backend; nsteps = 1)
    t = T0
    DTc = ctx.runconfig.timestep
    for _ in 1:nsteps
        if backend === :threads
            Operators.advance_active_parcels!(st, ctx, t, INDT)
        elseif backend === :ka_cpu
            Operators.advance_active_parcels_ka!(st, ctx, t, INDT;
                                                 backend = KernelAbstractions.CPU(; static = true))
        else
            Operators.advance_active_parcels_ka!(st, ctx, t, INDT; backend = backend)
        end
        t += DTc
    end
    return t
end
active(st) = findall(==(0x01), st.status)

function case_constant_w(backend, label)
    println("="^88)
    println("[1] $label — constant vertical wind, NoTurb, ctl=$(Int(CTL)): parity with threads")
    println("="^88)
    N = 2000; nsteps = 4
    ctx = mk_ctx(Ty.NoTurb(), levdot_w)
    a = mk_state(N, 20.0); b = mk_state(N, 20.0)
    tend = run!(a, ctx, :threads; nsteps); run!(b, ctx, backend; nsteps)
    ia = active(a)
    d_thr = mean(a.indl[ia] .- 20.0f0); d_ka = mean(b.indl[ia] .- 20.0f0)
    expect = W0 * 900.0 * nsteps
    @printf("  mean Δindl  threads %.4f   %s %.4f   analytic %.4f  (ratio %s/threads = %.3f)\n",
            d_thr, label, d_ka, expect, label, d_ka / d_thr)
    check("threads mean Δindl = W·Δt·nsteps (±1%)", isapprox(d_thr, expect; rtol = 0.01))
    check("$label mean Δindl within 3% of threads (was 3x before A1)",
          abs(d_ka - d_thr) / abs(d_thr) < 0.03, info = @sprintf("ratio %.4f", d_ka / d_thr))
    check("$label per-parcel indl identical to threads (max|Δ| < 1e-3 lvl)",
          maximum(abs.(a.indl[ia] .- b.indl[ia])) < 1e-3,
          info = @sprintf("max|Δ| %.2e", maximum(abs.(a.indl[ia] .- b.indl[ia]))))
    check("$label per-parcel indx/indy identical to threads (max|Δ| < 1e-3)",
          max(maximum(abs.(a.indx[ia] .- b.indx[ia])),
              maximum(abs.(a.indy[ia] .- b.indy[ia]))) < 1e-3)
    check("$label Dtime advanced like threads", all(b.Dtime[ia] .== tend) && all(a.Dtime[ia] .== tend))
    check("no NaN in $label positions", all(isfinite, b.indx[ia]) && all(isfinite, b.indl[ia]))
end

function case_substep_count(backend, label)
    println(); println("="^88)
    println("[1b] $label — sub-step count: Δt=3 s, ctl=9 (Float32 accumulation ran 10 sub-steps)")
    println("="^88)
    nsteps = 10
    ctx = mk_ctx(Ty.NoTurb(), levdot_w; ctl = 9, dt = Millisecond(3000))
    a = mk_state(200, 20.0); b = mk_state(200, 20.0)
    run!(a, ctx, :threads; nsteps); run!(b, ctx, backend; nsteps)
    ia = active(a)
    d_thr = mean(a.indl[ia] .- 20.0f0); d_ka = mean(b.indl[ia] .- 20.0f0)
    expect = W0 * 3.0 * nsteps
    @printf("  mean Δindl  threads %.5f   %s %.5f   analytic %.5f  (ratio %s/analytic = %.4f)\n",
            d_thr, label, d_ka, expect, label, d_ka / expect)
    check("$label Δindl = W·Δt·nsteps within 0.5% (was 10/9 with ctl+1 sub-steps)",
          isapprox(d_ka, expect; rtol = 0.005), info = @sprintf("ratio %.4f", d_ka / expect))
    check("$label per-parcel indl identical to threads (max|Δ| < 1e-4 lvl)",
          maximum(abs.(a.indl[ia] .- b.indl[ia])) < 1e-4,
          info = @sprintf("max|Δ| %.2e", maximum(abs.(a.indl[ia] .- b.indl[ia]))))
end

function case_entry_guard(backend, label)
    println(); println("="^88)
    println("[2] $label — entry guard: OOB parcels skipped, Dtime not advanced (CPU `continue`)")
    println("="^88)
    ctx = mk_ctx(Ty.NoTurb(), levdot_w)
    # (indx, indy) of 4 parcels outside the CPU entry test: indx<1, indx>=Ni, indy>=Nj, indy<1
    oob = [(0.5f0, 8.0f0), (Float32(NI), 8.0f0), (8.0f0, Float32(NJ) + 0.2f0), (8.0f0, 0.2f0)]
    for backendsym in (:threads, backend)
        st = mk_state(4, 20.0)
        for r in 1:4; st.indx[r] = oob[r][1]; st.indy[r] = oob[r][2]; end
        st.indx[5] = 8.0f0; st.indy[5] = 8.0f0; st.status[5] = 0x01; st.nest[5] = Int32(1)
        st.indl[5] = 20.0f0; st.Dtime[5] = T0; st.dt[5] = 1.0f0     # control (in range)
        run!(st, ctx, backendsym)
        lb = backendsym === :threads ? "threads" : label
        check("$lb: 4 OOB parcels untouched (indl, indx unchanged)",
              all(st.indl[1:4] .== 20.0f0) && all(st.indx[r] == oob[r][1] for r in 1:4))
        check("$lb: OOB parcels' Dtime stays at currentdate", all(st.Dtime[1:4] .== T0),
              info = string(unique(st.Dtime[1:4])))
        check("$lb: control parcel advanced", st.Dtime[5] == T0 + DT && st.indl[5] > 20.5f0)
    end
end

function case_nan_met(backend, label)
    println(); println("="^88)
    println("[3] $label — NaN levdot (masked cells) → zero vertical displacement, finite positions")
    println("="^88)
    ctx = mk_ctx(Ty.NoTurb(), levdot_nan)
    a = mk_state(2, 20.0); b = mk_state(2, 20.0)
    for st in (a, b)
        st.indx[1] = 13.5f0; st.indy[1] = 8.0f0     # inside the NaN block
        st.indx[2] = 8.0f0;  st.indy[2] = 8.0f0     # clean column
    end
    run!(a, ctx, :threads); run!(b, ctx, backend)
    @printf("  NaN-cell parcel Δindl: threads %.4f  %s %.4f ; clean parcel: threads %.4f  %s %.4f\n",
            a.indl[1] - 20, label, b.indl[1] - 20, a.indl[2] - 20, label, b.indl[2] - 20)
    check("$label NaN-cell parcel finite", isfinite(b.indl[1]) && isfinite(b.indx[1]))
    check("$label NaN-cell parcel: zero vertical displacement like threads",
          abs(b.indl[1] - a.indl[1]) < 1e-3 && abs(a.indl[1] - 20.0f0) < 1e-3)
    check("$label clean parcel identical to threads", abs(b.indl[2] - a.indl[2]) < 1e-3)
end

function case_hanna_turb(backend, label)
    println(); println("="^88)
    println("[4] $label — Hanna(false) turbulence in a 1 km PBL: vertical spread vs threads")
    println("="^88)
    N = 4000; nsteps = 8
    ctx = mk_ctx(Ty.Hanna(false), levdot_0)
    a = mk_state(N; indl_lo = 1.5, indl_hi = 4.5); b = mk_state(N; indl_lo = 1.5, indl_hi = 4.5)
    d0 = copy(a.indl)
    run!(a, ctx, :threads; nsteps); run!(b, ctx, backend; nsteps)
    ia = active(a); ib = active(b)
    s_thr = std(a.indl[ia] .- d0[ia]); s_ka = std(b.indl[ib] .- d0[ib])
    z_thr = mean(a.indl[ia] .- 1) * DZ; z_ka = mean(b.indl[ib] .- 1) * DZ
    @printf("  std(Δindl) over %d steps: threads %.3f lvl (%.0f m)   %s %.3f lvl (%.0f m)   ratio %.3f\n",
            nsteps, s_thr, s_thr * DZ, label, s_ka, s_ka * DZ, s_ka / s_thr)
    @printf("  mean height after %d steps: threads %.0f m   %s %.0f m   (PBL %.0f m)\n",
            nsteps, z_thr, label, z_ka, PBLH)
    check("threads Hanna(false) produces vertical spread (std Δindl > 0.3 lvl)", s_thr > 0.3)
    check("$label Hanna(false) produces vertical spread (was exactly 0 before A2)", s_ka > 0.3,
          info = @sprintf("std %.3f", s_ka))
    check("$label spread comparable to threads (ratio in [0.6, 1.7])",
          0.6 < s_ka / s_thr < 1.7, info = @sprintf("ratio %.3f", s_ka / s_thr))
    check("$label mean height comparable to threads (within 25% of PBL depth)",
          abs(z_ka - z_thr) < 0.25 * PBLH)
    check("$label parcels alive and finite", length(ib) == N && all(isfinite, b.indl[ib]))
end

function case_rng_seed()
    println(); println("="^88)
    println("[5] A7 — per-parcel seed depends on (slot, first-advance time, run seed)")
    println("="^88)
    rs = Operators.RNG_RUN_SEED[]
    s1 = Operators.rng_seed(17, Int64(Dates.value(T0)), rs)
    s2 = Operators.rng_seed(17, Int64(Dates.value(T0 + DT)), rs)
    s3 = Operators.rng_seed(17, Int64(Dates.value(T0)), rs ⊻ UInt64(1))
    s4 = Operators.rng_seed(18, Int64(Dates.value(T0)), rs)
    check("same slot, later injection → different seed", s1 != s2)
    check("same slot+time, different run seed → different seed", s1 != s3)
    check("neighbouring slots differ", s1 != s4)
    check("seeds never 0 (xorshift absorbing state)", all(!=(UInt64(0)), (s1, s2, s3, s4)))
    # the KA wrapper path: zeroed slot (boundary reuse) re-seeded from Dtime
    ctx = mk_ctx(Ty.NoTurb(), levdot_w)
    st = mk_state(1, 20.0)
    Operators.advance_active_parcels_ka!(st, ctx, T0, INDT; backend = KernelAbstractions.CPU(; static = true))
    r1 = st.rng_state[1]
    st.rng_state[1] = UInt64(0); st.Dtime[1] = T0 + Day(1)     # slot reused a day later
    Operators.advance_active_parcels_ka!(st, ctx, T0 + Day(1), INDT; backend = KernelAbstractions.CPU(; static = true))
    check("_ensure_rng_states!: reused slot gets a different stream", st.rng_state[1] != r1)
end

function case_dispatch_guard()
    println(); println("="^88)
    println("[7] A3/2026-09-21 — advance_dispatch!: nothing is threads-only any more; the wind mask still is,")
    println("    and an EXPLICIT backend request refuses to downgrade (error) while the auto default warns once")
    println("="^88)
    ENV["EASLEE_ADVANCE"] = "ka_cpu"            # read lazily by _advance_mode()
    Operators.ADVANCE_MODE[] = :unset; Operators._ADVANCE_REQUESTED[] = false
    Operators._KA_PHYS_WARNED[] = false; Operators._WIND_WARNED[] = false
    st = mk_state(50, 20.0)
    Operators.advance_dispatch!(st, mk_ctx(Ty.NoTurb(), levdot_w), T0, INDT)
    check("EASLEE_ADVANCE=ka_cpu selected for a KA-supported config",
          Operators.ADVANCE_MODE[] == :ka_cpu && !Operators._KA_PHYS_WARNED[])
    check("supported config: no override", Operators._ka_unsupported_switch(mk_ctx(Ty.Hanna(false), levdot_w)) === nothing)
    ctxm = mk_ctx(Ty.NoTurb(), levdot_w; mesoscale = true)
    check("mesoscale=true is supported on the KA path (2026-09-21 port)", Operators._ka_unsupported_switch(ctxm) === nothing)
    check("Hanna on-the-fly is supported on the KA path (turbulence_hanna_4d)",
          Operators._ka_unsupported_switch(mk_ctx(Ty.Hanna(true), levdot_w)) === nothing)
    # the one remaining threads-only feature: a frozen wind dimension (2-D / 1-D runs)
    ctx2d = mk_ctx(Ty.NoTurb(), levdot_w); ctx2d.nests[1].nestrunconfig.wind_mask = (true, true, false)
    st2 = mk_state(50, 20.0)
    threw = try; Operators.advance_dispatch!(st2, ctx2d, T0, INDT); false; catch e; occursin("Refusing to downgrade", sprint(showerror, e)); end
    check("explicit EASLEE_ADVANCE=ka_cpu + 2-D wind mask -> ERROR, not a silent threads run", threw)
    Operators.ADVANCE_MODE[] = :unset; delete!(ENV, "EASLEE_ADVANCE"); Operators._ADVANCE_REQUESTED[] = false
    Operators._WIND_WARNED[] = false
    st3 = mk_state(50, 20.0)
    Operators.advance_dispatch!(st3, ctx2d, T0, INDT)
    check("auto default (no EASLEE_ADVANCE) + 2-D wind mask -> threads, warned once, parcels advanced",
          Operators.ADVANCE_MODE[] == :threads && all(st3.Dtime[active(st3)] .== T0 + DT))
    Operators.ADVANCE_MODE[] = :unset; delete!(ENV, "EASLEE_ADVANCE")
end

# ── run ──────────────────────────────────────────────────────────────────────
case_constant_w(:ka_cpu, "ka_cpu")
case_substep_count(:ka_cpu, "ka_cpu")
case_entry_guard(:ka_cpu, "ka_cpu")
case_nan_met(:ka_cpu, "ka_cpu")
case_hanna_turb(:ka_cpu, "ka_cpu")
case_rng_seed()
case_dispatch_guard()

if get(ENV, "EASLEE_GPU", "0") == "1"
    if Operators.GPU_ENABLED && Operators._gpu_functional()
        println(); println("#"^88); println("GPU (CUDABackend) cases"); println("#"^88)
        gb = Operators.CUDA.CUDABackend()
        case_constant_w(gb, "ka_gpu")
        case_substep_count(gb, "ka_gpu")
        case_nan_met(gb, "ka_gpu")
        case_hanna_turb(gb, "ka_gpu")
        # A5: met-change detection — same buffers, changed levdot must be re-uploaded
        println(); println("[6] A5 — resident met buffers follow an in-place met change")
        ctx1 = mk_ctx(Ty.NoTurb(), copy(levdot_w))
        st = mk_state(200, 20.0)
        Operators.advance_active_parcels_ka!(st, ctx1, T0, INDT; backend = gb)
        d1 = mean(st.indl[active(st)] .- 20)
        ctx1.nests[1].metfields3D.levdot .*= 2.0f0        # in-place, t_index unchanged
        for i in active(st); st.indl[i] = 20.0f0; st.Dtime[i] = T0; end
        Operators.advance_active_parcels_ka!(st, ctx1, T0, INDT; backend = gb)
        d2 = mean(st.indl[active(st)] .- 20)
        @printf("  Δindl before %.4f  after doubling levdot in place %.4f\n", d1, d2)
        check("ka_gpu: in-place met change detected (Δ doubles)", isapprox(d2, 2 * d1; rtol = 0.02))
        check("ka_gpu: one buffer set per nest", length(Operators._KAGPU) == 1 && haskey(Operators._KAGPU, 1))
    else
        println("\nEASLEE_GPU=1 but CUDA not functional — GPU cases skipped")
    end
end

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
