# test_conv_entrain_detrain.jl ── ONE-STEP Monte-Carlo validation of the convective
# redistribution against the ARCHIVED ERA5/IFS convection fields.
#
# The well-mixed test (well_mixed.jl MODE=convection) proves a uniform tracer stays
# uniform, but says nothing about WHICH transport produced that fixed point.  This
# test closes that gap: parcels are seeded uniform in air mass in the strongest
# convective column, ONE convective step is applied, and the parcel statistics are
# compared against the archived fields — including avg_udr/avg_ddr, which the
# mass-flux scheme NEVER READS (it uses only the umf/dmf flux divergence), so the
# detrainment comparison is against an independent archived diagnostic.
#
#   massflux (CFG=configs/era5_conv.toml):
#     • F_up(k)  = m_pp/Δt · #{parcels with k0 ≤ k < k1}  must equal the archived
#       interface updraught flux  Mu_if(k) = ½(umf[k]+umf[k+1]);
#     • F_dn(k)  likewise must equal Mu_if(k) (downdraught Md_if + compensating
#       subsidence Mu_if−Md_if) — so up- and down-transport cancel exactly;
#     • detrainment: destination histogram of upward movers vs the flux-divergence
#       prediction (Mu_if[l−1]−Mu_if[l])⁺ and vs the ARCHIVED udr·Δz;
#     • entrainment: origin histogram vs (Mu_if[k]−Mu_if[k−1])⁺ and vs the archived-
#       implied  E = D_archived + ΔMu_if  (IFS archives D, not E: E = ∂Mu/∂z + D).
#
#   cape (CFG=configs/era5_cape.toml):
#     • the CAPE/CIN gate fires at the max-CAPE cell;
#     • moved parcels land UNIFORM IN AIR MASS below the diagnosed EL
#       (destination mass histogram ∝ dm(k), flat χ);
#     • parcels above the EL are untouched.
#
#   CFG=configs/era5_conv.toml NP=2000000 DT=900 julia -t 8 --project=easlee \
#     test/test_conv_entrain_detrain.jl
#
# Output: results/convtest_<name>.tsv (per-level profiles) + printed summary.
#
# Vertical layout (stacked, as the readers deliver it): zlev on the half levels, rho,
# umf, dmf, udr, ddr, T, puv, qv on the full levels; layer k (k ≥ 2) is the full level
# between zlev[k-1] and zlev[k], air mass rho[k]·Δz, parcels with indl ∈ [k-1, k);
# interface k sits at zlev[k] with flux ½(umf[k]+umf[k+1]) (0 at the ground, k = 1).

using Dates, Printf, Statistics, Random
global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
include(joinpath(REPO, "easlee", "runtime.jl"))
include(joinpath(REPO, "easlee", "build_run.jl"))
using .EASLEE.Operators: convective_redistribution!, diagnostic_convection!
const Op = EASLEE.Operators

const CFG = get(ENV, "CFG", joinpath(REPO, "configs", "era5_conv.toml"))
const NP  = parse(Int, get(ENV, "NP", "2000000"))
const DTS = parse(Float64, get(ENV, "DT", "900"))     # one convective step [s]
const ZCAP = 20000.0f0                                # seed the column up to here

spec = read_runspec(CFG)
for r in spec.release; r.particles = 1; r.time = [spec.stop + Hour(1), spec.stop + Hour(1)]; end
rt, obj = build_run(spec; prefetch=false, extra_slots=2NP)
ctx = rt.ctx; state = rt.state; rc = ctx.runconfig; nest = ctx.nests[1]
Ni = Int(nest.gridconfig.Ni); Nj = Int(nest.gridconfig.Nj)
Nk = Int(nest.levelconfig.Nk)

loopstart = rc.direction == -1 ? rc.enddate : rc.startdate
read_met_step!(rt, loopstart)
t_array = EASLEE.MetAdapter.met_time_index_array(loopstart, ctx.nests)
indt = t_array[1]

massflux = obj.convfields !== nothing
capemode = !massflux && obj.diagconv !== nothing
massflux || capemode || error("$(spec.name): neither convfields nor diagconv built — check [physics.convection].scheme")
scheme = massflux ? "massflux" : "cape"
@printf("=== conv entrain/detrain test [%s]  cfg=%s  NP=%d  Δt=%.0fs ===\n",
        scheme, basename(CFG), NP, DTS)

# ── strongest convective column (raw array scan at the integer time slice) ────
ti = clamp(Base.unsafe_trunc(Int, indt), 1, 2)
ic, jc = 0, 0; bestv = -Inf32
if massflux
    A = obj.convfields[1].umf
    for i in 3:Ni-2, j in 3:Nj-2
        m = 0.0f0
        @inbounds for k in 1:Nk; v = A[i,j,k,ti]; (isfinite(v) && v > m) && (m = v); end
        m > bestv && (global bestv = m; global ic = i; global jc = j)
    end
    @printf("strongest umf column: (i=%d, j=%d)  max umf=%.4f kg m⁻² s⁻¹\n", ic, jc, bestv)
else
    A = obj.diagconv[1].cape
    for i in 3:Ni-2, j in 3:Nj-2
        v = A[i,j,ti]
        (isfinite(v) && v > bestv) && (global bestv = v; global ic = i; global jc = j)
    end
    cin_here = obj.diagconv[1].cin[ic,jc,ti]
    @printf("max-CAPE cell: (i=%d, j=%d)  CAPE=%.0f J/kg  CIN=%.0f J/kg  (gate: cape≥%.0f, cin≤%.0f)\n",
            ic, jc, bestv, cin_here, obj.diagconvparams.cape_min, obj.diagconvparams.cin_max)
end

# ── column profiles, extracted EXACTLY as the scheme sees them ────────────────
fx = Float32(ic); fy = Float32(jc)
prof(Af) = (b = Vector{Float32}(undef, Nk); Op.get_vertical_profile!(b, Af, fx, fy, indt, Nk); b)
z   = prof(nest.metfields3D.zlev)
rho = prof(nest.metfields3D.rho)

# layer air mass dm[k] = rho[k]·(zlev[k] − zlev[k-1]), same floors as _build_conv_column!
dm = zeros(Float32, Nk)
for k in 2:Nk
    dz = max(z[k] - z[k-1], 1.0f0)
    dm[k] = max(rho[k], 1f-4) * dz
end
keep = [k for k in 2:Nk if isfinite(dm[k]) && dm[k] > 0 && isfinite(z[k-1]) && z[k-1] < ZCAP]

# ── seed NP parcels uniform in air mass, pinned to the column ─────────────────
w = Float64[dm[k] for k in keep]; w ./= sum(w)
cnt = round.(Int, w .* NP); cnt[argmax(w)] += NP - sum(cnt)
mass_pp = sum(dm[k] for k in keep) / NP            # each parcel = this much column mass [kg/m²]
rng = MersenneTwister(11); fill!(state.status, 0x00); s = 0
for (kk, k) in enumerate(keep), _ in 1:cnt[kk]
    global s += 1; s > length(state.status) && break
    @inbounds begin
        state.status[s] = 0x01; state.nest[s] = Int8(1); state.relbox[s] = Int32(1)
        state.indx[s] = fx; state.indy[s] = fy
        state.indl[s] = Float32(k - 1) + rand(rng, Float32)      # layer k: indl ∈ [k-1, k)
        state.Dtime[s] = loopstart; state.Dtimeinit[s] = loopstart
        state.is_new[s] = false; state.dt[s] = 1f0; state.mass[s] = Float32(mass_pp)
        state.uturb[s] = 0f0; state.vturb[s] = 0f0; state.wturb[s] = 0f0
    end
end
np = count(==(0x01), state.status); EASLEE.State.sync_dtime_ns!(state)
@printf("seeded %d parcels uniform-in-air-mass over levels %d…%d (z<%.0f m), m_pp=%.4g kg/m²\n",
        np, keep[1], keep[end], ZCAP, mass_pp)

layer_of(x) = clamp(Base.unsafe_trunc(Int, x) + 1, 2, Nk)
k0  = [layer_of(state.indl[i]) for i in 1:np]
zp0 = [Op._z_at(z, state.indl[i], Nk) for i in 1:np]    # parcel heights before the step

# ── ONE convective step ───────────────────────────────────────────────────────
if massflux
    convective_redistribution!(state, ctx, obj.convfields, DTS, t_array)
else
    diagnostic_convection!(state, ctx, obj.diagconv, DTS, t_array; params=obj.diagconvparams)
end
k1 = [layer_of(state.indl[i]) for i in 1:np]

nmoved = count(k1 .!= k0)
@printf("moved: %d parcels (%.2f%%)\n", nmoved, 100nmoved/np)

# ── χ after one step (mass-conservation sanity, both schemes) ─────────────────
cdest = zeros(Int, Nk); for x in k1; cdest[x] += 1; end
chiN = [cdest[k] * mass_pp / dm[k] for k in keep]
@printf("χ after 1 step: mean=%.4f  std=%.4f  min=%.3f  max=%.3f  (uniform ⇒ 1±noise)\n",
        mean(chiN), std(chiN), minimum(chiN), maximum(chiN))

mkpath(joinpath(REPO, "results"))
out = joinpath(REPO, "results", "convtest_$(spec.name).tsv")

if massflux
    cf = obj.convfields[1]
    mu  = prof(cf.umf);  mu  .= max.(mu, 0f0)
    mdr = prof(cf.dmf);  md  = max.(-mdr, 0f0)          # magnitude of downdraught
    udr = prof(cf.udr);  udr .= max.(udr, 0f0)
    ddr = prof(cf.ddr);  ddr .= max.(ddr, 0f0)

    muif(k) = (k < 2 || k >= Nk) ? 0f0 : 0.5f0*(mu[k]+mu[k+1])   # interface fluxes at zlev[k]
    mdif(k) = (k < 2 || k >= Nk) ? 0f0 : 0.5f0*(md[k]+md[k+1])

    # parcel-derived interface crossings (net displacement)
    fup = zeros(Float64, Nk); fdn = zeros(Float64, Nk)
    for i in 1:np
        a, b = k0[i], k1[i]
        if b > a;     for k in a:b-1;  fup[k] += 1; end
        elseif b < a; for k in b:a-1;  fdn[k] += 1; end
        end
    end
    fup .*= mass_pp / DTS; fdn .*= mass_pp / DTS

    # entrainment (origins) / detrainment (destinations) of NET-UPWARD movers
    epar = zeros(Float64, Nk); dpar = zeros(Float64, Nk)
    for i in 1:np
        k1[i] > k0[i] || continue
        epar[k0[i]] += 1; dpar[k1[i]] += 1
    end
    epar .*= mass_pp / DTS; dpar .*= mass_pp / DTS

    dflux = [max(muif(k-1) - muif(k), 0f0) for k in 1:Nk]   # flux-divergence detrainment
    eflux = [max(muif(k) - muif(k-1), 0f0) for k in 1:Nk]   # flux-divergence entrainment
    dz    = [k > 1 ? max(z[k]-z[k-1], 1f0) : 1f0 for k in 1:Nk]  # layer k thickness
    darch = [udr[k] * dz[k] for k in 1:Nk]                  # archived detrainment per layer
    earch = [darch[k] + (muif(k) - muif(k-1)) for k in 1:Nk] # archived-implied entrainment

    # summary correlations/errors over the convective layers (Mu_if > eps)
    act = [k for k in keep if muif(k) > 1f-4]
    r2(a,b) = (c = cor(Float64.(a), Float64.(b)); c*c)
    nrmse(a,b) = sqrt(mean((Float64.(a).-Float64.(b)).^2)) / max(maximum(Float64.(b)), eps())
    @printf("\n── interface fluxes over %d active interfaces (Mu_if > 1e-4) ──\n", length(act))
    @printf("F_up  vs Mu_if : R²=%.4f  NRMSE=%.3f   (Σ F_up=%.4f  Σ Mu_if=%.4f)\n",
            r2(fup[act], [muif(k) for k in act]), nrmse(fup[act], [muif(k) for k in act]),
            sum(fup[act]), sum(muif(k) for k in act))
    @printf("F_dn  vs Mu_if : R²=%.4f  NRMSE=%.3f   (Σ F_dn=%.4f)\n",
            r2(fdn[act], [muif(k) for k in act]), nrmse(fdn[act], [muif(k) for k in act]), sum(fdn[act]))
    da = [k for k in keep if dflux[k] > 1f-5]
    @printf("detrain(parcel) vs flux-div : R²=%.4f    vs archived udr·Δz : R²=%.4f\n",
            r2(dpar[da], dflux[da]), r2(dpar[da], darch[da]))
    ea = [k for k in keep if eflux[k] > 1f-5]
    @printf("entrain(parcel) vs flux-div : R²=%.4f    vs archived E=D+ΔMu : R²=%.4f\n",
            r2(epar[ea], eflux[ea]), r2(epar[ea], max.(earch[ea], 0.0)))

    open(out, "w") do io
        println(io, "# one-step conv test  scheme=massflux  cfg=$(basename(CFG))  NP=$np  dt=$(DTS)s  cell=($ic,$jc)  umfmax=$bestv")
        println(io, "# layer k = full level k between zbot = zlev[k-1] and ztop = zlev[k]; muif, mdif, fup, fdn at ztop")
        println(io, "k\tzbot\tztop\trho\tdm\tmu\tmd\tmuif\tmdif\tfup\tfdn\tepar\tdpar\teflux\tdflux\tearch\tdarch\tchi")
        for k in keep
            @printf(io, "%d\t%.1f\t%.1f\t%.5g\t%.3f\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.4f\n",
                    k, z[k-1], z[k], rho[k], dm[k], mu[k], md[k], muif(k), mdif(k), fup[k], fdn[k],
                    epar[k], dpar[k], eflux[k], dflux[k], max(earch[k],0), darch[k],
                    cdest[k]*mass_pp/dm[k])
        end
    end
else
    # ── CAPE scheme: recompute the EL exactly as the scheme does, then check the
    #    moved parcels are uniform in air mass below it.
    T  = prof(nest.metfields3D.T)
    pp = prof(nest.metfields3D.puv)
    q2 = if obj.diagconv[1].qv !== nothing
        prof(obj.diagconv[1].qv)[2]           # archived humidity (the scheme's path)
    else
        (pp[2] / (287.05f0 * max(rho[2],1f-4) * max(T[2],150f0)) - 1f0) / 0.608f0
    end
    the_p = Op._theta_e_bolton(T[2], pp[2], max(q2, 0f0))
    kel = 0; zel = NaN32                      # EL height = full level kel, ½(zlev[kel-1]+zlev[kel])
    for k in Nk:-1:3
        zf = 0.5f0 * (z[k-1] + z[k])
        (isfinite(zf) && zf <= 17000f0) || continue
        thes = Op._theta_e_bolton(T[k], max(pp[k],100f0), Op._qsat_bolton(T[k], max(pp[k],100f0)))
        if the_p > thes; global kel = k; global zel = zf; break; end
    end
    @printf("diagnosed EL: level %d  z=%.0f m  (θe_parcel=%.1f K, q_src=%.4g)\n",
            kel, zel, the_p, q2)

    moved_up   = count(i -> k1[i] != k0[i], 1:np)
    below      = kel > 0 ? [i for i in 1:np if k0[i] < kel] : Int[]
    movedbelow = count(i -> k1[i] != k0[i], below)
    above      = kel > 0 ? [i for i in 1:np if zp0[i] > zel] : Int[]
    touchedabove = count(i -> k1[i] != k0[i], above)
    @printf("below-EL parcels moved: %d/%d (%.1f%%)   above-EL touched: %d (must be 0)\n",
            movedbelow, length(below), 100movedbelow/max(length(below),1), touchedabove)

    # destination χ of moved parcels: uniform in air mass over [0, z_EL] ⇒ flat
    cmv = zeros(Int, Nk); for i in 1:np; k1[i] != k0[i] && (cmv[k1[i]] += 1); end
    kk = kel > 0 ? [k for k in keep if k < kel] : keep
    chimv = [cmv[k] / dm[k] for k in kk]; chimv ./= mean(chimv)
    @printf("moved-parcel destination χ over [0,z_EL): mean=1  std=%.4f  min=%.3f  max=%.3f\n",
            std(chimv), minimum(chimv), maximum(chimv))
    open(out, "w") do io
        println(io, "# one-step conv test  scheme=cape  cfg=$(basename(CFG))  NP=$np  dt=$(DTS)s  cell=($ic,$jc)  cape=$bestv  kel=$kel")
        println(io, "k\tzbot\tztop\tdm\tcmoved\tchi_moved\tchi_all")
        for k in keep
            cm = cmv[k]/dm[k]; cnorm = isempty(kk) ? NaN : cm/mean([cmv[q]/dm[q] for q in kk])
            @printf(io, "%d\t%.1f\t%.1f\t%.3f\t%d\t%.4f\t%.4f\n", k, z[k-1], z[k], dm[k], cmv[k],
                    k in kk ? cnorm : NaN, cdest[k]*mass_pp/dm[k])
        end
    end
end
@printf("\n✏  wrote %s\n", out)
