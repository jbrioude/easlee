# test_aerosol.jl — aerosol transport validation.
#
# Part 1 (no met data): unit physics —
#   settling velocities vs textbook values, Mie Qext limits, dry-deposition
#   velocity shape (minimum in the accumulation mode), lognormal bin setup.
# Part 2 (ERA5+CLD): 3 h integration with settling + dry dep + size-dependent
#   wet dep — mass budget closure, settling-rate check in a quiet column,
#   AOD map sanity.
#
# Usage:  julia --threads=8 --project=../easlee test_aerosol.jl [config.toml]

using Dates, Printf, Statistics
include("../easlee/src/io/initialize.jl")
include("../easlee/src/io/tools.jl")
include("../easlee/src/io/output.jl")
include("../easlee/src/config/EasleeConfig.jl"); using .EasleeConfig
include("../easlee/EASLEE.jl")
using .EASLEE
using .EASLEE.MetAdapter
using .EASLEE.State
using .EASLEE.Types: NoTurb as EASLEENoTurb
include("../easlee/runtime.jl")
include("../easlee/build_run.jl")

global gpu = false
npass = 0; nfail = 0
function check(label, ok, detail)
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-48s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end

println("── Part 1: unit physics ───────────────────────────────────────────")

# Settling at STP (T=293 K, p=101325 Pa), unit density spheres — textbook:
# 1 µm ≈ 3.5e-5, 10 µm ≈ 3.1e-3, 100 µm ≈ 0.25–0.31 m/s (transition regime)
vs1   = EASLEE.settling_velocity(1.0f-6,  1000.0f0, 293.15f0, 101325.0f0)
vs10  = EASLEE.settling_velocity(1.0f-5,  1000.0f0, 293.15f0, 101325.0f0)
vs100 = EASLEE.settling_velocity(1.0f-4,  1000.0f0, 293.15f0, 101325.0f0)
check("v_settle 1 µm",   2.8f-5 < vs1   < 4.2f-5, @sprintf("%.2e m/s (ref 3.5e-5)", vs1))
check("v_settle 10 µm",  2.6f-3 < vs10  < 3.6f-3, @sprintf("%.2e m/s (ref 3.1e-3)", vs10))
check("v_settle 100 µm", 0.20f0 < vs100 < 0.35f0, @sprintf("%.2f m/s (ref 0.25–0.31, transition)", vs100))
# Altitude dependence: slip correction → faster at low pressure
vs1_hi = EASLEE.settling_velocity(1.0f-6, 1000.0f0, 220.0f0, 5000.0f0)
check("slip correction at 5 kPa", vs1_hi > 3.0f0 * vs1, @sprintf("%.2e vs %.2e m/s", vs1_hi, vs1))

# Mie: large-sphere limit Qext → 2; strong-interaction peak ~3–4.5 near x≈4–6
qbig  = EASLEE.mie_qext(300.0, 1.53 - 0.003im)
qpeak = maximum(EASLEE.mie_qext(x, 1.33 - 1e-9im) for x in 4.0:0.25:7.0)
qray  = EASLEE.mie_qext(0.05, 1.53 - 0.003im)
check("Mie Qext large-x limit", 1.9 < qbig < 2.2, @sprintf("Qext(x=300)=%.3f", qbig))
check("Mie interference peak",  3.0 < qpeak < 4.6, @sprintf("max Qext=%.2f", qpeak))
check("Mie Rayleigh regime",    0.0 < qray < 0.01, @sprintf("Qext(x=0.05)=%.1e", qray))

# Dry deposition: vd minimum in the accumulation mode, large for coarse
vd(d) = EASLEE.deposition_velocity(Float32(d), 2000.0f0, 290.0f0, 101000.0f0,
                                   0.4f0, -100.0f0, 30.0f0)
vd_fine, vd_acc, vd_coarse = vd(1e-8), vd(3e-7), vd(2e-5)
check("vd(10 nm) > vd(0.3 µm)", vd_fine > 2vd_acc, @sprintf("%.1e vs %.1e m/s", vd_fine, vd_acc))
check("vd(20 µm) > vd(0.3 µm)", vd_coarse > 10vd_acc, @sprintf("%.1e vs %.1e m/s", vd_coarse, vd_acc))

# Lognormal bins: mass fractions sum to 1, Qext per bin finite
bins = EASLEE.AerosolBins(dg=2.0e-6, sigma_g=2.0, nbins=8, rho_p=2600.0, species=:dust)
check("lognormal bins", isapprox(sum(bins.mass_frac), 1.0, atol=1e-5) &&
      all(isfinite, bins.qext) && length(bins.d) == 8,
      @sprintf("d = %.2g–%.2g m, Qext = %.2f–%.2f", bins.d[1], bins.d[end],
               minimum(bins.qext), maximum(bins.qext)))
ebins = EASLEE.AerosolBins([0.1e-6, 1e-6, 10e-6, 100e-6]; species=:ash, rho_p=2400.0)
check("explicit bins", length(ebins.d) == 3 && isapprox(sum(ebins.mass_frac), 1.0, atol=1e-5),
      @sprintf("%d bins", length(ebins.d)))

# ── Hygroscopic growth (κ-Köhler) ─────────────────────────────────────────────
gf_sulf = EASLEE.growth_factor(0.61f0, 0.90f0)
gf_salt = EASLEE.growth_factor(1.28f0, 0.90f0)
gf_dust = EASLEE.growth_factor(0.03f0, 0.90f0)
check("GF sulfate @90% RH", isapprox(gf_sulf, (1 + 0.61 * 9)^(1/3); rtol=1e-3),
      @sprintf("%.3f (analytic 1.866)", gf_sulf))
check("GF ordering seasalt > sulfate > dust",
      gf_salt > gf_sulf > gf_dust > 1.0,
      @sprintf("%.2f > %.2f > %.3f", gf_salt, gf_sulf, gf_dust))
# Optics humidification: per-bin wet cross-section ratio at RH=0.9 vs dry
sbins = EASLEE.AerosolBins(dg=0.3e-6, sigma_g=1.8, nbins=8, rho_p=1770.0, species=:sulfate)
jr90  = findfirst(==(0.90f0), EASLEE.AEROSOL_RH_GRID)
bmid  = 4
fRH   = (sbins.qext[bmid, jr90] * gf_sulf^2) / max(sbins.qext[bmid, 1], 1e-6)
check("optics humidification f(90%) for sulfate", 1.5 < fRH < 12.0,
      @sprintf("σ_ext ratio %.2f at d=%.2g m", fRH, sbins.d[bmid]))
# Settling speeds up when wet (d², partially offset by lower density)
vs_dry = EASLEE.settling_velocity(1.0f-6, 1770.0f0, 290.0f0, 101000.0f0)
vs_wet = EASLEE.settling_velocity(1.0f-6 * gf_sulf, 1770.0f0/gf_sulf^3 + 1000.0f0*(1-1/gf_sulf^3),
                                  290.0f0, 101000.0f0)
check("wet settling > dry settling (1 µm sulfate @90%)", vs_wet > 1.5f0 * vs_dry,
      @sprintf("%.2e vs %.2e m/s", vs_wet, vs_dry))

println("\n── Part 2: 3 h ERA5 integration ───────────────────────────────────")
cfg = isempty(ARGS) ? joinpath(@__DIR__,"..","configs","example_dust.toml") : ARGS[1]
spec = read_runspec(cfg)
rt, _obj = build_run(spec; prefetch=false, extra_slots=10_000)
ctx   = rt.ctx
state = rt.state
rc    = ctx.runconfig
for n in 1:length(ctx.nests)
    ctx.nests[n].nestrunconfig.turbulence = EASLEENoTurb()
end
convfields = [EASLEE.ConvFields(rt.metfields3D[n]) for n in 1:rt.nest]
aero = EASLEE.AerosolModel(length(state.status), ctx, bins)

read_met_step!(rt, rc.startdate)

# Settling check in a quiet column: 60 µm particles at ~5 km must fall ≈ v_t·t
t_array = MetAdapter.met_time_index_array(rc.startdate, ctx.nests)
m3 = ctx.nests[1].metfields3D
nz = size(m3.zlev, 3)
zprof_q = [m3.zlev[100, 240, k, 1] for k in 1:nz]
k5 = findfirst(k -> zprof_q[k] > 5000.0f0, 1:nz)
big = EASLEE.AerosolBins([5.0e-5, 7.0e-5]; species=:dust, rho_p=2600.0)  # one ~60 µm bin
aero_big = EASLEE.AerosolModel(length(state.status), ctx, big)
ntest = 200
for p in 1:ntest
    state.status[p] = 0x01; state.nest[p] = 1
    state.indx[p] = 100.0f0; state.indy[p] = 240.0f0
    state.indl[p] = Float32(k5); state.mass[p] = 1.0f0
    aero_big.bin[p] = 1
end
for s in 1:4
    EASLEE.apply_settling!(aero_big, state, ctx, 900.0, t_array)
end
zend = let k = clamp(floor(Int, state.indl[1]), 1, nz-1)
    zprof_q[k] + (state.indl[1]-k)*(zprof_q[k+1]-zprof_q[k])
end
Tq = m3.T[100, 240, k5, 1]; pq = m3.puv[100, 240, k5, 1]
vt_ref = EASLEE.settling_velocity(big.d[1], big.rho_p, Tq, pq)
drop = zprof_q[k5] - zend
check("settling drop ≈ v_t·t (60 µm, 1 h)",
      0.7 * vt_ref * 3600 < drop < 1.3 * vt_ref * 3600,
      @sprintf("drop %.0f m vs v_t·t = %.0f m", drop, vt_ref * 3600))
for p in 1:ntest
    state.status[p] = 0x00; state.mass[p] = 0.0f0; aero_big.bin[p] = 0
end

# Full run: release (from RELEASES), advect + settle + dry dep + wet dep
ts = rc.timestep
nsteps = round(Int, abs(Dates.value(rc.enddate - rc.startdate)) / abs(Dates.value(ts)))
for itime in 0:nsteps
    currentdate = rc.startdate + itime * ts
    read_met_step!(rt, currentdate)
    EASLEE.step_model!(state, ctx, currentdate;
                       convfields = convfields, aerosol = aero)
    State.sync_dtime_ns!(state)
end

active_mass = sum(state.mass[i] for i in eachindex(state.mass) if state.status[i] == 0x01; init=0.0)
dd = sum(sum.(aero.drydep))
wd = sum(sum.(aero.wetdep))
@printf("\n  final active %.4f + drydep %.4f + wetdep %.4f = %.5f\n",
        active_mass, dd, wd, active_mass + dd + wd)
check("mass budget closes", isapprox(active_mass + dd + wd, 1.0, atol=2e-3),
      @sprintf("sum = %.5f", active_mass + dd + wd))
check("wet deposition active", wd > 1e-4, @sprintf("wetdep = %.4f", wd))

t_end = MetAdapter.met_time_index_array(rc.enddate, ctx.nests)
aod = EASLEE.compute_aod(aero, state, ctx, t_end)
check("AOD finite & positive", all(isfinite, aod) && maximum(aod) > 0,
      @sprintf("max AOD(550nm) = %.3g, %d cells > 0", maximum(aod), count(>(0), aod)))

# Per-bin occupancy ≈ mass fractions
cnt = zeros(Int, length(bins.d))
for i in eachindex(state.status)
    state.status[i] == 0x01 && aero.bin[i] > 0 && (cnt[aero.bin[i]] += 1)
end
occ = cnt ./ max(sum(cnt), 1)
check("bin occupancy ≈ mass fractions", maximum(abs.(occ .- bins.mass_frac)) < 0.02,
      @sprintf("max |Δ| = %.3f", maximum(abs.(occ .- bins.mass_frac))))

println()
@printf("=== %d passed, %d failed ===\n", npass, nfail)
exit(nfail == 0 ? 0 : 1)
