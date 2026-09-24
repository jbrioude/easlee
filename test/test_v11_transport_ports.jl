# test_v11_transport_ports.jl — pins the transport / met-reading fixes ported
# from the developer tree (v1.1) into this distribution on 2026-09-13.  The
# chemistry-side changes of v1.1 are NOT part of this distribution and are not
# tested here.
#   [1] tke_partition_h: Hanna-fill closure yields no horizontal variance;
#       native-e closure partitions e so that σ_w² + 2σ_h² = 2e EXACTLY and
#       τ_Lh is clamped to [TLW_FLOOR, 600 s].
#   [2] the horizontal OU step reproduces var(u') = σ_h² (stationary variance).
#   [3] esat_mk / esat_liq_mk: Murphy–Koop over liquid ≈ 611 Pa at 0 °C, ice
#       branch below 235 K, and the Bolton fit the plume rise used to carry sits
#       1.7–1.9× above e_si at 210 K (the documented reason for the swap).
#   [4] free-slot cursor `_first_unused` + `_release_nsteps` (share review B2c).
#   [5] [physics.turbulence].horizontal_tke and [output.trajectory].format parse,
#       reach the RunSpec, and survive a TOML round trip.
#   [6] NestRunConfig carries horiz_tke / vmix_*, and the CPU-only switch names
#       horizontal_tke (so a KA/GPU request cannot silently drop it).
#   [7] traj_open: "full" writes the header, "demovor" writes none.
# Synthetic data only — no met archive needed.

const REPO = abspath(joinpath(@__DIR__, ".."))

using Dates, Random, Statistics, TOML

include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
using .EasleeConfig.ConfigTOML: unknown_config_keys
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
const OP = EASLEE.Operators

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name, isempty(info) ? "" : "  ($info)")
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
println("="^78); println("[1] tke_partition_h"); println("="^78)
hanna  = OP.VMixClosure(false, NaN32, NaN32)
native = OP.VMixClosure(true, 0.4f0, 0.8f0)     # r_w = 0.4 ⇒ r_h = 0.3 each
e = 1.5f0; yl = 20.0f0
tlh0, sig0 = OP.tke_partition_h(hanna, e, yl)
check("Hanna-fill driver → (0, 0): no horizontal variance to partition", tlh0 == 0.0f0 && sig0 == 0.0f0)
tlh, sigh = OP.tke_partition_h(native, e, yl)
sigw = sqrt(2.0f0 * e * native.r_w)
check("native e: σ_w² + 2σ_h² = 2e exactly (r_w + 2r_h = 1)",
      isapprox(sigw^2 + 2sigh^2, 2e; rtol=1e-5), info="$(sigw^2 + 2sigh^2) vs $(2e)")
check("σ_h = √(2·e·r_h)", isapprox(sigh, sqrt(2e * 0.5f0 * (1 - native.r_w)); rtol=1e-6))
kh = native.c_k * yl * sqrt(e)
check("τ_Lh = K_h/σ_h² with K_h = C_K·ℓ·√e", isapprox(tlh, kh / sigh^2; rtol=1e-5), info="τ_Lh=$(tlh)s")
tlh_big, _ = OP.tke_partition_h(native, 1.0f-4, 1.0f6)     # K_h huge, σ_h tiny
check("τ_Lh capped at 600 s", tlh_big == 600.0f0, info="$(tlh_big)")
tlh_small, _ = OP.tke_partition_h(native, 50.0f0, 1.0f-3)  # K_h tiny
check("τ_Lh floored at TLW_FLOOR", tlh_small == OP.TLW_FLOOR, info="$(tlh_small)")

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[2] horizontal OU step: stationary variance = σ_h²"); println("="^78)
Random.seed!(11)
dt = 30.0f0; rr = exp(-dt / tlh); sq = sqrt(max(1.0f0 - rr * rr, 0.0f0))
n = 200_000; u = 0.0f0; acc = Vector{Float32}(undef, n)
for k in 1:n
    global u = u * rr + sigh * sq * randn(Float32)
    acc[k] = u
end
v = var(acc[20_000:end])
check("var(u') → σ_h² after spin-up (within 3 %)", isapprox(v, sigh^2; rtol=0.03), info="var=$(round(v,digits=4)) σ_h²=$(round(sigh^2,digits=4))")
check("mean(u') → 0", abs(mean(acc[20_000:end])) < 0.02 * sigh)

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[3] esat_mk / esat_liq_mk"); println("="^78)
check("liquid MK at 273.15 K ≈ 611 Pa", isapprox(OP.esat_liq_mk(273.15), 611.2; rtol=0.01), info="$(round(OP.esat_liq_mk(273.15),digits=1))")
check("esat_mk below 235 K is the ice law", OP.esat_mk(220.0) == OP.esat_ice(220.0))
check("esat_mk above 235 K is the liquid law", OP.esat_mk(250.0) == OP.esat_liq_mk(250.0))
r235 = OP.esat_liq_mk(235.0) / OP.esat_ice(235.0)
check("at the 235 K switch the two laws differ by the physical liquid/ice ratio (1.3–1.6)", 1.3 < r235 < 1.6, info="ratio $(round(r235,digits=3))")
bolton(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
rb = bolton(210.0) / OP.esat_ice(210.0)
check("Bolton liquid fit sits 1.7–1.9× above e_si at 210 K (why plume_rise swapped)", 1.7 < rb < 1.9, info="ratio $(round(rb,digits=2))")
check("esat_mk is monotone increasing 180–300 K", all(diff([OP.esat_mk(T) for T in 180.0:1.0:300.0]) .> 0))

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[4] free-slot cursor + release step count"); println("="^78)
st = EASLEE.State.ParcelState(10)
OP._FREE_CURSOR[] = 1
st.status[1:5] .= 0x01
a = OP._first_unused(st, 2); st.status[a] .= 0x01
check("first request from a fresh cursor → [6,7]", a == [6, 7] && OP._FREE_CURSOR[] == 8)
b = OP._first_unused(st, 2); st.status[b] .= 0x01
check("cursor continues → [8,9] (no rescan from 1)", b == [8, 9])
st.status[2] = 0x00
c = OP._first_unused(st, 2); st.status[c] .= 0x01
check("tail exhausted → wrap finds the slot freed below the cursor", c == [10, 2])
d = OP._first_unused(st, 1)
check("pool full → empty, cursor reset to 1", isempty(d) && OP._FREE_CURSOR[] == 1)
st.status[4] = 0x00
check("after the reset a freed slot is found again", OP._first_unused(st, 1) == [4])
OP._FREE_CURSOR[] = 99
check("cursor beyond the pool is clamped", OP._first_unused(st, 1) == [4])
OP._FREE_CURSOR[] = 1
tloops = DateTime(2022,1,1):Minute(15):DateTime(2022,1,2)
Random.seed!(7); okn = true
for _ in 1:50
    t1 = DateTime(2022,1,1) + Minute(rand(-60:1500)); t2 = t1 + Minute(rand(0:600))
    global okn &= OP._release_nsteps(tloops, t1, t2) == sum(t1 .<= tloops .<= t2)
end
check("_release_nsteps == sum(t1 .<= tloops .<= t2) on 50 random windows", okn)

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[5] config keys reach the RunSpec and round-trip"); println("="^78)
tmp = mktempdir()
cfg = joinpath(tmp, "t.toml")
write(cfg, """
name = "ports"
[run]
start = 2022-01-15T04:00:00
stop  = 2022-01-15T06:00:00
timestep = "15min"
[[meteorology]]
model = "ERA5"
path = "/nonexistent"
variables = "ERA5"
[[species]]
id = "TRACER"
kind = "gas"
[[release]]
species = "TRACER"
lon = [0.0, 1.0]
lat = [0.0, 1.0]
z = ["100m", "200m"]
time = [2022-01-15T04:00:00, 2022-01-15T05:00:00]
particles = 10
mass = "1kg"
[output]
path = "$(tmp)/out/"
interval = "1h"
sample = "15min"
  [output.grid]
  lon0 = -10.0
  lat0 = -10.0
  nx = 10
  ny = 10
  dx = 1.0
  dy = 1.0
  heights = [0, 1000]
  [output.trajectory]
  enabled = true
  format  = "demovor"
[physics]
  [physics.turbulence]
  scheme = "thomson"
  horizontal_tke = true
""")
spec = read_runspec(cfg; check_paths=false)
check("horizontal_tke parses to the RunSpec", spec.physics.turbulence.horizontal_tke === true)
check("trajectory.format parses to the RunSpec", spec.output.trajectory.format == "demovor")
tv = runspec_to_dict(spec)
check("horizontal_tke survives the TOML round trip", get(get(get(tv, "physics", Dict()), "turbulence", Dict()), "horizontal_tke", false) === true)
check("both keys are KNOWN to the schema (strict_keys would not reject them)",
      isempty(filter(p -> occursin("horizontal_tke", p[1]) || occursin("trajectory.format", p[1]),
                     unknown_config_keys(TOML.parsefile(cfg)))))

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[6] NestRunConfig fields + CPU-only switch"); println("="^78)
NRC = EASLEE.Types.NestRunConfig
check("NestRunConfig has horiz_tke", :horiz_tke in fieldnames(NRC))
check("NestRunConfig has the per-nest closure (vmix_native, vmix_rw, vmix_ck)",
      all(f -> f in fieldnames(NRC), (:vmix_native, :vmix_rw, :vmix_ck)))
# 2026-09-23: the KA kernel carries horizontal_tke; what must remain is the explicit-backend refusal
src = read(joinpath(REPO, "easlee", "src", "core", "operators.jl"), String)
check("advance_dispatch! refuses an explicitly requested backend it cannot honour (2026-09-21: horizontal_tke runs on the KA path, so it is no longer listed as threads-only)", occursin("_ADVANCE_REQUESTED", src) && occursin("Refusing to downgrade", src) && !occursin("horizontal_tke (TKE-driven horizontal K_h)", src))
bsrc = read(joinpath(REPO, "easlee", "build_run.jl"), String)
check("build_run refuses horizontal_tke with [physics.mixing], scheme=none or a Hanna-fill met (the KA refusal is retired: the kernel carries the term since 2026-09-21)",
      all(s -> occursin(s, bsrc), ("horizontal_tke and [physics.mixing] are both armed",
                                   "needs a turbulence scheme", "requires a native-e driver")) &&
      !occursin("horizontal_tke is implemented on the CPU advance only", bsrc))

# ═══════════════════════════════════════════════════════════════════════════════
println(); println("="^78); println("[7] traj_open header by format"); println("="^78)
st2 = EASLEE.State.ParcelState(20); st2.status[1:10] .= 0x01
tw = EASLEE.traj_open(spec, st2, 300.0)
EASLEE.traj_close(tw)
p = joinpath(spec.output.path, "traj_ports.csv")
check("demovor: no header line written", isfile(p) && filesize(p) == 0)
spec.output.trajectory.format = "full"
tw = EASLEE.traj_open(spec, st2, 300.0); EASLEE.traj_close(tw)
hdr = readline(p)
check("full: header is time,id,lon,lat,height_m,nest,mass_kg,age_s,rho_kgm3,pblh_m",
      hdr == "time,id,lon,lat,height_m,nest,mass_kg,age_s,rho_kgm3,pblh_m", info=hdr)

println(); println("$npass passed, $nfail failed")
nfail == 0 || exit(1)
