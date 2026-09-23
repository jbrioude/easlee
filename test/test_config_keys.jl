#!/usr/bin/env julia
# test_config_keys.jl — configuration layer (review 2026-09-04 A7):
#   [a] a misspelled key is refused; the message names the key and suggests the nearest schema key,
#   [b] strict_keys = false downgrades it to a warning and the config loads,
#   [c] every shipped configs/*.toml has zero unknown keys,
#   [d] [species.settling] and the legacy [species.fall] alias both populate settle_r_eff,
#   [e] species.molar_mass parses (auto convert mass_ratio no longer falls back to 1),
#   [f] physics.ctl = 2.5 is refused, ctl = 5 accepted,
#   [g] two aerosol species with different bins are refused (one shared AerosolBins),
#   [h] refractive_index parses / is validated.
#
#   julia --project=easlee test/test_config_keys.jl
using Printf, Dates, Logging
import TOML

const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
using .EasleeConfig.ConfigTOML: unknown_config_keys, check_config_keys
using .EasleeConfig.RunSpecs: parse_refractive_index

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

const TMP = mktempdir()
const BASE = """
name = "keys_test"
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
path = "$(TMP)/out/"
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
"""
function writecfg(name, extra; head="")
    p = joinpath(TMP, name * ".toml")
    write(p, head * BASE * extra)
    return p
end
errmsg(f) = try f(); "" catch e; sprint(showerror, e) end

println("="^78); println("[a] misspelled key is refused with the key + a suggestion"); println("="^78)
pa = writecfg("bad_kh", "[physics.mixing]\nscheme = \"pairwise\"\nKh = 500.0\n")
m = errmsg(() -> read_runspec(pa; check_paths=false))
check("read_runspec errors", !isempty(m))
check("message names physics.mixing.Kh", occursin("physics.mixing.Kh", m), info=first(m, 0))
check("message suggests physics.mixing.K_h", occursin("physics.mixing.K_h", m))
u = unknown_config_keys(TOML.parsefile(pa))
check("unknown_config_keys returns exactly the bad key", length(u) == 1 && u[1][1] == "physics.mixing.Kh" && u[1][2] == "physics.mixing.K_h", info=string(u))
pb = writecfg("bad_species", "[[species]]\nid = \"X\"\nkind = \"gas\"\nmolar_mas = 64.0\n[species.convert]\nto = \"TRACER\"\ntau = \"1d\"\n")
m = errmsg(() -> read_runspec(pb; check_paths=false))
check("array-of-tables key: species.molar_mas → species.molar_mass",
      occursin("species.molar_mas ", m) && occursin("species.molar_mass", m))

println("="^78); println("[b] strict_keys = false → warning, config loads"); println("="^78)
pw = writecfg("warn_kh", "[physics.mixing]\nscheme = \"pairwise\"\nKh = 500.0\n"; head="strict_keys = false\n")
spec = nothing
logged = let io = IOBuffer()          # capture the warning through a temporary logger
    lg = Logging.SimpleLogger(io, Logging.Warn)
    global spec = Logging.with_logger(lg) do
        read_runspec(pw; check_paths=false)
    end
    String(take!(io))
end
check("loads with strict_keys = false", spec !== nothing && spec.name == "keys_test")
check("a warning naming the key was emitted", occursin("physics.mixing.Kh", logged))
check("strict_keys itself is not an unknown key", isempty(unknown_config_keys(Dict("strict_keys" => false))))

println("="^78); println("[c] every shipped config has zero unknown keys"); println("="^78)
nbad = 0
for f in sort(readdir(joinpath(REPO, "configs"); join=true))
    endswith(f, ".toml") || continue
    uu = unknown_config_keys(TOML.parsefile(f))
    isempty(uu) || (global nbad += 1; println("    unknown in ", basename(f), ": ", uu))
end
check("configs/*.toml: zero unknown keys", nbad == 0, info="$(nbad) bad")

println("="^78); println("[d] [species.settling] and [species.fall] both populate settle_r_eff"); println("="^78)
ps = writecfg("settling", "[[species]]\nid = \"ICE\"\nkind = \"gas\"\n[species.settling]\nr_eff = \"5um\"\nrho = 900.0\n")
sp = read_runspec(ps; check_paths=false)
ice = sp.species[findfirst(s -> s.id == "ICE", sp.species)]
check("[species.settling] r_eff → 5e-6 m, rho 900", isapprox(ice.settle_r_eff, 5e-6) && ice.settle_rho == 900.0)
check("gas-kind species: settling bool stays true (table form)", ice.settling)
pf = writecfg("fall", "[[species]]\nid = \"ICE\"\nkind = \"gas\"\n[species.fall]\nr_eff = \"5um\"\n")
sp = read_runspec(pf; check_paths=false)
ice = sp.species[findfirst(s -> s.id == "ICE", sp.species)]
check("[species.fall] alias → same settle_r_eff, default rho 917", isapprox(ice.settle_r_eff, 5e-6) && ice.settle_rho == 917.0)
d = runspec_to_dict(sp)
e = d["species"][findfirst(x -> x["id"] == "ICE", d["species"])]
check("writer emits [species.settling] (not fall)", haskey(e, "settling") && e["settling"] isa AbstractDict && isapprox(e["settling"]["r_eff"], 5e-6))
rt = dict_to_runspec(d)
check("round-trip keeps settle_r_eff", isapprox(rt.species[findfirst(s -> s.id == "ICE", rt.species)].settle_r_eff, 5e-6))
pa2 = writecfg("aero_settle_table", "[[species]]\nid = \"ASH\"\nkind = \"aerosol\"\n[species.settling]\nenabled = false\nr_eff = \"50um\"\n")
sp = read_runspec(pa2; check_paths=false)
ash = sp.species[findfirst(s -> s.id == "ASH", sp.species)]
check("aerosol: settling table with enabled=false → settling false, r_eff set", !ash.settling && isapprox(ash.settle_r_eff, 50e-6))

println("="^78); println("[e] molar_mass parses"); println("="^78)
pm = writecfg("molar", "[[species]]\nid = \"SO2x\"\nkind = \"gas\"\nmolar_mass = 64.06\n[species.convert]\nto = \"SO4x\"\ntau = \"25d\"\n[[species]]\nid = \"SO4x\"\nkind = \"aerosol\"\naerosol_type = \"sulfate\"\nmolar_mass = 96.06\n")
sp = read_runspec(pm; check_paths=false)
so2 = sp.species[findfirst(s -> s.id == "SO2x", sp.species)]
so4 = sp.species[findfirst(s -> s.id == "SO4x", sp.species)]
check("molar_mass 64.06 / 96.06 parsed (ids not in the species DB)", so2.molar_mass == 64.06 && so4.molar_mass == 96.06)
check("auto mass_ratio M_to/M_from available", so2.convert_mass_ratio == 0.0 && isapprox(so4.molar_mass / so2.molar_mass, 1.4995; atol=1e-3))
check("writer emits molar_mass", any(x -> get(x, "molar_mass", 0.0) == 64.06, runspec_to_dict(sp)["species"]))
m = errmsg(() -> read_runspec(writecfg("molar_neg", "[[species]]\nid = \"Y\"\nkind = \"gas\"\nmolar_mass = -1.0\n"); check_paths=false))
check("negative molar_mass refused", occursin("molar_mass", m))

println("="^78); println("[f] physics.ctl integer check"); println("="^78)
m = errmsg(() -> read_runspec(writecfg("ctl25", "[physics]\nctl = 2.5\n"); check_paths=false))
check("ctl = 2.5 refused", occursin("physics.ctl", m) && occursin("2.5", m))
sp = read_runspec(writecfg("ctl5", "[physics]\nctl = 5\n"); check_paths=false)
check("ctl = 5 accepted", sp.physics.ctl == 5.0)
sp = read_runspec(writecfg("ctl10f", "[physics]\nctl = 10.0\n"); check_paths=false)
check("ctl = 10.0 accepted (integer-valued float)", sp.physics.ctl == 10.0)
m = errmsg(() -> read_runspec(writecfg("ctl0", "[physics]\nctl = 0\n"); check_paths=false))
check("ctl = 0 refused", occursin("physics.ctl", m))

println("="^78); println("[g] two aerosol species with different bins are refused"); println("="^78)
two = "[[species]]\nid = \"A1\"\nkind = \"aerosol\"\naerosol_type = \"dust\"\ndensity = 2600.0\n[species.size]\ndg = \"1um\"\nsigma_g = 2.0\nnbins = 8\n" *
      "[[species]]\nid = \"A2\"\nkind = \"aerosol\"\naerosol_type = \"dust\"\ndensity = 2600.0\n[species.size]\ndg = \"4um\"\nsigma_g = 2.0\nnbins = 8\n"
m = errmsg(() -> read_runspec(writecfg("two_aero_diff", two); check_paths=false))
check("different dg → refused, message names both ids and AerosolBins", occursin("\"A1\"", m) && occursin("\"A2\"", m) && occursin("AerosolBins", m))
same = replace(two, "dg = \"4um\"" => "dg = \"1um\"")
sp = read_runspec(writecfg("two_aero_same", same); check_paths=false)
check("identical bins/density/type/kappa → accepted", count(s -> s.kind == "aerosol", sp.species) == 2)
m = errmsg(() -> read_runspec(writecfg("two_aero_type", replace(same, "id = \"A2\"\nkind = \"aerosol\"\naerosol_type = \"dust\"" => "id = \"A2\"\nkind = \"aerosol\"\naerosol_type = \"ash\"")); check_paths=false))
check("different aerosol_type → refused", occursin("AerosolBins", m))

println("="^78); println("[h] refractive_index"); println("="^78)
check("parse \"1.53-0.003i\" → 1.53 - 0.003im", parse_refractive_index("1.53-0.003i") == ComplexF64(1.53, -0.003))
check("parse \"1.43-0.0i\"", parse_refractive_index("1.43-0.0i") == ComplexF64(1.43, 0.0))
check("parse \"auto\" → nothing", parse_refractive_index("auto") === nothing)
check("parse garbage → nothing", parse_refractive_index("glass") === nothing)
m = errmsg(() -> read_runspec(writecfg("ri_bad", "[[species]]\nid = \"A1\"\nkind = \"aerosol\"\nrefractive_index = \"glass\"\n"); check_paths=false))
check("bad refractive_index refused at validation", occursin("refractive_index", m))

println("="^78); println("[i] TOML-only spellings are known (wildcards, tables, scalar sugar)"); println("="^78)
sugar = """
[physics]
turbulence = "hanna"
[physics.nudge]
z = [18000.0, 40000.0]
[output.grid]
origin = [-10.0, -10.0]
spacing = [1.0, 1.0]
size = [10, 10]
[output.receptors]
[[output.receptors.station]]
id = "S1"
lon = 1.0
lat = 1.0
[[species]]
id = "H2O_v"
kind = "gas"
ice = true
[[species]]
id = "ASH"
kind = "aerosol"
[species.settling]
r_eff = "50um"
[[release]]
species = "H2O_v"
kind = "volcanic"
lon = 1.0
lat = 1.0
time = [2022-01-15T04:00:00, 2022-01-15T05:00:00]
[release.volcanic]
phases = [[2022-01-15T04:00:00, 2022-01-15T05:00:00, "20km", "30km", 1.0]]
[release.plume_rise]
scheme = "volcanic"
mode = "parcels"
[release.plume_rise.vent_mix]
magma_t = 1300.0
[release.chem_init]
so2_ppb = 100.0
HCl = 5.0
h2o_from_plume = true
[chemistry]
enabled = true
[chemistry.init]
default = "cams"
[chemistry.init.sources.cams]
kind = "grib"
dir = "/x"
[chemistry.init.sources.cams.vars]
O3 = "go3"
[[chemistry.init.species]]
name = "O3"
from = "cams"
[chemistry.aerosols]
enabled = true
"""
u = unknown_config_keys(TOML.parse(sugar))
check("wildcard / table / vector / sugar spellings: zero unknown keys", isempty(u), info=string(u))
u = unknown_config_keys(TOML.parse(replace(sugar, "dir = \"/x\"\n" => "dir = \"/x\"\nunit = \"ppb\"\n")))
check("…but a typo inside a keyed sub-table is still caught", length(u) == 1 && u[1][1] == "chemistry.init.sources.cams.unit", info=string(u))

println("\n", npass, " passed, ", nfail, " failed")
nfail == 0 || (println("FAIL"); exit(1))
