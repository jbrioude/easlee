# test_dust_kok.jl — Kok et al. (2014) dust emission + Kok (2011) size distribution:
#   [1] Shao & Lu fluid threshold ≈ 0.21 m/s at 1.225 kg/m³, D = 127 µm; decreases with ρ_a
#   [2] Fécan moisture factor: 1 below w′, > 1 above, increasing
#   [3] flux: 0 at/below threshold, positive above, increasing in u*; ∝ clay and bare;
#       Kok 2014 magnitude at u* = 0.5 m/s in the 1e-6–1e-4 kg m⁻² s⁻¹ range (their Fig. 4)
#   [4] Kok 2011 distribution: fractions sum to 1, peak of dV/dlnD between 2 and 10 µm,
#       mass fraction below 2.5 µm small (Kok 2011: ~ a few %), coarse half dominant
#   [5] AerosolBins with mass_frac = kok fractions; "kok" species through the RunSpec;
#       example config parses with kind = "dust" and the [release.dust] keys at RunSpec level.
#   julia --project=easlee test/test_dust_kok.jl
using Printf, Dates
global gpu = false; global chem = false
const REPO = joinpath(@__DIR__, "..")
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
const O = EASLEE.Operators
npass = 0; nfail = 0
function check(label, ok, detail="")
    global npass, nfail; ok ? (npass += 1) : (nfail += 1)
    @printf("  %-66s %s   %s\n", label, ok ? "PASS" : "FAIL", detail)
end
ut0 = O.dust_threshold_dry(1.225)
check("[1] fluid threshold 0.19–0.23 m/s at ρ0 (Kok 2014: ≈0.2)", 0.19 < ut0 < 0.23, @sprintf("%.3f", ut0))
check("[1] threshold decreases with air density", O.dust_threshold_dry(1.0) > ut0 > O.dust_threshold_dry(1.4))
check("[2] Fécan: 1 below w′ (clay 0.2 → w′ = 3.96 %)", O.fecan_moisture_factor(3.0, 0.2) == 1.0)
check("[2] Fécan: increasing above w′", 1.0 < O.fecan_moisture_factor(6.0, 0.2) < O.fecan_moisture_factor(12.0, 0.2))
check("[3] flux 0 at and below threshold", O.dust_flux_kok(ut0, 1.225) == 0.0 && O.dust_flux_kok(0.1, 1.225) == 0.0)
F3 = O.dust_flux_kok(0.3, 1.225); F5 = O.dust_flux_kok(0.5, 1.225); F7 = O.dust_flux_kok(0.7, 1.225)
check("[3] flux increasing in u*", 0 < F3 < F5 < F7, @sprintf("%.2e %.2e %.2e", F3, F5, F7))
check("[3] flux at u* = 0.5 (clay 0.2) in 1e-6–1e-4 kg m⁻² s⁻¹", 1e-6 < F5 < 1e-4, @sprintf("%.2e", F5))
check("[3] flux ∝ clay and bare", abs(O.dust_flux_kok(0.5, 1.225; clay=0.1) / F5 - 0.5) < 1e-9 && abs(O.dust_flux_kok(0.5, 1.225; bare=0.25) / F5 - 0.25) < 1e-9)
check("[3] soil moisture reduces the flux", O.dust_flux_kok(0.5, 1.225; w_pct=10.0) < F5)
check("[3] the dust-coefficient term: flux ratio F(0.7)/F(0.5) > (0.7/0.5)² (u*_st-dependent exponent)", F7 / F5 > (0.7 / 0.5)^2)
edges = exp.(range(log(0.2e-6), log(20e-6), length=9)); mf = O.kok_mass_fractions(edges)
check("[4] Kok 2011 fractions sum to 1, all ≥ 0", abs(sum(mf) - 1) < 1e-12 && all(mf .>= 0))
Ds = exp.(range(log(0.2e-6), log(20e-6), length=400)); pk = Ds[argmax(O.kok_dvdlnd.(Ds))]
check("[4] peak of dV/dlnD between 2 and 10 µm", 2e-6 < pk < 10e-6, @sprintf("%.1f µm", pk*1e6))
f25 = O.kok_mass_fractions([0.2e-6, 2.5e-6, 20e-6])[1]
check("[4] mass fraction below 2.5 µm between 2 and 12 %", 0.02 < f25 < 0.12, @sprintf("%.1f %%", 100*f25))
check("[4] coarse half (>2 µm) carries most of the mass", O.kok_mass_fractions([0.2e-6, 2e-6, 20e-6])[2] > 0.8)
bins = EASLEE.AerosolBins(edges; mass_frac=mf, rho_p=2650.0, species=:dust)
check("[5] AerosolBins accepts the Kok fractions (8 bins)", length(bins.d) == 8 && abs(sum(bins.mass_frac) - 1) < 1e-6)
s = read_runspec(joinpath(REPO, "configs", "example_dust_kok.toml"))
r = s.release[1]; sp = s.species[1]
check("[5] example config: kind dust, keys at RunSpec level", r.kind == "dust" && r.dust_clay == 0.15 && r.dust_w_pct == 0.5 && r.mass_per_parcel_kg == 20000.0 && r.particles == 400000, "clay=$(r.dust_clay) w=$(r.dust_w_pct) q=$(r.mass_per_parcel_kg) cap=$(r.particles)")
check("[5] example config: species size.distribution = kok", sp.size_distribution == "kok" && sp.nbins == 8)
errs = try; validate_runspec(s; check_paths=false); String[]; catch e; [sprint(showerror, e)]; end
check("[5] example config validates", isempty(errs), isempty(errs) ? "" : errs[1][1:min(120, end)])
d2 = runspec_to_dict(s); check("[5] round-trip keeps [release.dust] and size kok", d2["release"][1]["dust"]["clay"] == 0.15 && d2["species"][1]["size"]["distribution"] == "kok")
println("\n$(npass) passed, $(nfail) failed"); nfail == 0 || exit(1)
