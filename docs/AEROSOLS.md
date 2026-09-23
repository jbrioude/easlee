# EASLEE — Size-resolved aerosol transport

Status: 2026-06-12.  Implementation in `easlee/src/physics/aerosol.jl`;
validation in `validation/test_aerosol.jl` (16/16).

## Population definition

`AerosolBins` discretises the population into size bins; each parcel carries
ONE bin (sampled ∝ bin mass fraction at release — `AerosolModel.bin`), so the
ensemble reproduces the size distribution in expectation.

```julia
# Lognormal MASS distribution (median diameter, geometric σ):
bins = AerosolBins(dg=2.0e-6, sigma_g=2.0, nbins=8, rho_p=2600.0, species=:dust)

# Explicit size bins (edges in m, optional per-bin mass fractions):
bins = AerosolBins([0.1e-6, 1e-6, 10e-6, 100e-6];
                   mass_frac=[0.2, 0.5, 0.3], species=:ash, rho_p=2400.0)
```

Species sets the 550 nm complex refractive index (`AEROSOL_REFRACTIVE`):
`:dust` 1.53−0.003i · `:sulfate` 1.43 · `:bc` 1.95−0.79i · `:oc` 1.53−0.006i ·
`:seasalt` 1.50 · `:ash` 1.52−0.0015i · `:water` 1.33 — or pass `m_ri=...`
and `lambda=...` explicitly.

## Processes (operator-split after advection)

| Process | Formulation | Validation |
|---|---|---|
| **Settling** | Stokes–Cunningham (Sutherland μ(T), slip Cc(Kn(T,p))) + iterative Clift–Gauvin drag in the transition regime | 1/10/100 µm at STP: 3.5e-5 / 3.06e-3 / 0.25 m/s vs textbook 3.5e-5 / 3.1e-3 / 0.25–0.31; column test: 934 m fall vs v_t·t = 938 m |
| **Dry deposition** | `vd = vs + 1/(ra + rb + ra·rb·vs)` for parcels below `href` (default 30 m); ra from MO similarity (ψ_h, uses the per-dataset ustar/L), rb = 1/(u*(Sc^-2/3 + 10^(-3/St))) (Slinn) | accumulation-mode minimum reproduced (vd(0.3 µm)=1.7e-4 ≪ vd(10 nm), vd(20 µm)) |
| **Wet deposition** | size-dependent, from ConvFields tpf/cc: in-cloud Λ = f_act(d)·3.36e-4·I^0.79 (activation fraction → 0.9 above ~0.1 µm); below-cloud Λ = 1.5·E(d)·I_SI/d_drop with Slinn-shaped E(d) (Greenfield gap ~1e-3 at 0.5–1 µm) | 3 h ERA5 storm run: 8.7 % scavenged (vs 21 % for the fully-soluble bulk tracer — the gap is the Greenfield minimum) |
| **AOD** | per-bin Mie Q_ext (Bohren–Huffman series; **matches the Wiscombe 1979 benchmark exactly**: Q_ext(x=10, m=1.5−i) = 2.4173); τ_cell = Σ m·[Q_ext(π/4)d²/(ρ_p(π/6)d³)]/A_cell | Q_ext → 2.05 at x=300, interference peak 3.98, Rayleigh +2.9e-4 |

Mass budget over a full run (transport + settling + both depositions):
active 0.9125 + dry 0.0004 + wet 0.0871 = **1.00000**.

## Usage

```julia
bins = AerosolBins(dg=2.0e-6, sigma_g=2.0, nbins=8, rho_p=2600.0, species=:dust)
aero = AerosolModel(length(state.status), ctx, bins; href=30.0)
convfields = [ConvFields(rt.metfields3D[n]) for n in 1:rt.nest]   # for wet dep

for date in timesteps
    read_met_step!(rt, date)
    step_model!(state, ctx, date; convfields=convfields, aerosol=aero)
end

aod = compute_aod(aero, state, ctx)        # Niout × Njout column AOD at λ
# deposited mass maps (met grid, per nest):
aero.drydep[1]; aero.wetdep[1]
```

## Hygroscopic growth (κ-Köhler)

Each species carries a hygroscopicity κ (Petters & Kreidenweis 2007;
`AEROSOL_KAPPA`: sulfate 0.61, sea salt 1.28, OC 0.10, dust 0.03, BC/ash 0;
override with `AerosolBins(...; kappa=...)`).  The growth factor on the dry
diameter,  GF(RH) = (1 + κ·RH/(1−RH))^(1/3)  (a_w ≈ RH, capped at 98 %),
enters every process:

- **settling & dry deposition**: wet diameter d·GF and volume-mixed wet
  density (water 1000 kg/m³);
- **wet scavenging**: collection efficiency E(d_wet), activation f_act(d_wet);
- **optics/AOD**: per-bin Mie Q_ext is precomputed on an RH grid
  (`AEROSOL_RH_GRID`, 0–0.98) with the wet diameter AND the volume-mixed
  refractive index (water 1.33); `compute_aod(aero, state, ctx, t_array)`
  evaluates the parcel RH (reconstructed from ρ/T/p via Tv — no extra met
  fields needed) and interpolates the table.  Particle NUMBER always comes
  from the DRY mass and dry diameter.

Validated: GF(sulfate, 90 %) = 1.866 (analytic), seasalt > sulfate > dust
ordering, σ_ext humidification factor for accumulation-mode sulfate within
the observed 2–6 range, wet settling enhancement.

Notes:
- bins are assigned automatically to newly released parcels at the first
  aerosol step; per-parcel `state.mass` stays bulk DRY aerosol mass [kg];
- don't pass `wetdepparams` (bulk scheme) together with `aerosol` — the
  size-dependent scavenging replaces it;
- deposition is skipped in backward runs; settling is applied in both
  directions (it is a deterministic drift, sign handled by |Δt|... forward
  descent — for backward source-receptor runs of settling species, revisit);
- Kelvin term neglected (relevant only below ~50 nm at high RH);
- AOD assumes externally mixed, spherical particles at one wavelength.

## Limitation: no Aitken → accumulation growth (audited 2026-07-03)

EASLEE's aerosol is a **static-size sectional tracer**: every parcel is
assigned ONE dry bin at release (∝ bin mass fractions) and the bin diameters
never change during a run.  The κ-Köhler `growth_factor` is EQUILIBRIUM WATER
UPTAKE — reversible with RH, used for settling/deposition/scavenging/optics —
it does not grow the dry size.  Consequently the model has **no inter-mode
growth**: an Aitken particle stays Aitken until it deposits.

Missing processes (in order of importance for Aitken→accumulation):
1. **Condensational growth** — H2SO4 / organic vapours.  Nothing condenses
   onto the AerosolModel bins: the bins' dry size is fixed at release.
2. **In-cloud (aqueous) processing** — SO2 oxidation in droplets adds sulfate
   mass to activated particles on each cloud cycle (the process that carves
   the Hoppel minimum); fastest mass-adder in the marine boundary layer.
3. **Brownian coagulation** — Aitken self- and inter-modal coagulation
   (number sink, slow mass mover).
4. **Nucleation** (new-particle formation feeding the Aitken mode).

Reference magnitudes (`test/test_aerosol_growth.jl`: 60-bin sectional box
model, Fuchs coagulation + Fuchs–Sutugin H2SO4 condensation at the canonical
marine daytime 1e7 cm⁻³, 72 h, Aitken 600 cm⁻³/40 nm + accumulation
200 cm⁻³/160 nm):
- daytime condensation GR ≈ 0.8 nm/h (pure H2SO4; observed total 1–5 nm/h
  with organics/bases), Aitken dg 40 → 64 nm in 72 h,
- N(>100 nm) +31 %, total N −20 % (coagulation sink),
- EASLEE for the same case: dg 40 → 40 nm by construction.

Implementation path (when growth matters): keep per-parcel bins but make the
bin index dynamic — per-parcel dry-mass growth from (a) condensation of a
per-parcel sulfate source (kinetic Fuchs–Sutugin flux, cheap: d and T/P known)
and (b) an in-cloud sulfate increment when cc > threshold; re-bin the parcel when
its dry volume crosses a bin edge.
Coagulation needs inter-parcel pairing (the mixing pairlist machinery could
serve) and matters mostly for number, not mass — defer.
