# Mineral-dust emission — Kok et al. (2014) scheme (`[[release]] kind = "dust"`)

Added 25 Sept 2026 (Jerome: "implement the model from Jasper Kok"). Code:
`easlee/src/physics/dust_kok.jl`; test `test/test_dust_kok.jl`; example
`configs/example_dust_kok.toml`.

## Physics

**Vertical dust flux** (Kok et al. 2014, ACP 14, 13023, Eq. 18), per land cell, for
u* > u*_t and zero otherwise:

    F_d = C_d · f_bare · f_clay · ρ_a · (u*² − u*_t²) / u*_st · (u*/u*_t)^(C_α (u*_st − u*_st0)/u*_st0)
    C_d = C_d0 · exp(−C_e (u*_st − u*_st0)/u*_st0),   u*_st = u*_t · sqrt(ρ_a/ρ_a0)

Constants (their Table 1): C_d0 = 4.4×10⁻⁵, C_e = 2.0, C_α = 2.7, u*_st0 = 0.16 m s⁻¹,
ρ_a0 = 1.225 kg m⁻³. The soil's erodibility enters through the standardised threshold
u*_st (a lower threshold means a more erodible soil, a larger C_d and a steeper wind
dependence), the clay fraction f_clay and the bare-soil fraction f_bare.

**Fluid threshold**: Shao & Lu (2000) for the most erodible grain size, D = 127 µm:
u*_t0 = sqrt(A_N (ρ_p g D/ρ_a + γ/(ρ_a D))), A_N = 0.0123, γ = 1.65×10⁻⁴ kg s⁻²,
ρ_p = 2650 kg m⁻³ → 0.21 m s⁻¹ at ρ_a0. **Soil moisture**: Fécan et al. (1999),
u*_t = u*_t0 · sqrt(1 + 1.21 (w − w′)^0.68) for gravimetric moisture w [%] above the
residual w′ = 0.0014 c² + 0.17 c (c = clay %).

**Emitted size distribution** (Kok 2011, PNAS 108, 1016, Eq. 6), volume per ln D over
0.2–20 µm: dV/dlnD = (D/c_V)[1 + erf(ln(D/D_s)/(√2 ln σ_s))] exp(−(D/λ)³), D_s = 3.4 µm,
σ_s = 3.0, λ = 12 µm. Independent of wind and soil (brittle fragmentation). Selected for
the species with `size.distribution = "kok"` (bins log-spaced over 0.2–20 µm, or the
given `size.edges`); the bin mass fractions are the integrals of dV/dlnD. With 8 bins:
7 % of the mass below 2.5 µm, the volume peak near 9 µm.

## Engine

Like the sea-spray source: the release box reserves a species/relbox and schedules no
parcels; `apply_dust!` runs every transport step over the land cells (lsm ≥ 0.5) of the
box, accumulates F_d · area · Δt per cell and spawns one parcel of `mass_per_parcel` at
the surface of the cell whenever the accumulator holds one quantum. `particles` is the
parcel-pool cap. Inputs per cell from the 2-D met dict: `ustar` (required; ERA5 `zust`,
GFS `ustar`), `sp`/`lnsp` and `2t` for the air density (else 1.2 kg m⁻³), `swvl1`
(volumetric soil water, converted with a bulk density of 1.3 g cm⁻³) when the namelist
maps it, else `[release.dust].soil_moisture`. Clay and bare fractions are constants or
NetCDF lon/lat maps (`clay_file`, `bare_file`, variables `lon|longitude`, `lat|latitude`,
`clay`/`bare`). Accumulators are checkpointed (`chk[:dust]`); the run-end
"dust (Kok 2014) budget" line reports the spawned parcels, emitted mass, land cells, the
largest cell flux and the unspawned remainder.

## Configuration

```toml
[[species]]
id = "dust"; kind = "aerosol"; aerosol_type = "dust"; density = 2650.0
settling = true; dry_deposition = true; wet_deposition = true
[species.size]
distribution = "kok"; nbins = 8

[[release]]
species = "dust"; kind = "dust"
lon = [-15.0, 35.0]; lat = [15.0, 35.0]; z = [0.0, 50.0]; zkind = "agl"
time = [2025-01-01T12:00:00, 2025-01-02T00:00:00]
mass_per_parcel = "2000kg"; particles = 400000
[release.dust]
clay = 0.15; bare = 1.0; soil_moisture = 0.5      # or clay_file / bare_file maps
# threshold_d = "127um"; cd0 = 4.4e-5; ce = 2.0; calpha = 2.7   (defaults)
```

## Not included (yet)

Drag partition for roughness elements (folded into `bare`), the Kok et al. (2014) Part 2
global tuning of C_d0 against observations, a snow/frozen-soil mask (add via `bare_file`),
and the soil texture database: the `clay` map is the user's (e.g. from the Zobler classes
of `/data/20CRv3/timeInvariantSI/soiltype.nc`, or FAO/HWSD). The constants carry the
uncertainties of the paper (C_d0 ±0.5×10⁻⁵, C_e ±0.3, C_α ±1.0); the scheme is not tuned here.

## First run (25 Sept 2026, `/tmp/smoke_dust_kok.toml` = the example on ERA5 of 20 Sept 2025)

Sahara box 15–35° N, 15° W–35° E, 12 h from 12 UTC: 14 186 land cells, 0.40 Mt emitted
before the 200 000-parcel pool filled (2 t quantum, filled by ~15 UTC), largest cell flux
3.8×10⁻⁶ kg m⁻² s⁻¹ — the order of the fluxes of Kok et al. (2014, Fig. 4) at u* ≈ 0.4–0.5
m s⁻¹. The example config now uses a 20 t quantum. Not validated against observations
(AERONET/MODIS dust AOD would be the next step, with a soil map instead of clay = 0.15).
