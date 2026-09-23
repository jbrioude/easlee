# EASLEE species database

Generated from `easlee/data/species.toml`.  FLEXPART removal data parsed
verbatim; further species carry a molar weight only.  **?** = wet/dry deposition
unknown.  Lifetime τ = 1/(k_OH·[OH]+J) at [OH]=1.6e6 molec/cm³, T=273 K, given
only where OH kinetics are known — it is the timescale the `[[species]].decay`
knob can represent, not a modelled reaction.

| species | M | phase | decay | dry dep | wet dep | τ_chem |
|---|--:|---|---|---|---|--:|
| ACO2 | 58.0 | gas | stable | ? | ? | — |
| AEROSOL | ? | aero | stable | ρ1400,0.25µm | crain=1 | — |
| AIRTRACER | 29.0 | gas | stable | ? | ? | — |
| ALD2 | 24.0 | gas | stable | ? | ? | — |
| Ba-140 | 140.0 | aero | 12.8d | ρ1900,0.65µm | crain=3.6 | — |
| BC | ? | aero | stable | ρ2000,0.02µm | crain=1 | — |
| Br | 79.9 | gas | stable | ? | ? | — |
| BrCl | 115.4 | gas | stable | ? | ? | — |
| BrO | 95.9 | gas | stable | ? | ? | — |
| BrONO2 | 141.9 | gas | stable | ? | ? | — |
| C10H16 | 136.2 | gas | stable | ? | ? | 14h |
| C2H4 | 28.1 | gas | stable | ? | ? | — |
| C2H5OH | 46.1 | gas | stable | ? | ? | 58h |
| C2H6 | 30.0 | gas | stable | D=1.3 H=0.0019 f=0.0 | A=8.4e-05,B=0.62 | 40d |
| C2O3 | 75.0 | gas | stable | ? | ? | — |
| C3H6 | 42.1 | gas | stable | ? | ? | — |
| C3H8 | 44.1 | gas | stable | ? | ? | 8d |
| CH-3-I-131 | 131.0 | gas | 8.0d | ? | A=8.0e-07,B=0.62 | — |
| CH2O | 30.0 | gas | stable | ? | ? | 32h |
| CH3COCH3 | 58.1 | gas | stable | ? | ? | 4d |
| CH3COCHO | 72.1 | gas | stable | ? | ? | — |
| CH3COOH | 60.0 | gas | stable | ? | ? | — |
| CH3O2 | 47.0 | gas | stable | ? | ? | — |
| CH3OH | 32.0 | gas | stable | ? | ? | 9d |
| CH3OOH | 48.0 | gas | stable | ? | ? | — |
| CH4 | 16.0 | gas | stable | D=0.95 H=0.0014 f=0.0 | A=9.3e-05,B=0.62 | 5.4yr |
| Cl | 35.5 | gas | stable | ? | ? | — |
| Cl2 | 70.9 | gas | stable | ? | ? | — |
| Cl2O2 | 102.9 | gas | stable | ? | ? | — |
| ClNO2 | 81.5 | gas | stable | ? | ? | — |
| ClO | 51.5 | gas | stable | ? | ? | — |
| ClONO2 | 97.5 | gas | stable | ? | ? | — |
| CO | 28.0 | gas | stable | D=1.3 H=0.00097 f=0.0 | A=8.6e-05,B=0.62 | 50d |
| Cs-137 | 137.0 | aero | 30yr | ρ1900,0.65µm | crain=3.6 | — |
| DMS | 62.1 | gas | stable | ? | ? | — |
| DUST-COARSE | ? | aero | stable | ρ2500,12µm | crain=1 | — |
| DUST-FINE | ? | aero | stable | ρ2500,0.22µm | crain=1 | — |
| G-HCH | 290.8 | gas | stable | D=4.04 H=714.0 f=0.1 | A=3.9e-05,B=0.62 | — |
| H2 | 2.0 | gas | stable | ? | ? | — |
| H2O2 | 34.0 | gas | stable | D=1.4 H=100000.0 f=1.0 | A=7.4e-05,B=0.62 | 4d |
| HBr | 80.9 | gas | stable | D=1.9 H=1.0e9 f=0.0 | A=5.8e-05,B=0.62 | — |
| HCHO | 30.0 | gas | stable | D=1.3 H=6000.0 f=0.0 | A=8.4e-05,B=0.62 | 32h |
| HCl | 36.5 | gas | stable | ? | ? | — |
| HCOOH | 46.0 | gas | stable | ? | ? | — |
| HNO2 | 47.0 | gas | stable | D=1.6 H=48.0 f=0.1 | A=7.0e-05,B=0.62 | — |
| HNO3 | 63.0 | gas | stable | D=1.9 H=210000.0 f=0.0 | A=5.8e-05,B=0.62 | — |
| HNO4 | 79.0 | gas | stable | ? | ? | — |
| HO2 | 33.0 | gas | stable | ? | ? | — |
| HO2NO2 | 79.0 | gas | stable | ? | ? | — |
| HOBr | 96.9 | gas | stable | ? | ? | — |
| HOCl | 52.5 | gas | stable | ? | ? | — |
| HYPROPO2 | 91.0 | gas | stable | ? | ? | — |
| I-131 | 131.0 | aero | 8.0d | ρ1900,0.5µm | crain=3.6 | — |
| I2-131 | 131.0 | gas | 8.0d | D=2.7 H=100000.0 f=0.1 | A=8.0e-05,B=0.62 | — |
| IC3H7O2 | 75.0 | gas | stable | ? | ? | — |
| ISOP | 68.1 | gas | stable | ? | ? | 6h |
| ISPD | 70.0 | gas | stable | ? | ? | — |
| Kr-85 | 85.0 | gas | 11yr | ? | ? | — |
| MP | ? | aero | stable | ρ1234,1µm | A=-9.0e-10,B=-9.9 | — |
| MSA | 96.1 | gas | stable | ? | ? | — |
| N2O | 44.0 | gas | stable | D=1.6 H=0.024 f=0.0 | A=7.2e-05,B=0.62 | — |
| N2O5 | 108.0 | gas | stable | ? | ? | — |
| NH2 | 16.0 | gas | stable | ? | ? | — |
| NH3 | 17.0 | gas | stable | D=1.0 H=59.0 f=0.0 | A=9.8e-05,B=0.62 | — |
| NO | 30.0 | gas | stable | D=1.3 H=0.002 f=0.0 | A=8.4e-05,B=0.62 | — |
| NO2 | 46.0 | gas | stable | D=1.6 H=0.01 f=0.1 | A=7.1e-05,B=0.62 | — |
| NO3 | 62.0 | gas | stable | ? | ? | — |
| NO3-aero | ? | aero | stable | ρ2000,0.2µm | crain=1 | — |
| O | 16.0 | gas | stable | ? | ? | — |
| O3 | 48.0 | gas | stable | D=1.6 H=0.01 f=1.0 | A=7.4e-05,B=0.62 | 133d |
| OH | 17.0 | gas | stable | ? | ? | — |
| OLE | 24.0 | gas | stable | ? | ? | — |
| ORGNTR | 77.0 | gas | stable | ? | ? | — |
| PAN | 121.0 | gas | stable | D=2.6 H=3.6 f=0.1 | A=5.0e-05,B=0.62 | — |
| PAR | 12.0 | gas | stable | ? | ? | 9d |
| Particulate-Radionuclide | ? | aero | stable | ρ1900,0.65µm | crain=3.6 | — |
| PCB28 | 257.5 | gas | stable | D=3.81 H=3.02 f=0.1 | A=4.0e-05,B=0.62 | — |
| ROOH | 47.0 | gas | stable | ? | ? | — |
| ROR | 28.0 | gas | stable | ? | ? | — |
| Ru-106 | 106.0 | aero | 1yr | ρ1900,0.65µm | crain=3.6 | — |
| RXPAR | 12.0 | gas | stable | ? | ? | — |
| SO2 | 64.0 | gas | stable | D=1.9 H=1.3 f=0.0 | A=6.3e-05,B=0.62 | — |
| SO4 | 96.1 | gas | stable | ? | ? | — |
| SO4-aero | ? | aero | stable | ρ2000,0.2µm | crain=1 | — |
| Sr-90 | 90.0 | aero | 29yr | ρ1900,0.65µm | crain=3.6 | — |
| Xe-131m | 131.0 | gas | 11.8d | ? | ? | — |
| Xe-133 | 133.0 | gas | 5.2d | ? | ? | — |
| Xe-133m | 133.0 | gas | 2.2d | ? | ? | — |
| Xe-135 | 135.0 | gas | 0.4d | ? | ? | — |
| XO2 | 44.0 | gas | stable | ? | ? | — |
| XO2N | 44.0 | gas | stable | ? | ? | — |
| Y-91 | 91.0 | aero | 0yr | ρ1900,0.65µm | crain=3.6 | — |
