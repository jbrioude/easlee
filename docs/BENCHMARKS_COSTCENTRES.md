# Cost-centre scaling — plume vs domainfill, 100K–20M parcels (zion, 2026-07-18)

Machine zion (Threadripper 3970X, **32 threads**), ERA5 May-2025 **crop="off"
(met global 60S–60N non croppé)**, fenêtre 3 h / pas 15 min (12–13 pas),
turbulence Hanna, mixing pairwise, dépôt gaz HNO3 (dry+wet), kernel adaptatif.
Postes = `bench_report` (`[run] benchmark=true`). Temps = secondes TOTALES sur
le run. metread ~constant (indépendant de N) ; il domine à 100K, marginal à 20M.

NOTE : le dépôt HUMIDE construit les convfields → la redistribution convective
mass-flux tourne (poste `convection`) même avec `scheme="off"` — couplage
intrinsèque wetdep↔convection dans EASLEE.

## PLUME (release ponctuel, parcelles CONCENTRÉES)

| poste \ N            | 100K  | 1M    | 10M    | 20M    | scaling 1M→20M (×20) |
|----------------------|------:|------:|-------:|-------:|---------------------:|
| advection+turbulence | 3.16  | 17.05 | 151.2  | 298.7  | 17.5× (~linéaire)    |
| mixing (pairwise)    | 1.06  | 12.01 | 120.4  | 243.3  | 20.3× (~linéaire*)   |
| kernel adaptatif out | 1.07  | 2.75  | 96.9   | 272.9  | 99×  (SUPER-linéaire)|
| déposition dry+wet   | 0.80  | 1.21  | 5.92   | 10.66  | 8.8× (sous-linéaire) |
| (convection couplée) | 1.37  | 5.38  | 40.0   | 76.4   | 14×                  |
| metread (global)     | 129   | 128   | 128    | 130    | constant             |
| wall total           | 199 s | 237 s | 692 s  | 1272 s |                      |

## DOMAINFILL (boîte -10..70E × -35..10N × 0-20km, parcelles DISPERSÉES)

| poste \ N            | 100K  | 1M    | 10M    | 20M    | scaling 1M→20M (×20) |
|----------------------|------:|------:|-------:|-------:|---------------------:|
| advection+turbulence | 6.03  | 48.15 | 434.0  | 920.0  | 19.1× (~linéaire)    |
| mixing (pairwise)    | 0.004 | 0.094 | 1.00   | 1.75   | ~0 (parcelles éparses)|
| kernel adaptatif out | 1.17  | 2.83  | 37.1   | 79.3   | 28×  (super-linéaire) |
| déposition dry+wet   | 0.51  | 2.14  | 18.03  | 36.44  | 17×  (linéaire)      |
| (convection couplée) | 7.08  | 60.2  | 573.4  | 1154   | 19×  (DOMINE)        |
| metread (global)     | 74.8  | 74.4  | 75.7   | 75.1   | constant             |
| wall total           | 209 s | 310 s | 1305 s | 2486 s |                      |

*plume mixing ~linéaire car max_neighbors plafonné à 30 → O(N·30).

## Lectures

1. **advection+turbulence : linéaire en N** (les deux types). Domainfill ~3×
   plus cher/parcelle que plume — les parcelles dispersées sur le globe cassent
   la localité de cache du stencil météo.
2. **mixing : piloté par la DENSITÉ, pas par N.** Plume (dense) = poste majeur
   (243 s à 20M) ; domainfill (épars) = quasi gratuit (1,7 s). Peu de voisins
   co-localisés → peu de paires.
3. **kernel adaptatif (kNN) : SUPER-linéaire**, surtout en plume dense (99×
   pour ×20 parcelles). La recherche kNN coûte d'autant plus que les
   trajectoires se serrent. C'est le poste qui « rattrape » l'advection au 20M.
4. **déposition : toujours < 2 %**, ~linéaire (domainfill) / sous-linéaire
   (plume — seules les parcelles proches surface déposent).
5. **convection couplée au wetdep** domine le domainfill (48 % à 20M) : un
   remplissage global touche des colonnes convectives partout.
