# download/ — get public meteorology for EASLEE

Three scripts (Python 3, standard library only) fetch public meteorological data into a
directory and write its `AVAILABLE` index (plus the `AVAILABLE.files` sidecar), so a run
pointed at that directory starts without scanning the files.  The web GUI's **Data** tab
runs the same scripts.

| script | data | account | typical call |
|---|---|---|---|
| `gfs.py` | NCEP GFS 0.25/0.5/1.0° forecasts (NOMADS, last ~10 days; NOAA archive on AWS, 2021 onward) | none | `python3 download/gfs.py --analysis 2026-06-03T00 --hours 48 --out ./data/GFS` |
| `merra2.py` | NASA MERRA-2 3-hourly model-level reanalysis (GES DISC) | free Earthdata account | `EARTHDATA_TOKEN=… python3 download/merra2.py --start 2022-01-18 --end 2022-01-29 --out ./data/MERRA2` |
| `era5_archive.py` | prepared ERA5 archives of the validation cases (ETEX-1, CAPTEX-83), links in `links.toml` | none | `python3 download/era5_archive.py --dataset etex --out ./data/ERA5_ETEX` |
| `easlee_available.py` | (re)write the index of an existing directory | — | `python3 download/easlee_available.py --model GFS --dir ./data/GFS` |

Then in the run configuration:

```toml
[[meteorology]]
model     = "GFS"          # or MERRA2, ERA5
path      = "./data/GFS"
variables = "GFS"          # listvariables/<name>.vars
available = "auto"         # uses <path>/AVAILABLE when it matches the directory; "scan" forces a rescan
```

Notes
- GFS steps must be ≥ 3 h: EASLEE reads the 3-h averaged surface fluxes, which the analysis
  file (f000) does not carry.  Files are ~500 MB (0.25°), ~130 MB (0.5°), ~45 MB (1.0°).
- MERRA-2 needs the Earthdata account authorised for "NASA GESDISC DATA ARCHIVE"
  (profile → Applications).  Either a bearer token (`--token`, `EARTHDATA_TOKEN`) or a
  `~/.netrc` line `machine urs.earthdata.nasa.gov login <user> password <pw>`.
  Default collections `tavg3_3d_asm_Nv` + `tavg3_3d_mst_Ne` are what `MERRA2.vars` reads.
- General ERA5 retrievals are not scripted here: EASLEE reads flex_extract-style GRIB
  (`IFS_EA_MET_*` + `IFS_EA_CLD_*`); see https://www.flexpart.eu/flex_extract.
- The index is also written by the model itself after its first scan, and refreshed
  automatically when files are added, removed or resized (the sidecar records them).
