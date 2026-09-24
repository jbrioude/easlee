# download/ — get public meteorology for EASLEE

Three scripts (Python 3, standard library only) fetch public meteorological data into a
directory and write its `AVAILABLE` index (plus the `AVAILABLE.files` sidecar), so a run
pointed at that directory starts without scanning the files.  The web GUI's **Data** tab
runs the same scripts.

| script | data | account | typical call |
|---|---|---|---|
| `gfs.py` | NCEP GFS 0.25/0.5/1.0° forecasts (NOMADS, last ~10 days; NOAA archive on AWS, 2021 onward) | none | `python3 download/gfs.py --analysis 2026-06-03T00 --hours 48 --out ./data/GFS` |
| `merra2.py` | NASA MERRA-2 3-hourly model-level reanalysis (GES DISC) | free Earthdata account | `EARTHDATA_TOKEN=… python3 download/merra2.py --start 2022-01-18 --end 2022-01-29 --out ./data/MERRA2` |
| `era5_archive.py` | prepared ERA5 archives of the validation cases (ETEX-1, CAPTEX-83), links in `links.toml` | none | `python3 download/era5_archive.py --dataset etex` |
| `easlee_available.py` | (re)write the index of an existing directory | — | `python3 download/easlee_available.py --model GFS --dir ./data/GFS` |
| `repack_grib.py` | shrink a GRIB archive ~3x with 12-bit JPEG packing, for publishing or archiving | none | `python3 download/repack_grib.py /data/ETEX --out ./ETEX_jpeg` |

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

## Publishing an archive

`repack_grib.py` re-packs a GRIB archive from `grid_simple`/16-bit to
`grid_jpeg`/12-bit, which is about **3x smaller on disk** (measured 3.33x on the
ETEX ERA5 archive: 1.49 GB -> 0.45 GB per day).  It mirrors the directory tree,
copies non-GRIB files unchanged, never touches the source, and checks that every
output file holds exactly as many GRIB messages as its input.

The packing is lossy: each field is quantised to 4096 levels across its own
range, so the largest error is range/4096, about 0.004 m/s on a wind field and
0.004 K on a temperature field.  EASLEE's ETEX score is unchanged by it
(`docs/BENCHMARKS.md` 0).  Constant and near-constant fields, which the JPEG
coder rejects, keep their original packing; that is why the gain is ~3x and not
16/12.

`era5_archive.py --dataset etex` unpacks into `data/ERA5_ETEX` and writes the
`AVAILABLE` index with the `ERA5_CAPTEX` namelist, which is what the `configs/etex_*`
and `configs/validate_*` files already expect, so a fresh clone needs no path editing.
The CAPTEX archive lands in `data/ERA5_CAPTEX` the same way.

Pack the result as a plain `.tar`, not `.tar.gz`.  JPEG-packed GRIB is already
compressed, so gzip buys almost nothing and costs time at both ends;
`era5_archive.py` unpacks `.tar`, `.tar.gz` and `.zip` alike.
