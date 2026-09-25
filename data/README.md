# Data

This directory holds the data a run needs but the repository does not ship:
meteorology you download, and the measurements the validation runs score
against. Everything inside it is gitignored except this file, so nothing here
can be committed by accident.

## Validation observations

The measurements used by the CAPTEX and ETEX validation runs go here. They ship
separately, distributed by their own archives rather than by us.

```
data/ETEX/      emit-t1.txt, meas-t1.txt, stations.txt
data/CAPTEX/    emit-t1.txt, meas-t1.txt, stations.txt
```

`stations.txt` lists the receptor sites as `id lat lon`; `emit-t1.txt` and
`meas-t1.txt` are the release and measurement records in DATEM format.

If you already keep these datasets elsewhere on your machine, don't copy them
— point at them instead. Either link them in:

```bash
ln -s /where/you/keep/ETEX   data/ETEX
ln -s /where/you/keep/CAPTEX data/CAPTEX
```

or set an environment variable, which the `validation/check_*.jl` scripts read:

```bash
export EASLEE_ETEX_DIR=/where/you/keep/ETEX
export EASLEE_CAPTEX_DIR=/where/you/keep/CAPTEX
```

The TOML configs in `configs/` refer to `data/ETEX` and `data/CAPTEX` as paths
relative to the top of the repository, so start runs from there — or edit
`[validation] obs_dir` and `[output.receptors] file` to give an absolute path
of your own.

## Data sources

**ETEX-1 and CAPTEX-83 observations.**  Both come from the NOAA Air Resources
Laboratory's *Data Archive of Tracer Experiments and Meteorology* (DATEM),
<https://www.arl.noaa.gov/hysplit/datem-tracer-evaluation/datem/>, which is where
the `emit-t1.txt` / `meas-t1.txt` DATEM format used here originates.  The files
download directly, with no registration.  ARL asks that work using them carry the
acknowledgement:

> The authors gratefully acknowledge the NOAA Air Resources Laboratory (ARL) for
> the provision of the Data Archive of Tracer Experiments and Meteorology (DATEM).

The ETEX-1 measurements were made by the EC Joint Research Centre, whose original
archive is at
<https://remon.jrc.ec.europa.eu/past_activities/etex/site/datasets.htm>.  That
archive asks that the version number of the tracer data be quoted; the release-1
PMCH set is **etex1_v1.1.960505**.  Note that the readme distributed with it in
1996 describes the data as restricted to the participants of the ATMES exercise.
The files are today served openly from that public directory, but if you intend to
redistribute them rather than simply use them, ask the JRC first.

**Meteorology.**  The prepared archives are ECMWF ERA5, obtained through the
Copernicus Climate Change Service and redistributed under the Copernicus licence,
which asks for attribution.  Neither archive contains observation data: they hold
ERA5 fields, the `AVAILABLE` index and the MARS retrieval requests, nothing else.

## Meteorology

The scripts in [`../download/`](../download/README.md) fetch public
meteorological data (GFS, MERRA-2, prepared ERA5 archives) and write it here by
default — `data/GFS`, `data/MERRA2`, and so on. Point `[meteorology] path` at
whichever directory you used. These files are large, which is the other reason
this directory is gitignored.

See [../HOWTO.md](../HOWTO.md) for how the two fit together in a run.
