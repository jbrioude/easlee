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

## Meteorology

The scripts in [`../download/`](../download/README.md) fetch public
meteorological data (GFS, MERRA-2, prepared ERA5 archives) and write it here by
default — `data/GFS`, `data/MERRA2`, and so on. Point `[meteorology] path` at
whichever directory you used. These files are large, which is the other reason
this directory is gitignored.

See [../HOWTO.md](../HOWTO.md) for how the two fit together in a run.
