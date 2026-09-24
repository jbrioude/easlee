# Observation data for the validation runs

This directory is where EASLEE looks for the measurements used by the CAPTEX
and ETEX validation runs. It ships empty: the tracer-experiment datasets are
distributed by their own archives, not by us.

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

Note that the observations are only half of what a validation run needs: the
meteorological archive that drives the run is separate again, and is set by
`[meteorology] path` in each config. See [../HOWTO.md](../HOWTO.md).
