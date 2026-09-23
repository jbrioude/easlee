#!/usr/bin/env python3
"""
verify_nc.py — verify and sanitize EASLEE NetCDF output (any layout: the
engine grid_<name>.nc,
or the old single-species files).

For each file it reports, per DATA variable, the count of finite / NaN / Inf
values and the min/max of the finite data, then (with --fix) rewrites the file
in place with every non-finite value replaced by NaN, so the file only ever
contains finite numbers or the NaN fill — never Inf.

Usage:
  python3 tools/verify_nc.py juliaout/<run>/grid_<name>.nc [more.nc] [--fix]
"""
import sys
import numpy as np
from netCDF4 import Dataset

paths = [a for a in sys.argv[1:] if not a.startswith("--")]
FIX   = "--fix" in sys.argv

COORDS = {"lon", "lat", "plev", "height", "time"}

for path in paths:
    print(f"\n=== {path} ===")
    bad = False
    d = Dataset(path, "a" if FIX else "r")
    for vn in d.variables:
        if vn in COORDS:
            continue
        v = d.variables[vn]
        a = np.ma.filled(v[:].astype(float), np.nan)
        ninf = int(np.isinf(a).sum())
        nnan = int(np.isnan(a).sum())
        nfin = int(np.isfinite(a).sum())
        fin  = a[np.isfinite(a)]
        rng  = (float(fin.min()), float(fin.max())) if fin.size else (None, None)
        flag = "  <-- HAS INF!" if ninf else ""
        print(f"  {vn:12s} finite={nfin:9d} NaN={nnan:9d} Inf={ninf:6d}  "
              f"range=[{rng[0]},{rng[1]}]{flag}")
        if ninf:
            bad = True
        if FIX and (ninf or not np.array_equal(np.isnan(a), v[:].mask if hasattr(v[:], 'mask') else np.isnan(a))):
            a[~np.isfinite(a)] = np.nan
            v[:] = a
    if FIX:
        d.close()
        print("  -> sanitized (non-finite → NaN), written in place" if bad else "  -> already clean")
    else:
        d.close()
        print("  (clean)" if not bad else "  (Inf present — rerun with --fix)")
