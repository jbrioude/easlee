#!/usr/bin/env python3
"""Write the EASLEE AVAILABLE index of a meteorology directory from the file names.

EASLEE's `[meteorology].available = "auto"` scans the directory and opens every file to
read its valid time; that costs minutes per hundred GRIB files.  The download scripts of
this directory call this module after each download so the index is ready before the first
run, and the model reuses it as long as it matches what is on disk.

    python3 easlee_available.py --model GFS    --dir ./GFS
    python3 easlee_available.py --model MERRA2 --dir ./MERRA2
    python3 easlee_available.py --model ERA5   --dir ./ETEX_ERA5 [--vars ERA5_CONV]

Format (FLEXPART lineage, read by initialize.jl `readavailable!`):
    line 1   absolute path of the variables namelist (listvariables/<MODEL>.vars)
    line 2-3 header
    then     YYYYMMDD HHMMSS      <path relative to the met dir>  ON DISC
One line per VALID TIME; a file holding several times (MERRA-2) appears once per time.
When two files give the same valid time the shorter forecast lead wins (as the scanner).
"""
import argparse
import datetime as dt
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
INDEX_NAME = "AVAILABLE"


def vars_path(name):
    """listvariables/<name>.vars of this EASLEE tree (or an explicit path)."""
    if os.path.isfile(name):
        return os.path.abspath(name)
    cand = os.path.join(REPO, "listvariables", name if name.endswith(".vars") else name + ".vars")
    if not os.path.isfile(cand):
        sys.exit(f"variables namelist not found: {name} (looked for {cand})")
    return cand


# ── per-model (valid time, lead hours, relative path) ─────────────────────────────
GFS_PAT = re.compile(r"^gfs\.t(\d{2})z\.pgrb2\.(0p25|0p50|1p00)\.f(\d{3})$")


def gfs_entries(metdir):
    """<dir>/YYYYMMDD/HH/gfs.tHHz.pgrb2.0p25.fSSS  ->  valid = date + HH + SSS h."""
    out = []
    for root, _, files in os.walk(metdir):
        for f in files:
            m = GFS_PAT.match(f)
            if not m:
                continue
            rel = os.path.relpath(os.path.join(root, f), metdir)
            parts = rel.split(os.sep)
            day = next((p for p in parts if re.fullmatch(r"\d{8}", p)), None)
            if day is None:
                continue
            base = dt.datetime.strptime(day, "%Y%m%d") + dt.timedelta(hours=int(m.group(1)))
            step = int(m.group(3))
            out.append((base + dt.timedelta(hours=step), step, rel))
    return out


MERRA_PAT = re.compile(r"^MERRA2_\d{3}\.(inst\d|tavg\d)_(\d)d_([a-z]+)_N([a-z])\.(\d{8})\.nc4$")


def merra2_entries(metdir, main_token="_asm_Nv."):
    """One entry per valid time of each MAIN file (the companion collections are read via
    the main file, listvariables/MERRA2.vars `companions`).  Times from the file when
    netCDF4 is installed, else from the collection convention (tavgN: N/2 h offset)."""
    out = []
    try:
        import netCDF4  # optional
    except ImportError:
        netCDF4 = None
    for root, _, files in os.walk(metdir):
        for f in sorted(files):
            m = MERRA_PAT.match(f)
            if not m or main_token not in f:
                continue
            path = os.path.join(root, f)
            rel = os.path.relpath(path, metdir)
            day = dt.datetime.strptime(m.group(5), "%Y%m%d")
            times = None
            if netCDF4 is not None:
                try:
                    with netCDF4.Dataset(path) as ds:
                        tv = ds.variables["time"]
                        mu = re.search(r"since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})", tv.units)
                        base = dt.datetime.strptime(mu.group(1) + "T" + mu.group(2), "%Y-%m-%dT%H:%M:%S")
                        times = [base + dt.timedelta(minutes=float(t)) for t in tv[:]]
                except Exception:
                    times = None
            if times is None:
                kind, n = m.group(1), int(m.group(1)[-1])
                off = dt.timedelta(hours=n / 2) if kind.startswith("tavg") else dt.timedelta(0)
                times = [day + off + dt.timedelta(hours=n * k) for k in range(24 // n)]
            out.extend((t, 0, rel) for t in times)
    return out


ERA5_PAT = re.compile(r"^IFS_[A-Z0-9]+(?:_[A-Z]+)?_(\d{8})_(\d{3,4})_(\d{1,2})\.grb$")
ERA5_EXCLUDE = (re.compile(r"IFS_PT_MET"), re.compile(r"FULL_PT"), re.compile(r"^IFS_SFC_"))


def main_tokens(varsfile):
    """`companions = A>B, C>D` of the namelist -> the MAIN-file tokens [A, C] (as the
    model's scanner does); empty when the namelist declares no companions."""
    toks = []
    try:
        for line in open(varsfile):
            m = re.match(r"^\s*companions\s*=\s*(.+?)\s*(?:#.*)?$", line, re.I)
            if m:
                for pat in m.group(1).split(","):
                    parts = pat.split(">")
                    if len(parts) == 2:
                        toks.append(parts[0].strip())
    except OSError:
        pass
    return toks


def era5_entries(metdir, tokens=()):
    """IFS_EA_MET_YYYYMMDD_HHMM_S.grb or IFS_EA_YYYYMMDD_HHMM_S.grb -> base + S h (base hour
    without leading zero).  Only MAIN files are listed (companions _CLD_/_SFC_ come with them)."""
    out = []
    for root, _, files in os.walk(metdir):
        for f in files:
            if any(x.search(f) for x in ERA5_EXCLUDE):
                continue
            if tokens and not any(t in f for t in tokens):
                continue
            m = ERA5_PAT.match(f)
            if not m:
                continue
            hhmm = m.group(2)
            base = dt.datetime.strptime(m.group(1), "%Y%m%d") + dt.timedelta(
                hours=int(hhmm[:-2] or 0), minutes=int(hhmm[-2:]))
            step = int(m.group(3))
            out.append((base + dt.timedelta(hours=step), step, os.path.relpath(os.path.join(root, f), metdir)))
    return out


ENTRIES = {"GFS": gfs_entries, "MERRA2": merra2_entries, "ERA5": era5_entries}
DEFAULT_VARS = {"GFS": "GFS", "MERRA2": "MERRA2", "ERA5": "ERA5_CONV"}


def write_available(metdir, model, varsfile=None, out=None, quiet=False):
    """Build and write <metdir>/AVAILABLE; returns (path, n_times, first, last)."""
    model = model.upper()
    if model not in ENTRIES:
        sys.exit(f"unknown model {model}; known: {', '.join(ENTRIES)}")
    vp = vars_path(varsfile or DEFAULT_VARS[model])
    toks = main_tokens(vp)
    if model == "ERA5":
        entries = era5_entries(metdir, toks)
    elif model == "MERRA2":
        entries = merra2_entries(metdir, toks[0] if toks else "_asm_Nv.")
    else:
        entries = ENTRIES[model](metdir)
    if not entries:
        sys.exit(f"no {model} files found under {metdir}")
    best = {}
    for valid, step, rel in entries:
        if valid not in best or step < best[valid][0]:
            best[valid] = (step, rel)
    out = out or os.path.join(metdir, INDEX_NAME)
    with open(out, "w") as f:
        f.write(vp + "\n")
        f.write("YYYYMMDD HHMISS\n")
        f.write("________ ______      __________      __________\n")
        for valid in sorted(best):
            f.write(f"{valid:%Y%m%d} {valid:%H%M%S}      {best[valid][1]}  ON DISC\n")
    # sidecar <index>.files: the MAIN files (relative path, size) — the model's validity check
    mains = sorted({rel for _, _, rel in entries})
    with open(out + ".files", "w") as f:
        for rel in mains:
            f.write(f"{rel}\t{os.path.getsize(os.path.join(metdir, rel))}\n")
    ts = sorted(best)
    if not quiet:
        print(f"{out}: {len(ts)} valid times, {ts[0]} -> {ts[-1]}  (namelist {os.path.basename(vp)})")
        gaps = [(a, b) for a, b in zip(ts, ts[1:]) if (b - a) > (ts[1] - ts[0])]
        if gaps:
            print(f"  note: {len(gaps)} gap(s) longer than the first interval, e.g. {gaps[0][0]} -> {gaps[0][1]}")
    return out, len(ts), ts[0], ts[-1]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", required=True, choices=sorted(ENTRIES))
    ap.add_argument("--dir", required=True, help="meteorology directory")
    ap.add_argument("--vars", default=None, help="variables namelist name or path (default per model)")
    ap.add_argument("--out", default=None, help="index path (default <dir>/AVAILABLE)")
    a = ap.parse_args(argv)
    write_available(a.dir, a.model, a.vars, a.out)


if __name__ == "__main__":
    main()
