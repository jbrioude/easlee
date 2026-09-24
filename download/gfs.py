#!/usr/bin/env python3
"""Download NCEP GFS 0.25/0.5/1.0-degree GRIB2 forecasts for EASLEE and write the AVAILABLE index.

Two ways to say what you want:

  one analysis + its forecast (what "test it now" needs):
      python3 gfs.py --analysis 2026-06-03T00 --hours 48 --out ./GFS
  a date range of analyses, cycles and steps (a reanalysis-like archive):
      python3 gfs.py --start 2026-06-03 --end 2026-06-08 --cycles 00,12 --steps 3,6,9,12,15 --out ./GFS

Files land in  <out>/YYYYMMDD/HH/gfs.tHHz.pgrb2.<res>.fSSS  (the layout listvariables/GFS.vars
expects; steps must be >= 3 h because the surface fluxes EASLEE reads are 3-h averages).
Sources: NOMADS (last ~10 days) then the public NOAA archive on AWS (2021 onward); no account.
Point EASLEE at it with
    [[meteorology]]  model = "GFS"  path = "<out>"  variables = "GFS"  available = "auto"
Options: --res 0p25|0p50|1p00 (default 0p25), --source auto|nomads|aws, --dry-run,
         --no-available.  Existing complete files are skipped.
"""
import argparse
import datetime as dt
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import easlee_available  # noqa: E402

NOMADS = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
AWS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"


def urls(date, cycle, step, res, source):
    fn = f"gfs.t{cycle}z.pgrb2.{res}.f{step:03d}"
    day = date.strftime("%Y%m%d")
    cands = []
    if source in ("auto", "nomads"):
        cands.append(f"{NOMADS}/gfs.{day}/{cycle}/atmos/{fn}")
    if source in ("auto", "aws"):
        cands.append(f"{AWS}/gfs.{day}/{cycle}/atmos/{fn}")
    return fn, cands


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", "0") or 0)


def fetch(url, dest, tries=3):
    tmp = dest + ".part"
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
                total = int(r.headers.get("Content-Length", "0") or 0)
                done = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    done += len(chunk)
                    if total:
                        print(f"\r  {done/1e6:8.1f} / {total/1e6:.1f} MB", end="", flush=True)
            print()
            os.replace(tmp, dest)
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            print(f"  HTTP {e.code} ({e.reason}), retry {k+1}/{tries}")
        except Exception as e:  # network hiccup
            print(f"  {type(e).__name__}: {e}, retry {k+1}/{tries}")
        time.sleep(5 * (k + 1))
    return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--analysis", help="analysis time YYYY-MM-DDTHH (one cycle + forecast)")
    ap.add_argument("--hours", type=int, default=48, help="forecast length with --analysis (default 48)")
    ap.add_argument("--start", help="first analysis day YYYY-MM-DD")
    ap.add_argument("--end", help="last analysis day YYYY-MM-DD (default = start)")
    ap.add_argument("--cycles", default="00,12", help="cycles with --start (default 00,12)")
    ap.add_argument("--steps", default="3,6,9,12,15", help="forecast steps [h] with --start (default 3,6,9,12,15)")
    ap.add_argument("--res", default="0p25", choices=["0p25", "0p50", "1p00"])
    ap.add_argument("--out", default="./GFS")
    ap.add_argument("--source", default="auto", choices=["auto", "nomads", "aws"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-available", action="store_true", help="do not write the AVAILABLE index")
    a = ap.parse_args(argv)

    jobs = []   # (date, cycle, step)
    if a.analysis:
        t0 = dt.datetime.strptime(a.analysis, "%Y-%m-%dT%H")
        cyc = f"{(t0.hour // 6) * 6:02d}"
        for s in range(3, a.hours + 1, 3):
            jobs.append((t0.date(), cyc, s))
    elif a.start:
        d0 = dt.datetime.strptime(a.start, "%Y-%m-%d").date()
        d1 = dt.datetime.strptime(a.end, "%Y-%m-%d").date() if a.end else d0
        steps = sorted({int(s) for s in a.steps.split(",") if s.strip()})
        cycles = [c.strip().zfill(2) for c in a.cycles.split(",") if c.strip()]
        d = d0
        while d <= d1:
            for c in cycles:
                for s in steps:
                    jobs.append((d, c, s))
            d += dt.timedelta(days=1)
    else:
        ap.error("give --analysis (one forecast) or --start [--end] (date range)")
    if any(s < 3 for _, _, s in jobs):
        print("note: steps below 3 h have no 3-h averaged surface fluxes; EASLEE needs steps >= 3", file=sys.stderr)

    n_ok = n_skip = n_miss = 0
    for date, cyc, step in jobs:
        fn, cands = urls(date, cyc, step, a.res, a.source)
        dest = os.path.join(a.out, date.strftime("%Y%m%d"), cyc, fn)
        if a.dry_run:
            print(cands[0]); continue
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        got = False
        for url in cands:
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                try:
                    if remote_size(url) in (0, os.path.getsize(dest)):
                        print(f"exists: {dest}"); n_skip += 1; got = True; break
                except Exception:
                    print(f"exists: {dest}"); n_skip += 1; got = True; break
            print(f"get {url}")
            if fetch(url, dest):
                n_ok += 1; got = True; break
        if not got:
            print(f"  missing: {fn} ({date} {cyc}z)")
            n_miss += 1
    print(f"done: {n_ok} downloaded, {n_skip} already present, {n_miss} missing")
    if not a.dry_run and not a.no_available and (n_ok + n_skip) > 0:
        easlee_available.write_available(a.out, "GFS")
    return 0 if n_miss == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
