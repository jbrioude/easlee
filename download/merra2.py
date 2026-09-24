#!/usr/bin/env python3
"""Download MERRA-2 3-hourly model-level files from NASA GES DISC for EASLEE and write AVAILABLE.

    python3 merra2.py --start 2022-01-18 --end 2022-01-29 --out ./MERRA2

Collections (default = what listvariables/MERRA2.vars reads): tavg3_3d_asm_Nv (M2T3NVASM,
winds, T, q, omega, pressure, heights, cloud) and its companion tavg3_3d_mst_Ne (M2T3NEMST,
convective mass flux and precipitation fluxes).  Add --collections to fetch others
(tavg3_3d_trb_Ne, tavg3_3d_cld_Nv, inst3_3d_asm_Nv, tavg1_2d_flx_Nx, tavg1_2d_slv_Nx).

Authentication (free NASA Earthdata account, https://urs.earthdata.nasa.gov):
  either  --token <Earthdata bearer token>   or  EARTHDATA_TOKEN in the environment,
  or      a ~/.netrc entry:  machine urs.earthdata.nasa.gov login <user> password <pw>
Files land in <out>/MERRA2_<stream>.<collection>.YYYYMMDD.nc4.  Point EASLEE at it with
    [[meteorology]]  model = "MERRA2"  path = "<out>"  variables = "MERRA2"  available = "auto"
"""
import argparse
import datetime as dt
import http.cookiejar
import netrc
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import easlee_available  # noqa: E402

BASE = "https://data.gesdisc.earthdata.nasa.gov/data/MERRA2"
PRODUCT = {
    "tavg3_3d_asm_Nv": "M2T3NVASM", "tavg3_3d_mst_Ne": "M2T3NEMST", "tavg3_3d_trb_Ne": "M2T3NETRB",
    "tavg3_3d_cld_Nv": "M2T3NVCLD", "inst3_3d_asm_Nv": "M2I3NVASM", "tavg1_2d_flx_Nx": "M2T1NXFLX",
    "tavg1_2d_slv_Nx": "M2T1NXSLV", "inst1_2d_asm_Nx": "M2I1NXASM",
}
DEFAULT_COLL = ["tavg3_3d_asm_Nv", "tavg3_3d_mst_Ne"]


def stream(year):
    return 100 if year <= 1991 else 200 if year <= 2000 else 300 if year <= 2010 else 400


def opener(token):
    cj = http.cookiejar.CookieJar()
    handlers = [urllib.request.HTTPCookieProcessor(cj)]
    if not token:
        try:
            auth = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        except (FileNotFoundError, netrc.NetrcParseError):
            auth = None
        if auth is None:
            sys.exit("no Earthdata credentials: give --token, set EARTHDATA_TOKEN, or add urs.earthdata.nasa.gov to ~/.netrc")
        pm = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        pm.add_password(None, "https://urs.earthdata.nasa.gov", auth[0], auth[2])
        handlers.append(urllib.request.HTTPBasicAuthHandler(pm))
    op = urllib.request.build_opener(*handlers)
    if token:
        op.addheaders = [("Authorization", f"Bearer {token}")]
    return op


def fetch(op, url, dest, tries=3):
    tmp = dest + ".part"
    for k in range(tries):
        try:
            with op.open(url, timeout=120) as r, open(tmp, "wb") as f:
                if "text/html" in (r.headers.get("Content-Type") or ""):
                    raise RuntimeError("got an HTML page instead of data: check the Earthdata login / GES DISC authorisation")
                total = int(r.headers.get("Content-Length", "0") or 0); done = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk); done += len(chunk)
                    if total:
                        print(f"\r  {done/1e6:8.1f} / {total/1e6:.1f} MB", end="", flush=True)
            print(); os.replace(tmp, dest); return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            if e.code in (401, 403):
                sys.exit(f"HTTP {e.code} from GES DISC: the Earthdata account must be authorised for 'NASA GESDISC DATA ARCHIVE' (Earthdata profile > Applications)")
            print(f"  HTTP {e.code}, retry {k+1}/{tries}")
        except Exception as e:
            print(f"  {type(e).__name__}: {e}, retry {k+1}/{tries}")
        time.sleep(5 * (k + 1))
    return False


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--start", required=True, help="YYYY-MM-DD")
    ap.add_argument("--end", help="YYYY-MM-DD (default = start)")
    ap.add_argument("--collections", default=",".join(DEFAULT_COLL))
    ap.add_argument("--out", default="./MERRA2")
    ap.add_argument("--token", default=os.environ.get("EARTHDATA_TOKEN", ""))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-available", action="store_true")
    a = ap.parse_args(argv)
    d0 = dt.datetime.strptime(a.start, "%Y-%m-%d").date()
    d1 = dt.datetime.strptime(a.end, "%Y-%m-%d").date() if a.end else d0
    colls = [c.strip() for c in a.collections.split(",") if c.strip()]
    for c in colls:
        if c not in PRODUCT:
            sys.exit(f"unknown collection {c}; known: {', '.join(PRODUCT)}")
    op = None if a.dry_run else opener(a.token)
    a.dry_run or os.makedirs(a.out, exist_ok=True)
    n_ok = n_skip = n_miss = 0
    d = d0
    while d <= d1:
        for c in colls:
            got = False
            for st in (stream(d.year), 401):          # 401 = reprocessed months (2020-2021)
                fn = f"MERRA2_{st}.{c}.{d:%Y%m%d}.nc4"
                url = f"{BASE}/{PRODUCT[c]}.5.12.4/{d:%Y}/{d:%m}/{fn}"
                dest = os.path.join(a.out, fn)
                if a.dry_run:
                    print(url); got = True; break
                if os.path.exists(dest) and os.path.getsize(dest) > 0:
                    print(f"exists: {dest}"); n_skip += 1; got = True; break
                print(f"get {url}")
                if fetch(op, url, dest):
                    n_ok += 1; got = True; break
            if not got and not a.dry_run:
                print(f"  missing: {c} {d}"); n_miss += 1
        d += dt.timedelta(days=1)
    print(f"done: {n_ok} downloaded, {n_skip} already present, {n_miss} missing")
    if not a.dry_run and not a.no_available and (n_ok + n_skip) > 0:
        easlee_available.write_available(a.out, "MERRA2")
    return 0 if n_miss == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
