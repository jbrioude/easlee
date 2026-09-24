#!/usr/bin/env python3
"""Fetch a ready-made ERA5 archive for EASLEE's validation cases (ETEX-1, CAPTEX-83) and
write its AVAILABLE index.

    python3 era5_archive.py --dataset etex   --out ./ETEX_ERA5
    python3 era5_archive.py --dataset captex --out ./CAPTEX_ERA5
    python3 era5_archive.py --url https://.../my_era5.tar.gz --out ./ERA5

The archives are hourly ERA5 model-level GRIB files in EASLEE's IFS_EA_MET_/IFS_EA_CLD_
layout, as used by validation/validate_captex.jl and configs/output_etex.toml.  Their
download links are read from links.toml next to this script (fill them in once the files
are hosted).  General ERA5 retrievals need flex_extract (https://www.flexpart.eu/flex_extract)
with EASLEE's listvariables/ERA5*.vars naming; this script only unpacks a prepared archive.
"""
import argparse
import os
import sys
import tarfile
import urllib.request
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import easlee_available  # noqa: E402

LINKS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "links.toml")


def read_links():
    links = {}
    if os.path.isfile(LINKS):
        for line in open(LINKS):
            line = line.split("#")[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                links[k.strip()] = v.strip().strip('"')
    return links


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--dataset", choices=["etex", "captex"], help="named archive from links.toml")
    ap.add_argument("--url", help="archive URL (.tar.gz / .tgz / .zip), overrides --dataset")
    ap.add_argument("--out", default=None, help="target directory (default ./<DATASET>_ERA5)")
    ap.add_argument("--vars", default="ERA5_CONV", help="variables namelist for the index")
    ap.add_argument("--keep-archive", action="store_true")
    a = ap.parse_args(argv)
    url = a.url or read_links().get(f"{a.dataset}_era5", "") if (a.url or a.dataset) else ""
    if not url:
        sys.exit("no URL: give --url, or --dataset etex|captex with the link filled in download/links.toml")
    out = a.out or f"./{(a.dataset or 'archive').upper()}_ERA5"
    os.makedirs(out, exist_ok=True)
    arc = os.path.join(out, os.path.basename(url.split("?")[0]) or "era5_archive.tar.gz")
    if not os.path.exists(arc):
        print(f"get {url}")
        def hook(n, bs, total):
            if total > 0:
                print(f"\r  {n*bs/1e6:8.1f} / {total/1e6:.1f} MB", end="", flush=True)
        urllib.request.urlretrieve(url, arc, hook); print()
    else:
        print(f"exists: {arc}")
    print(f"unpacking into {out}")
    if zipfile.is_zipfile(arc):
        with zipfile.ZipFile(arc) as z:
            z.extractall(out)
    else:
        with tarfile.open(arc) as t:
            t.extractall(out)
    a.keep_archive or os.remove(arc)
    easlee_available.write_available(out, "ERA5", a.vars)
    print(f"[[meteorology]]  model = \"ERA5\"  path = \"{os.path.abspath(out)}\"  variables = \"{a.vars}\"  available = \"auto\"")


if __name__ == "__main__":
    main()
