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
import shutil
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import easlee_available  # noqa: E402

LINKS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "links.toml")

# The prepared archives are flex_extract IFS_EA_* files with IFS_SFC_* companions,
# which is what ERA5_CAPTEX.vars declares (ERA5_CONV.vars expects _MET_/_CLD_).
DEFAULT_VARS = {"etex": "ERA5_CAPTEX", "captex": "ERA5_CAPTEX"}


def read_links():
    links = {}
    if os.path.isfile(LINKS):
        for line in open(LINKS):
            line = line.split("#")[0].strip()
            if "=" in line:
                k, v = line.split("=", 1)
                links[k.strip()] = v.strip().strip('"')
    return links


def strip_top_level(out, keep):
    """If the archive wrapped everything in one directory, lift its contents up.

    `tar -cf a.tar DIR` stores DIR/... while `tar -cf a.tar -C DIR .` does not,
    and we want the meteorology directly under `out` either way.  `keep` is the
    downloaded archive, which must not be mistaken for the wrapper.
    """
    entries = [e for e in os.listdir(out) if os.path.join(out, e) != keep]
    if len(entries) != 1:
        return None
    inner = os.path.join(out, entries[0])
    if not os.path.isdir(inner):
        return None
    for name in os.listdir(inner):
        shutil.move(os.path.join(inner, name), os.path.join(out, name))
    os.rmdir(inner)
    return entries[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--dataset", choices=["etex", "captex"], help="named archive from links.toml")
    ap.add_argument("--url", help="archive URL (.tar.gz / .tgz / .zip), overrides --dataset")
    ap.add_argument("--out", default=None, help="target directory (default ./<DATASET>_ERA5)")
    ap.add_argument("--vars", default=None,
                    help="variables namelist for the index (default: ERA5_CAPTEX for the prepared archives)")
    ap.add_argument("--keep-archive", action="store_true")
    a = ap.parse_args(argv)
    url = a.url or read_links().get(f"{a.dataset}_era5", "") if (a.url or a.dataset) else ""
    if not url:
        sys.exit("no URL: give --url, or --dataset etex|captex with the link filled in download/links.toml")
    vars_ = a.vars or DEFAULT_VARS.get(a.dataset, "ERA5_CONV")
    out = a.out or f"./data/ERA5_{(a.dataset or 'archive').upper()}"
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
    wrapper = strip_top_level(out, arc)
    if wrapper:
        print(f"lifted the contents of {wrapper}/ up into {out}")
    a.keep_archive or os.remove(arc)
    easlee_available.write_available(out, "ERA5", vars_)
    print(f"[[meteorology]]  model = \"ERA5\"  path = \"{out}\"  variables = \"{vars_}\"  available = \"auto\"")


if __name__ == "__main__":
    main()
