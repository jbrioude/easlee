#!/usr/bin/env python3
"""repack_grib.py — shrink a GRIB meteorological archive with 12-bit JPEG packing.

ERA5 archives as retrieved are packed with `grid_simple` at 16 bits per value.
Repacking them to `grid_jpeg` at 12 bits makes them roughly 3x smaller on disk,
which is what turns a 33 GB archive into something one can publish.  The cost is
a quantisation of each field to 4096 levels across its own range, i.e. a maximum
absolute error of range/4096: about 0.004 m/s on a wind field and 0.004 K on a
temperature field.  EASLEE scores are unchanged by it (docs/BENCHMARKS.md 0).

    python3 download/repack_grib.py /data/ETEX --out ./ETEX_jpeg
    python3 download/repack_grib.py /data/ETEX --out ./ETEX_jpeg --jobs 8
    python3 download/repack_grib.py /data/ETEX --dry-run

The directory tree is mirrored: every GRIB file is repacked, everything else
(AVAILABLE, station lists, notes) is copied unchanged.  The source is never
modified.

Messages that cannot be JPEG-packed -- constant and near-constant fields, which
the JPEG coder rejects -- keep their original packing.  That is expected and is
why the saving is ~3x rather than the 16/12 ratio.  Every output file is checked
to hold exactly as many GRIB messages as its input; if one does not, the original
is copied instead and the file is reported.

Requires the ecCodes command-line tools (grib_set, grib_count), which ship with
the eccodes package that EASLEE already depends on for reading GRIB.
"""

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys

TOOLS = ("grib_set", "grib_count")


def have_tools():
    return [t for t in TOOLS if shutil.which(t) is None]


def is_grib(path):
    """True if the file starts with the GRIB magic, whatever its name."""
    try:
        with open(path, "rb") as f:
            return f.read(4) == b"GRIB"
    except OSError:
        return False


def count_messages(path):
    try:
        out = subprocess.run(["grib_count", path], capture_output=True, text=True, check=True)
        return int(out.stdout.strip())
    except (subprocess.CalledProcessError, ValueError):
        return -1


def repack_one(src, dst, bits):
    """Repack one GRIB file.  Returns (src, ok, note, size_in, size_out)."""
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    size_in = os.path.getsize(src)
    cmd = ["grib_set", "-f", "-r", "-s", f"packingType=grid_jpeg,bitsPerValue={bits}", src, dst]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
        return (src, False, f"grib_set produced nothing ({proc.returncode})", size_in, size_in)

    n_in, n_out = count_messages(src), count_messages(dst)
    if n_in < 0 or n_out < 0 or n_in != n_out:
        shutil.copy2(src, dst)
        return (src, False, f"message count {n_in} -> {n_out}, original kept", size_in, size_in)

    size_out = os.path.getsize(dst)
    if size_out >= size_in:          # already packed, or packing made it worse
        shutil.copy2(src, dst)
        return (src, True, "no gain, original kept", size_in, size_in)

    skipped = proc.stderr.count("ECCODES ERROR")
    note = f"{skipped} message(s) left as-is" if skipped else ""
    return (src, True, note, size_in, size_out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("src", help="meteorological archive directory to repack")
    ap.add_argument("--out", help="target directory (required unless --dry-run)")
    ap.add_argument("--bits", type=int, default=12, help="bits per value (default 12)")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 4) // 2),
                    help="files to repack in parallel (default: half the cores)")
    ap.add_argument("--dry-run", action="store_true", help="list what would be done, change nothing")
    a = ap.parse_args(argv)

    missing = have_tools()
    if missing and not a.dry_run:
        sys.exit(f"missing ecCodes tool(s): {', '.join(missing)}. "
                 "Install eccodes (Debian/Ubuntu: apt install libeccodes-tools).")
    if not os.path.isdir(a.src):
        sys.exit(f"not a directory: {a.src}")
    if not a.out and not a.dry_run:
        sys.exit("give --out <dir> (the source is never modified), or --dry-run")

    gribs, others = [], []
    for root, _, files in os.walk(a.src):
        for f in sorted(files):
            p = os.path.join(root, f)
            (gribs if is_grib(p) else others).append(p)

    if not gribs:
        sys.exit(f"no GRIB files found under {a.src}")

    total_in = sum(os.path.getsize(p) for p in gribs)
    print(f"{len(gribs)} GRIB files ({total_in/1e9:.2f} GB), {len(others)} other files")
    if a.dry_run:
        print(f"would repack to grid_jpeg / {a.bits} bits, expect roughly "
              f"{total_in/3e9:.2f} GB (~3x smaller)")
        return 0

    rel = lambda p: os.path.relpath(p, a.src)
    done_in = done_out = 0
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(repack_one, p, os.path.join(a.out, rel(p)), a.bits): p for p in gribs}
        for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            src, ok, note, s_in, s_out = fut.result()
            done_in += s_in
            done_out += s_out
            if not ok:
                failures.append((rel(src), note))
            print(f"\r  {i}/{len(gribs)}  {done_in/1e9:.2f} -> {done_out/1e9:.2f} GB", end="", flush=True)
    print()

    for p in others:
        dst = os.path.join(a.out, rel(p))
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        shutil.copy2(p, dst)

    ratio = done_in / done_out if done_out else 1.0
    print(f"\n{a.src} -> {a.out}")
    print(f"  {done_in/1e9:.2f} GB -> {done_out/1e9:.2f} GB  ({ratio:.2f}x smaller)")
    print(f"  {len(others)} non-GRIB file(s) copied unchanged")
    if failures:
        print(f"  {len(failures)} file(s) kept at their original packing:")
        for name, note in failures[:10]:
            print(f"    {name}: {note}")
        if len(failures) > 10:
            print(f"    ... and {len(failures)-10} more")
    print("\nCheck it before publishing:")
    print(f"  julia --project=easlee scripts/easlee.jl <a config pointed at {a.out}>")
    return 0


if __name__ == "__main__":
    sys.exit(main())
