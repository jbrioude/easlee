#!/bin/sh
# Build the distribution tarball.  Run from the parent directory of this tree.
#   sh easlee_v1_1_shared/make_tarball.sh
# Manifest.toml is EXCLUDED on purpose: one resolved on Julia 1.12 fails to
# precompile on 1.10 (reported by Yue, 2026-09-04).  Julia recreates them as
# soon as anyone runs in the tree, so they are excluded here rather than
# deleted once.  results/ and juliaout/ are run output, not distribution.
set -e
D=$(basename "$(cd "$(dirname "$0")" && pwd)")
cd "$(cd "$(dirname "$0")/.." && pwd)"
OUT="${D}_$(date +%Y%m%d).tar.gz"
tar czf "$OUT" \
    --exclude="$D/easlee/Manifest.toml" \
    --exclude="$D/webgui/Manifest.toml" \
    --exclude="$D/results/*" \
    --exclude="$D/juliaout/*" \
    --exclude='*.rej' --exclude='*.orig' --exclude='*.mem' \
    "$D"
echo "wrote $OUT"
n=$(tar tzf "$OUT" | grep -c 'Manifest.toml' || true)
[ "$n" -eq 0 ] || { echo "ERROR: $n Manifest.toml in the tarball"; exit 1; }
echo "check: no Manifest.toml shipped"
