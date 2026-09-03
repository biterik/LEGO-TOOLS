#!/usr/bin/env bash
#
# lego-dislo benchmark
# ====================
# Generates a large synthetic bcc sample (default ~1e8 atoms; override with
# BENCH_ATOMS, e.g. BENCH_ATOMS=1e7 for a quicker run), then times the
# field evaluation at 1, 4, and all available threads.  The "field time"
# printed by lego-dislo excludes file I/O; total wall time is reported too.
#
# The sample file is a few GB at full size — it is written to $BENCH_TMP
# (default: a fresh directory under $TMPDIR) and removed afterwards.
#
# Usage: bench.sh [BENCH_ATOMS=1e8] [BENCH_TMP=dir]
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
# Funding: NFDI-MatWerk

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
DISLO="$HERE/../lego-dislo"
ATOMS="${BENCH_ATOMS:-1e8}"

if [ ! -x "$DISLO" ]; then
    echo "ERROR: $DISLO not built — run 'make' at the repo root" >&2
    exit 2
fi

TMP="${BENCH_TMP:-$(mktemp -d -t lego-dislo-bench.XXXXXX)}"
CLEAN=""
[ -z "${BENCH_TMP:-}" ] && CLEAN="$TMP"
trap '[ -n "$CLEAN" ] && rm -rf "$CLEAN"' EXIT

# nx*ny*nz*2 atoms; pick a roughly cubic cell grid for the requested count.
NCELL=$(python3 -c "import math; n=float('$ATOMS'); print(int(round((n/2)**(1/3))))")
NATOMS=$((2 * NCELL * NCELL * NCELL))

MAXT=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

echo "lego-dislo benchmark"
echo "  cells:   ${NCELL}^3  ->  $NATOMS atoms"
echo "  tmp:     $TMP"
echo "  threads: 1, 4, $MAXT"
echo

cc -O2 -o "$TMP/gen_bcc" "$HERE/gen_bcc.c"
"$TMP/gen_bcc" "$NCELL" "$NCELL" "$NCELL" 3.1652 "$TMP/bench.data"

cat > "$TMP/card.yaml" <<'EOF'
material:
  name: W
  lattice_parameter: 3.1652
  cij:
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
box_orientation:
  x: [1, 0, 0]
  y: [0, 1, 0]
  z: [0, 0, 1]
dislocations:
  - point: [0.1, 0.05, 0.0]
    line: [0, 0, 1]
    burgers_direction: [1, 0, 0]
    burgers_magnitude: 1.0
    glide_normal: [0, 1, 0]
EOF

for t in 1 4 "$MAXT"; do
    echo "--- threads = $t ---"
    /usr/bin/time "$DISLO" --threads="$t" \
        "$TMP/card.yaml" "$TMP/bench.data" "$TMP/bench_out.data" \
        > /dev/null
    rm -f "$TMP/bench_out.data"
done
