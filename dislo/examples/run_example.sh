#!/usr/bin/env bash
#
# run_example.sh — one-command runner for the lego-dislo examples.
#
#   ./run_example.sh edge     W 1/2<111>{110} edge dislocation
#   ./run_example.sh screw    W 1/2<111> screw dislocation
#   ./run_example.sh dipole   two opposite edge dislocations
#
# Generates the small W crystal with make_crystal.py, runs lego-dislo with
# the corresponding card, and compares the result against the committed
# expected output (tolerance 1e-6 A, to allow for libm differences across
# platforms).  Outputs land in ./out/.
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
# Funding: NFDI-MatWerk

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
DISLO="$HERE/../lego-dislo"
CASE="${1:-edge}"

case "$CASE" in
    edge|dipole) X=1,1,1;  Y=-1,0,1; Z=1,-2,1 ;;
    screw)       X=-1,-1,2; Y=1,-1,0; Z=1,1,1 ;;
    *) echo "usage: $0 edge|screw|dipole" >&2; exit 1 ;;
esac

if [ ! -x "$DISLO" ]; then
    echo "ERROR: $DISLO not built — run 'make' at the repo root" >&2
    exit 2
fi

mkdir -p "$HERE/out"
IN="$HERE/out/${CASE}_in.data"
OUT="$HERE/out/${CASE}_out.data"

python3 "$HERE/make_crystal.py" --a0=3.1652 --x="$X" --y="$Y" --z="$Z" \
        --size=60,60,8.2 --center --out="$IN"

"$DISLO" --report="$HERE/out/${CASE}_report.txt" \
         "$HERE/$CASE.yaml" "$IN" "$OUT"

echo
echo "wrote $OUT (and ${CASE}_report.txt)"

EXPECTED="$HERE/expected/${CASE}_out.data"
if [ -f "$EXPECTED" ]; then
    python3 - "$OUT" "$EXPECTED" <<'EOF'
import sys

def read(path):
    pos, in_atoms = {}, False
    for line in open(path):
        s = line.split("#")[0].strip()
        if not s:
            continue
        if s.startswith("Atoms"):
            in_atoms = True
            continue
        if in_atoms:
            t = s.split()
            if not t[0].lstrip("-").isdigit():
                in_atoms = False
                continue
            pos[int(t[0])] = tuple(float(v) for v in t[2:5])
    return pos

a, b = read(sys.argv[1]), read(sys.argv[2])
assert set(a) == set(b), "atom id mismatch vs expected output"
mx = max(max(abs(pa[k] - b[i][k]) for k in range(3)) for i, pa in a.items())
print(f"comparison vs expected: {len(a)} atoms, max deviation {mx:.2e} A")
assert mx < 1e-6, "DEVIATES from the expected output!"
print("OK — matches the expected output")
EOF
else
    echo "(no expected/ reference for '$CASE' — skipping comparison)"
fi
