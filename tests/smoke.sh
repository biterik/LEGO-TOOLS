#!/usr/bin/env bash
#
# LEGO-TOOLS smoke test
# =====================
# Exercises afc + a representative subset of lego-tools binaries on the
# bundled Ni test structure, and checks that:
#   1. afc roundtrips xyz -> lammps -> cel -> xyz preserving atom count
#   2. lego-analyze --json reports the expected atom count
#   3. lego-shift, lego-pbc-wrap, lego-remove-per-atom all run without
#      error and preserve the atom count
#
# Exit code: 0 on success, non-zero on the first failure.
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
# Funding: NFDI-MatWerk

set -euo pipefail

# -- Locate repo root regardless of where the script is invoked from ----------
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
AFC="$ROOT/afc/afc"
LEGO="$ROOT/lego-tools"

TEST_DATA="$ROOT/afc/test-data/relax_d0_Ni_Mishin04_x001_y110_z-110_16x4x16.xyz.gz"
EXPECTED_ATOMS=92160

TMP="$(mktemp -d -t lego-smoke.XXXXXX)"
trap 'rm -rf "$TMP"' EXIT

PASS=0
FAIL=0
fail() {
    echo "  FAIL: $*" >&2
    FAIL=$((FAIL + 1))
}
ok() {
    echo "  ok  : $*"
    PASS=$((PASS + 1))
}

need() {
    if [ ! -x "$1" ]; then
        echo "ERROR: required binary not found or not executable: $1" >&2
        echo "       run 'make' at the repo root first" >&2
        exit 2
    fi
}

# -- Preconditions -----------------------------------------------------------
need "$AFC"
need "$LEGO/lego-analyze"
need "$LEGO/lego-shift"
need "$LEGO/lego-pbc-wrap"
need "$LEGO/lego-remove-per-atom"

if [ ! -f "$TEST_DATA" ]; then
    echo "ERROR: test data not found: $TEST_DATA" >&2
    exit 2
fi

echo "LEGO-TOOLS smoke test"
echo "  repo:       $ROOT"
echo "  tmp:        $TMP"
echo "  test data:  $(basename "$TEST_DATA") ($EXPECTED_ATOMS atoms expected)"
echo

cp "$TEST_DATA" "$TMP/input.xyz.gz"

# -- 1) afc roundtrip: xyz -> lammps -> cel -> xyz ---------------------------
echo "[1] afc roundtrip xyz -> lammps -> cel -> xyz"

"$AFC" "$TMP/input.xyz.gz" "$TMP/a.lmp"       >/dev/null 2>&1 \
    && ok "xyz -> lammps" || fail "afc xyz -> lammps"

"$AFC" -t Ni "$TMP/a.lmp" "$TMP/a.cel"        >/dev/null 2>&1 \
    && ok "lammps -> cel" || fail "afc lammps -> cel"

"$AFC" "$TMP/a.cel" "$TMP/a.xyz"              >/dev/null 2>&1 \
    && ok "cel -> xyz"   || fail "afc cel -> xyz"

# Count atoms in the final xyz: first line is N.
rt_n=$(head -n 1 "$TMP/a.xyz" | tr -d '[:space:]')
if [ "$rt_n" = "$EXPECTED_ATOMS" ]; then
    ok "roundtrip preserves atom count ($rt_n)"
else
    fail "roundtrip atom count: got $rt_n, expected $EXPECTED_ATOMS"
fi

# -- 2) lego-analyze --json reports expected count ---------------------------
echo "[2] lego-analyze --json on the lammps-data file"

if "$LEGO/lego-analyze" --json "$TMP/a.lmp" > "$TMP/analyze.json" 2>&1; then
    ok "lego-analyze --json ran"
else
    fail "lego-analyze --json exited non-zero"
    cat "$TMP/analyze.json" >&2 || true
fi

# Extract "natoms": NNN, robust to whitespace, no jq dependency.
an_n=$(tr -d '[:space:]' < "$TMP/analyze.json" \
         | grep -oE '"natoms":[0-9]+' \
         | head -n1 \
         | grep -oE '[0-9]+' || true)
if [ "$an_n" = "$EXPECTED_ATOMS" ]; then
    ok "lego-analyze natoms = $an_n"
else
    fail "lego-analyze natoms: got '${an_n:-<missing>}', expected $EXPECTED_ATOMS"
fi

# -- 3) lego-shift preserves atom count --------------------------------------
echo "[3] lego-shift 1.5 -0.25 3.0"

"$LEGO/lego-shift" 1.5 -0.25 3.0 "$TMP/a.lmp" "$TMP/shift.lmp" >/dev/null 2>&1 \
    && ok "lego-shift ran" || fail "lego-shift exited non-zero"

"$LEGO/lego-analyze" --json "$TMP/shift.lmp" > "$TMP/shift.json" 2>&1 || true
sh_n=$(tr -d '[:space:]' < "$TMP/shift.json" \
         | grep -oE '"natoms":[0-9]+' \
         | head -n1 \
         | grep -oE '[0-9]+' || true)
if [ "$sh_n" = "$EXPECTED_ATOMS" ]; then
    ok "lego-shift preserves atom count ($sh_n)"
else
    fail "lego-shift natoms: got '${sh_n:-<missing>}', expected $EXPECTED_ATOMS"
fi

# -- 4) lego-pbc-wrap preserves atom count -----------------------------------
echo "[4] lego-pbc-wrap"

"$LEGO/lego-pbc-wrap" "$TMP/shift.lmp" "$TMP/wrap.lmp" >/dev/null 2>&1 \
    && ok "lego-pbc-wrap ran" || fail "lego-pbc-wrap exited non-zero"

"$LEGO/lego-analyze" --json "$TMP/wrap.lmp" > "$TMP/wrap.json" 2>&1 || true
wr_n=$(tr -d '[:space:]' < "$TMP/wrap.json" \
         | grep -oE '"natoms":[0-9]+' \
         | head -n1 \
         | grep -oE '[0-9]+' || true)
if [ "$wr_n" = "$EXPECTED_ATOMS" ]; then
    ok "lego-pbc-wrap preserves atom count ($wr_n)"
else
    fail "lego-pbc-wrap natoms: got '${wr_n:-<missing>}', expected $EXPECTED_ATOMS"
fi

# -- 5) lego-remove-per-atom strips to id/type/x/y/z -------------------------
echo "[5] lego-remove-per-atom"

"$LEGO/lego-remove-per-atom" "$TMP/a.lmp" "$TMP/strip.lmp" >/dev/null 2>&1 \
    && ok "lego-remove-per-atom ran" || fail "lego-remove-per-atom exited non-zero"

"$LEGO/lego-analyze" --json "$TMP/strip.lmp" > "$TMP/strip.json" 2>&1 || true
st_n=$(tr -d '[:space:]' < "$TMP/strip.json" \
         | grep -oE '"natoms":[0-9]+' \
         | head -n1 \
         | grep -oE '[0-9]+' || true)
if [ "$st_n" = "$EXPECTED_ATOMS" ]; then
    ok "lego-remove-per-atom preserves atom count ($st_n)"
else
    fail "lego-remove-per-atom natoms: got '${st_n:-<missing>}', expected $EXPECTED_ATOMS"
fi

# -- 6) lego-dislo inserts a dislocation -------------------------------------
echo "[6] lego-dislo (W 1/2<111> edge in a small crystal)"

DISLO="$ROOT/dislo/lego-dislo"
need "$DISLO"

python3 "$ROOT/dislo/examples/make_crystal.py" \
    --a0=3.1652 --x=1,1,1 --y=-1,0,1 --z=1,-2,1 --size=30,30,8 --center \
    --out="$TMP/dislo_in.data" >/dev/null 2>&1 \
    && ok "make_crystal.py generated a sample" \
    || fail "make_crystal.py exited non-zero"

cat > "$TMP/dislo_card.yaml" <<'EOF'
material:
  name: W
  lattice_parameter: 3.1652
  cij:
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
box_orientation:
  x: [1, 1, 1]
  y: [-1, 0, 1]
  z: [1, -2, 1]
dislocations:
  - point: [0.1, 0.05, 0.0]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
EOF

di_in=$("$LEGO/lego-analyze" --json "$TMP/dislo_in.data" 2>/dev/null \
         | tr -d '[:space:]' | grep -oE '"natoms":[0-9]+' | head -n1 \
         | grep -oE '[0-9]+' || true)

if "$DISLO" "$TMP/dislo_card.yaml" "$TMP/dislo_in.data" "$TMP/dislo_out.data" \
        > "$TMP/dislo_report.txt" 2>/dev/null; then
    ok "lego-dislo ran"
else
    fail "lego-dislo exited non-zero"
fi

grep -q "K_GPa" "$TMP/dislo_report.txt" \
    && ok "report contains the K factor" \
    || fail "report is missing the K factor"

di_n=$("$LEGO/lego-analyze" --json "$TMP/dislo_out.data" 2>/dev/null \
         | tr -d '[:space:]' | grep -oE '"natoms":[0-9]+' | head -n1 \
         | grep -oE '[0-9]+' || true)
if [ -n "$di_n" ] && [ "$di_n" = "$di_in" ]; then
    ok "lego-dislo preserves atom count ($di_n)"
else
    fail "lego-dislo natoms: got '${di_n:-<missing>}', expected '${di_in:-?}'"
fi

echo
echo "----------------------------------------------------------------------"
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ] || exit 1
