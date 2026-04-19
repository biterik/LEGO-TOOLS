#!/usr/bin/env bash
#
# LEGO-TOOLS large-file regression test
# =====================================
# The read/write paths in atomio.c call zlib's gzread / gzwrite, whose
# `unsigned int` length parameter cannot transfer more than ~2 GiB in one
# call.  atomio.c clamps every such call to GZ_CHUNK_MAX bytes (default
# 1 GiB) and loops.  This regression test exercises that loop by rebuilding
# atomio with a *very* small chunk (4 KiB), which forces hundreds of gzread
# / gzwrite calls per file, and then asserts that the final output is
# byte-identical to a reference run with the default chunk size.
#
# Also exercises LEGO_DEBUG=1 to prove the diagnostic output is wired up.
#
# The test restores the default-chunk build on exit so subsequent `make`
# calls see a clean state.

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
cd "$ROOT"

TMP="$(mktemp -d -t lego-large.XXXXXX)"
trap 'rm -rf "$TMP"; echo "[cleanup] rebuilding lego-tools with default chunk size"; make -C lego-tools clean >/dev/null 2>&1 || true; make -C lego-tools >/dev/null 2>&1' EXIT

TEST_DATA="$ROOT/afc/test-data/relax_d0_Ni_Mishin04_x001_y110_z-110_16x4x16.xyz.gz"
EXPECTED_ATOMS=92160

PASS=0
FAIL=0
ok()   { echo "  ok  : $*"; PASS=$((PASS + 1)); }
fail() { echo "  FAIL: $*" >&2; FAIL=$((FAIL + 1)); }

# ---------------------------------------------------------------------------
# 1) Produce a reference output with the default (1 GiB) chunk size.
# ---------------------------------------------------------------------------
echo "[1] reference run with default GZ_CHUNK_MAX"
make -C lego-tools clean >/dev/null 2>&1 || true
make -C lego-tools       >/dev/null 2>&1

./afc/afc "$TEST_DATA" "$TMP/ref.lmp" >/dev/null 2>&1 \
    && ok "afc produced reference lammps file" \
    || fail "afc failed on reference build"

./lego-tools/lego-shift 1.0 2.0 3.0 "$TMP/ref.lmp" "$TMP/ref.shift.lmp" \
    >/dev/null 2>&1 && ok "lego-shift produced reference output" \
                    || fail "lego-shift failed on reference build"

# ---------------------------------------------------------------------------
# 2) Rebuild with tiny chunk and reproduce the same output.
# ---------------------------------------------------------------------------
echo "[2] tiny-chunk run with GZ_CHUNK_MAX=4096"
make -C lego-tools clean >/dev/null 2>&1 || true
make -C lego-tools \
    CFLAGS="-O3 -Wall -Wextra -std=c11 -g -fopenmp -DGZ_CHUNK_MAX=4096" \
    >/dev/null 2>&1

./afc/afc "$TEST_DATA" "$TMP/tiny.lmp" >/dev/null 2>&1 \
    && ok "afc ran (tiny chunk in *readers*)" \
    || fail "afc failed with tiny chunks"

./lego-tools/lego-shift 1.0 2.0 3.0 "$TMP/tiny.lmp" "$TMP/tiny.shift.lmp" \
    >/dev/null 2>&1 && ok "lego-shift ran (tiny chunks)" \
                    || fail "lego-shift failed with tiny chunks"

# Byte-identical outputs prove the chunking loop reassembles correctly.
if cmp -s "$TMP/ref.shift.lmp" "$TMP/tiny.shift.lmp"; then
    ok "output with 4 KiB chunks byte-matches output with 1 GiB chunks"
else
    fail "output differs between default and tiny chunk builds"
    # Help debugging by showing the first differing byte range.
    cmp "$TMP/ref.shift.lmp" "$TMP/tiny.shift.lmp" >&2 || true
fi

# ---------------------------------------------------------------------------
# 3) LEGO_DEBUG=1 emits the expected diagnostic lines.
# ---------------------------------------------------------------------------
echo "[3] LEGO_DEBUG=1 diagnostic output"
LEGO_DEBUG=1 ./lego-tools/lego-analyze "$TMP/tiny.lmp" > /dev/null 2>"$TMP/dbg.err"
for needle in "file size:" "finished reading" "detected format:" "parsed "; do
    if grep -qF "$needle" "$TMP/dbg.err"; then
        ok "LEGO_DEBUG output contains '$needle'"
    else
        fail "LEGO_DEBUG output missing '$needle'"
    fi
done

echo
echo "----------------------------------------------------------------------"
echo "Results: $PASS passed, $FAIL failed"
[ "$FAIL" -eq 0 ] || exit 1
