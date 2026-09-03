#!/usr/bin/env bash
#
# lego-dislo test suite
# =====================
# Runs, in order:
#   1. unit_tests.py        — self-contained (parser, frames, superposition,
#                             rotation equivalence); always runs
#   2. golden_validate.py   — validation against the legacy Disloelast
#                             golden files; auto-skips with a notice when
#                             the reference folder (~/DEVEL/Disloelast, or
#                             $DISLOELAST_DIR) is absent
#
# Exit code: 0 on success (skips count as success), non-zero on failure.
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
# Funding: NFDI-MatWerk

set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"

if [ ! -x "$HERE/../lego-dislo" ]; then
    echo "ERROR: $HERE/../lego-dislo not built — run 'make' at the repo root" >&2
    exit 2
fi

echo "lego-dislo tests"
echo "================"
echo
echo "--- unit tests -------------------------------------------------------"
python3 "$HERE/unit_tests.py"
echo
echo "--- golden validation ------------------------------------------------"
python3 "$HERE/golden_validate.py"
