#!/usr/bin/env bash
# Single a/2<110>{111} dislocation (screw, 30, 60, edge) in an fcc Ni slab, LAMMPS 'p p s':
# LEGO -> dcreator -> lego-cut -> gate -> LAMMPS relaxation (box as built AND corrected)
# -> full stress-tensor maps and far-field profiles.
#
# Needs: lego (github.com/biterik/LEGO), dcreator (github.com/biterik/dcreator), lego-cut (this
# repo), Python 3 with numpy + matplotlib, and LAMMPS with MANYBODY (Python module 'lammps' or
# $LMP pointing to an lmp binary), plus an eam/alloy file containing Ni.
#
#   POT=/path/to/Cu_Ni_Fischer_2018.eam.alloy ./run_example.sh            # all four characters
#   POT=... ./run_example.sh edge 60                                      # a subset
#
# Run time (2 cores, LAMMPS serial per relaxation): ~3 min per relaxation, 8 relaxations.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
: "${POT:?set POT to an eam/alloy file containing Ni (e.g. Fischer 2018 Cu-Ni from the NIST IPR)}"
ELEM="${ELEM:-Ni}"
A0="${A0:-3.52}"                   # 0 K lattice constant of the potential (Fischer 2018 Ni: 3.5200)
OUT="${OUT:-$HERE/out}"
CHARS=("${@:-screw 30 60 edge}")
read -r -a CHARS <<< "${CHARS[*]}"
export LEGO="${LEGO:-$(command -v lego)}" DCREATOR="${DCREATOR:-$(command -v dcreator)}"
export LEGO_CUT="${LEGO_CUT:-$HERE/../../lego-tools/lego-cut}"

names=()
for c in "${CHARS[@]}"; do
  python3 "$HERE/scripts/build_single_dislocation_cell.py" "$c" "$OUT/cells" --a0 "$A0"
  case $c in screw|60) o=x-211_y0-11_z111;; *) o=x-101_y1-21_z111;; esac
  suffix=""; if [[ $c =~ ^[0-9]+$ ]]; then suffix="deg"; fi
  n="${ELEM}_${c}${suffix}_${o}"; names+=("$n")
  for corr in 0 1; do
    python3 "$HERE/scripts/run_lammps_relax.py" "$n" "$OUT/cells" "$OUT/relaxed" "$POT" "$ELEM" "$corr"
  done
done
python3 "$HERE/scripts/plot_stress_tensor_maps.py" "$OUT/cells" "$OUT/relaxed" "$OUT/figures" "${names[@]}" --a0 "$A0"
