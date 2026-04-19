#!/usr/bin/env python3
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk
"""
shift_lmps.py — shift a LAMMPS data file by a vector (dx,dy,dz),
updating both the box bounds (xlo/xhi, ylo/yhi, zlo/zhi) and atom coordinates.

Example:
  ./shift_lmps.py -1 -1 20 slab.data shifted_slab.data

Notes:
- Works for orthogonal and triclinic boxes (keeps xy/xz/yz tilt factors unchanged).
- By default assumes the "Atoms" section has x,y,z in columns 3–5 (1-based),
  i.e. tokens[2], tokens[3], tokens[4] (0-based), which matches "Atoms # atomic"
  like: id type x y z [ix iy iz]
- If your atom style differs, use --xyz-cols to specify 0-based indices.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Tuple, Optional


RE_BOUNDS = re.compile(
    r"""^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+
         ([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+
         (xlo\s+xhi|ylo\s+yhi|zlo\s+zhi)\s*$""",
    re.VERBOSE,
)

RE_TILT = re.compile(
    r"""^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+
         ([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+
         ([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s+
         (xy\s+xz\s+yz)\s*$""",
    re.VERBOSE,
)

RE_SECTION_HEADER = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\b")


def fmt_float(x: float, spec: str) -> str:
    # Use a stable float format; user can override with --float-format
    return format(x, spec)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Shift LAMMPS data file (box + atoms) by dx dy dz.")
    p.add_argument("dx", type=float, help="Shift in x")
    p.add_argument("dy", type=float, help="Shift in y")
    p.add_argument("dz", type=float, help="Shift in z")
    p.add_argument("input", type=Path, help="Input LAMMPS data file")
    p.add_argument("output", type=Path, help="Output (shifted) LAMMPS data file")

    p.add_argument(
        "--xyz-cols",
        nargs=3,
        type=int,
        default=[2, 3, 4],
        metavar=("XCOL", "YCOL", "ZCOL"),
        help=(
            "0-based token indices for x,y,z in the Atoms section line. "
            "Default: 2 3 4 (i.e. columns 3–5 in 1-based counting)."
        ),
    )
    p.add_argument(
        "--float-format",
        default=".16g",
        help="Python format spec for floats written back (default: .16g).",
    )
    return p.parse_args()


def is_new_section_line(line: str) -> bool:
    # Any line starting with a word (letters) is considered a new section header,
    # e.g. 'Velocities', 'Bonds', 'Masses', etc.
    return bool(RE_SECTION_HEADER.match(line))


def shift_bounds_line(line: str, dx: float, dy: float, dz: float, float_spec: str) -> Optional[str]:
    m = RE_BOUNDS.match(line)
    if not m:
        return None
    lo = float(m.group(1))
    hi = float(m.group(2))
    kind = m.group(3).strip()

    if kind == "xlo xhi":
        lo += dx
        hi += dx
    elif kind == "ylo yhi":
        lo += dy
        hi += dy
    elif kind == "zlo zhi":
        lo += dz
        hi += dz
    else:
        return None

    return f"{fmt_float(lo, float_spec)} {fmt_float(hi, float_spec)} {kind}\n"


def try_shift_atom_line(
    line: str,
    dx: float,
    dy: float,
    dz: float,
    xyz_cols: Tuple[int, int, int],
    float_spec: str,
) -> Optional[str]:
    s = line.strip()
    if not s or s.startswith("#"):
        return None

    parts = s.split()
    xcol, ycol, zcol = xyz_cols

    # Must have enough tokens
    if max(xcol, ycol, zcol) >= len(parts):
        return None

    # Try parsing x,y,z as floats
    try:
        x = float(parts[xcol]) + dx
        y = float(parts[ycol]) + dy
        z = float(parts[zcol]) + dz
    except ValueError:
        return None

    parts[xcol] = fmt_float(x, float_spec)
    parts[ycol] = fmt_float(y, float_spec)
    parts[zcol] = fmt_float(z, float_spec)

    # Preserve a single-space separation; LAMMPS is fine with this.
    return " ".join(parts) + "\n"


def main() -> int:
    args = parse_args()
    dx, dy, dz = args.dx, args.dy, args.dz
    xyz_cols = (args.xyz_cols[0], args.xyz_cols[1], args.xyz_cols[2])
    float_spec = args.float_format

    lines = args.input.read_text(errors="replace").splitlines(keepends=True)

    out: List[str] = []
    in_atoms = False
    atoms_header_seen = False

    i = 0
    while i < len(lines):
        line = lines[i]

        # Update box bounds
        shifted = shift_bounds_line(line, dx, dy, dz, float_spec)
        if shifted is not None:
            out.append(shifted)
            i += 1
            continue

        # Keep triclinic tilt line as-is (xy xz yz)
        if RE_TILT.match(line):
            out.append(line)
            i += 1
            continue

        # Detect start of Atoms section
        if line.strip().startswith("Atoms"):
            out.append(line)
            in_atoms = True
            atoms_header_seen = True
            i += 1
            continue

        # If we are in Atoms section, shift coordinate lines until next section header
        if in_atoms:
            # LAMMPS data files typically have a blank line after the section title; keep it
            if not line.strip():
                out.append(line)
                i += 1
                continue

            # Stop if we hit a new section header (e.g. "Velocities", "Bonds", ...)
            # IMPORTANT: atom lines start with an integer id, not a word.
            if is_new_section_line(line) and not line.lstrip().startswith(("+", "-", ".", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")):
                in_atoms = False
                out.append(line)
                i += 1
                continue

            # Otherwise try to shift the atom line
            new_line = try_shift_atom_line(line, dx, dy, dz, xyz_cols, float_spec)
            out.append(new_line if new_line is not None else line)
            i += 1
            continue

        # Default: passthrough
        out.append(line)
        i += 1

    if not atoms_header_seen:
        raise SystemExit("ERROR: Did not find an 'Atoms' section header in the input file.")

    args.output.write_text("".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
