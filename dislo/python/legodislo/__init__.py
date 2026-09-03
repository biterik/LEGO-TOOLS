"""legodislo — Python driver for the lego-dislo tool (LEGO-TOOLS).

Inserts dislocations into atomistic configurations using the anisotropic
elastic solution (sextic method, Hirth & Lothe).  This module is a thin,
standard-library-only wrapper: it writes the input card, invokes the
compiled ``lego-dislo`` binary via subprocess, checks the return code, and
parses the solver report (per-dislocation P/F/G and K factors) into the
return value.

Example
-------
    from legodislo import Material, Dislocation, insert_dislocations

    result = insert_dislocations(
        input="W.data", output="W_edge.data",
        material=Material(name="W", a0=3.1652,
                          cubic=(522.4, 204.4, 160.6)),      # GPa
        box_orientation=((1, 1, 1), (-1, 0, 1), (1, -2, 1)),
        dislocations=[Dislocation(point=(0, 0, 0), line=(1, -2, 1),
                                  b_dir=(-1, -1, -1), b_mag=0.8663,
                                  glide_normal=(-1, 0, 1))])
    print(result.dislocations[0].K_GPa)

License: GPL-3.0
Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
         Max-Planck-Institut fuer Nachhaltige Materialien, Duesseldorf
Funding: NFDI-MatWerk
"""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

__all__ = ["Material", "Dislocation", "DisloResult", "Result",
           "insert_dislocations", "solve", "make_card", "find_binary",
           "LegoDisloError"]


class LegoDisloError(RuntimeError):
    """Raised when the lego-dislo binary cannot be found, fails, or its
    report cannot be parsed."""


Vec3 = Sequence[float]


@dataclass
class Material:
    """Elastic description of the crystal.

    Give either ``cubic=(c11, c12, c44)`` or a full 6x6 Voigt matrix
    ``cij`` (row-major, crystal frame).  Constants in GPa; ``a0`` in
    Angstrom.
    """
    a0: float
    name: str = "material"
    cubic: Optional[Tuple[float, float, float]] = None
    cij: Optional[Sequence[Sequence[float]]] = None

    def __post_init__(self):
        if (self.cubic is None) == (self.cij is None):
            raise ValueError("Material needs exactly one of cubic=(c11,c12,"
                             "c44) or cij=6x6 matrix")
        if self.cij is not None:
            if len(self.cij) != 6 or any(len(r) != 6 for r in self.cij):
                raise ValueError("Material.cij must be a 6x6 matrix")
        if self.a0 <= 0:
            raise ValueError("Material.a0 must be positive (Angstrom)")


@dataclass
class Dislocation:
    """One dislocation.

    point         a point on the line, box coordinates, Angstrom
    line          line direction, crystal coordinates (local Z axis)
    b_dir         Burgers vector direction, crystal coordinates
    b_mag         Burgers vector magnitude, lattice-parameter units
                  (length in Angstrom = b_mag * a0)
    glide_normal  glide-plane normal, crystal coordinates (local Y axis);
                  must be perpendicular to line
    cut_angle     branch-cut angle in degrees from the local +X axis
    """
    point: Vec3
    line: Vec3
    b_dir: Vec3
    b_mag: float
    glide_normal: Vec3
    cut_angle: float = 0.0


@dataclass
class DisloResult:
    """Parsed solver output for one dislocation."""
    P: List[Tuple[float, float]] = field(default_factory=list)
    F: List[List[float]] = field(default_factory=list)
    G: List[List[float]] = field(default_factory=list)
    burgers_vector_lp: Optional[List[float]] = None
    K_b2_over_4pi_GPa_lp2: Optional[float] = None
    K_b2_over_4pi_J_per_m: Optional[float] = None
    K_GPa: Optional[float] = None


@dataclass
class Result:
    """Return value of insert_dislocations / solve."""
    dislocations: List[DisloResult]
    report: str                  # full report text
    stderr: str                  # diagnostic output of the binary


def find_binary(binary: Optional[str] = None) -> str:
    """Locate the lego-dislo binary.

    Search order: explicit argument, $LEGO_DISLO, the build location
    relative to this package (dislo/lego-dislo), then $PATH.
    """
    candidates = []
    if binary:
        candidates.append(binary)
    if os.environ.get("LEGO_DISLO"):
        candidates.append(os.environ["LEGO_DISLO"])
    here = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(here, "..", "..", "lego-dislo"))
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return os.path.abspath(c)
    onpath = shutil.which("lego-dislo")
    if onpath:
        return onpath
    raise LegoDisloError(
        "cannot locate the lego-dislo binary; build it with 'make' in the "
        "LEGO-TOOLS repository, put it on $PATH, set $LEGO_DISLO, or pass "
        "binary=...")


def _fmt_vec(v: Vec3) -> str:
    return "[" + ", ".join(repr(float(c)) for c in v) + "]"


def make_card(material: Material,
              box_orientation: Sequence[Vec3],
              dislocations: Sequence[Dislocation]) -> str:
    """Render the input card (the exact flat YAML subset lego-dislo parses)."""
    if len(box_orientation) != 3:
        raise ValueError("box_orientation must be three direction vectors "
                         "(x, y, z)")
    if not dislocations:
        raise ValueError("at least one Dislocation is required")
    lines = ["material:",
             f"  name: {material.name}",
             f"  lattice_parameter: {material.a0!r}",
             "  cij:"]
    if material.cubic is not None:
        c11, c12, c44 = material.cubic
        lines.append("    cubic: {c11: %r, c12: %r, c44: %r}"
                     % (float(c11), float(c12), float(c44)))
    else:
        lines.append("    full:")
        for row in material.cij:
            lines.append("      - " + _fmt_vec(row))
    lines.append("box_orientation:")
    for name, v in zip(("x", "y", "z"), box_orientation):
        lines.append(f"  {name}: {_fmt_vec(v)}")
    lines.append("dislocations:")
    for d in dislocations:
        lines.append(f"  - point: {_fmt_vec(d.point)}")
        lines.append(f"    line: {_fmt_vec(d.line)}")
        lines.append(f"    burgers_direction: {_fmt_vec(d.b_dir)}")
        lines.append(f"    burgers_magnitude: {float(d.b_mag)!r}")
        lines.append(f"    glide_normal: {_fmt_vec(d.glide_normal)}")
        lines.append(f"    cut_angle: {float(d.cut_angle)!r}")
    return "\n".join(lines) + "\n"


def parse_report(text: str) -> List[DisloResult]:
    """Parse the lego-dislo report text into DisloResult objects."""
    out: List[DisloResult] = []
    cur: Optional[DisloResult] = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("dislocation ") and ":" not in s:
            cur = DisloResult()
            out.append(cur)
        elif cur is not None and ":" in s:
            key, val = s.split(":", 1)
            try:
                vals = [float(t) for t in val.split()]
            except ValueError:
                continue
            if not vals:
                continue
            if key.startswith("P_"):
                cur.P.append((vals[0], vals[1]))
            elif key.startswith("F_"):
                cur.F.append(vals)
            elif key.startswith("G_"):
                cur.G.append(vals)
            elif key == "burgers_vector_lp":
                cur.burgers_vector_lp = vals
            elif key in ("K_b2_over_4pi_GPa_lp2", "K_b2_over_4pi_J_per_m",
                         "K_GPa"):
                setattr(cur, key, vals[0])
    return out


def _invoke(options: List[str], tail: List[str], binary: Optional[str],
            material: Material, box_orientation: Sequence[Vec3],
            dislocations: Sequence[Dislocation]) -> Result:
    """Write the card, run 'lego-dislo OPTIONS CARD TAIL...', parse."""
    exe = find_binary(binary)
    card_text = make_card(material, box_orientation, dislocations)
    fd, card_path = tempfile.mkstemp(suffix=".yaml", prefix="legodislo.")
    try:
        with os.fdopen(fd, "w") as f:
            f.write(card_text)
        cmd = [exe] + options + [card_path] + tail
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise LegoDisloError(
                f"lego-dislo failed (exit {proc.returncode}):\n"
                f"{proc.stderr.strip()}\n--- card ---\n{card_text}")
        dislos = parse_report(proc.stdout)
        if len(dislos) != len(dislocations):
            raise LegoDisloError(
                f"could not parse the lego-dislo report (expected "
                f"{len(dislocations)} dislocations, found {len(dislos)}):\n"
                f"{proc.stdout}")
        return Result(dislocations=dislos, report=proc.stdout,
                      stderr=proc.stderr)
    finally:
        os.unlink(card_path)


def insert_dislocations(input: str, output: str, material: Material,
                        box_orientation: Sequence[Vec3],
                        dislocations: Sequence[Dislocation],
                        binary: Optional[str] = None,
                        out_format: Optional[str] = None,
                        threads: Optional[int] = None) -> Result:
    """Insert dislocations into the configuration `input`, write `output`.

    out_format: 'lammps' (default), 'lammps-dump', or 'imd'.
    Returns a Result with the parsed per-dislocation solver output.
    Raises LegoDisloError on any failure.
    """
    if not os.path.exists(input):
        raise LegoDisloError(f"input file not found: {input}")
    options = []
    if out_format:
        options.append(f"--out-format={out_format}")
    if threads:
        options.append(f"--threads={int(threads)}")
    return _invoke(options, [input, output], binary, material,
                   box_orientation, dislocations)


def solve(material: Material, box_orientation: Sequence[Vec3],
          dislocations: Sequence[Dislocation],
          binary: Optional[str] = None) -> Result:
    """Run the solver only (no atoms): returns P/F/G and K factors."""
    return _invoke(["--solve-only"], [], binary, material,
                   box_orientation, dislocations)
