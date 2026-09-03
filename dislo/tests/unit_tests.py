#!/usr/bin/env python3
"""unit_tests.py — self-contained tests for lego-dislo (no legacy data).

Covers:
  - cubic shorthand vs full 6x6 Voigt matrix give identical results
  - frame validation errors (non-orthogonal box, glide normal not
    perpendicular to the line, left-handed box)
  - card parser edge cases (missing keys, bad lists, tabs, unknown keys,
    asymmetric full matrix, comments)
  - superposition: a zero-Burgers second dislocation leaves the output
    byte-identical to the single-dislocation run (the multi-dislocation
    path adds exactly zero), and a dipole of opposite dislocations has
    decaying far-field displacements
  - rotation equivalence: the same physical dislocation inserted through
    two different box_orientation frames produces displacement fields that
    match after rotation, to 1e-8

Standard library only.  Exit code 0 = all pass, 1 = failure.

License: GPL-3.0
Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
         Max-Planck-Institut fuer Nachhaltige Materialien, Duesseldorf
Funding: NFDI-MatWerk
"""

import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "..", "lego-dislo")

FAILURES = []
TMP = tempfile.mkdtemp(prefix="lego-dislo-unit.")


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}: {name}" + (f"  ({detail})" if detail
                                                     else ""))
    if not ok:
        FAILURES.append(name)


def write(path, text):
    p = os.path.join(TMP, path)
    with open(p, "w") as f:
        f.write(text)
    return p


def run(args):
    return subprocess.run([BINARY] + args, capture_output=True, text=True)


CARD_HEAD = """\
material:
  name: W
  lattice_parameter: 3.1652
  cij:
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
"""

BOX_EDGE = """\
box_orientation:
  x: [1, 1, 1]
  y: [-1, 0, 1]
  z: [1, -2, 1]
"""

DISLO_EDGE = """\
dislocations:
  - point: [0.1234, 0.0567, 0.0]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
"""


# ---------------------------------------------------------------------------
# Atom-file helpers (LAMMPS data files, written/read directly)
# ---------------------------------------------------------------------------

def write_points(path, points, lo=-60.0, hi=60.0):
    p = os.path.join(TMP, path)
    with open(p, "w") as f:
        f.write("test points\n\n")
        f.write(f"{len(points)} atoms\n1 atom types\n\n")
        f.write(f"{lo} {hi} xlo xhi\n{lo} {hi} ylo yhi\n{lo} {hi} zlo zhi\n\n")
        f.write("Masses\n\n1 183.84\n\nAtoms # atomic\n\n")
        for i, (x, y, z) in enumerate(points, 1):
            f.write(f"{i} 1 {x:.12f} {y:.12f} {z:.12f}\n")
    return p


def read_positions(path):
    """Read id -> (x,y,z) from a LAMMPS data file."""
    pos = {}
    in_atoms = False
    with open(path) as f:
        for line in f:
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
                pos[int(t[0])] = (float(t[2]), float(t[3]), float(t[4]))
    return pos


def lcg_points(n, span, seed=12345):
    """Deterministic pseudo-random points in [-span, span]^3."""
    state = seed
    pts = []

    def rnd():
        nonlocal state
        state = (6364136223846793005 * state + 1442695040888963407) % 2**64
        return state / 2**64

    while len(pts) < n:
        p = tuple(span * (2 * rnd() - 1) for _ in range(3))
        pts.append(p)
    return pts


# ---------------------------------------------------------------------------
# 1. cubic vs full 6x6
# ---------------------------------------------------------------------------

def test_cubic_vs_full():
    print("[cubic vs full 6x6]")
    c11, c12, c44 = 522.4, 204.4, 160.6
    rows = []
    for i in range(6):
        row = [0.0] * 6
        if i < 3:
            for j in range(3):
                row[j] = c11 if i == j else c12
        else:
            row[i] = c44
        rows.append("      - [" + ", ".join(f"{v}" for v in row) + "]")
    card_full = ("material:\n  name: W\n  lattice_parameter: 3.1652\n"
                 "  cij:\n    full:\n" + "\n".join(rows) + "\n"
                 + BOX_EDGE + DISLO_EDGE)
    card_cubic = CARD_HEAD + BOX_EDGE + DISLO_EDGE

    r1 = run(["--solve-only", write("c_cubic.yaml", card_cubic)])
    r2 = run(["--solve-only", write("c_full.yaml", card_full)])
    check("both cards solve", r1.returncode == 0 and r2.returncode == 0,
          (r1.stderr + r2.stderr).strip())
    check("reports are identical", r1.stdout == r2.stdout)


# ---------------------------------------------------------------------------
# 2. validation errors
# ---------------------------------------------------------------------------

def test_validation_errors():
    print("[frame validation and parser errors]")
    cases = [
        ("box axes not orthogonal",
         CARD_HEAD + "box_orientation:\n  x: [1, 1, 1]\n  y: [1, 0, 1]\n"
         "  z: [1, -2, 1]\n" + DISLO_EDGE, "orthogonal"),
        ("box not right-handed",
         CARD_HEAD + "box_orientation:\n  x: [1, 1, 1]\n  y: [-1, 0, 1]\n"
         "  z: [-1, 2, -1]\n" + DISLO_EDGE, "right-handed"),
        ("glide normal not perpendicular to line",
         CARD_HEAD + BOX_EDGE +
         "dislocations:\n  - point: [0, 0, 0]\n    line: [1, -2, 1]\n"
         "    burgers_direction: [-1, -1, -1]\n"
         "    burgers_magnitude: 0.8663\n    glide_normal: [1, 0, 1]\n",
         "perpendicular"),
        ("missing required dislocation key",
         CARD_HEAD + BOX_EDGE +
         "dislocations:\n  - point: [0, 0, 0]\n    line: [1, -2, 1]\n"
         "    burgers_magnitude: 0.8663\n    glide_normal: [-1, 0, 1]\n",
         "missing"),
        ("missing lattice parameter",
         "material:\n  name: W\n  cij:\n"
         "    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}\n"
         + BOX_EDGE + DISLO_EDGE, "lattice_parameter"),
        ("bad list", CARD_HEAD +
         "box_orientation:\n  x: [1, 1]\n  y: [-1, 0, 1]\n  z: [1, -2, 1]\n"
         + DISLO_EDGE, ""),
        ("tab indentation", CARD_HEAD.replace("  name", "\tname")
         + BOX_EDGE + DISLO_EDGE, "tab"),
        ("unknown key", CARD_HEAD + "box_orientation:\n  x: [1, 1, 1]\n"
         "  y: [-1, 0, 1]\n  z: [1, -2, 1]\n  w: [1, 0, 0]\n" + DISLO_EDGE,
         "unknown"),
        ("asymmetric full matrix",
         "material:\n  name: X\n  lattice_parameter: 3.0\n  cij:\n    full:\n"
         "      - [100, 50, 50, 0, 0, 0]\n      - [51, 100, 50, 0, 0, 0]\n"
         "      - [50, 50, 100, 0, 0, 0]\n      - [0, 0, 0, 30, 0, 0]\n"
         "      - [0, 0, 0, 0, 30, 0]\n      - [0, 0, 0, 0, 0, 30]\n"
         + BOX_EDGE + DISLO_EDGE, "symmetric"),
    ]
    for i, (name, card, needle) in enumerate(cases):
        r = run(["--solve-only", write(f"bad{i}.yaml", card)])
        ok = (r.returncode != 0 and
              (needle.lower() in r.stderr.lower() if needle else True))
        check(name + " is rejected", ok, r.stderr.strip()[:90])

    # Comments and blank lines are fine.
    r = run(["--solve-only",
             write("comments.yaml",
                   "# a comment\n" + CARD_HEAD + "\n" + BOX_EDGE
                   + DISLO_EDGE + "    cut_angle: 0.0  # trailing comment\n")])
    check("comments and blank lines are accepted", r.returncode == 0,
          r.stderr.strip()[:90])


# ---------------------------------------------------------------------------
# 3. superposition
# ---------------------------------------------------------------------------

def test_superposition():
    print("[superposition]")
    pts = lcg_points(400, 55.0)
    fin = write_points("sup_in.data", pts)

    card1 = CARD_HEAD + BOX_EDGE + DISLO_EDGE
    card2 = CARD_HEAD + BOX_EDGE + DISLO_EDGE + """\
  - point: [200.0, 100.0, 0.0]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.0
    glide_normal: [-1, 0, 1]
"""
    out1 = os.path.join(TMP, "sup_out1.data")
    out2 = os.path.join(TMP, "sup_out2.data")
    r1 = run([write("sup1.yaml", card1), fin, out1])
    r2 = run([write("sup2.yaml", card2), fin, out2])
    check("runs succeed", r1.returncode == 0 and r2.returncode == 0,
          (r1.stderr + r2.stderr).strip()[:90])
    with open(out1, "rb") as f1, open(out2, "rb") as f2:
        check("zero-b second dislocation leaves output byte-identical",
              f1.read() == f2.read())

    # Dipole far-field decay: opposite dislocations 10 A apart vs a single
    # one.  Displacements are measured on the same input points.
    card_dip = CARD_HEAD + BOX_EDGE + """\
dislocations:
  - point: [-5.0, 0.03, 0.0]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
  - point: [5.0, 0.03, 0.0]
    line: [1, -2, 1]
    burgers_direction: [1, 1, 1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
"""
    outd = os.path.join(TMP, "sup_outd.data")
    rd = run([write("supd.yaml", card_dip), fin, outd])
    check("dipole run succeeds", rd.returncode == 0, rd.stderr.strip()[:90])

    pin = read_positions(fin)
    psgl = read_positions(out1)
    pdip = read_positions(outd)

    def u(pos, i):
        return tuple(pos[i][k] - pin[i][k] for k in range(3))

    far = [i for i in pin
           if math.hypot(pin[i][0], pin[i][1]) > 40.0]
    near = [i for i in pin
            if 15.0 < math.hypot(pin[i][0], pin[i][1]) < 25.0]
    check("test geometry has near and far atoms",
          len(far) > 20 and len(near) > 20,
          f"near={len(near)} far={len(far)}")

    def maxu(pos, ids):
        return max(math.sqrt(sum(c * c for c in u(pos, i))) for i in ids)

    mf_dip, mn_dip = maxu(pdip, far), maxu(pdip, near)
    mf_sgl = maxu(psgl, far)
    check("dipole far field decays (far < near)", mf_dip < 0.8 * mn_dip,
          f"near={mn_dip:.3f} A, far={mf_dip:.3f} A")
    check("dipole far field is small vs single dislocation",
          mf_dip < 0.5 * mf_sgl,
          f"dipole={mf_dip:.3f} A, single={mf_sgl:.3f} A")


# ---------------------------------------------------------------------------
# 4. rotation equivalence
# ---------------------------------------------------------------------------

def norm(v):
    s = math.sqrt(sum(c * c for c in v))
    return [c / s for c in v]


def matvec(m, v):
    return [sum(m[i][k] * v[k] for k in range(3)) for i in range(3)]


def test_rotation_equivalence():
    print("[rotation equivalence]")
    # Frame A: the edge frame; frame B: rotated (same crystal, same physical
    # dislocation, different box axes).
    A = [norm([1, 1, 1]), norm([-1, 0, 1]), norm([1, -2, 1])]
    B = [norm([-1, 0, 1]), norm([1, -2, 1]), norm([1, 1, 1])]
    # Q maps A-frame coordinates to B-frame coordinates: Q = R_B * R_A^T.
    Q = [[sum(B[i][k] * A[j][k] for k in range(3)) for j in range(3)]
         for i in range(3)]

    ptsA = []
    core_A = [0.1234, 0.0567, 0.0]
    for p in lcg_points(500, 45.0, seed=999):
        # avoid the branch-cut half-plane y=0, x>0 (A frame == dislo frame)
        if abs(p[1] - core_A[1]) < 1.0 and p[0] - core_A[0] > -1.0:
            continue
        # and the core region
        if math.hypot(p[0] - core_A[0], p[1] - core_A[1]) < 3.0:
            continue
        ptsA.append(p)
    ptsB = [matvec(Q, p) for p in ptsA]
    core_B = matvec(Q, core_A)

    dislo = """\
dislocations:
  - point: [{p[0]:.12f}, {p[1]:.12f}, {p[2]:.12f}]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
"""
    cardA = CARD_HEAD + BOX_EDGE + dislo.format(p=core_A)
    cardB = (CARD_HEAD + "box_orientation:\n  x: [-1, 0, 1]\n"
             "  y: [1, -2, 1]\n  z: [1, 1, 1]\n" + dislo.format(p=core_B))

    finA = write_points("rotA_in.data", ptsA)
    finB = write_points("rotB_in.data", ptsB)
    outA = os.path.join(TMP, "rotA_out.data")
    outB = os.path.join(TMP, "rotB_out.data")
    rA = run([write("rotA.yaml", cardA), finA, outA])
    rB = run([write("rotB.yaml", cardB), finB, outB])
    check("both frames run", rA.returncode == 0 and rB.returncode == 0,
          (rA.stderr + rB.stderr).strip()[:90])

    pA = read_positions(outA)
    pB = read_positions(outB)
    mx = 0.0
    for i, p0 in enumerate(ptsA, 1):
        uA = [pA[i][k] - p0[k] for k in range(3)]
        uB = [pB[i][k] - ptsB[i - 1][k] for k in range(3)]
        uAr = matvec(Q, uA)
        d = math.sqrt(sum((uB[k] - uAr[k]) ** 2 for k in range(3)))
        mx = max(mx, d)
    check("rotated displacement fields match to 1e-8", mx < 1e-8,
          f"n={len(ptsA)} max deviation = {mx:.3e} A")


def main():
    if not os.path.exists(BINARY):
        print(f"unit_tests.py: ERROR — {BINARY} not built (run make)")
        return 1
    try:
        test_cubic_vs_full()
        test_validation_errors()
        test_superposition()
        test_rotation_equivalence()
    finally:
        shutil.rmtree(TMP, ignore_errors=True)
    if FAILURES:
        print(f"\nunit_tests.py: {len(FAILURES)} FAILURE(S)")
        return 1
    print("\nunit_tests.py: all tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
