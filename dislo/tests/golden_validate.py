#!/usr/bin/env python3
"""golden_validate.py — validate lego-dislo against the Disloelast goldens.

Runs the four golden-reference validations from the port specification:

  1. Solver, edge : P/F/G of W_edge_results, K*b^2/4pi = 0.1327641e+11
  2. Solver, screw: P/F/G of W_screw_results
  3. Field, edge  : example/W_45x_45y_1.5z.s.bcc vs
                    example/anisod_W_45x_45y_1.5z.s.bcc.bz2   (197200 atoms)
  4. Field, screw : W_x-12-1_y-101_z111.sml.bcc vs
                    anisoscrew_W_x-12-1_y-101_z111.sml.bcc    (330000 atoms)

The legacy reference data lives OUTSIDE the repository (it is too large to
commit); this script auto-skips with a notice when that folder is absent.

Root-order note: the three sextic roots are compared after nearest-root
matching.  The 2005 x87 binary's Laguerre deflation order is not bit-
reproducible on modern hardware; for the screw case the roots come out in
a cyclically permuted order.  The summed displacement field is invariant
under this permutation, which validations 3-4 prove directly.

Frame note (validation 4): the golden screw file was produced by the legacy
workflow applying the field of input.w_111screw (solver frame X=[-1,-1,2],
Y=[1,-1,0], Z=[1,1,1]) directly in the SAMPLE's coordinates, whose axes are
x=[-1,2,-1], y=[-1,0,1], z=[1,1,1].  Because a 120-degree rotation about
[111] is a cubic symmetry operation, the P/F/G of both frames coincide, and
the equivalent card uses the sample's own axes as the dislocation frame
(glide_normal = [-1,0,1]).

Standard library only.  Exit code 0 = all pass (or skipped), 1 = failure.

License: GPL-3.0
Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
         Max-Planck-Institut fuer Nachhaltige Materialien, Duesseldorf
Funding: NFDI-MatWerk
"""

import bz2
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(HERE, "..", "lego-dislo")
GOLDEN_DIR = os.environ.get("DISLOELAST_DIR",
                            os.path.expanduser("~/DEVEL/Disloelast"))

FAILURES = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}: {name}" + (f"  ({detail})" if detail
                                                     else ""))
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# Legacy-format readers
# ---------------------------------------------------------------------------

def read_legacy_results(path):
    """Read a legacy channel-21 results file (W_edge_results et al.)."""
    with open(path) as f:
        lines = [ln for ln in f if ln.strip()]
    # 5 header lines, then core/ANG, then 3 P, 3 F, 3 G
    vals = [[float(t) for t in ln.split()] for ln in lines[6:15]]
    P = [(vals[i][0], vals[i][1]) for i in range(3)]
    F = [vals[3 + k] for k in range(3)]
    G = [vals[6 + k] for k in range(3)]
    return P, F, G


def read_legacy_atoms(path):
    """Tolerant reader for legacy raw-IMD atom files: plain lines
    'id type mass x y z', with the IMD #-header block sometimes appended
    at the END of the file.  Returns dict id -> (type, mass, x, y, z)."""
    opener = bz2.open if path.endswith(".bz2") else open
    atoms = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split()
            if len(t) < 6:
                continue
            atoms[int(t[0])] = (int(t[1]), float(t[2]),
                                float(t[3]), float(t[4]), float(t[5]))
    return atoms


def write_imd(path, atoms):
    """Write a proper IMD file (header first) that atomio can read."""
    with open(path, "w") as f:
        f.write("#F A 1 1 1 3 0 0\n#C number type mass x y z\n")
        f.write("#X 1000 0 0\n#Y 0 1000 0\n#Z 0 0 1000\n#E\n")
        for i in sorted(atoms):
            t, m, x, y, z = atoms[i]
            f.write(f"{i} {t} {m:.15g} {x:.15g} {y:.15g} {z:.15g}\n")


# ---------------------------------------------------------------------------
# Report parsing and comparison helpers
# ---------------------------------------------------------------------------

def parse_report(text):
    """Parse the lego-dislo report into a list of per-dislocation dicts."""
    dislos = []
    cur = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("dislocation ") and ":" not in s:
            cur = {"P": [], "F": [], "G": []}
            dislos.append(cur)
        elif cur is not None and ":" in s:
            key, val = s.split(":", 1)
            vals = [float(t) for t in val.split()]
            if key.startswith("P_"):
                cur["P"].append((vals[0], vals[1]))
            elif key.startswith("F_"):
                cur["F"].append(vals)
            elif key.startswith("G_"):
                cur["G"].append(vals)
            else:
                cur[key] = vals if len(vals) > 1 else vals[0]
    return dislos


def run_solver(card_text, label):
    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as f:
        f.write(card_text)
        card = f.name
    try:
        r = subprocess.run([BINARY, "--solve-only", card],
                           capture_output=True, text=True)
        if r.returncode != 0:
            check(f"{label}: solver run", False, r.stderr.strip())
            return None
        return parse_report(r.stdout)[0]
    finally:
        os.unlink(card)


def rel_dev(a, b):
    scale = max(abs(a), abs(b), 1e-3)   # legacy prints 7 digits; entries
    return abs(a - b) / scale           # below ~1e-3 are noise-dominated


def compare_solution(label, mine, path_golden, tol=2e-5):
    P_g, F_g, G_g = read_legacy_results(path_golden)
    P_m, F_m, G_m = mine["P"], mine["F"], mine["G"]

    # Nearest-root bijective matching golden index -> our index.
    perm = []
    for pg in P_g:
        dists = [math.hypot(pg[0] - pm[0], pg[1] - pm[1]) for pm in P_m]
        perm.append(dists.index(min(dists)))
    if sorted(perm) != [0, 1, 2]:
        check(f"{label}: root matching is bijective", False, str(perm))
        return
    check(f"{label}: root matching", True,
          "identity order" if perm == [0, 1, 2]
          else f"permuted order {perm} (deflation-order difference)")

    dev = 0.0
    for n in range(3):
        dev = max(dev, rel_dev(P_g[n][0], P_m[perm[n]][0]),
                  rel_dev(P_g[n][1], P_m[perm[n]][1]))
        for k in range(3):
            dev = max(dev, rel_dev(F_g[k][n], F_m[k][perm[n]]),
                      rel_dev(G_g[k][n], G_m[k][perm[n]]))
    check(f"{label}: P/F/G max relative deviation < {tol:g}", dev < tol,
          f"max rel dev = {dev:.2e}")


def compare_fields(label, mine, golden, tol=1e-3):
    assert len(mine) == len(golden), (len(mine), len(golden))
    mx, ss, mxid = 0.0, 0.0, None
    for i, (_, _, x, y, z) in mine.items():
        _, _, gx, gy, gz = golden[i]
        d = math.sqrt((x - gx) ** 2 + (y - gy) ** 2 + (z - gz) ** 2)
        ss += d * d
        if d > mx:
            mx, mxid = d, i
    rms = math.sqrt(ss / len(mine))
    check(f"{label}: max|dr| < {tol:g} A", mx < tol,
          f"n={len(mine)} max|dr|={mx:.3e} A (atom {mxid}), rms={rms:.3e} A")
    return mx, rms


def run_field(card_text, in_atoms, label):
    tmpdir = tempfile.mkdtemp(prefix="lego-dislo-golden.")
    card = os.path.join(tmpdir, "card.yaml")
    fin = os.path.join(tmpdir, "in.imd")
    fout = os.path.join(tmpdir, "out.imd")
    with open(card, "w") as f:
        f.write(card_text)
    write_imd(fin, in_atoms)
    r = subprocess.run([BINARY, "--out-format=imd", card, fin, fout],
                       capture_output=True, text=True)
    if r.returncode != 0:
        check(f"{label}: field run", False, r.stderr.strip())
        return None
    out = read_legacy_atoms(fout)
    for p in (card, fin, fout):
        os.unlink(p)
    os.rmdir(tmpdir)
    return out


# ---------------------------------------------------------------------------
# Cards
# ---------------------------------------------------------------------------

CARD_HEAD = """\
material:
  name: W
  lattice_parameter: 3.1652
  cij:
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
"""

CARD_EDGE = CARD_HEAD + """\
box_orientation:
  x: [1, 1, 1]
  y: [-1, 0, 1]
  z: [1, -2, 1]
dislocations:
  - point: [0.0, 0.0, 0.0]
    line: [1, -2, 1]
    burgers_direction: [-1, -1, -1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
"""

CARD_SCREW_SOLVER = CARD_HEAD + """\
box_orientation:
  x: [-1, -1, 2]
  y: [1, -1, 0]
  z: [1, 1, 1]
dislocations:
  - point: [0.0, 0.0, 0.0]
    line: [1, 1, 1]
    burgers_direction: [1, 1, 1]
    burgers_magnitude: 0.8663
    glide_normal: [1, -1, 0]
"""

# See the frame note in the module docstring.
CARD_SCREW_FIELD = CARD_HEAD + """\
box_orientation:
  x: [-1, 2, -1]
  y: [-1, 0, 1]
  z: [1, 1, 1]
dislocations:
  - point: [0.0, 0.0, 0.0]
    line: [1, 1, 1]
    burgers_direction: [1, 1, 1]
    burgers_magnitude: 0.8663
    glide_normal: [-1, 0, 1]
"""


def main():
    if not os.path.isdir(GOLDEN_DIR):
        print(f"golden_validate.py: SKIP — legacy reference folder not found "
              f"at {GOLDEN_DIR}\n(set DISLOELAST_DIR to override)")
        return 0
    if not os.path.exists(BINARY):
        print(f"golden_validate.py: ERROR — {BINARY} not built (run make)")
        return 1

    print(f"golden data: {GOLDEN_DIR}")

    # 1. Solver, edge -------------------------------------------------------
    print("[1] solver, W 1/2<111> edge")
    sol = run_solver(CARD_EDGE, "edge solver")
    if sol:
        compare_solution("edge solver", sol,
                         os.path.join(GOLDEN_DIR, "W_edge_results"))
        # K*b^2/4pi: golden 0.1327641e+11 in units of 1e11 N/m^2 * latpar^2;
        # the report gives GPa*latpar^2 -> multiply by 1e9 for N/m^2.
        k_si = sol["K_b2_over_4pi_GPa_lp2"] * 1e9
        check("edge solver: K*b^2/4pi = 0.1327641e+11",
              rel_dev(k_si, 0.1327641e11) < 1e-6,
              f"got {k_si:.7e}")

    # 2. Solver, screw ------------------------------------------------------
    print("[2] solver, W 1/2<111> screw")
    sol = run_solver(CARD_SCREW_SOLVER, "screw solver")
    if sol:
        compare_solution("screw solver", sol,
                         os.path.join(GOLDEN_DIR, "W_screw_results"))

    # 3. Field, edge --------------------------------------------------------
    print("[3] field, edge (197200 atoms)")
    fin = os.path.join(GOLDEN_DIR, "example", "W_45x_45y_1.5z.s.bcc")
    fgold = os.path.join(GOLDEN_DIR, "example",
                         "anisod_W_45x_45y_1.5z.s.bcc.bz2")
    if os.path.exists(fin) and os.path.exists(fgold):
        out = run_field(CARD_EDGE, read_legacy_atoms(fin), "edge field")
        if out:
            compare_fields("edge field", out, read_legacy_atoms(fgold))
    else:
        print("  SKIP: edge sample/golden files not found")

    # 4. Field, screw -------------------------------------------------------
    print("[4] field, screw (330000 atoms)")
    fin = os.path.join(GOLDEN_DIR, "W_x-12-1_y-101_z111.sml.bcc")
    fgold = os.path.join(GOLDEN_DIR,
                         "anisoscrew_W_x-12-1_y-101_z111.sml.bcc")
    if os.path.exists(fin) and os.path.exists(fgold):
        out = run_field(CARD_SCREW_FIELD, read_legacy_atoms(fin),
                        "screw field")
        if out:
            compare_fields("screw field", out, read_legacy_atoms(fgold))
    else:
        print("  SKIP: screw sample/golden files not found")

    if FAILURES:
        print(f"\ngolden_validate.py: {len(FAILURES)} FAILURE(S)")
        return 1
    print("\ngolden_validate.py: all golden validations passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
