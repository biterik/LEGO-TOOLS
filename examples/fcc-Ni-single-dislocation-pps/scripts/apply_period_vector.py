#!/usr/bin/env python3
"""Impose the corrected x-period vector (Lx', delta) on a cut single-dislocation slab and write a
LAMMPS triclinic data file.

Why a rotation: the corrected period vector along the glide direction has a component ALONG THE
LINE, a = (Lx', delta, 0) with delta = +-b_s/2 (Rodney 2004: "shift of b/2 in Y for atoms leaving
through X"). LAMMPS only tilts the SECOND box vector (b = (xy, ly, 0)), i.e. it can shift x across
the y boundary but not y across the x boundary. 'change_box all xy final delta' would therefore shear
the crystal by delta/Ly instead of delta/Lx -- wrong axis, ~Lx/Ly times too much (tested: screw cell,
sigma_xy went from -0.7 to +2.7 GPa). So the cell is written in the LAMMPS frame
    X = line direction (old y),  Y = -glide direction (old -x, shifted into [0, Lx')),  Z = old z
in which the old a becomes the LAMMPS b vector: (xy, ly) = (-delta, Lx').

Steps: affine x -> x*Lx'/Lx, y -> y + delta*x'/Lx' (keeps the seam a lattice match), then rotate.
Stress components back in the old (glide, line, normal) frame:
    s_xx = S_YY, s_yy = S_XX, s_zz = S_ZZ, s_xy = -S_XY, s_xz = -S_YZ, s_yz = S_XZ
usage: apply_period_vector.py CUT.lmp BOX-CORRECTION.lmp OUT.lmp
"""
import re, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from gate_single_dislocation_pps_cell import load_lammps_data   # noqa: E402

cut, corr, out = sys.argv[1:4]
v = {m.group(1): float(m.group(2)) for m in
     (re.match(r'variable\s+(\S+)\s+equal\s+(\S+)', l) for l in open(corr)) if m}
box, ids, r = load_lammps_data(cut)
Lx = box['x'][1] - box['x'][0]; Ly = box['y'][1] - box['y'][0]
Lxp, delta = v['Lx_corr'], v['xy_corr']
x1 = (r[:, 0] - box['x'][0]) * Lxp / Lx
y1 = (r[:, 1] - box['y'][0]) + delta * x1 / Lxp
X, Y, Z = y1, Lxp - x1, r[:, 2]
xy = -delta
# wrap into the parallelepiped (fractional coordinates along a = (Ly,0,0), b = (xy, Lxp, 0))
# (fractional coordinates from the UNWRAPPED positions first, then wrap both: wrapping Y alone
#  without the matching xy shift of X breaks the seam -- caught by the coordination check below)
sb = Y / Lxp
sa = (X - xy * sb) / Ly
sa, sb = np.mod(sa, 1.0), np.mod(sb, 1.0)
X, Y = sa * Ly + xy * sb, sb * Lxp
with open(out, 'w') as fh:
    fh.write(f'# {cut} with period vector (Lx\'={Lxp:.10f}, delta={delta:+.10f}); '
             f'LAMMPS frame X=line, Y=-glide, Z=normal\n\n{len(ids)} atoms\n1 atom types\n\n'
             f'0.0 {Ly:.10f} xlo xhi\n0.0 {Lxp:.10f} ylo yhi\n{box["z"][0]:.10f} {box["z"][1]:.10f} zlo zhi\n'
             f'{xy:.10f} 0.0 0.0 xy xz yz\n\nAtoms # atomic\n\n')
    for i, a, b, c in zip(ids, X, Y, Z):
        fh.write(f'{i} 1 {a:.10f} {b:.10f} {c:.10f}\n')
print(f'wrote {out}: lx {Ly:.4f} (line)  ly {Lxp:.4f} (glide, corrected)  xy {xy:+.4f}  atoms {len(ids)}')
