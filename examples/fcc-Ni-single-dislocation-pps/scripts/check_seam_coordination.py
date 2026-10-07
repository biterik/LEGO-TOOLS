#!/usr/bin/env python3
"""Independent check of a written LAMMPS data file (orthogonal or triclinic, 'p p s'): every atom
within 3 A of the periodic boundary along the GLIDE direction and more than 6 A away from the free
z surfaces must have exactly 12 neighbours closer than 3.0 A (fcc, a0 ~ 3.5 A).
usage: check_seam_coordination.py DATA.lmp glide-axis(x|y)      exit status 1 on failure"""
import sys
import numpy as np

f, ax = sys.argv[1], sys.argv[2]
L = open(f).read().split('\n')
box = {}; xy = 0.0
for l in L[:40]:
    for k in 'xyz':
        if f'{k}lo {k}hi' in l: box[k] = list(map(float, l.split()[:2]))
    if l.endswith('xy xz yz'): xy = float(l.split()[0])
i = next(k for k, l in enumerate(L) if l.strip().startswith('Atoms'))
r = np.array([l.split()[2:5] for l in L[i + 2:] if l.strip()], float)
lx, ly = box['x'][1] - box['x'][0], box['y'][1] - box['y'][0]
a, b = np.array([lx, 0, 0]), np.array([xy, ly, 0])
g = 0 if ax == 'x' else 1
s = r[:, g] - box['xy'[g]][0]; Lg = lx if g == 0 else ly
zlo, zhi = r[:, 2].min(), r[:, 2].max()
seam = (r[:, 2] > zlo + 6) & (r[:, 2] < zhi - 6) & ((s < 3) | (s > Lg - 3))
near = (s < 7) | (s > Lg - 7)
S, A = r[seam], r[near]
imgs = [m * a + k * b for m in (-1, 0, 1) for k in (-1, 0, 1)]
cnt = np.empty(len(S), int)
for j in range(len(S)):
    d = A - S[j]
    dmin = np.min([np.linalg.norm(d + im, axis=1) for im in imgs], axis=0)
    cnt[j] = ((dmin > 1e-6) & (dmin < 3.0)).sum()
h = dict(zip(*np.unique(cnt, return_counts=True)))
ok = set(h) == {12}
print(f'{f}: seam coordination {{neighbours: atoms}} = { {int(k): int(v) for k, v in h.items()} } -> {"PASS" if ok else "FAIL"}')
sys.exit(0 if ok else 1)
