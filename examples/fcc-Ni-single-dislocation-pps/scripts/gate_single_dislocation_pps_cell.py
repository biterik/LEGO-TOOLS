#!/usr/bin/env python3
"""Gate for a single straight dislocation in a slab periodic along the glide direction x and the
line y, free along the glide-plane normal z (LAMMPS 'p p s').

Rule: the box x-period vector must be the AVERAGE of the natural period vectors of the two
half-crystals above and below the glide plane,
    Lx' = (L_top + L_bot)/2,   L_half = N_half / (atoms per Angstrom of x of that half, perfect slab)
    delta = < U_y(low-x end) - U_y(high-x end) >  over the two halves  (dcreator rigid shift)
i.e. the x-period vector is a = (Lx', delta, 0): a y-offset ALONG THE LINE across the x boundary.
LAMMPS cannot express that with its xy tilt (which offsets x across the y boundary); use
apply_period_vector.py, which puts the line along LAMMPS x.
Afterwards each half carries the unavoidable +-b_e/(2Lx) normal and +-b_s/(2Lx) shear strain of
the periodic row of dislocations, with zero mean. Rodney, PRB 61, 8714 (2000): Lx + b/2 for the
edge; Rodney, Acta Mater 52, 607 (2004): shifted periodicity +b/2 along Y for atoms leaving through X.

Checks:  edge Burgers content L_bot - L_top (must be +-b_e) and the seam coordination after the
affine remap to (Lx', delta) (every interior atom within 3 A of the x seam must have 12 neighbours).

usage: gate_single_dislocation_pps_cell.py PERFECT INSERTED CUT z_glide x_line
"""
import sys
import numpy as np


def load_lammps_data(f):
    box, rows, in_atoms = {}, [], False
    with open(f) as fh:
        lines = fh.read().split('\n')
    for k, l in enumerate(lines):
        for ax in 'xyz':
            if f'{ax}lo {ax}hi' in l:
                box[ax] = list(map(float, l.split()[:2]))
        if l.strip().startswith('Atoms'):
            start = k + 2
            break
    for l in lines[start:]:
        s = l.split()
        if not s:
            if rows:
                break
            continue
        if s[0][0].isalpha():
            break
        rows.append(s[:5])
    a = np.array(rows, float)
    return box, a[:, 0].astype(int), a[:, 2:5]


def gate(perfect, inserted, cut, z_glide, x_line, seam_cut=3.0):
    bp, ip, rp = load_lammps_data(perfect)
    bi, ii, ri = load_lammps_data(inserted)
    bc, ic, rc = load_lammps_data(cut)
    L0 = bp['x'][1] - bp['x'][0]
    Ly = bp['y'][1] - bp['y'][0]
    op, oi = np.argsort(ip), np.argsort(ii)
    assert (ip[op] == ii[oi]).all(), 'atom ids of PERFECT and INSERTED differ'
    rp, ri = rp[op], ri[oi]
    u, x0, z0 = ri - rp, rp[:, 0], rp[:, 2]
    margin = 0.25 * L0
    dU = []
    for half in (z0 > z_glide + 10, z0 < z_glide - 10):
        lo = half & (x0 < x_line - margin)
        hi = half & (x0 > x_line + margin)
        dU.append(u[lo].mean(0) - u[hi].mean(0))
    xy_req = 0.5 * (dU[0][1] + dU[1][1])

    rho_t = (rp[:, 2] > z_glide).sum() / L0
    rho_b = (rp[:, 2] < z_glide).sum() / L0
    L_top = (rc[:, 2] > z_glide).sum() / rho_t
    L_bot = (rc[:, 2] < z_glide).sum() / rho_b
    Lx_req = 0.5 * (L_top + L_bot)
    Lx_now = bc['x'][1] - bc['x'][0]

    # seam coordination after the affine remap x -> x*Lx'/Lx, y -> y + xy'*x/Lx'
    x = (rc[:, 0] - bc['x'][0]) * Lx_req / Lx_now
    y = rc[:, 1] + xy_req * x / Lx_req
    x = np.mod(x, Lx_req)
    s = (x < 2 * seam_cut) | (x > Lx_req - 2 * seam_cut)
    q = np.c_[x[s], y[s], rc[s, 2]]
    zlo, zhi = rc[:, 2].min(), rc[:, 2].max()
    inner = (q[:, 2] > zlo + 6) & (q[:, 2] < zhi - 6) & ((q[:, 0] < seam_cut) | (q[:, 0] > Lx_req - seam_cut))
    zb = np.floor(q[:, 2] / 3.5).astype(int)
    hist = {}
    for b in np.unique(zb[inner]):
        A = q[(zb >= b - 1) & (zb <= b + 1)]
        B = q[inner & (zb == b)]
        d = B[:, None, :] - A[None, :, :]
        img = np.round(d[..., 0] / Lx_req)
        d[..., 0] -= Lx_req * img
        d[..., 1] -= xy_req * img
        d[..., 1] -= Ly * np.round(d[..., 1] / Ly)
        r = np.sqrt((d ** 2).sum(-1))
        n = ((r > 1e-6) & (r < 3.0)).sum(1)
        for v in n:
            hist[int(v)] = hist.get(int(v), 0) + 1
    ok = (set(hist) == {12})
    text = (f'atoms: perfect {len(ip)}  inserted {len(ii)}  after cut {len(ic)}  (deleted {len(ii) - len(ic)})\n'
            f'rigid-shift difference, low-x minus high-x end: top {np.round(dU[0], 4)}  bottom {np.round(dU[1], 4)}\n'
            f'natural x length: top {L_top:.4f}  bottom {L_bot:.4f}  -> edge content L_bot - L_top = {L_bot - L_top:+.4f} A\n'
            f'REQUIRED Lx = {Lx_req:.4f} A   (box as built {Lx_now:.4f}; difference {Lx_now - Lx_req:+.4f} A = {100 * (Lx_now / Lx_req - 1):+.3f} %)\n'
            f'REQUIRED delta = {xy_req:+.4f} A  (y-offset of the x-period vector; box as built 0)\n'
            f'seam coordination after remap (neighbours < 3.0 A : atoms): {dict(sorted(hist.items()))} -> {"PASS" if ok else "FAIL"}\n')
    return dict(Lx_req=Lx_req, xy_req=xy_req, Lx_now=Lx_now, pass_=ok, **{'pass': ok}, text=text)


if __name__ == '__main__':
    p, i, c, zg, xl = sys.argv[1:6]
    res = gate(p, i, c, float(zg), float(xl))
    print(res['text'], end='')
    sys.exit(0 if res['pass'] else 1)
