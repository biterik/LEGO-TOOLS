#!/usr/bin/env python3
"""Build one straight a/2<110>{111} dislocation in an fcc slab that is periodic along the glide
direction x and the line y, with free (111) surfaces in z (LAMMPS 'p p s').

Chain:  LEGO (perfect slab) -> dcreator (upper_half rigid shift) -> lego-cut (drop the slab the
rigid shift pushed out of the box) -> gate (box correction + seam check).

The gate computes the corrected x-period vector a = (Lx', delta, 0):
    Lx'   = (L_top + L_bot)/2  (= L0 - b_e/2)      delta = <dU_y> over the two halves (= +-b_s/2)
and writes it to box-correction_<tag>.lmp. apply_period_vector.py then imposes it by ONE affine map
and writes a triclinic LAMMPS cell with the LINE along LAMMPS x (needed because LAMMPS can only tilt
its second box vector; 'change_box xy' would shear along the wrong axis).

usage: build_single_dislocation_cell.py CHARACTER OUTDIR [--a0 3.52] [--Lx 150 --Ly 30 --Lz 150]
       CHARACTER: screw | 30 | 60 | edge
Binaries: $LEGO, $DCREATOR, $LEGO_CUT (default: found on $PATH as lego, dcreator, lego-cut).
"""
import argparse, os, shutil, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gate_single_dislocation_pps_cell import gate, load_lammps_data   # noqa: E402

# crystal frame, slab orientation and the dislocation; b in units of a0/2
ORIENT_A = dict(new_x=(-2, 1, 1), new_z=(1, 1, 1))     # x=[-211] y=[0-11] z=[111]
ORIENT_B = dict(new_x=(-1, 0, 1), new_z=(1, 1, 1))     # x=[-101] y=[1-21] z=[111]
CHARACTERS = {
    'screw': dict(orient=ORIENT_A, tag='x-211_y0-11_z111', b=(0, -1, 1), line=(0, -1, 1)),
    '60':    dict(orient=ORIENT_A, tag='x-211_y0-11_z111', b=(-1, 1, 0), line=(0, -1, 1)),
    '30':    dict(orient=ORIENT_B, tag='x-101_y1-21_z111', b=(0, -1, 1), line=(1, -2, 1)),
    'edge':  dict(orient=ORIENT_B, tag='x-101_y1-21_z111', b=(-1, 0, 1), line=(1, -2, 1)),
}
GLIDE_NORMAL = (1, 1, 1)


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def rotation(new_x, new_z):
    x, z = unit(new_x), unit(new_z)
    return np.array([x, np.cross(z, x), z])          # rows: sample axes in crystal coords


def run(cmd, log):
    with open(log, 'a') as fh:
        fh.write('$ ' + ' '.join(cmd) + '\n'); fh.flush()
        subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('character', choices=CHARACTERS)
    ap.add_argument('outdir')
    ap.add_argument('--a0', type=float, default=3.52)
    ap.add_argument('--Lx', type=float, default=150.0)
    ap.add_argument('--Ly', type=float, default=30.0)
    ap.add_argument('--Lz', type=float, default=150.0)
    ap.add_argument('--stripe', type=float, default=30.0)
    ap.add_argument('--element', default='Ni')
    ap.add_argument('--mass', type=float, default=58.6934)
    a = ap.parse_args()
    lego = os.environ.get('LEGO', shutil.which('lego') or 'lego')
    dcreator = os.environ.get('DCREATOR', shutil.which('dcreator') or 'dcreator')
    lego_cut = os.environ.get('LEGO_CUT', shutil.which('lego-cut') or 'lego-cut')

    c = CHARACTERS[a.character]
    os.makedirs(a.outdir, exist_ok=True)
    name = f'{a.element}_{a.character}{"deg" if a.character.isdigit() else ""}_{c["tag"]}'
    log = os.path.join(a.outdir, f'build-log_{name}.txt')
    open(log, 'w').close()
    f_perf = os.path.join(a.outdir, f'{name}_perfect-slab.lmp')
    f_ins = os.path.join(a.outdir, f'{name}_dcreator-inserted.lmp')
    f_cut = os.path.join(a.outdir, f'{name}_cut.lmp')

    # 1. LEGO perfect slab
    p_lego = os.path.join(a.outdir, f'lego_{name}.param')
    with open(p_lego, 'w') as fh:
        fh.write(f"""structure      fcc
lattice_const  {a.a0}
mass           {a.mass}
new_x          {' '.join(map(str, c['orient']['new_x']))}
new_z          {' '.join(map(str, c['orient']['new_z']))}
box_x          {a.Lx}
box_y          {a.Ly}
box_z          {a.Lz}
box_mode       periodic
outfile        {f_perf}
output_format  lammps
""")
    run([lego, p_lego], log)
    box, ids, r = load_lammps_data(f_perf)
    L0 = box['x'][1] - box['x'][0]

    # 2. glide plane: the (111) interplanar gap closest to mid-height; dcreator point exactly mid-gap
    zs = np.unique(np.round(r[:, 2], 4))
    gaps = 0.5 * (zs[1:] + zs[:-1])
    zg = gaps[np.argmin(np.abs(gaps - 0.5 * (zs[0] + zs[-1])))]
    x_line = 0.5 * L0

    # 3. rigid side: the edge part of the rigid shift must point OUT of the box through the face the
    #    rigid region reaches (overlap -> delete). Pointing in opens a b_e gap = opposite dislocation.
    R = rotation(c['orient']['new_x'], c['orient']['new_z'])
    b_crys = 0.5 * a.a0 * np.asarray(c['b'], float)
    b_s = R @ b_crys
    l_s, n_s = R @ unit(c['line']), R @ unit(GLIDE_NORMAL)
    x_d = np.cross(n_s, l_s)                                   # dcreator: x_d = y_d x z_d
    if abs(b_s[0]) < 1e-8:
        side = 'right'
    else:
        side = 'right' if b_s[0] * x_d[0] > 0 else 'left'

    p_dc = os.path.join(a.outdir, f'dcreator_{name}.dissparam')
    with open(p_dc, 'w') as fh:
        fh.write(f"""configfile         {f_perf}
outfile            {f_ins}
new_x              {' '.join(map(str, c['orient']['new_x']))}
new_z              {' '.join(map(str, c['orient']['new_z']))}
burgersv           {' '.join(f'{v:.6f}' for v in b_crys)}
linev              {' '.join(map(str, c['line']))}
glideplane         {' '.join(map(str, GLIDE_NORMAL))}
point              {x_line:.6f} 0.0 {zg:.6f}
stripe             {a.stripe}
displacement_mode  upper_half
rigid_side         {side}
input_format       lammps_data
output_format      lammps_data
""")
    run([dcreator, p_dc], log)

    # 4. delete exactly the atoms the rigid shift pushed out of [0, L0) in x
    run([lego_cut, '0', f'{L0 - 1e-9:.10f}', '-1e9', '1e9', '-1e9', '1e9', f_ins, f_cut], log)

    # 5. gate: required period vector + seam coordination after the remap
    res = gate(f_perf, f_ins, f_cut, zg, x_line)
    b_e, b_scr = abs(b_s[0]), abs(b_s[1])
    rep = os.path.join(a.outdir, f'gate-report_{name}.txt')
    with open(rep, 'w') as fh:
        fh.write(f'character {a.character}  b(sample) = {np.round(b_s, 4)}  |b_e| = {b_e:.4f}  |b_s| = {b_scr:.4f}\n')
        fh.write(f'dcreator rigid_side = {side}  (x_d = {np.round(x_d, 3)})  glide gap z = {zg:.4f}  line x = {x_line:.4f}\n')
        fh.write(res['text'])
    print(open(rep).read())
    if not res['pass']:
        sys.exit(f'GATE FAILED for {name} -- see {rep}')

    zmin, zmax = r[:, 2].min(), r[:, 2].max()
    d111 = a.a0 / np.sqrt(3.0)
    with open(os.path.join(a.outdir, f'box-correction_{name}.lmp'), 'w') as fh:
        fh.write(f"""# written by build_single_dislocation_cell.py for {name}
variable Lx_built  equal {res['Lx_now']:.10f}
variable Lx_corr   equal {res['Lx_req']:.10f}   # = (L_top + L_bot)/2
variable xy_corr   equal {res['xy_req']:.10f}   # = delta = <dU_y> over the halves: y-offset of the x-period vector
variable z_glide   equal {zg:.10f}
variable x_line    equal {x_line:.10f}
variable z_fix_lo  equal {zmin + 3 * d111 - 0.3:.6f}   # 3 (111) layers per face held in z
variable z_fix_hi  equal {zmax - 3 * d111 + 0.3:.6f}
variable b_edge    equal {b_e:.10f}
variable b_screw   equal {b_scr:.10f}
""")

    subprocess.run([sys.executable, os.path.join(HERE, 'apply_period_vector.py'), f_cut,
                    os.path.join(a.outdir, f'box-correction_{name}.lmp'),
                    os.path.join(a.outdir, f'{name}_corrected-period-vector.lmp')], check=True)
    # independent check of the file LAMMPS will read (glide direction is LAMMPS y there)
    subprocess.run([sys.executable, os.path.join(HERE, 'check_seam_coordination.py'),
                    os.path.join(a.outdir, f'{name}_corrected-period-vector.lmp'), 'y'], check=True)


if __name__ == '__main__':
    main()
