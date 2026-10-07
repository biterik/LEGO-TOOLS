#!/usr/bin/env python3
"""Run relax_single_dislocation_pps.lmp through the LAMMPS Python module (or set $LMP to a binary).
usage: run_lammps_relax.py NAME CELLDIR OUTDIR POTFILE ELEMENT CORR(0|1)"""
import os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
name, celldir, outdir, pot, elem, corr = sys.argv[1:7]
os.makedirs(outdir, exist_ok=True)
variant = 'corrected' if corr == '1' else 'uncorrected'
log = os.path.join(outdir, f'log_relax_{name}_{variant}.lammps')
args = ['-var', 'name', name, '-var', 'celldir', celldir, '-var', 'outdir', outdir,
        '-var', 'pot', pot, '-var', 'elem', elem, '-var', 'corr', corr, '-log', log, '-screen', 'none',
        '-in', os.path.join(here, 'relax_single_dislocation_pps.lmp')]
if os.environ.get('LMP'):
    subprocess.run([os.environ['LMP']] + args, check=True)
else:
    from lammps import lammps
    l = lammps(cmdargs=args[:-2])
    l.file(args[-1])
    l.close()
print(open(log).read().split('RESULT')[-1].split('\n')[0] if 'RESULT' in open(log).read() else 'NO RESULT')
