"""python -m legodislo — command-line equivalent of the legodislo API.

Runs the compiled lego-dislo binary (auto-located; override with
--binary or $LEGO_DISLO) on an existing card file:

  python -m legodislo CARD INPUT OUTPUT [--out-format=FMT] [--threads=N]
                                        [--report=FILE] [--binary=PATH]
  python -m legodislo --solve-only CARD [--binary=PATH]

License: GPL-3.0
Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
         Max-Planck-Institut fuer Nachhaltige Materialien, Duesseldorf
Funding: NFDI-MatWerk
"""

import argparse
import subprocess
import sys

from . import LegoDisloError, find_binary


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m legodislo",
        description="Insert dislocations (anisotropic elastic field) into "
                    "an atomistic configuration via the lego-dislo binary.")
    ap.add_argument("card", help="input card file (see dislo/README.md)")
    ap.add_argument("input", nargs="?",
                    help="input configuration (omit with --solve-only)")
    ap.add_argument("output", nargs="?",
                    help="output configuration (omit with --solve-only)")
    ap.add_argument("--solve-only", action="store_true",
                    help="solve and report only; no atoms")
    ap.add_argument("--out-format", choices=["lammps", "lammps-dump", "imd"],
                    help="output format (default: lammps)")
    ap.add_argument("--threads", type=int, help="number of OpenMP threads")
    ap.add_argument("--report", metavar="FILE",
                    help="write the solver report to FILE")
    ap.add_argument("--binary", metavar="PATH",
                    help="path to the lego-dislo binary")
    args = ap.parse_args(argv)

    if args.solve_only:
        if args.input or args.output:
            ap.error("--solve-only takes only CARD")
    elif not (args.input and args.output):
        ap.error("INPUT and OUTPUT are required (or use --solve-only)")

    try:
        exe = find_binary(args.binary)
    except LegoDisloError as e:
        print(f"legodislo: {e}", file=sys.stderr)
        return 1

    cmd = [exe]
    if args.solve_only:
        cmd.append("--solve-only")
    if args.out_format:
        cmd.append(f"--out-format={args.out_format}")
    if args.threads:
        cmd.append(f"--threads={args.threads}")
    if args.report:
        cmd.append(f"--report={args.report}")
    cmd.append(args.card)
    if not args.solve_only:
        cmd += [args.input, args.output]
    return subprocess.call(cmd)


if __name__ == "__main__":
    sys.exit(main())
