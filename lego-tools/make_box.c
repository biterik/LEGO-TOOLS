/* lego-make-box — ensure all atoms are inside the simulation box.
 *
 * Strategy:
 *   1. For periodic axes, wrap atoms back via PBC (only LAMMPS, orthogonal).
 *      For non-periodic axes (or non-LAMMPS inputs), this step is skipped.
 *   2. Compute the per-axis (min,max) of the resulting atom coordinates and,
 *      if any atom still lies outside the box, expand the box bounds so that
 *      every atom is contained. The box is never shrunk.
 *
 * Usage:
 *   lego-make-box [options] INPUT OUTPUT
 *
 * Options:
 *   --out-format=FMT   lammps|lammps-dump|imd  (default: lammps)
 *   --no-wrap          skip the PBC wrapping step
 *   --pad=EPS          extra padding added on each expanded side (default: 0)
 *   -h, --help
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "atomio.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-make-box [options] INPUT OUTPUT\n"
        "\n"
        "Wraps atoms with PBC (LAMMPS orthogonal only) and then expands the\n"
        "box bounds so that every atom lies inside.\n"
        "\n"
        "Options:\n"
        "  --out-format=FMT   lammps|lammps-dump|imd  (default: lammps)\n"
        "  --no-wrap          skip PBC wrapping step\n"
        "  --pad=EPS          extra padding on each expanded side (default 0)\n"
        "  -h, --help\n");
}

int main(int argc, char **argv) {
    const char *out_fmt_s = "lammps";
    int do_wrap = 1;
    double pad = 0.0;
    const char *pos[8];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage();
            return 0;
        } else if (!strncmp(a, "--out-format=", 13)) {
            out_fmt_s = a + 13;
        } else if (!strcmp(a, "--no-wrap")) {
            do_wrap = 0;
        } else if (!strncmp(a, "--pad=", 6)) {
            pad = strtod(a + 6, NULL);
        } else if (!strcmp(a, "--")) {
            /* end of options */
        } else if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-make-box: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 8) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if (npos != 2) { usage(); return 1; }
    const char *input  = pos[0];
    const char *output = pos[1];

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-make-box: unknown output format '%s'\n", out_fmt_s);
        return 1;
    }

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-make-box: input has no x/y/z columns\n");
        config_free(&c);
        return 1;
    }
    if (c.natoms == 0) {
        fprintf(stderr, "lego-make-box: no atoms in input\n");
        config_free(&c);
        return 1;
    }

    size_t N = c.natoms;
    int nc = c.ncols;
    int ix = c.idx_x, iy = c.idx_y, iz = c.idx_z;
    int iix = c.idx_ix, iiy = c.idx_iy, iiz = c.idx_iz;
    double *D = c.data;

    /* Step 1: PBC wrap (LAMMPS orthogonal only). */
    int wrap_ok = do_wrap
                && (c.format == FMT_LAMMPS_DATA || c.format == FMT_LAMMPS_DUMP)
                && !c.triclinic;
    int wrap_x = wrap_ok && c.pbc[0] && (c.xhi - c.xlo) > 0.0;
    int wrap_y = wrap_ok && c.pbc[1] && (c.yhi - c.ylo) > 0.0;
    int wrap_z = wrap_ok && c.pbc[2] && (c.zhi - c.zlo) > 0.0;

    size_t wrapped = 0;
    if (wrap_x || wrap_y || wrap_z) {
        double Lx = c.xhi - c.xlo;
        double Ly = c.yhi - c.ylo;
        double Lz = c.zhi - c.zlo;
        double xlo = c.xlo, ylo = c.ylo, zlo = c.zlo;

#pragma omp parallel for reduction(+:wrapped)
        for (size_t i = 0; i < N; i++) {
            double *row = D + i * (size_t)nc;
            int touched = 0;
            if (wrap_x) {
                double k = floor((row[ix] - xlo) / Lx);
                if (k != 0.0) {
                    row[ix] -= k * Lx;
                    if (iix >= 0) row[iix] += k;
                    touched = 1;
                }
            }
            if (wrap_y) {
                double k = floor((row[iy] - ylo) / Ly);
                if (k != 0.0) {
                    row[iy] -= k * Ly;
                    if (iiy >= 0) row[iiy] += k;
                    touched = 1;
                }
            }
            if (wrap_z) {
                double k = floor((row[iz] - zlo) / Lz);
                if (k != 0.0) {
                    row[iz] -= k * Lz;
                    if (iiz >= 0) row[iiz] += k;
                    touched = 1;
                }
            }
            if (touched) wrapped++;
        }
    }

    /* Step 2: parallel min/max reduction over coordinates. */
    double mnx =  DBL_MAX, mny =  DBL_MAX, mnz =  DBL_MAX;
    double mxx = -DBL_MAX, mxy = -DBL_MAX, mxz = -DBL_MAX;

#pragma omp parallel
    {
        double lmnx =  DBL_MAX, lmny =  DBL_MAX, lmnz =  DBL_MAX;
        double lmxx = -DBL_MAX, lmxy = -DBL_MAX, lmxz = -DBL_MAX;
#pragma omp for nowait
        for (size_t i = 0; i < N; i++) {
            const double *row = D + i * (size_t)nc;
            double x = row[ix], y = row[iy], z = row[iz];
            if (x < lmnx) lmnx = x;
            if (y < lmny) lmny = y;
            if (z < lmnz) lmnz = z;
            if (x > lmxx) lmxx = x;
            if (y > lmxy) lmxy = y;
            if (z > lmxz) lmxz = z;
        }
#pragma omp critical
        {
            if (lmnx < mnx) mnx = lmnx;
            if (lmny < mny) mny = lmny;
            if (lmnz < mnz) mnz = lmnz;
            if (lmxx > mxx) mxx = lmxx;
            if (lmxy > mxy) mxy = lmxy;
            if (lmxz > mxz) mxz = lmxz;
        }
    }

    /* Step 3: expand box only where necessary (never shrink). */
    double oxlo = c.xlo, oxhi = c.xhi;
    double oylo = c.ylo, oyhi = c.yhi;
    double ozlo = c.zlo, ozhi = c.zhi;
    int grew = 0;

    if (mnx < c.xlo) { c.xlo = mnx - pad; grew = 1; }
    if (mxx > c.xhi) { c.xhi = mxx + pad; grew = 1; }
    if (mny < c.ylo) { c.ylo = mny - pad; grew = 1; }
    if (mxy > c.yhi) { c.yhi = mxy + pad; grew = 1; }
    if (mnz < c.zlo) { c.zlo = mnz - pad; grew = 1; }
    if (mxz > c.zhi) { c.zhi = mxz + pad; grew = 1; }

    fprintf(stderr,
        "lego-make-box: wrapped %zu / %zu atoms; box %s\n"
        "  x: [%.6g %.6g] -> [%.6g %.6g]\n"
        "  y: [%.6g %.6g] -> [%.6g %.6g]\n"
        "  z: [%.6g %.6g] -> [%.6g %.6g]\n",
        wrapped, N, grew ? "grew" : "unchanged",
        oxlo, oxhi, c.xlo, c.xhi,
        oylo, oyhi, c.ylo, c.yhi,
        ozlo, ozhi, c.zlo, c.zhi);

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    return rc;
}
