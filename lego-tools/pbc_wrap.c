/* lego-pbc-wrap — wrap atoms back into the simulation box via PBC.
 *
 * LAMMPS only (data or dump). Triclinic and IMD are rejected.
 *
 * If image-flag columns (ix/iy/iz) are present they are updated so that
 * the unwrapped position x_real = x + ix*Lx (etc.) is preserved.
 *
 * Usage:
 *   lego-pbc-wrap [options] INPUT OUTPUT
 *
 * Options:
 *   --out-format=FMT   lammps|lammps-dump  (default: lammps)
 *   -x, -y, -z         wrap only this axis  (default: all three)
 *   -h, --help
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "atomio.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-pbc-wrap [options] INPUT OUTPUT\n"
        "\n"
        "Wraps atoms outside the simulation box back inside using PBC.\n"
        "Works on LAMMPS data and dump files only (orthogonal boxes).\n"
        "Image flags (ix/iy/iz) are updated if present.\n"
        "\n"
        "Options:\n"
        "  --out-format=FMT   lammps|lammps-dump  (default: lammps)\n"
        "  -x                 wrap only x axis\n"
        "  -y                 wrap only y axis\n"
        "  -z                 wrap only z axis\n"
        "                     (default: all three axes)\n"
        "  -h, --help\n");
}

int main(int argc, char **argv) {
    const char *out_fmt_s = "lammps";
    int do_x = 0, do_y = 0, do_z = 0;
    const char *pos[8];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage();
            return 0;
        } else if (!strncmp(a, "--out-format=", 13)) {
            out_fmt_s = a + 13;
        } else if (!strcmp(a, "-x")) {
            do_x = 1;
        } else if (!strcmp(a, "-y")) {
            do_y = 1;
        } else if (!strcmp(a, "-z")) {
            do_z = 1;
        } else if (!strcmp(a, "--")) {
            /* end of options */
        } else if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-pbc-wrap: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 8) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if (npos != 2) { usage(); return 1; }
    if (!do_x && !do_y && !do_z) { do_x = do_y = do_z = 1; }

    const char *input  = pos[0];
    const char *output = pos[1];

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt != FMT_LAMMPS_DATA && out_fmt != FMT_LAMMPS_DUMP) {
        fprintf(stderr, "lego-pbc-wrap: output format must be lammps "
                        "or lammps-dump\n");
        return 1;
    }

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.format != FMT_LAMMPS_DATA && c.format != FMT_LAMMPS_DUMP) {
        fprintf(stderr, "lego-pbc-wrap: input is not a LAMMPS file "
                        "(detected %s)\n", format_name(c.format));
        config_free(&c);
        return 1;
    }
    if (c.triclinic) {
        fprintf(stderr, "lego-pbc-wrap: triclinic boxes are not supported\n");
        config_free(&c);
        return 1;
    }
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-pbc-wrap: input has no x/y/z columns\n");
        config_free(&c);
        return 1;
    }

    double Lx = c.xhi - c.xlo;
    double Ly = c.yhi - c.ylo;
    double Lz = c.zhi - c.zlo;
    if (do_x && Lx <= 0.0) {
        fprintf(stderr, "lego-pbc-wrap: zero-length x box\n");
        config_free(&c); return 1;
    }
    if (do_y && Ly <= 0.0) {
        fprintf(stderr, "lego-pbc-wrap: zero-length y box\n");
        config_free(&c); return 1;
    }
    if (do_z && Lz <= 0.0) {
        fprintf(stderr, "lego-pbc-wrap: zero-length z box\n");
        config_free(&c); return 1;
    }

    size_t N = c.natoms;
    int nc = c.ncols;
    int ix = c.idx_x, iy = c.idx_y, iz = c.idx_z;
    int iix = c.idx_ix, iiy = c.idx_iy, iiz = c.idx_iz;
    double xlo = c.xlo, ylo = c.ylo, zlo = c.zlo;
    double *D = c.data;

    size_t wrapped = 0;

#pragma omp parallel for reduction(+:wrapped)
    for (size_t i = 0; i < N; i++) {
        double *row = D + i * (size_t)nc;
        int touched = 0;

        if (do_x) {
            double x = row[ix];
            double s = (x - xlo) / Lx;
            double k = floor(s);
            if (k != 0.0) {
                row[ix] = x - k * Lx;
                if (iix >= 0) row[iix] += k;
                touched = 1;
            }
        }
        if (do_y) {
            double y = row[iy];
            double s = (y - ylo) / Ly;
            double k = floor(s);
            if (k != 0.0) {
                row[iy] = y - k * Ly;
                if (iiy >= 0) row[iiy] += k;
                touched = 1;
            }
        }
        if (do_z) {
            double z = row[iz];
            double s = (z - zlo) / Lz;
            double k = floor(s);
            if (k != 0.0) {
                row[iz] = z - k * Lz;
                if (iiz >= 0) row[iiz] += k;
                touched = 1;
            }
        }
        if (touched) wrapped++;
    }

    fprintf(stderr, "lego-pbc-wrap: wrapped %zu / %zu atoms\n", wrapped, N);

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    return rc;
}
