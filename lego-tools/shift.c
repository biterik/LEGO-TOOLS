/* lego-shift — shift all atom coordinates and the box bounds by (dx,dy,dz).
 *
 * Usage:
 *   lego-shift [--out-format=...] dx dy dz INPUT OUTPUT
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "atomio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-shift [options] dx dy dz INPUT OUTPUT\n"
        "\n"
        "Options:\n"
        "  --out-format=FMT   lammps|lammps-dump|imd  (default: lammps)\n"
        "  -h, --help\n");
}

int main(int argc, char **argv) {
    const char *out_fmt_s = "lammps";
    const char *pos[16];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage();
            return 0;
        } else if (!strncmp(a, "--out-format=", 13)) {
            out_fmt_s = a + 13;
        } else if (!strcmp(a, "--")) {
            /* -- ends options */
        } else if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-shift: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 16) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if (npos != 5) { usage(); return 1; }

    double dx = strtod(pos[0], NULL);
    double dy = strtod(pos[1], NULL);
    double dz = strtod(pos[2], NULL);
    const char *input  = pos[3];
    const char *output = pos[4];

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-shift: unknown output format '%s'\n", out_fmt_s);
        return 1;
    }

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-shift: input has no x/y/z columns\n");
        return 1;
    }

    c.xlo += dx; c.xhi += dx;
    c.ylo += dy; c.yhi += dy;
    c.zlo += dz; c.zhi += dz;

    size_t N = c.natoms;
    int nc = c.ncols;
    int ix = c.idx_x, iy = c.idx_y, iz = c.idx_z;
    double *D = c.data;

#pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
        D[i * (size_t)nc + ix] += dx;
        D[i * (size_t)nc + iy] += dy;
        D[i * (size_t)nc + iz] += dz;
    }

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    return rc;
}
