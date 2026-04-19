/* lego-cut — keep atoms whose x,y,z lie inside a [min,max]^3 box.
 *
 * Usage:
 *   lego-cut [options] xmin xmax ymin ymax zmin zmax INPUT OUTPUT
 *
 *   --in-format=lammps|lammps-dump|imd     (default: auto-detect)
 *   --out-format=lammps|lammps-dump|imd    (default: lammps)
 *
 * Box bounds are preserved as in the input (per user choice).
 * Atoms outside the range are dropped and the atom count updated.
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
#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-cut [options] xmin xmax ymin ymax zmin zmax INPUT OUTPUT\n"
        "\n"
        "Keep atoms with xmin<=x<=xmax, ymin<=y<=ymax, zmin<=z<=zmax.\n"
        "Box bounds are NOT modified.\n"
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
        } else if (a[0] == '-' && a[1] == '-' && a[2] == 0) {
            /* -- ends options */
        } else if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-cut: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 16) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if (npos != 8) { usage(); return 1; }

    double xmin = strtod(pos[0], NULL);
    double xmax = strtod(pos[1], NULL);
    double ymin = strtod(pos[2], NULL);
    double ymax = strtod(pos[3], NULL);
    double zmin = strtod(pos[4], NULL);
    double zmax = strtod(pos[5], NULL);
    const char *input  = pos[6];
    const char *output = pos[7];

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-cut: unknown output format '%s'\n", out_fmt_s);
        return 1;
    }

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-cut: input has no x/y/z columns\n");
        return 1;
    }

    size_t N = c.natoms;
    int nc = c.ncols;
    unsigned char *keep = (unsigned char *)malloc(N);
    if (!keep) { fprintf(stderr, "lego-cut: oom\n"); return 1; }

    size_t nkeep = 0;
#pragma omp parallel for reduction(+:nkeep)
    for (size_t i = 0; i < N; i++) {
        double x = c.data[i * (size_t)nc + c.idx_x];
        double y = c.data[i * (size_t)nc + c.idx_y];
        double z = c.data[i * (size_t)nc + c.idx_z];
        int k = (x >= xmin && x <= xmax &&
                 y >= ymin && y <= ymax &&
                 z >= zmin && z <= zmax);
        keep[i] = (unsigned char)k;
        if (k) nkeep++;
    }

    /* Compact in place (sequential — writes are cheap vs the compare). */
    double *D = c.data;
    size_t w = 0;
    for (size_t i = 0; i < N; i++) {
        if (!keep[i]) continue;
        if (w != i) {
            memcpy(D + w * (size_t)nc, D + i * (size_t)nc,
                   (size_t)nc * sizeof(double));
        }
        w++;
    }
    free(keep);
    c.natoms = nkeep;

    fprintf(stderr, "lego-cut: kept %zu / %zu atoms\n", nkeep, N);

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    return rc;
}
