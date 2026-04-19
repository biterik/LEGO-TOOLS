/* lego-remove-per-atom — strip every per-atom column except
 * id, type, x, y, z (and drop a Velocities section, forces, charges, …).
 *
 * Usage:
 *   lego-remove-per-atom [--out-format=...] INPUT OUTPUT
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

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-remove-per-atom [options] INPUT OUTPUT\n"
        "\n"
        "Keeps only id, type, x, y, z columns. Drops velocities,\n"
        "forces, charges, energies, image flags, and any other\n"
        "per-atom quantities present.\n"
        "\n"
        "Options:\n"
        "  --out-format=FMT   lammps|lammps-dump|imd  (default: lammps)\n"
        "  -h, --help\n");
}

int main(int argc, char **argv) {
    const char *out_fmt_s = "lammps";

    int ai = 1;
    while (ai < argc && argv[ai][0] == '-') {
        if (!strcmp(argv[ai], "-h") || !strcmp(argv[ai], "--help")) {
            usage();
            return 0;
        } else if (!strncmp(argv[ai], "--out-format=", 13)) {
            out_fmt_s = argv[ai] + 13;
            ai++;
        } else {
            fprintf(stderr, "lego-remove-per-atom: unknown option %s\n", argv[ai]);
            return 1;
        }
    }
    if (argc - ai != 2) { usage(); return 1; }
    const char *input  = argv[ai++];
    const char *output = argv[ai++];

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-remove-per-atom: unknown output format '%s'\n",
                out_fmt_s);
        return 1;
    }

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-remove-per-atom: input has no x/y/z columns\n");
        return 1;
    }

    /* Build the stripped config by materialising a fresh, minimal layout. */
    int want_id   = c.idx_id   >= 0;
    int want_type = c.idx_type >= 0;
    int new_nc = (want_id ? 1 : 0) + (want_type ? 1 : 0) + 3;

    double *nd = (double *)malloc(c.natoms * (size_t)new_nc * sizeof(double));
    if (!nd) { fprintf(stderr, "lego-remove-per-atom: oom\n"); return 1; }

    int old_nc = c.ncols;
    /* Source column indices for each column in the new minimal layout.
     * Integer-ness of the new columns is fixed by the layout itself
     * (id/type int, x/y/z float) and is set below when we rewrite the
     * column metadata, so we don't need a parallel dst_is_int[] array. */
    int src[8];
    int ncopy = 0;
    if (want_id)   { src[ncopy++] = c.idx_id;   }
    if (want_type) { src[ncopy++] = c.idx_type; }
    src[ncopy++] = c.idx_x;
    src[ncopy++] = c.idx_y;
    src[ncopy++] = c.idx_z;

    size_t N = c.natoms;
#pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
        const double *s = c.data + i * (size_t)old_nc;
        double *d = nd + i * (size_t)new_nc;
        for (int k = 0; k < ncopy; k++) d[k] = s[src[k]];
    }

    /* Swap in the new layout. */
    free(c.data);
    c.data = nd;
    c.ncols = new_nc;

    /* Reset column metadata. */
    int k = 0;
    if (want_id)   { strcpy(c.col_names[k], "id");   c.col_is_int[k++] = 1; }
    if (want_type) { strcpy(c.col_names[k], "type"); c.col_is_int[k++] = 1; }
    strcpy(c.col_names[k], "x"); c.col_is_int[k++] = 0;
    strcpy(c.col_names[k], "y"); c.col_is_int[k++] = 0;
    strcpy(c.col_names[k], "z"); c.col_is_int[k++] = 0;
    for (int j = k; j < MAX_COLS; j++) {
        c.col_names[j][0] = 0;
        c.col_is_int[j] = 0;
    }
    c.has_velocities = 0;
    strcpy(c.atom_style, "atomic");
    config_refresh_std_indices(&c);

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    return rc;
}
