/* lego-analyze — per-column min/max, atom counts by type, box info,
 *                velocity presence check, and format detection.
 *
 * Usage:
 *   lego-analyze [--json] INPUT
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "atomio.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-analyze [--json] INPUT\n"
        "\n"
        "Reports format, atom count, types, box, per-column min/max,\n"
        "and whether velocities are present.\n");
}

typedef struct {
    double vmin, vmax;
} Range;

static void compute_ranges(const Config *c, Range *r) {
    size_t N = c->natoms;
    int nc = c->ncols;
    const double *D = c->data;

    for (int k = 0; k < nc; k++) {
        r[k].vmin =  DBL_MAX;
        r[k].vmax = -DBL_MAX;
    }
    if (N == 0) return;

#ifdef _OPENMP
    int nth = omp_get_max_threads();
    if (nth < 1) nth = 1;
    Range *tl = (Range *)malloc((size_t)nth * (size_t)nc * sizeof(Range));
    for (int t = 0; t < nth; t++)
        for (int k = 0; k < nc; k++) {
            tl[t * nc + k].vmin =  DBL_MAX;
            tl[t * nc + k].vmax = -DBL_MAX;
        }
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        Range *loc = &tl[tid * nc];
#pragma omp for
        for (size_t i = 0; i < N; i++) {
            const double *row = D + i * (size_t)nc;
            for (int k = 0; k < nc; k++) {
                double v = row[k];
                if (v < loc[k].vmin) loc[k].vmin = v;
                if (v > loc[k].vmax) loc[k].vmax = v;
            }
        }
    }
    for (int t = 0; t < nth; t++)
        for (int k = 0; k < nc; k++) {
            if (tl[t * nc + k].vmin < r[k].vmin) r[k].vmin = tl[t * nc + k].vmin;
            if (tl[t * nc + k].vmax > r[k].vmax) r[k].vmax = tl[t * nc + k].vmax;
        }
    free(tl);
#else
    for (size_t i = 0; i < N; i++) {
        const double *row = D + i * (size_t)nc;
        for (int k = 0; k < nc; k++) {
            double v = row[k];
            if (v < r[k].vmin) r[k].vmin = v;
            if (v > r[k].vmax) r[k].vmax = v;
        }
    }
#endif
}

static void count_by_type(const Config *c, size_t **counts_out, int *ntypes_out) {
    *counts_out = NULL;
    *ntypes_out = 0;
    if (c->idx_type < 0) return;
    int nc = c->ncols;
    size_t N = c->natoms;

    /* First pass: find max type. */
    int maxt = 0;
    for (size_t i = 0; i < N; i++) {
        int t = (int)c->data[i * (size_t)nc + c->idx_type];
        if (t > maxt) maxt = t;
    }
    if (maxt <= 0) return;
    size_t *cnt = (size_t *)calloc((size_t)maxt + 1, sizeof(size_t));
    for (size_t i = 0; i < N; i++) {
        int t = (int)c->data[i * (size_t)nc + c->idx_type];
        if (t >= 0 && t <= maxt) cnt[t]++;
    }
    *counts_out = cnt;
    *ntypes_out = maxt;
}

/* Rough validity checks. */
static int validate(const Config *c, char *reason, size_t reason_sz) {
    if (c->format == FMT_LAMMPS_DATA || c->format == FMT_LAMMPS_DUMP) {
        if (c->idx_id < 0 || c->idx_type < 0 ||
            c->idx_x < 0 || c->idx_y < 0 || c->idx_z < 0) {
            snprintf(reason, reason_sz,
                     "LAMMPS requires id, type, x, y, z columns");
            return 0;
        }
    }
    if (c->format == FMT_IMD) {
        if (c->idx_x < 0 || c->idx_y < 0 || c->idx_z < 0) {
            snprintf(reason, reason_sz, "IMD requires x, y, z columns");
            return 0;
        }
    }
    if (c->xhi < c->xlo || c->yhi < c->ylo || c->zhi < c->zlo) {
        snprintf(reason, reason_sz, "box bounds inverted");
        return 0;
    }
    reason[0] = 0;
    return 1;
}

int main(int argc, char **argv) {
    int json = 0;
    int ai = 1;
    while (ai < argc && argv[ai][0] == '-') {
        if (!strcmp(argv[ai], "--json")) { json = 1; ai++; }
        else if (!strcmp(argv[ai], "-h") || !strcmp(argv[ai], "--help")) {
            usage();
            return 0;
        } else {
            fprintf(stderr, "lego-analyze: unknown option %s\n", argv[ai]);
            return 1;
        }
    }
    if (argc - ai != 1) { usage(); return 1; }
    const char *input = argv[ai];

    Config c;
    if (config_read(&c, input) != 0) return 1;

    Range *ranges = (Range *)calloc((size_t)c.ncols, sizeof(Range));
    compute_ranges(&c, ranges);

    size_t *type_counts = NULL;
    int max_type = 0;
    count_by_type(&c, &type_counts, &max_type);

    char validity[256];
    int valid = validate(&c, validity, sizeof(validity));

    if (!json) {
        printf("file:           %s\n", input);
        printf("format:         %s\n", format_name(c.format));
        if (c.format == FMT_LAMMPS_DATA && c.atom_style[0])
            printf("atom_style:     %s\n", c.atom_style);
        printf("valid:          %s%s%s\n",
               valid ? "yes" : "no (",
               valid ? "" : validity,
               valid ? "" : ")");
        printf("natoms:         %zu\n", c.natoms);
        printf("velocities:     %s\n", c.has_velocities ? "yes" : "no");
        printf("triclinic:      %s\n", c.triclinic ? "yes" : "no");
        printf("box:\n");
        printf("  xlo xhi     %.16g %.16g   (Lx = %.16g)\n",
               c.xlo, c.xhi, c.xhi - c.xlo);
        printf("  ylo yhi     %.16g %.16g   (Ly = %.16g)\n",
               c.ylo, c.yhi, c.yhi - c.ylo);
        printf("  zlo zhi     %.16g %.16g   (Lz = %.16g)\n",
               c.zlo, c.zhi, c.zhi - c.zlo);
        if (c.triclinic)
            printf("  xy xz yz    %.16g %.16g %.16g\n", c.xy, c.xz, c.yz);

        if (type_counts) {
            printf("types (%d):\n", max_type);
            for (int t = 1; t <= max_type; t++)
                if (type_counts[t] > 0)
                    printf("  type %d:   %zu atoms\n", t, type_counts[t]);
        }

        printf("columns (%d):\n", c.ncols);
        printf("  %-16s %-24s %-24s\n", "name", "min", "max");
        for (int k = 0; k < c.ncols; k++) {
            if (c.col_is_int[k])
                printf("  %-16s %-24lld %-24lld\n",
                       c.col_names[k],
                       (long long)ranges[k].vmin,
                       (long long)ranges[k].vmax);
            else
                printf("  %-16s %-24.16g %-24.16g\n",
                       c.col_names[k], ranges[k].vmin, ranges[k].vmax);
        }
    } else {
        printf("{\n");
        printf("  \"file\": \"%s\",\n", input);
        printf("  \"format\": \"%s\",\n", format_name(c.format));
        printf("  \"atom_style\": \"%s\",\n", c.atom_style);
        printf("  \"valid\": %s,\n", valid ? "true" : "false");
        if (!valid) printf("  \"validity_reason\": \"%s\",\n", validity);
        printf("  \"natoms\": %zu,\n", c.natoms);
        printf("  \"velocities\": %s,\n", c.has_velocities ? "true" : "false");
        printf("  \"triclinic\": %s,\n", c.triclinic ? "true" : "false");
        printf("  \"box\": {\n");
        printf("    \"xlo\": %.16g, \"xhi\": %.16g,\n", c.xlo, c.xhi);
        printf("    \"ylo\": %.16g, \"yhi\": %.16g,\n", c.ylo, c.yhi);
        printf("    \"zlo\": %.16g, \"zhi\": %.16g", c.zlo, c.zhi);
        if (c.triclinic)
            printf(",\n    \"xy\": %.16g, \"xz\": %.16g, \"yz\": %.16g\n",
                   c.xy, c.xz, c.yz);
        else
            printf("\n");
        printf("  },\n");

        printf("  \"type_counts\": {");
        if (type_counts) {
            int first = 1;
            for (int t = 1; t <= max_type; t++) {
                if (type_counts[t] == 0) continue;
                printf("%s\"%d\": %zu", first ? "" : ", ", t, type_counts[t]);
                first = 0;
            }
        }
        printf("},\n");

        printf("  \"columns\": [\n");
        for (int k = 0; k < c.ncols; k++) {
            printf("    {\"name\": \"%s\", \"min\": %.16g, \"max\": %.16g, \"int\": %s}%s\n",
                   c.col_names[k], ranges[k].vmin, ranges[k].vmax,
                   c.col_is_int[k] ? "true" : "false",
                   k + 1 == c.ncols ? "" : ",");
        }
        printf("  ]\n");
        printf("}\n");
    }

    free(ranges);
    free(type_counts);
    config_free(&c);
    return 0;
}
