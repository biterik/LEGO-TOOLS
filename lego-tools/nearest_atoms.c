/* lego-nearest-atoms — find the k nearest atoms to a given point (x,y,z).
 *
 * Usage:
 *   lego-nearest-atoms [options] x y z INPUT
 *
 * Options:
 *   -k N            report N nearest atoms (default: 5)
 *   --tol=EPS       distance <= EPS counts as "sitting there" (default: 1e-6)
 *   --no-pbc        ignore periodic boundary conditions (default: use PBC
 *                   via minimum-image convention for orthogonal boxes)
 *   --json          machine-readable output
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
        "Usage: lego-nearest-atoms [options] x y z INPUT\n"
        "\n"
        "Finds the k nearest atoms to the point (x,y,z).\n"
        "\n"
        "Options:\n"
        "  -k N           number of neighbours to report (default: 5)\n"
        "  --tol=EPS      distance <= EPS counts as \"sitting there\"\n"
        "                 (default: 1e-6)\n"
        "  --no-pbc       ignore periodic boundary conditions\n"
        "  --json         machine-readable output\n"
        "  -h, --help\n");
}

typedef struct {
    double d2;       /* squared distance */
    size_t row;      /* row in c.data */
} Hit;

/* Max-heap of size up to k (largest d2 at index 0) so we can evict the
 * worst candidate in O(log k) when we find a better one. */
static void heap_up(Hit *h, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h[p].d2 >= h[i].d2) break;
        Hit t = h[p]; h[p] = h[i]; h[i] = t;
        i = p;
    }
}
static void heap_down(Hit *h, int n, int i) {
    for (;;) {
        int l = 2 * i + 1, r = 2 * i + 2, m = i;
        if (l < n && h[l].d2 > h[m].d2) m = l;
        if (r < n && h[r].d2 > h[m].d2) m = r;
        if (m == i) break;
        Hit t = h[m]; h[m] = h[i]; h[i] = t;
        i = m;
    }
}
static void heap_push(Hit *h, int *n, int k, Hit v) {
    if (*n < k) {
        h[*n] = v;
        heap_up(h, *n);
        (*n)++;
    } else if (v.d2 < h[0].d2) {
        h[0] = v;
        heap_down(h, *n, 0);
    }
}

static int cmp_hit(const void *a, const void *b) {
    double da = ((const Hit *)a)->d2;
    double db = ((const Hit *)b)->d2;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

int main(int argc, char **argv) {
    int k = 5;
    double tol = 1e-6;
    int use_pbc = 1;
    int json = 0;

    const char *pos[8];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage();
            return 0;
        } else if (!strcmp(a, "-k")) {
            if (ai + 1 >= argc) { usage(); return 1; }
            k = atoi(argv[++ai]);
        } else if (!strncmp(a, "-k=", 3)) {
            k = atoi(a + 3);
        } else if (!strncmp(a, "--tol=", 6)) {
            tol = strtod(a + 6, NULL);
        } else if (!strcmp(a, "--no-pbc")) {
            use_pbc = 0;
        } else if (!strcmp(a, "--json")) {
            json = 1;
        } else if (!strcmp(a, "--")) {
            /* -- ends options */
        } else if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-nearest-atoms: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 8) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if (npos != 4) { usage(); return 1; }
    if (k < 1) {
        fprintf(stderr, "lego-nearest-atoms: k must be >= 1\n");
        return 1;
    }

    double px = strtod(pos[0], NULL);
    double py = strtod(pos[1], NULL);
    double pz = strtod(pos[2], NULL);
    const char *input = pos[3];

    Config c;
    if (config_read(&c, input) != 0) return 1;
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-nearest-atoms: input has no x/y/z columns\n");
        config_free(&c);
        return 1;
    }
    if (c.natoms == 0) {
        fprintf(stderr, "lego-nearest-atoms: no atoms in input\n");
        config_free(&c);
        return 1;
    }

    double Lx = c.xhi - c.xlo;
    double Ly = c.yhi - c.ylo;
    double Lz = c.zhi - c.zlo;
    int pbc_x = use_pbc && Lx > 0.0;
    int pbc_y = use_pbc && Ly > 0.0;
    int pbc_z = use_pbc && Lz > 0.0;
    if (c.triclinic) {
        /* Minimum-image for triclinic is non-trivial; disable PBC in that
         * case and warn. */
        pbc_x = pbc_y = pbc_z = 0;
        if (use_pbc)
            fprintf(stderr, "lego-nearest-atoms: triclinic box, PBC disabled\n");
    }

    int k_cap = (int)((size_t)k < c.natoms ? (size_t)k : c.natoms);

    size_t N = c.natoms;
    int nc = c.ncols;
    int ix = c.idx_x, iy = c.idx_y, iz = c.idx_z;
    const double *D = c.data;

#ifdef _OPENMP
    int nth = omp_get_max_threads();
    if (nth < 1) nth = 1;
#else
    int nth = 1;
#endif
    Hit *tl_heaps = (Hit *)malloc((size_t)nth * (size_t)k_cap * sizeof(Hit));
    int *tl_n = (int *)calloc((size_t)nth, sizeof(int));
    if (!tl_heaps || !tl_n) {
        fprintf(stderr, "lego-nearest-atoms: oom\n");
        free(tl_heaps); free(tl_n);
        config_free(&c);
        return 1;
    }

#pragma omp parallel
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        Hit *h = &tl_heaps[(size_t)tid * (size_t)k_cap];
        int n = 0;

#pragma omp for nowait
        for (size_t i = 0; i < N; i++) {
            const double *row = D + i * (size_t)nc;
            double dx = row[ix] - px;
            double dy = row[iy] - py;
            double dz = row[iz] - pz;
            if (pbc_x) { dx -= Lx * rint(dx / Lx); }
            if (pbc_y) { dy -= Ly * rint(dy / Ly); }
            if (pbc_z) { dz -= Lz * rint(dz / Lz); }
            double d2 = dx*dx + dy*dy + dz*dz;
            Hit v = { d2, i };
            heap_push(h, &n, k_cap, v);
        }
        tl_n[tid] = n;
    }

    /* Merge thread-local heaps into a single global heap. */
    Hit *gh = (Hit *)malloc((size_t)k_cap * sizeof(Hit));
    if (!gh) {
        fprintf(stderr, "lego-nearest-atoms: oom\n");
        free(tl_heaps); free(tl_n);
        config_free(&c);
        return 1;
    }
    int gn = 0;
    for (int t = 0; t < nth; t++) {
        Hit *h = &tl_heaps[(size_t)t * (size_t)k_cap];
        for (int i = 0; i < tl_n[t]; i++)
            heap_push(gh, &gn, k_cap, h[i]);
    }
    free(tl_heaps);
    free(tl_n);

    /* Sort ascending by distance. */
    qsort(gh, (size_t)gn, sizeof(Hit), cmp_hit);

    double tol2 = tol * tol;
    int exact_hit = (gn > 0 && gh[0].d2 <= tol2);

    if (!json) {
        printf("query:   %.16g %.16g %.16g\n", px, py, pz);
        printf("file:    %s\n", input);
        printf("natoms:  %zu\n", c.natoms);
        printf("pbc:     %s%s%s%s\n",
               pbc_x ? "x" : "",
               pbc_y ? "y" : "",
               pbc_z ? "z" : "",
               (!pbc_x && !pbc_y && !pbc_z) ? "none" : "");
        printf("tol:     %.16g\n", tol);
        printf("exact:   %s", exact_hit ? "yes" : "no");
        if (exact_hit && c.idx_id >= 0) {
            long long id = (long long)D[gh[0].row * (size_t)nc + c.idx_id];
            printf(" (atom id %lld)", id);
        }
        printf("\n");

        printf("neighbours (%d):\n", gn);
        printf("  %-5s %-12s %-6s %-22s %-22s %-22s %-22s\n",
               "rank", "id", "type", "x", "y", "z", "distance");
        for (int i = 0; i < gn; i++) {
            const double *row = D + gh[i].row * (size_t)nc;
            long long id = c.idx_id   >= 0 ? (long long)row[c.idx_id]   : -1;
            long long tp = c.idx_type >= 0 ? (long long)row[c.idx_type] : -1;
            double x = row[ix], y = row[iy], z = row[iz];
            double dist = sqrt(gh[i].d2);
            printf("  %-5d %-12lld %-6lld %-22.16g %-22.16g %-22.16g %-22.16g\n",
                   i + 1, id, tp, x, y, z, dist);
        }
    } else {
        printf("{\n");
        printf("  \"file\": \"%s\",\n", input);
        printf("  \"query\": [%.16g, %.16g, %.16g],\n", px, py, pz);
        printf("  \"natoms\": %zu,\n", c.natoms);
        printf("  \"pbc\": [%s, %s, %s],\n",
               pbc_x ? "true" : "false",
               pbc_y ? "true" : "false",
               pbc_z ? "true" : "false");
        printf("  \"tol\": %.16g,\n", tol);
        printf("  \"exact\": %s,\n", exact_hit ? "true" : "false");
        printf("  \"neighbours\": [\n");
        for (int i = 0; i < gn; i++) {
            const double *row = D + gh[i].row * (size_t)nc;
            long long id = c.idx_id   >= 0 ? (long long)row[c.idx_id]   : -1;
            long long tp = c.idx_type >= 0 ? (long long)row[c.idx_type] : -1;
            double x = row[ix], y = row[iy], z = row[iz];
            double dist = sqrt(gh[i].d2);
            printf("    {\"rank\": %d, \"id\": %lld, \"type\": %lld, "
                   "\"x\": %.16g, \"y\": %.16g, \"z\": %.16g, "
                   "\"distance\": %.16g}%s\n",
                   i + 1, id, tp, x, y, z, dist,
                   i + 1 == gn ? "" : ",");
        }
        printf("  ]\n");
        printf("}\n");
    }

    free(gh);
    config_free(&c);
    return 0;
}
