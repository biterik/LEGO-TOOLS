/* lego-dislo — insert dislocations into atomistic configurations using the
 * anisotropic elastic solution (sextic method, Hirth & Lothe).
 *
 * Usage:
 *   lego-dislo [options] CARD INPUT OUTPUT
 *   lego-dislo --solve-only [options] CARD
 *
 *   --in-format=lammps|lammps-dump|imd   (default: auto-detect, via atomio)
 *   --out-format=lammps|lammps-dump|imd  (default: lammps)
 *   --threads=N                          (default: OMP_NUM_THREADS)
 *   --report=FILE                        (write K-factors, P/F/G per disloc.)
 *   --solve-only                         (solve + report only, no atoms)
 *   -h, --help
 *
 * For every atom, the displacement field of each dislocation is evaluated
 * in that dislocation's own frame and the results are superposed (plain
 * linear superposition, no periodic-image corrections).  Box bounds are
 * left unchanged — cut/shift afterwards with lego-cut / lego-shift.
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "../lego-tools/atomio.h"
#include "card.h"
#include "sextic.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    SexticSolution sol;
    double phicut[3];
    double M[3][3];        /* box frame -> dislocation frame */
    double point[3];       /* core point, box coordinates, Angstrom */
} SolvedDislo;

static void usage(void) {
    fprintf(stderr,
        "Usage: lego-dislo [options] CARD INPUT OUTPUT\n"
        "       lego-dislo --solve-only [options] CARD\n"
        "\n"
        "Insert dislocations (anisotropic elastic displacement field, sextic\n"
        "method) into an atomistic configuration.  Box bounds are left\n"
        "unchanged; cut/shift afterwards with lego-cut / lego-shift.\n"
        "\n"
        "Options:\n"
        "  --in-format=FMT    lammps|lammps-dump|imd  (default: auto-detect)\n"
        "  --out-format=FMT   lammps|lammps-dump|imd  (default: lammps)\n"
        "  --threads=N        number of OpenMP threads\n"
        "  --report=FILE      write the solver report (P/F/G, K) to FILE\n"
        "  --solve-only       solve and report only; no INPUT/OUTPUT\n"
        "  -h, --help\n");
}

static void norm3(double v[3]) {
    double s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= s; v[1] /= s; v[2] /= s;
}

static double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void cross3(const double a[3], const double b[3], double r[3]) {
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

static int zero3(const double v[3]) {
    return v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0;
}

#define ORTHO_TOL 1e-8

/* Build the box rotation R rows = unit box axes in crystal coordinates.
 * Errors out unless the axes are orthogonal and right-handed to 1e-8. */
static int build_box_rotation(const Card *card, double R[3][3]) {
    if (zero3(card->box_x) || zero3(card->box_y) || zero3(card->box_z)) {
        fprintf(stderr, "lego-dislo: box_orientation axis is a zero vector\n");
        return 1;
    }
    memcpy(R[0], card->box_x, sizeof(R[0]));
    memcpy(R[1], card->box_y, sizeof(R[1]));
    memcpy(R[2], card->box_z, sizeof(R[2]));
    for (int i = 0; i < 3; i++) norm3(R[i]);
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            if (fabs(dot3(R[i], R[j])) > ORTHO_TOL) {
                fprintf(stderr, "lego-dislo: box_orientation axes %c and %c "
                        "are not orthogonal (|cos| = %.3e)\n",
                        'x' + i, 'x' + j, fabs(dot3(R[i], R[j])));
                return 1;
            }
    double zc[3];
    cross3(R[0], R[1], zc);
    if (fabs(zc[0] - R[2][0]) > ORTHO_TOL ||
        fabs(zc[1] - R[2][1]) > ORTHO_TOL ||
        fabs(zc[2] - R[2][2]) > ORTHO_TOL) {
        fprintf(stderr, "lego-dislo: box_orientation is not right-handed "
                "(z != x cross y)\n");
        return 1;
    }
    return 0;
}

/* Solve one dislocation: frame setup, sextic solution, frame composition. */
static int solve_dislo(const Card *card, int idx, const double Rbc[3][3],
                       SolvedDislo *sd) {
    const CardDislo *d = &card->dislo[idx];

    if (zero3(d->line) || zero3(d->glide_normal) || zero3(d->b_dir)) {
        fprintf(stderr, "lego-dislo: dislocation %d: zero direction vector\n",
                idx + 1);
        return 1;
    }

    /* Dislocation frame in crystal coordinates: Z = line, Y = glide-plane
     * normal, X = Y x Z (right-handed). */
    double T[3][3];
    memcpy(T[2], d->line, sizeof(T[2]));
    memcpy(T[1], d->glide_normal, sizeof(T[1]));
    norm3(T[2]);
    norm3(T[1]);
    if (fabs(dot3(T[1], T[2])) > ORTHO_TOL) {
        fprintf(stderr, "lego-dislo: dislocation %d: glide_normal is not "
                "perpendicular to line (|cos| = %.3e)\n",
                idx + 1, fabs(dot3(T[1], T[2])));
        return 1;
    }
    cross3(T[1], T[2], T[0]);
    norm3(T[0]);

    /* Burgers vector components in the dislocation frame (latpar units). */
    double bu[3];
    memcpy(bu, d->b_dir, sizeof(bu));
    norm3(bu);
    double bur[3];
    for (int i = 0; i < 3; i++) bur[i] = d->b_mag * dot3(bu, T[i]);

    char err[256];
    if (sextic_solve(card->C6, (const double (*)[3])T, bur, &sd->sol,
                     err, sizeof(err))) {
        fprintf(stderr, "lego-dislo: dislocation %d: %s\n", idx + 1, err);
        return 1;
    }
    sextic_phicut(&sd->sol, d->cut_angle, sd->phicut);

    /* Composition box -> crystal -> dislocation: M = T * Rbc^T. */
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            sd->M[i][j] = T[i][0]*Rbc[j][0] + T[i][1]*Rbc[j][1]
                        + T[i][2]*Rbc[j][2];
    memcpy(sd->point, d->point, sizeof(sd->point));
    return 0;
}

static void report_write(FILE *fp, const Card *card, const SolvedDislo *sd) {
    fprintf(fp, "lego-dislo report\n");
    fprintf(fp, "material: %s\n", card->name[0] ? card->name : "(unnamed)");
    fprintf(fp, "lattice_parameter_A: %.10g\n", card->a0);
    fprintf(fp, "dislocations: %d\n", card->ndislo);
    for (int i = 0; i < card->ndislo; i++) {
        const CardDislo *d = &card->dislo[i];
        const SexticSolution *s = &sd[i].sol;
        double b2 = s->burgers[0]*s->burgers[0] + s->burgers[1]*s->burgers[1]
                  + s->burgers[2]*s->burgers[2];
        fprintf(fp, "\ndislocation %d\n", i + 1);
        fprintf(fp, "  point_A: %.10g %.10g %.10g\n",
                d->point[0], d->point[1], d->point[2]);
        fprintf(fp, "  line: %g %g %g\n",
                d->line[0], d->line[1], d->line[2]);
        fprintf(fp, "  glide_normal: %g %g %g\n",
                d->glide_normal[0], d->glide_normal[1], d->glide_normal[2]);
        fprintf(fp, "  burgers_vector_lp: %.10g %.10g %.10g\n",
                s->burgers[0], s->burgers[1], s->burgers[2]);
        fprintf(fp, "  burgers_magnitude_lp: %.10g\n", d->b_mag);
        fprintf(fp, "  cut_angle_deg: %.10g\n", d->cut_angle);
        for (int n = 0; n < 3; n++)
            fprintf(fp, "  P_%d: %15.7e %15.7e\n", n + 1,
                    s->pr[n], s->pp[n]);
        for (int k = 0; k < 3; k++)
            fprintf(fp, "  F_%d: %15.7e %15.7e %15.7e\n", k + 1,
                    s->F[k][0], s->F[k][1], s->F[k][2]);
        for (int k = 0; k < 3; k++)
            fprintf(fp, "  G_%d: %15.7e %15.7e %15.7e\n", k + 1,
                    s->G[k][0], s->G[k][1], s->G[k][2]);
        /* K_total is in GPa * latpar^2 (C in GPa, b in latpar units).
         * J/m:  * 1e9 Pa/GPa * (a0 * 1e-10 m)^2
         * K in GPa: * 4pi / |b_lp|^2 */
        double a0m = card->a0 * 1e-10;
        fprintf(fp, "  K_b2_over_4pi_GPa_lp2: %.7e\n", s->K_total);
        fprintf(fp, "  K_b2_over_4pi_J_per_m: %.7e\n",
                s->K_total * 1e9 * a0m * a0m);
        if (b2 > 0.0)
            fprintf(fp, "  K_GPa: %.7e\n", s->K_total * 4.0 * M_PI / b2);
    }
}

int main(int argc, char **argv) {
    const char *in_fmt_s = NULL;
    const char *out_fmt_s = "lammps";
    const char *report_file = NULL;
    int solve_only = 0;
    long threads = 0;
    const char *pos[8];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage();
            return 0;
        } else if (!strncmp(a, "--in-format=", 12)) {
            in_fmt_s = a + 12;
        } else if (!strncmp(a, "--out-format=", 13)) {
            out_fmt_s = a + 13;
        } else if (!strncmp(a, "--threads=", 10)) {
            threads = strtol(a + 10, NULL, 10);
        } else if (!strncmp(a, "--report=", 9)) {
            report_file = a + 9;
        } else if (!strcmp(a, "--solve-only")) {
            solve_only = 1;
        } else if (!strcmp(a, "--")) {
            /* -- ends options */
        } else if (a[0] == '-' && a[1] != 0) {
            fprintf(stderr, "lego-dislo: unknown option %s\n", a);
            return 1;
        } else {
            if (npos >= 8) { usage(); return 1; }
            pos[npos++] = a;
        }
    }
    if ((solve_only && npos != 1) || (!solve_only && npos != 3)) {
        usage();
        return 1;
    }
    const char *card_file = pos[0];

#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads((int)threads);
#else
    (void)threads;
#endif

    Format out_fmt = format_from_string(out_fmt_s);
    if (out_fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-dislo: unknown output format '%s'\n", out_fmt_s);
        return 1;
    }
    if (in_fmt_s && format_from_string(in_fmt_s) == FMT_UNKNOWN) {
        fprintf(stderr, "lego-dislo: unknown input format '%s'\n", in_fmt_s);
        return 1;
    }

    Card card;
    char err[256];
    if (card_parse(card_file, &card, err, sizeof(err))) {
        fprintf(stderr, "lego-dislo: %s\n", err);
        return 1;
    }

    double Rbc[3][3];
    if (build_box_rotation(&card, Rbc)) return 1;

    SolvedDislo *sd = (SolvedDislo *)malloc((size_t)card.ndislo * sizeof(*sd));
    if (!sd) { fprintf(stderr, "lego-dislo: oom\n"); return 1; }
    for (int i = 0; i < card.ndislo; i++)
        if (solve_dislo(&card, i, (const double (*)[3])Rbc, &sd[i])) {
            free(sd);
            return 1;
        }

    report_write(stdout, &card, sd);
    if (report_file) {
        FILE *fp = fopen(report_file, "w");
        if (!fp) {
            fprintf(stderr, "lego-dislo: cannot write report '%s'\n",
                    report_file);
            free(sd);
            return 1;
        }
        report_write(fp, &card, sd);
        fclose(fp);
    }
    if (solve_only) { free(sd); return 0; }

    const char *input  = pos[1];
    const char *output = pos[2];

    if (in_fmt_s) {
        Format want = format_from_string(in_fmt_s);
        Format got  = detect_format(input);
        if (got != FMT_UNKNOWN && got != want)
            fprintf(stderr, "lego-dislo: warning: --in-format=%s but '%s' "
                    "looks like %s (atomio auto-detects; proceeding)\n",
                    in_fmt_s, input, format_name(got));
    }

    Config c;
    if (config_read(&c, input) != 0) { free(sd); return 1; }
    if (c.idx_x < 0 || c.idx_y < 0 || c.idx_z < 0) {
        fprintf(stderr, "lego-dislo: input has no x/y/z columns\n");
        config_free(&c);
        free(sd);
        return 1;
    }

    size_t N = c.natoms;
    int nc = c.ncols;
    int ix = c.idx_x, iy = c.idx_y, iz = c.idx_z;
    double *D = c.data;
    double a0 = card.a0;
    int nd = card.ndislo;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < N; i++) {
        double *row = D + i * (size_t)nc;
        double rx = row[ix], ry = row[iy], rz = row[iz];
        double ax = 0.0, ay = 0.0, az = 0.0;   /* accumulated displacement */
        for (int d = 0; d < nd; d++) {
            const SolvedDislo *s = &sd[d];
            double vx = rx - s->point[0];
            double vy = ry - s->point[1];
            double vz = rz - s->point[2];
            /* Into the dislocation frame; only the in-plane coordinates
             * matter (the line is along the local Z axis). */
            double x = (s->M[0][0]*vx + s->M[0][1]*vy + s->M[0][2]*vz) / a0;
            double y = (s->M[1][0]*vx + s->M[1][1]*vy + s->M[1][2]*vz) / a0;
            double u[3];
            sextic_displacement(&s->sol, s->phicut, x, y, u);
            u[0] *= a0; u[1] *= a0; u[2] *= a0;
            /* Back to the box frame (M is orthogonal: inverse = transpose) */
            ax += s->M[0][0]*u[0] + s->M[1][0]*u[1] + s->M[2][0]*u[2];
            ay += s->M[0][1]*u[0] + s->M[1][1]*u[1] + s->M[2][1]*u[2];
            az += s->M[0][2]*u[0] + s->M[1][2]*u[1] + s->M[2][2]*u[2];
        }
        row[ix] = rx + ax;
        row[iy] = ry + ay;
        row[iz] = rz + az;
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    double dt = (double)(t1.tv_sec - t0.tv_sec)
              + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec);
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    fprintf(stderr, "lego-dislo: displaced %zu atoms (%d dislocation%s) "
            "in %.3f s field time, %d thread%s\n",
            N, nd, nd == 1 ? "" : "s", dt, nthreads,
            nthreads == 1 ? "" : "s");

    int rc = config_write(&c, output, out_fmt);
    config_free(&c);
    free(sd);
    return rc;
}
