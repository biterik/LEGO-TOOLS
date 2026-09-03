/* gen_bcc.c — fast synthetic bcc sample generator for benchmarking.
 *
 * Writes a LAMMPS data file with nx*ny*nz bcc unit cells (2 atoms/cell)
 * along the cube axes.  Used only by bench.sh; kept deliberately trivial
 * (the field evaluation does not care about the orientation, only about
 * the atom count).
 *
 * Usage: gen_bcc nx ny nz a0 OUTPUT
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "Usage: gen_bcc nx ny nz a0 OUTPUT\n");
        return 1;
    }
    long nx = atol(argv[1]), ny = atol(argv[2]), nz = atol(argv[3]);
    double a0 = atof(argv[4]);
    const char *out = argv[5];
    if (nx <= 0 || ny <= 0 || nz <= 0 || a0 <= 0.0) {
        fprintf(stderr, "gen_bcc: bad arguments\n");
        return 1;
    }

    FILE *fp = fopen(out, "w");
    if (!fp) { fprintf(stderr, "gen_bcc: cannot write %s\n", out); return 1; }
    setvbuf(fp, NULL, _IOFBF, 1 << 22);

    long n = 2 * nx * ny * nz;
    fprintf(fp, "bcc block %ldx%ldx%ld a0=%g (gen_bcc)\n\n", nx, ny, nz, a0);
    fprintf(fp, "%ld atoms\n1 atom types\n\n", n);
    fprintf(fp, "%.10g %.10g xlo xhi\n", -0.5 * nx * a0, 0.5 * nx * a0);
    fprintf(fp, "%.10g %.10g ylo yhi\n", -0.5 * ny * a0, 0.5 * ny * a0);
    fprintf(fp, "%.10g %.10g zlo zhi\n\n", -0.5 * nz * a0, 0.5 * nz * a0);
    fprintf(fp, "Masses\n\n1 183.84\n\nAtoms # atomic\n\n");

    long id = 0;
    double ox = -0.5 * nx * a0, oy = -0.5 * ny * a0, oz = -0.5 * nz * a0;
    for (long i = 0; i < nx; i++)
        for (long j = 0; j < ny; j++)
            for (long k = 0; k < nz; k++) {
                double x = ox + i * a0, y = oy + j * a0, z = oz + k * a0;
                fprintf(fp, "%ld 1 %.6f %.6f %.6f\n", ++id, x, y, z);
                fprintf(fp, "%ld 1 %.6f %.6f %.6f\n", ++id,
                        x + 0.5 * a0, y + 0.5 * a0, z + 0.5 * a0);
            }

    if (fclose(fp) != 0) { fprintf(stderr, "gen_bcc: write failed\n");
                           return 1; }
    fprintf(stderr, "gen_bcc: wrote %ld atoms to %s\n", n, out);
    return 0;
}
