/* card.h — input card parser for lego-dislo.
 *
 * The card uses a YAML-look syntax restricted to a small, flat subset so
 * the parser stays dependency-free.  See dislo/README.md for the exact
 * grammar; the Python driver only ever emits this subset.
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#ifndef CARD_H
#define CARD_H

#define CARD_MAX_DISLO 64

typedef struct {
    double point[3];          /* point on the line, box coordinates, A     */
    double line[3];           /* line direction, crystal coordinates       */
    double b_dir[3];          /* Burgers direction, crystal coordinates    */
    double b_mag;             /* |b| in lattice-parameter units            */
    double glide_normal[3];   /* glide-plane normal, crystal coordinates   */
    double cut_angle;         /* branch-cut angle, degrees (default 0)     */
} CardDislo;

typedef struct {
    char   name[64];          /* material name (informational)             */
    double a0;                /* lattice parameter, Angstrom               */
    double C6[6][6];          /* stiffness, Voigt, crystal frame, GPa      */

    double box_x[3];          /* crystallographic direction of box x axis  */
    double box_y[3];
    double box_z[3];

    int       ndislo;
    CardDislo dislo[CARD_MAX_DISLO];
} Card;

/* Parse the card file.  Returns 0 on success; on failure returns nonzero
 * and writes a message (with line number) into err. */
int card_parse(const char *path, Card *card, char *err, int errlen);

#endif /* CARD_H */
