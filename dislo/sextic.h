/* sextic.h — anisotropic elastic dislocation displacement field
 *            (sextic method, Hirth & Lothe, "Theory of Dislocations").
 *
 * Faithful C port of the 2005-era Fortran toolkit "Disloelast":
 *   elast-cubic.f  (solver; lineage: Alexey Girshick / U. Penn,
 *                   Bart Pestman 1988, Erik Bitzek 2005)
 *   FindRoots.f    (Laguerre root finder, Numerical Recipes style)
 *   anisodisloc.f  (field evaluation with branch-cut handling)
 *
 * The solver takes the full 6x6 stiffness matrix (Voigt notation) in the
 * crystal frame plus the dislocation coordinate frame, solves the sextic
 * eigenproblem and returns the three complex roots p(n), the coefficient
 * matrices F(k,n) (log terms) and G(k,n) (arctan terms), and the energy
 * prefactor K*b^2/4pi.
 *
 * The displacement of a point at in-plane position (x,y) — measured in the
 * dislocation frame, in units of the lattice parameter, relative to the
 * core — is
 *
 *   u_k = -(1/2pi) SUM_n [ F(k,n)*ln((x+PR(n)*y)^2 + (PP(n)*y)^2)
 *                         + G(k,n)*arctg_n(x,y) ]
 *
 * where arctg_n is the argument of (x+PR(n)*y) + i*PP(n)*y, branch-corrected
 * into [phicut(n), phicut(n)+2pi) so that the single 2pi discontinuity lies
 * along the direction ANG (the cut angle) in real (x,y) space.
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#ifndef SEXTIC_H
#define SEXTIC_H

#include <complex.h>

typedef struct {
    /* Roots p(n) = pr(n) + i*pp(n) of the sextic, after the legacy
     * degeneracy-breaking nudge and root selection. */
    double pr[3], pp[3];

    /* Coefficient matrices, indexed [k][n] (displacement component k,
     * root n) exactly as in the legacy channel-21 output. */
    double F[3][3];             /* log terms    */
    double G[3][3];             /* arctan terms */

    /* Burgers vector in the dislocation frame (lattice-parameter units),
     * copied from the input. */
    double burgers[3];

    /* Energy prefactor K*b^2/4pi, per component of b and total, in
     * (units of C) * (units of b)^2.  With C in GPa and b in lattice-
     * parameter units, multiply by 1e9*(a0*1e-10)^2 for J/m. */
    double K_comp[3];
    double K_total;
} SexticSolution;

/* Find the 6 complex roots of the real degree-6 polynomial
 *   a[0] + a[1]*p + ... + a[6]*p^6.
 * Faithful port of FindRoots.f (Laguerre with deflation + polishing);
 * the ROOT ORDER matters downstream and matches the legacy code.
 * Returns 0 on success, nonzero on failure. */
int findroots6(const double a[7], double _Complex root[6]);

/* Solve the sextic eigenproblem.
 *   C6      6x6 stiffness matrix, Voigt notation, crystal frame (any units;
 *           the displacement field is invariant to the overall scale of C).
 *   T       rows = unit vectors of the dislocation frame axes X, Y, Z
 *           expressed in crystal coordinates (X = glide dir, Y = glide-plane
 *           normal, Z = line direction; right-handed, X = Y x Z).
 *   burgers Burgers vector components in the dislocation frame
 *           (lattice-parameter units).
 * Returns 0 on success; on failure returns nonzero and writes a message
 * into err (if err != NULL, up to errlen bytes). */
int sextic_solve(const double C6[6][6], const double T[3][3],
                 const double burgers[3], SexticSolution *sol,
                 char *err, int errlen);

/* Precompute the branch-cut angles phicut[n] in each eta-plane for a cut
 * along direction ang_deg (degrees from +X in the real (x,y) plane). */
void sextic_phicut(const SexticSolution *sol, double ang_deg,
                   double phicut[3]);

/* Evaluate the displacement u[3] (lattice-parameter units) at in-plane
 * position (x,y) (lattice-parameter units, relative to the core), using
 * the legacy branch-cut logic. */
void sextic_displacement(const SexticSolution *sol, const double phicut[3],
                         double x, double y, double u[3]);

#endif /* SEXTIC_H */
