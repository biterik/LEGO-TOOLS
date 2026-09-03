/* findroots.c — roots of a real polynomial by Laguerre's method.
 *
 * Faithful C port of Disloelast's FindRoots.f (Numerical Recipes zroots/
 * laguer).  The deflation order determines which member of each complex-
 * conjugate root pair lands in which slot, and the downstream solver in
 * sextic.c selects roots by position — so the algorithm, constants, and
 * iteration order are kept exactly as in the Fortran original.
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#include "sextic.h"

#include <complex.h>
#include <math.h>
#include <stdio.h>

#define FR_EPS   1.0e-9         /* convergence criterion (FindRoots.f EPS)  */
#define FR_EPSS  1.0e-14        /* round-off estimate    (LAGUER EPSS)      */
#define FR_MAXIT 100

/* Find one root of the complex polynomial a[0] + a[1]x + ... + a[m]x^m,
 * starting from *x.  Returns 0 on success, 1 if MAXIT is exceeded. */
static int laguer(const double _Complex *a, int m, double _Complex *x,
                  int polish) {
    double _Complex xx = *x;

    for (int iter = 1; iter <= FR_MAXIT; iter++) {
        double _Complex b = a[m];
        double err = cabs(b);
        double _Complex d = 0.0, f = 0.0;
        double abx = cabs(xx);
        for (int j = m - 1; j >= 0; j--) {
            f = xx * f + d;
            d = xx * d + b;
            b = xx * b + a[j];
            err = cabs(b) + abx * err;
        }
        err *= FR_EPSS;
        if (cabs(b) <= err) {
            *x = xx;
            return 0;
        }
        double _Complex g  = d / b;
        double _Complex g2 = g * g;
        double _Complex h  = g2 - 2.0 * f / b;
        double _Complex sq = csqrt((double)(m - 1) * ((double)m * h - g2));
        double _Complex gp = g + sq;
        double _Complex gm = g - sq;
        if (cabs(gp) < cabs(gm)) gp = gm;
        double _Complex dx = (double)m / gp;
        double _Complex x1 = xx - dx;
        if (xx == x1) {
            *x = xx;
            return 0;
        }
        xx = x1;
        if (!polish && cabs(dx) <= FR_EPS * cabs(xx)) {
            *x = xx;
            return 0;
        }
    }
    fprintf(stderr, "lego-dislo: laguer: too many iterations\n");
    *x = xx;
    return 1;
}

int findroots6(const double a[7], double _Complex root[6]) {
    const int m = 6;
    double _Complex ad[7];
    int ifail;

    /* First estimates: deflate from degree m down to 1, starting each
     * search at x = 0. */
    for (int j = 0; j <= m; j++) ad[j] = a[j];

    for (int j = m; j >= 1; j--) {
        double _Complex x = 0.0;
        ifail = laguer(ad, j, &x, 0);
        if (ifail) return ifail;
        if (fabs(cimag(x)) <= 2.0 * (FR_EPS * FR_EPS) * fabs(creal(x)))
            x = creal(x);
        root[j - 1] = x;
        /* Synthetic division: deflate the found root. */
        double _Complex b = ad[j];
        for (int jj = j - 1; jj >= 0; jj--) {
            double _Complex c = ad[jj];
            ad[jj] = b;
            b = x * b + c;
        }
    }

    /* Polish each root against the undeflated polynomial. */
    for (int j = 0; j <= m; j++) ad[j] = a[j];
    for (int j = 0; j < m; j++) {
        ifail = laguer(ad, m, &root[j], 1);
        if (ifail) return ifail;
    }
    return 0;
}
