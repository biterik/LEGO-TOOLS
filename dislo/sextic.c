/* sextic.c — anisotropic elastic dislocation solver + displacement field.
 *
 * Faithful C port of Disloelast (see sextic.h for the lineage).  The
 * algorithm — including the deliberate small-offset regularization of
 * near-degenerate sextic roots, the positional root selection, and the
 * branch-cut handling of the arctan terms — is kept exactly as in the
 * Fortran original; only the arithmetic is uniformly double precision
 * (the original mixed in single-precision literals such as pi).
 *
 * Two documented deviations from the literal Fortran (neither is reachable
 * for well-conditioned inputs; both fix latent bugs in dead branches):
 *   1. In the eigenvector fallback ladder (Fortran loop 37), branch 3 of
 *      the original tests AR13^2+AI23^2 (a typo mixing two subdeterminants);
 *      here the intended |D13|^2 is used.
 *   2. Fortran branch 2 (label 32/33) never assigns AI(2,N), silently
 *      reusing a stale value; here the mathematically consistent imaginary
 *      part is computed.
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
#include <string.h>

/* Fortran SUFIX mapping: 9-index I -> tensor pair (M,N), 1-based:
 *   1:(1,1) 2:(2,2) 3:(3,3) 4:(2,3) 5:(3,1) 6:(1,2) 7:(3,2) 8:(1,3) 9:(2,1)
 * Stored 0-based below. */
static const int SUF_M[9] = {0, 1, 2, 1, 2, 0, 2, 0, 1};
static const int SUF_N[9] = {0, 1, 2, 2, 0, 1, 1, 2, 0};

/* Voigt index (0-based) for each 9-index: pairs (2,3),(3,1),(1,2) and their
 * transposes share the Voigt slots 4,5,6 by the minor symmetry of C. */
static const int VOIGT9[9] = {0, 1, 2, 3, 4, 5, 3, 4, 5};

static void mat9_mul(const double A[9][9], const double B[9][9],
                     double R[9][9]) {
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++) {
            double s = 0.0;
            for (int k = 0; k < 9; k++) s += A[i][k] * B[k][j];
            R[i][j] = s;
        }
}

/* Solve the 6x6 real system A x = b by Gaussian elimination with partial
 * pivoting (equivalent to the bundled LAPACK dgesv of the original).
 * Returns 0 on success, nonzero if singular. */
static int solve6(double A[6][6], double b[6]) {
    int perm[6] = {0, 1, 2, 3, 4, 5};
    for (int col = 0; col < 6; col++) {
        int piv = col;
        double amax = fabs(A[perm[col]][col]);
        for (int r = col + 1; r < 6; r++) {
            double v = fabs(A[perm[r]][col]);
            if (v > amax) { amax = v; piv = r; }
        }
        if (amax == 0.0) return 1;
        int t = perm[col]; perm[col] = perm[piv]; perm[piv] = t;
        int pr = perm[col];
        for (int r = col + 1; r < 6; r++) {
            int rr = perm[r];
            double f = A[rr][col] / A[pr][col];
            A[rr][col] = 0.0;
            for (int c = col + 1; c < 6; c++) A[rr][c] -= f * A[pr][c];
            b[rr] -= f * b[pr];
        }
    }
    double x[6];
    for (int row = 5; row >= 0; row--) {
        int rr = perm[row];
        double s = b[rr];
        for (int c = row + 1; c < 6; c++) s -= A[rr][c] * x[c];
        x[row] = s / A[rr][row];
    }
    memcpy(b, x, sizeof(x));
    return 0;
}

static void set_err(char *err, int errlen, const char *msg) {
    if (err && errlen > 0) snprintf(err, (size_t)errlen, "%s", msg);
}

int sextic_solve(const double C6[6][6], const double T[3][3],
                 const double burgers[3], SexticSolution *sol,
                 char *err, int errlen) {
    /* ---- 9x9 elastic constant matrix in crystal axes -------------------- */
    double C9[9][9], Q[9][9], QT[9][9], R9[9][9], CP[9][9];
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            C9[i][j] = C6[VOIGT9[i]][VOIGT9[j]];

    /* Rotation into the dislocation frame: Q(I,J) = T(K,M)*T(L,N) with
     * (M,N) = sufix(I), (K,L) = sufix(J);  CP = Q^T * C9 * Q. */
    for (int i = 0; i < 9; i++) {
        int m = SUF_M[i], n = SUF_N[i];
        for (int j = 0; j < 9; j++) {
            int k = SUF_M[j], l = SUF_N[j];
            Q[i][j] = T[k][m] * T[l][n];
        }
    }
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            QT[i][j] = Q[j][i];
    mat9_mul(C9, Q, R9);
    mat9_mul(QT, R9, CP);

    /* CP is indexed 1-based in the comments below to match the Fortran. */
#define CPV(i, j) CP[(i) - 1][(j) - 1]

    /* ---- coefficients of the quadratics a_ik(p) ------------------------- *
     * AA[L][0..2]: a_ik(p) = AA0 + AA1*p + AA2*p^2, L = 3*(i-1)+(k-1). */
    double AA[9][3] = {
        {CPV(1,1), 2.0*CPV(1,6),          CPV(6,6)},
        {CPV(1,6), CPV(1,2) + CPV(6,6),   CPV(2,6)},
        {CPV(1,5), CPV(1,4) + CPV(5,6),   CPV(4,6)},
        {CPV(1,6), CPV(6,6) + CPV(1,2),   CPV(2,6)},
        {CPV(6,6), 2.0*CPV(2,6),          CPV(2,2)},
        {CPV(5,6), CPV(4,6) + CPV(2,5),   CPV(2,4)},
        {CPV(1,5), CPV(5,6) + CPV(1,4),   CPV(4,6)},
        {CPV(5,6), CPV(2,5) + CPV(4,6),   CPV(2,4)},
        {CPV(5,5), 2.0*CPV(4,5),          CPV(4,4)},
    };

    /* ---- sextic polynomial: det|a_ik(p)| = 0 ---------------------------- *
     * Six triple products with alternating signs (rows are 1-based L). */
    static const int DET_ROWS[6][3] = {
        {1,5,9}, {1,6,8}, {2,6,7}, {2,4,9}, {3,4,8}, {3,5,7}
    };
    static const double DET_SGN[6] = {+1, -1, +1, -1, +1, -1};

    double a[7] = {0};
    for (int t = 0; t < 6; t++) {
        const double *as = AA[DET_ROWS[t][0] - 1];
        const double *bs = AA[DET_ROWS[t][1] - 1];
        const double *cs = AA[DET_ROWS[t][2] - 1];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    a[i + j + k] += DET_SGN[t] * as[i] * bs[j] * cs[k];
    }

    double _Complex z[6];
    if (findroots6(a, z) != 0) {
        set_err(err, errlen, "sextic root finding failed (Laguerre did not "
                             "converge)");
        return 1;
    }

    /* ---- legacy degeneracy-breaking nudge ------------------------------- *
     * Real or imaginary parts below 1e-7 are set to the small offset 1e-4
     * (this is the 0.1000000E-03 visible in the legacy P output). */
    double X[6], Y[6];
    for (int i = 0; i < 6; i++) {
        X[i] = creal(z[i]);
        if (fabs(X[i]) < 1.0e-7) X[i] = 1.0e-4;
        Y[i] = cimag(z[i]);
        if (fabs(Y[i]) < 1.0e-7) Y[i] = 1.0e-4;
    }

    /* ---- positional root selection (slots 1,3,5 of the deflation order) - */
    X[1] = X[2]; Y[1] = Y[2];
    X[2] = X[4]; Y[2] = Y[4];
    if (X[1] < X[0]) {
        double tx = X[0], ty = Y[0];
        X[0] = X[1]; Y[0] = Y[1];
        X[1] = tx;   Y[1] = ty;
    }

    for (int n = 0; n < 3; n++) {
        sol->pr[n] = X[n];
        sol->pp[n] = Y[n];
    }
    memcpy(sol->burgers, burgers, 3 * sizeof(double));

    /* ---- eigenvectors A(k,n): a_ik(p_n) A_k = 0 ------------------------- */
    double _Complex Avec[3][3];         /* Avec[k][n] */
    for (int n = 0; n < 3; n++) {
        double _Complex av[3][3];       /* av[i][k] = a_ik(p_n) */
        for (int i = 0; i < 3; i++)
            for (int k = 0; k < 3; k++) {
                const double *q = AA[3 * i + k];
                double sa = q[0] + q[1] * X[n]
                          + q[2] * (X[n] * X[n] - Y[n] * Y[n]);
                double sb = q[1] * Y[n] + 2.0 * q[2] * X[n] * Y[n];
                av[i][k] = sa + sb * I;
            }
        /* Subdeterminants of rows 1,2 (Fortran DAN). */
        double _Complex d23 = av[0][1] * av[1][2] - av[0][2] * av[1][1];
        double _Complex d13 = av[0][0] * av[1][2] - av[0][2] * av[1][0];
        double _Complex d12 = av[0][0] * av[1][1] - av[0][1] * av[1][0];

        double m12 = creal(d12) * creal(d12) + cimag(d12) * cimag(d12);
        double m23 = creal(d23) * creal(d23) + cimag(d23) * cimag(d23);
        double m13 = creal(d13) * creal(d13) + cimag(d13) * cimag(d13);

        if (m12 >= 1.0e-8) {
            Avec[0][n] = d23 / d12;
            Avec[1][n] = -d13 / d12;
            Avec[2][n] = 1.0;
        } else if (m23 >= 1.0e-8) {
            Avec[0][n] = 1.0;
            Avec[1][n] = -d13 / d23;
            Avec[2][n] = d12 / d23;
        } else if (m13 >= 1.0e-8) {
            Avec[0][n] = -d23 / d13;
            Avec[1][n] = 1.0;
            Avec[2][n] = -d12 / d13;
        } else {
            /* Block-diagonal special case: a13,a23,a31,a32 all ~0 (e.g. a
             * pure screw with full decoupling); then A3 = 0 and
             * A2 = -(a11/a12) A1. */
            int all0 = 0;
            static const int CHK[4] = {2, 5, 6, 7};  /* rows a13,a23,a31,a32 */
            for (int r = 0; r < 4; r++)
                for (int q = 0; q < 3; q++)
                    if (fabs(AA[CHK[r]][q]) > 1.0e-4) all0 = 1;
            double m_a12 = creal(av[0][1]) * creal(av[0][1])
                         + cimag(av[0][1]) * cimag(av[0][1]);
            if (all0 || m_a12 < 1.0e-8) {
                set_err(err, errlen, "eq. 13-87 in Hirth & Lothe cannot be "
                                     "solved (degenerate eigensystem)");
                return 1;
            }
            Avec[0][n] = 1.0;
            Avec[1][n] = -av[0][0] / av[0][1];
            Avec[2][n] = 0.0;
        }
    }

    double AR[3][3], AI[3][3];          /* [k][n] */
    for (int n = 0; n < 3; n++)
        for (int k = 0; k < 3; k++) {
            AR[k][n] = creal(Avec[k][n]);
            AI[k][n] = cimag(Avec[k][n]);
        }

    double SGN[3];
    for (int n = 0; n < 3; n++) SGN[n] = (Y[n] < 0.0) ? -1.0 : 1.0;

    /* ---- B(k,n) = SGN(n) * conj( sum_j bb_kj(p_n) A_j(n) ) -------------- *
     * bb_kj(p) = bb0[k][j] + bb1[k][j]*p  (Fortran BBR/BBI). */
    double bb0[3][3] = {
        {CPV(1,6), CPV(6,6), CPV(5,6)},
        {CPV(1,2), CPV(2,6), CPV(2,5)},
        {CPV(1,4), CPV(4,6), CPV(4,5)},
    };
    double bb1[3][3] = {
        {CPV(6,6), CPV(2,6), CPV(4,6)},
        {CPV(2,6), CPV(2,2), CPV(2,4)},
        {CPV(4,6), CPV(2,4), CPV(4,4)},
    };

    double BRm[3][3], BIm[3][3];        /* [k][n] */
    for (int n = 0; n < 3; n++) {
        double _Complex p = X[n] + Y[n] * I;
        for (int k = 0; k < 3; k++) {
            double _Complex s = 0.0;
            for (int j = 0; j < 3; j++)
                s += (bb0[k][j] + bb1[k][j] * p) * Avec[j][n];
            BRm[k][n] = SGN[n] * creal(s);
            BIm[k][n] = -SGN[n] * cimag(s);
        }
    }

    /* ---- 6x6 linear system for D(n) ------------------------------------- */
    double ACUL[6][6], DSOL[6];
    for (int i = 0; i < 6; i++) DSOL[i] = 0.0;
    for (int i = 0; i < 3; i++) DSOL[i] = burgers[i];
    for (int k = 0; k < 3; k++)
        for (int n = 0; n < 3; n++) {
            ACUL[k][n]         = SGN[n] * AR[k][n];
            ACUL[k][n + 3]     = -SGN[n] * AI[k][n];
            ACUL[k + 3][n]     = BRm[k][n];
            ACUL[k + 3][n + 3] = BIm[k][n];
        }
    if (solve6(ACUL, DSOL) != 0) {
        set_err(err, errlen, "singular 6x6 system for the coefficients D(n)");
        return 1;
    }
    double DR[3], DI[3];
    for (int n = 0; n < 3; n++) {
        DR[n] = DSOL[n];
        DI[n] = DSOL[n + 3];
    }

    /* ---- F (log) and G (arctan) coefficient matrices -------------------- */
    for (int k = 0; k < 3; k++)
        for (int n = 0; n < 3; n++) {
            sol->F[k][n] = -0.5 * (AI[k][n] * DR[n] + AR[k][n] * DI[n]);
            sol->G[k][n] = -(AR[k][n] * DR[n] - AI[k][n] * DI[n]);
        }

    /* ---- energy prefactor K*b^2/4pi ------------------------------------- */
    sol->K_total = 0.0;
    for (int i = 0; i < 3; i++) {
        double esum = 0.0;
        for (int k = 0; k < 3; k++)
            for (int n = 0; n < 3; n++) {
                double bbr = bb0[i][k] + bb1[i][k] * X[n];
                double bbi = bb1[i][k] * Y[n];
                esum += bbr * (AR[k][n] * DI[n] + AI[k][n] * DR[n])
                      + bbi * (AR[k][n] * DR[n] - AI[k][n] * DI[n]);
            }
        sol->K_comp[i] = burgers[i] * esum / (4.0 * M_PI);
        sol->K_total += sol->K_comp[i];
    }

#undef CPV
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Field evaluation (port of anisodisloc.f)                               */
/* ---------------------------------------------------------------------- */

void sextic_phicut(const SexticSolution *sol, double ang_deg,
                   double phicut[3]) {
    /* Normalize the cut angle into [0,360) as the legacy code did. */
    double ang = ang_deg;
    int nang = (int)(ang / 360.0);
    if (ang >= 360.0) ang -= nang * 360.0;
    if (ang < 0.0)    ang -= (nang - 1) * 360.0;
    ang *= M_PI / 180.0;

    /* Map the cut direction (cos ANG, sin ANG) into each eta-plane and
     * take its argument in [0, 2pi). */
    for (int n = 0; n < 3; n++) {
        double etar = cos(ang) + sol->pr[n] * sin(ang);
        double etai = sol->pp[n] * sin(ang);
        if (etar == 0.0) etar += 1.0e-7;
        double ph = atan(etai / etar);
        if (etar < 0.0) ph += M_PI;
        if (etar > 0.0 && etai < 0.0) ph += 2.0 * M_PI;
        phicut[n] = ph;
    }
}

void sextic_displacement(const SexticSolution *sol, const double phicut[3],
                         double x, double y, double u[3]) {
    /* xs persists across the k/n loops: the legacy code nudges XSHIFT in
     * place when x + PR(n)*y hits exactly zero, and every later term of the
     * same atom sees the nudged value.  Kept verbatim. */
    double xs = x;
    for (int k = 0; k < 3; k++) {
        double uk = 0.0;
        for (int n = 0; n < 3; n++) {
            double arg1 = xs + sol->pr[n] * y;
            if (arg1 == 0.0) {
                xs += 1.0e-7;
                arg1 = xs + sol->pr[n] * y;
            }
            double arg2 = sol->pp[n] * y;

            /* Argument of (arg1 + i*arg2) with the single 2pi discontinuity
             * along the +x axis of the eta-plane, i.e. in [0, 2pi)... */
            double arct = atan(arg2 / arg1);
            if (arg1 < 0.0) arct += M_PI;
            if (arg1 > 0.0 && arg2 < 0.0) arct += 2.0 * M_PI;
            /* ...then moved so the discontinuity lies on the chosen cut. */
            if (arct < phicut[n]) arct += 2.0 * M_PI;

            uk += sol->F[k][n] * log(arg1 * arg1 + arg2 * arg2)
                + sol->G[k][n] * arct;
        }
        u[k] = -uk / (2.0 * M_PI);
    }
}
