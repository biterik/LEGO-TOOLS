/* card.c — input card parser for lego-dislo (see card.h).
 *
 * Grammar (a strict subset of YAML; two-space indentation):
 *
 *   material:
 *     name: W
 *     lattice_parameter: 3.1652
 *     cij:
 *       cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
 *       # or:
 *       full:
 *         - [c11, c12, c13, c14, c15, c16]
 *         - ... (6 rows)
 *   box_orientation:
 *     x: [1, 1, 1]
 *     y: [-1, 0, 1]
 *     z: [1, -2, 1]
 *   dislocations:
 *     - point: [0.0, 0.0, 0.0]
 *       line: [1, -2, 1]
 *       burgers_direction: [-1, -1, -1]
 *       burgers_magnitude: 0.8663
 *       glide_normal: [-1, 0, 1]
 *       cut_angle: 0.0        # optional
 *
 * '#' starts a comment; blank lines are ignored; tabs are rejected.
 *
 * License: GPL-3.0
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#include "card.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char err[256];
    int  lineno;
    int  failed;
} PErr;

static void perr(PErr *pe, const char *fmt, const char *arg) {
    if (pe->failed) return;
    pe->failed = 1;
    char msg[192];
    snprintf(msg, sizeof(msg), fmt, arg ? arg : "");
    snprintf(pe->err, sizeof(pe->err), "card line %d: %s", pe->lineno, msg);
}

/* Strip comment + trailing whitespace in place; return indent (spaces). */
static int clean_line(char *s, PErr *pe) {
    char *h = strchr(s, '#');
    if (h) *h = 0;
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == ' ' || s[n-1] == '\t' ||
                     s[n-1] == '\r' || s[n-1] == '\n'))
        s[--n] = 0;
    int ind = 0;
    while (s[ind] == ' ') ind++;
    if (s[ind] == '\t') perr(pe, "tab indentation is not supported", NULL);
    return ind;
}

/* Parse "[a, b, ...]" into out[0..n-1]; requires exactly n numbers. */
static int parse_list(const char *s, double *out, int n, PErr *pe) {
    while (*s == ' ') s++;
    if (*s != '[') { perr(pe, "expected a list like [1, 2, 3]", NULL); return 1; }
    s++;
    for (int i = 0; i < n; i++) {
        char *end;
        double v = strtod(s, &end);
        if (end == s) { perr(pe, "expected %s numbers in list",
                              n == 3 ? "3" : "6"); return 1; }
        out[i] = v;
        s = end;
        while (*s == ' ') s++;
        if (i < n - 1) {
            if (*s != ',') { perr(pe, "expected ',' in list", NULL); return 1; }
            s++;
        }
    }
    while (*s == ' ') s++;
    if (*s != ']') { perr(pe, "expected ']' closing the list", NULL); return 1; }
    s++;
    while (*s == ' ') s++;
    if (*s) { perr(pe, "unexpected text after list", NULL); return 1; }
    return 0;
}

/* Parse "{k1: v1, k2: v2, ...}"; returns number of pairs or -1. */
#define MAX_PAIRS 8
static int parse_map(const char *s, char keys[MAX_PAIRS][16],
                     double vals[MAX_PAIRS], PErr *pe) {
    while (*s == ' ') s++;
    if (*s != '{') { perr(pe, "expected a map like {c11: 522.4, ...}", NULL);
                     return -1; }
    s++;
    int np = 0;
    for (;;) {
        while (*s == ' ') s++;
        if (*s == '}') { s++; break; }
        if (np >= MAX_PAIRS) { perr(pe, "too many entries in map", NULL);
                               return -1; }
        int kl = 0;
        while (*s && *s != ':' && *s != ' ' && kl < 15) keys[np][kl++] = *s++;
        keys[np][kl] = 0;
        while (*s == ' ') s++;
        if (*s != ':') { perr(pe, "expected ':' in map", NULL); return -1; }
        s++;
        char *end;
        vals[np] = strtod(s, &end);
        if (end == s) { perr(pe, "expected a number for map key '%s'",
                              keys[np]); return -1; }
        s = end;
        np++;
        while (*s == ' ') s++;
        if (*s == ',') { s++; continue; }
    }
    while (*s == ' ') s++;
    if (*s) { perr(pe, "unexpected text after map", NULL); return -1; }
    return np;
}

static double parse_scalar(const char *s, PErr *pe) {
    char *end;
    double v = strtod(s, &end);
    if (end == s) { perr(pe, "expected a number, got '%s'", s); return 0.0; }
    while (*end == ' ') end++;
    if (*end) { perr(pe, "unexpected text after number: '%s'", end); return 0.0; }
    return v;
}

enum Section { S_NONE, S_MATERIAL, S_CIJ, S_CIJ_FULL, S_BOX, S_DISLOS };

int card_parse(const char *path, Card *card, char *err, int errlen) {
    FILE *fp = fopen(path, "r");
    if (!fp) {
        snprintf(err, (size_t)errlen, "cannot open card file '%s'", path);
        return 1;
    }

    memset(card, 0, sizeof(*card));
    card->a0 = 0.0;

    PErr pe = {{0}, 0, 0};
    enum Section sec = S_NONE;
    int have_cij = 0, have_a0 = 0;
    int have_bx = 0, have_by = 0, have_bz = 0;
    int full_rows = 0;
    unsigned dislo_keys = 0;      /* bitmask of required keys of current item */
#define DK_POINT 1u
#define DK_LINE 2u
#define DK_BDIR 4u
#define DK_BMAG 8u
#define DK_GLIDE 16u

    char line[1024];
    while (!pe.failed && fgets(line, sizeof(line), fp)) {
        pe.lineno++;
        int ind = clean_line(line, &pe);
        if (pe.failed) break;
        char *s = line + ind;
        if (!*s) continue;

        /* List row under cij: full: */
        if (sec == S_CIJ_FULL && s[0] == '-' && s[1] != 0 &&
            strchr(s, ':') == NULL) {
            if (full_rows >= 6) { perr(&pe, "more than 6 rows under 'full:'",
                                       NULL); break; }
            parse_list(s + 1, card->C6[full_rows], 6, &pe);
            full_rows++;
            if (full_rows == 6) { have_cij = 1; sec = S_CIJ; }
            continue;
        }

        /* Dislocation list item: "- key: value" */
        int is_item = (s[0] == '-' && s[1] == ' ');
        if (is_item) {
            if (sec != S_DISLOS) { perr(&pe, "list item outside a list "
                                        "section", NULL); break; }
            if (card->ndislo > 0 &&
                (dislo_keys & (DK_POINT|DK_LINE|DK_BDIR|DK_BMAG|DK_GLIDE))
                    != (DK_POINT|DK_LINE|DK_BDIR|DK_BMAG|DK_GLIDE)) {
                perr(&pe, "previous dislocation is missing a required key",
                     NULL);
                break;
            }
            if (card->ndislo >= CARD_MAX_DISLO) {
                perr(&pe, "too many dislocations", NULL);
                break;
            }
            card->ndislo++;
            dislo_keys = 0;
            card->dislo[card->ndislo - 1].cut_angle = 0.0;
            s += 2;
            while (*s == ' ') s++;
        }

        /* Split "key: value" */
        char *colon = strchr(s, ':');
        if (!colon) { perr(&pe, "expected 'key: value', got '%s'", s); break; }
        *colon = 0;
        char *key = s;
        char *val = colon + 1;
        while (*val == ' ') val++;

        if (ind == 0 && !is_item) {
            if (!strcmp(key, "material"))             sec = S_MATERIAL;
            else if (!strcmp(key, "box_orientation")) sec = S_BOX;
            else if (!strcmp(key, "dislocations"))    sec = S_DISLOS;
            else { perr(&pe, "unknown top-level key '%s'", key); break; }
            if (*val) { perr(&pe, "unexpected value after '%s:'", key); break; }
            continue;
        }

        switch (sec) {
        case S_MATERIAL:
        case S_CIJ:
        case S_CIJ_FULL:
            if (!strcmp(key, "name")) {
                snprintf(card->name, sizeof(card->name), "%s", val);
            } else if (!strcmp(key, "lattice_parameter")) {
                card->a0 = parse_scalar(val, &pe);
                have_a0 = 1;
            } else if (!strcmp(key, "cij")) {
                if (*val) { perr(&pe, "unexpected value after 'cij:'", NULL);
                            break; }
                sec = S_CIJ;
            } else if (!strcmp(key, "cubic")) {
                if (sec != S_CIJ && sec != S_CIJ_FULL) {
                    perr(&pe, "'cubic' outside 'cij:'", NULL); break;
                }
                char keys[MAX_PAIRS][16];
                double vals[MAX_PAIRS];
                int np = parse_map(val, keys, vals, &pe);
                if (np < 0) break;
                double c11 = 0, c12 = 0, c44 = 0;
                int h11 = 0, h12 = 0, h44 = 0;
                for (int i = 0; i < np; i++) {
                    if (!strcmp(keys[i], "c11")) { c11 = vals[i]; h11 = 1; }
                    else if (!strcmp(keys[i], "c12")) { c12 = vals[i]; h12 = 1; }
                    else if (!strcmp(keys[i], "c44")) { c44 = vals[i]; h44 = 1; }
                    else { perr(&pe, "unknown cubic constant '%s'", keys[i]);
                           break; }
                }
                if (pe.failed) break;
                if (!h11 || !h12 || !h44) {
                    perr(&pe, "cubic needs c11, c12 and c44", NULL); break;
                }
                memset(card->C6, 0, sizeof(card->C6));
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        card->C6[i][j] = (i == j) ? c11 : c12;
                for (int i = 3; i < 6; i++) card->C6[i][i] = c44;
                have_cij = 1;
            } else if (!strcmp(key, "full")) {
                if (sec != S_CIJ) { perr(&pe, "'full' outside 'cij:'", NULL);
                                    break; }
                if (*val) { perr(&pe, "unexpected value after 'full:'", NULL);
                            break; }
                memset(card->C6, 0, sizeof(card->C6));
                full_rows = 0;
                sec = S_CIJ_FULL;
            } else {
                perr(&pe, "unknown material key '%s'", key);
            }
            break;

        case S_BOX:
            if (!strcmp(key, "x")) { parse_list(val, card->box_x, 3, &pe);
                                     have_bx = 1; }
            else if (!strcmp(key, "y")) { parse_list(val, card->box_y, 3, &pe);
                                          have_by = 1; }
            else if (!strcmp(key, "z")) { parse_list(val, card->box_z, 3, &pe);
                                          have_bz = 1; }
            else perr(&pe, "unknown box_orientation key '%s'", key);
            break;

        case S_DISLOS: {
            if (card->ndislo == 0) {
                perr(&pe, "key before the first '- ' list item", NULL);
                break;
            }
            CardDislo *d = &card->dislo[card->ndislo - 1];
            if (!strcmp(key, "point")) {
                parse_list(val, d->point, 3, &pe); dislo_keys |= DK_POINT;
            } else if (!strcmp(key, "line")) {
                parse_list(val, d->line, 3, &pe); dislo_keys |= DK_LINE;
            } else if (!strcmp(key, "burgers_direction")) {
                parse_list(val, d->b_dir, 3, &pe); dislo_keys |= DK_BDIR;
            } else if (!strcmp(key, "burgers_magnitude")) {
                d->b_mag = parse_scalar(val, &pe); dislo_keys |= DK_BMAG;
            } else if (!strcmp(key, "glide_normal")) {
                parse_list(val, d->glide_normal, 3, &pe);
                dislo_keys |= DK_GLIDE;
            } else if (!strcmp(key, "cut_angle")) {
                d->cut_angle = parse_scalar(val, &pe);
            } else {
                perr(&pe, "unknown dislocation key '%s'", key);
            }
            break;
        }

        default:
            perr(&pe, "'%s' outside any section", key);
        }
    }
    fclose(fp);

    /* Final completeness checks. */
    if (!pe.failed) {
        pe.lineno = 0;  /* file-level errors */
        if (sec == S_CIJ_FULL && full_rows < 6)
            perr(&pe, "'cij: full:' has fewer than 6 rows", NULL);
        else if (!have_cij)
            perr(&pe, "missing elastic constants (material: cij:)", NULL);
        else if (!have_a0 || card->a0 <= 0.0)
            perr(&pe, "missing or non-positive material: lattice_parameter",
                 NULL);
        else if (!have_bx || !have_by || !have_bz)
            perr(&pe, "box_orientation needs x, y and z", NULL);
        else if (card->ndislo == 0)
            perr(&pe, "no dislocations given", NULL);
        else if ((dislo_keys & (DK_POINT|DK_LINE|DK_BDIR|DK_BMAG|DK_GLIDE))
                     != (DK_POINT|DK_LINE|DK_BDIR|DK_BMAG|DK_GLIDE))
            perr(&pe, "last dislocation is missing a required key "
                      "(point, line, burgers_direction, burgers_magnitude, "
                      "glide_normal)", NULL);
    }

    /* Symmetry check on the stiffness matrix. */
    if (!pe.failed) {
        double cmax = 0.0;
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++)
                if (fabs(card->C6[i][j]) > cmax) cmax = fabs(card->C6[i][j]);
        if (cmax <= 0.0) {
            perr(&pe, "stiffness matrix is zero", NULL);
        } else {
            for (int i = 0; i < 6 && !pe.failed; i++)
                for (int j = i + 1; j < 6 && !pe.failed; j++)
                    if (fabs(card->C6[i][j] - card->C6[j][i]) > 1e-6 * cmax)
                        perr(&pe, "stiffness matrix is not symmetric", NULL);
        }
    }

    if (pe.failed) {
        if (pe.lineno > 0)
            snprintf(err, (size_t)errlen, "%s", pe.err);
        else
            snprintf(err, (size_t)errlen, "card: %s",
                     strchr(pe.err, ':') ? strchr(pe.err, ':') + 2 : pe.err);
        return 1;
    }
    return 0;
}
