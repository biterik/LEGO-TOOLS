/* atomio.c — shared I/O library for atomistic configurations.
 *
 * See atomio.h for API overview. This file implements readers/writers
 * for LAMMPS data files, LAMMPS dump files, and IMD ASCII checkpoint
 * files, with optional gzip support (zlib) and OpenMP-parallel parsing
 * of the per-atom rows.
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "atomio.h"

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* -------------------------------------------------------------------------- */
/* small helpers                                                              */
/* -------------------------------------------------------------------------- */

static void die(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "atomio: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

static char *xstrndup(const char *s, size_t n) {
    char *r = (char *)malloc(n + 1);
    memcpy(r, s, n);
    r[n] = 0;
    return r;
}

static int str_has_suffix(const char *s, const char *suf) {
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

static inline const char *skip_ws(const char *p, const char *end) {
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    return p;
}

/* -------------------------------------------------------------------------- */
/* Config book-keeping                                                        */
/* -------------------------------------------------------------------------- */

void config_init(Config *c) {
    memset(c, 0, sizeof(*c));
    c->format = FMT_UNKNOWN;
    c->pbc[0] = c->pbc[1] = c->pbc[2] = 1;
    c->idx_id = c->idx_type = c->idx_mol = c->idx_q = -1;
    c->idx_x  = c->idx_y  = c->idx_z  = -1;
    c->idx_vx = c->idx_vy = c->idx_vz = -1;
    c->idx_ix = c->idx_iy = c->idx_iz = -1;
    c->ntypes = 0;
    c->atom_style[0] = 0;
    c->imd_ndata = 0;
    c->imd_has_mass = 0;
}

void config_free(Config *c) {
    free(c->data);
    free(c->title);
    free(c->extra_sections);
    config_init(c);
}

int config_find_col(const Config *c, const char *name) {
    for (int i = 0; i < c->ncols; i++)
        if (strcmp(c->col_names[i], name) == 0) return i;
    return -1;
}

void config_refresh_std_indices(Config *c) {
    c->idx_id   = config_find_col(c, "id");
    c->idx_type = config_find_col(c, "type");
    c->idx_mol  = config_find_col(c, "mol");
    c->idx_q    = config_find_col(c, "q");
    c->idx_x    = config_find_col(c, "x");
    c->idx_y    = config_find_col(c, "y");
    c->idx_z    = config_find_col(c, "z");
    c->idx_vx   = config_find_col(c, "vx");
    c->idx_vy   = config_find_col(c, "vy");
    c->idx_vz   = config_find_col(c, "vz");
    c->idx_ix   = config_find_col(c, "ix");
    c->idx_iy   = config_find_col(c, "iy");
    c->idx_iz   = config_find_col(c, "iz");
    c->has_velocities =
        (c->idx_vx >= 0 && c->idx_vy >= 0 && c->idx_vz >= 0);
}

static int add_col(Config *c, const char *name, int is_int) {
    if (c->ncols >= MAX_COLS) die("too many columns (MAX_COLS=%d)", MAX_COLS);
    int idx = c->ncols++;
    /* Manual copy with explicit length clamp.  Using strncpy here trips
     * -Wstringop-truncation, and snprintf("%s", ...) trips
     * -Wformat-truncation, because callers may pass very long raw tokens.
     * Column names that don't fit in MAX_COL_NAME bytes are silently
     * truncated (they're internal lookup keys, not user output). */
    size_t n = strlen(name);
    if (n > (size_t)(MAX_COL_NAME - 1)) n = MAX_COL_NAME - 1;
    memcpy(c->col_names[idx], name, n);
    c->col_names[idx][n] = '\0';
    c->col_is_int[idx] = is_int;
    return idx;
}

/* -------------------------------------------------------------------------- */
/* format enum helpers                                                        */
/* -------------------------------------------------------------------------- */

const char *format_name(Format f) {
    switch (f) {
    case FMT_LAMMPS_DATA: return "lammps-data";
    case FMT_LAMMPS_DUMP: return "lammps-dump";
    case FMT_IMD:         return "imd";
    default:              return "unknown";
    }
}

Format format_from_string(const char *s) {
    if (!s) return FMT_UNKNOWN;
    if (!strcmp(s, "lammps") || !strcmp(s, "lammps-data") ||
        !strcmp(s, "data"))
        return FMT_LAMMPS_DATA;
    if (!strcmp(s, "lammps-dump") || !strcmp(s, "dump"))
        return FMT_LAMMPS_DUMP;
    if (!strcmp(s, "imd"))
        return FMT_IMD;
    return FMT_UNKNOWN;
}

/* -------------------------------------------------------------------------- */
/* file I/O (gzip-aware)                                                      */
/* -------------------------------------------------------------------------- */

/* Read an entire (possibly gzipped) file into a single NUL-terminated buffer. */
static int read_whole_file(const char *path, char **out, size_t *out_len) {
    gzFile gf = gzopen(path, "rb");
    if (!gf) {
        fprintf(stderr, "atomio: cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }
#if ZLIB_VERNUM >= 0x1235
    gzbuffer(gf, 1 << 20);
#endif
    size_t cap = 1 << 20, len = 0;
    char *buf = (char *)malloc(cap);
    if (!buf) { gzclose(gf); return -1; }
    for (;;) {
        if (len == cap) {
            cap *= 2;
            char *nb = (char *)realloc(buf, cap);
            if (!nb) { free(buf); gzclose(gf); return -1; }
            buf = nb;
        }
        int n = gzread(gf, buf + len, (unsigned)(cap - len));
        if (n < 0) {
            fprintf(stderr, "atomio: read error on %s\n", path);
            free(buf); gzclose(gf);
            return -1;
        }
        if (n == 0) break;
        len += (size_t)n;
    }
    gzclose(gf);
    char *nb = (char *)realloc(buf, len + 1);
    if (nb) buf = nb;
    buf[len] = 0;
    *out = buf;
    *out_len = len;
    return 0;
}

Format detect_format(const char *filename) {
    gzFile gf = gzopen(filename, "rb");
    if (!gf) return FMT_UNKNOWN;
    char head[512];
    int n = gzread(gf, head, (unsigned)sizeof(head) - 1);
    gzclose(gf);
    if (n <= 0) return FMT_UNKNOWN;
    head[n] = 0;
    const char *p = head;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    if (strncmp(p, "ITEM:", 5) == 0) return FMT_LAMMPS_DUMP;
    if (*p == '#' && (p[1] == 'F' || p[1] == 'C' || p[1] == 'X' || p[1] == 'E'))
        return FMT_IMD;
    /* Any file whose first non-blank line looks like LAMMPS data title. */
    return FMT_LAMMPS_DATA;
}

/* -------------------------------------------------------------------------- */
/* Output buffer with optional gzip backend                                   */
/* -------------------------------------------------------------------------- */

typedef struct {
    char   *buf;
    size_t  len;
    size_t  cap;
    gzFile  gf;
    FILE   *fp;
    int     use_gz;
} WBuf;

static int wbuf_open(WBuf *w, const char *path) {
    memset(w, 0, sizeof(*w));
    w->cap = 1 << 16;
    w->buf = (char *)malloc(w->cap);
    if (!w->buf) return -1;
    w->use_gz = str_has_suffix(path, ".gz");
    if (w->use_gz) {
        w->gf = gzopen(path, "wb");
        if (!w->gf) return -1;
#if ZLIB_VERNUM >= 0x1235
        gzbuffer(w->gf, 1 << 20);
#endif
    } else {
        w->fp = fopen(path, "wb");
        if (!w->fp) return -1;
        setvbuf(w->fp, NULL, _IOFBF, 1 << 20);
    }
    return 0;
}

static int wbuf_flush(WBuf *w) {
    if (w->len == 0) return 0;
    int rc = 0;
    if (w->use_gz) {
        int n = gzwrite(w->gf, w->buf, (unsigned)w->len);
        if (n != (int)w->len) rc = -1;
    } else {
        if (fwrite(w->buf, 1, w->len, w->fp) != w->len) rc = -1;
    }
    w->len = 0;
    return rc;
}

static int wbuf_close(WBuf *w) {
    int rc = wbuf_flush(w);
    if (w->use_gz && w->gf) gzclose(w->gf);
    if (!w->use_gz && w->fp) fclose(w->fp);
    free(w->buf);
    memset(w, 0, sizeof(*w));
    return rc;
}

static void wbuf_reserve(WBuf *w, size_t need) {
    if (w->len + need <= w->cap) return;
    if (w->cap < 1 << 20) {
        while (w->len + need > w->cap) w->cap *= 2;
        char *nb = (char *)realloc(w->buf, w->cap);
        if (!nb) die("oom");
        w->buf = nb;
    } else {
        wbuf_flush(w);
        if (need > w->cap) {
            w->cap = need * 2;
            char *nb = (char *)realloc(w->buf, w->cap);
            if (!nb) die("oom");
            w->buf = nb;
        }
    }
}

static void wbuf_write(WBuf *w, const char *s, size_t n) {
    wbuf_reserve(w, n);
    memcpy(w->buf + w->len, s, n);
    w->len += n;
}

static void wbuf_puts(WBuf *w, const char *s) {
    wbuf_write(w, s, strlen(s));
}

static void wbuf_printf(WBuf *w, const char *fmt, ...) {
    va_list ap;
    char tmp[512];
    va_start(ap, fmt);
    int n = vsnprintf(tmp, sizeof(tmp), fmt, ap);
    va_end(ap);
    if (n < 0) return;
    if ((size_t)n < sizeof(tmp)) {
        wbuf_write(w, tmp, (size_t)n);
        return;
    }
    wbuf_reserve(w, (size_t)n + 1);
    va_start(ap, fmt);
    vsnprintf(w->buf + w->len, (size_t)n + 1, fmt, ap);
    va_end(ap);
    w->len += (size_t)n;
}

/* -------------------------------------------------------------------------- */
/* LAMMPS data reader                                                         */
/* -------------------------------------------------------------------------- */

static const char *lmps_sections[] = {
    "Atoms", "Velocities", "Masses", "Ellipsoids", "Lines", "Triangles",
    "Bodies", "Bonds", "Angles", "Dihedrals", "Impropers",
    "Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs",
    "Dihedral Coeffs", "Improper Coeffs",
    "BondBond Coeffs", "BondAngle Coeffs", "MiddleBondTorsion Coeffs",
    "EndBondTorsion Coeffs", "AngleTorsion Coeffs",
    "AngleAngleTorsion Coeffs", "BondBond13 Coeffs", "AngleAngle Coeffs",
    "CMAP", "AtomFile",
    NULL
};

static int line_starts_section(const char *line, size_t len, char *out_name) {
    const char *p = line, *e = line + len;
    while (p < e && (*p == ' ' || *p == '\t')) p++;
    if (p >= e || !isupper((unsigned char)*p)) return 0;
    int best = -1, best_len = 0;
    for (int i = 0; lmps_sections[i]; i++) {
        int slen = (int)strlen(lmps_sections[i]);
        if (p + slen > e) continue;
        if (strncmp(p, lmps_sections[i], (size_t)slen) != 0) continue;
        char next = (p + slen < e) ? p[slen] : '\n';
        if (next == 0 || next == '\n' || next == '\r' || next == ' ' ||
            next == '\t' || next == '#') {
            if (slen > best_len) { best_len = slen; best = i; }
        }
    }
    if (best < 0) return 0;
    strcpy(out_name, lmps_sections[best]);
    return 1;
}

/* Describe a LAMMPS atom style as a sequence of columns.
 * Format in the string: space-separated tokens. Integer tokens are
 * prefixed with '@'. All styles get image flags "ix iy iz" appended
 * optionally (detected at read time).
 */
static const struct {
    const char *style;
    const char *cols;
} atom_style_cols[] = {
    {"atomic",    "@id @type x y z"},
    {"charge",    "@id @type q x y z"},
    {"molecular", "@id @mol @type x y z"},
    {"full",      "@id @mol @type q x y z"},
    {"bond",      "@id @mol @type x y z"},
    {"angle",     "@id @mol @type x y z"},
    {"sphere",    "@id @type diameter density x y z"},
    {"dipole",    "@id @type q x y z mux muy muz"},
    {"peri",      "@id @type volume density x y z"},
    {"meso",      "@id @type rho e cv x y z"},
    {NULL, NULL}
};

static const char *lookup_atom_style_cols(const char *style) {
    for (int i = 0; atom_style_cols[i].style; i++)
        if (strcmp(style, atom_style_cols[i].style) == 0)
            return atom_style_cols[i].cols;
    return NULL;
}

/* Count whitespace-separated tokens on a line. */
static int count_tokens(const char *p, const char *e) {
    int n = 0;
    while (p < e && *p != '#') {
        while (p < e && (*p == ' ' || *p == '\t')) p++;
        if (p >= e || *p == '#' || *p == '\n' || *p == '\r') break;
        n++;
        while (p < e && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' &&
               *p != '#')
            p++;
    }
    return n;
}

/* Locate the first non-blank, non-comment line in a section body. */
static const char *first_data_line(const char *p, const char *e,
                                   size_t *out_len) {
    while (p < e) {
        const char *l = p;
        const char *lend = memchr(p, '\n', (size_t)(e - p));
        if (!lend) lend = e;
        const char *q = skip_ws(l, lend);
        if (q < lend && *q != '#' && *q != '\r') {
            *out_len = (size_t)(lend - l);
            return l;
        }
        p = (lend < e) ? lend + 1 : e;
    }
    return NULL;
}

/* Parse a single atom row given a pre-built token→column map. */
static inline void parse_row_tokens(double *row, const int *col_map,
                                    int ntok_expected, const char *line,
                                    const char *end) {
    const char *p = line;
    for (int t = 0; t < ntok_expected; t++) {
        while (p < end && (*p == ' ' || *p == '\t')) p++;
        if (p >= end || *p == '\n' || *p == '\r' || *p == '#') break;
        char *nx;
        double v = strtod(p, &nx);
        if (nx == p) {
            /* non-numeric token; skip until whitespace */
            while (p < end && *p != ' ' && *p != '\t' && *p != '\n' &&
                   *p != '\r')
                p++;
        } else {
            p = nx;
        }
        int tgt = col_map[t];
        if (tgt >= 0) row[tgt] = v;
    }
}

/* Parse an entire section of rows in parallel. */
static void parse_rows_par(double *data, int ncols_out,
                           const int *col_map, int ntok_expected,
                           const char *text, size_t len, size_t nrows) {
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
    if (nth > 64) nth = 64;
#endif
    if (nrows < 10000 || nth < 2) {
        const char *p = text, *e = text + len;
        size_t row = 0;
        while (p < e && row < nrows) {
            const char *l = p;
            const char *lend = memchr(p, '\n', (size_t)(e - p));
            if (!lend) lend = e;
            const char *q = skip_ws(l, lend);
            if (q < lend && *q != '#' && *q != '\r') {
                parse_row_tokens(data + row * ncols_out, col_map, ntok_expected,
                                 q, lend);
                row++;
            }
            p = (lend < e) ? lend + 1 : e;
        }
        return;
    }

    /* Partition text into roughly equal byte ranges, snap to line starts. */
    size_t *starts  = (size_t *)calloc((size_t)nth + 1, sizeof(size_t));
    size_t *counts  = (size_t *)calloc((size_t)nth + 1, sizeof(size_t));
    size_t *offsets = (size_t *)calloc((size_t)nth + 1, sizeof(size_t));

    for (int t = 0; t < nth; t++) {
        size_t s = (len * (size_t)t) / (size_t)nth;
        if (t > 0) {
            while (s < len && text[s - 1] != '\n') s++;
        }
        starts[t] = s;
    }
    starts[nth] = len;

#pragma omp parallel for
    for (int t = 0; t < nth; t++) {
        const char *p = text + starts[t];
        const char *e = text + starts[t + 1];
        size_t n = 0;
        while (p < e) {
            const char *lend = memchr(p, '\n', (size_t)(e - p));
            if (!lend) lend = e;
            const char *q = skip_ws(p, lend);
            if (q < lend && *q != '#' && *q != '\r') n++;
            p = (lend < e) ? lend + 1 : e;
        }
        counts[t] = n;
    }

    offsets[0] = 0;
    for (int t = 0; t < nth; t++) offsets[t + 1] = offsets[t] + counts[t];

    if (offsets[nth] != nrows) {
        fprintf(stderr,
                "atomio: warning: parsed %zu rows, expected %zu\n",
                offsets[nth], nrows);
    }

#pragma omp parallel for
    for (int t = 0; t < nth; t++) {
        const char *p = text + starts[t];
        const char *e = text + starts[t + 1];
        size_t row = offsets[t];
        while (p < e && row < nrows) {
            const char *lend = memchr(p, '\n', (size_t)(e - p));
            if (!lend) lend = e;
            const char *q = skip_ws(p, lend);
            if (q < lend && *q != '#' && *q != '\r') {
                parse_row_tokens(data + row * (size_t)ncols_out, col_map,
                                 ntok_expected, q, lend);
                row++;
            }
            p = (lend < e) ? lend + 1 : e;
        }
    }

    free(starts);
    free(counts);
    free(offsets);
}

/* -------------------------------------------------------------------------- */
/* LAMMPS data reader (main)                                                  */
/* -------------------------------------------------------------------------- */

typedef struct {
    char   name[64];
    size_t start; /* offset in buf where section data begins */
    size_t end;   /* offset where section data ends */
    char   atom_style[32]; /* for Atoms section */
    size_t header_start;   /* offset of the section header line */
} SectionRef;

static int read_lammps_data(Config *c, char *buf, size_t buf_len) {
    c->format = FMT_LAMMPS_DATA;

    char *p = buf, *e = buf + buf_len;

    /* Title line */
    char *le = memchr(p, '\n', (size_t)(e - p));
    if (!le) le = e;
    c->title = xstrndup(p, (size_t)(le - p));
    /* strip trailing \r */
    size_t tl = strlen(c->title);
    if (tl && c->title[tl - 1] == '\r') c->title[tl - 1] = 0;
    p = (le < e) ? le + 1 : e;

    SectionRef secs[64];
    int nsec = 0;
    int past_header = 0;

    while (p < e) {
        char *lstart = p;
        le = memchr(p, '\n', (size_t)(e - p));
        if (!le) le = e;

        char name[64];
        if (line_starts_section(lstart, (size_t)(le - lstart), name)) {
            if (nsec > 0) secs[nsec - 1].end = (size_t)(lstart - buf);
            if (nsec >= 64) die("too many LAMMPS data sections");
            memset(&secs[nsec], 0, sizeof(secs[nsec]));
            strcpy(secs[nsec].name, name);
            secs[nsec].header_start = (size_t)(lstart - buf);
            secs[nsec].start = (size_t)((le < e ? le + 1 : e) - buf);
            secs[nsec].end   = (size_t)(e - buf);
            /* Extract atom_style comment for Atoms */
            if (strcmp(name, "Atoms") == 0) {
                char *hash = memchr(lstart, '#', (size_t)(le - lstart));
                if (hash) {
                    hash++;
                    while (hash < le && (*hash == ' ' || *hash == '\t')) hash++;
                    int j = 0;
                    while (hash < le && *hash != ' ' && *hash != '\t' &&
                           *hash != '\r' && *hash != '\n' && j < 31)
                        secs[nsec].atom_style[j++] = *hash++;
                    secs[nsec].atom_style[j] = 0;
                }
            }
            nsec++;
            past_header = 1;
            p = (le < e) ? le + 1 : e;
            continue;
        }

        if (!past_header) {
            /* Parse header line: counts and box. */
            char tmp[256];
            size_t ll = (size_t)(le - lstart);
            if (ll >= sizeof(tmp)) ll = sizeof(tmp) - 1;
            memcpy(tmp, lstart, ll);
            tmp[ll] = 0;
            char *tr = tmp;
            while (*tr == ' ' || *tr == '\t') tr++;
            if (*tr && *tr != '#' && *tr != '\r' && *tr != '\n') {
                long long ln = 0;
                int ni = 0;
                double v1, v2, v3;
                int matched = -1;
                if (sscanf(tr, "%lld atoms%n", &ln, &matched) >= 1 &&
                    matched > 0) {
                    c->natoms = (size_t)ln;
                } else if ((matched = -1,
                            sscanf(tr, "%d atom types%n", &ni, &matched) >= 1)
                           && matched > 0) {
                    c->ntypes = ni;
                } else if ((matched = -1,
                            sscanf(tr, "%lf %lf xlo xhi%n", &v1, &v2, &matched)
                            >= 2) && matched > 0) {
                    c->xlo = v1; c->xhi = v2;
                } else if ((matched = -1,
                            sscanf(tr, "%lf %lf ylo yhi%n", &v1, &v2, &matched)
                            >= 2) && matched > 0) {
                    c->ylo = v1; c->yhi = v2;
                } else if ((matched = -1,
                            sscanf(tr, "%lf %lf zlo zhi%n", &v1, &v2, &matched)
                            >= 2) && matched > 0) {
                    c->zlo = v1; c->zhi = v2;
                } else if ((matched = -1,
                            sscanf(tr, "%lf %lf %lf xy xz yz%n",
                                   &v1, &v2, &v3, &matched) >= 3) &&
                           matched > 0) {
                    c->xy = v1; c->xz = v2; c->yz = v3;
                    c->triclinic = 1;
                }
            }
        }
        p = (le < e) ? le + 1 : e;
    }

    /* Find sections of interest. */
    int iatoms = -1, ivel = -1;
    for (int i = 0; i < nsec; i++) {
        if (!strcmp(secs[i].name, "Atoms"))      iatoms = i;
        if (!strcmp(secs[i].name, "Velocities")) ivel   = i;
    }
    if (iatoms < 0) { fprintf(stderr, "atomio: no Atoms section\n"); return -1; }

    const char *as = buf + secs[iatoms].start;
    const char *ae = buf + secs[iatoms].end;

    /* atom_style */
    if (secs[iatoms].atom_style[0])
        strcpy(c->atom_style, secs[iatoms].atom_style);
    else
        strcpy(c->atom_style, "atomic");

    /* Figure out number of tokens per Atoms line. */
    size_t first_len = 0;
    const char *first = first_data_line(as, ae, &first_len);
    if (!first) { fprintf(stderr, "atomio: empty Atoms section\n"); return -1; }
    int ntok = count_tokens(first, first + first_len);

    /* Build column layout from atom_style. */
    const char *spec = lookup_atom_style_cols(c->atom_style);
    int base_ncols = 0;
    int base_is_int[MAX_COLS] = {0};
    char base_names[MAX_COLS][MAX_COL_NAME];
    if (!spec) {
        /* Unknown style → fall back to: id type x y z extra... */
        strcpy(base_names[0], "id");    base_is_int[0] = 1;
        strcpy(base_names[1], "type");  base_is_int[1] = 1;
        strcpy(base_names[2], "x");
        strcpy(base_names[3], "y");
        strcpy(base_names[4], "z");
        base_ncols = 5;
    } else {
        const char *q = spec;
        while (*q) {
            while (*q == ' ') q++;
            if (!*q) break;
            int is_int = 0;
            if (*q == '@') { is_int = 1; q++; }
            const char *start = q;
            while (*q && *q != ' ') q++;
            int len = (int)(q - start);
            if (len >= MAX_COL_NAME) len = MAX_COL_NAME - 1;
            memcpy(base_names[base_ncols], start, (size_t)len);
            base_names[base_ncols][len] = 0;
            base_is_int[base_ncols] = is_int;
            base_ncols++;
        }
    }

    /* Detect optional image flags (3 extra integer tokens). */
    int has_image = 0;
    int extra_tokens = ntok - base_ncols;
    if (extra_tokens == 3) has_image = 1;
    /* Any other extras → generic "colN" passthrough. */

    /* Register columns in Config. */
    for (int i = 0; i < base_ncols; i++)
        add_col(c, base_names[i], base_is_int[i]);
    if (has_image) {
        add_col(c, "ix", 1);
        add_col(c, "iy", 1);
        add_col(c, "iz", 1);
    } else {
        for (int k = 0; k < extra_tokens; k++) {
            char nm[32];
            snprintf(nm, sizeof(nm), "col%d", base_ncols + k + 1);
            add_col(c, nm, 0);
        }
    }
    config_refresh_std_indices(c);

    /* Build token→column map for the Atoms section. */
    int col_map[MAX_COLS];
    for (int i = 0; i < ntok && i < MAX_COLS; i++) col_map[i] = i;

    /* Allocate data. */
    c->data = (double *)calloc(c->natoms * (size_t)c->ncols, sizeof(double));
    if (!c->data) die("oom allocating atoms (%zu * %d)", c->natoms, c->ncols);

    parse_rows_par(c->data, c->ncols, col_map, ntok,
                   as, (size_t)(ae - as), c->natoms);

    /* Velocities section (optional). */
    if (ivel >= 0) {
        const char *vs = buf + secs[ivel].start;
        const char *ve_ = buf + secs[ivel].end;
        size_t vfl = 0;
        const char *vf = first_data_line(vs, ve_, &vfl);
        if (vf) {
            /* Add vx/vy/vz columns. */
            int ivx = add_col(c, "vx", 0);
            int ivy = add_col(c, "vy", 0);
            int ivz = add_col(c, "vz", 0);

            /* Grow data: insert 3 new cols per row. */
            int new_ncols = c->ncols;
            int old_ncols = new_ncols - 3;
            double *nd = (double *)malloc(c->natoms * (size_t)new_ncols *
                                           sizeof(double));
            if (!nd) die("oom velocities");
            for (size_t i = 0; i < c->natoms; i++) {
                memcpy(nd + i * (size_t)new_ncols,
                       c->data + i * (size_t)old_ncols,
                       (size_t)old_ncols * sizeof(double));
                nd[i * (size_t)new_ncols + ivx] = 0.0;
                nd[i * (size_t)new_ncols + ivy] = 0.0;
                nd[i * (size_t)new_ncols + ivz] = 0.0;
            }
            free(c->data);
            c->data = nd;
            config_refresh_std_indices(c);

            /* Velocity section rows: id vx vy vz. Assume same order as Atoms.
             * Verify via id match at the end. */
            int vel_ntok = count_tokens(vf, vf + vfl);
            if (vel_ntok < 4) {
                fprintf(stderr,
                        "atomio: warning: Velocities section has %d tokens\n",
                        vel_ntok);
            }

            /* Parse into a temporary [id, vx, vy, vz] array. */
            double *tmp =
                (double *)malloc(c->natoms * 4 * sizeof(double));
            int tmap[64];
            for (int i = 0; i < vel_ntok && i < 64; i++) tmap[i] = i < 4 ? i : -1;
            parse_rows_par(tmp, 4, tmap, vel_ntok,
                           vs, (size_t)(ve_ - vs), c->natoms);

            /* Attempt to assign in id-order using a direct id→row map. */
            long long min_id = (long long)(tmp[0]);
            long long max_id = min_id;
            for (size_t i = 1; i < c->natoms; i++) {
                long long v = (long long)(tmp[i * 4]);
                if (v < min_id) min_id = v;
                if (v > max_id) max_id = v;
            }
            size_t range = (size_t)(max_id - min_id + 1);

            size_t *id2row = NULL;
            if (c->idx_id >= 0 && range <= (size_t)(c->natoms) * 4 + 1024 &&
                range > 0) {
                id2row = (size_t *)malloc(range * sizeof(size_t));
                for (size_t i = 0; i < range; i++) id2row[i] = SIZE_MAX;
                for (size_t i = 0; i < c->natoms; i++) {
                    long long v =
                        (long long)(c->data[i * (size_t)c->ncols + c->idx_id]);
                    if (v - min_id >= 0 && (size_t)(v - min_id) < range)
                        id2row[v - min_id] = i;
                }
            }

            for (size_t i = 0; i < c->natoms; i++) {
                long long vid = (long long)(tmp[i * 4]);
                size_t row;
                if (id2row && vid - min_id >= 0 &&
                    (size_t)(vid - min_id) < range &&
                    id2row[vid - min_id] != SIZE_MAX) {
                    row = id2row[vid - min_id];
                } else {
                    row = i; /* fall back to same order */
                }
                c->data[row * (size_t)c->ncols + ivx] = tmp[i * 4 + 1];
                c->data[row * (size_t)c->ncols + ivy] = tmp[i * 4 + 2];
                c->data[row * (size_t)c->ncols + ivz] = tmp[i * 4 + 3];
            }
            free(id2row);
            free(tmp);
        }
    }

    /* Capture extra pass-through sections as raw text. */
    size_t extra_cap = 0;
    for (int i = 0; i < nsec; i++) {
        if (i == iatoms || i == ivel) continue;
        extra_cap += secs[i].end - secs[i].header_start + 2;
    }
    if (extra_cap) {
        c->extra_sections = (char *)malloc(extra_cap + 1);
        c->extra_sections[0] = 0;
        size_t off = 0;
        for (int i = 0; i < nsec; i++) {
            if (i == iatoms || i == ivel) continue;
            size_t n = secs[i].end - secs[i].header_start;
            memcpy(c->extra_sections + off, buf + secs[i].header_start, n);
            off += n;
            if (off && c->extra_sections[off - 1] != '\n')
                c->extra_sections[off++] = '\n';
        }
        c->extra_sections[off] = 0;
    }

    config_refresh_std_indices(c);
    return 0;
}

/* -------------------------------------------------------------------------- */
/* LAMMPS dump reader                                                         */
/* -------------------------------------------------------------------------- */

static int read_lammps_dump(Config *c, char *buf, size_t buf_len) {
    c->format = FMT_LAMMPS_DUMP;
    char *p = buf, *e = buf + buf_len;

    /* Walk line by line; handle ITEM: sections. */
    int have_natoms = 0;
    int have_box    = 0;
    int tri_box     = 0;
    int in_atoms    = 0;
    const char *atoms_start = NULL;
    int col_map[MAX_COLS];
    int ntok = 0;

    while (p < e) {
        char *le = memchr(p, '\n', (size_t)(e - p));
        if (!le) le = e;
        size_t ll = (size_t)(le - p);

        if (ll >= 5 && strncmp(p, "ITEM:", 5) == 0) {
            const char *tag = p + 5;
            while (tag < le && (*tag == ' ' || *tag == '\t')) tag++;
            if (!strncmp(tag, "TIMESTEP", 8)) {
                /* next line is timestep, skip */
                p = (le < e) ? le + 1 : e;
                le = memchr(p, '\n', (size_t)(e - p));
                if (!le) le = e;
            } else if (!strncmp(tag, "NUMBER OF ATOMS", 15)) {
                p = (le < e) ? le + 1 : e;
                le = memchr(p, '\n', (size_t)(e - p));
                if (!le) le = e;
                c->natoms = (size_t)strtoull(p, NULL, 10);
                have_natoms = 1;
            } else if (!strncmp(tag, "BOX BOUNDS", 10)) {
                if (strstr(tag, "xy")) tri_box = 1;
                for (int axis = 0; axis < 3; axis++) {
                    p = (le < e) ? le + 1 : e;
                    le = memchr(p, '\n', (size_t)(e - p));
                    if (!le) le = e;
                    double lo = 0, hi = 0, ti = 0;
                    if (tri_box) {
                        sscanf(p, "%lf %lf %lf", &lo, &hi, &ti);
                    } else {
                        sscanf(p, "%lf %lf", &lo, &hi);
                    }
                    if (axis == 0) { c->xlo = lo; c->xhi = hi; if (tri_box) c->xy = ti; }
                    if (axis == 1) { c->ylo = lo; c->yhi = hi; if (tri_box) c->xz = ti; }
                    if (axis == 2) { c->zlo = lo; c->zhi = hi; if (tri_box) c->yz = ti; }
                }
                if (tri_box) c->triclinic = 1;
                have_box = 1;
            } else if (!strncmp(tag, "ATOMS", 5)) {
                /* Column names follow on this same line. */
                const char *q = tag + 5;
                while (q < le && (*q == ' ' || *q == '\t')) q++;
                while (q < le) {
                    while (q < le && (*q == ' ' || *q == '\t')) q++;
                    if (q >= le || *q == '\r' || *q == '\n') break;
                    const char *nm = q;
                    while (q < le && *q != ' ' && *q != '\t' && *q != '\r' &&
                           *q != '\n') q++;
                    int nlen = (int)(q - nm);
                    if (nlen >= MAX_COL_NAME) nlen = MAX_COL_NAME - 1;
                    char name[MAX_COL_NAME] = {0};
                    memcpy(name, nm, (size_t)nlen);
                    /* Map dump aliases. */
                    int is_int = 0;
                    if (!strcmp(name, "id")   || !strcmp(name, "type") ||
                        !strcmp(name, "mol")  || !strcmp(name, "ix")   ||
                        !strcmp(name, "iy")   || !strcmp(name, "iz")   ||
                        !strcmp(name, "proc") || !strcmp(name, "procp1"))
                        is_int = 1;
                    /* unwrap variants → primary names */
                    if (!strcmp(name, "xs")  || !strcmp(name, "xu") ||
                        !strcmp(name, "xsu")) strcpy(name, "x");
                    if (!strcmp(name, "ys")  || !strcmp(name, "yu") ||
                        !strcmp(name, "ysu")) strcpy(name, "y");
                    if (!strcmp(name, "zs")  || !strcmp(name, "zu") ||
                        !strcmp(name, "zsu")) strcpy(name, "z");
                    add_col(c, name, is_int);
                }
                ntok = c->ncols;
                for (int i = 0; i < ntok; i++) col_map[i] = i;
                config_refresh_std_indices(c);

                /* Data lines start after this. */
                p = (le < e) ? le + 1 : e;
                in_atoms = 1;
                atoms_start = p;
                /* Seek to end of atoms (either EOF or next ITEM:) */
                const char *scan = atoms_start;
                while (scan < e) {
                    const char *nl = memchr(scan, '\n', (size_t)(e - scan));
                    if (!nl) nl = e;
                    if (nl - scan >= 5 && !strncmp(scan, "ITEM:", 5)) break;
                    scan = (nl < e) ? nl + 1 : e;
                }
                const char *atoms_end = scan;

                if (!have_natoms) {
                    fprintf(stderr, "atomio: dump has no NUMBER OF ATOMS\n");
                    return -1;
                }
                c->data = (double *)calloc(c->natoms * (size_t)c->ncols,
                                            sizeof(double));
                if (!c->data) die("oom atoms");
                parse_rows_par(c->data, c->ncols, col_map, ntok,
                               atoms_start, (size_t)(atoms_end - atoms_start),
                               c->natoms);
                p = (char *)atoms_end;
                in_atoms = 0;
                le = memchr(p, '\n', (size_t)(e - p));
                if (!le) le = e;
                continue;
            }
        }
        p = (le < e) ? le + 1 : e;
    }

    if (!have_box) {
        fprintf(stderr, "atomio: dump has no BOX BOUNDS\n");
        return -1;
    }
    (void)in_atoms;
    return 0;
}

/* -------------------------------------------------------------------------- */
/* IMD ASCII reader                                                           */
/* -------------------------------------------------------------------------- */

static int read_imd(Config *c, char *buf, size_t buf_len) {
    c->format = FMT_IMD;
    char *p = buf, *e = buf + buf_len;

    int n_id = 0, n_type = 0, n_mass = 0, n_pos = 0, n_vel = 0, n_data = 0;
    int parsed_F = 0;
    char col_names[MAX_COLS][MAX_COL_NAME];
    int  ncols_decl = 0;

    while (p < e) {
        char *le = memchr(p, '\n', (size_t)(e - p));
        if (!le) le = e;
        if (*p != '#') break;
        size_t ll = (size_t)(le - p);
        char line[512];
        if (ll >= sizeof(line)) ll = sizeof(line) - 1;
        memcpy(line, p, ll);
        line[ll] = 0;

        if (line[1] == 'F') {
            /* #F A n_id n_type n_mass n_pos n_vel n_data */
            char fmt = 0;
            int rc = sscanf(line + 2, " %c %d %d %d %d %d %d",
                            &fmt, &n_id, &n_type, &n_mass, &n_pos, &n_vel,
                            &n_data);
            if (rc < 5 || (fmt != 'A' && fmt != 'a')) {
                fprintf(stderr, "atomio: only ASCII IMD (#F A ...) supported\n");
                return -1;
            }
            if (rc < 7) n_data = 0;
            if (rc < 6) n_vel = 0;
            parsed_F = 1;
        } else if (line[1] == 'C') {
            /* Column names list. Everything after "#C". */
            const char *q = line + 2;
            while (*q == ' ' || *q == '\t') q++;
            while (*q) {
                const char *nm = q;
                while (*q && *q != ' ' && *q != '\t') q++;
                int nlen = (int)(q - nm);
                if (nlen > 0 && ncols_decl < MAX_COLS) {
                    if (nlen >= MAX_COL_NAME) nlen = MAX_COL_NAME - 1;
                    memcpy(col_names[ncols_decl], nm, (size_t)nlen);
                    col_names[ncols_decl][nlen] = 0;
                    ncols_decl++;
                }
                while (*q == ' ' || *q == '\t') q++;
            }
        } else if (line[1] == 'X') {
            double a, b, d;
            if (sscanf(line + 2, " %lf %lf %lf", &a, &b, &d) >= 1) {
                c->xlo = 0; c->xhi = a;
                c->xy = b; c->xz = d;
            }
        } else if (line[1] == 'Y') {
            double a, b, d;
            if (sscanf(line + 2, " %lf %lf %lf", &a, &b, &d) >= 2) {
                c->ylo = 0; c->yhi = b;
            }
        } else if (line[1] == 'Z') {
            double a, b, d;
            if (sscanf(line + 2, " %lf %lf %lf", &a, &b, &d) >= 3) {
                c->zlo = 0; c->zhi = d;
            }
        } else if (line[1] == 'E') {
            /* End of header. */
            p = (le < e) ? le + 1 : e;
            break;
        }
        p = (le < e) ? le + 1 : e;
    }

    if (!parsed_F) {
        fprintf(stderr, "atomio: no #F line in IMD file\n");
        return -1;
    }

    c->imd_has_mass = (n_mass > 0);
    c->imd_ndata    = n_data;

    /* Build columns in declared order. */
    if (ncols_decl > 0) {
        for (int i = 0; i < ncols_decl; i++) {
            int is_int =
                (!strcmp(col_names[i], "number") ||
                 !strcmp(col_names[i], "id")     ||
                 !strcmp(col_names[i], "type"));
            /* IMD "number" → rename to "id" for uniformity */
            if (!strcmp(col_names[i], "number")) strcpy(col_names[i], "id");
            add_col(c, col_names[i], is_int);
        }
    } else {
        /* Fall back to derived names based on counts. */
        if (n_id)   add_col(c, "id",   1);
        if (n_type) add_col(c, "type", 1);
        if (n_mass) add_col(c, "mass", 0);
        if (n_pos >= 1) add_col(c, "x", 0);
        if (n_pos >= 2) add_col(c, "y", 0);
        if (n_pos >= 3) add_col(c, "z", 0);
        if (n_vel >= 1) add_col(c, "vx", 0);
        if (n_vel >= 2) add_col(c, "vy", 0);
        if (n_vel >= 3) add_col(c, "vz", 0);
        for (int i = 0; i < n_data; i++) {
            char nm[32];
            snprintf(nm, sizeof(nm), "data%d", i + 1);
            add_col(c, nm, 0);
        }
    }

    int ntok = n_id + n_type + n_mass + n_pos + n_vel + n_data;
    if (c->ncols != ntok && ncols_decl > 0) {
        /* Trust #C if it disagrees with #F totals. */
        ntok = c->ncols;
    }
    config_refresh_std_indices(c);

    /* Count atoms in the body. */
    size_t body_off = (size_t)(p - buf);
    size_t nat = 0;
    {
        const char *q = p;
        while (q < e) {
            const char *nl = memchr(q, '\n', (size_t)(e - q));
            if (!nl) nl = e;
            const char *r = skip_ws(q, nl);
            if (r < nl && *r != '#' && *r != '\r') nat++;
            q = (nl < e) ? nl + 1 : e;
        }
    }
    c->natoms = nat;
    c->data = (double *)calloc(nat * (size_t)c->ncols, sizeof(double));
    if (!c->data) die("oom imd atoms");

    int col_map[MAX_COLS];
    for (int i = 0; i < ntok && i < MAX_COLS; i++) col_map[i] = i;
    parse_rows_par(c->data, c->ncols, col_map, ntok,
                   buf + body_off, buf_len - body_off, nat);

    /* Determine ntypes from type column if present. */
    if (c->idx_type >= 0) {
        int maxt = 0;
        for (size_t i = 0; i < c->natoms; i++) {
            int t = (int)c->data[i * (size_t)c->ncols + c->idx_type];
            if (t > maxt) maxt = t;
        }
        c->ntypes = maxt;
    }
    return 0;
}

/* -------------------------------------------------------------------------- */
/* Top-level dispatch                                                         */
/* -------------------------------------------------------------------------- */

int config_read(Config *c, const char *filename) {
    config_init(c);
    char *buf = NULL;
    size_t len = 0;
    if (read_whole_file(filename, &buf, &len) != 0) return -1;
    Format f = detect_format(filename);
    int rc = -1;
    switch (f) {
    case FMT_LAMMPS_DATA: rc = read_lammps_data(c, buf, len); break;
    case FMT_LAMMPS_DUMP: rc = read_lammps_dump(c, buf, len); break;
    case FMT_IMD:         rc = read_imd(c, buf, len);         break;
    default:
        fprintf(stderr, "atomio: unknown format for %s\n", filename);
        rc = -1;
    }
    free(buf);
    return rc;
}

/* -------------------------------------------------------------------------- */
/* LAMMPS data writer                                                         */
/* -------------------------------------------------------------------------- */

/* Pick an atom style for output based on which columns we have. */
static const char *pick_output_style(const Config *c) {
    if (c->atom_style[0] && lookup_atom_style_cols(c->atom_style))
        return c->atom_style;
    if (c->idx_mol >= 0 && c->idx_q >= 0) return "full";
    if (c->idx_mol >= 0)                  return "molecular";
    if (c->idx_q >= 0)                    return "charge";
    return "atomic";
}

static void fmt_double(char *out, size_t out_sz, double v) {
    snprintf(out, out_sz, "%.16g", v);
}

static int write_lammps_data(const Config *c, const char *path) {
    WBuf w;
    if (wbuf_open(&w, path) != 0) {
        fprintf(stderr, "atomio: cannot write %s\n", path);
        return -1;
    }

    const char *title = c->title ? c->title : "LAMMPS data file";
    wbuf_printf(&w, "%s\n\n", title);
    wbuf_printf(&w, "%zu atoms\n", c->natoms);
    if (c->ntypes > 0) wbuf_printf(&w, "%d atom types\n", c->ntypes);
    wbuf_puts(&w, "\n");
    wbuf_printf(&w, "%.16g %.16g xlo xhi\n", c->xlo, c->xhi);
    wbuf_printf(&w, "%.16g %.16g ylo yhi\n", c->ylo, c->yhi);
    wbuf_printf(&w, "%.16g %.16g zlo zhi\n", c->zlo, c->zhi);
    if (c->triclinic)
        wbuf_printf(&w, "%.16g %.16g %.16g xy xz yz\n", c->xy, c->xz, c->yz);
    wbuf_puts(&w, "\n");

    /* Emit pass-through sections (Masses, Bonds, ...). They already contain
     * their own section headers and blank lines. */
    if (c->extra_sections && c->extra_sections[0])
        wbuf_puts(&w, c->extra_sections);

    /* Atoms section. */
    const char *style = pick_output_style(c);
    wbuf_printf(&w, "\nAtoms # %s\n\n", style);

    /* Resolve output column order from style. */
    const char *spec = lookup_atom_style_cols(style);
    int out_idx[MAX_COLS];
    int out_is_int[MAX_COLS];
    int out_n = 0;
    if (spec) {
        const char *q = spec;
        while (*q) {
            while (*q == ' ') q++;
            if (!*q) break;
            int is_int = 0;
            if (*q == '@') { is_int = 1; q++; }
            const char *s = q;
            while (*q && *q != ' ') q++;
            char name[MAX_COL_NAME] = {0};
            int nlen = (int)(q - s);
            if (nlen >= MAX_COL_NAME) nlen = MAX_COL_NAME - 1;
            memcpy(name, s, (size_t)nlen);
            int idx = config_find_col(c, name);
            if (idx < 0) {
                /* Missing column: insert 0. Keep as-is, we signal via -1. */
                out_idx[out_n] = -1;
            } else {
                out_idx[out_n] = idx;
            }
            out_is_int[out_n] = is_int;
            out_n++;
        }
    }

    int write_image =
        (c->idx_ix >= 0 && c->idx_iy >= 0 && c->idx_iz >= 0);

    size_t N = c->natoms;
    int nc = c->ncols;
    const double *D = c->data;
    char tmp[64];

    for (size_t i = 0; i < N; i++) {
        for (int k = 0; k < out_n; k++) {
            if (k) wbuf_puts(&w, " ");
            double v = out_idx[k] >= 0 ? D[i * (size_t)nc + out_idx[k]] : 0.0;
            if (out_is_int[k]) {
                wbuf_printf(&w, "%lld", (long long)v);
            } else {
                fmt_double(tmp, sizeof(tmp), v);
                wbuf_puts(&w, tmp);
            }
        }
        if (write_image) {
            wbuf_printf(&w, " %lld %lld %lld",
                        (long long)D[i * (size_t)nc + c->idx_ix],
                        (long long)D[i * (size_t)nc + c->idx_iy],
                        (long long)D[i * (size_t)nc + c->idx_iz]);
        }
        wbuf_puts(&w, "\n");
    }

    if (c->has_velocities) {
        wbuf_puts(&w, "\nVelocities\n\n");
        int idv[3] = {c->idx_vx, c->idx_vy, c->idx_vz};
        for (size_t i = 0; i < N; i++) {
            long long id = c->idx_id >= 0
                ? (long long)D[i * (size_t)nc + c->idx_id]
                : (long long)(i + 1);
            wbuf_printf(&w, "%lld", id);
            for (int k = 0; k < 3; k++) {
                fmt_double(tmp, sizeof(tmp), D[i * (size_t)nc + idv[k]]);
                wbuf_puts(&w, " ");
                wbuf_puts(&w, tmp);
            }
            wbuf_puts(&w, "\n");
        }
    }

    return wbuf_close(&w);
}

/* -------------------------------------------------------------------------- */
/* LAMMPS dump writer                                                         */
/* -------------------------------------------------------------------------- */

static int write_lammps_dump(const Config *c, const char *path) {
    WBuf w;
    if (wbuf_open(&w, path) != 0) return -1;

    wbuf_puts(&w, "ITEM: TIMESTEP\n0\n");
    wbuf_printf(&w, "ITEM: NUMBER OF ATOMS\n%zu\n", c->natoms);
    if (c->triclinic) {
        wbuf_puts(&w, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n");
        wbuf_printf(&w, "%.16g %.16g %.16g\n", c->xlo, c->xhi, c->xy);
        wbuf_printf(&w, "%.16g %.16g %.16g\n", c->ylo, c->yhi, c->xz);
        wbuf_printf(&w, "%.16g %.16g %.16g\n", c->zlo, c->zhi, c->yz);
    } else {
        wbuf_puts(&w, "ITEM: BOX BOUNDS pp pp pp\n");
        wbuf_printf(&w, "%.16g %.16g\n", c->xlo, c->xhi);
        wbuf_printf(&w, "%.16g %.16g\n", c->ylo, c->yhi);
        wbuf_printf(&w, "%.16g %.16g\n", c->zlo, c->zhi);
    }

    wbuf_puts(&w, "ITEM: ATOMS");
    for (int k = 0; k < c->ncols; k++) wbuf_printf(&w, " %s", c->col_names[k]);
    wbuf_puts(&w, "\n");

    size_t N = c->natoms;
    int nc = c->ncols;
    const double *D = c->data;
    char tmp[64];
    for (size_t i = 0; i < N; i++) {
        for (int k = 0; k < nc; k++) {
            if (k) wbuf_puts(&w, " ");
            double v = D[i * (size_t)nc + k];
            if (c->col_is_int[k]) {
                wbuf_printf(&w, "%lld", (long long)v);
            } else {
                fmt_double(tmp, sizeof(tmp), v);
                wbuf_puts(&w, tmp);
            }
        }
        wbuf_puts(&w, "\n");
    }
    return wbuf_close(&w);
}

/* -------------------------------------------------------------------------- */
/* IMD writer                                                                 */
/* -------------------------------------------------------------------------- */

static int write_imd(const Config *c, const char *path) {
    WBuf w;
    if (wbuf_open(&w, path) != 0) return -1;

    int has_id   = c->idx_id   >= 0;
    int has_type = c->idx_type >= 0;
    int has_mass = config_find_col(c, "mass") >= 0;
    int has_pos  = c->idx_x >= 0 && c->idx_y >= 0 && c->idx_z >= 0;
    int has_vel  = c->has_velocities;

    if (!has_pos) {
        fprintf(stderr, "atomio: cannot write IMD without x/y/z columns\n");
        return -1;
    }

    /* Gather extra data columns (anything not in the fixed set). */
    int data_cols[MAX_COLS];
    int n_data = 0;
    for (int k = 0; k < c->ncols; k++) {
        const char *nm = c->col_names[k];
        if (!strcmp(nm, "id") || !strcmp(nm, "type") || !strcmp(nm, "mass") ||
            !strcmp(nm, "x")  || !strcmp(nm, "y")    || !strcmp(nm, "z")  ||
            !strcmp(nm, "vx") || !strcmp(nm, "vy")   || !strcmp(nm, "vz") ||
            !strcmp(nm, "ix") || !strcmp(nm, "iy")   || !strcmp(nm, "iz"))
            continue;
        data_cols[n_data++] = k;
    }

    wbuf_printf(&w, "#F A %d %d %d 3 %d %d\n",
                has_id ? 1 : 0, has_type ? 1 : 0, has_mass ? 1 : 0,
                has_vel ? 3 : 0, n_data);

    wbuf_puts(&w, "#C");
    if (has_id)   wbuf_puts(&w, " number");
    if (has_type) wbuf_puts(&w, " type");
    if (has_mass) wbuf_puts(&w, " mass");
    wbuf_puts(&w, " x y z");
    if (has_vel)  wbuf_puts(&w, " vx vy vz");
    for (int k = 0; k < n_data; k++)
        wbuf_printf(&w, " %s", c->col_names[data_cols[k]]);
    wbuf_puts(&w, "\n");

    double Lx = c->xhi - c->xlo;
    double Ly = c->yhi - c->ylo;
    double Lz = c->zhi - c->zlo;
    wbuf_printf(&w, "#X %.16g 0 0\n", Lx);
    wbuf_printf(&w, "#Y 0 %.16g 0\n", Ly);
    wbuf_printf(&w, "#Z 0 0 %.16g\n", Lz);
    wbuf_puts(&w, "#E\n");

    int mass_col = config_find_col(c, "mass");
    const double *D = c->data;
    size_t N = c->natoms;
    int nc = c->ncols;
    char tmp[64];

    for (size_t i = 0; i < N; i++) {
        int first = 1;
#define EMIT_INT(v) do {                                      \
            if (!first) wbuf_puts(&w, " ");                   \
            first = 0;                                        \
            wbuf_printf(&w, "%lld", (long long)(v));          \
        } while (0)
#define EMIT_DBL(v) do {                                      \
            if (!first) wbuf_puts(&w, " ");                   \
            first = 0;                                        \
            fmt_double(tmp, sizeof(tmp), (v));                \
            wbuf_puts(&w, tmp);                               \
        } while (0)

        if (has_id)   EMIT_INT(D[i * (size_t)nc + c->idx_id]);
        if (has_type) EMIT_INT(D[i * (size_t)nc + c->idx_type]);
        if (has_mass) EMIT_DBL(D[i * (size_t)nc + mass_col]);
        EMIT_DBL(D[i * (size_t)nc + c->idx_x] - c->xlo);
        EMIT_DBL(D[i * (size_t)nc + c->idx_y] - c->ylo);
        EMIT_DBL(D[i * (size_t)nc + c->idx_z] - c->zlo);
        if (has_vel) {
            EMIT_DBL(D[i * (size_t)nc + c->idx_vx]);
            EMIT_DBL(D[i * (size_t)nc + c->idx_vy]);
            EMIT_DBL(D[i * (size_t)nc + c->idx_vz]);
        }
        for (int k = 0; k < n_data; k++)
            EMIT_DBL(D[i * (size_t)nc + data_cols[k]]);
        wbuf_puts(&w, "\n");
#undef EMIT_INT
#undef EMIT_DBL
    }
    return wbuf_close(&w);
}

/* -------------------------------------------------------------------------- */
/* Top-level write                                                            */
/* -------------------------------------------------------------------------- */

int config_write(const Config *c, const char *filename, Format fmt) {
    if (fmt == FMT_UNKNOWN) fmt = FMT_LAMMPS_DATA;
    switch (fmt) {
    case FMT_LAMMPS_DATA: return write_lammps_data(c, filename);
    case FMT_LAMMPS_DUMP: return write_lammps_dump(c, filename);
    case FMT_IMD:         return write_imd(c, filename);
    default:
        fprintf(stderr, "atomio: unsupported output format\n");
        return -1;
    }
}
