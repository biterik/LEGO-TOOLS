/* lego-change-box — modify the simulation box bounds (orthogonal only).
 *
 * Two modes, chosen by argument count:
 *
 *   lego-change-box dx dy dz INPUT OUTPUT
 *       Add dx/dy/dz to the current xhi/yhi/zhi (negative shrinks).
 *       xlo/ylo/zlo are left unchanged.
 *
 *   lego-change-box xlo xhi ylo yhi zlo zhi INPUT OUTPUT
 *       Replace the bounds with explicit values.
 *
 * Streams the file line-by-line and rewrites only the box-bound lines.
 * The output format always matches the input format, and atoms are
 * copied through untouched — no full parse of the Atoms section, so
 * runtime is dominated by gzip de/compression, not atom processing.
 *
 * Supported inputs:
 *   - LAMMPS data   (lines "<lo> <hi> xlo xhi" / "ylo yhi" / "zlo zhi")
 *   - LAMMPS dump   (the 3 lines after "ITEM: BOX BOUNDS ...")
 *   - IMD ASCII     (#X / #Y / #Z header lines)
 *
 * Triclinic inputs are refused (cubic/orthogonal boxes only for now).
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <zlib.h>

static void usage(void) {
    fprintf(stderr,
        "Usage:\n"
        "  lego-change-box dx dy dz INPUT OUTPUT\n"
        "      Add dx dy dz to xhi yhi zhi (lo bounds unchanged).\n"
        "\n"
        "  lego-change-box xlo xhi ylo yhi zlo zhi INPUT OUTPUT\n"
        "      Replace the bounds with explicit values.\n"
        "\n"
        "Atoms are streamed through unchanged. Output format matches\n"
        "input format. Both files may have a .gz suffix. Orthogonal\n"
        "(cubic) boxes only.\n");
}

typedef struct {
    int    delta_mode;   /* 1 = add dx/dy/dz to hi; 0 = absolute */
    double dx, dy, dz;   /* deltas (delta mode) */
    double xlo, xhi;     /* absolute (!delta mode) */
    double ylo, yhi;
    double zlo, zhi;
} Params;

/* -------------------------------------------------------------------------- */
/* Streaming output                                                           */
/* -------------------------------------------------------------------------- */

typedef struct {
    gzFile gf;
    FILE  *fp;
    int    use_gz;
} Out;

static int out_open(Out *o, const char *path) {
    size_t n = strlen(path);
    o->use_gz = n > 3 && strcmp(path + n - 3, ".gz") == 0;
    if (o->use_gz) {
        o->gf = gzopen(path, "wb");
        o->fp = NULL;
#if ZLIB_VERNUM >= 0x1235
        if (o->gf) gzbuffer(o->gf, 1 << 20);
#endif
        return o->gf ? 0 : -1;
    }
    o->fp = fopen(path, "wb");
    o->gf = NULL;
    if (o->fp) setvbuf(o->fp, NULL, _IOFBF, 1 << 20);
    return o->fp ? 0 : -1;
}

static int out_write(Out *o, const char *data, size_t n) {
    if (n == 0) return 0;
    if (o->use_gz) {
        int r = gzwrite(o->gf, data, (unsigned)n);
        return r == (int)n ? 0 : -1;
    }
    return fwrite(data, 1, n, o->fp) == n ? 0 : -1;
}

static int out_close(Out *o) {
    int rc = 0;
    if (o->use_gz) {
        if (gzclose(o->gf) != Z_OK) rc = -1;
    } else if (o->fp) {
        if (fclose(o->fp) != 0) rc = -1;
    }
    return rc;
}

/* -------------------------------------------------------------------------- */
/* Line reader                                                                */
/* -------------------------------------------------------------------------- */

/* Growable line reader. Returns 0 on success, 1 on EOF, -1 on error.
 * out_line is reallocated as needed; *out_len is the length including any
 * trailing newline characters. */
static int read_line(gzFile gf, char **out_line, size_t *out_cap,
                     size_t *out_len) {
    size_t len = 0;
    if (*out_cap < 256) {
        *out_cap = 256;
        *out_line = (char *)realloc(*out_line, *out_cap);
        if (!*out_line) return -1;
    }
    for (;;) {
        if (gzgets(gf, *out_line + len, (int)(*out_cap - len)) == NULL) {
            if (len == 0) return 1;
            *out_len = len;
            return 0;
        }
        len += strlen(*out_line + len);
        if (len > 0 && (*out_line)[len - 1] == '\n') {
            *out_len = len;
            return 0;
        }
        /* No newline yet — grow and keep reading. */
        *out_cap *= 2;
        char *nb = (char *)realloc(*out_line, *out_cap);
        if (!nb) return -1;
        *out_line = nb;
    }
}

/* -------------------------------------------------------------------------- */
/* Format detection (from the first non-blank line)                           */
/* -------------------------------------------------------------------------- */

typedef enum { FMT_UNKNOWN, FMT_DATA, FMT_DUMP, FMT_IMD } Fmt;

static Fmt detect_format(const char *path) {
    gzFile gf = gzopen(path, "rb");
    if (!gf) return FMT_UNKNOWN;
    char head[1024];
    int n = gzread(gf, head, (int)sizeof(head) - 1);
    gzclose(gf);
    if (n <= 0) return FMT_UNKNOWN;
    head[n] = 0;
    const char *p = head;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    if (strncmp(p, "ITEM:", 5) == 0) return FMT_DUMP;
    if (*p == '#' && (p[1] == 'F' || p[1] == 'C' || p[1] == 'X' || p[1] == 'E'))
        return FMT_IMD;
    return FMT_DATA;
}

/* -------------------------------------------------------------------------- */
/* Rewriters                                                                  */
/* -------------------------------------------------------------------------- */

/* Strip trailing \n / \r\n and report what was stripped. */
static size_t line_content_len(const char *line, size_t len, const char **eol) {
    *eol = "";
    if (len >= 2 && line[len - 2] == '\r' && line[len - 1] == '\n') {
        *eol = "\r\n";
        return len - 2;
    }
    if (len >= 1 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
        *eol = line[len - 1] == '\n' ? "\n" : "\r";
        return len - 1;
    }
    return len;
}

/* Returns 1 if the trimmed line looks like "<num> <num> [<num>] SUFFIX". */
static int line_matches_suffix(const char *line, size_t clen,
                               const char *suffix) {
    size_t slen = strlen(suffix);
    if (clen < slen) return 0;
    /* Walk back over whitespace. */
    size_t e = clen;
    while (e > 0 && (line[e - 1] == ' ' || line[e - 1] == '\t')) e--;
    if (e < slen) return 0;
    if (memcmp(line + e - slen, suffix, slen) != 0) return 0;
    /* The char before suffix (if any) must be whitespace. */
    if (e - slen > 0) {
        char c = line[e - slen - 1];
        if (c != ' ' && c != '\t') return 0;
    }
    return 1;
}

/* Compute new (lo, hi) for one axis from the existing values and params. */
static void apply_axis(const Params *p, int axis, double cur_lo, double cur_hi,
                       double *new_lo, double *new_hi) {
    if (p->delta_mode) {
        *new_lo = cur_lo;
        if (axis == 0)      *new_hi = cur_hi + p->dx;
        else if (axis == 1) *new_hi = cur_hi + p->dy;
        else                *new_hi = cur_hi + p->dz;
    } else {
        if (axis == 0)      { *new_lo = p->xlo; *new_hi = p->xhi; }
        else if (axis == 1) { *new_lo = p->ylo; *new_hi = p->yhi; }
        else                { *new_lo = p->zlo; *new_hi = p->zhi; }
    }
}

/* Process a LAMMPS data file. */
static int rewrite_lammps_data(gzFile in, Out *out, const Params *p) {
    char   *line = NULL;
    size_t  cap  = 0, len = 0;
    int     rc;
    int     in_header = 1; /* box lines only appear in header */

    while ((rc = read_line(in, &line, &cap, &len)) == 0) {
        const char *eol;
        size_t clen = line_content_len(line, len, &eol);

        /* Detect section header (first uppercase word) to leave "header" mode. */
        if (in_header) {
            size_t i = 0;
            while (i < clen && (line[i] == ' ' || line[i] == '\t')) i++;
            if (i < clen && isupper((unsigned char)line[i])) {
                static const char *keywords[] = {
                    "Atoms", "Velocities", "Masses", "Bonds", "Angles",
                    "Dihedrals", "Impropers", "Pair", "PairIJ", "Bond",
                    "Angle", "Dihedral", "Improper", "Ellipsoids", "Lines",
                    "Triangles", "Bodies", "CMAP", "AtomFile", NULL
                };
                for (int k = 0; keywords[k]; k++) {
                    size_t kl = strlen(keywords[k]);
                    if (i + kl <= clen &&
                        memcmp(line + i, keywords[k], kl) == 0) {
                        char next = (i + kl < clen) ? line[i + kl] : '\n';
                        if (next == ' ' || next == '\t' || next == '\r' ||
                            next == '\n' || next == '#' || next == 0) {
                            in_header = 0;
                            break;
                        }
                    }
                }
            }
        }

        if (in_header) {
            if (line_matches_suffix(line, clen, "xy xz yz")) {
                fprintf(stderr,
                        "lego-change-box: input has xy xz yz tilt factors "
                        "(triclinic); refusing. Use orthogonal input.\n");
                free(line);
                return 1;
            }
            int axis = -1;
            const char *suffix = NULL;
            if      (line_matches_suffix(line, clen, "xlo xhi")) { axis = 0; suffix = "xlo xhi"; }
            else if (line_matches_suffix(line, clen, "ylo yhi")) { axis = 1; suffix = "ylo yhi"; }
            else if (line_matches_suffix(line, clen, "zlo zhi")) { axis = 2; suffix = "zlo zhi"; }
            if (axis >= 0) {
                double cur_lo = 0, cur_hi = 0;
                if (sscanf(line, "%lf %lf", &cur_lo, &cur_hi) != 2) {
                    fprintf(stderr,
                            "lego-change-box: cannot parse %s line\n", suffix);
                    free(line);
                    return 1;
                }
                double nlo, nhi;
                apply_axis(p, axis, cur_lo, cur_hi, &nlo, &nhi);
                char buf[160];
                int n = snprintf(buf, sizeof(buf), "%.16g %.16g %s%s",
                                 nlo, nhi, suffix, eol);
                if (out_write(out, buf, (size_t)n) != 0) { free(line); return 1; }
                continue;
            }
        }
        if (out_write(out, line, len) != 0) { free(line); return 1; }
    }
    free(line);
    return rc < 0 ? 1 : 0;
}

/* Process a LAMMPS dump file. */
static int rewrite_lammps_dump(gzFile in, Out *out, const Params *p) {
    char   *line = NULL;
    size_t  cap  = 0, len = 0;
    int     rc;
    int     box_lines_left = 0; /* after "ITEM: BOX BOUNDS" */
    int     axis = 0;

    while ((rc = read_line(in, &line, &cap, &len)) == 0) {
        const char *eol;
        size_t clen = line_content_len(line, len, &eol);

        if (clen >= 16 && memcmp(line, "ITEM: BOX BOUNDS", 16) == 0) {
            if (strstr(line, "xy") != NULL) {
                fprintf(stderr,
                        "lego-change-box: input dump has xy/xz/yz tilt "
                        "(triclinic); refusing.\n");
                free(line);
                return 1;
            }
            box_lines_left = 3;
            axis = 0;
            if (out_write(out, line, len) != 0) { free(line); return 1; }
            continue;
        }

        if (box_lines_left > 0) {
            double cur_lo = 0, cur_hi = 0;
            if (sscanf(line, "%lf %lf", &cur_lo, &cur_hi) != 2) {
                fprintf(stderr,
                        "lego-change-box: cannot parse dump BOX BOUNDS line\n");
                free(line);
                return 1;
            }
            double nlo, nhi;
            apply_axis(p, axis, cur_lo, cur_hi, &nlo, &nhi);
            char buf[128];
            int n = snprintf(buf, sizeof(buf), "%.16g %.16g%s", nlo, nhi, eol);
            if (out_write(out, buf, (size_t)n) != 0) { free(line); return 1; }
            box_lines_left--;
            axis++;
            continue;
        }

        if (out_write(out, line, len) != 0) { free(line); return 1; }
    }
    free(line);
    return rc < 0 ? 1 : 0;
}

/* Process an IMD ASCII checkpoint. IMD body coordinates are referenced
 * to the origin (0,0,0); the #X/#Y/#Z lines carry the box vectors, and
 * the box length on each axis equals the diagonal component. */
static int rewrite_imd(gzFile in, Out *out, const Params *p) {
    char   *line = NULL;
    size_t  cap  = 0, len = 0;
    int     rc;
    int     in_header = 1;

    while ((rc = read_line(in, &line, &cap, &len)) == 0) {
        const char *eol;
        size_t clen = line_content_len(line, len, &eol);

        if (in_header && clen >= 2 && line[0] == '#') {
            if (line[1] == 'E') {
                in_header = 0;
            } else if (line[1] == 'X' || line[1] == 'Y' || line[1] == 'Z') {
                int axis = line[1] == 'X' ? 0 : (line[1] == 'Y' ? 1 : 2);
                double a = 0, b = 0, d = 0;
                sscanf(line + 2, " %lf %lf %lf", &a, &b, &d);
                double off0 = (axis == 0) ? b : (axis == 1 ? a : a);
                double off1 = (axis == 2) ? b : (axis == 0 ? d : d);
                if (off0 != 0.0 || off1 != 0.0) {
                    fprintf(stderr,
                            "lego-change-box: IMD #%c has tilt components; "
                            "refusing (cubic boxes only).\n", line[1]);
                    free(line);
                    return 1;
                }
                /* current axis length */
                double cur_len =
                    axis == 0 ? a : (axis == 1 ? b : d);
                double new_len;
                if (p->delta_mode) {
                    double delta = axis == 0 ? p->dx :
                                   axis == 1 ? p->dy : p->dz;
                    new_len = cur_len + delta;
                } else {
                    new_len = axis == 0 ? (p->xhi - p->xlo) :
                              axis == 1 ? (p->yhi - p->ylo) :
                                          (p->zhi - p->zlo);
                }
                char buf[160];
                int n = snprintf(buf, sizeof(buf),
                                 "#%c %s%s",
                                 line[1],
                                 axis == 0 ? "" : (axis == 1 ? "0 " : "0 0 "),
                                 "");
                /* Rebuild properly. */
                if (axis == 0)
                    n = snprintf(buf, sizeof(buf), "#X %.16g 0 0%s",
                                 new_len, eol);
                else if (axis == 1)
                    n = snprintf(buf, sizeof(buf), "#Y 0 %.16g 0%s",
                                 new_len, eol);
                else
                    n = snprintf(buf, sizeof(buf), "#Z 0 0 %.16g%s",
                                 new_len, eol);
                if (out_write(out, buf, (size_t)n) != 0) { free(line); return 1; }
                continue;
            }
        }

        if (out_write(out, line, len) != 0) { free(line); return 1; }
    }
    free(line);
    return rc < 0 ? 1 : 0;
}

/* -------------------------------------------------------------------------- */
/* main                                                                       */
/* -------------------------------------------------------------------------- */

int main(int argc, char **argv) {
    const char *pos[16];
    int npos = 0;

    for (int ai = 1; ai < argc; ai++) {
        const char *a = argv[ai];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(); return 0; }
        if (!strcmp(a, "--")) continue;
        if (a[0] == '-' && !(a[1] >= '0' && a[1] <= '9') && a[1] != '.') {
            fprintf(stderr, "lego-change-box: unknown option %s\n", a);
            return 1;
        }
        if (npos >= 16) { usage(); return 1; }
        pos[npos++] = a;
    }

    Params p;
    memset(&p, 0, sizeof(p));
    const char *input, *output;

    if (npos == 5) {
        p.delta_mode = 1;
        p.dx = strtod(pos[0], NULL);
        p.dy = strtod(pos[1], NULL);
        p.dz = strtod(pos[2], NULL);
        input  = pos[3];
        output = pos[4];
    } else if (npos == 8) {
        p.delta_mode = 0;
        p.xlo = strtod(pos[0], NULL);
        p.xhi = strtod(pos[1], NULL);
        p.ylo = strtod(pos[2], NULL);
        p.yhi = strtod(pos[3], NULL);
        p.zlo = strtod(pos[4], NULL);
        p.zhi = strtod(pos[5], NULL);
        input  = pos[6];
        output = pos[7];
        if (p.xhi <= p.xlo || p.yhi <= p.ylo || p.zhi <= p.zlo) {
            fprintf(stderr,
                    "lego-change-box: each lo must be strictly less than hi\n");
            return 1;
        }
    } else {
        usage();
        return 1;
    }

    Fmt fmt = detect_format(input);
    if (fmt == FMT_UNKNOWN) {
        fprintf(stderr, "lego-change-box: cannot detect format of %s\n",
                input);
        return 1;
    }

    gzFile in = gzopen(input, "rb");
    if (!in) {
        fprintf(stderr, "lego-change-box: cannot open %s: %s\n",
                input, strerror(errno));
        return 1;
    }
#if ZLIB_VERNUM >= 0x1235
    gzbuffer(in, 1 << 20);
#endif

    Out out;
    if (out_open(&out, output) != 0) {
        fprintf(stderr, "lego-change-box: cannot open output %s\n", output);
        gzclose(in);
        return 1;
    }

    int rc = 1;
    switch (fmt) {
    case FMT_DATA: rc = rewrite_lammps_data(in, &out, &p); break;
    case FMT_DUMP: rc = rewrite_lammps_dump(in, &out, &p); break;
    case FMT_IMD:  rc = rewrite_imd       (in, &out, &p); break;
    default: break;
    }

    gzclose(in);
    if (out_close(&out) != 0) rc = 1;
    return rc;
}
