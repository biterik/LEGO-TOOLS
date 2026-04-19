/* converter.c — Atomic Format Converter implementation.
 *
 * Readers/writers for Extended XYZ, LAMMPS data (atomic style), and
 * CEL (Dr Probe supercell).  gz-transparent I/O via zlib.
 * OpenMP-parallel parsing of per-atom data rows.
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#define _POSIX_C_SOURCE 200809L
#include "converter.h"
#include "elements.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <zlib.h>

/* -------------------------------------------------------------------------- */
/* Large-file I/O chunk size                                                  */
/* -------------------------------------------------------------------------- */
/*
 * zlib's legacy gzread/gzwrite take an `unsigned int` length and return an
 * `int`, so any single call transferring more than ~2 GiB silently
 * truncates (older zlib) or fails with -1 (zlib >= 1.2.9).  Every read /
 * write call below is clamped to GZ_CHUNK_MAX bytes so that multi-GiB
 * files work regardless of whether the underlying file is gzipped --
 * gzopen transparently passes uncompressed data through, with the same
 * unsigned-int size ceiling on each gzread call.
 *
 * Default GZ_CHUNK_MAX is 1 GiB, well below INT_MAX.  Override at compile
 * time with e.g. -DGZ_CHUNK_MAX=4096 to exercise the chunking loop in
 * tests.  Mirrors the identical knob in lego-tools/atomio.c.
 */
#ifndef GZ_CHUNK_MAX
#define GZ_CHUNK_MAX ((size_t)1 << 30)   /* 1 GiB */
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================== */
/* Small helpers                                                       */
/* ================================================================== */

static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "afc: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(1);
}

static int str_has_suffix(const char *s, const char *suf)
{
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

static inline const char *skip_ws(const char *p, const char *end)
{
    while (p < end && (*p == ' ' || *p == '\t')) p++;
    return p;
}

static inline const char *skip_to_eol(const char *p, const char *end)
{
    while (p < end && *p != '\n') p++;
    return p;
}

static inline const char *next_line(const char *p, const char *end)
{
    p = skip_to_eol(p, end);
    return (p < end) ? p + 1 : end;
}

/* ================================================================== */
/* Config lifecycle                                                    */
/* ================================================================== */

Config config_create(void)
{
    Config cfg;
    memset(&cfg, 0, sizeof(cfg));
    return cfg;
}

void config_free(Config *cfg)
{
    free(cfg->atoms.id);
    free(cfg->atoms.type);
    free(cfg->atoms.x);
    free(cfg->atoms.y);
    free(cfg->atoms.z);
    memset(cfg, 0, sizeof(*cfg));
}

void atoms_reserve(Atoms *a, int capacity)
{
    if (capacity <= a->capacity) return;
    a->id   = (int *)   realloc(a->id,   (size_t)capacity * sizeof(int));
    a->type = (int *)   realloc(a->type,  (size_t)capacity * sizeof(int));
    a->x    = (double *)realloc(a->x,     (size_t)capacity * sizeof(double));
    a->y    = (double *)realloc(a->y,     (size_t)capacity * sizeof(double));
    a->z    = (double *)realloc(a->z,     (size_t)capacity * sizeof(double));
    if (!a->id || !a->type || !a->x || !a->y || !a->z)
        die("out of memory allocating %d atoms", capacity);
    a->capacity = capacity;
}

/* ================================================================== */
/* Format helpers                                                      */
/* ================================================================== */

Format format_from_extension(const char *path)
{
    if (str_has_suffix(path, ".xyz.gz") || str_has_suffix(path, ".xyz"))
        return FMT_XYZ;
    if (str_has_suffix(path, ".lmp")  || str_has_suffix(path, ".lmp.gz") ||
        str_has_suffix(path, ".data") || str_has_suffix(path, ".data.gz") ||
        str_has_suffix(path, ".lammps") || str_has_suffix(path, ".lammps.gz"))
        return FMT_LAMMPS;
    if (str_has_suffix(path, ".cel") || str_has_suffix(path, ".cel.gz"))
        return FMT_CEL;
    return FMT_UNKNOWN;
}

Format format_from_string(const char *s)
{
    if (!s) return FMT_UNKNOWN;
    if (!strcmp(s, "xyz"))    return FMT_XYZ;
    if (!strcmp(s, "lammps") || !strcmp(s, "lmp") || !strcmp(s, "data"))
        return FMT_LAMMPS;
    if (!strcmp(s, "cel"))    return FMT_CEL;
    return FMT_UNKNOWN;
}

const char *format_name(Format f)
{
    switch (f) {
    case FMT_XYZ:    return "xyz";
    case FMT_LAMMPS: return "lammps";
    case FMT_CEL:    return "cel";
    default:         return "unknown";
    }
}

/* ================================================================== */
/* gz-transparent file I/O                                             */
/* ================================================================== */

typedef struct {
    char   *data;
    size_t  len;
} FileBuf;

/* Is afc debug output requested via LEGO_DEBUG=1 ?  (Shared env var with
 * lego-tools so the whole suite has one debug knob.) */
static int afc_debug(void) {
    static int cached = -1;
    if (cached < 0) {
        const char *e = getenv("LEGO_DEBUG");
        cached = (e && *e && *e != '0') ? 1 : 0;
    }
    return cached;
}

/* Pretty-print a byte size into a caller-supplied buffer. */
static void afc_human_size(char *out, size_t outsz, size_t n) {
    static const char *units[] = {"B", "KiB", "MiB", "GiB", "TiB"};
    double v = (double)n;
    int u = 0;
    while (v >= 1024.0 && u < 4) { v /= 1024.0; u++; }
    snprintf(out, outsz, "%.2f %s (%zu bytes)", v, units[u], n);
}

static int file_read_all(const char *path, FileBuf *fb)
{
    int dbg = afc_debug();

    off_t file_size = -1;
    {
        struct stat st;
        if (stat(path, &st) == 0) file_size = st.st_size;
    }

    gzFile gf = gzopen(path, "rb");
    if (!gf) {
        fprintf(stderr, "afc: cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }
#if ZLIB_VERNUM >= 0x1235
    gzbuffer(gf, 1 << 20);
#endif
    if (dbg) {
        char sz[64] = "unknown";
        if (file_size >= 0) afc_human_size(sz, sizeof(sz), (size_t)file_size);
        fprintf(stderr, "afc: opening %s (file size: %s, chunk: %zu)\n",
                path, sz, (size_t)GZ_CHUNK_MAX);
    }
    size_t cap = 1 << 20, len = 0;
    char *buf = (char *)malloc(cap);
    if (!buf) {
        fprintf(stderr, "afc: out of memory allocating %zu-byte read buffer\n", cap);
        gzclose(gf);
        return -1;
    }
    size_t next_report = (size_t)1 << 30;   /* debug: report every 1 GiB read */
    for (;;) {
        if (len + (1 << 16) > cap) {
            size_t new_cap = cap * 2;
            if (new_cap <= cap) {
                fprintf(stderr,
                        "afc: read buffer size would overflow size_t at %zu bytes\n",
                        cap);
                free(buf); gzclose(gf);
                return -1;
            }
            char *nb = (char *)realloc(buf, new_cap);
            if (!nb) {
                fprintf(stderr,
                        "afc: out of memory growing read buffer to %zu bytes "
                        "(already read %zu)\n", new_cap, len);
                free(buf); gzclose(gf);
                return -1;
            }
            buf = nb;
            cap = new_cap;
        }
        size_t want = cap - len;
        if (want > GZ_CHUNK_MAX) want = GZ_CHUNK_MAX;
        int n = gzread(gf, buf + len, (unsigned)want);
        if (n < 0) {
            int gzerrno = 0;
            const char *gzmsg = gzerror(gf, &gzerrno);
            char sread[64], swant[64];
            afc_human_size(sread, sizeof(sread), len);
            afc_human_size(swant, sizeof(swant), want);
            fprintf(stderr,
                    "afc: read error on %s\n"
                    "     after reading %s, requested chunk %s\n"
                    "     zlib error: %s (errno=%d)\n"
                    "     (set LEGO_DEBUG=1 for verbose progress; "
                    "rebuild with -DGZ_CHUNK_MAX=<bytes> to shrink chunks)\n",
                    path, sread, swant, gzmsg ? gzmsg : "(none)", gzerrno);
            free(buf); gzclose(gf);
            return -1;
        }
        if (n == 0) break;
        len += (size_t)n;
        if (dbg && len >= next_report) {
            char sread[64];
            afc_human_size(sread, sizeof(sread), len);
            fprintf(stderr, "afc: read progress on %s: %s\n", path, sread);
            next_report += (size_t)1 << 30;
        }
    }
    gzclose(gf);
    char *nb = (char *)realloc(buf, len + 1);
    if (nb) buf = nb;
    buf[len] = 0;
    fb->data = buf;
    fb->len  = len;
    if (dbg) {
        char sread[64];
        afc_human_size(sread, sizeof(sread), len);
        fprintf(stderr, "afc: finished reading %s: %s\n", path, sread);
    }
    return 0;
}

/* ---- Buffered writer (gz-aware) ---- */

typedef struct {
    char   *buf;
    size_t  len, cap;
    gzFile  gf;
    FILE   *fp;
    int     use_gz;
} WBuf;

static int wbuf_open(WBuf *w, const char *path)
{
    memset(w, 0, sizeof(*w));
    w->cap = 1 << 16;
    w->buf = (char *)malloc(w->cap);
    if (!w->buf) return -1;
    w->use_gz = str_has_suffix(path, ".gz");
    if (w->use_gz) {
        w->gf = gzopen(path, "wb");
        if (!w->gf) { free(w->buf); return -1; }
#if ZLIB_VERNUM >= 0x1235
        gzbuffer(w->gf, 1 << 20);
#endif
    } else {
        w->fp = fopen(path, "wb");
        if (!w->fp) { free(w->buf); return -1; }
        setvbuf(w->fp, NULL, _IOFBF, 1 << 20);
    }
    return 0;
}

/* Flush w->buf.  Chunks every gzwrite() call to at most GZ_CHUNK_MAX bytes,
 * mirroring lego-tools/atomio.c.  Normal flushes happen at ~1 MiB so
 * chunking never triggers in practice, but this makes the path safe if a
 * caller ever buffers a very large amount at once. */
static int wbuf_flush(WBuf *w)
{
    if (w->len == 0) return 0;
    int rc = 0;
    if (w->use_gz) {
        size_t off = 0;
        while (off < w->len) {
            size_t want = w->len - off;
            if (want > GZ_CHUNK_MAX) want = GZ_CHUNK_MAX;
            int n = gzwrite(w->gf, w->buf + off, (unsigned)want);
            if (n <= 0 || (size_t)n != want) {
                int gzerrno = 0;
                const char *gzmsg = gzerror(w->gf, &gzerrno);
                char swant[64], sdone[64];
                afc_human_size(swant, sizeof(swant), want);
                afc_human_size(sdone, sizeof(sdone), off);
                fprintf(stderr,
                        "afc: gzwrite failed: asked %s after %s, "
                        "gzwrite returned %d, zlib error: %s (errno=%d)\n",
                        swant, sdone, n, gzmsg ? gzmsg : "(none)", gzerrno);
                rc = -1;
                break;
            }
            off += (size_t)n;
        }
    } else {
        size_t got = fwrite(w->buf, 1, w->len, w->fp);
        if (got != w->len) {
            char sdone[64], swant[64];
            afc_human_size(sdone, sizeof(sdone), got);
            afc_human_size(swant, sizeof(swant), w->len);
            fprintf(stderr,
                    "afc: short fwrite: wrote %s of %s (errno=%d: %s)\n",
                    sdone, swant, errno, strerror(errno));
            rc = -1;
        }
    }
    w->len = 0;
    return rc;
}

static int wbuf_close(WBuf *w)
{
    int rc = wbuf_flush(w);
    if (w->use_gz && w->gf) gzclose(w->gf);
    if (!w->use_gz && w->fp) fclose(w->fp);
    free(w->buf);
    memset(w, 0, sizeof(*w));
    return rc;
}

static void wbuf_ensure(WBuf *w, size_t need)
{
    if (w->len + need <= w->cap) return;
    if (w->cap >= (1 << 20)) {
        wbuf_flush(w);
        if (need > w->cap) {
            w->cap = need * 2;
            char *nb = (char *)realloc(w->buf, w->cap);
            if (!nb) die("oom");
            w->buf = nb;
        }
    } else {
        while (w->len + need > w->cap) w->cap *= 2;
        char *nb = (char *)realloc(w->buf, w->cap);
        if (!nb) die("oom");
        w->buf = nb;
    }
}

static void wbuf_write(WBuf *w, const char *s, size_t n)
{
    wbuf_ensure(w, n);
    memcpy(w->buf + w->len, s, n);
    w->len += n;
}

static void wbuf_puts(WBuf *w, const char *s) { wbuf_write(w, s, strlen(s)); }

static void wbuf_printf(WBuf *w, const char *fmt, ...)
{
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
    wbuf_ensure(w, (size_t)n + 1);
    va_start(ap, fmt);
    vsnprintf(w->buf + w->len, (size_t)n + 1, fmt, ap);
    va_end(ap);
    w->len += (size_t)n;
}

/* ================================================================== */
/* Box utilities                                                       */
/* ================================================================== */

void box_to_lammps(const Box *box, double lo[3], double hi[3], double tilt[3])
{
    /* Assumes upper-triangular cell:
     *   a = (cell[0][0], 0, 0)
     *   b = (cell[1][0], cell[1][1], 0)
     *   c = (cell[2][0], cell[2][1], cell[2][2])
     */
    lo[0] = box->origin[0];
    lo[1] = box->origin[1];
    lo[2] = box->origin[2];
    hi[0] = box->origin[0] + box->cell[0][0];
    hi[1] = box->origin[1] + box->cell[1][1];
    hi[2] = box->origin[2] + box->cell[2][2];
    tilt[0] = box->cell[1][0]; /* xy */
    tilt[1] = box->cell[2][0]; /* xz */
    tilt[2] = box->cell[2][1]; /* yz */
}

void box_from_lammps(Box *box, const double lo[3], const double hi[3],
                     const double tilt[3])
{
    memset(box, 0, sizeof(*box));
    box->cell[0][0] = hi[0] - lo[0];
    box->cell[1][0] = tilt[0]; /* xy */
    box->cell[1][1] = hi[1] - lo[1];
    box->cell[2][0] = tilt[1]; /* xz */
    box->cell[2][1] = tilt[2]; /* yz */
    box->cell[2][2] = hi[2] - lo[2];
    box->origin[0] = lo[0];
    box->origin[1] = lo[1];
    box->origin[2] = lo[2];
}

void box_to_cellparams(const Box *box, double *a, double *b, double *c,
                       double *alpha, double *beta, double *gamma)
{
    const double *va = box->cell[0];
    const double *vb = box->cell[1];
    const double *vc = box->cell[2];

    *a = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);
    *b = sqrt(vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2]);
    *c = sqrt(vc[0]*vc[0] + vc[1]*vc[1] + vc[2]*vc[2]);

    double ab = *a * *b, ac = *a * *c, bc = *b * *c;

    double cos_alpha = (vb[0]*vc[0] + vb[1]*vc[1] + vb[2]*vc[2]) / bc;
    double cos_beta  = (va[0]*vc[0] + va[1]*vc[1] + va[2]*vc[2]) / ac;
    double cos_gamma = (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]) / ab;

    /* clamp to [-1,1] for safety */
    if (cos_alpha >  1.0) cos_alpha =  1.0;
    if (cos_alpha < -1.0) cos_alpha = -1.0;
    if (cos_beta  >  1.0) cos_beta  =  1.0;
    if (cos_beta  < -1.0) cos_beta  = -1.0;
    if (cos_gamma >  1.0) cos_gamma =  1.0;
    if (cos_gamma < -1.0) cos_gamma = -1.0;

    *alpha = acos(cos_alpha) * 180.0 / M_PI;
    *beta  = acos(cos_beta)  * 180.0 / M_PI;
    *gamma = acos(cos_gamma) * 180.0 / M_PI;
}

void box_from_cellparams(Box *box, double a, double b, double c,
                         double alpha, double beta, double gamma)
{
    memset(box, 0, sizeof(*box));
    double ar = alpha * M_PI / 180.0;
    double br = beta  * M_PI / 180.0;
    double gr = gamma * M_PI / 180.0;

    box->cell[0][0] = a;
    box->cell[1][0] = b * cos(gr);
    box->cell[1][1] = b * sin(gr);
    box->cell[2][0] = c * cos(br);
    box->cell[2][1] = c * (cos(ar) - cos(br) * cos(gr)) / sin(gr);
    double cx = box->cell[2][0], cy = box->cell[2][1];
    box->cell[2][2] = sqrt(c * c - cx * cx - cy * cy);
}

int box_invert(const Box *box, double inv[3][3])
{
    const double (*m)[3] = box->cell;

    double det = m[0][0] * (m[1][1]*m[2][2] - m[1][2]*m[2][1])
               - m[0][1] * (m[1][0]*m[2][2] - m[1][2]*m[2][0])
               + m[0][2] * (m[1][0]*m[2][1] - m[1][1]*m[2][0]);

    if (fabs(det) < 1e-30) {
        fprintf(stderr, "afc: singular cell matrix (det = %g)\n", det);
        return -1;
    }
    double idet = 1.0 / det;

    inv[0][0] =  (m[1][1]*m[2][2] - m[1][2]*m[2][1]) * idet;
    inv[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1]) * idet;
    inv[0][2] =  (m[0][1]*m[1][2] - m[0][2]*m[1][1]) * idet;
    inv[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0]) * idet;
    inv[1][1] =  (m[0][0]*m[2][2] - m[0][2]*m[2][0]) * idet;
    inv[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0]) * idet;
    inv[2][0] =  (m[1][0]*m[2][1] - m[1][1]*m[2][0]) * idet;
    inv[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0]) * idet;
    inv[2][2] =  (m[0][0]*m[1][1] - m[0][1]*m[1][0]) * idet;
    return 0;
}

/* ================================================================== */
/* Coordinate transforms                                               */
/* ================================================================== */

void cart_to_frac(const Box *box, int n,
                  const double *cx, const double *cy, const double *cz,
                  double *fx, double *fy, double *fz)
{
    double inv[3][3];
    if (box_invert(box, inv) != 0) die("cannot invert cell for cart_to_frac");

    double ox = box->origin[0], oy = box->origin[1], oz = box->origin[2];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n > 10000)
#endif
    for (int i = 0; i < n; i++) {
        double dx = cx[i] - ox;
        double dy = cy[i] - oy;
        double dz = cz[i] - oz;
        fx[i] = inv[0][0]*dx + inv[0][1]*dy + inv[0][2]*dz;
        fy[i] = inv[1][0]*dx + inv[1][1]*dy + inv[1][2]*dz;
        fz[i] = inv[2][0]*dx + inv[2][1]*dy + inv[2][2]*dz;
    }
}

void frac_to_cart(const Box *box, int n,
                  const double *fx, const double *fy, const double *fz,
                  double *cx, double *cy, double *cz)
{
    const double (*c)[3] = box->cell;
    double ox = box->origin[0], oy = box->origin[1], oz = box->origin[2];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n > 10000)
#endif
    for (int i = 0; i < n; i++) {
        cx[i] = ox + fx[i]*c[0][0] + fy[i]*c[1][0] + fz[i]*c[2][0];
        cy[i] = oy + fx[i]*c[0][1] + fy[i]*c[1][1] + fz[i]*c[2][1];
        cz[i] = oz + fx[i]*c[0][2] + fy[i]*c[1][2] + fz[i]*c[2][2];
    }
}

/* ================================================================== */
/* Line-index helper for parallel row parsing                          */
/* ================================================================== */

/* Build an array of byte offsets to the start of each non-blank,
 * non-comment line in [text, text+len).  Returns the count.           */
static int index_data_lines(const char *text, size_t len,
                            size_t **offsets_out, int max_lines)
{
    size_t *off = (size_t *)malloc((size_t)(max_lines + 1) * sizeof(size_t));
    if (!off) die("oom in index_data_lines");
    int n = 0;
    const char *p = text, *end = text + len;
    while (p < end && n < max_lines) {
        const char *eol = skip_to_eol(p, end);
        const char *q = skip_ws(p, eol);
        if (q < eol && *q != '#' && *q != '\r' && *q != '\n') {
            off[n++] = (size_t)(p - text);
        }
        p = (eol < end) ? eol + 1 : end;
    }
    *offsets_out = off;
    return n;
}

/* ================================================================== */
/* Extended XYZ reader                                                 */
/* ================================================================== */

/* Parse a Lattice="..." field from the comment line. */
static int parse_lattice(const char *start, const char *end, double cell[3][3])
{
    const char *p = start;
    while (p < end - 8) {
        if (strncmp(p, "Lattice=", 8) == 0) { p += 8; break; }
        p++;
    }
    if (p >= end - 8) return -1;

    /* skip opening quote */
    if (*p == '"') p++;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            while (p < end && (*p == ' ' || *p == '\t')) p++;
            char *nx;
            cell[i][j] = strtod(p, &nx);
            if (nx == p) return -1;
            p = nx;
        }
    }
    return 0;
}

/* Parse an Origin="..." field. */
static int parse_origin(const char *start, const char *end, double origin[3])
{
    const char *p = start;
    while (p < end - 7) {
        if (strncmp(p, "Origin=", 7) == 0) { p += 7; break; }
        p++;
    }
    if (p >= end - 7) {
        origin[0] = origin[1] = origin[2] = 0.0;
        return 0; /* Origin is optional */
    }
    if (*p == '"') p++;
    for (int i = 0; i < 3; i++) {
        while (p < end && (*p == ' ' || *p == '\t')) p++;
        char *nx;
        origin[i] = strtod(p, &nx);
        if (nx == p) return -1;
        p = nx;
    }
    return 0;
}

int read_xyz(Config *cfg, const char *path)
{
    FileBuf fb;
    if (file_read_all(path, &fb) != 0) return -1;

    const char *p = fb.data, *end = fb.data + fb.len;

    /* Line 1: atom count */
    char *nx;
    long natoms = strtol(p, &nx, 10);
    if (nx == p || natoms <= 0) {
        fprintf(stderr, "afc: invalid atom count in %s\n", path);
        free(fb.data); return -1;
    }
    p = next_line(nx, end);

    /* Line 2: comment / metadata */
    const char *line2 = p;
    const char *line2_end = skip_to_eol(p, end);

    if (parse_lattice(line2, line2_end, cfg->box.cell) != 0) {
        fprintf(stderr, "afc: cannot parse Lattice= in %s\n", path);
        free(fb.data); return -1;
    }
    parse_origin(line2, line2_end, cfg->box.origin);

    /* Save comment */
    size_t clen = (size_t)(line2_end - line2);
    if (clen >= sizeof(cfg->comment)) clen = sizeof(cfg->comment) - 1;
    memcpy(cfg->comment, line2, clen);
    cfg->comment[clen] = 0;

    p = next_line(line2_end, end);

    /* Allocate atoms */
    atoms_reserve(&cfg->atoms, (int)natoms);
    cfg->atoms.count = (int)natoms;

    /* Index data lines */
    size_t *line_off;
    int nlines = index_data_lines(p, (size_t)(end - p), &line_off, (int)natoms);
    if (nlines < (int)natoms) {
        fprintf(stderr, "afc: only %d data lines, expected %ld atoms in %s\n",
                nlines, natoms, path);
        free(line_off); free(fb.data); return -1;
    }

    /* Parallel parse */
    const char *base = p;
    int max_type = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_type) if(natoms > 10000)
#endif
    for (int i = 0; i < (int)natoms; i++) {
        const char *lp = base + line_off[i];
        const char *le = skip_to_eol(lp, end);
        lp = skip_ws(lp, le);

        char *ep;
        int id   = (int)strtol(lp, &ep, 10); lp = ep;
        int type = (int)strtol(lp, &ep, 10); lp = ep;
        double x = strtod(lp, &ep);           lp = ep;
        double y = strtod(lp, &ep);           lp = ep;
        double z = strtod(lp, &ep);

        cfg->atoms.id[i]   = id;
        cfg->atoms.type[i] = type;
        cfg->atoms.x[i]    = x;
        cfg->atoms.y[i]    = y;
        cfg->atoms.z[i]    = z;
        if (type > max_type) max_type = type;
    }

    cfg->typemap.ntypes = max_type;

    free(line_off);
    free(fb.data);
    return 0;
}

/* ================================================================== */
/* Extended XYZ writer                                                 */
/* ================================================================== */

int write_xyz(const Config *cfg, const char *path)
{
    WBuf w;
    if (wbuf_open(&w, path) != 0) {
        fprintf(stderr, "afc: cannot open %s for writing\n", path);
        return -1;
    }

    int n = cfg->atoms.count;
    const Box *box = &cfg->box;

    /* Line 1 */
    wbuf_printf(&w, "%d\n", n);

    /* Line 2 */
    wbuf_printf(&w,
        "Lattice=\"%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\" "
        "Origin=\"%.10g %.10g %.10g\" "
        "Properties=id:I:1:species:S:1:pos:R:3\n",
        box->cell[0][0], box->cell[0][1], box->cell[0][2],
        box->cell[1][0], box->cell[1][1], box->cell[1][2],
        box->cell[2][0], box->cell[2][1], box->cell[2][2],
        box->origin[0], box->origin[1], box->origin[2]);

    /* Determine if we have element symbols for output */
    int have_symbols = 0;
    for (int t = 1; t <= cfg->typemap.ntypes; t++) {
        if (cfg->typemap.symbols[t][0]) { have_symbols = 1; break; }
    }

    /* Atom lines — parallel format into per-thread buffers */
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
    if (nth > 64) nth = 64;
#endif

    if (n < 10000 || nth < 2) {
        /* Serial path */
        for (int i = 0; i < n; i++) {
            if (have_symbols && cfg->atoms.type[i] >= 1 &&
                cfg->atoms.type[i] <= cfg->typemap.ntypes) {
                wbuf_printf(&w, "%d %s %.10g %.10g %.10g\n",
                    cfg->atoms.id[i],
                    cfg->typemap.symbols[cfg->atoms.type[i]],
                    cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
            } else {
                wbuf_printf(&w, "%d %d %.10g %.10g %.10g\n",
                    cfg->atoms.id[i], cfg->atoms.type[i],
                    cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
            }
        }
    } else {
        /* Parallel: each thread fills a local buffer, then we write them in order */
        size_t est_per_atom = 80;
        char **tbufs = (char **)calloc((size_t)nth, sizeof(char *));
        size_t *tlens = (size_t *)calloc((size_t)nth, sizeof(size_t));

        for (int t = 0; t < nth; t++) {
            int lo = (n * t) / nth;
            int hi = (n * (t + 1)) / nth;
            tbufs[t] = (char *)malloc((size_t)(hi - lo) * est_per_atom);
            if (!tbufs[t]) die("oom");
        }

#pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            int lo = (n * tid) / nth;
            int hi = (n * (tid + 1)) / nth;
            char *bp = tbufs[tid];
            size_t blen = 0;
            size_t bcap = (size_t)(hi - lo) * est_per_atom;

            for (int i = lo; i < hi; i++) {
                if (blen + 256 > bcap) {
                    bcap *= 2;
                    size_t offset = (size_t)(bp - tbufs[tid]);
                    tbufs[tid] = (char *)realloc(tbufs[tid], bcap);
                    if (!tbufs[tid]) die("oom");
                    bp = tbufs[tid] + offset;
                }
                int wrote;
                if (have_symbols && cfg->atoms.type[i] >= 1 &&
                    cfg->atoms.type[i] <= cfg->typemap.ntypes) {
                    wrote = snprintf(bp + blen, bcap - blen,
                        "%d %s %.10g %.10g %.10g\n",
                        cfg->atoms.id[i],
                        cfg->typemap.symbols[cfg->atoms.type[i]],
                        cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
                } else {
                    wrote = snprintf(bp + blen, bcap - blen,
                        "%d %d %.10g %.10g %.10g\n",
                        cfg->atoms.id[i], cfg->atoms.type[i],
                        cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
                }
                blen += (size_t)wrote;
            }
            tlens[tid] = blen;
        }

        for (int t = 0; t < nth; t++) {
            wbuf_write(&w, tbufs[t], tlens[t]);
            free(tbufs[t]);
        }
        free(tbufs);
        free(tlens);
    }

    return wbuf_close(&w);
}

/* ================================================================== */
/* LAMMPS data reader (atomic style)                                   */
/* ================================================================== */

int read_lammps(Config *cfg, const char *path)
{
    FileBuf fb;
    if (file_read_all(path, &fb) != 0) return -1;

    const char *p = fb.data, *end = fb.data + fb.len;
    long natoms = 0;
    int ntypes = 0;
    double lo[3] = {0}, hi[3] = {0}, tilt[3] = {0};
    int has_tilt = 0;

    /* Line 1: title */
    const char *title_end = skip_to_eol(p, end);
    size_t tlen = (size_t)(title_end - p);
    if (tlen >= sizeof(cfg->comment)) tlen = sizeof(cfg->comment) - 1;
    memcpy(cfg->comment, p, tlen);
    cfg->comment[tlen] = 0;
    /* strip trailing \r */
    if (tlen > 0 && cfg->comment[tlen - 1] == '\r')
        cfg->comment[--tlen] = 0;
    p = next_line(title_end, end);

    /* Parse header lines until we hit "Atoms" */
    const char *atoms_section = NULL;
    while (p < end) {
        const char *lp = p;
        const char *le = skip_to_eol(p, end);
        const char *q = skip_ws(lp, le);

        /* Check for "Atoms" section header */
        if (le - q >= 5 && strncmp(q, "Atoms", 5) == 0) {
            atoms_section = next_line(le, end);
            /* skip blank line after "Atoms ..." header */
            const char *tmp = atoms_section;
            const char *tmp_end = skip_to_eol(tmp, end);
            const char *tmp_q = skip_ws(tmp, tmp_end);
            if (tmp_q >= tmp_end || *tmp_q == '\n' || *tmp_q == '\r')
                atoms_section = next_line(tmp_end, end);
            break;
        }

        /* Skip blank / comment lines */
        if (q >= le || *q == '#' || *q == '\r' || *q == '\n') {
            p = next_line(le, end);
            continue;
        }

        /* Try to parse header keywords */
        char tmp[256];
        size_t ll = (size_t)(le - q);
        if (ll >= sizeof(tmp)) ll = sizeof(tmp) - 1;
        memcpy(tmp, q, ll);
        tmp[ll] = 0;

        int matched = -1;
        long long ln = 0;
        int ni = 0;
        double v1, v2, v3;

        if (sscanf(tmp, "%lld atoms%n", &ln, &matched) >= 1 && matched > 0) {
            natoms = (long)ln;
        } else if ((matched = -1, sscanf(tmp, "%d atom types%n", &ni, &matched) >= 1) && matched > 0) {
            ntypes = ni;
        } else if ((matched = -1, sscanf(tmp, "%lf %lf xlo xhi%n", &v1, &v2, &matched) >= 2) && matched > 0) {
            lo[0] = v1; hi[0] = v2;
        } else if ((matched = -1, sscanf(tmp, "%lf %lf ylo yhi%n", &v1, &v2, &matched) >= 2) && matched > 0) {
            lo[1] = v1; hi[1] = v2;
        } else if ((matched = -1, sscanf(tmp, "%lf %lf zlo zhi%n", &v1, &v2, &matched) >= 2) && matched > 0) {
            lo[2] = v1; hi[2] = v2;
        } else if ((matched = -1, sscanf(tmp, "%lf %lf %lf xy xz yz%n", &v1, &v2, &v3, &matched) >= 3) && matched > 0) {
            tilt[0] = v1; tilt[1] = v2; tilt[2] = v3;
            has_tilt = 1;
        }

        p = next_line(le, end);
    }

    if (!atoms_section || natoms <= 0) {
        fprintf(stderr, "afc: no Atoms section or invalid atom count in %s\n", path);
        free(fb.data); return -1;
    }

    (void)has_tilt;
    box_from_lammps(&cfg->box, lo, hi, tilt);
    cfg->typemap.ntypes = ntypes;

    /* Allocate */
    atoms_reserve(&cfg->atoms, (int)natoms);
    cfg->atoms.count = (int)natoms;

    /* Index data lines */
    size_t *line_off;
    int nlines = index_data_lines(atoms_section, (size_t)(end - atoms_section),
                                  &line_off, (int)natoms);
    if (nlines < (int)natoms) {
        fprintf(stderr, "afc: only %d atom lines, expected %ld in %s\n",
                nlines, natoms, path);
        free(line_off); free(fb.data); return -1;
    }

    /* Parallel parse: id type x y z [ix iy iz] */
    const char *base = atoms_section;
    double Lx = cfg->box.cell[0][0];
    double Ly = cfg->box.cell[1][1];
    double Lz = cfg->box.cell[2][2];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(natoms > 10000)
#endif
    for (int i = 0; i < (int)natoms; i++) {
        const char *lp = base + line_off[i];
        const char *le = skip_to_eol(lp, end);
        lp = skip_ws(lp, le);

        char *ep;
        int id   = (int)strtol(lp, &ep, 10); lp = ep;
        int type = (int)strtol(lp, &ep, 10); lp = ep;
        double x = strtod(lp, &ep);           lp = ep;
        double y = strtod(lp, &ep);           lp = ep;
        double z = strtod(lp, &ep);           lp = ep;

        /* Check for image flags */
        while (lp < le && (*lp == ' ' || *lp == '\t')) lp++;
        if (lp < le && (*lp == '-' || (*lp >= '0' && *lp <= '9'))) {
            int ix = (int)strtol(lp, &ep, 10); lp = ep;
            int iy = (int)strtol(lp, &ep, 10); lp = ep;
            int iz = (int)strtol(lp, &ep, 10);
            /* unwrap */
            x += ix * Lx;
            y += iy * Ly;
            z += iz * Lz;
        }

        cfg->atoms.id[i]   = id;
        cfg->atoms.type[i] = type;
        cfg->atoms.x[i]    = x;
        cfg->atoms.y[i]    = y;
        cfg->atoms.z[i]    = z;
    }

    free(line_off);
    free(fb.data);
    return 0;
}

/* ================================================================== */
/* LAMMPS data writer (atomic style)                                   */
/* ================================================================== */

int write_lammps(const Config *cfg, const char *path)
{
    WBuf w;
    if (wbuf_open(&w, path) != 0) {
        fprintf(stderr, "afc: cannot open %s for writing\n", path);
        return -1;
    }

    int n = cfg->atoms.count;
    double lo[3], hi[3], tilt[3];
    box_to_lammps(&cfg->box, lo, hi, tilt);
    int triclinic = (fabs(tilt[0]) > 1e-12 ||
                     fabs(tilt[1]) > 1e-12 ||
                     fabs(tilt[2]) > 1e-12);

    /* Title */
    if (cfg->comment[0])
        wbuf_printf(&w, "%s\n", cfg->comment);
    else
        wbuf_puts(&w, "LAMMPS data file via afc\n");

    wbuf_printf(&w, "\n%d atoms\n", n);
    wbuf_printf(&w, "%d atom types\n", cfg->typemap.ntypes > 0 ?
                cfg->typemap.ntypes : 1);

    wbuf_printf(&w, "\n%.15g %.15g xlo xhi\n", lo[0], hi[0]);
    wbuf_printf(&w, "%.15g %.15g ylo yhi\n",   lo[1], hi[1]);
    wbuf_printf(&w, "%.15g %.15g zlo zhi\n",   lo[2], hi[2]);
    if (triclinic)
        wbuf_printf(&w, "%.15g %.15g %.15g xy xz yz\n",
                    tilt[0], tilt[1], tilt[2]);

    wbuf_puts(&w, "\nAtoms # atomic\n\n");

    /* Atom lines */
    int nth = 1;
#ifdef _OPENMP
    nth = omp_get_max_threads();
    if (nth > 64) nth = 64;
#endif

    if (n < 10000 || nth < 2) {
        for (int i = 0; i < n; i++) {
            wbuf_printf(&w, "%d %d %.15g %.15g %.15g\n",
                cfg->atoms.id[i], cfg->atoms.type[i],
                cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
        }
    } else {
        size_t est = 100;
        char **tbufs = (char **)calloc((size_t)nth, sizeof(char *));
        size_t *tlens = (size_t *)calloc((size_t)nth, sizeof(size_t));

        for (int t = 0; t < nth; t++) {
            int lo_i = (n * t) / nth;
            int hi_i = (n * (t + 1)) / nth;
            tbufs[t] = (char *)malloc((size_t)(hi_i - lo_i) * est);
            if (!tbufs[t]) die("oom");
        }

#pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            int lo_i = (n * tid) / nth;
            int hi_i = (n * (tid + 1)) / nth;
            size_t blen = 0;
            size_t bcap = (size_t)(hi_i - lo_i) * est;

            for (int i = lo_i; i < hi_i; i++) {
                if (blen + 256 > bcap) {
                    bcap *= 2;
                    tbufs[tid] = (char *)realloc(tbufs[tid], bcap);
                    if (!tbufs[tid]) die("oom");
                }
                int wrote = snprintf(tbufs[tid] + blen, bcap - blen,
                    "%d %d %.15g %.15g %.15g\n",
                    cfg->atoms.id[i], cfg->atoms.type[i],
                    cfg->atoms.x[i], cfg->atoms.y[i], cfg->atoms.z[i]);
                blen += (size_t)wrote;
            }
            tlens[tid] = blen;
        }

        for (int t = 0; t < nth; t++) {
            wbuf_write(&w, tbufs[t], tlens[t]);
            free(tbufs[t]);
        }
        free(tbufs);
        free(tlens);
    }

    return wbuf_close(&w);
}

/* ================================================================== */
/* CEL reader                                                          */
/* ================================================================== */

int read_cel(Config *cfg, const char *path)
{
    FileBuf fb;
    if (file_read_all(path, &fb) != 0) return -1;

    const char *p = fb.data, *end = fb.data + fb.len;

    /* Line 1: comment */
    const char *le = skip_to_eol(p, end);
    size_t clen = (size_t)(le - p);
    if (clen >= sizeof(cfg->comment)) clen = sizeof(cfg->comment) - 1;
    memcpy(cfg->comment, p, clen);
    cfg->comment[clen] = 0;
    if (clen > 0 && cfg->comment[clen - 1] == '\r')
        cfg->comment[--clen] = 0;
    /* Strip leading "# " from comment */
    p = next_line(le, end);

    /* Line 2: 0 a b c alpha beta gamma (nm, degrees) */
    le = skip_to_eol(p, end);
    {
        const char *q = skip_ws(p, le);
        char *ep;
        /* skip the leading 0 */
        strtol(q, &ep, 10); q = ep;
        double a_nm     = strtod(q, &ep); q = ep;
        double b_nm     = strtod(q, &ep); q = ep;
        double c_nm     = strtod(q, &ep); q = ep;
        double alpha    = strtod(q, &ep); q = ep;
        double beta     = strtod(q, &ep); q = ep;
        double gamma_d  = strtod(q, &ep);

        /* nm -> Angstroms */
        box_from_cellparams(&cfg->box, a_nm * 10.0, b_nm * 10.0, c_nm * 10.0,
                            alpha, beta, gamma_d);
    }
    p = next_line(le, end);

    /* Count atoms (lines until '*') */
    int natoms = 0;
    {
        const char *scan = p;
        while (scan < end) {
            const char *sl = skip_ws(scan, end);
            if (sl < end && *sl == '*') break;
            const char *sle = skip_to_eol(scan, end);
            const char *sq = skip_ws(scan, sle);
            if (sq < sle && *sq != '#' && *sq != '\r' && *sq != '\n')
                natoms++;
            scan = (sle < end) ? sle + 1 : end;
        }
    }

    atoms_reserve(&cfg->atoms, natoms);
    cfg->atoms.count = natoms;

    /* Temp arrays for fractional coords */
    double *fx = (double *)malloc((size_t)natoms * sizeof(double));
    double *fy = (double *)malloc((size_t)natoms * sizeof(double));
    double *fz = (double *)malloc((size_t)natoms * sizeof(double));
    if (!fx || !fy || !fz) die("oom");

    /* First pass: read symbols and build type map, store fractional coords */
    cfg->typemap.ntypes = 0;
    int idx = 0;
    while (p < end && idx < natoms) {
        le = skip_to_eol(p, end);
        const char *q = skip_ws(p, le);
        if (q >= le || *q == '*' || *q == '#' || *q == '\r' || *q == '\n') {
            if (q < le && *q == '*') break;
            p = next_line(le, end);
            continue;
        }

        /* Parse symbol */
        char sym[AFC_MAX_SYMBOL] = {0};
        int si = 0;
        while (q < le && !(*q == ' ' || *q == '\t') && si < AFC_MAX_SYMBOL - 1)
            sym[si++] = *q++;
        sym[si] = 0;

        /* Look up or assign type */
        int type = 0;
        for (int t = 1; t <= cfg->typemap.ntypes; t++) {
            if (strcmp(cfg->typemap.symbols[t], sym) == 0) { type = t; break; }
        }
        if (type == 0) {
            type = ++cfg->typemap.ntypes;
            if (type >= AFC_MAX_TYPES) die("too many atom types");
            strcpy(cfg->typemap.symbols[type], sym);
        }

        /* Parse fractional coords */
        char *ep;
        q = skip_ws(q, le);
        double f0 = strtod(q, &ep); q = ep;
        double f1 = strtod(q, &ep); q = ep;
        double f2 = strtod(q, &ep);
        /* remaining fields (occ, Biso, 0,0,0) are ignored on read */

        cfg->atoms.id[idx]   = idx + 1;
        cfg->atoms.type[idx] = type;
        fx[idx] = f0;
        fy[idx] = f1;
        fz[idx] = f2;
        idx++;

        p = next_line(le, end);
    }
    cfg->atoms.count = idx;

    /* Convert fractional -> Cartesian */
    frac_to_cart(&cfg->box, idx, fx, fy, fz,
                 cfg->atoms.x, cfg->atoms.y, cfg->atoms.z);

    free(fx); free(fy); free(fz);
    free(fb.data);
    return 0;
}

/* ================================================================== */
/* CEL writer                                                          */
/* ================================================================== */

int write_cel(const Config *cfg, const char *path)
{
    /* Check that we have element symbols */
    for (int t = 1; t <= cfg->typemap.ntypes; t++) {
        if (!cfg->typemap.symbols[t][0]) {
            fprintf(stderr,
                "afc: type %d has no element symbol. "
                "Use -t to specify a type map (e.g. -t 1=Ni,2=Al).\n", t);
            return -1;
        }
    }

    WBuf w;
    if (wbuf_open(&w, path) != 0) {
        fprintf(stderr, "afc: cannot open %s for writing\n", path);
        return -1;
    }

    int n = cfg->atoms.count;

    /* Comment line */
    if (cfg->comment[0])
        wbuf_printf(&w, "# %s\n", cfg->comment);
    else
        wbuf_puts(&w, "# Generated by afc\n");

    /* Cell parameters (Angstroms -> nm) */
    double a, b, c, alpha, beta, gamma;
    box_to_cellparams(&cfg->box, &a, &b, &c, &alpha, &beta, &gamma);
    wbuf_printf(&w, "     0 %10.4f %10.4f %10.4f %8.4f %8.4f %8.4f\n",
                a / 10.0, b / 10.0, c / 10.0, alpha, beta, gamma);

    /* Convert Cartesian -> fractional */
    double *fx = (double *)malloc((size_t)n * sizeof(double));
    double *fy = (double *)malloc((size_t)n * sizeof(double));
    double *fz = (double *)malloc((size_t)n * sizeof(double));
    if (!fx || !fy || !fz) die("oom");

    cart_to_frac(&cfg->box, n, cfg->atoms.x, cfg->atoms.y, cfg->atoms.z,
                 fx, fy, fz);

    /* Write atom lines */
    for (int i = 0; i < n; i++) {
        const char *sym = cfg->typemap.symbols[cfg->atoms.type[i]];
        wbuf_printf(&w, "  %-4s %10.6f %10.6f %10.6f %10.6f %10.6f "
                    "%10.6f %10.6f %10.6f\n",
                    sym, fx[i], fy[i], fz[i],
                    1.0, 0.0, 0.0, 0.0, 0.0);
    }

    wbuf_puts(&w, "   *\n");

    free(fx); free(fy); free(fz);
    return wbuf_close(&w);
}

/* ================================================================== */
/* Dispatch                                                            */
/* ================================================================== */

int config_read(Config *cfg, const char *path, Format fmt)
{
    switch (fmt) {
    case FMT_XYZ:    return read_xyz(cfg, path);
    case FMT_LAMMPS: return read_lammps(cfg, path);
    case FMT_CEL:    return read_cel(cfg, path);
    default:
        fprintf(stderr, "afc: unknown input format\n");
        return -1;
    }
}

int config_write(const Config *cfg, const char *path, Format fmt)
{
    switch (fmt) {
    case FMT_XYZ:    return write_xyz(cfg, path);
    case FMT_LAMMPS: return write_lammps(cfg, path);
    case FMT_CEL:    return write_cel(cfg, path);
    default:
        fprintf(stderr, "afc: unknown output format\n");
        return -1;
    }
}
