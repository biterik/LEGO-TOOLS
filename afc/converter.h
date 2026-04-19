/* converter.h — Atomic Format Converter public API.
 *
 * Converts between Extended XYZ, LAMMPS data (atomic style), and
 * CEL (Dr Probe supercell) formats.  Transparent gzip support via zlib.
 * OpenMP-parallel parsing for large files.
 *
 * Internal unit: Angstroms, Cartesian coordinates.
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#ifndef CONVERTER_H
#define CONVERTER_H

#include <stddef.h>

/* ------------------------------------------------------------------ */
/* Constants                                                           */
/* ------------------------------------------------------------------ */

#define AFC_MAX_SYMBOL  4   /* "Ni\0" fits in 4 bytes */
#define AFC_MAX_TYPES 128

/* ------------------------------------------------------------------ */
/* Format enum                                                         */
/* ------------------------------------------------------------------ */

typedef enum {
    FMT_UNKNOWN = 0,
    FMT_XYZ,
    FMT_LAMMPS,
    FMT_CEL
} Format;

/* ------------------------------------------------------------------ */
/* Box — 3x3 lattice matrix (row-major) + origin                       */
/* ------------------------------------------------------------------ */

typedef struct {
    double cell[3][3];   /* cell[i] = i-th lattice vector */
    double origin[3];    /* box origin offset */
} Box;

/* ------------------------------------------------------------------ */
/* Atoms — Structure-of-Arrays for cache-friendly access               */
/* ------------------------------------------------------------------ */

typedef struct {
    int    *id;          /* 1-based atom IDs */
    int    *type;        /* 1-based numeric type */
    double *x, *y, *z;  /* Cartesian positions (Angstroms) */
    int     count;
    int     capacity;
} Atoms;

/* ------------------------------------------------------------------ */
/* TypeMap — numeric type ↔ element symbol                             */
/* ------------------------------------------------------------------ */

typedef struct {
    int  ntypes;
    char symbols[AFC_MAX_TYPES][AFC_MAX_SYMBOL]; /* 1-indexed */
} TypeMap;

/* ------------------------------------------------------------------ */
/* Config — top-level container                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    Box     box;
    Atoms   atoms;
    TypeMap typemap;
    char    comment[256];
} Config;

/* ------------------------------------------------------------------ */
/* Lifecycle                                                           */
/* ------------------------------------------------------------------ */

Config config_create(void);
void   config_free(Config *cfg);
void   atoms_reserve(Atoms *a, int capacity);

/* ------------------------------------------------------------------ */
/* Format helpers                                                      */
/* ------------------------------------------------------------------ */

Format      format_from_extension(const char *path);
Format      format_from_string(const char *name);
const char *format_name(Format f);

/* ------------------------------------------------------------------ */
/* I/O — gz-transparent.  Returns 0 on success.                        */
/* ------------------------------------------------------------------ */

int config_read (Config *cfg, const char *path, Format fmt);
int config_write(const Config *cfg, const char *path, Format fmt);

/* Individual readers / writers */
int read_xyz   (Config *cfg, const char *path);
int read_lammps(Config *cfg, const char *path);
int read_cel   (Config *cfg, const char *path);

int write_xyz   (const Config *cfg, const char *path);
int write_lammps(const Config *cfg, const char *path);
int write_cel   (const Config *cfg, const char *path);

/* ------------------------------------------------------------------ */
/* Box utilities                                                       */
/* ------------------------------------------------------------------ */

void box_to_lammps(const Box *box, double lo[3], double hi[3], double tilt[3]);
void box_from_lammps(Box *box, const double lo[3], const double hi[3],
                     const double tilt[3]);

void box_to_cellparams(const Box *box, double *a, double *b, double *c,
                       double *alpha, double *beta, double *gamma);
void box_from_cellparams(Box *box, double a, double b, double c,
                         double alpha, double beta, double gamma);

int  box_invert(const Box *box, double inv[3][3]);

/* ------------------------------------------------------------------ */
/* Coordinate transforms (OpenMP-parallel over atoms)                  */
/* ------------------------------------------------------------------ */

void cart_to_frac(const Box *box, int n,
                  const double *cx, const double *cy, const double *cz,
                  double *fx, double *fy, double *fz);

void frac_to_cart(const Box *box, int n,
                  const double *fx, const double *fy, const double *fz,
                  double *cx, double *cy, double *cz);

#endif /* CONVERTER_H */
