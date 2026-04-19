/* atomio.h — shared I/O library for atomistic configurations.
 *
 * Supports:
 *   - LAMMPS data files  (header + Masses/Atoms/Velocities/... sections)
 *   - LAMMPS dump files  (ITEM: TIMESTEP / ATOMS / ...)
 *   - IMD ASCII checkpoint files (#F A ... / #C ... / #X / #Y / #Z / #E)
 *
 * Reading auto-detects the format. Writing takes an explicit target format.
 * Gzipped inputs and outputs are transparently supported via zlib (suffix .gz).
 *
 * Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
 *          Max-Planck-Institut fuer Nachhaltige Materialien,
 *          Duesseldorf, Germany
 * Funding: NFDI-MatWerk
 */
#ifndef ATOMIO_H
#define ATOMIO_H

#include <stddef.h>

typedef enum {
    FMT_UNKNOWN = 0,
    FMT_LAMMPS_DATA,
    FMT_LAMMPS_DUMP,
    FMT_IMD
} Format;

#define MAX_COLS 64
#define MAX_COL_NAME 32

typedef struct {
    Format format;

    /* Box */
    double xlo, xhi, ylo, yhi, zlo, zhi;
    int    triclinic;
    double xy, xz, yz;
    int    pbc[3];                    /* 1 = periodic (default) */

    /* Counts */
    size_t natoms;
    int    ntypes;

    /* Per-atom columns. Row-major: data[i*ncols + c]. */
    int     ncols;
    char    col_names[MAX_COLS][MAX_COL_NAME];
    int     col_is_int[MAX_COLS];     /* 1 => print as integer */
    double *data;

    /* Quick indices (-1 if not present). */
    int idx_id, idx_type, idx_mol, idx_q;
    int idx_x,  idx_y,   idx_z;
    int idx_vx, idx_vy,  idx_vz;
    int idx_ix, idx_iy,  idx_iz;

    int has_velocities;

    /* LAMMPS-specific context. */
    char atom_style[32];              /* e.g. "atomic", "charge", "full" */
    char *title;                      /* first line of a LAMMPS data file */
    char *extra_sections;             /* raw text of pass-through sections
                                       * (Masses, Bonds, Coeffs, …)        */

    /* IMD-specific context: count of "data" columns beyond vel in #F line. */
    int imd_ndata;
    int imd_has_mass;                 /* whether input IMD had a mass column */
} Config;

void   config_init (Config *c);
void   config_free (Config *c);

Format detect_format     (const char *filename);
const char *format_name  (Format f);
Format format_from_string(const char *s);

/* Read from file (format auto-detected). Returns 0 on success. */
int config_read (Config *c, const char *filename);

/* Write to file in the requested format. Returns 0 on success. */
int config_write(const Config *c, const char *filename, Format fmt);

/* Column helpers. */
int  config_find_col            (const Config *c, const char *name);
void config_refresh_std_indices (Config *c);

#endif /* ATOMIO_H */
