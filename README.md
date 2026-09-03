# LEGO-TOOLS

A collection of command-line tools for manipulating and converting atomistic simulation data.
The suite consists of three components: **lego-tools**, a set of utilities for filtering,
transforming, and inspecting atomic configurations; **afc** (Atomic Format Converter),
a tool for converting between common atomistic file formats; and **lego-dislo**, a tool
that inserts dislocations via the anisotropic elastic displacement field (sextic method).

All tools are written in C (C11) with optional OpenMP parallelism and transparent gzip
support via zlib.  Scripts in `lego-tools/scripts/` provide supplementary functionality
in Python 3 and AWK.

---

## Author and Acknowledgements

**Erik Bitzek**  
Max-Planck-Institut für Nachhaltige Materialien, Düsseldorf, Germany  
<e.bitzek@mpi-susmat.de>

Funded by **NFDI-MatWerk** (National Research Data Infrastructure for Materials Science
and Engineering Workflows).

---

## License

This software is released under the **GNU General Public License v3.0** (GPL-3.0).
You are free to use, study, and redistribute it. Modifications and derivative works
must also be distributed under GPL-3.0, ensuring that improvements flow back to the
community. See [LICENSE](LICENSE) for the full text.

When using LEGO-TOOLS in published work, please cite this repository and give credit
to the original author.

---

## Repository Layout

```
LEGO-TOOLS/
├── Makefile                  Top-level build entry point
├── README.md
├── LICENSE
├── lego-tools/               Atomistic manipulation utilities
│   ├── Makefile
│   ├── atomio.c / atomio.h   Shared I/O library (LAMMPS data/dump, IMD)
│   ├── cut.c                 → lego-cut
│   ├── shift.c               → lego-shift
│   ├── analyze.c             → lego-analyze
│   ├── remove_per_atom.c     → lego-remove-per-atom
│   ├── change_box.c          → lego-change-box
│   ├── nearest_atoms.c       → lego-nearest-atoms
│   ├── pbc_wrap.c            → lego-pbc-wrap
│   ├── make_box.c            → lego-make-box
│   └── scripts/
│       ├── cut.awk           AWK filter for IMD-format coordinate cutting
│       ├── lmps_get_max.awk  AWK script to report coordinate extremes (LAMMPS dump)
│       └── shift_lmps.py     Python script to shift LAMMPS data files
├── afc/                      Atomic Format Converter
│   ├── Makefile
│   ├── main.c                CLI entry point → afc
│   ├── converter.c / .h      Core conversion library
│   └── elements.h            Element symbol table (Z = 1–118)
└── dislo/                    Anisotropic elastic dislocation insertion
    ├── Makefile
    ├── README.md             Full documentation (theory, card format, API)
    ├── TOOLCARD.md           One-page reference (for humans and LLM agents)
    ├── main.c                CLI entry point → lego-dislo
    ├── card.c / .h           Input-card parser
    ├── sextic.c / .h         Sextic solver + displacement field
    ├── findroots.c           Laguerre polynomial root finder
    ├── python/legodislo/     Python driver (stdlib-only)
    ├── examples/             Runnable edge/screw/dipole examples
    └── tests/                Unit tests, golden validation, benchmark
```

---

## Dependencies

| Dependency | Purpose | Notes |
|---|---|---|
| C11 compiler (GCC or Clang) | Compilation | GCC ≥ 7 or Clang ≥ 5 recommended |
| zlib (`-lz`) | Transparent gzip I/O | Usually pre-installed; `apt install zlib1g-dev` on Debian/Ubuntu |
| OpenMP | Parallel parsing (optional) | `libgomp` on Linux; `brew install libomp` on macOS |
| Python 3 | `shift_lmps.py` script | Standard library only |
| AWK | AWK scripts | Any POSIX-compliant AWK |

---

## Compiling

### Build everything (recommended)

```bash
cd LEGO-TOOLS
make
```

Binaries are placed inside their respective subdirectories (`lego-tools/` and `afc/`).
Copy or symlink them to a directory on your `$PATH` as needed.

### Build without OpenMP

If OpenMP is not available on your system:

```bash
make no-omp
```

### Build individual components

```bash
make lego-tools   # build only the lego-tools suite
make afc          # build only the Atomic Format Converter
make dislo        # build only lego-dislo (also builds atomio)
```

### macOS (Apple Silicon / Homebrew)

Install `libomp` before building:

```bash
brew install libomp
make
```

The Makefile auto-detects the Homebrew prefix for libomp on Darwin.

### Clean build artefacts

```bash
make clean
```

---

## Tool Reference

### lego-tools suite

All lego-tools utilities read LAMMPS data files, LAMMPS dump files, and IMD ASCII
checkpoint files.  Input format is auto-detected; output format defaults to LAMMPS data
and can be overridden with `--out-format=lammps|lammps-dump|imd`.
Gzipped inputs and outputs (`.gz`) are handled transparently.

---

#### `lego-cut`

Keep only atoms whose coordinates fall inside a specified axis-aligned bounding box.

```
lego-cut [options] xmin xmax ymin ymax zmin zmax INPUT OUTPUT

Options:
  --in-format=lammps|lammps-dump|imd   (default: auto-detect)
  --out-format=lammps|lammps-dump|imd  (default: lammps)
  -h, --help
```

The simulation box bounds are preserved unchanged; only the atom list is filtered.

---

#### `lego-shift`

Translate all atom coordinates and the simulation box by a vector (dx, dy, dz).

```
lego-shift [options] dx dy dz INPUT OUTPUT

Options:
  --out-format=lammps|lammps-dump|imd  (default: lammps)
  -h, --help
```

---

#### `lego-analyze`

Print a summary of a configuration file: detected format, atom count, type counts,
box geometry, per-column min/max values, and whether velocities are present.

```
lego-analyze [--json] INPUT
```

The optional `--json` flag produces machine-readable output.

---

#### `lego-remove-per-atom`

Strip all per-atom data columns except `id`, `type`, `x`, `y`, `z`.
Useful for removing velocities, forces, charges, image flags, etc. before passing
a file to tools that expect a minimal atom style.

```
lego-remove-per-atom [options] INPUT OUTPUT

Options:
  --out-format=lammps|lammps-dump|imd  (default: lammps)
  -h, --help
```

---

#### `lego-change-box`

Modify the simulation box bounds without parsing or rewriting the atom data.
Two modes:

```
# Expand/shrink the upper bound by (dx, dy, dz):
lego-change-box dx dy dz INPUT OUTPUT

# Set explicit bounds:
lego-change-box xlo xhi ylo yhi zlo zhi INPUT OUTPUT
```

Orthogonal boxes only. Streams the file line-by-line — fast even on large gzipped files.

---

#### `lego-nearest-atoms`

Find the k atoms nearest to a query point, with optional periodic boundary conditions.

```
lego-nearest-atoms [options] x y z INPUT

Options:
  -k N         number of nearest atoms to report (default: 5)
  --tol=EPS    distance ≤ EPS is treated as "on site" (default: 1e-6)
  --no-pbc     ignore PBC (default: minimum-image convention)
  --json       machine-readable output
  -h, --help
```

---

#### `lego-pbc-wrap`

Wrap atoms back into the primary simulation box using periodic boundary conditions.
If image-flag columns (ix, iy, iz) are present they are updated consistently so that
the unwrapped position is preserved.

```
lego-pbc-wrap [options] INPUT OUTPUT

Options:
  --out-format=lammps|lammps-dump  (default: lammps)
  -x, -y, -z                      wrap only the specified axis (default: all three)
  -h, --help
```

LAMMPS format only (orthogonal and triclinic boxes rejected for triclinic).

---

#### `lego-make-box`

Ensure every atom lies inside the simulation box.  First wraps atoms via PBC for
periodic LAMMPS configurations; then, if any atom is still outside, expands the box.
The box is never shrunk.

```
lego-make-box [options] INPUT OUTPUT

Options:
  --out-format=lammps|lammps-dump|imd  (default: lammps)
  --no-wrap                            skip the PBC wrapping step
  --pad=EPS                            extra padding added when expanding (default: 0)
  -h, --help
```

---

### Scripts

#### `scripts/shift_lmps.py`

Python alternative to `lego-shift` that works directly on raw LAMMPS data files
(plain text, not gzipped) and supports triclinic boxes.

```
python3 shift_lmps.py [options] dx dy dz INPUT OUTPUT

Options:
  --xyz-cols XCOL YCOL ZCOL   0-based column indices for x,y,z (default: 2 3 4)
  --float-format FMT          Python format spec for output floats (default: .16g)
```

#### `scripts/cut.awk`

AWK script for cutting atoms by coordinate range from IMD-format files.  Takes bounds
as positional arguments:

```
awk -f cut.awk xmin xmax ymin ymax zmin zmax INPUT
```

#### `scripts/lmps_get_max.awk`

AWK script that reports coordinate extremes from a LAMMPS dump file
(format: `id type x y z vx vy vz`).

```
awk -f lmps_get_max.awk DUMP_FILE
```

---

### afc — Atomic Format Converter

Convert atomistic configuration files between Extended XYZ, LAMMPS data (atomic style),
and CEL (Dr Probe supercell) formats.  Gzipped files are handled transparently.

```
afc [options] INPUT OUTPUT

Options:
  -i FMT   Input format:  xyz, lammps, cel  (auto-detected from extension)
  -o FMT   Output format: xyz, lammps, cel  (auto-detected from extension)
  -t MAP   Type map, e.g. "1=Ni,2=Al" or "Ni,Al"
  -h       Show help

Supported file extensions:
  .xyz[.gz]              Extended XYZ
  .lmp[.gz] .data[.gz]  LAMMPS data (atomic style)
  .cel[.gz]              CEL (Dr Probe supercell)
```

The `-t` type map is required when converting to CEL from formats that use
numeric type IDs (LAMMPS).  Examples:

```bash
afc structure.lmp structure.xyz
afc -t Ni grain.lmp grain.cel
afc -t "1=Ni,2=Al" alloy.data alloy.xyz
afc input.xyz.gz output.lmp
```

---

### `lego-dislo` — anisotropic elastic dislocation insertion

Insert one or many dislocations into an atomistic configuration (up to ~10⁸ atoms)
by adding the anisotropic elastic displacement field of each dislocation (sextic
method, Hirth & Lothe, *Theory of Dislocations*). A faithful, OpenMP-parallel port
of the 2005 Fortran toolkit *Disloelast*, driven by a single small card file
(elastic constants, box orientation, and per-dislocation line direction, Burgers
vector, and glide plane, all as crystallographic directions).

```
lego-dislo [options] CARD INPUT OUTPUT
lego-dislo --solve-only CARD             # solver report only, no atoms

Options:
  --in-format=lammps|lammps-dump|imd   (default: auto-detect)
  --out-format=lammps|lammps-dump|imd  (default: lammps)
  --threads=N                          (default: OMP_NUM_THREADS)
  --report=FILE                        write P/F/G and K factors to FILE
  -h, --help
```

Multiple dislocations superpose linearly (no periodic-image corrections); box
bounds are left unchanged — cut/shift afterwards with `lego-cut` / `lego-shift`.
A stdlib-only Python driver lives in `dislo/python/legodislo`. See
[`dislo/README.md`](dislo/README.md) for the card format, theory, validation
against the original 2005 golden outputs, and worked examples
(`dislo/examples/`); [`dislo/TOOLCARD.md`](dislo/TOOLCARD.md) is a one-page
compact reference.

---

## Contributing

Contributions are welcome — bug fixes, new tools, additional format support, and
documentation improvements.  Please:

1. Fork the repository on GitHub.
2. Create a feature branch (`git checkout -b feature/my-improvement`).
3. Commit your changes with a clear message.
4. Open a pull request describing what you changed and why.

All contributed code must be compatible with GPL-3.0.

---

## Citation

If LEGO-TOOLS contributes to published research, please cite it as:

> Erik Bitzek, *LEGO-TOOLS*, Max-Planck-Institut für Nachhaltige Materialien,
> Düsseldorf, Germany. https://github.com/biterik/LEGO-TOOLS

Funding acknowledgement: *This work was supported by NFDI-MatWerk.*
