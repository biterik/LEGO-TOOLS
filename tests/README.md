# tests/

End-to-end smoke tests for the merged LEGO-TOOLS build.

## `smoke.sh`

Exercises `afc` and a representative subset of `lego-tools` binaries on the
bundled Ni test structure (`afc/test-data/relax_d0_Ni_Mishin04_...xyz.gz`,
92,160 atoms) and checks that every stage preserves the expected atom count.

Run from the repo root:

```bash
make test
```

or directly:

```bash
./tests/smoke.sh
```

The script requires all binaries to have been built already (`make` at the
repo root). It exits non-zero on the first failure; individual checks are
printed with an `ok`/`FAIL` prefix.

Checks performed:

1. `afc` roundtrip: `xyz.gz -> lammps-data -> cel -> xyz`.
2. `lego-analyze --json` reports 92,160 atoms on the lammps-data file.
3. `lego-shift` preserves atom count.
4. `lego-pbc-wrap` preserves atom count.
5. `lego-remove-per-atom` preserves atom count.

## `large-file.sh`

Regression test for the >4 GiB file path in `atomio.c` **and**
`afc/converter.c`.  Zlib's legacy `gzread`/`gzwrite` can't transfer more
than ~2 GiB in a single call, so both I/O libraries chunk every call to
`GZ_CHUNK_MAX` bytes (default 1 GiB).  This script rebuilds *both*
`lego-tools` and `afc` with a tiny `GZ_CHUNK_MAX=4096`, which forces
*hundreds* of gzread/gzwrite calls per file, and asserts that the output
is byte-identical to a reference build at the default chunk size.  Also
checks that `LEGO_DEBUG=1` emits the expected diagnostic lines on both
code paths.

Run from the repo root:

```bash
make test-large
```

The script restores a clean default-chunk build on exit.

## Debug output

All `lego-tools` binaries honour the `LEGO_DEBUG=1` environment variable
and print extra progress/diagnostic lines on stderr, e.g.:

```text
atomio: opening input.fcc (file size: 7.12 GiB (7648102912 bytes), chunk: 1073741824)
atomio: read progress on input.fcc: 1.00 GiB (1073741824 bytes)
atomio: read progress on input.fcc: 2.00 GiB (2147483648 bytes)
...
atomio: finished reading input.fcc: 7.12 GiB (7648102912 bytes)
atomio: detected format: imd
atomio: parsed 227829632 atoms, 5 columns, ntypes=1
```

Set `LEGO_DEBUG=1` when diagnosing failures on unusually large files.
