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
