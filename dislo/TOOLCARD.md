# lego-dislo — tool card

Insert one or more dislocations into an atomistic configuration (LAMMPS
data, LAMMPS dump, or IMD; gzip transparent) by adding the anisotropic
elastic displacement field (sextic method, Hirth & Lothe). Multiple
dislocations superpose linearly; box bounds are left unchanged (cut/shift
afterwards with `lego-cut` / `lego-shift`). No periodic-image corrections.

## CLI

    lego-dislo [options] CARD INPUT OUTPUT
    lego-dislo --solve-only CARD             # report only, no atoms

    --in-format=lammps|lammps-dump|imd   default: auto-detect
    --out-format=lammps|lammps-dump|imd  default: lammps
    --threads=N                          default: OMP_NUM_THREADS
    --report=FILE                        also write the report to FILE

The report (stdout and `--report`) lists per dislocation the sextic roots
`P_n`, coefficient matrices `F_k`/`G_k`, and the energy prefactor:
`K_b2_over_4pi_J_per_m`, and `K_GPa` (energy factor K).

## Card file (complete example; this is the entire schema)

```yaml
material:
  name: W
  lattice_parameter: 3.1652          # Angstrom
  cij:                               # GPa; EITHER cubic shorthand ...
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
    # ... OR the full 6x6 Voigt matrix in the crystal frame:
    # full:
    #   - [522.4, 204.4, 204.4, 0, 0, 0]
    #   - ... (6 rows: c[1][1..6] .. c[6][1..6])
box_orientation:                     # crystallographic direction of each
  x: [1, 1, 1]                       # box axis; must be mutually orthogonal
  y: [-1, 0, 1]                      # and right-handed (z = x cross y)
  z: [1, -2, 1]
dislocations:                        # one or more entries
  - point: [0.0, 0.0, 0.0]           # point on the line, box coords, A
    line: [1, -2, 1]                 # line direction, crystal coords
    burgers_direction: [-1, -1, -1]  # crystal coords
    burgers_magnitude: 0.8663        # units of lattice_parameter
                                     # (length in A = value * a0)
    glide_normal: [-1, 0, 1]         # crystal coords; must be
                                     # perpendicular to line
    cut_angle: 0.0                   # optional, degrees, default 0
```

Only this flat two-level YAML subset is understood (`#` comments, 2-space
indents, inline `[lists]` / `{maps}`; no tabs, no anchors, no nesting
beyond what is shown). The displacement discontinuity (branch cut) lies in
the half-plane at `cut_angle` degrees from the local +X = `glide_normal x
line` direction; atoms exactly on the cut belong to its lower side.

## Python (stdlib-only; module at dislo/python/legodislo)

```python
from legodislo import Material, Dislocation, insert_dislocations
result = insert_dislocations(
    input="W.data", output="W_edge.data",
    material=Material(name="W", a0=3.1652, cubic=(522.4, 204.4, 160.6)),
    box_orientation=((1, 1, 1), (-1, 0, 1), (1, -2, 1)),
    dislocations=[Dislocation(point=(0, 0, 0), line=(1, -2, 1),
                              b_dir=(-1, -1, -1), b_mag=0.8663,
                              glide_normal=(-1, 0, 1))],
    binary=None)                     # auto-locates lego-dislo; or a path
print(result.dislocations[0].K_GPa)  # also .P/.F/.G, .K_b2_over_4pi_J_per_m
```

Raises `LegoDisloError` on failure. `solve(...)` runs the solver without
atoms. CLI twin: `python -m legodislo CARD INPUT OUTPUT`.

## Gotchas

- `glide_normal . line = 0` is enforced (1e-8); box axes likewise.
- Mixed (edge/screw/any) character is fine: `burgers_direction` is free.
- ~10^8 atoms x 1 dislocation: seconds of field time (OpenMP); total run
  is dominated by file I/O.
