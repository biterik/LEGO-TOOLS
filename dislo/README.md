# lego-dislo — anisotropic elastic dislocation insertion

`lego-dislo` inserts one or many dislocations into atomistic
configurations (up to ~10⁸ atoms) by adding the **anisotropic elastic
displacement field** of each dislocation to the atom positions. It is a
modern, faithful port of the 2005 Fortran toolkit *Disloelast*
(`elast-cubic.f` / `elastbart.f` + `anisodisloc.f`).

**Lineage:** the sextic solver goes back to Alexey Girshick (Department of
Materials Science and Engineering, University of Pennsylvania), with the
cubic-system extensions by Bart Pestman (1988) and the displacement
pipeline assembled by Erik Bitzek (2005). The method is described in
J. P. Hirth and J. Lothe, *Theory of Dislocations* (2nd ed., Wiley 1982),
ch. 13.

---

## Theory in one paragraph

For a straight dislocation along the local Z axis in a general anisotropic
linear-elastic medium, the equilibrium equations admit displacement
solutions of the form u ∝ f(x + p·y) with complex p. Requiring
non-trivial solutions leads to a **sextic** (degree-6) polynomial
det|a_ik(p)| = 0 whose six roots come in three complex-conjugate pairs
p(n) = PR(n) ± i·PP(n) (n = 1..3). The boundary conditions (Burgers
vector b, no net traction) fix complex coefficients D(n), and the
displacement of a point at in-plane position (x, y) relative to the core,
in lattice-parameter units, is

    u_k = -(1/2π) Σ_{n=1..3} [ F(k,n)·ln((x + PR(n)·y)² + (PP(n)·y)²)
                              + G(k,n)·arctg_n(x,y) ]

where arctg_n is the argument of the complex number
(x + PR(n)·y) + i·PP(n)·y, branch-corrected so that its single 2π jump —
the **cut**, where the displacement jumps by b — lies along the chosen
direction (`cut_angle`) in the (x, y) plane. The elastic energy per
length is (K·b²/4π)·ln(R/r₀); the prefactor K·b²/4π is reported. The
field is scale-invariant in the elastic constants (only K depends on the
absolute scale) and depends only on the in-plane coordinates.

## Frame conventions

Three frames are involved; all rotations are built from crystallographic
directions given in the card:

```
  crystal frame ── box_orientation ──▶ box frame        (atom coordinates)
       │
       └── line / glide_normal ──▶ dislocation frame:
               Z = line direction
               Y = glide-plane normal      (Y ⊥ Z enforced)
               X = Y × Z                   (right-handed)

                 Y (glide normal)
                 │      cut (cut_angle = 0): along +X
                 │    ─────────▶
                 ●───────────── X (glide direction)
                Z (line, out of the plane)
```

For every atom: position → shifted by `point` → rotated
box → crystal → dislocation frame → (x, y) in lattice-parameter units →
u evaluated → rotated back → **added** to the position. With several
dislocations the fields are evaluated independently, each in its own
frame, and superposed linearly.

## Input card

One small file drives a run (YAML-look syntax; see the grammar note
below):

```yaml
material:
  name: W
  lattice_parameter: 3.1652          # Angstrom
  cij:                               # GPa — either the cubic shorthand:
    cubic: {c11: 522.4, c12: 204.4, c44: 160.6}
    # ... or the full 6x6 Voigt matrix in the crystal frame:
    # full:
    #   - [522.4, 204.4, 204.4, 0, 0, 0]
    #   - [204.4, 522.4, 204.4, 0, 0, 0]
    #   - [204.4, 204.4, 522.4, 0, 0, 0]
    #   - [0, 0, 0, 160.6, 0, 0]
    #   - [0, 0, 0, 0, 160.6, 0]
    #   - [0, 0, 0, 0, 0, 160.6]
box_orientation:                     # crystallographic direction of each
  x: [1, 1, 1]                       # box axis; must be orthogonal and
  y: [-1, 0, 1]                      # right-handed (z = x × y), checked
  z: [1, -2, 1]                      # to 1e-8
dislocations:
  - point: [0.0, 0.0, 0.0]           # point on the line, box coords, A
    line: [1, -2, 1]                 # crystal coords
    burgers_direction: [-1, -1, -1]  # crystal coords (any character)
    burgers_magnitude: 0.8663        # in units of lattice_parameter
    glide_normal: [-1, 0, 1]         # crystal coords, ⊥ line (checked)
    cut_angle: 0.0                   # optional, degrees from local +X
```

**Units policy.** Coordinates and the lattice parameter are in Angstrom.
Elastic constants are in GPa (the displacement field only needs them
relatively; the K factor is reported in SI as well). The Burgers-vector
magnitude is in lattice-parameter units, exactly like the legacy inputs:
its length in Angstrom is `burgers_magnitude × lattice_parameter`
(e.g. ½⟨111⟩ in bcc: √3/2 = 0.8663).

**Card grammar** (hand-written parser, no external dependencies): the
flat two-level structure shown above is the *whole* language — top-level
sections `material` / `box_orientation` / `dislocations`, two-space
indentation, `key: value` pairs, inline lists `[a, b, c]`, one inline map
for `cubic:`, `-` list items, and `#` comments. Tabs, anchors, multi-line
strings, or deeper nesting are rejected. The Python driver only ever
emits this subset.

## CLI

```
lego-dislo [options] CARD INPUT OUTPUT
lego-dislo --solve-only [options] CARD

  --in-format=lammps|lammps-dump|imd    (default: auto-detect, via atomio)
  --out-format=lammps|lammps-dump|imd   (default: lammps)
  --threads=N                           (default: OMP_NUM_THREADS)
  --report=FILE                         (write the solver report to FILE)
  -h, --help
```

Formats and transparent gzip come from the shared `atomio` library
(LAMMPS data, LAMMPS dump, IMD ASCII). The solver report — per
dislocation the roots `P_n`, the coefficient matrices `F_k` (log terms)
and `G_k` (arctan terms), and the energy prefactor — is always printed to
stdout; `--report` additionally writes it to a file:

```
  K_b2_over_4pi_GPa_lp2: 1.3276411e+01     # K·b²/4π, GPa·(b in a0 units)²
  K_b2_over_4pi_J_per_m: 1.3300961e-09     # the same in J/m
  K_GPa: 2.2230741e+02                     # energy factor K
```

**Box bounds are left unchanged** (legacy behavior): the field displaces
atoms, some of which may end up outside the box. Cut or shift afterwards
with `lego-cut` / `lego-shift` — exactly like the 2005 workflow did with
its awk cutters.

## Python API

Standard-library-only driver in `python/legodislo` (add
`dislo/python` to `PYTHONPATH` or copy the package):

```python
from legodislo import Material, Dislocation, insert_dislocations

result = insert_dislocations(
    input="W.data", output="W_edge.data",
    material=Material(name="W", a0=3.1652, cubic=(522.4, 204.4, 160.6)),
    box_orientation=((1, 1, 1), (-1, 0, 1), (1, -2, 1)),
    dislocations=[Dislocation(point=(0, 0, 0), line=(1, -2, 1),
                              b_dir=(-1, -1, -1), b_mag=0.8663,
                              glide_normal=(-1, 0, 1))],
    binary=None)          # auto-locates lego-dislo; or pass a path

d = result.dislocations[0]
print(d.K_GPa, d.K_b2_over_4pi_J_per_m)   # plus d.P, d.F, d.G
```

It writes the card to a temporary file, invokes the binary, checks the
return code, parses the report, and raises `LegoDisloError` with the full
diagnostics on any failure. `solve(...)` runs the solver without touching
atoms; `python -m legodislo CARD INPUT OUTPUT` is the CLI twin. The full
6×6 matrix goes in as `Material(a0=..., cij=[[...6 rows...]])`.

## Examples

See [`examples/`](examples/README.md): a W ½⟨111⟩ edge, a W ½⟨111⟩
screw, and a two-dislocation dipole, each one command
(`./run_example.sh edge|screw|dipole`) with committed expected outputs.

## Validation

The port is validated against the stored golden outputs of the original
2005 Fortran binaries (`dislo/tests/`, auto-skipped when the legacy
reference folder is absent):

- solver: every printed digit of the legacy `W_edge_results` /
  `W_screw_results` P/F/G tables (≲4e-7 relative), and
  K·b²/4π = 0.1327641e+11 (legacy units);
- field: max atom-position deviation 2.0e-6 Å over the 197200-atom edge
  reference and 2.7e-7 Å over the 330000-atom screw reference
  (the remaining difference stems from single-precision constants, e.g.
  π = 3.141592654, in the legacy code — this port is uniformly double
  precision);
- self-tests: rotation-equivalence across box frames to 1e-12 Å,
  superposition identities, cubic-vs-full-6×6 identity, parser and frame
  validation errors (`tests/run_tests.sh`).

## Performance

OpenMP-parallel over atoms; work is O(natoms × ndislocations) and the
solver cost is negligible. On an Apple-silicon Mac (M5, 18 cores),
10⁸ atoms × 1 dislocation take ~4.5 s of field time on one thread and
~0.4 s on 18 threads; total wall time is dominated by reading/writing the
multi-GB configuration files. `tests/bench.sh` reproduces the numbers.

## Limitations

- **Linear elasticity**: the field diverges at the core; positions very
  near the core (and on the cut) are as unphysical as they were in 2005.
  Relax the configuration afterwards.
- **No periodic-image corrections**: multiple dislocations superpose
  naively; fields do not satisfy PBC. This matches the original scope.
- The displacement discontinuity (cut) lies along `cut_angle` from the
  local +X axis; atoms exactly on the cut get the lower-side value.
- Box bounds are not updated (see above).

## License and attribution

GPL-3.0. Erik Bitzek, Max-Planck-Institut für Nachhaltige Materialien,
Düsseldorf. Funded by NFDI-MatWerk. Please cite this repository and
Hirth & Lothe when publishing results obtained with lego-dislo.
