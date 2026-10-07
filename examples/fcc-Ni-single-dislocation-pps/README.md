# Example: one straight dislocation in an fcc Ni slab — screw, 30°, 60°, edge

The chain **LEGO → dcreator → lego-cut → box correction → LAMMPS** for a single a/2⟨110⟩{111}
dislocation in a slab that is periodic along the glide direction and the line, with free (111)
surfaces (LAMMPS `boundary p p s`). Every character is relaxed twice, with the box left as built
and with the corrected period vector, and the full stress tensor of both is plotted.

Frame everywhere: **x = glide direction** (periodic), **y = line direction** (periodic),
**z = glide-plane normal** (free surfaces; the 3 outermost (111) layers per face move only in
x and y, the 2D-dynamic boundary of Rodney & Martin 2000).

## The rule

dcreator shifts one half-crystal by the full Burgers vector **b** in the region between the line
and one x-face, and leaves the rest untouched. The upper and lower halves therefore want
different periods along x. The box's x-period vector has to be their **average**:

| | formula | screw | 30° | 60° | edge |
|---|---|---|---|---|---|
| length along the glide direction | Lx′ = L0 − b_e/2 | L0 | L0 − 0.622 Å | L0 − 1.078 Å | L0 − 1.245 Å |
| offset along the line | δ = ±b_s/2 | ±1.245 Å | ±1.078 Å | ±0.622 Å | 0 |

Here b_e = b sinθ and b_s = b cosθ are the edge and screw components (b = a0/√2 = 2.489 Å for
a0 = 3.52 Å). The signs depend on the construction; the gate script computes both numbers from
the actual displacement field.

After the correction the two halves carry **equal and opposite** strains of ±b_e/(2Lx) (normal)
and ±b_s/(2Lx) (shear). That is the unavoidable field of the periodic row of dislocations; its
mean is zero. This is Rodney's construction: Lx + b/2 for an edge built by adding half planes
(Rodney & Martin, PRB 61, 8714 (2000)), and "a shift of +b/2 in Y for atoms leaving through X"
for a screw (Rodney, Acta Mater. 52, 607 (2004)).

### Three implementation traps (all hit while making this example)

1. **The offset is along the line, across the x boundary.** LAMMPS can only tilt its second box
   vector, b = (xy, ly, 0): `xy` offsets *x* across the *y* boundary. So
   `change_box all xy final δ` shears the crystal by δ/Ly instead of δ/Lx. Here that was
   ~5 times too much and of the wrong kind: in the screw cell, σ_xy went from −0.64 to +2.7 GPa,
   and the energy rose. `apply_period_vector.py` therefore writes the corrected cell with the
   **line along LAMMPS x** (X = line, Y = −glide, Z = normal), where the period vector becomes
   the LAMMPS b vector (xy, ly) = (−δ, Lx′). The plotting script rotates the stresses back.
2. **Wrap triclinic coordinates through fractional coordinates.** Wrapping the y coordinate
   alone, without the matching xy shift of x, breaks the seam. `check_seam_coordination.py`
   checks the written file independently and caught exactly this.
3. **The rigid shift must push material out of the box, not in.** `rigid_side` is defined in
   the dislocation frame (x_d = n × l), not the sample frame. If the edge part of the shift
   points *into* the box, a gap of b_e opens at the seam. That gap is an opposite dislocation,
   and no box change repairs it. `build_single_dislocation_cell.py` picks the side so that the
   shifted atoms leave the box; `lego-cut 0 L0` then deletes exactly them, and the seam is a
   perfect lattice match. There is no cut value to choose.

`lego-change-box` edits only the header (upper bound, no remap, no tilt). It must not be used
for this step.

## Run it

```bash
# binaries: lego (github.com/biterik/LEGO), dcreator (github.com/biterik/dcreator),
#           lego-cut (this repo); LAMMPS with MANYBODY as Python module 'lammps' or $LMP
export LEGO=/path/to/LEGO/lego  LEGO_STRUCT_DIR=/path/to/LEGO/structures
export DCREATOR=/path/to/dcreator/build/dcreator
POT=/path/to/Cu_Ni_Fischer_2018.eam.alloy ./run_example.sh            # all four characters
POT=/path/to/Cu_Ni_Fischer_2018.eam.alloy ./run_example.sh edge 60    # a subset
```

The potential is not shipped. Use the Fischer, Schmitz & Eich (2019) Cu–Ni EAM from the NIST
Interatomic Potentials Repository (entry 2019--Fischer-F-Schmitz-G-Eich-S-M--Cu-Ni); its Ni has
a0 = 3.5200 Å at 0 K, so `A0=3.52` is the default. Any eam/alloy file with Ni works: set `A0`
to its 0 K lattice constant.

Cells are ~150 × 30 × 146 Å (51 000–59 000 atoms). FIRE minimisation to fnorm < 1e-4 eV/Å
takes ~3 min per relaxation on one core with the LAMMPS Python wheel (2025-07-22), so 8
relaxations take ~25 min.

| step | script | output (`out/`) |
|---|---|---|
| 1. perfect slab | `lego` | `cells/<name>_perfect-slab.lmp`, `cells/lego_<name>.param` |
| 2. insert the dislocation | `dcreator`, `upper_half`, `point` exactly mid-gap between two (111) planes | `cells/<name>_dcreator-inserted.lmp`, `cells/dcreator_<name>.dissparam` |
| 3. delete the pushed-out slab | `lego-cut 0 L0 …` | `cells/<name>_cut.lmp` |
| 4. gate: edge content, Lx′, δ, seam | `scripts/gate_single_dislocation_pps_cell.py` | `cells/gate-report_<name>.txt`, `cells/box-correction_<name>.lmp` |
| 5. impose (Lx′, δ), line along LAMMPS x | `scripts/apply_period_vector.py` | `cells/<name>_corrected-period-vector.lmp` |
| 6. independent seam check | `scripts/check_seam_coordination.py` | (exit status) |
| 7. relax, as built and corrected | `scripts/relax_single_dislocation_pps.lmp` via `scripts/run_lammps_relax.py` | `relaxed/log_relax_<name>_<variant>.lammps`, `relaxed/relaxed_<name>_<variant>_stress-tensor.dump` |
| 8. stress tensor maps and far-field profiles | `scripts/plot_stress_tensor_maps.py` | `figures/` |

Steps 1–6 are done by `scripts/build_single_dislocation_cell.py`; `run_example.sh` drives
everything. Atom data files and dumps are git-ignored; parameter files, gate reports, box
corrections, logs, figures and the summary table are kept as the reference output.

## Results (Fischer 2019 EAM, Ni, a0 = 3.52 Å)

Far-field mean stress (GPa), averaged over |x − x_line| > Lx/4, 15 Å or more from the glide
plane and outside the 2D layers. "top/bot" = upper/lower half-crystal, "mean" = their average.
Positive = tension. Full table: `out/figures/stress-tensor-farfield-summary.dat`.

| character | box | σ_xx top / bot / **mean** | σ_xy top / bot / **mean** | energy vs corrected | line moved during relaxation |
|---|---|---|---|---|---|
| screw | as built | −0.06 / −0.02 / **−0.04** | −1.17 / −0.12 / **−0.64** | +11.4 eV | +14.9 Å |
| screw | corrected | −0.02 / −0.03 / **−0.03** | −0.49 / +0.49 / **+0.00** | 0 | +2.5 Å |
| 30° | as built | +2.30 / +0.13 / **+1.22** | +1.01 / +0.14 / **+0.57** | +18.7 eV | +12.5 Å |
| 30° | corrected | +1.07 / −1.10 / **−0.02** | +0.41 / −0.41 / **+0.00** | 0 | +2.1 Å |
| 60° | as built | +4.03 / +0.28 / **+2.16** | −0.53 / −0.09 / **−0.31** | +38.1 eV | −5.8 Å |
| 60° | corrected | +1.89 / −1.97 / **−0.04** | −0.23 / +0.27 / **+0.02** | 0 | +0.7 Å |
| edge | as built | +4.49 / +0.34 / **+2.41** | 0.00 / −0.01 / **0.00** | +38.9 eV | +1.2 Å |
| edge | corrected | +2.12 / −2.26 / **−0.07** | 0.00 / −0.01 / **0.00** | 0 | +1.2 Å |

Line positions are the centre of the stacking-fault ribbon (hcp atoms within ±3 Å of the glide
plane) relative to where dcreator put the line.

What the example shows:
- With the box as built, the **mean** far-field stress is the error. Missing Lx′ leaves up to
  +2.4 GPa of uniform tension (σ_yy, σ_zz follow through Poisson coupling). Missing δ leaves a
  mean σ_xy ≈ μ·b_s/(2Lx), which on the screw component acts normal to the glide plane, i.e. it
  drives cross-slip.
- The as-built screw and 30° lines **glided 12–15 Å** during the 0 K relaxation. The corrected
  cells stay put within the ±2.5 Å rounding of the construction.
- After the correction every mean is below 0.07 GPa, and the halves are antisymmetric. The
  correction lowers the energy by 11–39 eV per cell (0.4–1.5 eV/Å of line).
- What remains is physical and scales as 1/Lx: ±2.1 GPa σ_xx in the edge cell at Lx = 148 Å,
  ±1.25 GPa at Lx = 250 Å. It biases solute segregation between the two halves. The
  antisymmetric σ_zz (±0.35 GPa) is the bending moment that the z-fixed surface layers take up.

Figures (`out/figures/`):
- `stress-tensor-maps_<name>.png`: all six components as y-averaged x–z maps (one (111) plane
  per z bin, 2.5 Å in x). Top row: box as built; bottom row: corrected. Each component has a
  symmetric colour scale shared by both rows, set by the interior away from the core. The black
  contour marks the stacking fault (hcp atoms); dotted lines mark the 2D layers.
- `stress-tensor-farfield-profiles.png`: σ_ij(z) far from the line, as built vs corrected, all
  characters.

## Orientations and Burgers vectors used

| character | x (glide) | y (line) | z (normal) | b (a0/2 units) |
|---|---|---|---|---|
| screw | [-211] | [0-11] | [111] | [0-11] |
| 60° | [-211] | [0-11] | [111] | [-110] |
| 30° | [-101] | [1-21] | [111] | [0-11] |
| edge | [-101] | [1-21] | [111] | [-101] |
