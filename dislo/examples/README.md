# lego-dislo examples

Three small, runnable examples (a couple of thousand atoms each), in the
spirit of the original 2005 `example_aniso_edge_howto.txt` workflow:
build a crystal, run the solver, add the displacement field, then
post-process with the other LEGO tools.

Build the tools first (`make` at the repository root). Then:

```sh
./run_example.sh edge      # W 1/2<111>{110} edge dislocation
./run_example.sh screw     # W 1/2<111> screw dislocation
./run_example.sh dipole    # two opposite edge dislocations (superposition)
```

Each run:

1. generates a small oriented W (bcc, a0 = 3.1652 A) crystal with
   `make_crystal.py`, centered on the origin, so the dislocation core
   (`point:` in the card, chosen slightly off-lattice) sits mid-crystal;
2. runs `lego-dislo` with the corresponding card (`edge.yaml`,
   `screw.yaml`, `dipole.yaml`), which prints the sextic solution
   (roots P, coefficient matrices F and G) and the energy prefactor
   K·b²/4π, and writes them to `out/<case>_report.txt`;
3. adds the anisotropic elastic displacement field to every atom and
   writes `out/<case>_out.data` (LAMMPS data format);
4. compares against the committed reference in `expected/` (tolerance
   1e-6 A — the files should agree to ~1e-12, the tolerance only absorbs
   libm differences between platforms).

For the edge case the reported factors reproduce the legacy Disloelast
`W_edge_results` reference: K·b²/4π = 1.3276411e+10 N/m² (in
lattice-parameter units of b), i.e. 1.3300961e-09 J/m, K = 222.3 GPa.

Typical follow-up steps, exactly like the 2005 workflow (the box bounds
are left unchanged by lego-dislo):

```sh
../../lego-tools/lego-cut  -- -80 80 -80 80 -4 4  out/edge_out.data cut.data
../../lego-tools/lego-shift 5 0 0 out/edge_out.data shifted.data
```

`make_crystal.py --help` documents the generator (any orthogonal
right-handed orientation, any box size).
