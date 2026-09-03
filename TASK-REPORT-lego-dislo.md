# TASK REPORT — lego-dislo (Disloelast modernization)

Branch: `feature/lego-dislo` (6 commits on top of `main`, not pushed, not
merged). Date: 2026-09-03. All validations pass on this machine (M5
MacBook Pro, macOS, Apple clang + Homebrew libomp, 18 cores, 64 GB).

## What was built

`dislo/` — a new LEGO-TOOLS component:

| File | Content |
|---|---|
| `sextic.c/.h` | Faithful C port of `elast-cubic.f`: 6×6→9×9 stiffness handling, tensor rotation into the dislocation frame, sextic determinant, root post-processing (legacy 1e-4 nudge + positional selection), eigenvectors, 6×6 solve (dgesv-equivalent), F/G coefficient matrices, K·b²/4π; plus the field evaluation ported from `anisodisloc.f` (branch-cut logic, PHICUT, in-place nudges) |
| `findroots.c` | Faithful port of `FindRoots.f` (Laguerre with deflation + polishing; order-preserving) |
| `card.c/.h` | Hand-written parser for the flat YAML-subset card (cubic shorthand or full 6×6 Voigt, GPa) |
| `main.c` | `lego-dislo` CLI on atomio (LAMMPS data/dump, IMD, gzip), OpenMP atom loop, frame composition box→crystal→dislocation, superposition, report writer, `--solve-only`, field timing |
| `python/legodislo/` | Stdlib-only driver: `Material`/`Dislocation` dataclasses, `insert_dislocations()`/`solve()`, report parsing into typed results, `LegoDisloError`, binary auto-location, `python -m legodislo` CLI twin |
| `examples/` | `make_crystal.py` (oriented bcc generator), edge/screw/dipole cards, one-command `run_example.sh` with committed `expected/` outputs (~130 KB each) |
| `tests/` | `unit_tests.py` (self-contained), `golden_validate.py` (auto-skips without `~/DEVEL/Disloelast`; `DISLOELAST_DIR` overrides), `run_tests.sh`, `bench.sh` + `gen_bcc.c` |
| docs | `dislo/README.md` (theory, frames, card grammar, API, limitations), `dislo/TOOLCARD.md` (1-page agent reference), top-level README section + layout, `tests/README.md` update, top-level Makefile + .gitignore wiring |

## Validation numbers

Golden references read from `/Users/e.bitzek/DEVEL/Disloelast` (harness
auto-skips with a notice if absent).

| # | Case | Result |
|---|---|---|
| 1 | Solver, W edge (X=[111], Y=[-101], Z=[1-21], b∥[-1-1-1], \|b\|=0.8663) | every printed digit of `W_edge_results` reproduced; **max rel dev 3.7e-7** on P/F/G; K·b²/4π = 1.3276411e10 N/m² vs golden 0.1327641e+11 (**rel dev < 1e-7**); in SI 1.3300961e-09 J/m, K = 222.3 GPa |
| 2 | Solver, W screw (X=[-1-12], Y=[1-10], Z=[111], b∥[111]) | `W_screw_results` reproduced; **max rel dev 2.0e-7** on P/F/G (roots matched by nearest-root pairing, see decision 2) |
| 3 | Field, edge — `example/W_45x_45y_1.5z.s.bcc` (197 200 atoms) vs `anisod_…bz2` | **max \|Δr\| = 1.97e-6 Å, RMS = 1.46e-6 Å** (tolerance 1e-3 Å) |
| 4 | Field, screw — `W_x-12-1_y-101_z111.sml.bcc` (330 000 atoms) vs `anisoscrew_…` | **max \|Δr\| = 2.73e-7 Å, RMS = 1.73e-7 Å** |
| 5 | Superposition | zero-b second dislocation → output **byte-identical** to single run; dipole (opposite b, 10 Å apart): max \|u\| 0.284 Å at r≈20 Å → 0.103 Å at r>40 Å (decays), vs 15.5 Å for the single dislocation at r>40 Å |
| 6 | Rotation equivalence (same physical dislocation through two different `box_orientation` frames, 492 points) | **max deviation 1.0e-12 Å** (tolerance 1e-8) |
| 7 | Unit tests | cubic vs full 6×6 → byte-identical reports; 9 rejection cases (non-orthogonal / left-handed box, glide normal ∦⊥ line, missing keys, bad lists, tabs, unknown keys, asymmetric C) all rejected with clear messages; comments accepted |

The residual ~1e-6 Å in case 3 is fully explained by two legacy
single-precision artifacts the double-precision port deliberately does
not reproduce: `anisodisloc.f` re-read P/F/G from the E14.7-formatted
results file (≈8 significant digits), and used pi = 3.141592654 assigned
through a single-precision literal (rel. error 2.7e-8 in the arctan
terms). No atom sits near enough to the branch cut for these to amplify.

**gfortran cross-check** (optional item): `elast-cubic.f` + bundled
LAPACK recompiled with Homebrew gfortran (`-std=legacy -w`, Accelerate
for BLAS) on this machine reproduces both golden results files **exactly,
digit for digit, including the root order** — confirming the golden files
and, transitively, the port.

## Benchmark (`dislo/tests/bench.sh`)

M5 MacBook Pro, 18 cores, 64 GB. "Field time" is the displacement loop
only (reported by the binary); "total" includes reading/writing the
LAMMPS data files.

| Atoms | Threads | Field time | Total wall |
|---|---|---|---|
| 10 000 422 | 1 | 0.454 s | 7.6 s |
| 10 000 422 | 4 | 0.121 s | 8.2 s |
| 10 000 422 | 18 | 0.041 s | 7.5 s |
| 99 672 064 | 1 | 4.483 s | 73.1 s |
| 99 672 064 | 4 | 1.266 s | 81.2 s |
| 99 672 064 | 18 | **0.372 s** | 75.2 s |

Field evaluation scales ~linearly (12× at 18 threads) and meets the
target (10⁸ atoms in seconds). Total wall time at 10⁸ atoms is dominated
by ~11 GB of text I/O in atomio (~5.4 GB in, ~6 GB out) and is
essentially thread-independent. Memory stays within atomio's in-memory
layout (~5 GB atom data + I/O buffers); no additional per-atom arrays are
allocated (displacements accumulate in registers per atom).

## Decisions taken (beyond those fixed in the task)

1. **Ported `elast-cubic.f`, not `elastbart.f`.** The diff shows
   `elastbart.f` is the older implicit-typing variant (default REAL,
   COMPLEX*8); `elast-cubic.f` is the explicit double-precision version
   the reference runs used (its output matches the goldens). The general
   6×6 capability is preserved by building the 9×9 matrix directly from
   the full Voigt matrix (`cij: full:`), which reduces exactly to the
   legacy cubic/tetragonal construction.
2. **Sextic root order is matched, not forced.** The deflation order of
   Laguerre's method starting from x=0 is sensitive at the last bit; the
   2005 x87 binary's order is not reproducible in principle (80-bit
   intermediates). The port keeps the exact legacy algorithm (order and
   all); the golden comparison pairs roots by nearest-neighbour matching
   and permutes F/G columns accordingly. On this machine the harness
   cards happen to give the identity order for both cases. The summed
   displacement field is invariant under the permutation — proven
   directly by validations 3–6. Consequence: the *report* columns may be
   permuted on other platforms; documented in the harness.
3. **Legacy regularization kept verbatim**: after root finding, real or
   imaginary parts with magnitude < 1e-7 are set to +1e-4 (the
   `0.1000000E-03` visible in the legacy P output) *before* everything
   downstream — this is the deliberate degeneracy-breaking the task
   flagged. Same for the field-side nudges (`x += 1e-7` when
   x + PR·y == 0, persisting across the term loops; ETAR += 1e-7 in
   PHICUT) and the [0,2π) + PHICUT branch correction, ported
   conditional-for-conditional. The dead "core cleaning" loop was not
   resurrected (single field evaluation per atom, as effective in 2005).
4. **Two documented micro-deviations in dead code paths** (commented in
   `sextic.c`): Fortran loop-37 branch 3 tests `AR13²+AI23²` (a typo
   mixing two subdeterminants) — the port uses the intended |D13|²; and
   branch 2 never assigns AI(2,N) (stale-value bug) — the port computes
   the consistent imaginary part. Neither branch is reachable for any
   validated input (all cases take branch 1); a faithful reproduction of
   the stale-memory behavior would be unreproducible anyway.
5. **Screw golden frame resolved** (task validation 4): the legacy run
   (`myrun_elastbart_Wscrew.bsh` + `input.w_111screw`) applied the solver
   output of frame X=[-1,-1,2], Y=[1,-1,0], Z=[1,1,1] directly in the
   sample's coordinates, whose axes are x=[-1,2,-1], y=[-1,0,1],
   z=[1,1,1]. A 120° rotation about [111] is a cubic symmetry operation,
   so both frames have identical rotated stiffness and identical P/F/G;
   the equivalent card therefore uses the sample's own axes as the
   dislocation frame (`glide_normal: [-1,0,1]`). Validated at 2.7e-7 Å;
   frame composition itself is independently covered by validation 6.
6. **π is double precision** (M_PI) throughout the new code — a "clean
   double-precision port" per design decision 4; see the error budget
   above.
7. **`--in-format` is accepted but advisory**: atomio's reader is
   auto-detecting by design; on a mismatch lego-dislo warns and proceeds.
   (Kept for CLI-spec compatibility rather than adding a forced-format
   path to atomio.)
8. **Legacy raw-IMD quirk kept out of atomio** (per spec): the golden
   files (plain `id type mass x y z` lines with the IMD header appended
   at the END) are read by a tolerant reader inside the test harness,
   which feeds proper IMD to the binary.
9. **Committed expected example outputs** (~130 KB each, well under the
   1 MB limit) instead of checksums — float text is not bit-portable
   across libms, so the runner compares numerically (1e-6 Å).
10. **`make no-omp` in `dislo/` forces an OpenMP-free `atomio.o`
    rebuild** — linking a stale OpenMP atomio.o without `-lomp` fails.
    (The same latent staleness exists in `lego-tools/no-omp` itself; not
    touched, as it is out of scope.)
11. **Report format**: plain `key: value` lines (greppable, parsed by the
    Python driver); K given three ways — GPa·(b in a0)², J/m, and the
    energy factor K in GPa (K = 4π·[K·b²/4π]/|b|², sanity-checked against
    the isotropic estimate μ/(1−ν) ≈ 223 GPa for W).

## How to run everything

```sh
make                       # builds lego-tools, afc, dislo/lego-dislo
make test                  # smoke test incl. a lego-dislo section
./dislo/tests/run_tests.sh # unit tests + golden validation (auto-skip)
BENCH_ATOMS=1e8 ./dislo/tests/bench.sh
cd dislo/examples && ./run_example.sh edge   # or screw | dipole
```

## Open questions

1. **Report column order across platforms** (decision 2): if bit-stable
   report output ever matters (e.g. for regression diffs on CI across
   architectures), a canonical root ordering could be *added* — but it
   would then deliberately deviate from the legacy positional selection.
   Left as-is pending Erik's preference.
2. **End-to-end wall time at 10⁸ atoms** is ~75 s, I/O-bound in atomio's
   text parsing/writing (the field itself takes 0.4 s). If "well under a
   minute" must hold end-to-end rather than for the field evaluation,
   atomio's writer would need parallelization or a binary format — out of
   scope here.
3. **Legacy type-0 atoms**: the Disloelast samples use IMD type 0, which
   LAMMPS data cannot represent (types are 1-based). The IMD path handles
   them fine; converting such files to LAMMPS data output would need a
   type remap (existing atomio behavior, unchanged).
4. `elast-tetra-ag.f` (tetragonal variant with a different input layout)
   and `strain.f` (stress evaluation) were not ported — the stress-field
   machinery (BSR/BSI arrays, the SGM loop) exists in `elast-cubic.f` but
   its output was discarded there too. The full 6×6 card input covers the
   tetragonal case functionally.
