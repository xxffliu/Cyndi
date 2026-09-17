# Cyndi v4.0 — correctness, performance, and a build that works

Cyndi generates diverse low-energy conformers with a multi-objective
evolutionary algorithm (Liu et al., *BMC Bioinformatics* **10**, 101, 2009).
This release is the result of an external code review on 2026-09-16 and the
work that followed it.

Everything below was measured on the full 329-molecule test set, MMFF94, fixed
seed 0.42, one core. The metric is the minimum Kabsch-aligned heavy-atom RMSD
between each input (bioactive) conformer and the generated set — lower is
better.

| build / settings | wall time | mean min-RMSD | median | ≤1.0 Å | ≤2.0 Å | confs/mol |
|---|---|---|---|---|---|---|
| v2 (previous release), pop 100 / gen 300 / opt-iter 50 | 451 s | 0.633 Å | 0.437 Å | 80.5% | 95.7% | 56.8 |
| + algorithm and crash fixes, same settings | 340 s | 0.639 Å | 0.431 Å | 79.3% | 96.4% | 65.8 |
| + optimiser and force-field work | 138 s | 0.596 Å | 0.365 Å | 79.6% | 97.0% | 68.6 |
| **v4.0 fast default** (mutation 0.4, gen 300) | **131 s** | **0.586 Å** | 0.351 Å | 81.5% | 97.0% | 70.6 |
| **v4.0 accurate** (mutation 0.4, gen 1000) | **344 s** | **0.508 Å** | **0.260 Å** | 81.8% | **98.5%** | 78.8 |
| gen 3000 (diminishing returns) | 723 s | 0.479 Å | 0.240 Å | 85.4% | 97.9% | 79.5 |

**3.4× faster at better accuracy**, or **24% lower mean min-RMSD in 25% less
wall time**, against the previous release.

One thing to be clear about: the bug fixes on their own bought correctness and
a 1.3× speed-up but **no measurable accuracy change** (0.633 → 0.639 Å is
inside the noise of a single seed). The accuracy came from spending the freed
time on a larger search budget.

## Fixed

**Polynomial mutation was dead code.** `real_mutate` computed `delta` with
integer division, so `delta` was always 0, `deltaq` was always 0, and the
operator never changed a gene: over a 300-generation run of `3ert` it visited
5400 genes and altered **2** of them. The entire search had been driven by SBX
crossover alone. Every previously published conclusion about
`MOGA_Mutation_Probability` was therefore measuring nothing but a change in
random-number consumption.

**Circular crossover undid itself.** The shortest-arc mapping ended with
`par2 = wrap_torsion(par1 + d)`, which folds the result straight back into
`[-72, 71]`. For 24.7% of parent pairs that restored the original value, and
SBX then spanned the long way round the circle.

**Segfault on inputs over ~580 molecules.** The per-molecule RNG seed
(`BasicSeed_ + 0.001 × index`) left the `(0,1)` interval the Knuth subtractive
generator requires, after which `randomperc()` returned values of order 1e10
and `(int)(rnds * PopSize_)` became an arbitrary array index.

**Chunked and whole-file runs disagreed** (the long-standing
`KNOWN_ISSUES.md` #1). `GlobalIdx` was incremented at the end of the loop body,
so the four `continue` paths skipped it and every molecule after the first
skipped one got a different seed in a chunked run. A second defect in
`MOL2IO::write` — stream state (`adjustfield`, `floatfield`, `precision`,
`showpoint`) inherited from the previous molecule's atom section — made the
first molecule of each chunk format differently. Both fixed and verified:
6339 conformer records byte-identical between one process and three chunks.

**The RMSD evaluator itself was wrong.** `bench/analyze.py` obtained the
largest eigenvalue of the quaternion matrix F by power iteration, which
converges to the eigenvalue of largest *magnitude*. F is symmetric and
traceless, so it always has negative eigenvalues, and whenever the most
negative one dominated the reported RMSD was far too large: **1084 of 3344
conformer pairs (32.4%) wrong by more than 1e-3 Å, worst error 14.7 Å**, and
39 of 327 molecules got the wrong minimum (1cle reported 16.4 Å, actually
2.73 Å). Any number published from an earlier version of that script needs
re-measuring.

Also: duplicate conformers accepted into the diversity archive (622
"conformers" of which 68 were distinct, for a one-rotatable-bond ligand);
`MOGA_SBX 15 20` parsing only the first value; `MOGA_Max_Opt_Gradient` parsed
and then never handed to the minimiser; `Conformer::rmsd` filled from
`-Rg`; `ForceField`/`MMFF94`/`TAFF` copy constructors that shallow-copied
component pointers every copy then deleted; `ForceField::add_unassigned_atom`
with an inverted membership test; RNG state declared `static` in a header;
and the output table's column widths, which glued the conformer count to the
timing.

## Changed

- **Retuned defaults**: `MOGA_Mutation_Probability` 0.4 (not the textbook
  `1/N` — with η_m = 20 the periodic polynomial operator is a ~16° local
  perturbation, so a high rate behaves like a local search layered on SBX),
  `MOGA_Max_Opt_Iteration` 20, `MOGA_Max_Opt_Gradient` 0.5.
- The CG line search is capped at 10 interpolation steps instead of 50. On the
  329 set that is 178 s → 138 s with mean/median min-RMSD unchanged.
- One square root instead of two in the MMFF94 non-bonded force kernels;
  `MMFF94.parm` parsed once per run instead of two to three times per molecule.

## Added

- **A CMake build** (`CMakeLists.txt`, `CMakePresets.json`) replacing the
  hand-written `build.sh` and the VS2008-era `Cyndi.sln`. On Windows,
  `build_windows.bat` installs `cmake` and `ninja` with winget if needed and
  builds through the MinGW preset. Verified with GCC 16.2.0, CMake 4.4.3 and
  Ninja 1.13.2.
- **Five CTest regression tests** for the invariants this release established —
  not a unit-test suite, but the things that have actually broken:
  chunk-vs-whole reproducibility, the >580-molecule seed crash, the RMSD
  evaluator (including the near-mirror-image case that defeated power
  iteration), a single-molecule smoke run, and CLI argument rejection. The
  first two fail against a pre-v4 binary, which is the point of having them.
- `bench/rmsd_core.py`, shared by `bench/analyze.py` and `bench/bench_rmsd.py`,
  using `numpy.linalg.eigvalsh` with a numpy-free fallback that power-iterates
  on `F + cI`.
- `bench/RETUNE_2026-09-16.md` with the full parameter sweep.

## Removed

The pre-CMake build paths (`build.sh`, `Cyndi.sln`, the `.vcproj`/`.vcxproj`
files), `Cyndi_batch.py` (Python 2, did not run on any current interpreter —
use `bench/run_parallel.py`), the benchmark result files produced by the broken
evaluator, the superseded 20-molecule tuning sweep, three sets of duplicated
data and parameter files, and the 2008-era sample output that shipped in the
source tree.

## Still open

- **CG minimisation aborts on a large fraction of conformers** — 1520 over the
  329 set, and on some inputs essentially every conformer. The message is now a
  counter rather than a per-conformer print, so it is visible; the underlying
  `ConicLineSearch` port from BALL has not been diagnosed. See
  `KNOWN_ISSUES.md` #5.
- **`MOL` copy semantics are incomplete** (`KNOWN_ISSUES.md` #6). Harmless
  today because nothing downstream calls `get_bond_list()` or `apply_rotor()`
  on a copy; a trap for anything new.
- ~750 `-Wall -Wextra` warnings, including 167 out-of-order member
  initialisers. `CYNDI_WERROR` and `CYNDI_WERROR_REORDER` exist and default to
  OFF for exactly that reason.
- The electrostatic term uses a buffered Coulomb potential with dielectric 1,
  which in vacuum over-rewards intramolecular salt bridges. A
  distance-dependent dielectric is worth an A/B on flexible charged ligands.
