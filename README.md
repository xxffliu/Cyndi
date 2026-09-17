# Cyndi

**Cyndi: a multi-objective evolution algorithm based method for bioactive molecular conformational generation**

Xiaofeng Liu, Fang Bai, Sisheng Ouyang, Xicheng Wang, Honglin Li, Hualiang Jiang

*BMC Bioinformatics* **10**, 101 (2009). DOI: [10.1186/1471-2105-10-101](https://link.springer.com/article/10.1186/1471-2105-10-101)

## Overview

Cyndi generates diverse, low-energy molecular conformers using a **multi-objective evolutionary algorithm (MOEA)**. Each rotatable dihedral angle is encoded as a gene; the population is evolved under competing objectives:

| Objective | Meaning |
|---|---|
| VDW energy | sterically favourable |
| Torsion energy | dihedral preferences |
| RMSD vs. input | conformational diversity |
| −gyration radius | molecular extension / compactness |

The Pareto archive (ε-domination based) yields conformers that are both energy-favoured and evenly scattered across conformational space. Reported performance: average minimum RMSD of 0.864 Å to bioactive conformations over 329 test molecules, ~0.49 s per molecule.

## Building

CMake is the supported build (same layout as SHAFTS: a `cyndi_core` static
library holding everything except `main()`, one executable, and the end-to-end
invariants registered with CTest).

### Any platform

```bash
cmake --preset default        # or: cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release
cmake --build --preset default -j
ctest --preset default
```

Requires CMake 3.15+ and a C++11 compiler. The build directory is
`build-cmake` (`build-mingw` for the Windows preset); `build/` is left free
because it was the scratch directory for the 2026 benchmark runs and is
gitignored. The five registered tests take about seven seconds; see "Tests"
below for what they cover.

### Windows (MinGW-w64)

```bat
build_windows.bat
```

That is the whole thing: it installs `cmake` and `ninja` with winget if they are
not already on PATH, then configures, builds and tests through the
`mingw-release` preset. `build_windows.bat debug` and `build_windows.bat linear`
select the other two presets. To drive the presets yourself:

```bat
cmake --preset mingw-release
cmake --build --preset mingw-release
ctest --preset mingw-release
```

`build_windows.bat` finds the toolchain from `%CYNDI_MINGW%`, then from `g++`
on PATH, then from `%USERPROFILE%\tools\mingw64`, `C:\msys64\mingw64` and
`C:\mingw64`. The preset itself pins nothing -- it takes `g++` from PATH. To
pin a specific compiler without touching a tracked file, put it in
`CMakeUserPresets.json` (gitignored):

```json
{ "version": 3, "configurePresets": [
  { "name": "mine", "inherits": "mingw-release",
    "cacheVariables": { "CMAKE_CXX_COMPILER": "C:/msys64/mingw64/bin/g++.exe" } } ] }
```

Two Windows-specific things the build handles that the old hand-written
`build.sh` used to handle by hand (it and the VS2008-era `Cyndi.sln` were
removed in the 2026-09-17 cleanup; CMake replaces both):

- **Sysroot headers.** This MinGW-w64 build does not put its own
  `<prefix>/include` on the default search path. CMake derives it from the
  compiler location and adds it with `-idirafter`, so it is searched last and
  never shadows `include/`.
- **Runtime DLLs.** `CYNDI_STATIC_RUNTIME` is ON by default on MinGW, so the
  executable is statically linked and needs no `libstdc++-6.dll` /
  `libgcc_s_seh-1.dll` / `libwinpthread-1.dll` beside it. The old script
  copied those three next to the binary instead.

`MMFF94.parm` and `TAFF.parm` are opened by bare filename, i.e. relative to the
working directory. The build copies both next to the executable so it runs from
the build directory as-is.

### Build options

| option | default | what it does |
|---|---|---|
| `CYNDI_BUILD_TESTS` | ON | register the CTest suite |
| `CYNDI_TORSION_LINEAR` | OFF | pre-v1 linear torsion encoding (clamp at ±180° instead of wrapping) — for A/B against the original behaviour |
| `CYNDI_STATIC_RUNTIME` | ON on MinGW | link libgcc/libstdc++ statically |
| `CYNDI_SANITIZE` | *(empty)* | `address`, `thread` or `undefined` |
| `CYNDI_WERROR` | OFF | `-Werror`. This tree still emits ~750 `-Wall -Wextra` warnings, so turning it on fails the build today. |
| `CYNDI_WERROR_REORDER` | OFF | `-Werror=reorder`. 167 out-of-order member initialisers remain, 108 of them in the single `MOGAParam` constructor in `include/ParamInput.h`. All harmless today (every initialiser is a literal), but that is luck rather than design — worth cleaning up so this can be turned on. |

### Tests

`ctest` registers the end-to-end invariants the 2026-09 review established, not
a unit-test suite (there isn't one):

| test | what it protects |
|---|---|
| `runs_single_molecule` | the program completes on a small ligand |
| `bad_option_rejected` | an unknown CLI option exits non-zero |
| `chunk_reproducibility` | a `-startidx`/`-maxmols` run reproduces the same molecule from a whole-file run byte for byte — `KNOWN_ISSUES.md` #1, which has broken twice for two different reasons |
| `seed_wraparound` | molecules past global index 580 do not segfault — `KNOWN_ISSUES.md` #4 |
| `rmsd_evaluator` | `bench/rmsd_core.py` returns the true minimum RMSD, including on the near-mirror-image pairs where the old power-iteration evaluator returned the wrong eigenvalue |

The last three are Python and skip themselves (CTest "Skipped") when Python 3 or
their input data is missing. The first two of those three fail against a
pre-2026-09-16 binary, which is the point of having them.


## Usage

```
Cyndi -input input.mol2 -output output.mol2 [-parm CyndiParam.in]
```

- Input: SYBYL MOL2 (single molecule or multi-molecule `@<TRIPOS>MOLECULE` blocks)
- Output: MOL2 with up to `MOGA_Max_Conformers` conformers per molecule
- Parameters: see `Cyndi/CyndiParam.in` (population size, generations, crossover/mutation probabilities, energy cutoff, RMSD clustering threshold, force field choice TAFF/MMFF94, etc.)

Batch mode:

```bash
python bench/run_parallel.py -input <file.mol2> -parm parameter/CyndiParam.in -n 4 -out all.mol2 -seed 0.42
```

### Test data

`bench/329_test_set.mol2` is the 329-molecule test set from the paper (PDB
ligands); `bench/subset20.mol2` is a 20-molecule slice of it for quick checks.
`Cyndi/` holds the four small example ligands (`1hnn`, `3ert`, `1a28`,
`EGFR_ah`) and the force-field parameter files.

## Layout

```
src/                C++ sources (MOGA core, force fields, minimizers, IO)
include/            headers
Cyndi/              force-field parameter files and the small example ligands
parameter/          the default CyndiParam.in
bench/              329-molecule benchmark set, RMSD evaluator, run drivers
tests/              the Python end-to-end tests CTest registers
CMakeLists.txt      the build (see Building)
CMakePresets.json   configure/build/test presets
build_windows.bat   one-shot Windows build through the MinGW preset
```

## v4 Correctness + Performance Pass (2026-09-16)

Follows an external code review (2026-09-16), working through its section 6
rollout order. **The measured numbers in the v1/v2 section below
are not trustworthy** and are kept only for history: the evaluator that
produced them (`bench/analyze.py`, `bench/bench_rmsd.py`) had a bug of its own
(see "Evaluation" below). Everything in this section was re-measured with the
corrected evaluator on the full 329-molecule test set, fixed seed 0.42, one
core, MMFF94.

### Where the code stood, and where it stands

| build / settings | wall time | mean min-RMSD | median | <=1.0 A | <=2.0 A | confs/mol |
|---|---|---|---|---|---|---|
| v2 code as-found, pop 100 / gen 300 / opt-iter 50 | 451 s | 0.633 A | 0.437 A | 80.5% | 95.7% | 56.8 |
| + algorithm & crash fixes (v4 step 1), same settings | 340 s | 0.639 A | 0.431 A | 79.3% | 96.4% | 65.8 |
| + optimiser & force-field work (v4 step 3) | 138 s | 0.596 A | 0.365 A | 79.6% | 97.0% | 68.6 |
| **+ retuned: mutation 0.4** (fast default) | **131 s** | **0.586 A** | 0.351 A | 81.5% | 97.0% | 70.6 |
| **+ gen 1000, mutation 0.4** (accurate default) | **344 s** | **0.508 A** | **0.260 A** | 81.8% | **98.5%** | 78.8 |
| + gen 3000, mutation 0.4 (diminishing returns) | 723 s | 0.479 A | 0.240 A | 85.4% | 97.9% | 79.5 |

So: **3.4x faster at better accuracy**, or **24% lower mean min-RMSD for 25%
less wall time** than the code this pass started from. The bug fixes alone
bought correctness and a 1.3x speed-up but no measurable accuracy change --
the accuracy comes from spending the freed time on more search.

### Retuning (all previous tuning rested on a dead mutation operator)

`real_mutate` computed `delta` with integer division, so `delta` was always 0,
`deltaq` was always 0, and polynomial mutation never changed a gene: over a
300-generation run of `3ert` the operator touched 5400 genes and altered **2**
of them (both from the `y == yl` fallback path). The whole search was driven by
SBX crossover alone. After the fix, 1006 of 1103 drawn genes actually move.
Every conclusion in the v1/v2 tuning sweep about `MOGA_Mutation_Probability`
was therefore measuring nothing but a change of random-number consumption.

Re-swept at pop 100 / gen 300 (mean min-RMSD / median / <=1 A / <=2 A):

| mutation | 0.05 | 0.10 | 0.20 | 0.30 | 0.40 | 0.60 |
|---|---|---|---|---|---|---|
| mean | 0.621 | 0.624 | 0.596 | 0.594 | **0.586** | 0.581 |
| median | 0.373 | 0.397 | 0.365 | 0.361 | 0.351 | 0.343 |
| <=1 A | 77.2% | 79.3% | 79.6% | 79.9% | **81.5%** | 81.5% |
| <=2 A | 95.4% | 96.0% | 97.0% | 96.4% | **97.0%** | 95.7% |

Higher mutation helps -- the opposite of the usual `1/N` guidance, because
with `MOGA_SBX`'s eta_m = 20 the periodic polynomial operator is a ~16 deg
local perturbation, not a restart. 0.4 is the recommended default; past 0.4
the mean keeps creeping down but the <=2 A rate starts to suffer.

Search budget matters more than population size. One "generation" in
`execuateMOGA` produces two children, so total energy evaluations are roughly
`PopSize + 2 * MaxNumGen` -- population only sets the initial sample:

| setting | wall time | mean | median | <=1 A | <=2 A |
|---|---|---|---|---|---|
| pop 100 / gen 300 | 138 s | 0.596 | 0.365 | 79.6% | 97.0% |
| pop 300 / gen 300 | 117 s | 0.609 | 0.385 | 79.6% | 96.0% |
| pop 100 / gen 1000 | 339 s | 0.521 | 0.268 | 83.6% | 97.9% |
| pop 100 / gen 1000, mut 0.4 | 344 s | 0.508 | 0.260 | 81.8% | 98.5% |
| pop 100 / gen 3000, mut 0.4 | 723 s | 0.479 | 0.240 | 85.4% | 97.9% |

Returns keep coming but flatten: tripling the budget again (gen 1000 -> 3000)
costs 2.1x the time for a further 6% on the mean. gen 1000 is the sensible
default; gen 3000 is worth it when a specific flexible ligand matters more
than throughput.

### Algorithm and crash fixes (step 1 -- the review's patch, plus what it missed)

- **Polynomial mutation was dead code** (integer division in `delta`, see
  above). The review's patch fixed the downstream rounding but not `delta`
  itself. On the circular encoding the bounded form has no meaning, so the
  periodic branch now uses Deb's unbounded operator (`delta = 1`, i.e.
  `xy = 0`) and wraps the child; `TORSION_LINEAR` keeps the bounded form with
  the correct `delta1`/`delta2` split.
- **Circular crossover undid itself.** The shortest-arc mapping ended with
  `par2 = wrap_torsion(par1 + d)`, which folds the result straight back into
  `[-72, 71]` -- for 24.7% of parent pairs that restored the original value and
  SBX spanned the long way round the circle. `par2 = par1 + d` is now left
  unwrapped and only the children are wrapped, with the unbounded SBX spread
  (`alpha = 2`).
- **SBX `beta` used integer division** and children were truncated toward zero
  (a systematic bias); both fixed.
- **`GlobalIdx` skipped molecules**, which was the root cause of
  `KNOWN_ISSUES.md` #1. See that file -- now closed and verified byte-identical.
- **Segfault on inputs over ~580 molecules**: the per-molecule seed left the
  `(0,1)` interval the Knuth generator requires. Reproduced and fixed.
- **Duplicate conformers**: the diversity archive accepted individuals that
  were exact duplicates (`dom_check` needs a strict improvement, so equal
  objective vectors are never dominated), and `con_update` handed it children
  the epsilon archive had rejected for being *over the energy cutoff*. For a
  one-rotatable-bond ligand this produced 622 "conformers" of which 68 were
  distinct, all of them minimised. 1hnn 1.82 s -> 0.25 s, 1a28 3.75 s ->
  0.28 s.
- **`Conformer` field mapping** still used the original 4-objective layout, so
  `rmsd` was filled from `-Rg`. Only visible with
  `MOGA_Optimize_Conformer N`.
- `MOGA_SBX 15 20` parsed only the first value, so `n_distribution_m` was 15.
- Archive-update index errors (`i--` after a deletion; restart from slot 0),
  and the output table's column widths, which glued the conformer count and
  the timing together ("100.144822" = 10 conformers + 0.144822 s).

### Optimiser and force field (step 3)

- **`MOGA_Max_Opt_Gradient` was parsed and then never used.** The minimiser ran
  against the built-in `MAX_GRADIENT` of 0.01 kJ/(mol A), which it essentially
  never reaches, so every conformer burned the full iteration budget. It is now
  passed through `cgm.setup(ff, option)`; the shipped parameter files use 0.5
  with `MOGA_Max_Opt_Iteration 20`. It also no longer depends on appearing
  after `MOGA_Optimize_Conformer` in the file.
- **Line search capped at 10 interpolation steps** (was 50). Each step costs a
  full energy *and* gradient evaluation. On the 329 set this is 178 s -> 138 s
  with mean/median min-RMSD unchanged at 0.596 / 0.365 A -- the extra 40
  interpolations never bought anything.
- **One square root instead of two** in `MMFF94VDW::update_forces` and
  `MMFF94Ele::update_forces` (`length()` followed by `normalize()`, which calls
  `length()` again), on the hottest loop in the program. ~9%.
- **`MMFF94.parm` is parsed once**, not on every force-field setup (260 KB,
  two to three times per molecule).
- Copying a `ForceField`, `MMFF94` or `TAFF` shallow-copied the component
  pointers that every copy's destructor then deleted. Now `= delete`.
- `ForceField::add_unassigned_atom` had its membership test inverted, so the
  list of atoms with unassigned force-field types was always empty.
- `individual::gyration_radius` was never initialised; `individual::operator=`
  dropped `cv`; `MOGA::clear()` reset `taff_` but not `mmff94_` (the default
  force field) or the diversity archive.
- The RNG state (`oldrand[55]` et al.) was `static` in `random.h`, i.e. one
  private copy per translation unit. Now `extern` + a single definition.
- `MOL::updateRMSD()` summed over heavy atoms, divided by all atoms, and never
  took a square root. (Nothing calls it; fixed rather than left as a trap.)
- The "No convergence and the step computation failed" message is now a counter
  reported once. It fires 1520 times over the 329 set and on essentially every
  conformer of some inputs -- see `KNOWN_ISSUES.md` #5.

### Evaluation

`bench/analyze.py` and `bench/bench_rmsd.py` found the largest eigenvalue of
the quaternion matrix F by power iteration, which converges to the eigenvalue
of largest *magnitude*. F is symmetric and traceless, so it always has negative
eigenvalues, and whenever the most negative one dominated the RMSD came out
far too large. Measured over one 329-molecule run: **1084 of 3344 conformer
pairs (32.4%) wrong by more than 1e-3 A, worst error 14.7 A**, and 39 of 327
molecules got the wrong *minimum* (1cle 16.4 A reported vs 2.73 A actual,
1hvl 11.8 vs 4.09, 4est 10.4 vs 2.94). Mean min-RMSD over the set: 1.175 A
reported vs 0.995 A actual.

Both scripts now share `bench/rmsd_core.py`, which uses
`numpy.linalg.eigvalsh` (with a numpy-free fallback that power-iterates on
`F + cI`, where the shift makes largest-magnitude and largest-algebraic
coincide). Heavy atoms are selected from the SYBYL type column rather than the
atom-name column. The C++ side (`MOL::minimizeRMSD`) uses a full Jacobi
decomposition and was never affected.

**Any number in the v1/v2 section below that came from the old evaluator needs
re-measuring before it is quoted**, including the per-version "improved on N /
regressed on M" counts.

### Batch reproducibility

`bench/run_parallel.py`'s bit-identity claim now actually holds. Besides the
`GlobalIdx` fix, `MOL2IO::write` was leaving stream state (`adjustfield`,
`floatfield`, `precision`, `showpoint`) set by the atom section, which the
*next* molecule's header inherited -- so the first molecule of every chunk
formatted differently from the same molecule in a whole-file run. Verified on
a 90-molecule subset: 6339 conformer records byte-identical between one
process and three chunks (the `# Created ... on <date>` line aside).

## Notes (2026 maintenance)

- Ported to build with modern GCC 16 (MinGW) — original code is VS2008-era
- Fixed: `dynamic_cast == false` (MMFF94StretchBend), `name == "****"` comparison-as-assignment (MOL2IO), CG minimizer never initialized in the post-MOGA optimization step
- Optimized: force-field setup hoisted out of the per-conformer loop (~3x speedup), `pow()` → multiplication in VDW kernels
- Clustering now uses Kabsch-aligned heavy-atom RMSD (was unaligned RMSD)

## v1 Algorithm Rework (2026)

> **Historical.** The min-RMSD figures in this section were produced with the
> broken RMSD evaluator described under "Evaluation" above and should not be
> quoted. The v4 section has re-measured numbers.


Rebuilt the MOGA objective design and torsion encoding. Validated on the full
**329-molecule test set** (min-RMSD of generated conformers to the input
/bioactive conformer, Kabsch-aligned heavy atoms, lower is better):

| Configuration | mean min-RMSD | median | ≤1.0 Å | ≤2.0 Å | mean confs | total time |
|---|---|---|---|---|---|---|
| original 4-obj (VDW, torsion, RMSD, -Rg) | 1.044 Å | 0.697 Å | 63.8% | 88.1% | 22.6 | 200 s |
| v1 3-obj (confEnergy, RMSD, -Rg) | 0.814 Å | 0.545 Å | 72.3% | 93.0% | 25.2 | 223 s |
| v2 3-obj + dual archive | 0.732 Å | 0.536 Å | 74.8% | 94.5% | 31.6 | 417 s |
| **v2 + tuned defaults — current** | **0.638 Å** | **0.457 Å** | **80.5%** | **97.0%** | 56.1 | 1065 s |

- **v1: 22% lower mean min-RMSD** vs original; improved on 164 molecules,
  regressed on 71, equal on 94. The largest gains are on flexible ligands
  (e.g. 1f0u 11.1→1.1 Å, 1hos 11.0→2.6 Å, 1r1h 4.6→1.6 Å).
- **v2: further 10% lower mean min-RMSD** vs v1 (0.814→0.732 Å) by adding an
  ε-NSGA-II style diversity archive; improved on 134 molecules, regressed on
  85. Notably fixes v1's worst failures (1rne 10.2→2.1 Å, 1mts 2.0→0.68 Å,
  1ian 2.1→1.1 Å). Runtime roughly doubles (the diversity archive roughly
  doubles the archive size; tuning `MOGA_Population_Size`/`MOGA_Max_Generation`
  can trade accuracy for speed).
- **Tuned defaults: further 13% lower** (0.732→0.638 Å) via a parameter sweep
  on a 20-molecule subset: 300 generations (vs 100) and mutation 0.2 (vs 0.1)
  gave the largest gains; a finer ε grid (2 0.2 0.1 2 vs 3 0.3 0.1 2) helped
  further; 500 generations gave no additional gain. ≤1.0 Å success rate
  74.8%→80.5%, ≤2.0 Å 94.5%→97.0%. Runtime ~1065 s for the full set (more
  conformers kept: 56.1/mol on average).
- All runs: 329/329 molecules processed, 0 failures, fixed seed 0.42.

Changes:
- **Objective f0 = conformation-dependent energy** (torsion + VDW +
  electrostatics). The original split VDW/torsion objectives were strongly
  correlated; electrostatics (H-bonds, charge distribution) was previously
  excluded from the search entirely. Full FF energy is NOT used because
  stretch/bend/oop terms are constant under dihedral rotation and flatten the
  Pareto front into a single epsilon box.
- **f1 = aligned RMSD** vs. input conformer (unchanged role, now the only
  diversity objective in 2-obj mode).
- **Circular torsion encoding**: torsions are periodic; crossover maps parents
  onto the shortest arc and children wrap around ±180° instead of clamping.
- **v2: dual archive (ε-NSGA-II style)**: the original ε-dominance grid
  archive guarantees convergence; a second maximin (farthest-neighbour)
  diversity archive keeps spread non-dominated solutions that the ε grid
  rejected, and same-box competition uses normalized crowding distance
  (`obj_distance`/`min_archive_distance`/`update_diversity_archive` in
  `ConGen.cpp`).
- **Fixed silent ε-parameter bug**: `MOGA_Epsilon_Quaternion` in the parameter
  file was never parsed (the parser only knew the four separate
  `MOGA_*_Epsilon` keys), so grid sizes silently fell back to code defaults.
- **Fixed polynomial mutation bugs**: `int val = ...` truncated the mutation
  step to 0/1 (destroying the distribution); the random redraw also assigned
  `randomperc()` ∈ [0,1) to an `int`, always yielding 0.
- Tuned defaults: `MOGA_Num_Objectives 3`,
  `MOGA_Max_Generation 300`, `MOGA_Mutation_Probability 0.2`,
  `MOGA_Epsilon_Quaternion 2 0.2 0.1 2`, `MOGA_Energy_Cutoff 60`
  (see `Cyndi/CyndiParam.in`). Tuning sweep on a 20-molecule subset:
  300 generations and higher mutation (0.2) cut mean min-RMSD 0.668→0.477 Å
  (−29%); a finer ε grid (2, 0.2) helped further; 500 generations gave no
  additional gain.
- Benchmark harness: `bench/runner.py`, `bench/analyze.py`. (The per-molecule
  `bench/analysis_*.txt` result files were removed in the 2026-09 cleanup --
  they came from the broken evaluator described above. `bench/RETUNE_2026-09-16.md`
  has the re-measured numbers.)
- Multiprocess batch driver: `bench/run_parallel.py` splits a multi-molecule
  mol2 across N independent Cyndi processes and merges the outputs in input
  order. With a fixed seed the merged output is bit-identical to a single
  run (per-molecule seeds are derived from the molecule's global index via
  the `-startidx`/`-maxmols` CLI options). Example:

      python bench/run_parallel.py -input 329_test_set.mol2 \\
          -parm Cyndi/CyndiParam.in -n 4 -out all.mol2 -seed 0.42

  Measured on a 20-molecule subset: 68.2 s (1 process) vs 35.4 s (4
  processes) on 4 logical cores, identical output.

## License

Released under the [MIT License](LICENSE).
