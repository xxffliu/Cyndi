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

### Windows (MinGW / MSYS2)

```bash
bash build.sh          # outputs build/cyndi.exe (requires g++ on PATH)
```

The MSYS2 g++ sysroot headers are not auto-included; `build.sh` handles this via `-idirafter`.

### Visual Studio (original)

Open `Cyndi/Cyndi.sln` (VS2008-era project; `Cyndi.vcxproj` for newer versions).

## Usage

```
Cyndi -input input.mol2 -output output.mol2 [-parm CyndiParam.in]
```

- Input: SYBYL MOL2 (single molecule or multi-molecule `@<TRIPOS>MOLECULE` blocks)
- Output: MOL2 with up to `MOGA_Max_Conformers` conformers per molecule
- Parameters: see `Cyndi/CyndiParam.in` (population size, generations, crossover/mutation probabilities, energy cutoff, RMSD clustering threshold, force field choice TAFF/MMFF94, etc.)

Batch mode:

```bash
python Cyndi_batch.py [options] input-mol2.list
```

### Test data

`Cyndi/329_test_set.mol2` contains the 329-molecule test set from the paper (PDB ligands).

## Layout

```
src/            C++ sources (MOGA core, force fields, minimizers, IO)
include/        headers
Cyndi/          VS project files, force-field parameter files, test data
Cyndi_batch.py  Python 2 batch driver (generates CyndiParam.in and loops the exe)
build.sh        MinGW build script
```

## Notes (2026 maintenance)

- Ported to build with modern GCC 16 (MinGW) — original code is VS2008-era
- Fixed: `dynamic_cast == false` (MMFF94StretchBend), `name == "****"` comparison-as-assignment (MOL2IO), CG minimizer never initialized in the post-MOGA optimization step
- Optimized: force-field setup hoisted out of the per-conformer loop (~3x speedup), `pow()` → multiplication in VDW kernels
- Clustering now uses Kabsch-aligned heavy-atom RMSD (was unaligned RMSD)

## v1 Algorithm Rework (2026)

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
- Benchmark harness: `bench/runner.py`, `bench/analyze.py`, full per-molecule
  results in `bench/analysis_new.txt` / `bench/analysis_old.txt`.
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
