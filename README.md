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

Rebuilt the MOGA objective design and torsion encoding, validated on a
20-molecule subset of the 329 test set (min-RMSD to the input/bioactive
conformer, lower is better):

| Configuration | mean min-RMSD | mean conformers | speed |
|---|---|---|---|
| original 4-obj (VDW, torsion, RMSD, -Rg) | 0.604 Å | 27.6 | 0.67 s/mol |
| v1 2-obj (confEnergy, RMSD) | 0.524 Å | 6.3 | 0.26 s/mol |
| **v1 3-obj (confEnergy, RMSD, -Rg) — default** | **0.503 Å** | 30.6 | 0.86 s/mol |

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
- **Fixed silent ε-parameter bug**: `MOGA_Epsilon_Quaternion` in the parameter
  file was never parsed (the parser only knew the four separate
  `MOGA_*_Epsilon` keys), so grid sizes silently fell back to code defaults.
- **Fixed polynomial mutation bugs**: `int val = ...` truncated the mutation
  step to 0/1 (destroying the distribution); the random redraw also assigned
  `randomperc()` ∈ [0,1) to an `int`, always yielding 0.
- Tuned defaults: `MOGA_Num_Objectives 3`, `MOGA_Epsilon_Quaternion 3 0.3 0.1 2`,
  `MOGA_Energy_Cutoff 60` (see `Cyndi/CyndiParam.in`).

## License

Released under the [MIT License](LICENSE).
