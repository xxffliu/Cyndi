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

## License

Released under the [MIT License](LICENSE).
