# Known Issues

## 1. [RESOLVED 2026-09-16] Chunked vs single-process reproducibility differs on some molecules

**Root cause:** `GlobalIdx` in `execuateCONGEN()` was incremented at the very
*end* of the loop body, but four paths `continue` past that point: zero
rotatable bonds, more than 30 rotatable bonds, `congen.setup()` failure, and
`num_conf == 0`. The `-startidx` skip path (near the top of the loop) *did*
increment it correctly. So every molecule after the first skipped one received
a global index — and therefore an RNG seed — one lower in a whole-file run than
in a chunked run.

In the 329-molecule test set `1byb` (molecule #34) has more than 30 rotatable
bonds and is skipped, so from `1glp` (#84) onward the two invocations used
different seeds. That is exactly the reported pattern: the 20-molecule subset
was reproducible (it contains no skipped molecule) while the full set diverged
from `1glp` on.

**Fix:** capture the index before the body runs and advance immediately:

```cpp
Counter += 1;
int ThisIdx = GlobalIdx;   // index of THIS molecule
GlobalIdx += 1;            // advance now, so every 'continue' path also counts
...
congen.set_global_idx(ThisIdx);
```

**Verified:** with `MOGA_Random_Seed 0.42`, `1glp` produces 11 conformers whose
`@<TRIPOS>ATOM` blocks are byte-identical between a whole-file run and
`-startidx 83 -maxmols 1`.

A second, smaller reproducibility defect was fixed at the same time in
`MOL2IO::write`: the atom-count line inherited the stream's `adjustfield`,
which the atom/bond writer leaves set to `left`. The first molecule written to
any file therefore got a right-adjusted count line and all later ones a
left-adjusted one — so the first molecule of every chunk formatted differently
from the same molecule inside a whole-file run. The writer now sets the
adjustfield explicitly.

The three hypotheses listed in the previous version of this file
(uninitialised memory, static force-field state, `MOL2IO` read state) can all
be closed. Two of them did point at real, separate defects, which are also
fixed now: `individual::gyration_radius` was never initialised and
`individual::operator=` dropped `cv` (`ConGen.cpp`), and `MOGA::clear()` reset
`taff_` but not `mmff94_` — the *default* force field — nor `div_archive`.

---

## 2. [RESOLVED] `individual::f[]` out-of-bounds read in archive dump

Fixed in `53c873d`. `execuateMOGA` read `f[2]`/`f[3]` unconditionally even
when `NumObjects_ < 4`; the `f` vector is sized `NumObjects_`, so this was
heap UB (read stale garbage). Single-threaded runs happened to read stable
zeros, so results were unaffected.

Superseded 2026-09-16: the `Conformer` field mapping was still the original
four-objective layout, so `Conformer::rmsd` was being filled from `f[2]`
(which under the v1 layout is `-gyration radius`, not the RMSD) and
`TorsionEnergy` from `f[1]` (the RMSD). That only mattered on the
`MOGA_Optimize_Conformer N` path, which writes `TorVec[i].rmsd` straight into
the output; it now reads `f[1]` and `individual::gyration_radius`.

## 3. [RESOLVED] Global RNG state leaked between molecules

Original `execuateMOGA` called `warmup_random(seed)` without zeroing
`oldrand[]` and without resetting `jrand`, so the stream for a given seed
depended on how many random numbers previous molecules consumed. Now uses
`randomize(seed)` (zeroes + warmup + jrand reset).

Also fixed 2026-09-16: `oldrand[55]`, `rndx1`, `rndx2` and `rndcalcflag` were
declared `static` **in `random.h`**, giving every translation unit that
includes it a private copy of the generator state. It happened to work only
because `randomize`, `warmup_random`, `advance_random` and `randomperc` all
live in `random.cpp` and so all touched that file's copy. They are now
`extern` in the header and defined once in `random.cpp`.

## 4. [RESOLVED 2026-09-16] Segfault on inputs with more than ~580 molecules

`eff_seed = BasicSeed_ + 0.001 * mo_global_idx_` walks out of the `(0,1)`
interval that the Knuth subtractive generator in `random.cpp` requires — it
only corrects differences below 0, never values above 1. With the documented
seed `0.42` this happens at molecule 581, after which `randomperc()` returns
values of order 1e10 and `(int)(rnds * PopSize_)` becomes an arbitrary index.
Reproduced by concatenating the 329-molecule set three times and running
`-startidx 650 -maxmols 3`: segmentation fault before the fix, clean run after.
The seed is now reduced with `fmod(..., 1.0)` (identical results for
`idx < 580`, so existing benchmarks are unaffected).

Longer term this generator should be replaced by `std::mt19937_64` seeded with
`splitmix64(global_index)`, which removes the whole class of problem — but that
changes every existing result, so it needs its own baseline.

## 5. [OPEN] CG minimisation aborts on a large fraction of conformers

Over the 329-molecule test set, 1520 conformer minimisations end with
`aborted_ = true` ("no convergence and the step computation failed") out of
roughly 22 000 — and on some inputs, such as `EGFR_ah.mol2`, essentially every
minimisation aborts. An aborted minimisation leaves the conformer at whatever
coordinates the line search reached, so it is not a crash and the conformer is
still energy-filtered afterwards, but it means the `ConicLineSearch` port from
BALL is failing to bracket a step on a large class of inputs.

This has been made visible rather than fixed: the per-conformer console message
is now a counter reported once at the end of the run. Diagnosing it properly
(and the obvious alternative — replacing the conic/cubic/quadratic
interpolation cascade with Armijo backtracking or L-BFGS) is the next piece of
optimiser work.

## 6. [OPEN] `MOL` copy semantics are incomplete

`MOL`'s copy constructor and assignment operator do not copy `bk_pos_`,
`rmsd_`, `rotation_map`, `rotor_id_map`, `_vring`, `charge_type_` or
`is_initialized_`, and `_vfragment` is a shallow copy of owning raw pointers
(the destructor deletes them, so a copy of a molecule with fragments is a
double free). `ATOM::operator=` copies `bond_` as raw pointers into the *source*
molecule's bonds.

Nothing downstream currently calls `get_bond_list()` or `apply_rotor()` on a
copy, so this does not bite today — `Cyndi.cpp`'s `FailedVector`, the
`priority_queue<MOL>` and `LeaderClustering` only read coordinates. Any new
code that does will. The fix is either a full deep copy with the adjacency
lists rebuilt afterwards, or making `MOL` non-copyable with an explicit
`clone()`.

`ForceField`, `MMFF94` and `TAFF` had the same shape of defect and have been
made non-copyable (`= delete`) as of 2026-09-16.

## 7. [RESOLVED 2026-09-17] `Cyndi_batch.py` was Python 2

`print` statements; it did not run on any current interpreter. Removed in the
2026-09-17 repository cleanup. `bench/run_parallel.py` does the same job
(and splits the work across processes).
