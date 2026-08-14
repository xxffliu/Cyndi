# Known Issues

## 1. [OPEN] Chunked vs single-process reproducibility differs on some molecules

**Status:** open -- the multiprocess driver (`bench/run_parallel.py`) is
correct for the 20-molecule subset (verified bit-identical) but still
diverges on some molecules of the full 329 set (e.g. `1glp`).

**Symptom.** With a fixed seed (`MOGA_Random_Seed 0.42`), the same molecule
yields different conformer sets depending on how the process reached it:

| invocation                                    | 1glp conformers (before filter) |
|-----------------------------------------------|--------------------------------|
| single run over the whole file                | 209                            |
| `-startidx 83 -maxmols 1` (chunked driver)    | 184 / 112                      |
| `-startidx 0 -maxmols 84`                     | 99 / 113                       |

The 20-molecule subset (`sub20_full.mol2`) is fully reproducible between a
single run and 4 chunked runs; the divergence appears only in the full 329
set, so it correlates with molecules processed earlier in the file.

**What has been established (2026-08):**

- Per-molecule seeds are derived deterministically:
  `seed = BasicSeed_ + 0.001 * mo_global_idx_`, where `mo_global_idx_` is the
  molecule's 0-based index in the whole file (`-startidx`/`-maxmols` only
  select the window; the global index is tracked inside Cyndi).
- `execuateMOGA` now calls `randomize(seed)` (zeroes `oldrand[55]` + warmup)
  instead of the original bare `warmup_random(seed)`, and `warmup_random`
  resets `jrand = 0`. This fixed the *stream*: debug builds print the same
  first random number (`V3DBG12 first=...`) for the same seed regardless of
  window, and even the same second random number (`V3DBG13`).
- **Yet the MOGA run itself still differs** (before-filter counts differ by
  up to 2x). So the RNG *stream* is identical, but something else in the
  per-molecule state still depends on processing history.
- `randomnormaldeviate`/`noise`/`rndcalcflag` are never called by Cyndi;
  `srand/rand` (via `TimeRandomSeed`) is not used when a fixed seed is set.
- FF (`MMFF94`/`TAFF`) and `MOL::initialize` do not call the RNG at all.

**Hypotheses not yet ruled out:**

1. Uninitialized memory in `individual` (constructor `f.resize(4)` but no
   value initialization) or elsewhere that happens to read deterministic
   garbage in one process layout and different garbage in another. The
   `f[]` out-of-bounds read was already fixed (`53c873d`) but other vectors
   (e.g. `box`) may still carry stale values across molecules.
2. Some static/global mutable state in a force-field component or in
   `MOL` that `clear()`/`setup()` does not fully reset.
3. `MOL2IO::read` internal state when skipping molecules (`-startidx` path
   reads and discards N molecules; the full-run path processes them).

**How to reproduce:**

```
# single window, whole file (slow, ~25 min on 4 cores)
python bench/run_parallel.py -input bench/329_test_set.mol2 \
    -parm bench/CyndiParam_tuned_seed.in -n 1 -seed 0.42 -out full.mol2

# chunked
python bench/run_parallel.py -input bench/329_test_set.mol2 \
    -parm bench/CyndiParam_tuned_seed.in -n 4 -seed 0.42 -out chunked.mol2

# compare molecule 1glp
python bench/analyze.py bench/329_test_set.mol2 full.mol2 full
python bench/analyze.py bench/329_test_set.mol2 chunked.mol2 chunked
```

Or fast single-molecule repro:

```
cyndi.exe -input 329_test_set.mol2 -output a.mol2 -parm seed.in -startidx 83 -maxmols 1
cyndi.exe -input 329_test_set.mol2 -output b.mol2 -parm seed.in -maxmols 84   # 1glp is #84
# a.mol2 and b.mol2 differ for 1glp even though seed/stream are identical
```

**Impact:** outputs of a chunked batch run may differ slightly from a single
run on a few molecules. Both are valid conformer sets; the difference is a
reproducibility/consistency issue, not a crash or a wrong-answer bug.

---

## 2. [RESOLVED] `individual::f[]` out-of-bounds read in archive dump

Fixed in `53c873d`. `execuateMOGA` read `f[2]`/`f[3]` unconditionally even
when `NumObjects_ < 4`; the `f` vector is sized `NumObjects_`, so this was
heap UB (read stale garbage). Single-threaded runs happened to read stable
zeros, so results were unaffected; guarded by `nobj >= 3` / `nobj >= 4` now.

## 3. [RESOLVED] Global RNG state leaked between molecules

Original `execuateMOGA` called `warmup_random(seed)` without zeroing
`oldrand[]` and without resetting `jrand`, so the stream for a given seed
depended on how many random numbers previous molecules consumed. Now uses
`randomize(seed)` (zeroes + warmup + jrand reset). This was required for
chunk-vs-full reproducibility and also makes every molecule's stream
independent, which is more correct.
