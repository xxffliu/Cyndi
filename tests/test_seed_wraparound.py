#!/usr/bin/env python3
"""Molecules past global index ~580 must not crash the run.

`KNOWN_ISSUES.md` #4: the per-molecule seed is `BasicSeed_ + 0.001 * index`,
and the Knuth subtractive generator in `random.cpp` only corrects differences
below 0 -- never values at or above 1. With the documented seed 0.42 the seed
leaves `(0,1)` at molecule 581, after which `randomperc()` returns values of
order 1e10 and `(int)(rnds * PopSize_)` is an arbitrary array index. Before the
fix this was a reliable segmentation fault.

The fix reduces the seed with `fmod(..., 1.0)`. Nothing in the normal test set
reaches index 580, so this test builds a long input by concatenating the given
file until it is long enough, then processes a few molecules from beyond the
boundary.

Usage: test_seed_wraparound.py <cyndi-exe> <input.mol2> <param.in>
Exit codes: 0 pass, 1 fail, 77 skipped (input not available).
"""
import os
import re
import subprocess
import sys
import tempfile

START_INDEX = 650      # comfortably past the 581 boundary
N_MOLECULES = 3


def main():
    exe, molfile, parm = sys.argv[1], sys.argv[2], sys.argv[3]
    if not os.path.exists(molfile):
        print("SKIP: %s not present" % molfile)
        sys.exit(77)

    with open(molfile, errors="replace") as fh:
        text = fh.read()
    blocks = [b for b in re.split(r"(?=@<TRIPOS>MOLECULE)", text) if b.strip()]
    if not blocks:
        print("SKIP: no molecules in %s" % molfile)
        sys.exit(77)

    need = START_INDEX + N_MOLECULES
    repeats = -(-need // len(blocks))          # ceil
    workdir = tempfile.mkdtemp(prefix="cyndi_seed_")
    long_input = os.path.join(workdir, "long.mol2")
    with open(long_input, "w") as fh:
        fh.write(text * repeats)
    print("built a %d-molecule input (%d x %d)"
          % (len(blocks) * repeats, len(blocks), repeats))

    proc = subprocess.run(
        [exe, "-input", long_input, "-output", os.path.join(workdir, "out.mol2"),
         "-parm", parm, "-startidx", str(START_INDEX), "-maxmols", str(N_MOLECULES)],
        cwd=os.getcwd(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    if proc.returncode != 0:
        # A negative return code on POSIX is the killing signal; 0xC0000005 on
        # Windows is an access violation. Either way, name it.
        print("FAIL cyndi exited %d starting at molecule %d -- this is the seed "
              "wraparound crash if it looks like a signal or an access violation"
              % (proc.returncode, START_INDEX))
        print(proc.stdout[-4000:])
        sys.exit(1)

    if "Number of Molecules Processed: %d" % N_MOLECULES not in proc.stdout:
        print("FAIL run completed but did not report %d molecules processed"
              % N_MOLECULES)
        print(proc.stdout[-4000:])
        sys.exit(1)

    print("ok  molecules %d-%d processed without crashing"
          % (START_INDEX, START_INDEX + N_MOLECULES - 1))


if __name__ == "__main__":
    main()
