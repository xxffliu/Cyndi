#!/usr/bin/env python3
"""A chunked run must reproduce a whole-file run, molecule for molecule.

This is `KNOWN_ISSUES.md` #1, and it has broken twice for different reasons:

  * `GlobalIdx` was incremented at the end of the loop body, so every `continue`
    path (0 rotatable bonds, >30 rotatable bonds, setup failure, no conformers)
    skipped it -- and each skipped molecule shifted the RNG seed of every later
    molecule by one in a whole-file run but not in a chunked one.
  * `MOL2IO::write` inherited stream state (`adjustfield`, `floatfield`,
    `precision`, `showpoint`) left set by the previous molecule's atom section,
    so the FIRST molecule of every chunk was formatted differently from the
    same molecule inside a whole-file run.

Both are invisible in a 20-molecule subset, which is why they survived. The
test therefore runs far enough into the file to cross a skipped molecule
(`1byb`, #34 of the 329-molecule set, has more than 30 rotatable bonds).

Compares the `@<TRIPOS>ATOM` blocks, which is where the science is; the
`# Created and modified by GAPS on <date>` line is excluded because it is a
timestamp.

Usage: test_chunk_reproducibility.py <cyndi-exe> <input.mol2> <param.in>
Exit codes: 0 pass, 1 fail, 77 skipped (input not available).
"""
import os
import subprocess
import sys
import tempfile

# How many molecules the whole-file leg processes. Must be past the first
# skipped molecule for the test to have any power; 90 clears 1byb (#34) with
# room to spare and keeps the run to a couple of minutes.
N_MOLECULES = 90
# Molecules checked individually, given as 0-based global indices. 83 is 1glp,
# the molecule the original bug report was filed against.
SPOT_CHECK_INDICES = [83]


def records(path):
    """{(index, name): (atom lines,)} for every record in a mol2 file."""
    with open(path, errors="replace") as fh:
        lines = fh.read().split("\n")
    starts = [i for i, l in enumerate(lines) if l.startswith("@<TRIPOS>MOLECULE")]
    out = []
    for k, i in enumerate(starts):
        end = starts[k + 1] if k + 1 < len(starts) else len(lines)
        block = lines[i:end]
        name = block[1].strip() if len(block) > 1 else ""
        try:
            a = block.index("@<TRIPOS>ATOM")
            b = block.index("@<TRIPOS>BOND")
            atoms = tuple(block[a + 1:b])
        except ValueError:
            atoms = ()
        out.append((name, atoms))
    return out


def run(exe, args, cwd):
    proc = subprocess.run([exe] + args, cwd=cwd,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        print("cyndi exited %d\n%s" % (proc.returncode, proc.stdout[-4000:]))
        sys.exit(1)


def main():
    exe, molfile, parm = sys.argv[1], sys.argv[2], sys.argv[3]
    if not os.path.exists(molfile):
        print("SKIP: %s not present" % molfile)
        sys.exit(77)

    workdir = tempfile.mkdtemp(prefix="cyndi_chunk_")
    # The force-field parameter files are opened relative to the working
    # directory, so run where ctest put them.
    cwd = os.getcwd()

    whole = os.path.join(workdir, "whole.mol2")
    run(exe, ["-input", molfile, "-output", whole, "-parm", parm,
              "-maxmols", str(N_MOLECULES)], cwd)
    whole_records = records(whole)

    failures = []
    for idx in SPOT_CHECK_INDICES:
        solo = os.path.join(workdir, "solo_%d.mol2" % idx)
        run(exe, ["-input", molfile, "-output", solo, "-parm", parm,
                  "-startidx", str(idx), "-maxmols", "1"], cwd)
        solo_records = records(solo)
        if not solo_records:
            failures.append("molecule %d produced no output in the chunked run" % idx)
            continue
        name = solo_records[0][0]
        from_whole = [r for r in whole_records if r[0] == name]
        if len(from_whole) != len(solo_records):
            failures.append("%s: %d conformers whole-file vs %d chunked"
                            % (name, len(from_whole), len(solo_records)))
            continue
        differing = sum(1 for a, b in zip(from_whole, solo_records) if a[1] != b[1])
        if differing:
            failures.append("%s: %d of %d conformers differ in their atom blocks"
                            % (name, differing, len(solo_records)))
        else:
            print("ok  %s: %d conformers identical (whole-file vs -startidx %d)"
                  % (name, len(solo_records), idx))

    if failures:
        for f in failures:
            print("FAIL " + f)
        sys.exit(1)
    print("chunked and whole-file runs agree")


if __name__ == "__main__":
    main()
