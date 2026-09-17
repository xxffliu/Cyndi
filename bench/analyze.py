#!/usr/bin/env python3
"""Analyze a Cyndi benchmark run: for each molecule in the input set, compute
the minimum aligned heavy-atom RMSD between the input (bioactive) conformer and
the generated conformer set.  Reports per-molecule and summary statistics.

Usage:  python analyze.py <input.mol2> <generated.mol2> <label> [-q]

The RMSD kernel lives in rmsd_core.py.  It replaces the power-iteration
eigenvalue solver this script used before the 2026-09-16 review: power
iteration returns the eigenvalue of largest MAGNITUDE, which for the quaternion
matrix F is frequently a negative one, and that inflated 32% of the pairwise
RMSDs (worst case 14.6 A).  See rmsd_core.py for the details.  Numbers produced
by older versions of this script are not comparable with these.
"""
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rmsd_core import min_rmsd_per_molecule, summarize  # noqa: E402


def analyze(input_fn, out_fn, label, verbose=True):
    rows = min_rmsd_per_molecule(input_fn, out_fn)
    if verbose:
        print('%-8s %6s %9s %7s' % ('mol', 'nconf', 'minRMSD', 'atoms'))
        for name, nconf, best, natoms in rows:
            print('%-8s %6d %9.3f %7d' % (name, nconf, best, natoms))
        print()
    summarize(rows, label)
    return rows


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if a != '-q']
    analyze(args[0], args[1], args[2] if len(args) > 2 else args[1],
            verbose='-q' not in sys.argv)
