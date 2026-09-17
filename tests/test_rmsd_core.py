#!/usr/bin/env python3
"""bench/rmsd_core.py must return the true minimum RMSD.

The evaluator this replaced took the largest eigenvalue of the quaternion
matrix F by power iteration, which converges to the eigenvalue of largest
MAGNITUDE. F is symmetric and traceless, so it always has negative eigenvalues,
and whenever the most negative one dominated the reported RMSD was far too
large -- 32% of conformer pairs in one 329-molecule run, worst case 14.7 A.

Three checks, none of which the old implementation passes:

  1. A rigidly rotated + translated copy of a structure has RMSD 0.
  2. The numpy-free fallback agrees with the numpy path.
  3. A near-mirror-image pair, the case that made power iteration pick the
     wrong eigenvalue, gives the same answer both ways and is bounded by the
     unaligned RMSD.

Exit codes: 0 pass, 1 fail, 77 skipped (rmsd_core.py or numpy unavailable).
"""
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "bench"))

try:
    import rmsd_core
except ImportError as exc:
    print("SKIP: cannot import rmsd_core (%s)" % exc)
    sys.exit(77)

if rmsd_core._np is None:
    print("SKIP: numpy not available, nothing to cross-check the fallback against")
    sys.exit(77)

np = rmsd_core._np


def rotation(ax, ay, az):
    def rx(t):
        return np.array([[1, 0, 0], [0, math.cos(t), -math.sin(t)], [0, math.sin(t), math.cos(t)]])

    def ry(t):
        return np.array([[math.cos(t), 0, math.sin(t)], [0, 1, 0], [-math.sin(t), 0, math.cos(t)]])

    def rz(t):
        return np.array([[math.cos(t), -math.sin(t), 0], [math.sin(t), math.cos(t), 0], [0, 0, 1]])

    return rz(az) @ ry(ay) @ rx(ax)


def fallback_rmsd(a, b):
    """kabsch_rmsd with the numpy path disabled."""
    saved = rmsd_core._np
    rmsd_core._np = None
    try:
        return rmsd_core.kabsch_rmsd(a, b)
    finally:
        rmsd_core._np = saved


def main():
    random.seed(20260917)
    failures = []

    # 1. rigid motion must give 0
    for trial in range(20):
        n = random.randint(6, 40)
        pts = np.array([[random.uniform(-8, 8) for _ in range(3)] for _ in range(n)])
        R = rotation(random.uniform(0, 6.28), random.uniform(0, 6.28), random.uniform(0, 6.28))
        moved = pts @ R.T + np.array([random.uniform(-20, 20) for _ in range(3)])
        v = rmsd_core.kabsch_rmsd([tuple(p) for p in pts], [tuple(q) for q in moved])
        if v > 1e-6:
            failures.append("rigid motion trial %d gave RMSD %.6f, expected 0" % (trial, v))

    # 2. fallback agrees with numpy
    for trial in range(20):
        n = random.randint(6, 30)
        a = [tuple(random.uniform(-6, 6) for _ in range(3)) for _ in range(n)]
        b = [tuple(random.uniform(-6, 6) for _ in range(3)) for _ in range(n)]
        va = rmsd_core.kabsch_rmsd(a, b)
        vb = fallback_rmsd(a, b)
        if abs(va - vb) > 1e-6:
            failures.append("trial %d: numpy %.6f vs numpy-free fallback %.6f" % (trial, va, vb))

    # 3. near-mirror images: F's most negative eigenvalue dominates here, which
    #    is exactly where power iteration returned a wrong (inflated) RMSD.
    for trial in range(20):
        n = random.randint(8, 30)
        pts = np.array([[random.uniform(-6, 6) for _ in range(3)] for _ in range(n)])
        mirrored = pts * np.array([1.0, 1.0, -1.0])
        mirrored = mirrored + np.array([[random.gauss(0, 0.05) for _ in range(3)] for _ in range(n)])
        a = [tuple(p) for p in pts]
        b = [tuple(q) for q in mirrored]
        va = rmsd_core.kabsch_rmsd(a, b)
        vb = fallback_rmsd(a, b)
        if abs(va - vb) > 1e-6:
            failures.append("mirror trial %d: numpy %.6f vs fallback %.6f" % (trial, va, vb))
        # The optimally superposed RMSD can never exceed the RMSD of the two
        # centred structures with no rotation at all.
        ca = pts - pts.mean(0)
        cb = mirrored - mirrored.mean(0)
        unaligned = math.sqrt(float(((ca - cb) ** 2).sum()) / n)
        if va > unaligned + 1e-9:
            failures.append("mirror trial %d: aligned RMSD %.4f exceeds unaligned %.4f"
                            % (trial, va, unaligned))

    if failures:
        for f in failures[:10]:
            print("FAIL " + f)
        print("%d failure(s)" % len(failures))
        sys.exit(1)
    print("rmsd_core: rigid motion, numpy/fallback agreement and mirror cases all pass")


if __name__ == "__main__":
    main()
