#!/usr/bin/env python3
"""Exact minimum heavy-atom Kabsch RMSD for Cyndi benchmarks.

Background (code review 2026-09-16, section 4)
----------------------------------------------
The previous `kabsch_rmsd` in bench/analyze.py and build/bench_rmsd.py obtained
the largest eigenvalue of the 4x4 quaternion matrix F by plain power iteration.
Power iteration converges to the eigenvalue of largest MAGNITUDE, not the
algebraically largest one.  F is symmetric with trace 0, so it always has
negative eigenvalues, and whenever the most negative one dominates (common when
two conformers differ a lot, or when the best superposition is close to a
mirror image) the routine returned a negative `lam` and a wildly inflated RMSD.
On an 84-molecule run, 268 of 844 conformer pairs were wrong by more than
1e-3 A, the worst by 14.6 A.

This module computes the algebraically largest eigenvalue exactly with
numpy.linalg.eigvalsh, and falls back -- when numpy is unavailable -- to power
iteration on the shifted matrix F + c*I (c = max abs row sum >= spectral
radius), which is positive semi-definite, so there "largest magnitude" and
"algebraically largest" coincide.

Everything here works on heavy atoms only, selected from the mol2 atom-type
column (field 6, e.g. "C.3"/"N.ar"/"H"), which is the element-bearing field --
the atom-name column can start with 'H' for non-hydrogens and vice versa.
"""
import math
import re

try:
    import numpy as _np
except ImportError:  # pragma: no cover - fallback path
    _np = None


def element_of(atom_type):
    """Element symbol from a mol2 SYBYL atom type ('C.3' -> 'C')."""
    return atom_type.split('.')[0]


def parse_multi_mol2(path, heavy_only=True):
    """Parse a (multi-)mol2 file into {name: [ (n_atoms, 3) coordinate lists ]}.

    Conformers of one molecule keep the file's atom order, so conformer i and
    conformer j are index-aligned and no atom matching is needed.
    """
    txt = open(path, errors='ignore').read()
    mols = {}
    for block in re.split(r'@<TRIPOS>MOLECULE', txt)[1:]:
        lines = block.strip().splitlines()
        if not lines:
            continue
        name = lines[0].strip()
        coords = []
        in_atom = False
        for ln in lines[1:]:
            if ln.startswith('@<TRIPOS>ATOM'):
                in_atom = True
                continue
            if ln.startswith('@<TRIPOS>'):
                in_atom = False
                continue
            if not in_atom:
                continue
            p = ln.split()
            if len(p) < 6:
                continue
            if heavy_only and element_of(p[5]).upper() == 'H':
                continue
            coords.append((float(p[2]), float(p[3]), float(p[4])))
        if coords:
            mols.setdefault(name, []).append(coords)
    return mols


def _max_eig_sym4(F):
    """Algebraically largest eigenvalue of a symmetric 4x4 matrix."""
    if _np is not None:
        return float(_np.linalg.eigvalsh(_np.asarray(F, dtype=float)).max())
    # numpy-free fallback: shift into the PSD cone so that the largest
    # magnitude eigenvalue IS the algebraically largest one.
    c = max(sum(abs(v) for v in row) for row in F)
    G = [[F[i][j] + (c if i == j else 0.0) for j in range(4)] for i in range(4)]
    q = [1.0, 0.5, 0.25, 0.125]   # generic start: never orthogonal to the top vector
    for _ in range(500):
        nq = [sum(G[i][j] * q[j] for j in range(4)) for i in range(4)]
        nm = math.sqrt(sum(x * x for x in nq))
        if nm < 1e-300:
            break
        q = [x / nm for x in nq]
    lam = sum(q[i] * sum(G[i][j] * q[j] for j in range(4)) for i in range(4))
    return lam - c


def kabsch_rmsd(a, b):
    """Minimum RMSD between two index-aligned coordinate lists after optimal
    translation + rotation (no reflection)."""
    n = min(len(a), len(b))
    if n == 0:
        return float('nan')
    if _np is not None:
        A = _np.asarray(a[:n], dtype=float)
        B = _np.asarray(b[:n], dtype=float)
        va = A - A.mean(0)
        vb = B - B.mean(0)
        R = va.T @ vb
        ss = float((va ** 2).sum() + (vb ** 2).sum())
    else:
        ca = [sum(p[i] for p in a[:n]) / n for i in range(3)]
        cb = [sum(p[i] for p in b[:n]) / n for i in range(3)]
        va = [[p[i] - ca[i] for i in range(3)] for p in a[:n]]
        vb = [[p[i] - cb[i] for i in range(3)] for p in b[:n]]
        R = [[sum(va[k][i] * vb[k][j] for k in range(n)) for j in range(3)]
             for i in range(3)]
        ss = sum(sum(x * x for x in p) for p in va) + \
             sum(sum(x * x for x in p) for p in vb)
    F = [[R[0][0] + R[1][1] + R[2][2], R[1][2] - R[2][1], R[2][0] - R[0][2], R[0][1] - R[1][0]],
         [R[1][2] - R[2][1], R[0][0] - R[1][1] - R[2][2], R[0][1] + R[1][0], R[0][2] + R[2][0]],
         [R[2][0] - R[0][2], R[0][1] + R[1][0], -R[0][0] + R[1][1] - R[2][2], R[1][2] + R[2][1]],
         [R[0][1] - R[1][0], R[0][2] + R[2][0], R[1][2] + R[2][1], -R[0][0] - R[1][1] + R[2][2]]]
    lam = _max_eig_sym4(F)
    return math.sqrt(max(0.0, (ss - 2.0 * lam) / n))


def min_rmsd_per_molecule(ref_path, gen_path):
    """[(name, n_conformers, min_rmsd, n_heavy_atoms), ...] sorted by name."""
    refs = parse_multi_mol2(ref_path)
    gens = parse_multi_mol2(gen_path)
    rows = []
    for name in sorted(refs):
        ref = refs[name][0]
        gen = gens.get(name, [])
        best = min((kabsch_rmsd(ref, g) for g in gen), default=float('nan'))
        rows.append((name, len(gen), best, len(ref)))
    return rows


def summarize(rows, label, stream=None):
    import sys
    out = stream or sys.stdout
    vals = [r[2] for r in rows if not math.isnan(r[2])]
    nconfs = [r[1] for r in rows]
    if not vals:
        print('[%s] nothing evaluated' % label, file=out)
        return {}
    vals_sorted = sorted(vals)
    stats = {
        'label': label,
        'n': len(rows),
        'evaluated': len(vals),
        'mean': sum(vals) / len(vals),
        'median': vals_sorted[len(vals) // 2],
        'le1': 100.0 * sum(1 for v in vals if v <= 1.0) / len(vals),
        'le2': 100.0 * sum(1 for v in vals if v <= 2.0) / len(vals),
        'mean_confs': sum(nconfs) / len(nconfs),
        'total_confs': sum(nconfs),
    }
    print('[%(label)s] n=%(n)d evaluated=%(evaluated)d mean=%(mean).3f A '
          'median=%(median).3f A <=1A=%(le1).1f%% <=2A=%(le2).1f%% '
          'mean confs=%(mean_confs).1f total confs=%(total_confs)d' % stats, file=out)
    return stats


if __name__ == '__main__':
    import sys
    ref, gen = sys.argv[1], sys.argv[2]
    label = sys.argv[3] if len(sys.argv) > 3 else gen
    rows = min_rmsd_per_molecule(ref, gen)
    if '-v' in sys.argv:
        print('%-8s %6s %9s %6s' % ('mol', 'nconf', 'minRMSD', 'atoms'))
        for r in rows:
            print('%-8s %6d %9.3f %6d' % r)
        print()
    summarize(rows, label)
