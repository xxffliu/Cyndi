#!/usr/bin/env python
"""Analyze Cyndi benchmark: for each molecule in the input set, compute the
minimum aligned RMSD between the input (bioactive) conformer and the generated
conformer set. Reports per-molecule and summary statistics."""
import re, sys, math, os

def parse_multi(fn):
    """Parse multi-mol2: return {name: [conformer_atom_lists, ...]}"""
    txt = open(fn, errors='ignore').read()
    blocks = re.split(r'@<TRIPOS>MOLECULE', txt)[1:]
    mols = {}
    for b in blocks:
        lines = b.strip().splitlines()
        if not lines:
            continue
        name = lines[0].strip()
        atoms = []
        in_atom = False
        for ln in lines[1:]:
            if ln.startswith('@<TRIPOS>ATOM'):
                in_atom = True; continue
            if ln.startswith('@<TRIPOS>'):
                in_atom = False; continue
            if in_atom:
                p = ln.split()
                if len(p) >= 6:
                    atoms.append((p[1], float(p[2]), float(p[3]), float(p[4])))
        mols.setdefault(name, []).append(atoms)
    return mols

def kabsch_rmsd(a, b):
    ha = [(x, y, z) for (n, x, y, z) in a if n[0] != 'H']
    hb = [(x, y, z) for (n, x, y, z) in b if n[0] != 'H']
    n = min(len(ha), len(hb))
    if n == 0:
        return float('nan')
    ca = tuple(sum(p[i] for p in ha) / n for i in range(3))
    cb = tuple(sum(p[i] for p in hb) / n for i in range(3))
    va = [(p[0]-ca[0], p[1]-ca[1], p[2]-ca[2]) for p in ha[:n]]
    vb = [(p[0]-cb[0], p[1]-cb[1], p[2]-cb[2]) for p in hb[:n]]
    R = [[0.0]*3 for _ in range(3)]
    for k in range(n):
        for i in range(3):
            for j in range(3):
                R[i][j] += va[k][i] * vb[k][j]
    F = [[R[0][0]+R[1][1]+R[2][2], R[1][2]-R[2][1], R[2][0]-R[0][2], R[0][1]-R[1][0]],
         [R[1][2]-R[2][1], R[0][0]-R[1][1]-R[2][2], R[0][1]+R[1][0], R[0][2]+R[2][0]],
         [R[2][0]-R[0][2], R[0][1]+R[1][0], -R[0][0]+R[1][1]-R[2][2], R[1][2]+R[2][1]],
         [R[0][1]-R[1][0], R[0][2]+R[2][0], R[1][2]+R[2][1], -R[0][0]-R[1][1]+R[2][2]]]
    q = [1.0, 0.0, 0.0, 0.0]
    for _ in range(200):
        nq = [sum(F[i][j]*q[j] for j in range(4)) for i in range(4)]
        nm = math.sqrt(sum(x*x for x in nq))
        if nm < 1e-300:
            break
        q = [x/nm for x in nq]
    lam = sum(q[i]*sum(F[i][j]*q[j] for j in range(4)) for i in range(4))
    ss = sum(sum(p[i]**2 for i in range(3)) for p in va) + \
         sum(sum(p[i]**2 for i in range(3)) for p in vb)
    return math.sqrt(max(0.0, (ss - 2*lam)/n))

def analyze(input_fn, out_fn, label):
    refs = parse_multi(input_fn)
    gens = parse_multi(out_fn)
    print('%-6s %6s %8s %8s' % ('mol', 'nconf', 'minRMSD', 'refAtoms'))
    rows = []
    for name, reflist in sorted(refs.items()):
        ref = reflist[0]
        gen = gens.get(name, [])
        if gen:
            best = min(kabsch_rmsd(ref, g) for g in gen)
        else:
            best = float('nan')
        rows.append((name, len(gen), best, len(ref)))
        print('%-6s %6d %8.3f %8d' % (name, len(gen), best, len(ref)))
    vals = [r[2] for r in rows if not math.isnan(r[2])]
    nconfs = [r[1] for r in rows]
    if vals:
        print()
        print('[%s] molecules: %d, evaluated: %d' % (label, len(rows), len(vals)))
        print('[%s] mean min-RMSD: %.3f A   median: %.3f A' % (label, sum(vals)/len(vals), sorted(vals)[len(vals)//2]))
        print('[%s] min-RMSD <= 1.0 A: %d/%d (%.1f%%)' % (label, sum(1 for v in vals if v <= 1.0), len(vals), 100.0*sum(1 for v in vals if v <= 1.0)/len(vals)))
        print('[%s] min-RMSD <= 2.0 A: %d/%d (%.1f%%)' % (label, sum(1 for v in vals if v <= 2.0), len(vals), 100.0*sum(1 for v in vals if v <= 2.0)/len(vals)))
        print('[%s] mean conformers/mol: %.1f   total conformers: %d' % (label, sum(nconfs)/len(nconfs), sum(nconfs)))
    return rows

if __name__ == '__main__':
    input_fn, out_fn, label = sys.argv[1], sys.argv[2], sys.argv[3]
    analyze(input_fn, out_fn, label)
