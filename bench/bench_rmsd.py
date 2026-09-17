#!/usr/bin/env python3
"""Per-molecule benchmark with min-RMSD evaluation: for each named molecule,
generate conformers with a given Cyndi build, then compute the minimum aligned
heavy-atom RMSD between the generated set and the INPUT conformer (the
bioactive reference).  Lower min-RMSD = better reproduction.

Usage:  python bench_rmsd.py <source.mol2> <name1,name2,...> <cyndi-exe> <param.in>

The RMSD kernel lives in rmsd_core.py; see that file for why the
previous power-iteration eigenvalue solver in this script was wrong.
"""
import os
import re
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
from rmsd_core import parse_multi_mol2, kabsch_rmsd  # noqa: E402

# Cyndi's per-molecule table row: name, rot bonds, raw confs, kept confs, time.
ROW_RE = re.compile(r'^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.eE+-]+)\s*$', re.M)


def extract(src, names, outdir):
    txt = open(src, errors='ignore').read()
    got = {}
    for b in re.split(r'(?=@<TRIPOS>MOLECULE)', txt):
        m = re.search(r'@<TRIPOS>MOLECULE\s*\n(\S+)', b)
        if m and m.group(1) in names:
            got[m.group(1)] = b
    os.makedirs(outdir, exist_ok=True)
    for n, b in got.items():
        open(os.path.join(outdir, n + '.mol2'), 'w').write(b)
    return list(got)


def conformer_count(stdout, molname):
    """Kept-conformer count for `molname` from Cyndi's stdout table."""
    for m in ROW_RE.finditer(stdout):
        if m.group(1) == molname:
            return int(m.group(4))
    return 0


def run(exe, molfile, outfile, parm):
    t0 = time.time()
    r = subprocess.run([exe, '-input', molfile, '-output', outfile, '-parm', parm],
                       capture_output=True, text=True, timeout=600)
    return time.time() - t0, r.stdout


def main():
    src, names, exe, parm = sys.argv[1], sys.argv[2].split(','), sys.argv[3], sys.argv[4]
    outdir = 'bm_' + os.path.splitext(os.path.basename(exe))[0]
    got = extract(src, names, outdir)
    print('%-10s %8s %8s %9s %s' % ('mol', 'conf', 'time_s', 'minRMSD', 'status'))
    tot_conf, tot_t = 0, 0.0
    for n in got:
        mf = os.path.join(outdir, n + '.mol2')
        of = os.path.join(outdir, n + '_out.mol2')
        dt, out = run(exe, mf, of, parm)
        nconf = conformer_count(out, n)
        refs = parse_multi_mol2(mf)
        gens = parse_multi_mol2(of) if os.path.exists(of) else {}
        ref = refs.get(n, [None])[0]
        gen = gens.get(n, [])
        best = min((kabsch_rmsd(ref, g) for g in gen), default=float('nan')) \
            if ref else float('nan')
        tot_conf += nconf
        tot_t += dt
        print('%-10s %8d %8.3f %9.3f %s' % (n, nconf, dt, best, 'OK' if nconf else 'FAIL'))
    print('TOTAL: %d confs, %.2f s, avg %.3f s/mol'
          % (tot_conf, tot_t, tot_t / max(1, len(got))))


if __name__ == '__main__':
    main()
