#!/usr/bin/env python
"""Tuning sweep: run each param variant on subset20, compute mean min-RMSD
and runtime, print a comparison table. Usage: python tune_sweep.py"""
import subprocess, os, re, sys, time, math, importlib.util

os.chdir(os.path.dirname(os.path.abspath(__file__)))
EXE = os.path.abspath('cyndi.exe')

# load analyze.py for kabsch + parsing
spec = importlib.util.spec_from_file_location('analyze', 'analyze.py')
an = importlib.util.module_from_spec(spec)
spec.loader.exec_module(an)

PARAMS = sorted(f for f in os.listdir('tune_params') if f.endswith('.in'))
print('%-16s %10s %10s %10s %8s' % ('variant', 'meanRMSD', 'median', 'conf', 'time_s'))

results = {}
for pf in PARAMS:
    name = pf[:-3]
    out = 'tune_out_%s.mol2' % name
    log = 'tune_log_%s.txt' % name
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run([EXE, '-input', 'subset20.mol2', '-output', out,
                            '-parm', os.path.join('tune_params', pf)],
                           stdout=fh, stderr=subprocess.STDOUT, text=True, timeout=900)
    dt = time.time() - t0
    # per-molecule min-RMSD vs input
    refs = an.parse_multi('subset20.mol2')
    gens = an.parse_multi(out) if os.path.exists(out) else {}
    vals, nconfs = [], 0
    for molname, reflist in refs.items():
        gen = gens.get(molname, [])
        nconfs += len(gen)
        if gen:
            vals.append(min(an.kabsch_rmsd(reflist[0], g) for g in gen))
    if vals:
        mean = sum(vals) / len(vals)
        med = sorted(vals)[len(vals) // 2]
    else:
        mean = med = float('nan')
    results[name] = (mean, med, nconfs, dt)
    print('%-16s %10.3f %10.3f %10d %8.1f' % (name, mean, med, nconfs, dt))

print()
best = min(results, key=lambda k: results[k][0])
print('BEST mean min-RMSD: %s (%.3f A)' % (best, results[best][0]))
