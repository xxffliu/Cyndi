#!/usr/bin/env python
"""Run Cyndi on the 329 test set with the given exe/params. Usage:
python runner.py <exe> <param> <out_mol2> <log>"""
import subprocess, sys, time, os

exe, parm, out, log = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
os.chdir(os.path.dirname(os.path.abspath(__file__)))
t0 = time.time()
with open(log, 'w') as fh:
    r = subprocess.run([os.path.abspath(exe), '-input', '329_test_set.mol2',
                        '-output', out, '-parm', parm],
                       stdout=fh, stderr=subprocess.STDOUT, text=True)
dt = time.time() - t0
print('exit=%d elapsed=%.1fs' % (r.returncode, dt))
