#!/usr/bin/env python3
"""Run Cyndi on the 329-molecule test set with the given exe and parameters.

Usage: python runner.py <exe> <param.in> <out.mol2> <log.txt>

Cyndi opens MMFF94.parm and TAFF.parm by bare filename, i.e. relative to the
working directory, so this copies them next to the input set before running.
The canonical copies live in Cyndi/ -- there is deliberately only one of each
in the repository.
"""
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PARM_SOURCE = os.path.join(HERE, os.pardir, 'Cyndi')


def stage_parm_files(dest):
    for name in ('MMFF94.parm', 'TAFF.parm'):
        src = os.path.join(PARM_SOURCE, name)
        dst = os.path.join(dest, name)
        if not os.path.exists(dst):
            if not os.path.exists(src):
                sys.exit('missing force-field parameter file: %s' % src)
            shutil.copyfile(src, dst)


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    exe, parm, out, log = sys.argv[1:5]
    exe = os.path.abspath(exe)
    parm = os.path.abspath(parm)
    os.chdir(HERE)
    stage_parm_files(HERE)
    t0 = time.time()
    with open(log, 'w') as fh:
        r = subprocess.run([exe, '-input', '329_test_set.mol2',
                            '-output', out, '-parm', parm],
                           stdout=fh, stderr=subprocess.STDOUT, text=True)
    print('exit=%d elapsed=%.1fs' % (r.returncode, time.time() - t0))
    return r.returncode


if __name__ == '__main__':
    sys.exit(main())
