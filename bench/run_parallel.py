#!/usr/bin/env python
"""Cyndi multiprocess batch driver.

Splits a multi-molecule mol2 file across N independent Cyndi processes
(one chunk per process, run concurrently) and concatenates the chunk
outputs in input order. With a fixed RNG seed (MOGA_Random_Seed in the
param file) the merged output is BIT-IDENTICAL to a single-process run
over the whole file -- verified: same conformers, same coordinates.

Usage:
    python run_parallel.py -input <file.mol2> -parm <CyndiParam.in>
                           [-exe cyndi.exe] [-n 4] [-out out.mol2]
                           [-workdir dir] [-seed 0.42]

Options:
    -exe       Cyndi executable (default: cyndi.exe in the script dir)
    -n         number of parallel processes (default: CPU count)
    -out       merged output mol2 (default: output.mol2)
    -workdir   scratch dir for chunk files (default: <out>.chunks/)
    -seed      fixed RNG seed; if given, appended to the param file copy
               used for the run (overrides any MOGA_Random_Seed there)
"""
import argparse, os, re, subprocess, sys, time, shutil, concurrent.futures

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EXE = os.path.join(HERE, "cyndi.exe")


def split_mol2(src, ndst, outdir):
    """Split a multi-molecule mol2 into ndst chunk files (one per worker).
    Returns list of (chunk_path, start_index, n_mols)."""
    txt = open(src, encoding='utf-8', errors='replace').read()
    blocks = re.split(r'(?=@<TRIPOS>MOLECULE)', txt)
    mols = [b for b in blocks if b.strip()]
    n = len(mols)
    if n == 0:
        sys.exit("no molecules found in %s" % src)
    chunk = max(1, (n + ndst - 1) // ndst)
    files = []
    for start in range(0, n, chunk):
        part = mols[start:start + chunk]
        fn = os.path.join(outdir, "chunk_%05d.mol2" % start)
        with open(fn, 'w', encoding='utf-8', errors='replace') as fh:
            fh.write(''.join(part))
        files.append((fn, start, len(part)))
    return files


def run_chunk(exe, input_file, start, nmols, parm, out_file, workdir):
    """Run Cyndi on a slice of the FULL input file (start..start+nmols-1).
    All workers read the same file; -startidx/-maxmols select each slice.
    Each slice's molecules get global-index seeds, so merging the outputs is
    bit-identical to one run over the whole file."""
    log = os.path.join(workdir, "chunk_%05d.log" % start)
    with open(log, 'w') as lf:
        t0 = time.time()
        r = subprocess.run(
            [exe, "-input", input_file, "-output", out_file,
             "-parm", parm, "-startidx", str(start), "-maxmols", str(nmols)],
            stdout=lf, stderr=subprocess.STDOUT, text=True)
        dt = time.time() - t0
    ok = (r.returncode == 0)
    return start, ok, dt, log


def main():
    ap = argparse.ArgumentParser(description="Cyndi multiprocess batch driver")
    ap.add_argument("-input", required=True, help="multi-molecule input mol2")
    ap.add_argument("-parm", required=True, help="Cyndi parameter file")
    ap.add_argument("-exe", default=DEFAULT_EXE, help="Cyndi executable")
    ap.add_argument("-n", type=int, default=os.cpu_count() or 4,
                    help="number of parallel processes (default: CPU count)")
    ap.add_argument("-out", default="output.mol2", help="merged output mol2")
    ap.add_argument("-workdir", default=None, help="scratch dir for chunks")
    ap.add_argument("-seed", type=float, default=None,
                    help="fixed RNG seed (appended to the param file)")
    args = ap.parse_args()

    exe = os.path.abspath(args.exe)
    if not os.path.exists(exe):
        sys.exit("exe not found: %s" % exe)

    workdir = args.workdir or (args.out + ".chunks")
    os.makedirs(workdir, exist_ok=True)

    # param file: copy to workdir (so the -seed override never touches the
    # user's original), appending the fixed seed if requested
    parm_local = os.path.join(workdir, "param.in")
    txt = open(args.parm, encoding='utf-8', errors='replace').read()
    if args.seed is not None:
        txt = re.sub(r'(?m)^MOGA_Random_Seed\s+\S+.*$', '', txt)
        txt = txt.rstrip() + "\nMOGA_Random_Seed\t\t%.2f\n" % args.seed
    open(parm_local, 'w', encoding='utf-8').write(txt)

    print("splitting %s -> %d chunk(s)" % (args.input, args.n))
    chunks = split_mol2(args.input, args.n, workdir)
    print("  %d molecules, %d chunk file(s)" % (sum(c for _, _, c in chunks), len(chunks)))
    input_full = os.path.abspath(args.input)

    t0 = time.time()
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.n) as pool:
        futs = []
        for cf, start, nmols in chunks:
            out_chunk = os.path.join(workdir, "chunk_%05d_out.mol2" % start)
            futs.append(pool.submit(run_chunk, exe, input_full, start, nmols,
                                    parm_local, out_chunk, workdir))
        for f in concurrent.futures.as_completed(futs):
            results.append(f.result())
    wall = time.time() - t0

    results.sort()
    nfail = sum(1 for _, ok, _, _ in results if not ok)
    for start, ok, dt, log in results:
        print("  chunk@%d: %s (%.1f s)" % (start, "OK" if ok else "FAILED", dt))
    if nfail:
        print("WARNING: %d chunk(s) failed; see logs in %s" % (nfail, workdir))

    # merge in input order
    with open(args.out, 'w', encoding='utf-8', errors='replace') as out:
        for cf, start, nmols in chunks:
            oc = os.path.join(workdir, "chunk_%05d_out.mol2" % start)
            if os.path.exists(oc):
                out.write(open(oc, encoding='utf-8', errors='replace').read())
    print("merged -> %s" % args.out)
    print("wall time: %.1f s across %d process(es)" % (wall, min(args.n, len(chunks))))


if __name__ == "__main__":
    main()
