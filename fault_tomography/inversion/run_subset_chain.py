#!/usr/bin/env python
"""Run the daily joint-inversion chain (solver -> transmission_T -> transmission_deconv) on a STATION
SUBSET, in its own results directory. Built 2026-09-02 for the split-network coherence test on the
47.2N deep decline (and the roster-stable / cohort-leave-out variants). NOTHING in res_catalog is
touched: the subset directory gets symlinks to G.npz + era_table.csv and a FILTERED copy of
family_days.parquet / families_daily.csv; all three chain scripts read everything from RESDIR.

Usage:
  python fault_tomography/inversion/run_subset_chain.py --name halfA --include stations_A.txt \
      [--exclude stations_X.txt] [--lam 0.464] [--workers 14] [--stages solve,T,deconv]
Station tags are the `tag` column of families_daily.csv (e.g. B018p90f40, WISH, NLWA); one per line,
case-sensitive. --include/--exclude may be combined (include first, then exclude).
Output dir: fault_tomography/inversion/splits/<name>/  (+ chain.log)
"""
import argparse, os, subprocess, sys, time
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(HERE, "res_catalog")


def read_list(fn):
    return [x.strip() for x in open(fn) if x.strip() and not x.startswith("#")]


def prepare(name, include, exclude, outroot=None, drop_fe=None, parity=None, min_sites=None):
    """drop_fe: CSV with a column `fe` (tag|fam|e<era>) [+ optional `flagged` bool] of family-eras to DROP.
    parity: 'odd'|'even' -> keep only odd/even family-eras within every station (sorted by fe id) = design (E)."""
    out = os.path.join(outroot or os.path.join(HERE, "splits"), name)
    os.makedirs(out, exist_ok=True)
    for f in ("G.npz", "era_table.csv", "cells.csv", "pairs.csv"):
        dst = os.path.join(out, f)
        if not os.path.exists(dst):
            os.symlink(os.path.join(SRC, f), dst)
    fams = pd.read_csv(f"{SRC}/families_daily.csv")
    tags = set(fams.tag.unique())
    if include:
        tags &= set(include)
    if exclude:
        tags -= set(exclude)
    tags = sorted(tags)
    fams = fams[fams.tag.isin(tags)]
    fd = pd.read_parquet(f"{SRC}/family_days.parquet")
    fd = fd[fd.tag.isin(tags)]
    fams["fe"] = fams.tag + "|" + fams.fam + "|e" + fams.era.astype(str)
    fd["fe"] = fd.tag + "|" + fd.fam + "|e" + fd.era.astype(str)
    n_drop = 0
    if drop_fe:
        dcsv = pd.read_csv(drop_fe)
        bad = set(dcsv.fe[dcsv.flagged.astype(bool)]) if "flagged" in dcsv else set(dcsv.fe)
        n_drop = int(fd.fe.isin(bad).sum())
        fd = fd[~fd.fe.isin(bad)]; fams = fams[~fams.fe.isin(bad)]
    if parity:
        fams = fams.sort_values("fe")
        fams["k"] = fams.groupby("tag").cumcount()
        keep_fe = set(fams.fe[(fams.k % 2) == (1 if parity == "odd" else 0)])
        fd = fd[fd.fe.isin(keep_fe)]; fams = fams[fams.fe.isin(keep_fe)]
    fd = fd.drop(columns=["fe"]); fams = fams.drop(columns=[c for c in ("fe", "k") if c in fams])
    with open(f"{out}/subset_spec.txt", "w") as fh:
        fh.write(f"include={include is not None} exclude={exclude is not None} drop_fe={drop_fe} "
                 f"parity={parity} rows_dropped_by_fe={n_drop} min_sites={min_sites}\n")
    fd.to_parquet(f"{out}/family_days.parquet", index=False)
    fams.to_csv(f"{out}/families_daily.csv", index=False)
    with open(f"{out}/stations.txt", "w") as fh:
        fh.write("\n".join(tags) + "\n")
    print(f"[{name}] stations {len(tags)} | family-eras {len(fams)} | family-days {len(fd):,} "
          f"| span {fd.date.min().date()} -> {fd.date.max().date()}", flush=True)
    return out


def run(cmd, env, log):
    with open(log, "a") as fh:
        fh.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} $ {' '.join(cmd)}\n"); fh.flush()
        rc = subprocess.call(cmd, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        raise SystemExit(f"stage failed rc={rc}: {' '.join(cmd)} (see {log})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--include", default=None)
    ap.add_argument("--exclude", default=None)
    ap.add_argument("--lam", type=float, default=0.464, help="fixed lam_f (0 = re-oracle)")
    ap.add_argument("--lam-ratio", type=float, default=0.125)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--max-sweeps", type=int, default=30)
    ap.add_argument("--stages", default="solve,T,deconv")
    ap.add_argument("--outroot", default=None, help="parent dir (default fault_tomography/inversion/splits)")
    ap.add_argument("--drop-fe", default=None, help="CSV of family-eras to drop (col fe, optional col flagged)")
    ap.add_argument("--parity", default=None, choices=["odd", "even"], help="design (E) family split")
    ap.add_argument("--min-sites", type=int, default=None, help="override solver MIN_SITES (env MIN_SITES)")
    args = ap.parse_args()
    inc = read_list(args.include) if args.include else None
    exc = read_list(args.exclude) if args.exclude else None
    out = prepare(args.name, inc, exc, args.outroot, args.drop_fe, args.parity, args.min_sites)
    env = dict(os.environ, RESDIR=os.path.relpath(out, ROOT), WORKERS=str(args.workers))
    if args.min_sites:
        env["MIN_SITES"] = str(args.min_sites)
    log = f"{out}/chain.log"
    stages = args.stages.split(",")
    py = sys.executable
    if "solve" in stages:
        run([py, f"{HERE}/invert_dvv_4d_daily.py", "--lam", str(args.lam), "--lam-ratio",
             str(args.lam_ratio), "--workers", str(args.workers),
             "--max-sweeps", str(args.max_sweeps)], env, log)
    if "T" in stages:
        run([py, f"{HERE}/transmission_T.py"], env, log)
    if "deconv" in stages:
        run([py, f"{HERE}/transmission_deconv.py"], env, log)
    print(f"[{args.name}] DONE -> {out}", flush=True)
