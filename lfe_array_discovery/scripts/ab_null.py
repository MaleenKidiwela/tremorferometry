#!/usr/bin/env python
"""STAGE 7 -- the WHOLE-CHAIN SITE-SHIFT NULL for the Royer reproduction.

`ab_common.NULL_K` (env `AB_NULL=k`) makes every stage read a station-specific circularly shifted
day WITHIN THE SAME EPISODE, with `t0` rewritten to the requested date so all stations stay on one
absolute clock.  Site-coherent noise, the diurnal cycle and each station's own gap pattern survive
untouched; only inter-station coherence -- what an LFE family is -- is destroyed.

The null gets the IDENTICAL chain and the IDENTICAL table (memory `rank-dont-discard`): the same
autocorrelation, the same cluster sweep, the same growth iterations, the same locator, the same
columns.  Evidence is the real-versus-null separation as a function of each criterion, never a
threshold chosen on the real run.

  python scripts/ab_null.py --sub NW --reals 3
-> results/ab_royer/<SUB>/null<k>/{auto_window,clusters,grow,locate}/...
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

PY = sys.executable


def sh(cmd, env):
    A.log("$ AB_NULL=%s %s" % (env.get("AB_NULL"), " ".join(cmd[1:5]) + " ..."))
    r = subprocess.run(cmd, cwd=os.path.dirname(HERE), env=env)
    if r.returncode:
        raise SystemExit(f"null stage failed ({r.returncode}): {' '.join(cmd)}")


def realisation(sub, k, grow_iters, boot):
    """The IDENTICAL literal chain under the site-shift null: same detector, the REAL run's
    clustering configuration (--headline-from, never re-selected on the null), same growth pool
    and scan days, same Hyp2000+hypoDD location, same secondary locator."""
    root = f"{A.RES}/{sub}/null{k}"
    os.makedirs(root, exist_ok=True)
    env = dict(os.environ, AB_NULL=str(k))
    sh([PY, f"{HERE}/ab_auto.py", "--sub", sub, "--cols", "window",
        "--out", f"{root}/auto_window"], env)
    if not os.path.exists(f"{root}/clusters/summary.json"):
        sh([PY, f"{HERE}/ab_cluster.py", "--sub", sub, "--auto", f"{root}/auto_window",
            "--out", f"{root}/clusters",
            "--headline-from", f"{A.RES}/{sub}/clusters/summary.json"], env)
    if not os.path.exists(f"{root}/grow/summary.json"):
        sh([PY, f"{HERE}/ab_grow.py", "--sub", sub, "--tpl", f"{root}/clusters/tpl_it0.npz",
            "--out", f"{root}/grow", "--max-iters", str(grow_iters), "--tol", "0.01"], env)
    its = sorted(int(f.split("it")[1].split(".")[0])
                 for f in os.listdir(f"{root}/grow") if f.startswith("det_it"))
    if not its:
        A.log(f"[null {k}] no detections to locate")
        return
    last = max(its)
    if not os.path.exists(f"{root}/hypo/loc_summary.csv"):
        sh([PY, f"{HERE}/ab_hypo.py", "--sub", sub, "--det", f"{root}/grow/det_it{last}.parquet",
            "--out", f"{root}/hypo", "--boot", str(boot)], env)
    if not os.path.exists(f"{root}/locate/depth_summary.csv"):
        sh([PY, f"{HERE}/ab_locate.py", "--sub", sub,
            "--det", f"{root}/grow/det_it{last}.parquet",
            "--tpl", f"{root}/grow/tpl_it{last}.npz",
            "--out", f"{root}/locate", "--boot", str(boot)], env)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
    ap.add_argument("--reals", type=int, default=3)
    ap.add_argument("--only", default="", help="run just this realisation index (e.g. --only 2)")
    ap.add_argument("--grow-iters", type=int, default=8)     # AUDIT FIX A9: = the real run
    ap.add_argument("--boot", type=int, default=60)          # AUDIT FIX A9: = the real run
    a = ap.parse_args()
    ks = [int(a.only)] if a.only else list(range(1, a.reals + 1))
    for k in ks:
        realisation(a.sub, k, a.grow_iters, a.boot)
    A.log(f"[null {a.sub}] realisations {ks} done")


if __name__ == "__main__":
    main()
