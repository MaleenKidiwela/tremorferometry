#!/usr/bin/env python
"""Stage 7 -- WHOLE-CHAIN SITE-DAY-SHIFT NULL on the Washington bed (>= 3 realisations).

Identical in design to `w_null_siteshift.py` (each SITE's day series circularly shifted by a
site-specific whole-day offset, members shifting together, so site-coherent noise is preserved
exactly and only true inter-site coherence is destroyed); it differs only in that every stage is
launched through `x_stage.py` / `x_locate_gates.py`, so the null runs on the WA bed.  Each
realisation is its own subprocess chain, so a shift can never leak into the real pass.

Accept (Merlin Q5(b)): the real gated count must be >= 5 x the null mean AND the null mean must be
<= 1 gated "family" per realisation; >= 3 per realisation means the FA rule is broken and the run
stops.

    python scripts/x_null_siteshift.py --out-root results/step5 --n 3
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import x_common as X                                                        # noqa: E402

PY = sys.executable


def sh(cmd):
    X.log("$ " + " ".join(os.path.basename(c) for c in cmd[:5]) + (" ..." if len(cmd) > 5 else ""))
    r = subprocess.run(cmd, cwd=os.path.dirname(HERE))
    if r.returncode:
        raise SystemExit(f"null stage failed ({r.returncode}): {' '.join(cmd)}")


def realisation(k, out_root, hours_file=""):
    W = X.activate()
    shift = W.site_shift(k)
    d = f"{out_root}/null{k}"
    os.makedirs(d, exist_ok=True)
    js = json.dumps(shift)
    json.dump(shift, open(f"{d}/shift.json", "w"), indent=1)
    X.log(f"[null{k}] site shifts {shift}")
    ST = [PY, f"{HERE}/x_stage.py"]
    hrs = (["--hours", "file", "--hours-file", hours_file] if hours_file
           else ["--hours", "tremor"])
    if not os.path.exists(f"{d}/ac.parquet"):
        sh(ST + ["w_autocorr"] + hrs + ["--out", f"{d}/ac", "--t-record",
                                        str(W.T_ENG_FLOOR), "--shift", js])
    if not os.path.exists(f"{d}/cl_seeds.csv"):
        sh(ST + ["w_cluster", "cluster", "--pairs", f"{d}/ac.parquet",
                 "--thr", f"{out_root}/T.json", "--out", f"{d}/cl"])
    sh(ST + ["w_grow", "handoff", "--clusters", f"{d}/cl_clusters.parquet",
             "--seeds", f"{d}/cl_seeds.csv", "--out", f"{d}/bank0.npz", "--shift", js])
    sh(ST + ["w_grow", "quiet", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0", "--shift", js])
    sh(ST + ["w_grow", "grow", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
             "--thr", f"{d}/p0_thresholds.csv", "--shift", js])
    sh(ST + ["w_grow", "scramble", "--clusters", f"{d}/cl_clusters.parquet",
             "--seeds", f"{d}/cl_seeds.csv", "--out", f"{d}/bank0scr.npz", "--shift", js])
    sh(ST + ["w_grow", "grow", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
             "--thr", f"{d}/p0_thresholds.csv", "--shift", js])
    sh(ST + ["w_grow", "survive", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
             "--seeds", f"{d}/cl_seeds.csv", "--scramble-hits", f"{d}/p0scr_hits.parquet",
             "--scramble-bank", f"{d}/bank0scr.npz"])
    prev, tag = f"{d}/bank0.npz", "p0"
    for it in range(1, W.N_ITERATIONS + 1):
        sh(ST + ["w_iterate", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
                 "--survivors", f"{d}/p0_survival.csv", "--out-bank", f"{d}/bank{it}.npz",
                 "--out-diag", f"{d}/iter{it}_diag.csv", "--shift", js])
        sh(ST + ["w_grow", "quiet", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
                 "--shift", js])
        sh(ST + ["w_grow", "grow", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
                 "--thr", f"{d}/p{it}_thresholds.csv", "--shift", js])
        sh(ST + ["w_grow", "survive", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
                 "--seeds", f"{d}/cl_seeds.csv"])
        prev, tag = f"{d}/bank{it}.npz", f"p{it}"
    sh(ST + ["w_dedupe", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
             "--out", f"{d}/dedupe.csv"])
    sh([PY, f"{HERE}/x_locate_gates.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--dedupe", f"{d}/dedupe.csv", "--survival", f"{d}/{tag}_survival.csv",
        "--out", f"{d}/lg", "--shift", js])
    j = json.load(open(f"{d}/lg.json"))
    j["realisation"] = k
    j["shift"] = shift
    json.dump(j, open(f"{d}/summary.json", "w"), indent=1)
    return j


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--only", type=int, default=-1)
    ap.add_argument("--hours-file", default="")
    a = ap.parse_args()
    ks = [a.only] if a.only >= 0 else list(range(1, a.n + 1))
    rows = []
    for k in ks:
        p = f"{a.out_root}/null{k}/summary.json"
        rows.append(json.load(open(p)) if os.path.exists(p)
                    else realisation(k, a.out_root, a.hours_file))
    R = pd.DataFrame(rows)
    p = f"{a.out_root}/null_summary.csv"
    if os.path.exists(p):
        R = pd.concat([pd.read_csv(p), R], ignore_index=True).drop_duplicates("realisation",
                                                                              keep="last")
    R.to_csv(p, index=False)
    X.log(f"[null] {len(R)} realisations; gated families per realisation: {list(R.n_gate_pass)}")


if __name__ == "__main__":
    main()
