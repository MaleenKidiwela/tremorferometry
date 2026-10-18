#!/usr/bin/env python
"""Stage 7 -- WHOLE-CHAIN SITE-DAY-SHIFT NULL (>= 3 realisations).

The ENTIRE pipeline (autocorrelation -> clusters -> handoff -> growth -> survival -> iterations ->
dedupe -> locate -> gates) is re-run with each SITE's day-series circularly shifted by a
site-specific whole-day offset.  Members of a site shift TOGETHER, so site-coherent noise is
preserved exactly and only true inter-site coherence is destroyed.  The count of null "gated
families" per realisation is the family-level false-alarm rate (accept: <= 1 per realisation;
>= 3 means the FA rule is broken and the run stops).

Each realisation runs as its own subprocess chain so the shift can never leak into the real pass.
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
import w_common as W                                                        # noqa: E402

PY = sys.executable


def sh(cmd):
    W.log("$ " + " ".join(cmd[:6]) + (" ..." if len(cmd) > 6 else ""))
    r = subprocess.run(cmd, cwd=os.path.dirname(HERE))
    if r.returncode:
        raise SystemExit(f"stage failed ({r.returncode}): {' '.join(cmd)}")


def realisation(k, T, out_root):
    shift = W.site_shift(k)
    d = f"{out_root}/null{k}"
    os.makedirs(d, exist_ok=True)
    js = json.dumps(shift)
    json.dump(shift, open(f"{d}/shift.json", "w"), indent=1)
    W.log(f"[null{k}] site shifts {shift}")
    sh([PY, f"{HERE}/w_autocorr.py", "--hours", "tremor", "--out", f"{d}/ac",
        "--t-record", str(W.T_ENG_FLOOR), "--shift", js])
    sh([PY, f"{HERE}/w_cluster.py", "cluster", "--pairs", f"{d}/ac.parquet",
        "--thr", f"{out_root}/T.json", "--out", f"{d}/cl"])
    sh([PY, f"{HERE}/w_grow.py", "handoff", "--clusters", f"{d}/cl_clusters.parquet",
        "--seeds", f"{d}/cl_seeds.csv", "--out", f"{d}/bank0.npz", "--shift", js])
    sh([PY, f"{HERE}/w_grow.py", "quiet", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--shift", js])
    sh([PY, f"{HERE}/w_grow.py", "grow", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--thr", f"{d}/p0_thresholds.csv", "--shift", js])
    sh([PY, f"{HERE}/w_grow.py", "scramble", "--clusters", f"{d}/cl_clusters.parquet",
        "--seeds", f"{d}/cl_seeds.csv", "--out", f"{d}/bank0scr.npz", "--shift", js])
    sh([PY, f"{HERE}/w_grow.py", "grow", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
        "--thr", f"{d}/p0_thresholds.csv", "--shift", js])
    sh([PY, f"{HERE}/w_grow.py", "survive", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--seeds", f"{d}/cl_seeds.csv", "--scramble-hits", f"{d}/p0scr_hits.parquet",
        "--scramble-bank", f"{d}/bank0scr.npz"])
    prev, tag = f"{d}/bank0.npz", "p0"
    for it in range(1, W.N_ITERATIONS + 1):
        sh([PY, f"{HERE}/w_iterate.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
            "--survivors", f"{d}/p0_survival.csv", "--out-bank", f"{d}/bank{it}.npz",
            "--out-diag", f"{d}/iter{it}_diag.csv", "--shift", js])
        sh([PY, f"{HERE}/w_grow.py", "quiet", "--bank", f"{d}/bank{it}.npz",
            "--out", f"{d}/p{it}", "--shift", js])
        sh([PY, f"{HERE}/w_grow.py", "grow", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
            "--thr", f"{d}/p{it}_thresholds.csv", "--shift", js])
        sh([PY, f"{HERE}/w_grow.py", "survive", "--bank", f"{d}/bank{it}.npz",
            "--out", f"{d}/p{it}", "--seeds", f"{d}/cl_seeds.csv"])
        prev, tag = f"{d}/bank{it}.npz", f"p{it}"
    sh([PY, f"{HERE}/w_dedupe.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--out", f"{d}/dedupe.csv"])
    sh([PY, f"{HERE}/w_locate_gates.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
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
    ap.add_argument("--n", type=int, default=W.N_SITESHIFT_REALISATIONS)
    ap.add_argument("--only", type=int, default=-1)
    a = ap.parse_args()
    T = json.load(open(f"{a.out_root}/T.json"))["T"]
    ks = [a.only] if a.only >= 0 else list(range(1, a.n + 1))
    rows = []
    for k in ks:
        rows.append(realisation(k, T, a.out_root))
    R = pd.DataFrame(rows)
    p = f"{a.out_root}/null_summary.csv"
    if os.path.exists(p):
        R = pd.concat([pd.read_csv(p), R], ignore_index=True).drop_duplicates("realisation",
                                                                              keep="last")
    R.to_csv(p, index=False)
    W.log(f"[null] {len(R)} realisations; gated families per realisation: "
          f"{list(R.n_gate_pass)}")


if __name__ == "__main__":
    main()
