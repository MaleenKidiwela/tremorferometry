#!/usr/bin/env python
"""STEP 1c' STAGE 9 -- WHOLE-CHAIN SITE-DAY-SHIFT NULL x3, counting GATED FAMILIES **and
SURVIVORS**.

Identical to `w_null_siteshift` (stages 4-8 re-run with each SITE's day series circularly shifted
by a site-specific whole-day offset, members shifting together, so site-coherent noise is preserved
and only true inter-site coherence is destroyed) with two changes Merlin's prescription 9 requires:

  * the per-realisation report carries the number of SURVIVORS as well as gated families -- with
    ~2,000 pair seeds per realisation the survivor count is the sensitive false-alarm readout;
  * every stage is guarded by its own output file, so a relaunched queue resumes instead of
    repeating GPU hours.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

PY = sys.executable


def sh(cmd, guard=""):
    if guard and os.path.exists(guard):
        W.log(f"[null] skip (guard {os.path.basename(guard)} exists)")
        return
    W.log("$ " + " ".join(os.path.basename(c) for c in cmd[:3]) + " ...")
    r = subprocess.run(cmd, cwd=os.path.dirname(HERE))
    if r.returncode:
        raise SystemExit(f"null stage failed ({r.returncode}): {' '.join(cmd)}")
    if guard and not os.path.exists(guard):
        raise SystemExit(f"null stage produced no {guard}: {' '.join(cmd)}")


def n_surv(p):
    S = W.read_csv_safe(p, W.SURV_COLS)
    return int(S.survives.astype(bool).sum()) if len(S) else 0


def realisation(k, out_root):
    shift = W.site_shift(k)
    d = f"{out_root}/null{k}"
    os.makedirs(d, exist_ok=True)
    js = json.dumps(shift)
    json.dump(shift, open(f"{d}/shift.json", "w"), indent=1)
    W.log(f"[null{k}] site shifts {shift}")
    G = f"{HERE}/w_grow.py"
    sh([PY, f"{HERE}/w_autocorr.py", "--hours", "tremor", "--out", f"{d}/ac",
        "--t-record", str(W.T_ENG_FLOOR), "--shift", js], f"{d}/ac.parquet")
    sh([PY, f"{HERE}/w_cluster.py", "cluster", "--pairs", f"{d}/ac.parquet",
        "--thr", f"{out_root}/T.json", "--out", f"{d}/cl"], f"{d}/cl_seeds.csv")
    sh([PY, G, "handoff", "--clusters", f"{d}/cl_clusters.parquet", "--seeds",
        f"{d}/cl_seeds.csv", "--out", f"{d}/bank0.npz", "--shift", js], f"{d}/bank0.npz")
    sh([PY, G, "quiet", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0", "--shift", js],
       f"{d}/p0_thresholds.csv")
    sh([PY, G, "grow", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--thr", f"{d}/p0_thresholds.csv", "--shift", js], f"{d}/p0_hits.parquet")
    sh([PY, G, "scramble", "--clusters", f"{d}/cl_clusters.parquet", "--seeds",
        f"{d}/cl_seeds.csv", "--out", f"{d}/bank0scr.npz", "--shift", js], f"{d}/bank0scr.npz")
    sh([PY, G, "grow", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
        "--thr", f"{d}/p0_thresholds.csv", "--shift", js], f"{d}/p0scr_hits.parquet")
    sh([PY, G, "survive", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--seeds", f"{d}/cl_seeds.csv", "--scramble-hits", f"{d}/p0scr_hits.parquet",
        "--scramble-bank", f"{d}/bank0scr.npz"], f"{d}/p0_survival.csv")
    prev, tag = f"{d}/bank0.npz", "p0"
    for it in range(1, W.N_ITERATIONS + 1):
        sh([PY, f"{HERE}/w_iterate.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
            "--survivors", f"{d}/p0_survival.csv", "--out-bank", f"{d}/bank{it}.npz",
            "--out-diag", f"{d}/iter{it}_diag.csv", "--shift", js], f"{d}/bank{it}.npz")
        sh([PY, G, "quiet", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}", "--shift", js],
           f"{d}/p{it}_thresholds.csv")
        sh([PY, G, "grow", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
            "--thr", f"{d}/p{it}_thresholds.csv", "--shift", js], f"{d}/p{it}_hits.parquet")
        sh([PY, G, "survive", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
            "--seeds", f"{d}/cl_seeds.csv"], f"{d}/p{it}_survival.csv")
        prev, tag = f"{d}/bank{it}.npz", f"p{it}"
    sh([PY, f"{HERE}/w_dedupe.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--out", f"{d}/dedupe.csv"], f"{d}/dedupe.csv")
    sh([PY, f"{HERE}/w_locate_gates.py", "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
        "--dedupe", f"{d}/dedupe.csv", "--survival", f"{d}/{tag}_survival.csv",
        "--out", f"{d}/lg", "--shift", js], f"{d}/lg.json")
    j = json.load(open(f"{d}/lg.json"))
    seeds = W.read_csv_safe(f"{d}/cl_seeds.csv")
    j.update(realisation=k, shift=shift, n_seeds=int(len(seeds)),
             n_survivors_p0=n_surv(f"{d}/p0_survival.csv"),
             n_survivors_final=n_surv(f"{d}/p{W.N_ITERATIONS}_survival.csv"),
             n_pairs=int(len(pd.read_parquet(f"{d}/ac.parquet")))
             if os.path.exists(f"{d}/ac.parquet") else 0)
    json.dump(j, open(f"{d}/summary.json", "w"), indent=1, default=str)
    W.log(f"[null{k}] seeds {j['n_seeds']}, survivors p0 {j['n_survivors_p0']} / final "
          f"{j['n_survivors_final']}, located {j.get('n_located')}, gated {j.get('n_gate_pass')}")
    return j


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--n", type=int, default=W.N_SITESHIFT_REALISATIONS)
    ap.add_argument("--only", type=int, default=-1)
    a = ap.parse_args()
    ks = [a.only] if a.only >= 0 else list(range(1, a.n + 1))
    rows = []
    for k in ks:
        p = f"{a.out_root}/null{k}/summary.json"
        rows.append(json.load(open(p)) if os.path.exists(p) else realisation(k, a.out_root))
    R = pd.DataFrame(rows)
    p = f"{a.out_root}/null_summary.csv"
    if os.path.exists(p):
        R = pd.concat([pd.read_csv(p), R], ignore_index=True).drop_duplicates("realisation",
                                                                              keep="last")
    R = R.sort_values("realisation")
    R.to_csv(p, index=False)
    W.log(f"[null] {len(R)} realisations; gated {list(R.n_gate_pass)}; "
          f"survivors p0 {list(R.n_survivors_p0)}")


if __name__ == "__main__":
    main()
