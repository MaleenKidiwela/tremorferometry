#!/usr/bin/env python
"""FINAL SPEC, STEP 5 (GPU) -- THE ONE BLIND RUN.

Stages 7-10 of 1c' under `results/freeze_step1c_prime_v3.json`:

  7   the frozen stage-7 seed set (every >= 2-node component of pairs >= T_pair = 8.0, top 2,000
      by max pair net_mad) -- built once, hashed into the v3 freeze, and NOT rebuilt here;
  8   the SELECTED order (frozen or staged) with ETS-day-only survival, dedupe, locate, Gates 0-2;
  9   the WHOLE-CHAIN SITE-DAY-SHIFT null x3, re-running the autocorrelation, the clustering and
      the SELECTED order under each shift, reporting GATED FAMILIES **and** SURVIVORS;
  10  `z_score.py` -- the frozen scorer, against oracle_reference_2005_v3.csv.

PASS -> the WA_GO conditions of Merlin's first audit apply with v3 in place of 1c'.
FAIL or FA_BROKEN -> STOP Arm A and WA; the study closes.
results/WA_GO is NEVER created here.
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
import q_common as Q                                                        # noqa: E402
import q_chain as QC                                                        # noqa: E402

PY = sys.executable
STUDY = os.path.dirname(HERE)


def null_realisation(k, run_dir, order, T_json):
    d = f"{run_dir}/null{k}"
    os.makedirs(d, exist_ok=True)
    if os.path.exists(f"{d}/summary.json"):
        return json.load(open(f"{d}/summary.json"))
    shift = W.site_shift(k)
    js = json.dumps(shift)
    json.dump(shift, open(f"{d}/shift.json", "w"), indent=1)
    Q.log(f"[null{k}] site shifts {shift}")
    QC.sh([PY, f"{HERE}/w_autocorr.py", "--hours", "tremor", "--out", f"{d}/ac",
           "--t-record", str(W.T_ENG_FLOOR), "--shift", js], f"{d}/ac.parquet")
    QC.sh([PY, f"{HERE}/w_cluster.py", "cluster", "--pairs", f"{d}/ac.parquet",
           "--thr", T_json, "--out", f"{d}/cl"], f"{d}/cl_seeds.csv")
    r = QC.run(d, f"{d}/cl_seeds.csv", f"{d}/cl_clusters.parquet", order=order, shift_json=js)
    j = Q.jload(f"{d}/lg.json", {}) or {}
    j.update(realisation=k, shift=shift, n_seeds=r["n_seeds"],
             n_survivors_p0=r["n_survivors_p0"], n_survivors_final=r["n_survivors_final"],
             n_pairs=int(len(pd.read_parquet(f"{d}/ac.parquet")))
             if os.path.exists(f"{d}/ac.parquet") else 0)
    Q.jdump(j, f"{d}/summary.json")
    Q.log(f"[null{k}] seeds {j['n_seeds']}, survivors p0 {j['n_survivors_p0']} / final "
          f"{j['n_survivors_final']}, located {j.get('n_located')}, gated {j.get('n_gate_pass')}")
    return j


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=f"{W.RES}/step1c_prime")
    ap.add_argument("--dir", default=f"{Q.OUT}/s5")
    ap.add_argument("--oracle", default=f"{Q.OUT}/oracle_reference_2005_v3.csv")
    ap.add_argument("--out", default=f"{W.RES}/final_score")
    ap.add_argument("--n-null", type=int, default=Q.S5_NULL_REALISATIONS)
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    FR = json.load(open(f"{W.RES}/freeze_step1c_prime_v3.json"))
    order = FR["selected_order"]
    Q.log(f"[s5] THE ONE BLIND RUN, order = {order}, seeds = {FR['seeds']['n']}")

    # ---- stage 7/8: the selected order on the frozen stage-7 seed set --------------------
    r = QC.run(a.dir, f"{a.run}/cl_seeds.csv", f"{a.run}/cl_clusters.parquet", order=order)

    # ---- stage 9: whole-chain site-shift null x3 -----------------------------------------
    rows = [null_realisation(k, a.dir, order, f"{a.run}/T.json")
            for k in range(1, a.n_null + 1)]
    N = pd.DataFrame(rows)
    N.sort_values("realisation").to_csv(f"{a.dir}/null_summary.csv", index=False)
    Q.log(f"[s5/null] gated {list(N.n_gate_pass)}; survivors p0 {list(N.n_survivors_p0)}; "
          f"final {list(N.n_survivors_final)}")

    # ---- stage 10: the frozen scorer -----------------------------------------------------
    rc = subprocess.run([PY, f"{HERE}/z_score.py", "--run", a.dir, "--oracle", a.oracle,
                         "--out", a.out], cwd=STUDY).returncode
    if rc:
        raise SystemExit(f"z_score failed ({rc})")
    S = json.load(open(a.out + ".json"))
    S["order"] = order
    S["chain"] = r
    Q.jdump(S, a.out + ".json")
    Q.log(f"[s5] {S['verdict_line']}")
    Q.marker(f"{Q.OUT}/S5_VERDICT.txt", S["verdict_line"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
