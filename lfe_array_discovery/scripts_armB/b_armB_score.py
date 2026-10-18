#!/usr/bin/env python
"""ARM B stage 6 -- scoring against the 10 isolated hidden Bostock families + the frozen bars.

Runs the frozen scorer `w_score_step1.py` unchanged (in process) and then applies Arm B's OWN
pre-registered ABSOLUTE bars from `results/freeze_armB.json`, which do not depend on any
reference detector:

    a hidden ISOLATED family is REDISCOVERED_ABS when a cluster lies <= 5 km from it AND its
    one-to-one recall at +- 2 s is >= 0.30 AND its purity under the 3-km cluster truth is >= 0.50.

The RELATIVE clauses of `w_score_step1` need a reference.  Merlin rejected the transferred 2012
template reference for 2005 (F4a).  If `results/step1c_prime/oracle_reference_2005.csv` exists
(agent H's Step-6 oracle ceiling) it is used; otherwise the run is marked
`relative_reference = "pending"` and ONLY the absolute bars are read as the Arm-B verdict.

Usage
  python scripts/b_armB_score.py --clusters results/armB/lg_clusters.parquet \
      --gates results/armB/lg_gates.csv --nulls results/armB/null_summary.csv \
      --out results/armB_score
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("W_BED", "c")
ORACLE = "/home/jovyan/tremorferometry/lfe_array_discovery/results/step1c_prime/" \
         "oracle_reference_2005.csv"


def resolve_reference(work):
    """(path, note).  The 2005 oracle ceiling if agent H has written it, else the 2012 transfer."""
    if os.path.exists(ORACLE):
        try:
            R = pd.read_csv(ORACLE)
        except Exception as e:                                              # noqa: BLE001
            return None, f"oracle unreadable ({type(e).__name__}); relative reference pending"
        ren = {}
        if "recall_ref" not in R.columns:
            for a, b in (("recall", "recall_ref"), ("purity", "purity_ref"),
                         ("C_iso_recall", "recall_ref"), ("C_iso_purity", "purity_ref")):
                if a in R.columns and b not in R.columns:
                    ren[a] = b
        R = R.rename(columns=ren)
        key = "bkey" if "bkey" in R.columns else ("hidden" if "hidden" in R.columns else None)
        if key is None or "recall_ref" not in R.columns or "purity_ref" not in R.columns:
            return None, ("oracle_reference_2005.csv present but its columns are not "
                          f"({list(R.columns)[:8]}); relative reference pending")
        p = os.path.join(work, "oracle_reference_2005_asref.csv")
        R.rename(columns={key: "bkey"}).to_csv(p, index=False)
        return p, f"2005 oracle ceiling ({ORACLE})"
    return None, "relative reference PENDING (results/step1c_prime/oracle_reference_2005.csv " \
                 "does not exist yet); absolute bars govern the Arm-B verdict"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--gates", default="")
    ap.add_argument("--nulls", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--freeze", default="results/freeze_armB.json")
    a = ap.parse_args()
    import w_common as W
    import w_score_step1 as SC

    work = os.path.dirname(a.out) or "."
    os.makedirs(work, exist_ok=True)
    ref, note = resolve_reference(work)
    if ref:
        SC.PICKREF = ref
    argv = [sys.argv[0], "--clusters", a.clusters, "--out", a.out, "--year", "2005",
            "--windows", f"{W.bed_cfg()['day0']},{W.bed_cfg()['day1']}"]
    if a.gates:
        argv += ["--gates", a.gates]
    if a.nulls and os.path.exists(a.nulls):
        argv += ["--nulls", a.nulls]
    old = sys.argv
    sys.argv = argv
    try:
        SC.main()
    finally:
        sys.argv = old

    F = json.load(open(a.freeze))
    bars = F["accept"]
    D = pd.read_csv(a.out + "_isolated_detail.csv")
    ok = ((D.dist_km <= bars["centroid_km"]) & (D.recall >= bars["recall_abs"])
          & (D.purity_cluster3km >= bars["purity_abs_cluster3km"]))
    D = D.assign(rediscovered_abs=ok.fillna(False))
    D.to_csv(a.out + "_isolated_detail.csv", index=False)
    n_abs = int(D.rediscovered_abs.sum())
    S = json.load(open(a.out + ".json"))
    nulls = S.get("null_gated_per_realisation") or []
    nullc = []
    for p in sorted(__import__("glob").glob(os.path.join(work, "null*/cl_meta.json"))):
        try:
            nullc.append(int(json.load(open(p)).get("n_clusters", -1)))
        except Exception:                                                   # noqa: BLE001
            pass
    verdict = ("PASS" if (n_abs >= bars["n_isolated"]
                          and (not nullc or max(nullc) <= bars["null_clusters_le"]))
               else ("FA_RULE_BROKEN" if nullc and max(nullc) > bars["null_clusters_le"]
                     and n_abs >= bars["n_isolated"]
                     else ("ARMB_INADEQUATE" if n_abs < 2 else "BELOW_BAR")))
    out = dict(arm="B (LFEDetect jms -> PyOcto -> DBSCAN -> frozen growth)",
               relative_reference=note,
               n_isolated_rediscovered_abs=n_abs,
               n_isolated_rediscovered_relative=S.get("n_isolated_rediscovered"),
               per_family=D.to_dict("records"),
               median_centroid_km=S.get("median_centroid_km"),
               median_recall=S.get("median_recall"),
               median_purity_cluster3km=S.get("median_purity_cluster3km"),
               spurious_split_rate=S.get("spurious_split_rate"),
               null_gated_per_realisation=nulls, null_clusters_per_realisation=nullc,
               bars=bars, verdict=verdict,
               stop_rule=("Merlin STOP: Arm B(jms) north < 2/10 -> STOP the study "
                          "('array-era discovery not reproducible from the public archive at "
                          "this SNR with these methods')"))
    json.dump(out, open(a.out + "_armB.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "per_family"}, indent=1, default=str))


if __name__ == "__main__":
    main()
