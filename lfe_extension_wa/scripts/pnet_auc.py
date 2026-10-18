#!/usr/bin/env python
"""Stage A/B-features: AUC of the PNET candidate features
  (i)  data vs scrambled-control candidates (all 4 controls pooled), and
  (ii) truth-matched (own family, and <= 3 km group) vs unmatched DATA candidates (days the family's array run scanned),
by evaluation set (fit_in_sample = build half V2008+V2010; V2010_test; V2008_test; A2011) x region.
Reference columns: the step-1 features n_z2 (best single feature) and the step-1 P(real).

  python scripts/pnet_auc.py --tag w05   -> results/pnet/stageA_<tag>_auc.csv
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

FE = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_mean", "pnet_rel_max"]
REF = ["n_z2", "P"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    F = pd.concat([pd.read_parquet(p, columns=["date", "fid", "series", "t"] + FE)
                   for p in sorted(glob.glob(f"{PC.OUT}/feat_{a.tag}/feat_*.parquet"))], ignore_index=True)
    C = pd.read_parquet(f"{W.RES}/step1_prob/candidates.parquet",
                        columns=["date", "fid", "series", "t", "region", "eval_set", "truth_ok", "m_own", "m_nbr3"] + REF)
    M = F.merge(C, on=["date", "fid", "series", "t"], how="inner")
    W.log(PC.LOG, f"[auc:{a.tag}] {len(F)} feature rows, {len(C)} candidates, merged {len(M)}")
    rows = []
    for es in ["fit_in_sample", "V2010_test", "V2008_test", "A2011", "ALL_heldout"]:
        ge = M[M.eval_set.isin(["V2010_test", "V2008_test", "A2011"])] if es == "ALL_heldout" else M[M.eval_set == es]
        for rg in ["ALL"] + W.REGIONS:
            g = ge if rg == "ALL" else ge[ge.region == rg]
            d = g[g.series == "data"]; c = g[g.series != "data"]
            t = d[d.truth_ok.astype(bool)]
            for lev in ["ctl", "own", "nbr3"]:
                if lev == "ctl":
                    pos, neg = d, c
                else:
                    mm = t[f"m_{lev}"].astype(bool)
                    pos, neg = t[mm], t[~mm]
                r = dict(tag=a.tag, eval_set=es, region=rg, contrast={"ctl": "data_vs_control", "own": "matched_own_vs_unmatched",
                                                                      "nbr3": "matched_3km_vs_unmatched"}[lev],
                         n_pos=len(pos), n_neg=len(neg))
                for f in FE + REF:
                    r[f] = PC.auc(pos[f], neg[f])
                rows.append(r)
    O = pd.DataFrame(rows)
    O.to_csv(f"{PC.OUT}/stageA_{a.tag}_auc.csv", index=False)
    W.log(PC.LOG, f"[auc:{a.tag}]\n" + O[O.region == "ALL"].round(3).to_string())


if __name__ == "__main__":
    main()
