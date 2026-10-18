#!/usr/bin/env python
"""SH_V3 step A (2026-09-29): re-score detections with Shelly (2025) Text S1 final-catalogue rules.

Variants (cumulative), each at k = 6..12 x daily MAD:
  v2      net_mad >= k only (must reproduce the SH_V2 tables -- a check of this re-implementation)
  +mcc    AND mean cc over live channels >= 0.12
  +floor  AND summed cc >= 2.75
  S1      AND 4 s declustering across all families (greedy by summed cc)   <- his final catalogue
  S1x     as S1 but the 4 s suppression only between DIFFERENT families
Scores (one-to-one +-2 s, t8_score functions, t8_score_aa truth-day filter), days-mode source and fam:
  family-specific  -- our family f vs his family f (the v2 metric)
  family-agnostic  -- all our detections vs all his events pooled, then same_fam = matches whose
                      families agree.
Declustering is applied to the full detection file (all families ran on all test days) BEFORE the
fam-mode day filter.

  python scripts/sh_v3_rescore.py --ours results/sh_v2/mf_it0_test.parquet --tag v2it0 --out results/sh_v3/rescore_v2it0.csv
"""
import os, sys, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v3_common import *
TH = [6, 7, 8, 9, 10, 11, 12]
VAR = {"v2": dict(meancc=False, floor=False, decl=False),
       "+mcc": dict(meancc=True, floor=False, decl=False),
       "+floor": dict(meancc=True, floor=True, decl=False),
       "S1": dict(meancc=True, floor=True, decl=True, same_family=True),
       "S1x": dict(meancc=True, floor=True, decl=True, same_family=False)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--ours", required=True)
    ap.add_argument("--tag", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--t0", default="2018-01-01"); ap.add_argument("--t1", default="2024-12-01")
    ap.add_argument("--variants", default=",".join(VAR))
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    O = pd.read_parquet(a.ours); O["fam"] = O.fam.astype(str)
    t0, t1 = pd.Timestamp(a.t0), pd.Timestamp(a.t1)
    O = O[(O.t >= t0) & (O.t < t1)].reset_index(drop=True)
    S = his_catalogue(); T = S.rename(columns={"meancc": "cc"})[["fam", "t", "cc"]]
    T = T[(T.t >= t0) & (T.t < t1)]
    TM = {m: truth_for(T, m) for m in ("source", "fam")}
    rows = []; fam_rows = []
    for v in a.variants.split(","):
        for k in TH:
            E = apply_rules(O, k, **VAR[v])
            for m in ("source", "fam"):
                Tm, DF = TM[m]
                Om = ours_for(E, m, DF).reset_index(drop=True)
                F, g = score_all(Om, Tm.reset_index(drop=True))
                rows.append(dict(tag=a.tag, variant=v, xmad=k, mode=m, **g,
                                 recall_fam=g["match_fam"] / max(g["his"], 1), prec_fam=g["match_fam"] / max(g["ours"], 1),
                                 recall_any=g["match_any"] / max(g["his"], 1), prec_any=g["match_any"] / max(g["ours"], 1),
                                 fam_ge1=int((F.n_match > 0).sum()), fam_p30=int((F.precision >= 0.3).sum())))
                F["tag"] = a.tag; F["variant"] = v; F["xmad"] = k; F["mode"] = m; fam_rows.append(F)
                r = rows[-1]
                log(f"[{a.tag} {v} {k}x {m}] ours {r['ours']} his {r['his']} | fam R {r['recall_fam']:.3f} P {r['prec_fam']:.3f}"
                    f" | any R {r['recall_any']:.3f} P {r['prec_any']:.3f} same_fam {r['same_fam']}")
    pd.DataFrame(rows).to_csv(a.out, index=False)
    pd.concat(fam_rows).to_csv(a.out.replace(".csv", "_fam.csv"), index=False)
if __name__ == "__main__": main()
