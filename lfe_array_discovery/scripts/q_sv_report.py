#!/usr/bin/env python
"""Build results/SCREEN_VALIDATION_NORTH.md from the measured tables.  Prefix `q_sv_`."""
from __future__ import annotations
import json, os
import numpy as np, pandas as pd

OUT = "results/q_screenval"
CLAUSES = [("c_a", "(a) quiet-day FA <= 0.1/tpl-day AND tremor rate >= 3x it"),
           ("c_b", "(b) tremor co-timing > same-day base rate"),
           ("c_c1", "(c1) located on >= 5 stations"),
           ("c_c2", "(c2) shuffled-station control passes"),
           ("c_c3", "(c3) INTERIOR to the wide grid"),
           ("c_c4", "(c4) |dz| <= 8 km of Slab2"),
           ("c_d", "(d) >= 3 distinct tremor days"),
           ("c_e", "(e) 2nd-half re-stack relocates within 10 km")]


def famfrac(S, key, col):
    f = S.groupby(key)[col].max()
    return int(f.sum()), len(f)


def main():
    V = pd.read_csv(f"{OUT}/screen_north.csv")
    for c, _ in CLAUSES:
        if c not in V.columns:
            V[c] = False
        V[c] = V[c].fillna(False).astype(bool)
    K = V[V.cls == "kept_family"]; H = V[V.cls == "hidden_family"]
    N = V[V.cls == "unlocatable"]; CN = V[V.cls == "candidate_new"]
    rows = []
    for c, nm in CLAUSES:
        kf, kn = famfrac(K, "bkey", c)
        hf, hn = famfrac(H, "hid_bkey", c)
        rows.append(dict(clause=c, name=nm,
                         kept_fam=f"{kf}/{kn} = {kf/kn:.0%}",
                         kept_fam_frac=kf / kn,
                         kept_tpl=f"{int(K[c].sum())}/{len(K)} = {K[c].mean():.0%}",
                         hid_fam=f"{hf}/{hn} = {hf/hn:.0%}",
                         hid_tpl=f"{int(H[c].sum())}/{len(H)} = {H[c].mean():.0%}",
                         cand=f"{int(CN[c].sum())}/{len(CN)}",
                         noise_pass=f"{int(N[c].sum())}/{len(N)} = {N[c].mean():.0%}",
                         noise_reject=1 - N[c].mean()))
    T = pd.DataFrame(rows)
    T.to_csv(f"{OUT}/clause_table.csv", index=False)
    print(T[["clause", "kept_fam", "kept_tpl", "hid_fam", "noise_pass"]].to_string(index=False))
    # verdicts
    def verdict(r):
        if r.kept_fam_frac < 0.85:
            return "UNJUSTIFIED as coded"
        return "JUSTIFIED" if r.noise_reject >= 0.50 else "weak (passes true families, little power)"
    T["verdict"] = [verdict(r) for r in T.itertuples()]
    print(T[["clause", "verdict"]].to_string(index=False))
    # cumulative
    cum = pd.Series(True, index=V.index)
    cr = []
    for c, nm in [x for x in CLAUSES]:
        cum = cum & V[c]
        kf = K.assign(x=cum.loc[K.index].values).groupby("bkey").x.max().sum()
        hf = H.assign(x=cum.loc[H.index].values).groupby("hid_bkey").x.max().sum()
        cr.append(dict(clause=c, kept_tpl=int(cum.loc[K.index].sum()), kept_fam=int(kf),
                       hid_tpl=int(cum.loc[H.index].sum()), hid_fam=int(hf),
                       cand=int(cum.loc[CN.index].sum()), noise=int(cum.loc[N.index].sum())))
    CU = pd.DataFrame(cr)
    CU.to_csv(f"{OUT}/cumulative_table.csv", index=False)
    print(CU.to_string(index=False))
    V.to_csv(f"{OUT}/screen_north.csv", index=False)


if __name__ == "__main__":
    main()
