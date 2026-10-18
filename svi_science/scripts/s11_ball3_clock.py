#!/usr/bin/env python
"""s11_ball3_clock -- S11 ballistic round 3, group clock terms and corrected pair tables (round3/PLAN.md §1).
Round-2 quartet logic generalised to co-located groups (SVI: PGC quartet; WA-NW: G1 = B005/B006/B007/LRIV,
G2 = B013/B943): per (family, bin, label) LS of within-group pair differences, median gauge, leave-family-out average
(s11_ball2_clock.aggregate). Writes pair_corrected (dS) and pairP_corrected (dP3, same K from S), clock_terms, and the
co-located floor before/after per bin.

    python s11_ball3_clock.py --region SVI|WA_NW [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L
import s11_ball2_clock as CK
from s11_ball_analysis import scatter


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--region", default="SVI"); ap.add_argument("--tag", default="")
    a = ap.parse_args(); reg, TB, _, _ = L.dirs(a.region)
    pr = L.load(reg, "lags_pair", a.tag)
    D = L.prep_pairs(pr)
    T = L.group_terms(D, reg)
    T, Kt = CK.aggregate(T)
    T.to_parquet(f"{TB}/clock_terms_fam{a.tag}.parquet", index=False); Kt.to_csv(f"{TB}/clock_terms{a.tag}.csv", index=False)
    Dc = L.correct(D, T, reg)
    Dc.to_parquet(f"{TB}/pair_corrected{a.tag}.parquet", index=False)
    keys = ["fam", "pair", "epA", "epB", "bin", "lab"]
    DP = L.prep_pairs(pr, "dP3", "sig_dP3", "railP3").merge(Dc[keys + ["kcorr", "has_corr"]], on=keys, how="left")
    DP["kcorr"] = DP.kcorr.fillna(0.0); DP["has_corr"] = DP.has_corr.fillna(False)
    DP["rc0"] = DP.r + DP.kcorr; DP["rc"] = DP.rc0 - DP.groupby(L.PKEY).rc0.transform("mean")
    DP.to_parquet(f"{TB}/pairP_corrected{a.tag}.parquet", index=False)
    # co-located floor before / after (within-group pairs), per bin
    ng = L.node_group(reg); rows = []
    W = Dc[Dc.A.map(ng).notna() & (Dc.A.map(ng) == Dc.B.map(ng))]
    for (bk, gname), g in W.groupby(["bin", W.A.map(ng)]):
        gc = g[g.has_corr]
        b0 = scatter(gc.assign(r=gc.r)) if len(gc) else (np.nan, np.nan)
        a0 = scatter(gc.assign(r=gc.rc)) if len(gc) else (np.nan, np.nan)
        rows.append(dict(bin=bk, group=gname, n_values=len(g), n_corr=len(gc), sd_before=b0[0], rob_before=b0[1],
                         sd_after=a0[0], rob_after=a0[1]))
        for pair, gp in gc.groupby("pair"):
            b1 = scatter(gp.assign(r=gp.r)); a1 = scatter(gp.assign(r=gp.rc))
            rows.append(dict(bin=bk, group=gname, pair=pair, n_values=len(gp), n_corr=len(gp), sd_before=b1[0],
                             rob_before=b1[1], sd_after=a1[0], rob_after=a1[1]))
    F = pd.DataFrame(rows); F.to_csv(f"{TB}/coloc_floor{a.tag}.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(F[F.pair.isna()].round(2).to_string() if "pair" in F else F.round(2).to_string())
        print(Kt[Kt.glitch & (Kt.bin == "EP")].sort_values("date")[["station", "epoch", "lab", "date", "K", "sigK", "n_fam",
                                                                       "sign_agree"]].round(2).to_string())
    open(f"{L.ROOT3}/CLOCK_{reg}{a.tag}_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
