#!/usr/bin/env python
"""s11_ball3_jointp -- S11 ballistic round 3, §4 joint multi-station P alignment (round3/PLAN.md).
Site-joint S / P3 (equal-weight mean of the members' cc curves on common rows) against single members on the same rows:
median bootstrap sigma_P, P cc, rail rate; best member = the instrument with the lowest median sigma_P in that group and
bin type (headline), per-stack minimum (column). Rule 4. Decade S - P precision (YR) joint vs members.

    python s11_ball3_jointp.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L
import s11_ball3_common as CM


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args(); tag = a.tag
    rng = np.random.default_rng(20261013)
    rows, dec = [], []
    for reg in ("SVI", "WA_NW"):
        p = f"{L.dirs(reg)[1]}/site_joint{tag}.parquet"
        if not os.path.exists(p):
            continue
        J = pd.read_parquet(p)
        for (grp, bk, kind), g in J.groupby(["group", "bin", "kind"]):
            mem = [c[2:-6] for c in g.columns if c.startswith("m_") and c.endswith("_sigP3") and g[c].notna().any()]
            medm = {m: float(g[f"m_{m}_sigP3"].median()) for m in mem}
            best = min(medm, key=medm.get) if medm else None
            ok = ~g.railP3
            r = dict(region=reg, group=grp, bin=bk, kind=kind, n_stacks=len(g), n_mem_med=float(g.n_mem.median()),
                     joint_sigP3=float(g.loc[ok, "sigP3"].median()), joint_ccP3=float(g.ccP3.median()), joint_rail=float(g.railP3.mean()),
                     joint_sigSmP3=float(g.loc[ok, "sigSmP3"].median()),
                     best_member=best, best_member_sigP3=medm.get(best, np.nan),
                     best_member_sigP3_same=float(g.loc[g[f"m_{best}_sigP3"].notna() & ok, f"m_{best}_sigP3"].median()) if best else np.nan,
                     best_member_ccP3=float(g[f"m_{best}_ccP3"].median()) if best else np.nan,
                     best_member_sigSmP3=float(g[f"m_{best}_sigSmP3"].median()) if best else np.nan,
                     perstack_min_sigP3=float(g.loc[ok, "best_sigP3"].median()))
            r["ratio_fixed_best"] = r["joint_sigP3"] / r["best_member_sigP3_same"] if best else np.nan
            r["ratio_perstack_min"] = r["joint_sigP3"] / r["perstack_min_sigP3"]
            for m in mem:
                r[f"med_sigP3_{m}"] = medm[m]
            rows.append(r)
        # decade S - P precision (YR): joint vs each member on the same rows
        d = J[(J.kind == "data") & (J.bin == "YR")].copy()
        for grp, g in d.groupby("group"):
            cand = {"joint": ("lagS", "lagP3", "sigSmP3")}
            for c in g.columns:
                if c.startswith("m_") and c.endswith("_lagSmP3"):
                    cand[c[2:-8]] = (c, None, c.replace("lagSmP3", "sigSmP3"))
            for nm, (c1, c2, sc) in cand.items():
                x = g.copy()
                x["v"] = (x[c1] - x[c2]) if c2 else x[c1]
                if c2:
                    x.loc[x.railS | x.railP3, "v"] = np.nan
                x = x[np.isfinite(x.v) & np.isfinite(x[sc])]
                x["K"] = x.groupby("fam").v.transform("size"); x = x[x.K >= 3]
                if x.fam.nunique() < 3:
                    continue
                x["r"] = x.v - x.groupby("fam").v.transform("mean"); x["sid"] = x.groupby("fam").ngroup()
                x["sig"] = np.maximum(x[sc], 0.5); x["t"] = x.lab.str[:4].astype(int)
                cc = CM.Contr(x.sid, x.r, x.sig ** 2, x.t); C, v = cc.contrast(); okc = np.isfinite(C)
                pr = L.pool_re(C[okc], np.sqrt(v[okc]))
                rob = 1.4826 * np.median(np.abs(x.r))
                dec.append(dict(region=reg, group=grp, estimator=nm, n_series=int(okc.sum()), decade_ms=pr["mean"],
                                decade_sig_ms=pr["sig"], chi2dof=pr["chi2dof"], year_scatter_rob_ms=rob))
    R = pd.DataFrame(rows); Dd = pd.DataFrame(dec)
    R.to_csv(f"{L.ROOT3}/tables/jointp_quality{tag}.csv", index=False); Dd.to_csv(f"{L.ROOT3}/tables/jointp_decade{tag}.csv", index=False)
    v = {}
    D = R[(R.kind == "data") & R.bin.isin(["YR", "EP"])]
    passing = D[D.ratio_fixed_best <= 0.7].group.unique().tolist()
    ctrl = R[(R.kind == "ctrl")].joint_ccP3.median() if (R.kind == "ctrl").any() else np.nan
    v = dict(groups_passing=passing, n_groups=int(D.group.nunique()), ctrl_joint_ccP3_median=float(ctrl),
             PASS=bool(len(passing) >= 2 and (not np.isfinite(ctrl) or ctrl < 0.3)))
    json.dump(v, open(f"{L.ROOT3}/tables/rule4{tag}.json", "w"), indent=1, default=float)
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(R[["region", "group", "bin", "kind", "n_stacks", "joint_sigP3", "best_member", "best_member_sigP3_same", "ratio_fixed_best",
                 "perstack_min_sigP3", "ratio_perstack_min", "joint_ccP3", "best_member_ccP3", "joint_sigSmP3", "best_member_sigSmP3"]].round(3).to_string())
        print(Dd.round(3).to_string()); print(v)


if __name__ == "__main__":
    main()
