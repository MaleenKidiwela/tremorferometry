#!/usr/bin/env python
"""WA FIXED-10 re-run: station liveness + 'strict' unit-days (2026-10-07).  Reads OUT/mf/prod (run_*.json, famday_*.parquet),
OUT/setup/fixed_sets.json and the post_v2 rate table of one network; writes only under that network's post directory.

liveness   per station x era: share of production days on which the station was loaded and passed the guards (mf run
           members), + every absence run >= 30 days  -> <post>/station_liveness_by_era.csv, <post>/station_outages_ge30d.csv
strict     a family-day is STRICT when every station of the family's fixed template set contributes that day:
             core (F9):  the family's PB-borehole template stations (one sensor + one frozen template set each, 2007-2026)
             all  (F10): the same + UW.GNW whenever GNW's sensor epoch of that date has a template for the family
                         (BH.50 epoch to 2010-11-30: 260 families; BH.40 / HH.100 from 2010-12-01: 52 families)
           Olympic A group-day strict = every observable member family strict.  -> <post>/unit_day_strict.parquet
           Rates are reported BOTH ways: all scanned unit-days (post_v2 tables as they are) and strict unit-days only
           (<post>/strict/rates_summary_v2.csv, same function and columns as post_v2's summary).

  WA_PROD_DIR=... WA_POST_NET=core WA_POST_V2_DIR=... WA_POST_HEAD_LAB=own python scripts/wa_f10_strict.py
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
import wa_prod_post_v2 as V2

NET = V2.NET
O2 = V2.O2


def liveness(FX):
    rows = []
    for p in sorted(glob.glob(f"{P.OUT}/mf/prod/run_*.json")):
        J = json.load(open(p))
        mem = set(J["members"].get("all", []))
        for s in FX["F10"]:
            rows.append((J["date"], J["era"], s, s in mem))
    D = pd.DataFrame(rows, columns=["date", "era", "netsta", "live"])
    T = D.groupby(["netsta", "era"]).live.agg(["mean", "sum", "size"]).reset_index()
    T.columns = ["netsta", "era", "live_frac", "live_days", "days"]
    T.to_csv(f"{O2}/station_liveness_by_era.csv", index=False)
    out = []
    for s, g in D.sort_values("date").groupby("netsta"):
        v = g.live.to_numpy(); d = g.date.to_numpy(); i = 0
        while i < len(v):
            if not v[i]:
                j = i
                while j < len(v) and not v[j]:
                    j += 1
                if j - i >= 30:
                    out.append((s, d[i], d[j - 1], j - i, P.era_of(d[i])))
                i = j
            else:
                i += 1
    pd.DataFrame(out, columns=["netsta", "first_day", "last_day", "n_days", "era_start"]).to_csv(
        f"{O2}/station_outages_ge30d.csv", index=False)
    V2.log(f"[strict] liveness:\n{T.pivot(index='netsta', columns='era', values='live_frac').round(3).to_string()}")


def strict_days(FX):
    fam = {int(k): v for k, v in FX["families"].items()}
    gnw_bh50 = {f for f, v in fam.items() if any("BH.50" in n for n in v["GNW_nodes"])}
    gnw_bh40 = {f for f, v in fam.items() if any("BH.40" in n for n in v["GNW_nodes"])}
    rows = []
    for p in sorted(glob.glob(f"{P.OUT}/mf/prod/famday_*.parquet")):
        X = pd.read_parquet(p, columns=["net", "fid", "date", "stations", "observable"])
        X = X[(X.net == NET) & X.observable.astype(bool)]
        if not len(X):
            continue
        d = str(X.date.iloc[0])
        gset = gnw_bh50 if d < "2010-12-01" else gnw_bh40
        for f, st in zip(X.fid.astype(int), X.stations.astype(str)):
            have = set(st.split(",")) if st else set()
            need = set(fam.get(f, {}).get("F9", []))
            if NET == "all" and f in gset:
                need.add("UW.GNW")
            rows.append((f, d, len(need), len(need & have), bool(need) and need <= have))
    S = pd.DataFrame(rows, columns=["fid", "date", "n_fixed", "n_fixed_live", "strict"])
    F = P.families()
    grp = dict(zip(F.fid, F.pos_group.fillna("")))
    S["pos_group"] = S.fid.map(grp)
    fam_u = S[S.pos_group == ""].assign(unit_kind="family", unit=lambda x: x.fid.astype(str))
    G = S[S.pos_group != ""].groupby(["pos_group", "date"]).agg(n_fixed=("n_fixed", "sum"), n_fixed_live=("n_fixed_live", "sum"),
                                                                 strict=("strict", "all")).reset_index()
    G = G.rename(columns={"pos_group": "unit"}).assign(unit_kind="olyA_group")
    U = pd.concat([fam_u[["unit_kind", "unit", "date", "n_fixed", "n_fixed_live", "strict"]],
                   G[["unit_kind", "unit", "date", "n_fixed", "n_fixed_live", "strict"]]], ignore_index=True)
    U.to_parquet(f"{O2}/unit_day_strict.parquet", index=False)
    return U


def main():
    FX = json.load(open(f"{P.OUT}/setup/fixed_sets.json"))
    liveness(FX)
    U = strict_days(FX)
    R = pd.read_parquet(f"{O2}/rates_daily_v2.parquet")
    R = R.merge(U[["unit_kind", "unit", "date", "strict"]], on=["unit_kind", "unit", "date"], how="left")
    R["strict"] = R.strict.fillna(False).astype(bool)
    sh = R.groupby(["region", "era"]).strict.mean().round(3)
    V2.log(f"[strict] net {NET}: strict share of scanned unit-days by region x era:\n{sh.unstack().to_string()}")
    V = pd.concat([V2.view_rows(R[R.strict], n) for n, _, _ in V2.VIEWS], ignore_index=True)
    os.makedirs(f"{O2}/strict", exist_ok=True)
    V2.O2 = f"{O2}/strict"                       # summary() writes to the module's O2
    V2.summary(V)
    V2.O2 = O2
    sh.unstack().to_csv(f"{O2}/strict_share_by_region_era.csv")
    V2.log(f"[strict] done -> {O2}/unit_day_strict.parquet, {O2}/strict/rates_summary_v2.csv")


if __name__ == "__main__":
    main()
