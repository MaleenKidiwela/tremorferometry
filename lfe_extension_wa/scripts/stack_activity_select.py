#!/usr/bin/env python
"""Pick the families + stations of the continuous 2012-06 -> 2014-12 activity series, and write its setup.

Rule (fixed before the validation results were seen; a compute/download budget, not a gate on any result):
  families  per-family AUC of A (data) for 'active day' (>= 5 array detections) at W = 1, pooled over the
            validation runs present (valA, valB, valC); families with >= 3 active and >= 3 inactive scanned days;
            the TOP 20 by AUC.
  stations  for those families, each station scores sum_f max(0, mean z_data - mean z_control1 on the family's
            active W = 1 days); the TOP 12 stations, at most 2 per 3-km site.
  -> results/stack_activity/cont_selection.json, cont_fids.csv, cont_stations.csv,
     results/stack_activity/setup_cont/{station_days,stations,families}.csv  (WA_SETUP for wa_fetch / wa_cache)
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
from stack_activity import OUT
from stack_activity_report import auc

N_FAM, N_STA, N_ACTIVE = 20, 12, 5
D0, D1 = "2012-06-01", "2014-12-31"
MB_PER_DAY = 12.3        # trial fetch: 26.3 GB / 2134 station-days


def main():
    FWs, SWs = [], []
    for r in ["valA", "valB", "valC"]:
        p = f"{OUT}/score_{r}/family_window.parquet"
        if os.path.exists(p):
            x = pd.read_parquet(p)
            x = x[(x.W == 1) & x.n_truth.notna() & (x.n_scanned_days > 0)].assign(run=r)
            FWs.append(x)
            s = pd.read_parquet(f"{OUT}/score_{r}/station_window.parquet",
                                columns=["W", "label", "block", "fid", "netsta", "z_data", "z_control1"])
            SWs.append(s[s.W == 1].assign(run=r))
    FW = pd.concat(FWs, ignore_index=True)
    SW = pd.concat(SWs, ignore_index=True)
    FW["active"] = FW.n_truth >= N_ACTIVE
    rows = []
    for f, g in FW.groupby("fid"):
        na = int(g.active.sum()); ni = int((~g.active).sum())
        rows.append(dict(fid=f, n_active=na, n_inactive=ni,
                         auc=auc(g.A_data.to_numpy(), g.active.to_numpy().astype(float)) if na >= 3 and ni >= 3 else np.nan))
    R = pd.DataFrame(rows).dropna(subset=["auc"]).sort_values("auc", ascending=False)
    top = R.head(N_FAM)
    act = FW[FW.fid.isin(top.fid) & FW.active][["run", "label", "block", "fid"]]
    s = SW.merge(act, on=["run", "label", "block", "fid"])
    st = s.groupby(["fid", "netsta"])[["z_data", "z_control1"]].mean()
    st["gain"] = (st.z_data - st.z_control1).clip(lower=0)
    sta = st.groupby("netsta").gain.sum().sort_values(ascending=False)
    site = dict(zip(pd.read_csv(f"{W.RES}/step4b/setup/stations.csv").netsta, pd.read_csv(f"{W.RES}/step4b/setup/stations.csv").site))
    stations, per_site = [], {}
    for ns_ in sta.index:                       # at most 2 stations per 3-km site (site diversity for the network)
        q = site.get(ns_, ns_)
        if per_site.get(q, 0) < 2 and len(stations) < N_STA:
            stations.append(ns_); per_site[q] = per_site.get(q, 0) + 1
    sta = sta[stations]
    # ---------- setup for fetch / cache
    sd = f"{OUT}/setup_cont"
    os.makedirs(sd, exist_ok=True)
    S4 = pd.read_csv(f"{W.RES}/step4b/setup/stations.csv")
    S4.to_csv(f"{sd}/stations.csv", index=False)
    pd.read_csv(f"{W.RES}/setup/families.csv").to_csv(f"{sd}/families.csv", index=False)
    N = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]])
    X = pd.read_csv(f"{W.RES}/step4b/setup/extra_nodes.csv")
    N = pd.concat([N[~(N.net + "." + N.sta).isin(set(X.net + "." + X.sta))], X]).drop_duplicates("node")
    out = []
    for d in pd.date_range(D0, D1, freq="D"):
        ds = d.strftime("%Y-%m-%d")
        for ns in stations:
            net, sta_ = ns.split(".")
            n = N[(N.net == net) & (N.sta == sta_) & (N.start <= ds) & (N.end > ds)]
            out.append(dict(window="CONT", date=ds, kind="cont", netsta=ns, meta_live=len(n) > 0,
                            meta_nodes=";".join(n.node), raw=W.raw_path(ns, ds) or ""))
    SD = pd.DataFrame(out)
    SD.to_csv(f"{sd}/station_days.csv", index=False)
    live = SD[SD.meta_live]
    nf = int((live.raw == "").sum())
    F = pd.read_csv(f"{W.RES}/setup/families.csv")
    uid = dict(zip(F.fid, F.family_uid))
    top.to_csv(f"{OUT}/cont_fids.csv", index=False)
    pd.DataFrame(dict(netsta=stations, gain=sta.head(N_STA).values)).to_csv(f"{OUT}/cont_stations.csv", index=False)
    J = dict(rule=__doc__, fids=[int(x) for x in top.fid], family_uids=[uid[int(x)] for x in top.fid],
             auc=[float(x) for x in top.auc], stations=stations, station_gain=[float(x) for x in sta.head(N_STA)],
             n_station_days=int(len(SD)), n_meta_live=int(len(live)), n_on_disk=int((live.raw != "").sum()),
             n_to_fetch=nf, est_GB=round(nf * MB_PER_DAY / 1000, 1), date_range=f"{D0}:{D1}")
    json.dump(J, open(f"{OUT}/cont_selection.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in J.items() if k != "rule"}, indent=1))


if __name__ == "__main__":
    main()
