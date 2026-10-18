#!/usr/bin/env python
"""ARM B stage 4 -- space-time clustering of the associated events into family candidates.

DBSCAN on the EPICENTRES (local equirectangular km projection), eps = 5 km, min_samples = 20,
over the WHOLE window (time is not a clustering coordinate: an LFE family repeats through the
episode, so a time radius would split it).  Clusters with fewer than 20 events are dropped.

The output is written in the EXACT format `w_grow.py handoff` consumes, so stages 5-6 are the
frozen Arm-A engine with no modification:
    <out>_clusters.parquet   cid, t          (cluster member ORIGIN times, the growth axis)
    <out>_seeds.csv          cid, date, n_members, n_in_hour, n_in_3h, n_bins, t0, t1, t_seed
plus <out>_centroids.csv with the PyOcto centroid of each cluster (diagnostic only -- the
location used for scoring comes from the frozen `w_locate_gates` stage, not from here).

Usage
  python scripts/b_armB_cluster.py --events results/armB/events_t30.csv --out results/armB/cl
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
import w_common as W                                                        # noqa: E402

EPS_KM = 5.0
MIN_EVENTS = 20
SEED_COLS = ["cid", "date", "n_members", "n_in_hour", "n_in_3h", "n_bins", "t0", "t1", "t_seed",
             "max_net_mad", "n_pairs"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--eps-km", type=float, default=EPS_KM)
    ap.add_argument("--min-events", type=int, default=MIN_EVENTS)
    a = ap.parse_args()
    from sklearn.cluster import DBSCAN

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    E = pd.read_csv(a.events) if os.path.exists(a.events) else pd.DataFrame()
    rows, mem, cen = [], [], []
    if len(E):
        lat0 = float(E.latitude.median())
        x = (E.longitude.to_numpy(float) - float(E.longitude.median())) * 111.32 * np.cos(
            np.deg2rad(lat0))
        y = (E.latitude.to_numpy(float) - lat0) * 110.57
        lab = DBSCAN(eps=a.eps_km, min_samples=a.min_events).fit(np.c_[x, y]).labels_
        E = E.assign(label=lab)
        cid = 0
        for k, g in E[E.label >= 0].groupby("label"):
            if len(g) < a.min_events:
                continue
            t = np.sort(g.time.to_numpy(float))
            cid += 1
            c = f"B{cid:05d}"
            hh = np.floor(t / 3600.0).astype(int)
            n_hour = max(int((hh == h).sum()) for h in np.unique(hh))
            n_3h = max(int(((hh >= h - 1) & (hh <= h + 1)).sum()) for h in np.unique(hh))
            nb = len(np.unique(np.floor(t / W.SEED_BIN_S).astype(int)))
            ts = float(np.median(t))
            rows.append(dict(cid=c, date=str(pd.Timestamp(ts, unit="s").date()),
                             n_members=int(len(t)), n_in_hour=n_hour, n_in_3h=n_3h, n_bins=nb,
                             t0=float(t.min()), t1=float(t.max()), t_seed=ts,
                             max_net_mad=float(g.mean_prob.max()) if "mean_prob" in g else np.nan,
                             n_pairs=int(len(t))))
            mem.append(pd.DataFrame(dict(cid=c, t=t)))
            cen.append(dict(cid=c, n_events=int(len(t)), lat=float(g.latitude.median()),
                            lon=float(g.longitude.median()), dep=float(g.depth.median()),
                            n_days=int(pd.to_datetime(g.time, unit="s").dt.date.nunique()),
                            rms_resid=float(g.rms_resid.median())
                            if "rms_resid" in g else np.nan))
    Cl = pd.DataFrame(rows, columns=SEED_COLS)
    if len(Cl):
        Cl = Cl.sort_values("n_members", ascending=False).head(W.SEED_BUDGET).reset_index(drop=True)
        keep = set(Cl.cid)
        M = pd.concat([m for m in mem if m.cid.iloc[0] in keep], ignore_index=True)
        C = pd.DataFrame(cen)
        C = C[C.cid.isin(keep)]
    else:
        M = pd.DataFrame(columns=["cid", "t"])
        C = pd.DataFrame(columns=["cid", "n_events", "lat", "lon", "dep", "n_days", "rms_resid"])
    Cl.to_csv(a.out + "_seeds.csv", index=False)
    M.to_parquet(a.out + "_clusters.parquet", index=False)
    C.to_csv(a.out + "_centroids.csv", index=False)
    json.dump(dict(n_events=int(len(E)), n_clusters=int(len(Cl)), eps_km=a.eps_km,
                   min_events=a.min_events, seed_budget=W.SEED_BUDGET,
                   n_noise=int((E.label < 0).sum()) if len(E) and "label" in E else 0,
                   median_cluster_size=float(Cl.n_members.median()) if len(Cl) else None),
              open(a.out + "_meta.json", "w"), indent=1)
    W.log(f"[cluster] {len(E)} events -> {len(Cl)} clusters >= {a.min_events} events "
          f"(eps {a.eps_km} km); median size "
          f"{float(Cl.n_members.median()) if len(Cl) else 0}")


if __name__ == "__main__":
    main()
