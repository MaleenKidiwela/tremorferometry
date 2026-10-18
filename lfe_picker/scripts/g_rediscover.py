#!/usr/bin/env python
"""P3 REDISCOVERY SCORER — did hidden families re-emerge as clusters, and were kept families split?

Merlin's criterion (plan review, P3): a hidden family counts as REDISCOVERED when some proposed
cluster
  (a) has a centroid within 5 km of it,
  (b) matches >= 30 % of that family's Bostock events ONE-TO-ONE (+-2 s) inside the evaluated windows,
  (c) has member purity >= 50 % against that family, and
  (d) is not > 50 % explained by a KEPT family's events.
Also reported: the SPURIOUS-SPLIT rate for kept families (one true family appearing as >= 2 clusters).

Detector-agnostic: clusters are (cluster id, member times, centroid lat/lon); families are
(id, lat, lon, event times).  Nothing here knows how the clusters were made.

    python scripts/g_rediscover.py --clusters X.parquet --events Y.parquet --families F.csv \
           --hidden 8,10,12 --windows 2012-08-30,2012-10-05
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
import g_common as G                                             # noqa: E402


def one2one(a, b, tol=2.0):
    """Greedy-by-time one-to-one matching of times a against sorted times b.
    Returns a boolean mask over a (matched) and the index in b it took (or -1).
    Same recipe as lfe_extension d1_score_epochs.one2one / a1_hf_veto.one2one."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    mask = np.zeros(len(a), bool)
    take = np.full(len(a), -1, int)
    if len(a) == 0 or len(b) == 0:
        return mask, take
    order = np.argsort(a, kind="stable")
    used = np.zeros(len(b), bool)
    j = 0
    for i in order:
        ti = a[i]
        while j < len(b) and b[j] < ti - tol:
            j += 1
        k = j
        while k < len(b) and b[k] <= ti + tol:
            if not used[k]:
                used[k] = True
                mask[i] = True
                take[i] = k
                break
            k += 1
    return mask, take


def in_windows(t, windows):
    t = np.asarray(t, float)
    if not windows:
        return np.ones(len(t), bool)
    m = np.zeros(len(t), bool)
    for a, b in windows:
        m |= (t >= a) & (t <= b)
    return m


def score(clusters, events, families, hidden, windows=None, tol=2.0, r_km=5.0,
          min_match=0.30, min_purity=0.50, max_explained=0.50, split_frac=0.20):
    """clusters: dict cid -> dict(times=array, lat=, lon=)
       events:   dict fid -> array of event times
       families: dict fid -> (lat, lon)
       hidden:   iterable of fid hidden from the detector's labels"""
    hidden = set(hidden)
    kept = [f for f in families if f not in hidden]
    ev = {f: np.sort(np.asarray(t, float)[in_windows(t, windows)]) for f, t in events.items()}
    cl = {c: dict(times=np.sort(np.asarray(v["times"], float)[in_windows(v["times"], windows)]),
                  lat=v["lat"], lon=v["lon"]) for c, v in clusters.items()}
    # cluster -> per-family matched counts (one-to-one, the CLUSTER side is the thing matched)
    expl = {}
    for c, v in cl.items():
        e = {}
        for f, t in ev.items():
            m, _ = one2one(v["times"], t, tol)
            e[f] = int(m.sum())
        expl[c] = e
    rows = []
    for h in sorted(hidden):
        if h not in families:
            continue
        best = None
        for c, v in cl.items():
            n = len(v["times"])
            if n == 0:
                continue
            d = float(G.epi_km(v["lat"], v["lon"], *families[h]))
            mh, _ = one2one(ev[h], v["times"], tol)             # family side: recall
            recall = float(mh.mean()) if len(ev[h]) else 0.0
            purity = expl[c][h] / n
            kbest = max((expl[c][k] / n for k in kept), default=0.0)
            kwho = max(kept, key=lambda k: expl[c][k], default=None) if kept else None
            ok = (d <= r_km and recall >= min_match and purity >= min_purity
                  and kbest <= max_explained)
            r = dict(hidden=h, cluster=c, dist_km=d, n_cluster=n, n_family=len(ev[h]),
                     recall=recall, purity=purity, kept_explained=kbest, kept_family=kwho,
                     rediscovered=bool(ok))
            if best is None or (ok, recall) > (best["rediscovered"], best["recall"]):
                best = r
        rows.append(best or dict(hidden=h, cluster=None, rediscovered=False))
    R = pd.DataFrame(rows)
    # ---- spurious splits of KEPT families
    srows = []
    for k in kept:
        n = len(ev[k])
        if n == 0:
            continue
        parts = []
        for c, v in cl.items():
            nc = len(v["times"])
            if nc == 0:
                continue
            mk, _ = one2one(ev[k], v["times"], tol)
            frac = float(mk.mean())
            pur = expl[c][k] / nc
            if frac >= split_frac and pur >= min_purity:
                parts.append((c, frac, pur))
        srows.append(dict(family=k, n_events=n, n_parts=len(parts),
                          split=bool(len(parts) >= 2),
                          parts=",".join(str(p[0]) for p in parts)))
    S = pd.DataFrame(srows)
    summ = dict(
        n_hidden=int(len(hidden)), n_kept=int(len(kept)), n_clusters=int(len(cl)),
        n_rediscovered=int(R.rediscovered.sum()) if len(R) else 0,
        rediscovery_rate=float(R.rediscovered.mean()) if len(R) else np.nan,
        n_kept_split=int(S.split.sum()) if len(S) else 0,
        spurious_split_rate=float(S.split.mean()) if len(S) else np.nan,
        criterion=dict(r_km=r_km, tol_s=tol, min_match=min_match, min_purity=min_purity,
                       max_explained=max_explained, split_frac=split_frac))
    return R, S, summ


def _demo():
    """Self-test on 2012 truth: 'clusters' = the frozen template detections of each family (the best
    a detector could hand P3), families = Bostock, hidden = the 30 holdout families.  Scores what the
    EXISTING detector would score, so the scorer is calibrated before any learned model uses it."""
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv").set_index("bkey")
    gt = pd.read_parquet(f"{G.REP}/results/e5r_gt_bostock2012.parquet")
    det = pd.read_parquet(f"{G.EXT}/results/e5r_detections_a1_mergeB_L12.parquet")
    det = det[(det.L == 12.0) & (det.net_mad >= 8.25)]
    hid = json.load(open(f"{G.RESULTS}/holdouts.json"))
    hidden = [int(h["bkey"]) for h in hid["hidden_families"]]
    clusters, events, families = {}, {}, {}
    for b in F.index:
        families[int(b)] = (float(F.loc[b, "lat"]), float(F.loc[b, "lon"]))
        events[int(b)] = gt.t_origin.values[gt.bkey.values == b]
    for b, g in det.groupby("fam"):
        b = int(b)
        if b not in families or len(g) < 20:
            continue
        clusters[b] = dict(times=g.t.to_numpy(float), lat=families[b][0], lon=families[b][1])
    R, S, summ = score(clusters, events, families, hidden)
    os.makedirs(G.GATES, exist_ok=True)
    R.to_csv(f"{G.GATES}/g_rediscovery_demo.csv", index=False)
    S.to_csv(f"{G.GATES}/g_rediscovery_demo_splits.csv", index=False)
    json.dump(summ, open(f"{G.GATES}/g_rediscovery_demo.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))
    print(R.to_string(index=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--clusters", default="")
    ap.add_argument("--events", default="")
    ap.add_argument("--families", default="")
    ap.add_argument("--hidden", default="")
    ap.add_argument("--windows", default="")
    ap.add_argument("--out", default=f"{G.GATES}/g_rediscovery")
    a = ap.parse_args()
    if a.demo or not a.clusters:
        return _demo()
    C = pd.read_parquet(a.clusters) if a.clusters.endswith("parquet") else pd.read_csv(a.clusters)
    E = pd.read_parquet(a.events) if a.events.endswith("parquet") else pd.read_csv(a.events)
    F = pd.read_csv(a.families).set_index("bkey")
    clusters = {c: dict(times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]), lon=float(g.lon.iloc[0]))
                for c, g in C.groupby("cid")}
    events = {int(b): g.t.to_numpy(float) for b, g in E.groupby("bkey")}
    families = {int(b): (float(F.loc[b, "lat"]), float(F.loc[b, "lon"])) for b in F.index}
    hidden = [int(x) for x in a.hidden.split(",") if x]
    windows = []
    for w in a.windows.split(";"):
        if w:
            x, y = w.split(",")
            windows.append((pd.Timestamp(x).timestamp(), pd.Timestamp(y).timestamp()))
    R, S, summ = score(clusters, events, families, hidden, windows or None)
    R.to_csv(a.out + ".csv", index=False)
    S.to_csv(a.out + "_splits.csv", index=False)
    json.dump(summ, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
