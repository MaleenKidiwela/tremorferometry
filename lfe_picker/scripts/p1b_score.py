#!/usr/bin/env python
"""P1b scoring — single-station S-peak recall vs peaks/station-day, from a stored scan.

Metric (Merlin P1b): threshold each station-day at the value of its R-th strongest declustered peak,
match the surviving peaks one-to-one to that station-day's labelled S arrivals within +-1 s, and
report recall = matched / labelled, per station, over the year.

THE ARITHMETIC CEILING is reported next to every number: with one CF series per station-day, R peaks
per day cannot match more than R events, so recall <= R x station-days / labelled events.  During an
ETS the north stations carry ~1,000-1,500 Bostock events per day, so 100 peaks/station-day caps
recall near 0.07 -- the template's 0.119 in P0-O is at 100 peaks per TEMPLATE-day (x116 families).

  python scripts/p1b_score.py --tag pnet_w05 --year 2012
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_scan as S

OUT = f"{P.RES}/p1"
RATES = [10, 30, 100, 300, 720, 1000, 3000, 6000, 11600]


def score(tag, year, tol=1.0, rates=RATES, npz=None, labels=None):
    z = np.load(npz or f"{OUT}/scan_{tag}_{year}.npz")
    days = pd.read_csv((npz or f"{OUT}/scan_{tag}_{year}.npz").replace(".npz", "_days.csv"))
    L = labels if labels is not None else pd.read_parquet(f"{P.RES}/labels_v0.parquet")
    ps = L[(L.kind == "pos_S") & (L.year == year)]
    hid, hsta = P.holdouts()
    rows = []
    for r in days.itertuples():
        key = f"{r.sta}_{r.doy}"
        if f"i_{key}" not in z:
            continue
        idx, val = z[f"i_{key}"], z[f"v_{key}"]
        pt = r.t0 + idx / P.FS
        g = ps[(ps.sta == r.sta) & (ps.doy == r.doy)]
        if not len(g):
            continue
        obs = g.t.to_numpy()
        o = np.argsort(obs)
        obs = obs[o]
        prov = g.provenance.to_numpy()[o]
        fam = g.fam.to_numpy()[o]
        mag = g.aux.to_numpy()[o]
        hidm = g.hidden_family.to_numpy()[o]
        frac = r.n / P.NDAY
        for R in rates:
            k = max(1, int(round(R * frac)))
            if len(val) == 0:
                hit = np.zeros(len(obs), bool)
            else:
                thr = np.sort(val)[::-1][min(k, len(val)) - 1]
                sel = pt[val >= thr][:k * 3]
                hit = S.one2one(sel, obs, tol=tol)
            rows.append(pd.DataFrame(dict(sta=r.sta, doy=r.doy, rate=R, hit=hit, prov=prov,
                                          fam=fam, mag=mag, hidden=hidm,
                                          npeak_avail=len(val), k=k, dayfrac=frac)))
    D = pd.concat(rows, ignore_index=True)
    D["held_out_station"] = D.sta.isin(hsta)
    return D, days


def summarize(D, days, rate=100):
    """The ceiling uses ONLY the station-days that carry labelled events -- those are the days the
    recall denominator comes from, so they are the days whose peak budget can match them."""
    d = D[D.rate == rate]
    per = d.groupby("sta").agg(n=("hit", "size"), rec=("hit", "mean")).reset_index()
    nd = d.drop_duplicates(["sta", "doy"]).groupby("sta").dayfrac.sum()
    per["station_days"] = per.sta.map(nd)
    per["ceiling"] = np.minimum(1.0, rate * per.station_days / per.n)
    per["rec_over_ceiling"] = per.rec / per.ceiling
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--tol", type=float, default=1.0)
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()
    D, days = score(a.tag, a.year, tol=a.tol)
    sfx = a.suffix or f"{a.tag}_{a.year}_tol{a.tol:g}"
    D.to_parquet(f"{OUT}/score_{sfx}.parquet", index=False)

    curve = D.groupby(["sta", "rate"]).hit.agg(["size", "mean"]).reset_index()
    curve.columns = ["sta", "rate", "n", "recall"]
    curve.to_csv(f"{OUT}/curve_{sfx}.csv", index=False)
    per = summarize(D, days, rate=100)
    per.to_csv(f"{OUT}/per_station_{sfx}.csv", index=False)
    print(f"\n== {a.tag} {a.year}  recall at 100 peaks/station-day (tol +-{a.tol} s) ==")
    print(per.to_string(index=False))
    for R in [720, 3000]:
        p = summarize(D, days, rate=R)
        print(f"\n-- at {R} peaks/station-day --")
        print(p[["sta", "rec", "ceiling"]].to_string(index=False))
    byp = D[D.rate == 100].groupby("prov").hit.agg(["size", "mean"])
    print("\nby provenance @100:\n", byp.to_string())
    byh = D[D.rate == 100].groupby("hidden").hit.agg(["size", "mean"])
    print("\nby hidden family @100:\n", byh.to_string())


if __name__ == "__main__":
    main()
