#!/usr/bin/env python
"""P3 STEP 0(a) — chance-correct every P1b single-station recall curve.

Merlin's flaw #4: the P1b headline "4.1x the template cc at the template's own peak budget" is NOT
chance corrected.  At 11,600 declustered peaks/station-day the +-1 s matching tolerance covers
>= 27 % of the day, so a large part of the recall is accidental.

Correction recipe = EXACTLY scripts/p2_score.py's: shift the DETECTOR's times by +-307 s, rematch
one-to-one against the same labels, average the two, and report

    recall_cc = (recall - chance) / (1 - chance).

Nothing is re-run and nothing is re-tuned: this reads the stored scans of P1b
(results/p1/scan_<tag>_<year>.npz) and re-scores them.  The 2011 evaluation was spent on
2026-09-17 and is recorded in the ledger; this is bookkeeping on its stored product, so no new
referee budget is consumed.

    python scripts/q_p1b_chance.py --out results/p3/p1b_chance_corrected.csv
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_scan as S

OUT = f"{P.RES}/p1"
RATES = [10, 30, 100, 300, 720, 1000, 3000, 6000, 11600]
SHIFT = 307.0


def score_tag(tag, year, labels, tol=1.0, rates=RATES, shift=SHIFT):
    npz = f"{OUT}/scan_{tag}_{year}.npz"
    if not os.path.exists(npz):
        return None
    z = np.load(npz)
    days = pd.read_csv(npz.replace(".npz", "_days.csv"))
    ps = labels[(labels.kind == "pos_S") & (labels.year == year)]
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
        obs = np.sort(g.t.to_numpy())
        frac = r.n / P.NDAY
        for R in rates:
            k = max(1, int(round(R * frac)))
            if len(val) == 0:
                h = hp = hm = np.zeros(len(obs), bool)
            else:
                thr = np.sort(val)[::-1][min(k, len(val)) - 1]
                sel = pt[val >= thr][:k * 3]
                h = S.one2one(sel, obs, tol=tol)
                hp = S.one2one(sel + shift, obs, tol=tol)
                hm = S.one2one(sel - shift, obs, tol=tol)
            rows.append(dict(sta=r.sta, doy=r.doy, rate=R, n=len(obs),
                             hit=int(h.sum()), hit_p=int(hp.sum()), hit_m=int(hm.sum()),
                             npeak=int(len(val)), k=k, dayfrac=frac))
    if not rows:
        return None
    D = pd.DataFrame(rows)
    g = D.groupby(["sta", "rate"]).agg(n=("n", "sum"), hit=("hit", "sum"),
                                       hit_p=("hit_p", "sum"), hit_m=("hit_m", "sum"),
                                       station_days=("dayfrac", "sum")).reset_index()
    g["tag"] = tag
    g["year"] = year
    g["recall"] = g.hit / g.n
    g["chance"] = 0.5 * (g.hit_p + g.hit_m) / g.n
    g["recall_cc"] = (g.recall - g.chance) / (1.0 - g.chance)
    g["ceiling"] = np.minimum(1.0, g.rate * g.station_days / g.n)
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", default="pnet_w05,pnet_w10,lin_w05")
    ap.add_argument("--years", default="2010,2012,2011")
    ap.add_argument("--tol", type=float, default=1.0)
    ap.add_argument("--out", default=f"{P.RES}/p3/p1b_chance_corrected.csv")
    a = ap.parse_args()
    L = pd.read_parquet(f"{P.RES}/labels_v0.parquet")
    out = []
    for tag in a.tags.split(","):
        for y in [int(x) for x in a.years.split(",")]:
            g = score_tag(tag, y, L, tol=a.tol)
            if g is None:
                print(f"[skip] {tag} {y}")
                continue
            out.append(g)
            r = g[(g.sta == "PGC")].set_index("rate")
            print(f"[{tag} {y}] PGC raw/chance/cc  "
                  + "  ".join(f"{R}:{r.recall[R]:.4f}/{r.chance[R]:.4f}/{r.recall_cc[R]:.4f}"
                              for R in [100, 720, 11600] if R in r.index), flush=True)
    G = pd.concat(out, ignore_index=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    G.to_csv(a.out, index=False)
    print(f"[done] -> {a.out}  ({len(G)} rows)")


if __name__ == "__main__":
    main()
