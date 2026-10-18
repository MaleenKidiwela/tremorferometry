#!/usr/bin/env python
"""QC / K choice for ab_mute.py burst muting (2026-09-26).

For the eye-classified families of a finished run (default run 1: NW/grow det_it14 + tpl_it14 stations,
classes from NW/hypo_28/review/eye_class.csv: 33 lfe / 136 noise), every (detection, template station)
pair defines the correlated window [t, t+15 s] on that station.  For each (K, qmax) we report
  loss_lfe   = fraction of real-LFE detection-station windows touched by a mute (the window would be
               invalid, exactly as a gap), and the fraction of real detections losing >= 1 / >= 3 stations,
  hit_noise  = the same fraction for the noise families,
  muted      = overall muted fraction of live 2 s windows (all pool stations, all scan days).
Also writes the distributions of R (station ratio) and Q (neighbour-site array ratio) for all windows,
real-LFE windows and noise windows (figure mute_ratio_distributions.png).

  python scripts/ab_mute_qc.py --out results/ab_royer/NW_mute/mute
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_mute as M                                                         # noqa: E402

KS = [10, 15, 20, 30, 50]
QS = [(1.5, np.inf), (2.0, np.inf), (1.5, 1.5), (2.0, 1.5), (2.0, 1.3), (1.5, 1.3), (3.0, 1.3)]
TPL = 15.0


def _day(args):
    date, out, det, stmap = args
    z = np.load(f"{out}/energy/{date}.npz")
    stas = list(z["stas"])
    R, Q, QL = M.ratios(z["e"], stas)
    live = np.isfinite(R)
    res = dict(date=date, live=live.sum(1))
    for K in KS:
        for q in QS:
            mu = M.mute_from(R, Q, K, q[0], QL=QL, qlong=q[1]) & live
            res[("muted", K, q)] = mu.sum(1)
    # per detection-station windows
    mid = pd.Timestamp(date).timestamp()
    recs = []
    if len(det):
        k0 = np.floor((det.t.values - mid) / M.WIN_S).astype(int)
        k1 = np.ceil((det.t.values + TPL - mid) / M.WIN_S).astype(int)
        MU = {(K, q): M.mute_from(R, Q, K, q[0], QL=QL, qlong=q[1]) for K in KS for q in QS}
        for j, (g, c) in enumerate(zip(det.grp.values, det.cls.values)):
            a, b = max(k0[j], 0), min(k1[j], M.NK)
            for s in stmap.get(int(g), []):
                i = stas.index(s)
                if not live[i, a:b].all():
                    continue
                rr = R[i, a:b]
                row = [int(g), c, j, s, float(np.nanmax(rr)), float(np.nanmin(Q[i, a:b]))]
                row += [bool(MU[kq][i, a:b].any()) for kq in MU]
                recs.append(row)
    # random sample of R, Q for the distribution figure
    rng = np.random.default_rng(abs(hash(date)) % 2**32)
    ii = np.flatnonzero(live.ravel())
    ii = rng.choice(ii, min(4000, ii.size), replace=False)
    res["sampleR"] = R.ravel()[ii].astype(np.float32)
    res["sampleQ"] = Q.ravel()[ii].astype(np.float32)
    return res, recs, stas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{A.RES}/NW_mute/mute")
    ap.add_argument("--run", default=f"{A.RES}/NW/grow")
    ap.add_argument("--it", default="14")
    ap.add_argument("--eye", default=f"{A.RES}/NW/hypo_28/review/eye_class.csv")
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    det = pd.read_parquet(f"{a.run}/det_it{a.it}.parquet")
    ec = pd.read_csv(a.eye)
    cls = dict(zip(ec.grp, ec["class"]))
    det = det[det.grp.isin(cls)].copy()
    det["cls"] = det.grp.map(cls)
    ts = pd.read_csv(f"{a.run}/tpl_it{a.it}_stations.csv")
    stmap = {int(r.cid): r.stations.split(",") for r in ts.itertuples()}
    dates = A.scan_days("NW")
    by = dict(tuple(det.groupby("date")))
    jobs = [(d, a.out, by.get(d, det.iloc[:0]), stmap) for d in dates]
    tot, recs, SR, SQ = {}, [], [], []
    cols = ["grp", "cls", "j", "sta", "Rmax", "Qmin"] + [f"m_K{K}_q{q[0]}_{q[1]}" for K in KS for q in QS]
    with ProcessPoolExecutor(a.workers) as ex:
        for res, rr, stas in ex.map(_day, jobs):
            for k, v in res.items():
                if k in ("date", "sampleR", "sampleQ"):
                    continue
                tot[k] = tot.get(k, 0) + v
            for r in rr:
                recs.append([r[0], r[1], f"{res['date']}_{r[2]}"] + r[3:])
            SR.append(res["sampleR"]); SQ.append(res["sampleQ"])
    df = pd.DataFrame(recs, columns=cols)
    df.to_parquet(f"{a.out}/qc_detwin.parquet")
    live = tot["live"]
    rows = []
    for K in KS:
        for q in QS:
            c = f"m_K{K}_q{q[0]}_{q[1]}"
            r = dict(K=K, qmax=q[0], qlong=q[1], muted_all=tot[("muted", K, q)].sum() / live.sum())
            for cl in ["lfe", "noise"]:
                x = df[df.cls == cl]
                r[f"win_{cl}"] = x[c].mean()
                per = x.groupby("j")[c].sum()
                r[f"det_ge1_{cl}"] = (per >= 1).mean()
                r[f"det_ge3_{cl}"] = (per >= 3).mean()
            rows.append(r)
    sw = pd.DataFrame(rows)
    sw.to_csv(f"{a.out}/qc_sweep.csv", index=False)
    print(sw.round(4).to_string())
    for K in KS:
        for q in QS:
            pd.Series(tot[("muted", K, q)] / np.maximum(live, 1), index=stas).to_csv(
                f"{a.out}/qc_station_frac_K{K}_q{q[0]}_{q[1]}.csv", header=["frac"])
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    SR, SQ = np.concatenate(SR), np.concatenate(SQ)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    bins = np.logspace(-1.5, 4, 90)
    for v, lab in [(SR, "all live 2 s windows (sample)"),
                   (df[df.cls == "lfe"].Rmax.values, "real-LFE det. windows (max over 15 s)"),
                   (df[df.cls == "noise"].Rmax.values, "noise-family det. windows (max over 15 s)")]:
        v = v[np.isfinite(v)]
        ax[0].hist(np.clip(v, bins[0], bins[-1]), bins, density=True, histtype="step", lw=1.5,
                   label=lab)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("R = window energy / channel daily median (max over components)")
    ax[0].legend(fontsize=8)
    bq = np.logspace(-1.3, 2, 80)
    for v, lab in [(SQ, "all live windows"), (SQ[SR > 20], "windows with R > 20"),
                   (df[(df.cls == "lfe") & (df.Rmax > 20)].Qmin.values, "real-LFE windows with R>20 (min Q)")]:
        v = v[np.isfinite(v)]
        ax[1].hist(np.clip(v, bq[0], bq[-1]), bq, density=True, histtype="step", lw=1.5, label=lab)
    ax[1].axvline(2, color="k", ls=":")
    ax[1].set_xscale("log"); ax[1].set_xlabel("Q = neighbour-site array ratio (median of 5 nearest other sites)")
    ax[1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{a.out}/mute_ratio_distributions.png", dpi=110)


if __name__ == "__main__":
    main()
