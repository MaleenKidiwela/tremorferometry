#!/usr/bin/env python
"""Validation plots for ab_mute.py on a few NW days (2026-09-26).

Per day: for W020, W060 and B013 (site ringing) the hour with the most muted time: 1-8 Hz trace (N),
muted spans shaded red, real-LFE detection windows of run 1 (33 eye-classified families, det_it14,
[t, t+15 s] on that station if it is a template station) as green bars; below, R (station ratio) and
Q (neighbour-site array ratio) with the thresholds.  Plus zooms: 3 muted bursts and 3 real-LFE
detections per station-day (before / after muting+notching).  And the per-day loss table.

  python scripts/ab_mute_validate.py --dates 2008-05-08,2008-05-15,2011-08-03
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_mute as M                                                         # noqa: E402

import matplotlib                                                            # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                              # noqa: E402

STAS = ["W020", "W060", "B013"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{A.RES}/NW_mute/mute")
    ap.add_argument("--dates", default="2008-05-08,2008-05-15,2011-08-03")
    ap.add_argument("--K", type=float, default=15.0)
    ap.add_argument("--qmax", type=float, default=2.0)
    ap.add_argument("--qlong", type=float, default=1.3)
    a = ap.parse_args()
    det = pd.read_parquet(f"{A.RES}/NW/grow/det_it14.parquet")
    ec = pd.read_csv(f"{A.RES}/NW/hypo_28/review/eye_class.csv")
    lfe = set(ec.grp[ec["class"] == "lfe"])
    ts = pd.read_csv(f"{A.RES}/NW/grow/tpl_it14_stations.csv")
    stmap = {int(r.cid): r.stations.split(",") for r in ts.itertuples()}
    det = det[det.grp.isin(lfe)]
    vd = f"{a.out}/validate"
    os.makedirs(vd, exist_ok=True)
    rows = []
    for date in a.dates.split(","):
        z = np.load(f"{a.out}/energy/{date}.npz")
        stas = list(z["stas"])
        R, Q, QL = M.ratios(z["e"], stas)
        MU = M.mute_from(R, Q, a.K, a.qmax, QL=QL, qlong=a.qlong)
        mid = pd.Timestamp(date).timestamp()
        dd = det[det.date == date]
        # loss table over ALL template stations for this day
        for r in dd.itertuples():
            for s in stmap.get(int(r.grp), []):
                i = stas.index(s)
                k0, k1 = int((r.t - mid) // 2), int(np.ceil((r.t + 15 - mid) / 2))
                if not np.isfinite(R[i, k0:k1]).all():
                    continue
                rows.append(dict(date=date, grp=r.grp, t=r.t, sta=s, muted=bool(MU[i, k0:k1].any())))
        fig, axs = plt.subplots(len(STAS) * 2, 1, figsize=(16, 4.2 * len(STAS)),
                                gridspec_kw=dict(height_ratios=[2, 1] * len(STAS)))
        zoom = []
        for si, s in enumerate(STAS):
            i = stas.index(s)
            hr = int(np.argmax([MU[i, h * 1800:(h + 1) * 1800].sum() for h in range(24)]))
            raw = M._raw_day(s, date)
            ax, bx = axs[2 * si], axs[2 * si + 1]
            if raw is None:
                ax.set_title(f"{s} {date}: no data"); continue
            t0, X, m = raw
            a0 = int(round((mid + hr * 3600 - t0) * A.FS))
            seg = X[1, a0:a0 + 144000]
            tt = np.arange(seg.size) / A.FS / 60
            ax.plot(tt, seg, "k", lw=0.3)
            kk = np.arange(hr * 1800, (hr + 1) * 1800)
            for k in kk[MU[i, kk]]:
                ax.axvspan((k - hr * 1800) * 2 / 60, (k - hr * 1800 + 1) * 2 / 60, color="r", alpha=0.3, lw=0)
            y = np.nanpercentile(np.abs(seg), 99.9) if seg.size else 1
            nL = nM = 0
            for r in dd.itertuples():
                if s not in stmap.get(int(r.grp), []):
                    continue
                x = (r.t - mid - hr * 3600) / 60
                if 0 <= x <= 60:
                    k0, k1 = int((r.t - mid) // 2), int(np.ceil((r.t + 15 - mid) / 2))
                    hit = MU[i, k0:k1].any()
                    ax.plot([x, x + 0.25], [y * 1.05] * 2, color="m" if hit else "g", lw=4)
                    nL += 1; nM += int(hit)
            ax.set_xlim(0, 60)
            ax.set_ylim(-y * 1.15, y * 1.15)
            ax.set_title(f"{s} N {date} {hr:02d}h (1-8 Hz). red = muted ({MU[i, kk].mean():.1%} of hour); "
                         f"green = real-LFE detection windows kept ({nL - nM}), magenta = muted ({nM})",
                         fontsize=9)
            bx.semilogy(np.arange(1800) * 2 / 60, R[i, kk], lw=0.5, label="R (station)")
            bx.semilogy(np.arange(1800) * 2 / 60, Q[i, kk], lw=0.8, label="Q (neighbour sites, 6 s)")
            bx.semilogy(np.arange(1800) * 2 / 60, QL[i, kk], lw=0.8, label="QL (neighbour sites, 60 s)")
            bx.axhline(a.qlong, color="g", ls=":")
            bx.axhline(a.K, color="r", ls=":"); bx.axhline(a.qmax, color="b", ls=":")
            bx.set_xlim(0, 60); bx.legend(fontsize=7, loc="upper right"); bx.set_xlabel("minutes")
            # zoom candidates
            mk = np.flatnonzero(MU[i] & np.isfinite(R[i]))
            if mk.size:
                # the strongest 3 distinct bursts
                order = mk[np.argsort(-R[i, mk])]
                picked = []
                for k in order:
                    if all(abs(k - p) > 60 for p in picked):
                        picked.append(k)
                    if len(picked) == 3:
                        break
                zoom += [(s, "muted burst", mid + k * 2) for k in picked]
            dl = [r.t for r in dd.itertuples() if s in stmap.get(int(r.grp), [])]
            if dl:
                rng = np.random.default_rng(0)
                zoom += [(s, "real-LFE det", t) for t in rng.choice(dl, min(3, len(dl)), replace=False)]
        plt.tight_layout()
        plt.savefig(f"{vd}/hour_{date}.png", dpi=80)
        plt.close()
        # zooms: before (grey) and after (black) mute + notch
        os.environ["AB_MUTE_DIR"] = a.out
        fig, axs = plt.subplots(len(zoom), 1, figsize=(12, 1.5 * len(zoom) + 1), squeeze=False)
        for ax, (s, lab, tc) in zip(axs[:, 0], zoom):
            b = M._raw_day(s, date)
            c = A.load_day(s, date)
            for (t0, X, m), col, lw in [(b, "0.6", 1.0), (c, "k", 0.6)]:
                i0 = int(round((tc - 20 - t0) * A.FS))
                i0 = max(i0, 0)
                seg = X[1, i0:i0 + int(60 * A.FS)]
                ax.plot(np.arange(seg.size) / A.FS - 20, seg, color=col, lw=lw)
            if lab == "real-LFE det":
                ax.axvspan(0, 15, color="g", alpha=0.1)
            ax.set_title(f"{s} {date} {pd.Timestamp(tc, unit='s'):%H:%M:%S} {lab} (grey raw, black muted+notched)",
                         fontsize=8)
        plt.tight_layout()
        plt.savefig(f"{vd}/zoom_{date}.png", dpi=80)
        plt.close()
        os.environ.pop("AB_MUTE_DIR", None)
    df = pd.DataFrame(rows)
    df.to_csv(f"{vd}/lfe_window_loss.csv", index=False)
    g = df.groupby("date").muted.agg(["size", "mean"])
    g["det_ge1"] = df.groupby(["date", "t"]).muted.any().groupby("date").mean()
    print(g.to_string())
    print("by station:", df.groupby("sta").muted.mean().round(3).sort_values().to_dict())


if __name__ == "__main__":
    main()
