#!/usr/bin/env python
"""Record-section review sheets.

  --mode xs     : cross-station detections (results/xs/<date>_det.parquet): the raw 1-10 Hz
                  window on every pool station, N/E/Z, sorted by distance from ALP0, S times of the
                  detecting stations marked.  Plus a detections-per-hour panel for the days.
  --mode stacks : family stacks npz (keys '<grp>_<STA>', (3, n) over LONG_PRE/LONG_LEN), one
                  figure per family, picks.csv overlaid if given.

  python scripts/xs_review.py --mode xs --xs results/xs --n 6 --out figures/xs
  python scripts/xs_review.py --mode stacks --stacks results/hypo/stacks_loc.npz --picks results/hypo/picks.csv --out figures/review
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import xs_common as C

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LONG_PRE, LONG_LEN = 12.0, 44.0


def dist_order(names, ref="XA.ALP0"):
    la0, lo0 = C.COORD[ref]
    key = lambda s: np.hypot((C.COORD[s][1] - lo0) * np.cos(np.radians(la0)), C.COORD[s][0] - la0) * 111.2
    return sorted(names, key=key), {s: key(s) for s in names}


def plot_xs(xs_dir, dates, n, out):
    os.makedirs(out, exist_ok=True)
    D = pd.concat([pd.read_parquet(f"{xs_dir}/{d}_det.parquet").assign(date=d) for d in dates], ignore_index=True)
    # detections per hour
    fig, ax = plt.subplots(figsize=(10, 2.6))
    h = ((D.t_win - pd.Timestamp(dates[0]).timestamp()) // 3600).astype(int)
    ax.bar(*np.unique(h, return_counts=True), width=0.9, color="0.3")
    ax.set_xlim(0, 24 * len(dates))
    ax.set_xlabel(f"hour since {dates[0]} 00:00 UTC")
    ax.set_ylabel("detections / h")
    ax.set_title(f"cross-station detections, {len(D)} total")
    fig.tight_layout()
    fig.savefig(f"{out}/xs_per_hour.png", dpi=110)
    plt.close(fig)
    stas = [s for s in C.STATIONS.netsta if os.path.isdir(f"{C.CACHE}/{s}")]
    order, dist = dist_order(stas)
    sel = D.sort_values("maxwin", ascending=False).head(n)
    for r in sel.itertuples():
        pre, ln = 6.0, 30.0
        fig, axs = plt.subplots(1, 3, figsize=(13, 9), sharey=True)
        for k, s in enumerate(order):
            z = C.load_day(s, r.date)
            if z is None:
                continue
            t0, X, m = z
            i0 = int(round((r.t_win - pre - t0) * C.FS))
            nn = int(ln * C.FS)
            if i0 < 0 or i0 + nn > X.shape[1] or not m[i0:i0 + nn].all():
                continue
            seg = X[:, i0:i0 + nn].astype(float)
            mx = np.abs(seg[1:]).max() if np.abs(seg[1:]).max() > 0 else np.abs(seg).max()
            tt = np.arange(nn) / C.FS - pre
            for c, ax in zip((1, 2, 0), axs):
                if mx > 0 and np.abs(seg[c]).max() > 0:
                    ax.plot(tt, k + 0.45 * seg[c] / mx, lw=0.5, color="k")
                col = f"tS_{s}"
                if col in D.columns and np.isfinite(D.at[r.Index, col]):
                    ax.plot([D.at[r.Index, col] - r.t_win] * 2, [k - 0.45, k + 0.45], color="r", lw=1.5)
        for ax, lab in zip(axs, ("N", "E", "Z")):
            ax.axvspan(0, 18, color="0.9", zorder=-1)
            ax.set_title(lab)
            ax.set_xlabel("s after detection-window start")
        axs[0].set_yticks(range(len(order)))
        axs[0].set_yticklabels([f"{s} {dist[s]:.0f}km" for s in order], fontsize=7)
        fig.suptitle(f"{pd.Timestamp(r.t_win, unit='s')}  maxwin {r.maxwin:.2f}  n_ch {r.n_ch}  "
                     f"sigr_max {r.sigr_max:.2f}  (grey = 18 s detection window; red = S times; "
                     f"each station scaled to its horizontal max)", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{out}/xs_{pd.Timestamp(r.t_win, unit='s'):%Y%m%dT%H%M%S}.png", dpi=100)
        plt.close(fig)


def plot_stacks(npz, picks, out, ndet=None):
    os.makedirs(out, exist_ok=True)
    z = np.load(npz)
    ST = {}
    for k in z.files:
        g, s = k.split("_", 1)
        ST.setdefault(int(g), {})[s] = z[k]
    P = pd.read_csv(picks) if picks and os.path.exists(picks) else pd.DataFrame()
    bare = {s.split(".")[1]: s for s in C.STATIONS.netsta}
    for g, T in sorted(ST.items()):
        names = [bare.get(s, s) for s in T]
        order, dist = dist_order([s for s in names if s in C.COORD])
        fig, axs = plt.subplots(1, 3, figsize=(13, 9), sharey=True)
        for k, s in enumerate(order):
            st = T[s.split(".")[1]] if s.split(".")[1] in T else T[s]
            tt = np.arange(st.shape[1]) / C.FS - LONG_PRE
            mx = np.abs(st).max()
            for c, ax in zip((1, 2, 0), axs):
                if mx > 0:
                    ax.plot(tt, k + 0.45 * st[c] / mx, lw=0.5, color="k")
            if len(P):
                q = P[(P.grp == g) & (P.sta == s.split(".")[1])]
                for r in q.itertuples():
                    for ph, col in (("tS", "r"), ("tP", "b")):
                        v = getattr(r, ph, np.nan)
                        if np.isfinite(v):
                            for ax in axs:
                                ax.plot([v - LONG_PRE] * 2, [k - 0.45, k + 0.45], color=col, lw=1.2)
        for ax, lab in zip(axs, ("N", "E", "Z")):
            ax.set_title(lab)
            ax.set_xlabel("s after template start")
        axs[0].set_yticks(range(len(order)))
        axs[0].set_yticklabels([f"{s} {dist[s]:.0f}km" for s in order], fontsize=7)
        fig.suptitle(f"family {g}" + (f"  n_det {ndet.get(g)}" if ndet else ""), fontsize=10)
        fig.tight_layout()
        fig.savefig(f"{out}/fam_{g:06d}.png", dpi=100)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["xs", "stacks"])
    ap.add_argument("--xs", default=f"{C.RES}/xs")
    ap.add_argument("--dates", default=",".join(C.DETECT_DAYS))
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--stacks")
    ap.add_argument("--picks")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if a.mode == "xs":
        plot_xs(a.xs, a.dates.split(","), a.n, a.out)
    else:
        plot_stacks(a.stacks, a.picks, a.out)


if __name__ == "__main__":
    main()
