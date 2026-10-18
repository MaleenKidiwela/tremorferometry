#!/usr/bin/env python
"""bw_burst2006_recsec -- eye-check record sections (2-8 Hz) of the 20 Hz caches around the early-2006 burst peaks.

Read-only on the caches and mf/feat.  Single process, light.  Per day: (a) 2-8 Hz 3-C envelope (10-s RMS, divided by the
station's whole-day median) for every roster station with a cache, over a few hours around the peak, with the frozen
and era data detections (net_mad >= 9) as ticks; (b) 90-s zooms of the 2-8 Hz N component around the 3 strongest frozen
detections, all stations, each trace normalised to its own max in the zoom.
Out: svi_science/results/backward/first_look/burst2006/fig_recsec_<date>.png, station_env_hourly.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
import bw_p4env_detect as E

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
OUT = f"{BW}/first_look/burst2006"
FEAT = f"{BW}/mf/feat"
T = 9.0
# date, window start / end (UTC hour), label
CASES = [("2006-03-28", 2, 16, "Mar 28 burst"), ("2006-03-29", 6, 18, "Mar 29 burst"),
         ("2006-01-16", 0, 8, "Jan 16 burst"), ("2006-05-26", 0, 15, "May 26 burst"),
         ("2006-03-20", 6, 16, "quiet ref")]
FS = 20.0
SOS = butter(4, [2 / (FS / 2), 8 / (FS / 2)], "band", output="sos")
C = E.coords()


def day_trace(sta, d):
    p = E.cache_path(sta, d)
    if p is None:
        return None
    z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"])
    day0 = pd.Timestamp(d).value / 1e9
    n = int(86400 * FS)
    X = np.zeros((3, n)); i0 = int(round((t0 - day0) * FS)); a, b = max(i0, 0), min(i0 + x.shape[1], n)
    if b <= a:
        return None
    X[:, a:b] = x[:, a - i0:b - i0]
    X[~np.isfinite(X)] = 0
    gap = ~np.any(X != 0, axis=0)
    Y = sosfiltfilt(SOS, X, axis=1)
    Y[:, gap] = np.nan
    return Y


def dets(run, d):
    p = f"{FEAT}/{run}/{d}.parquet"
    if not os.path.exists(p):
        return pd.DataFrame(columns=["series", "t", "net_mad", "fam"])
    A = pd.read_parquet(p, columns=["series", "t", "net_mad", "fam"])
    return A[A.net_mad >= T]


def main():
    rows = []
    for d, h0, h1, lab in CASES:
        stas = [s for s in E.ROSTER if E.cache_path(s, d)]
        # order: by latitude (N on top)
        stas = sorted(stas, key=lambda s: C[s][0])
        Ys = {}
        for s in stas:
            Y = day_trace(s, d)
            if Y is not None:
                Ys[s] = Y
        stas = [s for s in stas if s in Ys]
        Df = dets("p5f_2006", d); De = dets("p5e_2006", d)
        day0 = pd.Timestamp(d).value / 1e9
        fig = plt.figure(figsize=(14, 11))
        ax = fig.add_axes([0.06, 0.40, 0.90, 0.55])
        k = int(10 * FS)
        for i, s in enumerate(stas):
            Y = Ys[s]
            env = np.sqrt(np.nanmean((Y ** 2)[:, : (Y.shape[1] // k) * k].reshape(3, -1, k), axis=(0, 2)))
            med = np.nanmedian(env)
            r = env / med
            for hh in range(24):
                seg = r[hh * 360:(hh + 1) * 360]
                rows.append(dict(date=d, sta=s, hour_utc=hh, env_rel_med=float(np.nanmedian(seg)) if np.isfinite(seg).any() else np.nan,
                                 env_rel_p90=float(np.nanpercentile(seg, 90)) if np.isfinite(seg).any() else np.nan))
            tt = (np.arange(len(r)) + 0.5) * 10 / 3600
            m = (tt >= h0) & (tt <= h1)
            ax.plot(tt[m], i + 0.12 * np.clip(r[m], 0, 8), lw=0.6, color="k")
        for _, q in Df[Df.series == "data"].iterrows():
            ax.plot([(q.t - day0) / 3600] * 2, [-1.6, -1.1], color="C0", lw=0.4)
        for _, q in Df[Df.series == "c1"].iterrows():
            ax.plot([(q.t - day0) / 3600] * 2, [-2.4, -1.9], color="C3", lw=0.4)
        for _, q in De[De.series == "data"].iterrows():
            ax.plot([(q.t - day0) / 3600] * 2, [-3.2, -2.7], color="C2", lw=0.4)
        ax.set_yticks(list(range(len(stas))) + [-1.35, -2.15, -2.95])
        ax.set_yticklabels([f"{s} ({C[s][0]:.2f}N)" for s in stas] + ["frozen data", "frozen ctrl", "era data"], fontsize=7)
        ax.set_xlim(h0, h1); ax.set_ylim(-3.5, len(stas) + 0.5)
        ax.set_xlabel("hour of UTC day (PST = UTC-8)")
        ax.set_title(f"{d} {lab}: 2-8 Hz 3-C envelope, 10-s RMS / station-day median (clipped at 8x), stations by latitude",
                     fontsize=9)
        # zooms
        top = Df[(Df.series == "data") & (Df.t >= day0 + h0 * 3600) & (Df.t <= day0 + h1 * 3600)]
        top = top.sort_values("net_mad", ascending=False)
        pick = []
        for _, q in top.iterrows():
            if all(abs(q.t - p.t) > 600 for p in pick):
                pick.append(q)
            if len(pick) == 3:
                break
        for j, q in enumerate(pick):
            a2 = fig.add_axes([0.06 + j * 0.31, 0.04, 0.28, 0.30])
            i0 = int((q.t - day0 - 30) * FS); i1 = i0 + int(90 * FS)
            for i, s in enumerate(stas):
                y = Ys[s][1, max(i0, 0):i1]
                if not np.isfinite(y).any():
                    continue
                y = y / np.nanmax(np.abs(y))
                a2.plot(np.arange(len(y)) / FS - 30, i + 0.45 * y, lw=0.4, color="k")
            a2.set_yticks(range(len(stas))); a2.set_yticklabels(stas, fontsize=6)
            a2.axvline(0, color="C0", lw=0.6)
            ts = pd.to_datetime(q.t, unit="s").strftime("%H:%M:%S")
            a2.set_title(f"N comp 2-8 Hz, fam {int(q.fam)}, net_mad {q.net_mad:.1f}, {ts} UTC", fontsize=7)
            a2.set_xlabel("s from detection (template start)", fontsize=7)
        fig.savefig(f"{OUT}/fig_recsec_{d}.png", dpi=100); plt.close(fig)
        print(d, "stations", stas, flush=True)
    pd.DataFrame(rows).to_csv(f"{OUT}/station_env_hourly.csv", index=False, float_format="%.3f")


if __name__ == "__main__":
    main()
