#!/usr/bin/env python
"""S12c pilot §4a: long stacks, yearly stability, LVL-peak test (P1). No time-change measurement here.
Usage: acf02_visibility.py --stations CN.PGC,CN.VGZ,CN.SNB --years 2018-2021 [--mask none|tremor]
Out: tables/visibility_<tag>.csv, figures/visibility_<STA>_<band>_<tag>.png
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import hilbert
from acf_common import *

ap = argparse.ArgumentParser()
ap.add_argument("--stations", default="CN.PGC,CN.VGZ,CN.SNB")
ap.add_argument("--years", default="2018-2021")
ap.add_argument("--mask", default="none")
ap.add_argument("--tag", default="pilot")
a = ap.parse_args()
y0, y1 = map(int, a.years.split("-")); years = list(range(y0, y1 + 1))
B = pd.read_csv(f"{WS}/tables/predicted_lag_bands.csv")
PAIRS = ["ZZ", "RR", "TT", "ZR+", "ZR-", "ZT+", "ZT-", "RT+", "RT-"]
rows = []


def agc(x, w=5.0):
    from scipy.ndimage import uniform_filter1d
    e = np.sqrt(uniform_filter1d(x ** 2, int(w * FS)))
    return x / (e + 1e-12)


for ns in a.stations.split(","):
    W = windows_for(ns)
    for band in ["b05_2", "b1_4"]:
        try:
            D, H, R = load(ns, years, band)
        except ValueError:
            print(f"no data {ns}", flush=True); continue
        ok = None
        ok = hour_mask(ns, D, R, a.mask)
        fig, axs = plt.subplots(len(PAIRS), 1, figsize=(11, 2.0 * len(PAIRS)), sharex=True)
        yrs = D.astype("datetime64[Y]").astype(int) + 1970
        for ip, pr in enumerate(PAIRS):
            base = pr.rstrip("+-")
            S, n = day_stack(H[base], ok)
            S = split_pair(S, pr)
            good = np.isfinite(S[:, 0])
            L = np.nanmean(S[good], 0)
            Y = {y: np.nanmean(S[good & (yrs == y)], 0) for y in years if (good & (yrs == y)).sum() > 30}
            ph = PAIR_PHASE.get(pr, "SS" if base == "RT" else "PS")
            env = np.abs(hilbert(L))
            ax = axs[ip]
            t = LAG_AUTO
            m = t >= 2
            sc = np.abs(agc(L)[m]).max()
            for y, ys in Y.items():
                ax.plot(t[m], agc(ys)[m] / sc, lw=0.6, alpha=0.6)
            ax.plot(t[m], agc(L)[m] / sc, "k", lw=1.0)
            for iface, col in [("top", "tab:blue"), ("bottom", "tab:red")]:
                rr = B[(B.station == ns) & (B.interface == iface)].iloc[0]
                ax.axvspan(rr[f"{ph}_lo"], rr[f"{ph}_hi"], color=col, alpha=0.15)
            for side, ls in [("pre", ":"), ("lvl", "--")]:
                for v in W[ph][side]:
                    ax.axvline(v, color="0.4", ls=ls, lw=0.7)
            ax.set_ylabel(pr, rotation=0, labelpad=18)
            # stability + peaks
            yl = list(Y.values())
            rec = dict(station=ns, band=band, pair=pr, phase=ph, n_days=int(good.sum()),
                       n_hours=int(np.isfinite(H[base][:, :, 0]).sum() if ok is None else (np.isfinite(H[base][:, :, 0]) & ok).sum()))
            for wn, (w0, w1) in list(W[ph].items()) + [("late", (27.0, 30.0))]:
                mm = (t >= w0) & (t < w1)
                ccs = [np.corrcoef(yl[i][mm], yl[j][mm])[0, 1] for i in range(len(yl)) for j in range(i + 1, len(yl))]
                rec[f"ccyr_{wn}"] = round(float(np.mean(ccs)), 3) if ccs else np.nan
                rec[f"rms_{wn}"] = float(np.sqrt(np.mean(L[mm] ** 2)))
            for iface in ["top", "bottom"]:
                rr = B[(B.station == ns) & (B.interface == iface)].iloc[0]
                lo, hi = rr[f"{ph}_lo"] - 0.3, rr[f"{ph}_hi"] + 0.3
                mm = (t >= lo) & (t <= hi)
                k = np.argmax(env[mm]); tp = t[mm][k]
                nb = (np.abs(t - tp) <= 4) & (t >= 2)
                ratio = env[mm][k] / np.median(env[nb])
                same = 0
                for ys in yl:
                    ey = np.abs(hilbert(ys)); same += int(abs(t[mm][np.argmax(ey[mm])] - tp) <= 0.3)
                rec[f"peak_{iface}_lag"] = round(float(tp), 2); rec[f"peak_{iface}_ratio"] = round(float(ratio), 2)
                rec[f"peak_{iface}_years"] = f"{same}/{len(yl)}"
                ax.plot(tp, 0.9, "v", color="k", ms=4)
            rows.append(rec)
            ax.set_title(f"{pr}: ccyr pre {rec['ccyr_pre']:.2f} lvl {rec['ccyr_lvl']:.2f} late {rec['ccyr_late']:.2f}; "
                         f"top pk {rec['peak_top_lag']} s ×{rec['peak_top_ratio']} ({rec['peak_top_years']}), "
                         f"bottom pk {rec['peak_bottom_lag']} s ×{rec['peak_bottom_ratio']} ({rec['peak_bottom_years']})",
                         fontsize=8, loc="left")
        axs[-1].set_xlabel("lag (s)"); axs[-1].set_xlim(2, 30)
        fig.suptitle(f"{ns} {band} {a.years} mask={a.mask}: long stack (black, AGC 5 s) and yearly stacks (colour); "
                     f"blue/red = predicted LVL top/bottom; dotted = pre window, dashed = LVL window", fontsize=9)
        fig.tight_layout(); fig.savefig(f"{WS}/figures/visibility_{ns}_{band}_{a.tag}_{a.mask}.png", dpi=110); plt.close(fig)
        print(f"{ns} {band} done", flush=True)
T = pd.DataFrame(rows)
T.to_csv(f"{WS}/tables/visibility_{a.tag}_{a.mask}.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
print(T[["station", "band", "pair", "n_days", "ccyr_pre", "ccyr_lvl", "ccyr_ctl_a", "ccyr_late", "peak_top_lag", "peak_top_ratio",
         "peak_top_years", "peak_bottom_lag", "peak_bottom_ratio", "peak_bottom_years"]].to_string(index=False))
