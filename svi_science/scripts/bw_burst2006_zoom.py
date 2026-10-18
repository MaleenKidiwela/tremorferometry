#!/usr/bin/env python
"""bw_burst2006_zoom -- 20-min close-up (2-8 Hz, 1-s RMS envelope, each station / its own 20-min median) on burst days,
with frozen data detections; read-only on the caches and mf/feat.  Out: burst2006/fig_zoom_<date>_<hhmm>.png"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import bw_burst2006_recsec as R
OUT = R.OUT
CASES = [("2006-03-29", 9.5), ("2006-03-28", 9.0), ("2006-03-20", 8.3), ("2006-01-16", 1.5)]
for d, h0 in CASES:
    stas = sorted([s for s in R.E.ROSTER if R.E.cache_path(s, d)], key=lambda s: R.C[s][0])
    Df = R.dets("p5f_2006", d); De = R.dets("p5e_2006", d); day0 = pd.Timestamp(d).value / 1e9
    a0, a1 = int(h0 * 3600 * R.FS), int((h0 * 3600 + 1200) * R.FS)
    fig, ax = plt.subplots(1, 1, figsize=(13, 8))
    i = 0; lab = []
    for s in stas:
        Y = R.day_trace(s, d)
        if Y is None: continue
        y = Y[:, a0:a1]; k = int(R.FS)
        e = np.sqrt(np.nanmean((y ** 2).reshape(3, -1, k), axis=(0, 2)))
        e = e / np.nanmedian(e)
        ax.plot(h0 * 60 + np.arange(len(e)) / 60, i + 0.25 * np.clip(e, 0, 6), lw=0.6, color="k"); lab.append(s); i += 1
    for df, y0, c in ((Df[Df.series == "data"], -1.0, "C0"), (De[De.series == "data"], -1.8, "C2")):
        tm = (df.t - day0) / 60
        tm = tm[(tm >= h0 * 60) & (tm <= h0 * 60 + 20)]
        ax.vlines(tm, y0 - 0.3, y0 + 0.3, color=c, lw=0.8)
    ax.set_yticks(list(range(len(lab))) + [-1.0, -1.8]); ax.set_yticklabels(lab + ["frozen data", "era data"], fontsize=7)
    ax.set_xlabel("minute of UTC day"); ax.set_title(f"{d} {int(h0):02d}:{int(h0 % 1 * 60):02d} UTC +20 min: 2-8 Hz 3-C 1-s RMS / 20-min median (clip 6x)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_zoom_{d}_{int(h0):02d}{int(h0 % 1 * 60):02d}.png", dpi=100); plt.close(fig)
    print(d, h0, flush=True)
