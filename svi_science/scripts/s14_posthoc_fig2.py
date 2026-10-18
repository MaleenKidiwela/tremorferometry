#!/usr/bin/env python
"""s14_posthoc_fig2 -- POST HOC: station-episode mean CC residual (after n + noise) beside station-episode mean dv/v
anomaly, per station, over time (SVI 7 stations; WA 6 sites with most cells). Uses tables/posthoc_station_episode_*.csv."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import s14_common as C
fig, axs = plt.subplots(4, 1, figsize=(14, 12), sharex=True)
for i, r in enumerate(("SVI", "WA_NW")):
    S = pd.read_csv(f"{C.TB}/posthoc_station_episode_{r}.csv")
    E = C.load_episode(r)[["episode", "date"]].groupby("episode").date.median()
    S["date"] = S.episode.map(E)
    U = C.units(r)
    stas = S.groupby("netsta").k.sum().sort_values(ascending=False).index[:7]
    a, b = axs[2 * i], axs[2 * i + 1]
    for ax in (a, b):
        for _, u in U[U.cls == "major"].iterrows():
            ax.axvspan(u.start, u.end + pd.Timedelta(days=1), color="#f4b183", alpha=0.3, lw=0)
    cols = plt.cm.tab10(np.arange(len(stas)))
    for c, ns in zip(cols, stas):
        q = S[S.netsta == ns].sort_values("date")
        a.plot(q.date, q.m, "-o", ms=2.5, lw=0.7, color=c, label=ns)
        b.plot(q.date, q.dv, "-o", ms=2.5, lw=0.7, color=c)
    a.axhline(0, color="k", lw=0.4); b.axhline(0, color="k", lw=0.4)
    a.set_ylabel("CC residual (Fisher z)"); b.set_ylabel("dv/v anomaly (%)"); b.set_ylim(-0.6, 0.6)
    a.set_title(f"{r} (POST HOC): station-episode mean over families (≥ 4) — CC after stack size + noise", fontsize=9)
    b.set_title(f"{r}: same station-episodes — mean dv/v anomaly (2–5 s)", fontsize=9)
    a.legend(fontsize=7, ncol=7, loc="lower left")
fig.savefig(f"{C.FIG}/figC2_posthoc_station_episode.png", dpi=120, bbox_inches="tight")
