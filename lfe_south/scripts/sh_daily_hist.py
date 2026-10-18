#!/usr/bin/env python
"""Daily event counts through Shelly (2025)'s catalogue 2018-2024, the days we ran, and our forward run to 2026-09.
Reads saved files only. 2026-10-02."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.dates as mdates
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_south/scripts")
from sh_v3_common import his_catalogue
R = "/home/jovyan/tremorferometry/lfe_south/results"
H = his_catalogue(); tc = [c for c in H.columns if c in ("t", "time", "datetime")][0]
H["date"] = pd.to_datetime(H[tc]).dt.tz_localize(None).dt.normalize() if pd.to_datetime(H[tc]).dt.tz is not None else pd.to_datetime(H[tc]).dt.normalize()
his = H.groupby("date").size()
D = pd.read_csv(f"{R}/sh_conf/day_table.csv", parse_dates=["date"])
v3b = D[D.run == "v3b"]; build = set(v3b[v3b.build_any].date); test = set(v3b[~v3b.build_any].date)
F = pd.read_csv(f"{R}/sh_kmrfix/t6_daily_per_family_v3.csv.gz", parse_dates=["date"])
fw = F.groupby("date")[["S1_v3_scaled", "sumP_v2b_scaled", "S1_v2"]].sum(min_count=1).dropna(how="all")
days = pd.date_range("2018-01-01", "2026-09-30", freq="D")
his = his.reindex(days, fill_value=0)

ovl = (fw.index >= "2023-12-01") & (fw.index < "2024-12-01"); ext = fw.index >= "2024-12-01"
jo = pd.concat([his.rename("his"), fw.sumP_v2b_scaled.rename("ours")], axis=1).loc["2023-12-01":"2024-11-30"]
r = jo.corr().iloc[0, 1]
fig = plt.figure(figsize=(17, 12)); gs = fig.add_gridspec(3, 1, height_ratios=[4, 0.7, 3], hspace=0.35)
a = fig.add_subplot(gs[0])
a.bar(his.index, his.values, width=1.0, color="0.35", lw=0, label=f"Shelly catalogue ({int(his.sum()):,} events, 2018-01 → 2024-11)")
a.bar(fw.index[ext], fw.sumP_v2b_scaled[ext], width=1.0, color="C0", lw=0, label="ours, expected repeats ΣP scaled to his (2024-12 → 2026-09)")
a.step(fw.index[ovl], fw.sumP_v2b_scaled[ovl], where="mid", color="C0", lw=0.6, alpha=0.9, label="ours in the overlap year")
for d, l in [("2024-12-01", "his catalogue ends"), ("2025-09-17", "KMR swap → 16 stations")]:
    a.axvline(pd.Timestamp(d), color="k", ls=":", lw=1); a.text(pd.Timestamp(d), 222, " " + l, fontsize=8)
a.set_ylim(0, 240); a.set_ylabel("events per day"); a.legend(loc="upper left", fontsize=9)
a.set_title("Daily LFE events, all 27 Mendocino families, 2018-01 → 2026-09 (daily bins)", loc="left", fontsize=12)
a.xaxis.set_major_locator(mdates.YearLocator()); a.xaxis.set_major_formatter(mdates.DateFormatter("%Y")); a.set_xlim(days[0], days[-1])
b = fig.add_subplot(gs[1], sharex=a)
for st, c, y, lab in [(build, "C1", 2, f"template-build days ({len(build)})"), (test, "C2", 1, f"test days scored ({len(test)})")]:
    x = sorted(st); b.scatter(x, [y] * len(x), marker="|", s=120, color=c, label=lab)
b.axvspan(pd.Timestamp("2023-12-01"), pd.Timestamp("2026-09-29"), color="C0", alpha=0.15, label="continuous forward run (every day)")
b.set_yticks([1, 2]); b.set_yticklabels(["test", "build"]); b.set_ylim(0.3, 2.7); b.legend(loc="upper left", fontsize=8, ncol=3)
b.set_title("Days our detector ran", loc="left", fontsize=11)
c = fig.add_subplot(gs[2])
c.bar(jo.index, jo.his, width=1.0, color="0.6", lw=0, label="his")
c.step(jo.index, jo.ours, where="mid", color="C0", lw=0.9, label=f"ours (ΣP, scaled) — daily r = {r:.2f}")
c.set_ylabel("events per day"); c.legend(loc="upper left", fontsize=9); c.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
c.set_title("Zoom: overlap year 2023-12 → 2024-11, his vs ours day by day", loc="left", fontsize=11)
out = f"{R}/sh_daily/daily_events_shelly.png"; fig.savefig(out, dpi=110, bbox_inches="tight"); print(out, "r", round(r, 3))
his.rename("his_events").to_csv(f"{R}/sh_daily/his_daily_counts.csv")
fw.to_csv(f"{R}/sh_daily/ours_daily_forward.csv")
print("his total", int(his.sum()), "days with events", int((his > 0).sum()), "max/day", int(his.max()), his.idxmax().date())
print("ran days in 2018-2023-11:", len([d for d in build | test if d < pd.Timestamp('2023-12-01')]))
