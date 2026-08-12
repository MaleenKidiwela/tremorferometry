"""HDW dv/v restricted to the EHZ ('previous sensor', pre-2022-09 swap) era,
with southern vs northern candidate patches broken out. Uses the per-era CSV
(EHZ referenced to its own mean), so it's a clean single-sensor record.
"""
from __future__ import annotations

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HLAT, HLON = 47.649, -123.053
SWAP = "2022-09-01"
SMOOTH = 60

s = pd.read_csv("data/hdw_densify_set.summary.csv")
dlat = s.lat - HLAT
dlon = (s.lon - HLON) * np.cos(np.radians(HLAT))
s["az"] = (np.degrees(np.arctan2(dlon, dlat))) % 360
south = set(s[(s.az >= 110) & (s.az <= 250)].family_id)
north = set(s.family_id) - south

d = pd.read_csv("data/daily_dvv_HDW_coda_1to4_perera.csv")
d["date"] = pd.to_datetime(d["date"])
d = d[d.date < SWAP]                                  # EHZ era only


def med(df):
    m = df.groupby("date")["dvv"].median() * 100
    return m.rolling(SMOOTH, center=True, min_periods=SMOOTH // 3).median()


fig, ax = plt.subplots(figsize=(13, 5))
ax.axhline(0, color="0.6", lw=0.6)
allm = med(d)
sm = med(d[d.patch.isin(south)])
nm = med(d[d.patch.isin(north)])
ax.plot(allm.index, allm.values, color="k", lw=1.8, zorder=5,
        label=f"all 42 patches")
ax.plot(sm.index, sm.values, color="#d1495b", lw=1.4, zorder=4,
        label=f"southern patches (n={len(south)})")
ax.plot(nm.index, nm.values, color="#2e6f95", lw=1.4, zorder=4,
        label=f"northern patches (n={len(north)})")

# minor same-sensor response epochs in the EHZ era
for dt in ["2008-04-01", "2020-12-04"]:
    ax.axvline(pd.Timestamp(dt), color="0.55", lw=1.0, ls=":", zorder=2)
    ax.text(pd.Timestamp(dt), ax.get_ylim()[1], " EHZ epoch", rotation=90,
            va="top", ha="left", fontsize=6.5, color="0.4")
ax.axvline(pd.Timestamp("2001-02-28"), color="purple", lw=1.1, ls="--", zorder=2)
ax.text(pd.Timestamp("2001-02-28"), ax.get_ylim()[0], " Nisqually", rotation=90,
        va="bottom", ha="right", fontsize=7, color="purple")

ax.set_ylabel("dv/v (%)")
ax.set_xlabel("date")
ax.xaxis.set_major_locator(mdates.YearLocator(2))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.set_title("HDW coda dv/v (1-4 s) -- EHZ 'previous sensor' era only (pre-2022-09 swap), "
             f"per-era ref, {SMOOTH}-d median", fontsize=10)
ax.legend(loc="lower left", fontsize=8, ncol=3)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig("figures/smoke_dvv_HDW_EHZ_south_north.png", dpi=150)
print("wrote figures/smoke_dvv_HDW_EHZ_south_north.png")
print(f"span {d.date.min().date()} .. {d.date.max().date()} | "
      f"south patches {len(south)}, north {len(north)}")
