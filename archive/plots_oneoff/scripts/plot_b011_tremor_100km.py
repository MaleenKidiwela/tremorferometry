#!/usr/bin/env python
"""Map all PNSN tremor within 100 km of B011 (master catalog) to verify coverage
is not clipped. Overlays station, 100 km ring, and the 609 picker-LFE families."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd

STA_LAT, STA_LON = 48.65, -123.448  # B011
R_KM = 100.0

t = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv",
                usecols=["time", "lat", "lon"])
dkm = np.sqrt(((t.lat - STA_LAT) * 111.0) ** 2 +
              ((t.lon - STA_LON) * 111.0 * np.cos(np.radians(STA_LAT))) ** 2)
near = t[dkm <= R_KM].copy()
near["year"] = pd.to_datetime(near.time).dt.year
print(f"tremor within {R_KM:.0f} km of B011: {len(near):,}")
print("lat span:", round(near.lat.min(), 2), "..", round(near.lat.max(), 2))
print("year span:", int(near.year.min()), "..", int(near.year.max()),
      "| distinct years:", near.year.nunique())

# 100 km ring
th = np.linspace(0, 2 * np.pi, 400)
ring_lat = STA_LAT + (R_KM / 111.0) * np.sin(th)
ring_lon = STA_LON + (R_KM / (111.0 * np.cos(np.radians(STA_LAT)))) * np.cos(th)

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(11, 10))
ax = plt.axes(projection=proj)
ax.set_extent([-125.6, -121.4, 47.4, 49.9], crs=proj)
ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#efece6", zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#dce6ee", zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, color="#7a7a7a", zorder=1)
gl = ax.gridlines(draw_labels=True, lw=0.3, color="#b8b8b8", alpha=0.6)
gl.top_labels = gl.right_labels = False

ax.scatter(near.lon, near.lat, s=3, c="#3a6ea5", alpha=0.18, lw=0, zorder=2,
           label=f"PNSN tremor <100 km (n={len(near):,})", transform=proj)
ax.plot(ring_lon, ring_lat, color="#c0392b", lw=1.6, ls="--", zorder=3,
        label="100 km radius", transform=proj)

# 609 picker-LFE families for comparison
try:
    pk = pd.read_csv("data/family_picker_p70_m3nyr.csv")
    lfe = pk[pk.pred == "LFE"]
    ax.scatter(lfe.lon, lfe.lat, s=14, c="#2c7d3f", edgecolor="white", lw=0.25,
               alpha=0.9, zorder=4, label=f"609 picker-LFE families", transform=proj)
except Exception as e:
    print("families overlay skipped:", e)

ax.scatter([STA_LON], [STA_LAT], s=360, marker="*", color="#1a1a1a",
           edgecolor="white", lw=1.0, zorder=6, label="B011 station", transform=proj)

ax.set_title(f"PNSN tremor within 100 km of B011 (master catalog)\n"
             f"coverage check — n={len(near):,}, lat {near.lat.min():.2f}–{near.lat.max():.2f}°N",
             fontsize=12.5, pad=12)
ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
plt.savefig("figures/b011_tremor_100km.png", dpi=170, bbox_inches="tight")
print("saved figures/b011_tremor_100km.png")
