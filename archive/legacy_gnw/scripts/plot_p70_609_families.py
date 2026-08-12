#!/usr/bin/env python
"""Plot locations of the 609 picker-LFE B011 families
(P(LFE)>=0.7 candidates, >=3 members, no year rule), colored by family-stack P(LFE)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd

STA_LAT, STA_LON = 48.65, -123.448  # B011 (co-located with PGC)

pk = pd.read_csv("data/family_picker_p70_m3nyr.csv")
lfe = pk[pk.pred == "LFE"].copy()
eq = pk[pk.pred == "EQ"].copy()

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(10.5, 10))
ax = plt.axes(projection=proj)
ax.set_extent([-125.4, -122.4, 47.5, 49.5], crs=proj)

ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#efece6", zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#dce6ee", zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, color="#7a7a7a", zorder=1)
gl = ax.gridlines(draw_labels=True, lw=0.3, color="#b8b8b8", alpha=0.6)
gl.top_labels = gl.right_labels = False

# 40 EQ multiplets, faint, for context
ax.scatter(eq.lon, eq.lat, s=22, marker="^", color="#c0392b", edgecolor="white",
           lw=0.3, alpha=0.55, zorder=3, label=f"picker EQ multiplet (n={len(eq)})",
           transform=proj)

# 609 picker-LFE families, colored by family-stack P(LFE)
sc = ax.scatter(lfe.lon, lfe.lat, c=lfe.p_lfe, cmap="viridis", vmin=0.5, vmax=1.0,
                s=26, edgecolor="white", lw=0.3, alpha=0.95, zorder=4,
                label=f"picker LFE (n={len(lfe)})", transform=proj)

ax.scatter([STA_LON], [STA_LAT], s=360, marker="*", color="#1a1a1a",
           edgecolor="white", lw=1.0, zorder=6, label="B011 station", transform=proj)

cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.02)
cb.set_label("family-stack P(LFE)", fontsize=10)

ax.set_title("B011 — 609 picker-LFE families\n"
             "P(LFE)≥0.7 candidates · ≥3 members · no year rule",
             fontsize=12.5, pad=12)
ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
plt.savefig("figures/b011_p70_609_family_locations.png", dpi=170, bbox_inches="tight")
print("saved figures/b011_p70_609_family_locations.png")
print(f"LFE={len(lfe)} EQ={len(eq)}")
