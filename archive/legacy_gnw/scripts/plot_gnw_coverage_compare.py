"""Map GNW: already-done families vs coverage-optimal candidates, highlighting
the coverage GAPS (candidates not yet densified). Shows whether the 75 done give
good azimuth/distance coverage within +/-100 km."""
from __future__ import annotations

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GLAT, GLON = 47.5641, -122.8250
fam = pd.read_csv("data/gnw_pnsn_families.summary.csv").set_index("family_id")
done = (set(pd.read_csv("data/mf_gnw_all.csv", usecols=["template"]).template.unique())
        | set(pd.read_csv("data/mf_gnwcircle_all.csv", usecols=["template"]).template.unique()))
cand = pd.read_csv("data/gnw_coverage_candidates.summary.csv")
done_ids = [f for f in done if f in fam.index]
done_xy = fam.loc[done_ids]
new = cand[~cand.family_id.isin(done)]          # coverage gaps
hit = cand[cand.family_id.isin(done)]           # candidates already done

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(10, 9))
ax = plt.axes(projection=proj)
lons = list(done_xy.lon) + list(cand.lon) + [GLON]
lats = list(done_xy.lat) + list(cand.lat) + [GLAT]
pad = 0.25
ax.set_extent([min(lons) - pad, max(lons) + pad, min(lats) - pad, max(lats) + pad], crs=proj)
ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#cfe2f3")
ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f5ead6")
ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#cfe2f3")
ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, edgecolor="#555")
ax.add_feature(cfeature.BORDERS.with_scale("10m"), lw=0.5, edgecolor="#888", linestyle=":")

# +/-100 km circle
th = np.linspace(0, 2 * np.pi, 200)
ax.plot(GLON + (100 / (111 * np.cos(np.radians(GLAT)))) * np.cos(th),
        GLAT + (100 / 111.0) * np.sin(th), color="red", lw=1.0, ls="--",
        transform=proj, zorder=4, label="100 km")

ax.scatter(done_xy.lon, done_xy.lat, s=22, c="0.55", marker="o", edgecolor="none",
           zorder=5, transform=proj, label=f"already done ({len(done_ids)})")
ax.scatter(hit.lon, hit.lat, s=70, facecolor="none", edgecolor="green", lw=1.4,
           marker="D", zorder=6, transform=proj, label=f"coverage cand. done ({len(hit)})")
sc = ax.scatter(new.lon, new.lat, c=new.snr, cmap="autumn_r", vmin=10, vmax=18,
                s=90, marker="D", edgecolor="k", lw=0.5, zorder=7, transform=proj,
                label=f"coverage GAP, not done ({len(new)})")
fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.02).set_label("SNR (gap candidates)")
ax.plot(GLON, GLAT, marker="^", color="red", markersize=15, markeredgecolor="k",
        linestyle="none", zorder=8, transform=proj, label="GNW")

gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray", alpha=0.4, linestyle="--")
gl.top_labels = gl.right_labels = False
ax.set_title("GNW: done families vs coverage-optimal candidates (red diamonds = "
             "coverage gaps not yet densified)", fontsize=9.5)
ax.legend(loc="upper left", fontsize=7.5, framealpha=0.9)
fig.savefig("figures/smoke_gnw_coverage_compare.png", dpi=150, bbox_inches="tight")
print("wrote figures/smoke_gnw_coverage_compare.png")
