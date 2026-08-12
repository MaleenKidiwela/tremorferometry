#!/usr/bin/env python
"""Plot the 527 B011 families: P(LFE)>0.8, >=3 members, 2010-2017, no year rule.
Colored by family-stack P(LFE), sized by seed member count."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd

STA_LAT, STA_LON = 48.65, -123.448  # B011

pk = pd.read_csv("data/family_picker_p70_2010_2017_m3.csv")
sm = pd.read_csv("data/b011_disc_p70_2010_2017_m3.summary.csv").set_index("family_id")
pk = pk.join(sm[["n_members"]], on="fam")
hi = pk[pk.p_lfe > 0.8].copy()

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(10.5, 10))
ax = plt.axes(projection=proj)
ax.set_extent([-125.4, -122.4, 47.5, 49.5], crs=proj)
ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#efece6", zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#dce6ee", zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, color="#7a7a7a", zorder=1)
gl = ax.gridlines(draw_labels=True, lw=0.3, color="#b8b8b8", alpha=0.6)
gl.top_labels = gl.right_labels = False

sc = ax.scatter(hi.lon, hi.lat, c=hi.p_lfe, cmap="viridis", vmin=0.8, vmax=1.0,
                s=20 + 16 * (hi.n_members - 3), edgecolor="white", lw=0.3,
                alpha=0.9, zorder=4, transform=proj)
ax.scatter([STA_LON], [STA_LAT], s=380, marker="*", color="#1a1a1a",
           edgecolor="white", lw=1.0, zorder=6, label="B011 station", transform=proj)
for m in (3, 5, 7):
    ax.scatter([], [], s=20 + 16 * (m - 3), c="#3b7a57", edgecolor="white",
               lw=0.3, label=f"{m} seed members")
cb = plt.colorbar(sc, ax=ax, shrink=0.55, pad=0.02)
cb.set_label("family-stack P(LFE)", fontsize=10)
ax.set_title(f"B011 — {len(hi)} LFE families\n"
             "P(LFE)>0.8 · ≥3 members · 2010–2017 · no year rule",
             fontsize=12.5, pad=12)
ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
plt.savefig("figures/b011_p70_1017_m3_p80.png", dpi=170, bbox_inches="tight")
print(f"saved. n={len(hi)}")
