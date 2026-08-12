#!/usr/bin/env python
"""Plot locations of the P>=0.7 / >=4-member / no-year-rule B011 families,
colored by Application-B picker class, with B011 station + Lin catalog LFEs overlaid."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd

STA_LAT, STA_LON = 48.65, -123.448  # B011 (co-located with PGC)

pk = pd.read_csv("data/family_picker_p70_m4nyr.csv")
sm = pd.read_csv("data/b011_disc_p70_m4nyr.summary.csv").set_index("family_id")
pk = pk.join(sm[["n_members"]], on="fam")
lfe = pk[pk.pred == "LFE"]
eq = pk[pk.pred == "EQ"]

lin = pd.read_csv("catalogs/lin_families_ets_2010_vi.csv")

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection=proj)
ax.set_extent([-125.2, -122.5, 47.6, 49.4], crs=proj)

ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#efece6", zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#dce6ee", zorder=0)
ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, color="#7a7a7a", zorder=1)
gl = ax.gridlines(draw_labels=True, lw=0.3, color="#b8b8b8", alpha=0.6)
gl.top_labels = gl.right_labels = False

# Lin catalog LFEs (ground-truth reference)
ax.scatter(lin.lon, lin.lat, s=42, marker="s", facecolor="none",
           edgecolor="#c8a020", lw=1.0, alpha=0.8, zorder=2,
           label=f"Lin (2023) catalog LFEs (n={len(lin)})", transform=proj)

# our families, sized by seed member count
ax.scatter(lfe.lon, lfe.lat, s=18 + 10 * (lfe.n_members - 4), marker="o",
           color="#2c7d3f", edgecolor="white", lw=0.4, alpha=0.9, zorder=4,
           label=f"picker LFE (n={len(lfe)})", transform=proj)
ax.scatter(eq.lon, eq.lat, s=18 + 10 * (eq.n_members - 4), marker="^",
           color="#c0392b", edgecolor="white", lw=0.4, alpha=0.9, zorder=4,
           label=f"picker EQ multiplet (n={len(eq)})", transform=proj)

ax.scatter([STA_LON], [STA_LAT], s=340, marker="*", color="#1a1a1a",
           edgecolor="white", lw=1.0, zorder=6, label="B011 station", transform=proj)

ax.set_title("B011 families  (P(LFE)≥0.7 candidates, ≥4 members, no year rule)\n"
             "colored by Application-B stack picker class",
             fontsize=12.5, pad=12)
ax.legend(loc="upper left", fontsize=9.5, framealpha=0.92)
plt.savefig("figures/b011_p70_m4_family_locations.png", dpi=170, bbox_inches="tight")
print("saved figures/b011_p70_m4_family_locations.png")
print(f"LFE={len(lfe)} EQ={len(eq)}  station B011 at {STA_LAT},{STA_LON}")
