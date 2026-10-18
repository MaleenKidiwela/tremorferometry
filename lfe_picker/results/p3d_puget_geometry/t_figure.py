"""T-study step 3: figures/t_puget_geometry.png -- interface cells coloured by the number of 3-C
SITES within 45 km hypocentral in 2024, three panels (North reference / Puget-Olympic / Oregon),
with the stations that serve them overlaid.  Sequential single-hue ramp (count = magnitude)."""
import json
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.lines import Line2D
import cartopy.crs as ccrs, cartopy.feature as cfeature

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
FIG = f"{ROOT}/lfe_picker/figures/t_puget_geometry.png"
YEAR = 2024

C = pd.read_csv(f"{OUT}/cells.csv")
C = C[(C.year == YEAR) & (C.pool == "all")]
S = pd.read_csv(f"{OUT}/sites_all.csv")
FAM = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv")
alive = np.array([any(str(a) <= f"{YEAR}-07-01" <= str(b) for a, b in json.loads(i)) for i in S.intervals])
S = S[alive]

INK, MUTED, GRID = "#1c1c1c", "#5b5b5b", "#d9d9d9"
RAMP = LinearSegmentedColormap.from_list("cnt", ["#eef3fa", "#c6d9ef", "#94b6dd", "#5d8fc4",
                                                 "#2f6aa6", "#134e85", "#06305c"])
BNDS = [0, 1, 2, 3, 5, 8, 12, 40]
NORM = BoundaryNorm(BNDS, RAMP.N)

PANELS = [("North_48_49", "North 48-49 N  (reference)", "",
           (-125.4, -122.4, 47.75, 49.15)),
          ("Puget_47_48", "Puget / Olympic 47-48 N", "", (-124.4, -121.6, 46.75, 48.15)),
          ("Oregon_43.5_46", "Oregon 43.5-46 N", "", (-124.8, -121.6, 43.25, 46.15))]

fig = plt.figure(figsize=(15.0, 6.8))
TITLES = []
axes = []
for i, (reg, title, sub_t, ext) in enumerate(PANELS):
    ax = fig.add_subplot(1, 3, i + 1, projection=ccrs.PlateCarree())
    ax.set_extent(ext, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#fbfaf7", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#f2f5f7", edgecolor="none", zorder=0)
    ax.coastlines("50m", linewidth=0.5, color="#9a9a9a", zorder=1)
    g = C[C.region == reg]
    ax.scatter(g.lon, g.lat, c=g.n45, cmap=RAMP, norm=NORM, s=17, marker="s",
               linewidths=0, zorder=3, transform=ccrs.PlateCarree())
    ax.scatter(g.lon[~g.tremor_active], g.lat[~g.tremor_active], facecolors="none",
               edgecolors="#ffffff", s=17, marker="s", linewidths=0.45, zorder=4,
               transform=ccrs.PlateCarree())          # hollow = little/no tremor within 20 km
    sub = S[(S.lat.between(ext[2] - 0.6, ext[3] + 0.6)) & (S.lon.between(ext[0] - 0.7, ext[1] + 0.7))]
    bh = sub.kind == "borehole"
    ax.scatter(sub.lon[~bh], sub.lat[~bh], marker="^", s=34, facecolor="#ffffff",
               edgecolor=INK, linewidths=0.8, zorder=6, transform=ccrs.PlateCarree())
    ax.scatter(sub.lon[bh], sub.lat[bh], marker="s", s=26, facecolor="#ffffff",
               edgecolor="#b5651d", linewidths=0.9, zorder=6, transform=ccrs.PlateCarree())
    if reg == "North_48_49":
        ax.scatter(FAM.lon, FAM.lat, marker="o", s=13, facecolor="none", edgecolor="#b5651d",
                   linewidths=0.8, zorder=7, transform=ccrs.PlateCarree())
    med = int(np.median(g.n45[g.tremor_active]))
    TITLES.append((ax, title, f"median {med} sites $\\leq$45 km over tremor-active cells"
                   + (f"  --  {sub_t}" if sub_t else "")))
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color=GRID, alpha=0.9, linestyle=":")
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = dict(size=8, color=MUTED)
    axes.append(ax)

cax = fig.add_axes([0.33, 0.115, 0.34, 0.022])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=RAMP), cax=cax, orientation="horizontal",
                  spacing="uniform", ticks=BNDS)
cb.set_label(f"3-C sites within 45 km hypocentral of the Slab2 interface cell ({YEAR})",
             fontsize=9, color=INK)
cb.ax.xaxis.set_label_position("top")
cb.ax.tick_params(labelsize=8, color=MUTED, labelcolor=MUTED); cb.outline.set_visible(False)

leg = [Line2D([], [], marker="^", color="none", markerfacecolor="#fff", markeredgecolor=INK,
              markersize=7, label="3-C surface site (alive 2024)"),
       Line2D([], [], marker="s", color="none", markerfacecolor="#fff", markeredgecolor="#b5651d",
              markersize=6, label="3-C borehole site"),
       Line2D([], [], marker="o", color="none", markerfacecolor="#fff", markeredgecolor="#b5651d",
              markersize=6, label="130 Bostock LFE families (north panel)"),
       Line2D([], [], marker="s", color="none", markerfacecolor="#fff", markeredgecolor="#9a9a9a",
              markersize=7, label="cell with <100 tremor epicentres within 20 km (white outline)")]
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=2, frameon=False,
           fontsize=8.5, labelcolor=INK, columnspacing=3.0)
fig.suptitle("Hypocentral station density on the 28-45 km Cascadia interface, 2024 "
             "(sites = stations merged within 1 km)", fontsize=13, color=INK, y=0.985)
fig.subplots_adjust(left=0.045, right=0.985, top=0.83, bottom=0.20, wspace=0.18)
fig.canvas.draw()                      # equal-aspect axes shrink only at draw time
for ax, t, sub_t in TITLES:             # titles on one common baseline, left-aligned to each map
    x0 = ax.get_position().x0
    fig.text(x0, 0.905, t, fontsize=11.5, color=INK, ha="left", va="bottom", fontweight="bold")
    fig.text(x0, 0.875, sub_t, fontsize=9.5, color=MUTED, ha="left", va="bottom")
fig.savefig(FIG, dpi=170, facecolor="white")
print("wrote", FIG)
