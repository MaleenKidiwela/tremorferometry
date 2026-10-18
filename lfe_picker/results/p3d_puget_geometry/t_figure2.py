"""T-study (revised brief) step 4: figures/t_stations_2019plus_map.png
Every 3-C seismometer station with data 2019+ in the box, coloured by sensor type, roster vs
non-roster marked, over PNSN 2019+ tremor density, with Slab2 30/40 km contours and the 130
Bostock families.  Three stacked panels (phone-legible): whole box, north zoom, Olympics/Puget zoom.
Categorical palette #2a78d6 / #4a3aa7 / #eda100 (validated: lightness, chroma, CVD, normal-vision;
contrast WARN on yellow relieved by marker shape + legend)."""
import json
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.lines import Line2D
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from scipy.ndimage import gaussian_filter
import cartopy.crs as ccrs, cartopy.feature as cfeature

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_picker/results/p3d_puget_geometry"
FIG = f"{ROOT}/lfe_picker/figures/t_stations_2019plus_map.png"
INK, MUTED, GRID = "#1c1c1c", "#56554f", "#d9d9d9"
COL = {"surface_broadband": "#2a78d6", "borehole": "#4a3aa7", "surface_short_period": "#eda100"}
MRK = {"surface_broadband": "^", "borehole": "s", "surface_short_period": "D"}
FAMC = "#eb6834"
TRC = LinearSegmentedColormap.from_list("tr", ["#ffffff", "#e8e6f2", "#cfcada", "#a9a2bd"])

ST = pd.read_csv(f"{OUT}/stations_all_2019plus.csv")
ST = ST[ST.detect_capable].reset_index(drop=True)
ST["alive25"] = [any(a <= "2025-07-01" <= b for a, b in json.loads(i)) for i in ST.intervals]
SI = pd.read_csv(f"{OUT}/sites.csv")
RK = pd.read_csv(f"{OUT}/station_ranking.csv")
FAM = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv")
TR = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
TR = TR[(TR.time >= "2019-01-01") & TR.lat.between(46.3, 49.7) & TR.lon.between(-125.7, -121.3)]
alive25 = np.array([any(s <= "2025-07-01" <= e for s, e in json.loads(i)) for i in SI.intervals])
SI["alive25"] = alive25

sl = pd.read_csv(f"{ROOT}/data/cas_slab2_input_04-18.csv", usecols=["lat", "lon", "depth", "etype"])
sl = sl[~sl.etype.isin(["BA", "TO"])]
sl = sl[(sl.depth >= 10) & (sl.depth <= 70)].dropna(subset=["lat", "lon", "depth"])
_L = LinearNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
_N = NearestNDInterpolator(sl[["lon", "lat"]].values, sl.depth.values)
gx, gy = np.meshgrid(np.linspace(-125.7, -121.3, 320), np.linspace(46.3, 49.7, 320))
gz = _L(np.column_stack([gx.ravel(), gy.ravel()]))
valid = np.isfinite(gz).reshape(gx.shape)      # inside the Slab2 control-point hull
gz = np.where(np.isfinite(gz), gz, _N(np.column_stack([gx.ravel(), gy.ravel()]))).reshape(gx.shape)
gz = gaussian_filter(gz, 9.0)           # the Slab2 INPUT file is a scattered control-point cloud;
                                        # smooth before contouring (point queries stay unsmoothed)
gz = np.ma.array(gz, mask=~valid)       # never contour the nearest-neighbour extrapolation

BA = pd.read_csv(f"{OUT}/bostock_array_stations.csv")     # the 21 template stations, 2003-2005
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]
TOPN = ["JRDRB", "LHLYB", "DAOB", "ALBHB", "PTRF", "ODUC", "ODEER"]
TOPP = list(RK[RK.greedy_rank_puget > 0].sort_values("greedy_rank_puget").site.head(10))
BASE = list(SI[SI.is_base_plus].site)

PANELS = [("a", "Whole box: every 3-C seismometer with data 2019+", (-125.5, -121.5, 46.5, 49.5), []),
          ("b", "North / Salish Sea: the 130 Bostock LFE families", (-124.9, -122.7, 47.9, 49.3), TOPN + ["SYMB", "BUCKS"]),
          ("c", "Olympics / Puget: tremor-weighted interface", (-124.7, -121.9, 46.6, 48.3), TOPP)]

fig = plt.figure(figsize=(8.6, 16.4))
for i, (tag, title, ext, labels) in enumerate(PANELS):
    ax = fig.add_subplot(3, 1, i + 1, projection=ccrs.PlateCarree())
    ax.set_extent(ext, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale("50m"), facecolor="#fbfaf7", edgecolor="none", zorder=0)
    ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="#f1f4f6", edgecolor="none", zorder=0)
    H, xe, ye = np.histogram2d(TR.lon, TR.lat, bins=[np.arange(ext[0], ext[1] + .02, .02),
                                                     np.arange(ext[2], ext[3] + .02, .02)])
    ax.pcolormesh(xe, ye, H.T, cmap=TRC, norm=PowerNorm(0.45, vmin=0, vmax=max(H.max(), 1)),
                  zorder=1, shading="flat", transform=ccrs.PlateCarree())
    cs = ax.contour(gx, gy, gz, levels=[30, 40], colors=["#3f5d78"], linewidths=1.2,
                    linestyles=["--", "-"], zorder=4, transform=ccrs.PlateCarree())
    ax.clabel(cs, fmt="%d km", fontsize=8.5, colors="#3f5d78")
    ax.coastlines("50m", linewidth=0.5, color="#b9b9b9", zorder=3)
    ax.scatter(FAM.lon, FAM.lat, s=17, marker="o", facecolor="none", edgecolor=FAMC,
               linewidths=0.9, zorder=5, transform=ccrs.PlateCarree())
    for t, g in ST.groupby("sensor_type"):
        for alv, (fc, ec, lw, sz) in [(True, (COL[t], INK, 0.7, 46)), (False, ("white", COL[t], 1.1, 34))]:
            gg = g[g.alive25 == alv]
            ax.scatter(gg.lon, gg.lat, s=sz, marker=MRK[t], facecolor=fc, edgecolor=ec,
                       linewidths=lw, zorder=6, transform=ccrs.PlateCarree())
            r = gg[gg.in_roster183]
            ax.scatter(r.lon, r.lat, s=sz * 3.1, marker=MRK[t], facecolor="none",
                       edgecolor="#4a4a4a", linewidths=0.8, zorder=7, transform=ccrs.PlateCarree())
    if tag == "b":                       # Bostock's own 2003-2005 template array
        tmp, per = BA[~BA.permanent], BA[BA.permanent]
        ax.scatter(tmp.lon, tmp.lat, s=54, marker="v", facecolor=FAMC, edgecolor="#7a2f10",
                   linewidths=0.7, zorder=7, transform=ccrs.PlateCarree())
        ax.scatter(per.lon, per.lat, s=64, marker="v", facecolor="none", edgecolor="#7a2f10",
                   linewidths=1.3, zorder=7, transform=ccrs.PlateCarree())
        for k, r in enumerate(BA.sort_values(["lat", "lon"]).itertuples()):
            dy, ha = ((-12, "center") if k % 3 == 0 else
                      (8, "center") if k % 3 == 1 else (-3, "right"))
            dx = 0 if ha == "center" else -7
            ax.annotate(r.sta, (r.lon, r.lat), xytext=(dx, dy), textcoords="offset points",
                        fontsize=6.8, color="#7a2f10", ha=ha, zorder=9,
                        path_effects=HALO, transform=ccrs.PlateCarree())
    for s in labels:
        r = SI[SI.site == s]
        if len(r) and ext[0] < float(r.lon.iloc[0]) < ext[1] and ext[2] < float(r.lat.iloc[0]) < ext[3]:
            dy = -11 if s in ("LHLYB", "ALBHB", "OHOH", "DOSE", "WINDI", "PAN4H") else 5
            ha, dx = ("right", -8) if s in ("ALBHB", "DAOB", "ODEER") else ("left", 7)
            ax.annotate(s, (float(r.lon.iloc[0]), float(r.lat.iloc[0])), xytext=(dx, dy),
                        textcoords="offset points", fontsize=9.0, color=INK, zorder=10, ha=ha,
                        fontweight="bold", path_effects=HALO, transform=ccrs.PlateCarree())
    ax.text(0.012, 0.965, f"({tag}) {title}", transform=ax.transAxes, fontsize=11.5, color=INK,
            va="top", ha="left", fontweight="bold",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=3.0), zorder=9)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color=GRID, linestyle=":", alpha=0.9)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = dict(size=9, color=MUTED)

leg = [Line2D([], [], marker="^", ls="none", mfc=COL["surface_broadband"], mec=INK, ms=8,
              label="surface broadband"),
       Line2D([], [], marker="s", ls="none", mfc=COL["borehole"], mec=INK, ms=7, label="borehole"),
       Line2D([], [], marker="D", ls="none", mfc=COL["surface_short_period"], mec=INK, ms=7,
              label="surface short-period"),
       Line2D([], [], marker="^", ls="none", mfc="white", mec=COL["surface_broadband"], ms=7,
              label="not operating in 2025 (hollow)"),
       Line2D([], [], marker="^", ls="none", mfc="none", mec="#4a4a4a", ms=13,
              markeredgewidth=0.8, label="outer ring = in the 183-station roster"),
       Line2D([], [], marker="o", ls="none", mfc="none", mec=FAMC, ms=7,
              label="130 Bostock LFE families"),
       Line2D([], [], marker="v", ls="none", mfc=FAMC, mec="#7a2f10", ms=8,
              label="Bostock 2003-05 array: 15 POLARIS temporaries"),
       Line2D([], [], marker="v", ls="none", mfc="none", mec="#7a2f10", ms=8, markeredgewidth=1.3,
              label="Bostock 2003-05 array: 6 permanent CN"),
       Line2D([], [], ls="--", color="#3f5d78", label="Slab2 30 km (dashed) / 40 km (solid)"),
       Line2D([], [], marker="s", ls="none", mfc="#a9a2bd", mec="none", ms=9,
              label="PNSN tremor density 2019-2026")]
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.004), ncol=2, frameon=False,
           fontsize=9.2, labelcolor=INK, columnspacing=2.0, handletextpad=0.7)
fig.suptitle("Seismic stations available 2019-2026, Cascadia 46.5-49.5 N\n"
             "3-C seismometers only (strong-motion and Z-only listed separately)",
             fontsize=13.5, color=INK, y=0.995)
fig.subplots_adjust(left=0.085, right=0.985, top=0.955, bottom=0.088, hspace=0.10)
fig.savefig(FIG, dpi=150, facecolor="white")
print("wrote", FIG)
