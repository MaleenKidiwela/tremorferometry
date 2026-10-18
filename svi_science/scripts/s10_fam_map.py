"""Map: one LFE family + the stations that recorded it (user request 2026-10-11)."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
R = "/home/jovyan/tremorferometry"
fam = pd.read_csv(f"{R}/data/catalogue/families_v1_3.csv")
fr = fam[(fam.region == "S_Vancouver_Is") & (fam.family_id.astype(str) == str(FAM))].iloc[0]
snr = {kv.split(":")[0]: float(kv.split(":")[1]) for kv in str(fr.clear_stations).split(";") if ":" in kv}
inv = pd.read_csv(f"{R}/svi_science/results/backward_scope/inventory_stations.csv").drop_duplicates("sta").set_index("sta")
polaris = ["KHVB", "MGCB", "SHDB", "SHVB", "TWBB", "KELB", "SILB", "SSIB", "TSJB", "TWGB", "TWKB", "LZB"]
other = ["B003", "B004", "B005", "B006", "B007", "LRIV"]
SPAN = {"PGC": "2003–2026", "VGZ": "2003–2026", "SNB": "2003–2021", "NLLB": "2003–2026", "B009": "2005–2026",
        "B010": "2005–2026", "B011": "2005–2026", "B926": "2007–2026", "B001": "2005–2026", "LZB": "2003–2014"}
fig = plt.figure(figsize=(12, 10)); ax = fig.add_axes([0.05, 0.06, 0.8, 0.88], projection=ccrs.PlateCarree())
ax.set_extent([-124.6, -122.8, 47.9, 49.35])
try:
    ax.add_feature(cf.LAND.with_scale("10m"), facecolor="#f2efe6"); ax.add_feature(cf.OCEAN.with_scale("10m"), facecolor="#dfeaf3")
    ax.add_feature(cf.COASTLINE.with_scale("10m"), lw=0.6); ax.add_feature(cf.BORDERS.with_scale("10m"), lw=0.5, ls="--")
except Exception as e:
    print("no coastlines:", e)
gl = ax.gridlines(draw_labels=True, lw=0.3, color="grey"); gl.top_labels = gl.right_labels = False
for s in other:
    if s in inv.index:
        r = inv.loc[s]; ax.plot(r.lon, r.lat, "v", ms=7, mfc="lightgrey", mec="grey", transform=ccrs.PlateCarree())
        ax.text(r.lon + 0.02, r.lat - 0.03, s, fontsize=7, color="grey", transform=ccrs.PlateCarree())
for s in polaris:
    if s in inv.index:
        r = inv.loc[s]; ax.plot(r.lon, r.lat, "^", ms=8, mfc="tab:green", mec="k", mew=0.5, transform=ccrs.PlateCarree())
        lab = "KELB/KLNB" if s == "KELB" else s
        ax.text(r.lon + 0.02, r.lat + 0.015, f"{lab}\n{SPAN.get(s, '2003–2006')}", fontsize=7, color="darkgreen", transform=ccrs.PlateCarree())
cm = plt.cm.plasma; vmax = max(snr.values())
for s, v in snr.items():
    if s not in inv.index: continue
    r = inv.loc[s]
    ax.plot(r.lon, r.lat, "s" if not s.startswith("B") else "o", ms=11, mfc=cm(v / vmax), mec="k", mew=0.8, transform=ccrs.PlateCarree(), zorder=5)
site = ["PGC", "B009", "B010", "B011"]
for s, v in snr.items():
    if s not in inv.index or s in site: continue
    r = inv.loc[s]
    ax.text(r.lon + 0.03, r.lat, f"{s}  SNR {v:.1f}  ({SPAN.get(s, '')})", fontsize=8, fontweight="bold", transform=ccrs.PlateCarree(), zorder=6)
r = inv.loc["PGC"]
txt = "PGC site (4 co-located):\n" + "\n".join(f"  {s}  SNR {snr[s]:.1f}  ({SPAN[s]})" for s in site if s in snr)
ax.annotate(txt, xy=(r.lon, r.lat), xytext=(-123.35, 48.93), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
            textcoords=ccrs.PlateCarree()._as_mpl_transform(ax), fontsize=8, fontweight="bold",
            bbox=dict(fc="white", ec="k", alpha=0.85), arrowprops=dict(arrowstyle="-", lw=0.8), zorder=8)
ax.plot(fr.lon, fr.lat, "*", ms=22, mfc="red", mec="k", transform=ccrs.PlateCarree(), zorder=7)
ax.text(fr.lon + 0.03, fr.lat - 0.06, f"family {FAM}\n{fr.lat:.2f}°N {abs(fr.lon):.2f}°W, {fr.depth_km:.0f} km deep", fontsize=9, color="red", fontweight="bold", transform=ccrs.PlateCarree(), zorder=7)
sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, vmax)); cax = fig.add_axes([0.88, 0.25, 0.02, 0.5]); cb = fig.colorbar(sm, cax=cax); cb.set_label("template SNR (clear stations)")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker="*", ls="", mfc="red", mec="k", ms=14, label="LFE family"),
                   Line2D([], [], marker="s", ls="", mfc="grey", mec="k", ms=9, label="Canadian station (clear, colour = SNR)"),
                   Line2D([], [], marker="o", ls="", mfc="grey", mec="k", ms=9, label="PB borehole (clear, colour = SNR)"),
                   Line2D([], [], marker="^", ls="", mfc="tab:green", mec="k", ms=9, label="POLARIS + LZB (era templates, 2003–2006)"),
                   Line2D([], [], marker="v", ls="", mfc="lightgrey", mec="grey", ms=8, label="other network stations (not clear)")],
          loc="lower left", fontsize=8, framealpha=0.9)
ax.set_title(f"SVI LFE family {FAM} and the stations that recorded it", fontsize=12)
out = f"{R}/svi_science/results/s10_dvv/pilot/figures/s10_map_fam{FAM}.png"; fig.savefig(out, dpi=120); print(out)
