import warnings, glob, os; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cfeature
def famll(pid):
    a = pid.split("__")[0].split("_"); return float(a[0]), float(a[1])
# all reliable (certified) families across all stations
fam = set()
for fp in glob.glob("data/*_causality_cert.csv"):
    try:
        c = pd.read_csv(fp); [fam.add(p) for p in c[c.reliable].fam]
    except Exception: pass
for fp in glob.glob("data/*_fwd_vs_rev_coda.csv"):
    stem = os.path.basename(fp).replace("_fwd_vs_rev_coda.csv", "")
    if not os.path.exists(f"data/{stem}_causality_cert.csv"):
        try:
            c = pd.read_csv(fp); [fam.add(p) for p in c[c.ratio > 1.5].fam]
        except Exception: pass
flat = np.array([famll(p) for p in fam]); flat, flon = flat[:, 0], flat[:, 1]
# station coords (from the inversion catalog)
pr = pd.read_csv("fault_tomography/inversion/res_catalog/pairs.csv")
st = pr.drop_duplicates("tag")[["sta_lat", "sta_lon"]].values
proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(9, 12)); a = plt.axes(projection=proj)
a.set_extent([-127.5, -121, 39.5, 50.5], crs=proj)
a.add_feature(cfeature.LAND, facecolor="#f2f2f2"); a.add_feature(cfeature.OCEAN, facecolor="#dbeafe")
a.coastlines(resolution="50m", lw=.7); a.add_feature(cfeature.BORDERS, lw=.4); a.add_feature(cfeature.STATES, lw=.3)
a.scatter(flon, flat, s=6, c="#c0392b", alpha=.35, transform=proj, label=f"certified LFE families (n={len(fam)})", edgecolors="none")
a.scatter(st[:, 1], st[:, 0], s=55, c="#1a4d8f", marker="^", edgecolors="k", lw=.5, transform=proj, label=f"stations (n={len(st)})", zorder=5)
gl = a.gridlines(draw_labels=True, lw=.2); gl.top_labels = gl.right_labels = False
a.legend(loc="upper right", fontsize=11); a.set_title(f"Cascadia LFE-coda dv/v network: {len(fam)} certified families + {len(st)} stations", fontsize=12)
fig.savefig("fault_tomography/inversion/res_catalog/families_stations_map.png", dpi=130, bbox_inches="tight")
print(f"wrote families_stations_map.png | families {len(fam)} | stations {len(st)}")
