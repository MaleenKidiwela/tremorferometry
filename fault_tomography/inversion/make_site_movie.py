"""Movie of the NEAR-SITE term from the daily joint inversion: SITES (true near-station dv/v,
one value per station per day — the ground under each station: shallow crust, weather, season),
drawn at station locations. Companion to the near-source (interface) movie. Full record (site
terms need no transmission correction — they are direct per-station estimates; the 2009-10
exclusion applies to the transmission-corrected INTERFACE product only)."""
import os, re
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cartopy.crs as ccrs, cartopy.feature as cfeature

RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
STRIDE = int(os.environ.get("STRIDE", "7")); FPS = int(os.environ.get("FPS", "12"))
VLIM = float(os.environ.get("VLIM", "1.0"))

z = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
S = np.asarray(z["SITES"], float)                    # (T, NS) true near-station dv/v (%)
days = pd.to_datetime(z["days"]); tags = [str(x) for x in z["sta_tags"]]
base = [re.sub(r"p\d+f\d+$", "", t) for t in tags]
coord = {r.sta: (r.lat, r.lon) for _, r in pd.read_csv("data/candidate_stations_post2020.csv").iterrows()}
have = [i for i, b in enumerate(base) if b in coord]
miss = sorted(set(b for b in base if b not in coord))
if miss:
    print(f"WARNING: {len(miss)} stations without coords, dropped: {miss[:10]}")
slat = np.array([coord[base[i]][0] for i in have]); slon = np.array([coord[base[i]][1] for i in have])
S = S[:, have]
ti = np.arange(0, len(days), STRIDE)

med = pd.Series(np.nanmedian(S, axis=1), index=days).rolling("90D", min_periods=30).mean()

proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(7.2, 9.6))
axm = fig.add_axes([0.06, 0.24, 0.88, 0.72], projection=proj)
axt = fig.add_axes([0.10, 0.05, 0.84, 0.13])
axm.set_extent([slon.min() - 0.7, slon.max() + 0.7, slat.min() - 0.5, slat.max() + 0.5], crs=proj)
axm.add_feature(cfeature.LAND, facecolor="#f0f0f0")
axm.add_feature(cfeature.OCEAN, facecolor="#dbeafe")
axm.coastlines(resolution="50m", lw=0.6)
cmap = plt.cm.RdBu_r.copy()
cmap.set_bad((0.85, 0.85, 0.85, 0.35))               # inactive station = translucent gray
sc = axm.scatter(slon, slat, c=np.zeros(len(slon)), s=42, marker="o", cmap=cmap,
                 vmin=-VLIM, vmax=VLIM, transform=proj, edgecolors="#00000033", lw=0.4, zorder=3)
cb = fig.colorbar(sc, ax=axm, shrink=0.6, pad=0.02)
cb.set_label("near-site dv/v (%, site term, rel. era gauge)")
ttl = axm.set_title("", fontsize=11)
axt.plot(med.index, med.values, lw=1.0, color="#334155")
axt.axhline(0, lw=0.5, color="#94a3b8")
axt.set_ylabel("median site term (%)", fontsize=8); axt.tick_params(labelsize=7)
lo, hi = np.nanpercentile(med, [1, 99]); pad = 0.2 * (hi - lo + 1e-9)
axt.set_ylim(lo - pad, hi + pad)
cursor = axt.axvline(days[0], color="#dc2626", lw=1.2)

def update(k):
    t = ti[k]
    v = S[t]
    good = np.isfinite(v)
    sc.set_array(np.ma.masked_invalid(v))            # array drives colors; masked -> bad color
    cursor.set_xdata([days[t], days[t]])
    ttl.set_text(f"Near-SITE term (station ground) — {days[t].date()}  ({good.sum()} stations)")
    return sc, cursor, ttl

anim = FuncAnimation(fig, update, frames=len(ti), blit=False)
out = f"{RD}/inversion_4d_daily_movie_sites.mp4"
anim.save(out, writer=FFMpegWriter(fps=FPS, bitrate=2500), dpi=110)
print(f"wrote {out} | {len(ti)} frames @ {FPS} fps | {len(have)} stations")
