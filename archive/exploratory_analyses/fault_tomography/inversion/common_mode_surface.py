"""Monthly SHALLOW (common-mode) dv/v SURFACE map.
Each station's COMMON MODE = mean dv/v across ITS families (the near-receiver, shallow signal shared by all
families at a site). Per month, interpolate the per-station common-mode ANOMALY across space -> a surface
dv/v field (the complement of the deep residual inversion). Era-split + per-era 2019-24 baseline demean so
instrument/metadata steps don't fake surface features. Trailing WIN-month mean per station; grid masked
beyond RAD km of the nearest station (no wild extrapolation).
Outputs {RESDIR}/common_mode_surface_movie.mp4 + common_mode_surface_by_year.png.
Reads {RESDIR}/{pairs.csv, pair_months.parquet, era_table.csv}."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cartopy.crs as ccrs, cartopy.feature as cfeature
RC = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
DESEASON = os.environ.get("DESEASON", "0") == "1"          # default: keep seasonal (it's a real shallow signal)
WIN = int(os.environ.get("WIN", "3")); STEP = 0.1; RAD = float(os.environ.get("RAD", "45"))  # mask radius (km)

pairs = pd.read_csv(f"{RC}/pairs.csv")
sta = pairs.drop_duplicates("tag")[["tag", "sta_lat", "sta_lon"]].reset_index(drop=True)
pm = pd.read_parquet(f"{RC}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["t"] = pd.PeriodIndex(pm.ym, freq="M")
# station COMMON MODE per month = n_days-weighted mean of dvv over the station's families (cells)
cm = pm.assign(wd=pm.dvv_month*pm.n_days).groupby(["tag", "t"]).agg(num=("wd", "sum"), den=("n_days", "sum")).reset_index()
cm["cm"] = cm.num/cm.den; cm["year"] = cm.t.dt.year
# instrument era-split
et = pd.read_csv(f"{RC}/era_table.csv"); bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
cm["eid"] = [int(sum(1 for b in bnd_of.get(tg, []) if pd.Timestamp(tt.to_timestamp()) >= b)) for tg, tt in zip(cm.tag, cm.t)]
# optional conditional deseason per (station, era)
def harm(x): return np.column_stack([np.sin(2*np.pi*x/12), np.cos(2*np.pi*x/12), np.sin(4*np.pi*x/12), np.cos(4*np.pi*x/12)])
def deseason(g):
    v = g.cm.values.astype(float); x = (g.t-g.t.min()).apply(lambda z: z.n).values.astype(float); n = len(v)
    if DESEASON and n >= 24:
        X = np.column_stack([np.ones_like(x), harm(x)]); b, *_ = np.linalg.lstsq(X, v, rcond=None)
        RSSf = ((v-X@b)**2).sum(); RSSc = ((v-v.mean())**2).sum(); F = ((RSSc-RSSf)/4)/(RSSf/(n-5)) if RSSf > 0 else 0
        if 1-fdist.cdf(F, 4, n-5) < 0.05: v = v-harm(x)@b[1:]
    return v
cm = cm.sort_values(["tag", "eid", "t"])
cm["cmd"] = np.concatenate([deseason(g) for _, g in cm.groupby(["tag", "eid"])])
# per (station, era) demean to 2019-24 baseline (instrument-step safe)
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); base = g.cmd.values[inw].mean() if inw.sum() >= 6 else g.cmd.values.mean()
    return g.cmd.values - base
cm["anom"] = np.concatenate([demean(g) for _, g in cm.groupby(["tag", "eid"])])

# trailing WIN-month mean per station
months = pd.period_range("2011-01", "2026-06", freq="M")
piv = cm.pivot_table(index="tag", columns="t", values="anom", aggfunc="mean").reindex(columns=months)
roll = piv.T.rolling(WIN, min_periods=1).mean().T
stll = sta.set_index("tag").reindex(piv.index); slat = stll.sta_lat.values; slon = stll.sta_lon.values

# grid + distance mask
glon = np.arange(slon.min()-0.3, slon.max()+0.3+1e-9, STEP); glat = np.arange(slat.min()-0.3, slat.max()+0.3+1e-9, STEP)
GLO, GLA = np.meshgrid(glon, glat)
lat0, lon0 = slat.mean(), slon.mean()
def xy(la, lo): return np.column_stack([(np.asarray(lo)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(la)-lat0)*111])
gxy = xy(GLA.ravel(), GLO.ravel())
def surface(vals, la, lo):
    Z = griddata(np.column_stack([lo, la]), vals, (GLO, GLA), method="linear")
    dmin, _ = cKDTree(xy(la, lo)).query(gxy); Z[(dmin.reshape(GLO.shape) > RAD)] = np.nan
    return Z
def frame(k):
    col = roll.values[:, k]; ok = np.isfinite(col)
    if ok.sum() < 8: return None, 0
    return surface(col[ok], slat[ok], slon[ok]), int(ok.sum())
surfs = [frame(k) for k in range(len(months))]
nst = [s[1] for s in surfs]; allv = np.concatenate([s[0][np.isfinite(s[0])] for s in surfs if s[0] is not None])
VL = float(np.nanpercentile(np.abs(allv), 97)); VL = min(max(VL, 0.1), 0.6)
print(f"months {len(months)} | stations/frame median {int(np.median([n for n in nst if n]))} | vlim ±{VL:.2f}%")

ext = [slon.min()-0.5, slon.max()+0.5, slat.min()-0.4, slat.max()+0.4]; proj = ccrs.PlateCarree()
# ---- monthly movie ----
fig = plt.figure(figsize=(6, 9)); ax = plt.axes(projection=proj); ax.set_extent(ext, crs=proj)
ax.add_feature(cfeature.LAND, facecolor="#efefef"); ax.add_feature(cfeature.OCEAN, facecolor="#dbeafe"); ax.coastlines(resolution="50m", lw=.6)
Z0 = next(s[0] for s in surfs if s[0] is not None)
pcm = ax.pcolormesh(glon, glat, np.ma.masked_invalid(Z0), cmap="RdBu_r", vmin=-VL, vmax=VL, transform=proj, shading="auto")
scat = ax.scatter(slon, slat, s=8, c="k", marker="^", transform=proj, zorder=5)
plt.colorbar(pcm, ax=ax, fraction=.04, label="shallow common-mode dv/v anomaly (%)")
ttl = ax.set_title("", fontsize=12)
def upd(k):
    Z = surfs[k][0]
    pcm.set_array(np.ma.masked_invalid(Z).ravel() if Z is not None else np.ma.masked_all(GLO.size))
    ttl.set_text(f"shallow common-mode dv/v  {months[k]}  ({surfs[k][1]} stns, {WIN}-mo)")
    return pcm, ttl
FuncAnimation(fig, upd, frames=len(months), blit=False).save(f"{RC}/common_mode_surface_movie.mp4", writer="ffmpeg", fps=6, dpi=110)
print(f"wrote {RC}/common_mode_surface_movie.mp4")
plt.close(fig)
# ---- annual-mean sample maps ----
yrs = list(range(2012, 2026)); yidx = {y: [k for k, m in enumerate(months) if m.year == y] for y in yrs}
fig2, axs = plt.subplots(2, 7, figsize=(22, 8), subplot_kw={"projection": proj})
for j, y in enumerate(yrs):
    a = axs.flat[j]; ks = [k for k in yidx[y] if surfs[k][0] is not None]
    a.set_extent(ext, crs=proj); a.add_feature(cfeature.LAND, facecolor="#efefef"); a.add_feature(cfeature.OCEAN, facecolor="#dbeafe"); a.coastlines(resolution="50m", lw=.5)
    if ks:
        Zy = np.nanmean([surfs[k][0] for k in ks], axis=0)
        pcm = a.pcolormesh(glon, glat, np.ma.masked_invalid(Zy), cmap="RdBu_r", vmin=-VL, vmax=VL, transform=proj, shading="auto")
    a.scatter(slon, slat, s=3, c="k", marker="^", transform=proj, zorder=5); a.set_title(f"{y}", fontsize=10)
fig2.suptitle(f"Shallow common-mode dv/v SURFACE (annual mean, interpolated from per-station common mode, ±{RAD:.0f} km mask, mirror-free)", fontsize=13, y=.99)
cax = fig2.add_axes([1.005, 0.3, 0.008, 0.4]); plt.colorbar(pcm, cax=cax, label="dv/v anomaly (%)")
fig2.tight_layout(); fig2.savefig(f"{RC}/common_mode_surface_by_year.png", dpi=120, bbox_inches="tight")
print(f"wrote {RC}/common_mode_surface_by_year.png")
