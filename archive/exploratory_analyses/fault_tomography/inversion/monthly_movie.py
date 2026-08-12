#!/usr/bin/env python
"""Monthly deep delta-beta/beta MOVIE on the 0.2 deg grid (common-mode removed, mirror-free). Inverts a TRAILING
3-month window per frame (monthly cadence, enough data/frame). Saves an animated gif of the deep-cell maps."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cartopy.crs as ccrs, cartopy.feature as cfeature
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20"); LAM = float(os.environ.get("LAM", "0.5")); STEP = float(os.environ.get("STEP", "0.2"))
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); clat = d["cell_lat"]; clon = d["cell_lon"]; depth = d["depth_km"]
cxy = d["cxy"]; ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; deep = depth > 30   # family-median depth
_, nbr = cKDTree(cxy).query(cxy, k=min(6, Mf)); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]: L[i, i] += 1; L[i, j] -= 1
LtL = L.T@L
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["pair"] = pm.tag+"|"+pm.cell; pm["year"] = pm.t.dt.year
et = pd.read_csv(f"{RD}/era_table.csv"); bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
ts = pm.t.dt.to_timestamp().values; eid = np.zeros(len(pm), int); strad = np.zeros(len(pm), bool)
for i, (tg, t) in enumerate(zip(pm.tag.values, ts)):
    bs = bnd_of.get(tg, [])
    if bs: eid[i] = int(sum(1 for b in bs if t >= np.datetime64(b))); strad[i] = any(abs((pd.Timestamp(t)-b).days) <= 35 for b in bs)
pm["stera"] = pm.tag+"__e"+eid.astype(str); pm = pm[~strad].copy()
def harm(t): yr = 12.; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def deseason(g):
    v = g.dvv_month.values.astype(float); t = (g.t-g.t.min()).apply(lambda x: x.n).values.astype(float); n = len(v)
    if n >= 24:
        X = np.column_stack([np.ones_like(t), harm(t)]); b, *_ = np.linalg.lstsq(X, v, rcond=None)
        RSSf = ((v-X@b)**2).sum(); RSSc = ((v-v.mean())**2).sum(); F = ((RSSc-RSSf)/4)/(RSSf/(n-5)) if RSSf > 0 else 0
        if 1-fdist.cdf(F, 4, n-5) < 0.05: v = v-harm(t)@b[1:]
    return v
pm = pm.sort_values(["pair", "t"]); pm["ds"] = np.concatenate([deseason(g) for _, g in pm.groupby("pair")])
pm["cmode"] = pm.groupby(["stera", "ym"]).ds.transform("mean"); pm["resid"] = pm.ds - pm.cmode
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); return g.resid.values - (g.resid.values[inw].mean() if inw.sum() >= 6 else g.resid.values.mean())
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")]); pm["mon"] = pm.t
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
def solve(sub):
    a = sub.groupby(["row", "pair"]).agg(dvv=("anom", "mean"), nm=("anom", "size")).reset_index()
    if len(a) < 20: return None
    sig = np.array([max(1e-3, sig_pair.get(p, TVM)) for p in a.pair]); w = 1/sig; Gr = Gc[a.row.values]
    AtA = (Gr*w[:, None]).T@(Gr*w[:, None]); Minv = inv(AtA+LAM**2*LtL+1e-4*np.eye(Mf))
    m = Minv@((Gr.T*w**2)@a.dvv.values); ill = np.isin(np.arange(Mf), a.row.values); m[~ill] = np.nan; return m
months = pd.period_range("2010-06", "2026-06", freq="M")
frames = []
for m in months:
    win = pm[(pm.t <= m) & (pm.t >= m-2)]                       # trailing 3-month window
    mm = solve(win); frames.append((str(m), mm if mm is not None else np.full(Mf, np.nan)))
print(f"monthly frames: {len(frames)}; illuminated cells (median): {int(np.median([np.isfinite(fr[1][deep]).sum() for fr in frames]))}")

# --- animation ---
di = np.where(deep)[0]; plon, plat = clon[di], clat[di]
ulon = np.round(np.arange(plon.min(), plon.max()+STEP/2, STEP), 3); ulat = np.round(np.arange(plat.min(), plat.max()+STEP/2, STEP), 3)
lox = {v: i for i, v in enumerate(ulon)}; lax = {v: i for i, v in enumerate(ulat)}
edx = np.append(ulon-STEP/2, ulon[-1]+STEP/2); edy = np.append(ulat-STEP/2, ulat[-1]+STEP/2); LON, LAT = np.meshgrid(edx, edy)
def zof(mm):
    Z = np.full((len(ulat), len(ulon)), np.nan)
    for ci in di:
        v = mm[ci]
        if np.isfinite(v):
            iy = lax.get(round(clat[ci], 3)); ix = lox.get(round(clon[ci], 3))
            if iy is not None and ix is not None: Z[iy, ix] = v
    return Z
proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(6, 9)); ax = plt.axes(projection=proj)
ax.set_extent([plon.min()-0.6, plon.max()+0.6, plat.min()-0.6, plat.max()+0.6], crs=proj)
ax.add_feature(cfeature.LAND, facecolor="#f0f0f0"); ax.add_feature(cfeature.OCEAN, facecolor="#dbeafe"); ax.coastlines(resolution="50m", lw=.6)
pc = ax.pcolormesh(LON, LAT, np.ma.masked_invalid(zof(frames[0][1])), cmap="RdBu_r", vmin=-0.4, vmax=0.4, transform=proj)
plt.colorbar(pc, ax=ax, fraction=.04, label="δβ/β (%, model)"); ttl = ax.set_title(frames[0][0], fontsize=13)
def upd(i):
    pc.set_array(np.ma.masked_invalid(zof(frames[i][1])).ravel()); ttl.set_text(f"deep δβ/β  {frames[i][0]}  (3-mo trailing)"); return pc, ttl
anim = FuncAnimation(fig, upd, frames=len(frames), blit=False)
anim.save(f"{RD}/dbb_monthly_movie.mp4", writer="ffmpeg", fps=6, dpi=110)
print(f"wrote {RD}/dbb_monthly_movie.mp4 ({len(frames)} frames)")
