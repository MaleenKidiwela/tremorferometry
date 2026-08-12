#!/usr/bin/env python
"""3-D VOLUME inversion of REAL dv/v (same grid as volume3d_checkerboard_new.png: 0.2 deg x 6 depth layers).
Common-mode removed per station (mean across its families); residual inverted for delta-beta/beta in the crustal
volume with a 3-D graph-Laplacian. Mirror-FREE. Per year. Output res_catalog/volume_cm.npz + per-layer maps."""
import os, sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
sys.path.insert(0, "fault_tomography/kernels")
from kernel import kernel_singlescatter
from scipy.stats import f as fdist
RC = "fault_tomography/inversion/res_catalog"
BETA, ELL, W1, W2 = 3.5, 40.0, 2.0, 4.0
HGRID, DEPTHS = 0.2, np.array([3., 12., 22., 32., 40., 48.])   # SAME grid as the checkerboard
pairs = pd.read_csv(f"{RC}/pairs.csv"); cells = pd.read_csv(f"{RC}/cells.csv").sort_values("cell").reset_index(drop=True)
cidx = {c: i for i, c in enumerate(cells.cell)}
lat0, lon0 = cells.cell_lat.mean(), cells.cell_lon.mean()
def proj(lat, lon): return np.column_stack([(np.asarray(lon)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(lat)-lat0)*111])
CXY = proj(cells.cell_lat.values, cells.cell_lon.values); CZ = cells.depth_km.values
# volume model grid
mlat = np.arange(cells.cell_lat.min()-0.3, cells.cell_lat.max()+0.3+1e-9, HGRID)
mlon = np.arange(cells.cell_lon.min()-0.3, cells.cell_lon.max()+0.3+1e-9, HGRID)
GLA, GLO = np.meshgrid(mlat, mlon, indexing="ij")
hlat = np.repeat(GLA.ravel(), len(DEPTHS)); hlon = np.repeat(GLO.ravel(), len(DEPTHS))
mz = np.tile(DEPTHS, GLA.size); mxy = proj(hlat, hlon); Mv = len(mz)
print(f"volume grid: {Mv} cells ({GLA.size} cols x {len(DEPTHS)} layers), {HGRID} deg horizontal")
# forward operator (source = interface cell, model = volume)
sxy = proj(pairs.sta_lat.values, pairs.sta_lon.values); G = np.zeros((len(pairs), Mv), np.float32)
for i, r in pairs.iterrows():
    ci = cidx[r.cell]; G[i] = -kernel_singlescatter(mxy, mz, CXY[ci], float(CZ[ci]), sxy[i], BETA, W1, W2, ELL)
hit = (np.abs(G) > 0).sum(0); ill = hit >= 3; G = G[:, ill]; mxy, mz, hlat, hlon = mxy[ill], mz[ill], hlat[ill], hlon[ill]
M = G.shape[1]; print(f"illuminated volume cells (>=3 rays): {M}")
prow = {(t, c): i for i, (t, c) in enumerate(zip(pairs.tag, pairs.cell))}
# 3-D graph Laplacian (depth stretched x3)
P3 = np.column_stack([mxy, mz*3.0]); _, nb = cKDTree(P3).query(P3, k=8); Lap = np.zeros((M, M))
for i in range(M):
    for j in nb[i, 1:]: Lap[i, i] += 1; Lap[i, j] -= 1
LtL = Lap.T@Lap

# --- real data: deseason + per-station common-mode removal + era-split ---
pm = pd.read_parquet(f"{RC}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["prow"] = [prow.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.prow >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag+"|"+pm.cell
et = pd.read_csv(f"{RC}/era_table.csv"); bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
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
pm["cmode"] = pm.groupby(["stera", "ym"]).ds.transform("mean")   # per-station COMMON MODE (mean across families)
pm["resid"] = pm.ds - pm.cmode
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); return g.resid.values - (g.resid.values[inw].mean() if inw.sum() >= 6 else g.resid.values.mean())
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))

LAM = float(os.environ.get("LAM", "0.5"))
def solve(sub):
    a = sub.groupby(["prow", "pair"]).agg(dvv=("anom", "mean"), nm=("anom", "size")).reset_index()
    if len(a) < 40: return None, None
    sig = np.array([max(1e-3, sig_pair.get(p, TVM))/np.sqrt(max(1., nm)) for p, nm in zip(a.pair, a.nm)])
    w = 1/sig; Gr = G[a.prow.values]; AtA = (Gr*w[:, None]).T@(Gr*w[:, None]); reg = LAM**2*LtL+1e-4*np.eye(M)
    Minv = inv(AtA+reg); m = Minv@((Gr.T*w**2)@a.dvv.values); return m, np.array([cid in set(a.prow.values) for cid in range(M)]) if False else m*0+1
wins = list(range(2010, 2026)); VOL = np.full((len(wins), M), np.nan)
for wi, y in enumerate(wins):
    m, _ = solve(pm[pm.year == y])
    if m is not None: VOL[wi] = m
np.savez_compressed(f"{RC}/volume_cm.npz", wins=np.array(wins), VOL=VOL, hlat=hlat, hlon=hlon, mz=mz, depths=DEPTHS, lam=LAM, hgrid=HGRID)
print(f"deep(>30km) layer volume-mean delta-beta/beta by year (%): " + " ".join(f"{y}:{np.nanmean(VOL[wi][mz>30]):+.2f}" for wi, y in enumerate(wins) if np.isfinite(VOL[wi]).any()))
print(f"wrote {RC}/volume_cm.npz ({M} illuminated cells x {len(wins)} years)")
