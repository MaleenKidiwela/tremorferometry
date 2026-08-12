#!/usr/bin/env python
"""Explicit COMMON-MODE subtraction + coarse-grid DEEP dv/v inversion (user's plan).
(1) DETERMINE the common mode = per-station median dv/v across that station's families each month (the shallow,
    near-receiver signal shared by all families at a station).
(2) SUBTRACT it -> family-specific residual = the deep, near-source signal.
(3) INVERT the residual on the COARSE grid (0.5 deg, ~50 km cells = one resolution element) for deep delta-beta/beta.
Corrections kept: conditional per-pair deseason; per-family baseline demean; instrument era-split (station x era
common mode); NO site terms (common mode already removed explicitly). Matched scrambled-year null for significance.
Outputs res_catalog_g50/coarse_cm.npz + prints; figures via plot_coarse_cm.py.  Mirror-FREE.
"""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g50")
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float)                                   # +K convention
f = d["captured"].astype(float); clat = d["cell_lat"]; clon = d["cell_lon"]; depth = d["depth_km"]
mixed = d["depth_mixed"] if "depth_mixed" in d.files else np.zeros(len(depth), bool)  # fuzzy near-30km cells
cxy = d["cxy"]; cellids = d["cell"]; ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; deep = depth > 30   # family-median depth
_, nbr = cKDTree(cxy).query(cxy, k=min(6, Mf)); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T@L
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag + "|" + pm.cell

# instrument era-split -> the common mode is per (station, era)
et = pd.read_csv(f"{RD}/era_table.csv") if os.path.exists(f"{RD}/era_table.csv") else pd.read_csv("fault_tomography/inversion/res_catalog/era_table.csv")
bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
ts = pm.t.dt.to_timestamp().values; eid = np.zeros(len(pm), int); strad = np.zeros(len(pm), bool)
for i, (tg, t) in enumerate(zip(pm.tag.values, ts)):
    bs = bnd_of.get(tg, [])
    if bs:
        eid[i] = int(sum(1 for b in bs if t >= np.datetime64(b))); strad[i] = any(abs((pd.Timestamp(t)-b).days) <= 35 for b in bs)
pm["stera"] = pm.tag + "__e" + eid.astype(str); pm = pm[~strad].copy()

# (1) conditional per-pair deseason (annual+semiannual, only if F-test significant)
from scipy.stats import f as fdist
def harm(t): yr = 12.; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def deseason(g):
    v = g.dvv_month.values.astype(float); t = (g.t-g.t.min()).apply(lambda x: x.n).values.astype(float); n = len(v)
    if n >= 24:
        X = np.column_stack([np.ones_like(t), harm(t)]); b, *_ = np.linalg.lstsq(X, v, rcond=None)
        RSSf = ((v-X@b)**2).sum(); RSSc = ((v-v.mean())**2).sum(); F = ((RSSc-RSSf)/4)/(RSSf/(n-5)) if RSSf > 0 else 0
        if 1-fdist.cdf(F, 4, n-5) < 0.05: v = v - harm(t)@b[1:]
    return v
pm = pm.sort_values(["pair", "t"]); pm["ds"] = np.concatenate([deseason(g) for _, g in pm.groupby("pair")])

# (2) COMMON MODE = per (station-era, month) median across families; SUBTRACT
pm["cmode"] = pm.groupby(["stera", "ym"]).ds.transform("mean")
pm["resid"] = pm.ds - pm.cmode
# per-family baseline demean (2019-2024 window) -> anomaly
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); base = g.resid.values[inw].mean() if inw.sum() >= 6 else g.resid.values.mean()
    return g.resid.values - base
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
rho = pm.groupby("pair").anom.apply(lambda a: float(np.clip(np.corrcoef(a.values[:-1], a.values[1:])[0, 1], 0, .95)) if len(a) > 8 and a.std() > 0 else 0.).to_dict()

# common-mode time series (network median of the per-station common modes, for display)
cm_net = pm.groupby("year").cmode.mean()

# (3) INVERT residual per year on coarse grid (NO site terms; common mode already removed)
def agg(sub):
    a = sub.groupby(["row", "pair"]).agg(dvv=("anom", "mean"), nm=("anom", "size")).reset_index()
    nf = lambda p: (1-rho.get(p, 0))/(1+rho.get(p, 0))
    a["sig"] = [max(1e-3, sig_pair.get(p, TVM))/np.sqrt(max(1., nm*nf(p))) for p, nm in zip(a.pair, a.nm)]
    return a
LAM = float(os.environ.get("LAM", "0.3"))
def solve(a, cov=False):
    w = 1/a.sig.values; Gr = Gc[a.row.values]; AtA = (Gr*w[:, None]).T@(Gr*w[:, None])
    reg = LAM**2*LtL + 1e-4*np.eye(Mf); Minv = inv(AtA+reg); m = Minv@((Gr.T*w**2)@a.dvv.values)
    if cov: return m, np.sqrt(np.clip(np.diag(Minv@AtA@Minv), 0, None))
    return m
wins = list(range(2010, 2026)); MODEL = np.full((len(wins), Mf), np.nan); SIGM = np.full((len(wins), Mf), np.nan); idx = np.full(len(wins), np.nan)
for wi, y in enumerate(wins):
    a = agg(pm[pm.year == y])
    if len(a) < 30 or a.row.nunique() < 10: continue
    m, sm = solve(a, cov=True); ill = np.isin(np.arange(Mf), a.row.values)
    m[~ill] = np.nan; sm[~ill] = np.nan; MODEL[wi] = m; SIGM[wi] = sm; idx[wi] = np.nanmean(m[deep])
# matched null
nidx = []
for s in range(50):
    rng = np.random.RandomState(s); pn = pm[pm.year < 2026].copy(); pn["year"] = pn.groupby("pair").year.transform(lambda x: rng.permutation(x.values))
    v = [np.nanmean(solve(agg(pn[pn.year == y]))[deep]) for y in wins[:-1] if len(agg(pn[pn.year == y])) >= 30]
    nidx.append(np.nanmax(np.abs(v)) if v else np.nan)
NULL = float(np.nanpercentile(nidx, 95)); realmax = float(np.nanmax(np.abs(idx))); pct = 100*float(np.mean(np.array(nidx) < realmax))

# ---- ETS-phase composite: deep dv/v in ETS (tremor-active) months minus inter-ETS months ----
# per-CELL ETS months = tremor within 100 km, monthly count z-scored WITHIN year, z>0.5 (async N/S ETS + rate growth removed)
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"], errors="coerce"); tr = tr.dropna(subset=["t", "lat", "lon"])
lat0, lon0 = clat.mean(), clon.mean()
def _xy(la, lo): return np.column_stack([(np.asarray(lo)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(la)-lat0)*111])
_tree = cKDTree(_xy(tr.lat.values, tr.lon.values)); tr_ym = tr.t.dt.to_period("M"); tr_yr = tr.t.dt.year.values
_nn = _tree.query_ball_point(_xy(clat, clon), r=100)
cell_ets = {}
for k, cid in enumerate(cellids):
    ev = tr_ym.values[_nn[k]]; yy = tr_yr[_nn[k]]
    if len(ev) == 0: cell_ets[cid] = set(); continue
    dfc = pd.DataFrame({"ym": ev, "yr": yy}); mc = dfc.groupby(["yr", "ym"]).size().reset_index(name="n")
    z = mc.groupby("yr").n.transform(lambda s: (s-s.mean())/(s.std() if s.std() > 0 else 1))
    cell_ets[cid] = set(mc.ym[z > 0.5].astype(str))
pm["is_ets"] = [ym in cell_ets.get(c, set()) for c, ym in zip(pm.cell, pm.ym)]
comp = []
for p, g in pm.groupby("pair"):
    e = g.anom[g.is_ets]; iq = g.anom[~g.is_ets]
    if len(e) >= 3 and len(iq) >= 3:
        comp.append(dict(row=int(g.row.iloc[0]), pair=p, dvv=float(e.mean()-iq.mean()), nm=int(min(len(e), len(iq)))))
comp = pd.DataFrame(comp)
ets_map = np.full(Mf, np.nan); ets_sig = np.full(Mf, np.nan); ets_deep_mean = np.nan
if len(comp) >= 20:
    nf = lambda p: (1-rho.get(p, 0))/(1+rho.get(p, 0))
    comp["sig"] = [max(1e-3, sig_pair.get(p, TVM))/np.sqrt(max(1., nm*nf(p))) for p, nm in zip(comp.pair, comp.nm)]
    m, sm = solve(comp, cov=True); illc = np.isin(np.arange(Mf), comp.row.values)
    ets_map = np.where(illc, m, np.nan); ets_sig = np.where(illc, sm, np.nan); ets_deep_mean = float(np.nanmean(ets_map[deep]))

np.savez_compressed(f"{RD}/coarse_cm.npz", wins=np.array(wins), MODEL=MODEL, SIGM=SIGM, idx=idx, clat=clat, clon=clon,
                    depth=depth, deep=deep, mixed=mixed, null=NULL, pct=pct, cm_net=cm_net.reindex(wins).values, lam=LAM, cellids=cellids,
                    ets_map=ets_map, ets_sig=ets_sig)
print(f"COARSE-GRID COMMON-MODE-SUBTRACTED DEEP INVERSION ({os.path.basename(RD)}, {int(deep.sum())} deep cells, mirror-FREE)")
print(f"common mode removed (network median dv/v by year, %): " + " ".join(f"{y}:{cm_net.get(y,np.nan):+.3f}" for y in wins if y in cm_net.index))
print(f"deep residual index by year (model %): " + " ".join(f"{y}:{idx[wi]:+.2f}" for wi, y in enumerate(wins) if np.isfinite(idx[wi])))
print(f"NULL 95pct {NULL:.3f}% | real max |idx| {realmax:.3f}% at {pct:.0f}th pctile -> {'not significant (<95th)' if pct < 95 else 'significant'}")
print(f"ETS composite (ETS - inter-ETS): {len(comp)} pairs | deep-cell mean {ets_deep_mean:+.4f}% (model)")
print(f"wrote {RD}/coarse_cm.npz")
