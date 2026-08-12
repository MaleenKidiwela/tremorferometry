"""Deep INDEX at MONTHLY cadence (single-panel time series) for the coarse common-mode grids
(family-centroid, mirror-free). Solves the deep inversion on a TRAILING 3-month window each month;
plots the deep network-mean delta-beta/beta anomaly per month + a scrambled-time null band.
Reuses monthly_movie.py's data prep exactly. Reads {RESDIR}/G.npz + pair_months.parquet + era_table.csv."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20")
LAM = float(os.environ.get("LAM", "0.5")); STEP = os.environ.get("STEP", "0.2"); WIN = int(os.environ.get("WIN", "3"))
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); depth = d["depth_km"]
cxy = d["cxy"]; ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; deep = depth > 30
_, nbr = cKDTree(cxy).query(cxy, k=min(6, Mf)); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
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
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
months = pd.period_range("2010-06", "2026-06", freq="M")
prow = {p: r for p, r in zip(pm.pair, pm.row)}                       # pair -> G row (1:1)
piv = pm.pivot_table(index="pair", columns="t", values="anom", aggfunc="mean").reindex(columns=months)
prows = np.array([prow[p] for p in piv.index]); psig = np.array([max(1e-3, sig_pair.get(p, TVM)) for p in piv.index])
regmat = LAM**2*LtL + 1e-4*np.eye(Mf); allrows = np.arange(Mf)
def deep_index(A):                                                  # A = pairs x months (trailing mean) -> monthly deep index
    out = np.full(A.shape[1], np.nan)
    for k in range(A.shape[1]):
        col = A[:, k]; ok = np.isfinite(col)
        if ok.sum() < 20: continue
        rows = prows[ok]; w = 1/psig[ok]; Gr = Gc[rows]
        AtA = (Gr*w[:, None]).T@(Gr*w[:, None]); m = inv(AtA+regmat)@((Gr.T*w**2)@col[ok])
        m[~np.isin(allrows, rows)] = np.nan; out[k] = np.nanmean(m[deep])
    return out
def trailing(P): return P.T.rolling(WIN, min_periods=1).mean().T.values   # trailing WIN-month mean per pair
idx = deep_index(trailing(piv))
# scrambled-time null: permute each pair's observed months, recompute the monthly index.
# The reported statistic is max|idx| over ALL months, so the null must be the distribution of the
# per-scramble max|idx| (NOT a pooled pointwise band) -- else the multiple-months search is uncounted.
base = piv.values; nvals = []
for s in range(40):
    A = base.copy()
    for r in range(A.shape[0]):
        fin = np.where(np.isfinite(A[r]))[0]
        if len(fin) > 1: A[r, fin] = A[r, fin][np.random.RandomState(1000*s+r).permutation(len(fin))]
    nvals.append(deep_index(pd.DataFrame(A, columns=months).T.rolling(WIN, min_periods=1).mean().T.values))
nvals = np.array(nvals)
nullmax = np.nanmax(np.abs(nvals), axis=1)                    # per-scramble max|idx| over months
NULLMAX95 = float(np.nanpercentile(nullmax, 95)); BAND = float(np.nanpercentile(np.abs(nvals), 95))
realmax = float(np.nanmax(np.abs(idx))); pctile = 100*float(np.mean(nullmax < realmax))
nfin = int(np.isfinite(idx).sum())
verdict = "significant (>95th of max-null)" if pctile >= 95 else "NOT significant (within max-null)"
print(f"monthly deep index: {nfin}/{len(months)} months | max|idx| {realmax:.3f}% at {pctile:.0f}th pctile of max-null (95th={NULLMAX95:.3f}%) -> {verdict}")

fig, a = plt.subplots(figsize=(11, 4.8))
tx = months.to_timestamp()
a.axhspan(-BAND, BAND, color="grey", alpha=.18, label=f"pointwise null (±{BAND:.2f}%)")
a.axhline(NULLMAX95, color="grey", ls="--", lw=1); a.axhline(-NULLMAX95, color="grey", ls="--", lw=1, label=f"max-over-months null 95th (±{NULLMAX95:.2f}%)")
a.plot(tx, idx, "-", color="#1a1a2e", lw=1.4, label=f"deep network-mean anomaly ({WIN}-mo trailing)")
a.axhline(0, color="k", lw=.6)
a.set_xlabel("date"); a.set_ylabel("deep δβ/β anomaly (%, model)")
a.set_title(f"Deep interface index — MONTHLY ({STEP}° grid, {int(deep.sum())} deep cells, common-mode, mirror-free)\n"
            f"max |idx| {realmax:.2f}% at {pctile:.0f}th pctile of max-over-months null → {verdict}", fontsize=10.5)
a.legend(fontsize=9, loc="upper right"); a.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{RD}/deep_index_monthly.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/deep_index_monthly.png")
