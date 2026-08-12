"""Illustrate the composite for ONE cell: is dv/v lower in tremor-active (ETS) months than tremor-quiet months?
Picks the best-DATA-COVERED deep cell (not the biggest effect), shows its monthly dv/v with ETS months marked
and the two means. Reuses the exact prep of ets_composite_confirm.py. Mirror-free."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g50")
d = np.load(f"{RD}/G.npz", allow_pickle=True)
clat = d["cell_lat"]; clon = d["cell_lon"]; depth = d["depth_km"]; cellids = d["cell"]; deep = depth > 30
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag+"|"+pm.cell
# era-split
et = pd.read_csv(f"{RD}/era_table.csv")
bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
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
# per-cell local ETS months
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"], errors="coerce"); tr = tr.dropna(subset=["t", "lat", "lon"])
lat0, lon0 = clat.mean(), clon.mean()
def xy(la, lo): return np.column_stack([(np.asarray(lo)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(la)-lat0)*111])
nn = cKDTree(xy(tr.lat.values, tr.lon.values)).query_ball_point(xy(clat, clon), r=100)
tr_ym = tr.t.dt.to_period("M"); tr_yr = tr.t.dt.year.values
cell_ets = {}
for k, cid in enumerate(cellids):
    ev = tr_ym.values[nn[k]]; yy = tr_yr[nn[k]]
    if len(ev) == 0: cell_ets[cid] = set(); continue
    mc = pd.DataFrame({"ym": ev, "yr": yy}).groupby(["yr", "ym"]).size().reset_index(name="n")
    z = mc.groupby("yr").n.transform(lambda s: (s-s.mean())/(s.std() if s.std() > 0 else 1))
    cell_ets[cid] = set(pd.PeriodIndex(mc.ym[z > 0.5]))
# cell-aggregate monthly dv/v (mean across the cell's pairs)
cm = pm.groupby(["cell", "t"]).anom.mean().reset_index()
# pick best-covered deep cell that has >=8 ETS and >=8 quiet months
best = None
for cid in cellids[deep]:
    g = cm[cm.cell == cid]; e = np.array([tt in cell_ets.get(cid, set()) for tt in g.t])
    if e.sum() >= 8 and (~e).sum() >= 8 and (best is None or len(g) > best[1]):
        best = (cid, len(g))
cid = best[0]; g = cm[cm.cell == cid].sort_values("t"); e = np.array([tt in cell_ets.get(cid, set()) for tt in g.t])
ci = list(cellids).index(cid); mE = g.anom[e].mean(); mQ = g.anom[~e].mean()
print(f"cell {cid} ({clat[ci]:.2f}N {clon[ci]:.2f}E, {depth[ci]:.0f} km) | {len(g)} months, {int(e.sum())} ETS / {int((~e).sum())} quiet")
print(f"mean dv/v ETS = {mE:+.4f}% | mean dv/v quiet = {mQ:+.4f}% | ETS - quiet = {mE-mQ:+.4f}%")

fig, a = plt.subplots(figsize=(12, 5))
tx = g.t.dt.to_timestamp()
a.axhline(0, color="k", lw=.5)
a.plot(tx, g.anom, "-", color="#999", lw=.8, zorder=1)
a.scatter(tx[~e], g.anom[~e], s=26, c="#2c7fb8", label=f"tremor-quiet months (n={int((~e).sum())})", zorder=3)
a.scatter(tx[e], g.anom[e], s=42, c="#d7301f", marker="D", label=f"tremor-ACTIVE (ETS) months (n={int(e.sum())})", zorder=4)
a.axhline(mQ, color="#2c7fb8", ls="--", lw=1.6, label=f"mean quiet = {mQ:+.3f}%")
a.axhline(mE, color="#d7301f", ls="--", lw=1.6, label=f"mean ETS = {mE:+.3f}%")
a.annotate("", xy=(tx.iloc[-1], mE), xytext=(tx.iloc[-1], mQ), arrowprops=dict(arrowstyle="<->", color="k"))
a.text(tx.iloc[-1], (mE+mQ)/2, f"  ETS−quiet\n  = {mE-mQ:+.3f}%", va="center", fontsize=10)
a.set_xlabel("date"); a.set_ylabel("dv/v anomaly (%, common-mode removed)")
a.set_title(f"One deep cell {clat[ci]:.2f}°N {clon[ci]:.2f}°E ({depth[ci]:.0f} km) — is dv/v lower in tremor-active months?\n"
            f"mean ETS {mE:+.3f}%  vs  mean quiet {mQ:+.3f}%   →  difference {mE-mQ:+.3f}% ({'lower in ETS' if mE<mQ else 'higher in ETS'})", fontsize=11)
a.legend(fontsize=8.5, loc="best", ncol=2); a.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{RD}/ets_one_cell.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/ets_one_cell.png")
