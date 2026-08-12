"""CURIOSITY run (does not touch production): collapse the whole interface to ONE shared delta-beta/beta value
per month, solved jointly with per-station site terms. Design per month: [interface col = each pair's TOTAL
interface coupling (Gc row-sum); site cols = station indicators]. Reads {RESDIR}/{G.npz, pair_months.parquet}.
Same conditional-deseason + 2019-24 demean as production. Writes interface_single_ts.png."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import lstsq
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20")
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); Mf = G.shape[1]
ptag = d["pair_tag"]; pcell = d["pair_cell"]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
iface_coup = Gc.sum(1)                                          # per-row total interface coupling (uniform-interface column)
FMED = float(np.median(f[keep]))
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag+"|"+pm.cell
# conditional per-pair deseason + 2019-24 demean (production-style) -> anom
def harm(t): return np.column_stack([np.sin(2*np.pi*t/12), np.cos(2*np.pi*t/12), np.sin(4*np.pi*t/12), np.cos(4*np.pi*t/12)])
def deseason(g):
    v = g.dvv_month.values.astype(float); t = (g.t-g.t.min()).apply(lambda x: x.n).values.astype(float); n = len(v)
    if n >= 24:
        X = np.column_stack([np.ones_like(t), harm(t)]); b, *_ = lstsq(X, v, rcond=None)
        RSSf = ((v-X@b)**2).sum(); RSSc = ((v-v.mean())**2).sum(); F = ((RSSc-RSSf)/4)/(RSSf/(n-5)) if RSSf > 0 else 0
        if 1-fdist.cdf(F, 4, n-5) < 0.05: v = v-harm(t)@b[1:]
    return v
pm = pm.sort_values(["pair", "t"]); pm["ds"] = np.concatenate([deseason(g) for _, g in pm.groupby("pair")])
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); return g.ds.values - (g.ds.values[inw].mean() if inw.sum() >= 6 else g.ds.values.mean())
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
RIDGE_S = 1e-3
def solve1(a):
    ut = pd.unique(a.tag.values); tc = {t: k for k, t in enumerate(ut)}; ns = len(ut)
    S = np.zeros((len(a), ns)); S[np.arange(len(a)), [tc[t] for t in a.tag]] = 1.0
    A = np.hstack([iface_coup[a.row.values][:, None], S])       # col 0 = interface, cols 1: = site terms
    w = 1/np.array([max(1e-3, sig_pair.get(p, TVM)) for p in a.pair]); Aw = A*w[:, None]
    reg = np.eye(1+ns)*RIDGE_S; reg[0, 0] = 1e-6                # tiny ridge; interface barely penalized
    m = np.linalg.inv(Aw.T@Aw + reg) @ ((A.T*w**2) @ a.dvv.values)
    return float(m[0])
def agg(sub):
    return sub.groupby(["row", "pair"]).agg(dvv=("anom", "mean"), tag=("tag", "first")).reset_index()
months = pd.period_range("2010-06", "2026-06", freq="M"); ts = np.full(len(months), np.nan)
for k, mth in enumerate(months):
    a = agg(pm[pm.t == mth])
    if len(a) >= 40 and a.tag.nunique() >= 8: ts[k] = solve1(a)
yrs = list(range(2010, 2027)); yts = np.full(len(yrs), np.nan)   # incl. 2026 (PARTIAL year)
for i, y in enumerate(yrs):
    a = agg(pm[pm.year == y])
    if len(a) >= 40 and a.tag.nunique() >= 8: yts[i] = solve1(a)
print(f"single-interface delta-beta/beta (model %): monthly std {np.nanstd(ts):.3f} | yearly range [{np.nanmin(yts):+.3f},{np.nanmax(yts):+.3f}] | data-space std {np.nanstd(ts)*FMED:.4f}%")

fig, ax = plt.subplots(figsize=(12, 5))
tx = months.to_timestamp()
ax.plot(tx, ts, "-", color="#888", lw=1, label="monthly")
ax.plot([pd.Timestamp(f"{y}-07-01") for y in yrs], yts, "o-", color="#c0392b", lw=2, ms=6, label="yearly")
ax.axhline(0, color="k", lw=.5)
ax.set_xlabel("date"); ax.set_ylabel("single-interface δβ/β (%, model)")
ax.set_title(f"Whole interface collapsed to ONE δβ/β value over time (+ per-station site terms, deseasoned)\n"
             f"[{os.path.basename(RD)}] monthly std {np.nanstd(ts):.2f}% model ≈ {np.nanstd(ts)*FMED:.3f}% data-space", fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{RD}/interface_single_ts.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/interface_single_ts.png")

# --- yearly-only clean figure ---
fig2, a2 = plt.subplots(figsize=(10, 5))
a2.plot(yrs[:-1], yts[:-1], "o-", color="#c0392b", lw=2.2, ms=8)
a2.plot([yrs[-1]], [yts[-1]], "D", color="#c0392b", mfc="white", ms=9, label="2026 (partial)"); a2.legend(fontsize=9)
for y, v in zip(yrs, yts):
    if np.isfinite(v): a2.annotate(f"{v:+.2f}", (y, v), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
a2.axhline(0, color="k", lw=.6)
a2.set_xlabel("year"); a2.set_ylabel("single-interface δβ/β (%, model)"); a2.set_xticks(yrs); a2.tick_params(axis="x", rotation=45)
a2.set_title(f"Whole interface = ONE δβ/β value per YEAR (+ site terms, deseasoned) [{os.path.basename(RD)}]\n"
             f"range [{np.nanmin(yts):+.2f},{np.nanmax(yts):+.2f}]% model ≈ [{np.nanmin(yts)*FMED:+.3f},{np.nanmax(yts)*FMED:+.3f}]% data-space", fontsize=11)
a2.grid(alpha=.3); fig2.tight_layout(); fig2.savefig(f"{RD}/interface_single_ts_yearly.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/interface_single_ts_yearly.png")
