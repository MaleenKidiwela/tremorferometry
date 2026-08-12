"""Do the SITE TERMS (near-station) carry the seasonal cycle? Solve the joint inversion at MONTHLY resolution
WITHOUT deseasoning, so the seasonal signal is present; then check whether it lands in the per-station site
terms (near-receiver, expected seasonal) vs the interface cells (near-source, expected NOT seasonal).
Reads {RESDIR}/{G.npz, pair_months.parquet}. Writes site_seasonality.png + prints amplitudes."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
from scipy.stats import f as fdist
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20")
LAM_F = float(os.environ.get("LAM", "0.231")); LAM_RATIO = 0.05/0.40
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); cxy = d["cxy"]; Mf = G.shape[1]
ptag = d["pair_tag"]; pcell = d["pair_cell"]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
FMED = float(np.median(f[keep]))     # model->data conversion (cells are model-space, sites data-space)
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
_, nbr = cKDTree(cxy).query(cxy, k=min(6, Mf)); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T@L
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag+"|"+pm.cell; pm["mon"] = pm.t.dt.month
# NO deseason -> keep the seasonal signal. per-pair demean (2019-24) only, so pairs are comparable anomalies.
pm["dv"] = pm.dvv_month.astype(float)
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); return g.dv.values - (g.dv.values[inw].mean() if inw.sum() >= 6 else g.dv.values.mean())
pm = pm.sort_values(["pair", "t"]); pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
reg0 = np.zeros((Mf, Mf)); reg0[:Mf, :Mf] = LAM_F**2*LtL; reg0[np.arange(Mf), np.arange(Mf)] += 1e-4
def solve_month(a):
    ut = pd.unique(a.tag.values); tc = {t: k for k, t in enumerate(ut)}; ns = len(ut)
    S = np.zeros((len(a), ns)); S[np.arange(len(a)), [tc[t] for t in a.tag]] = 1.0
    A = np.hstack([Gc[a.row.values], S]); w = 1/np.array([max(1e-3, sig_pair.get(p, TVM)) for p in a.pair])
    Aw = A*w[:, None]; reg = np.zeros((Mf+ns, Mf+ns)); reg[:Mf, :Mf] = reg0
    si = np.arange(Mf, Mf+ns); reg[si, si] += (LAM_F*LAM_RATIO)**2
    m = inv(Aw.T@Aw + reg) @ ((A.T*w**2) @ a.dvv.values)
    return m[:Mf], {t: m[Mf+tc[t]] for t in ut}                    # cells, {station: site term}
months = pd.period_range("2010-06", "2026-06", freq="M")
site_rows = []; cell_series = np.full((len(months), Mf), np.nan)
for k, mth in enumerate(months):
    sub = pm[pm.t == mth]
    a = sub.groupby(["row", "pair"]).agg(dvv=("anom", "mean"), tag=("tag", "first")).reset_index()
    if len(a) < 40 or a.tag.nunique() < 8: continue
    cells, sites = solve_month(a)
    cell_series[k] = cells
    for st, v in sites.items(): site_rows.append((mth.month, st, float(v)))
sdf = pd.DataFrame(site_rows, columns=["mon", "sta", "val"])

# --- seasonality: fit annual+semiannual per station site-term series; and per interface cell ---
def seasonal_amp_frac(series_by_key):
    amps = []; nseas = 0; ntot = 0
    for key, g in series_by_key:
        v = g["val"].values if hasattr(g, "columns") else g
        mo = g["mon"].values if hasattr(g, "columns") else None
        if mo is None or len(v) < 24: continue
        t = np.arange(len(v))
        X = np.column_stack([np.ones_like(t, float), np.sin(2*np.pi*mo/12), np.cos(2*np.pi*mo/12), np.sin(4*np.pi*mo/12), np.cos(4*np.pi*mo/12)])
        b, *_ = np.linalg.lstsq(X, v, rcond=None); rssf = ((v-X@b)**2).sum(); rssc = ((v-v.mean())**2).sum()
        F = ((rssc-rssf)/4)/(rssf/(len(v)-5)) if rssf > 0 else 0; p = 1-fdist.cdf(F, 4, len(v)-5)
        amp = np.hypot(b[1], b[2]); amps.append(amp); ntot += 1; nseas += (p < 0.05)
    return np.array(amps), nseas, ntot
site_amps, site_ns, site_nt = seasonal_amp_frac(sdf.groupby("sta"))
# interface cells: build (mon, cell) frame
cm = []
for k, mth in enumerate(months):
    for ci in range(Mf):
        if np.isfinite(cell_series[k, ci]): cm.append((mth.month, ci, cell_series[k, ci]))
cdf = pd.DataFrame(cm, columns=["mon", "cell", "val"])
cell_amps, cell_ns, cell_nt = seasonal_amp_frac(cdf.groupby("cell"))
cell_amps = cell_amps * FMED                                     # model-space -> data-space (fair vs site terms)
# calendar-month composite (network mean); cells converted to data-space
site_cyc = sdf.groupby("mon").val.mean(); cell_cyc = cdf.groupby("mon").val.mean() * FMED
print(f"[data-space, FMED={FMED:.3f}]")
print(f"SITE terms (near-station): {site_ns}/{site_nt} stations seasonal (p<0.05) | median seasonal amp {np.nanmedian(site_amps):.3f}%")
print(f"INTERFACE cells (near-source): {cell_ns}/{cell_nt} cells seasonal (p<0.05) | median seasonal amp {np.nanmedian(cell_amps):.4f}% (model {np.nanmedian(cell_amps)/FMED:.3f}%)")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
mo = np.arange(1, 13)
ax[0].plot(site_cyc.index, site_cyc.values, "o-", color="#b8860b", lw=2, label=f"SITE terms (near-station), amp {np.nanmedian(site_amps):.2f}%")
ax[0].plot(cell_cyc.index, cell_cyc.values, "s-", color="#1a4d8f", lw=2, label=f"interface cells (near-source), amp {np.nanmedian(cell_amps):.2f}%")
ax[0].axhline(0, color="k", lw=.5); ax[0].set_xlabel("calendar month"); ax[0].set_ylabel("dv/v anomaly (%)")
ax[0].set_title("Seasonal cycle (network mean by calendar month)"); ax[0].legend(fontsize=9); ax[0].grid(alpha=.3); ax[0].set_xticks(mo)
ax[1].hist(site_amps, bins=25, alpha=.6, color="#b8860b", label=f"site terms ({100*site_ns/max(1,site_nt):.0f}% seasonal)")
ax[1].hist(cell_amps, bins=25, alpha=.6, color="#1a4d8f", label=f"interface cells ({100*cell_ns/max(1,cell_nt):.0f}% seasonal)")
ax[1].set_xlabel("per-series seasonal amplitude (%)"); ax[1].set_ylabel("count"); ax[1].set_title("Seasonal amplitude distribution"); ax[1].legend(fontsize=9); ax[1].grid(alpha=.3)
fig.suptitle("Do SITE TERMS carry the seasonal cycle? (monthly joint solve, NOT deseasoned)", fontsize=13)
fig.tight_layout(); fig.savefig(f"{RD}/site_seasonality.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/site_seasonality.png")
