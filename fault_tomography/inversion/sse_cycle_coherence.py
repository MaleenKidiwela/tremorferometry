#!/usr/bin/env python
"""Two continuous-record lenses on 'does coda dv/v track modeled slow slip' at the 125-km deep-patch scale (N+C+S):
 (1) SAWTOOTH: correlation of deseasoned patch-mean dv/v with the model's cumulative thrust slip (relative to the
     long-term trend = the inter-episode loading sawtooth) and with the trailing slip-in-window; null = rigid circular
     shift of the model series (2000 draws). Products: interface (production + gateD0.3), family TOTAL, family
     site-removed, SITES of patch stations.
 (2) COHERENCE at the recurrence band: magnitude-squared coherence between patch-mean dv/v and slip-in-window
     (Welch, 2-yr segments) integrated over 0.6-1.2 cycles/yr, vs the same shift null.
Outputs: res_catalog_sse/cycle_coherence.csv, figures/sse_cycle_coherence.png
"""
import numpy as np, pandas as pd
from scipy.signal import coherence
SD = "fault_tomography/invers ion/res_catalog_sse".replace(" ", ""); RD = "fault_tomography/inversion/res_catalog"
rng = np.random.RandomState(5); ND = 2000
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells)
clat = zc["cell_lat"]; dep = zc["iface_depth_km"]; pdist = zc["patch_dist_km"]; gi = pd.Series(np.arange(G), index=days)
cum = zc["slip_thrust_cum"].astype(float); inc = zc["slip_thrust_inc30"].astype(float)
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True); PAT = np.asarray(zt["patterns"], float); band = zt["patch_band"]; plat = zt["patch_lat"]
deepq = [q for q in np.argsort(plat) if band[q] == "deep" and 40.5 < plat[q] < 49.5]
def load_prod(path):
    z = np.load(path, allow_pickle=True); return pd.DataFrame(np.asarray(z["MODEL"], float), index=pd.to_datetime(z["days"])).reindex(days).values, z
Y, zp = load_prod(f"{RD}/inversion_4d_daily.npz"); Yg, _ = load_prod("fault_tomography/inversion/splits/gateD0.3/inversion_4d_daily.npz")
SITES = pd.DataFrame(np.asarray(zp["SITES"], float), index=pd.to_datetime(zp["days"]), columns=list(zp["sta_tags"])).reindex(days)
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "cell", "date", "dvv"]); cpos = {c: i for i, c in enumerate(cells)}
fd = fd[fd.cell.isin(cpos)]; fd["ci"] = fd.cell.map(cpos); fd["v"] = -fd.dvv; fd["dw"] = fd.v - fd.groupby(["tag", "date"]).v.transform("mean")
cd = fd.groupby(["ci", "date"]).agg(tot=("v", "mean"), dw=("dw", "mean")).reset_index(); okd = cd.date.isin(days)
TOT = np.full((G, M), np.nan); DW = np.full((G, M), np.nan)
TOT[gi[cd.date[okd]].values, cd.ci[okd].values] = cd.tot[okd].values; DW[gi[cd.date[okd]].values, cd.ci[okd].values] = cd.dw[okd].values
sta_rows = fd.groupby(["ci", "tag"]).size().reset_index(name="n"); del fd, cd
t0, t1 = pd.Timestamp("2011-01-01"), pd.Timestamp("2026-07-20"); win = (days >= t0) & (days <= t1)
def patch_mean(X, q):
    w = PAT[q]; fin = np.isfinite(X); num = np.nansum(np.where(fin, X, 0) * w[None, :], 1); den = (fin * w[None, :]).sum(1)
    out = num / np.where(den > 0, den, np.nan); out[(fin * w[None, :]).sum(1) / w.sum() < 0.5] = np.nan; return out
def deseason(y):
    s = pd.Series(y, index=days); return (s - s.groupby(s.index.month).transform("mean")).values
def corr(a, b):
    m = np.isfinite(a) & np.isfinite(b) & win; return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 500 else np.nan
def coh_band(a, b):
    m = np.isfinite(a) & np.isfinite(b) & win
    if m.sum() < 1500: return np.nan
    aa = pd.Series(a).interpolate(limit_direction="both").values[win]; bb = pd.Series(b).interpolate(limit_direction="both").values[win]
    f, C = coherence(aa, bb, fs=365.25, nperseg=730); sel = (f >= 0.6) & (f <= 1.2); return float(C[sel].mean())
rows = []
for q in deepq:
    w = PAT[q]; cs = np.where(w > 0.3)[0]
    sr = sta_rows[sta_rows.ci.isin(cs)].groupby("tag").n.sum(); sr = sr / sr.sum(); stas = [t for t in sr[sr >= 0.05].index if t in SITES.columns]
    prods = {"interface": deseason(patch_mean(Y, q)), "interface gateD0.3": deseason(patch_mean(Yg, q)), "family TOTAL": deseason(patch_mean(TOT, q)),
             "family site-removed": deseason(patch_mean(DW, q)), "SITES": deseason(SITES[stas].mean(axis=1).values) if stas else np.full(G, np.nan)}
    cum_q = patch_mean(cum, q); inc_q = patch_mean(inc, q)
    for pname, y in prods.items():
        r_cum, r_inc, ch = corr(y, cum_q), corr(y, inc_q), coh_band(y, inc_q)
        nc, ni, nch = [], [], []
        for _ in range(ND):
            s = rng.randint(365, G - 365); cs_ = np.roll(cum_q, s); is_ = np.roll(inc_q, s); nc.append(corr(y, cs_)); ni.append(corr(y, is_))
            if _ < 300: nch.append(coh_band(y, is_))
        nc, ni, nch = np.array(nc), np.array(ni), np.array(nch)
        rows.append(dict(patch_lat=round(plat[q], 1), product=pname, r_cumslip=r_cum, p_cum=float(np.nanmean(np.abs(nc) >= abs(r_cum))), r_slipwin=r_inc,
                         p_inc=float(np.nanmean(np.abs(ni) >= abs(r_inc))), coh_recur=ch, p_coh=float(np.nanmean(nch >= ch)) if np.isfinite(ch) else np.nan, null_sd_r=float(np.nanstd(nc))))
    print(f"patch {plat[q]:.1f}N done", flush=True)
tab = pd.DataFrame(rows); tab.to_csv(f"{SD}/cycle_coherence.csv", index=False)
print(tab.round(3).to_string(index=False))
for pname in tab["product"].unique():
    t = tab[tab["product"] == pname]
    print(f"{pname:20s}: median r(cum slip) {t.r_cumslip.median():+.3f} (patches with p<0.05: {(t.p_cum < 0.05).sum()}/{len(t)}) | median r(slip-in-window) {t.r_slipwin.median():+.3f} (p<0.05: {(t.p_inc < 0.05).sum()}) | coherence 0.6-1.2 cpy median {t.coh_recur.median():.3f} (p<0.05: {(t.p_coh < 0.05).sum()})")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 3, figsize=(16, 4))
for k, (col, ttl) in enumerate((("r_cumslip", "corr with cumulative slip (loading sawtooth)"), ("r_slipwin", "corr with slip in trailing window"), ("coh_recur", "coherence 0.6-1.2 cycles/yr"))):
    for pname, mk in zip(tab["product"].unique(), "osd^v"):
        t = tab[tab["product"] == pname]; ax[k].plot(t.patch_lat, t[col], marker=mk, lw=.8, label=pname)
    ax[k].axhline(0, color="k", lw=.5); ax[k].set_title(ttl, fontsize=9); ax[k].set_xlabel("deep patch latitude")
ax[0].legend(fontsize=7); fig.suptitle("Continuous-record lenses: patch-mean coda dv/v vs modeled slow slip (deseasoned; nulls = circular shift)", fontsize=10)
fig.savefig("figures/sse_cycle_coherence.png", dpi=110, bbox_inches="tight"); print("wrote figures/sse_cycle_coherence.png, cycle_coherence.csv")
