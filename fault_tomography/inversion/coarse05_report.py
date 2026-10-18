#!/usr/bin/env python
"""Post-processing of the COARSE (0.5-deg band x depth band) daily inversion in fault_tomography/inversion/coarse05:
 (1) per-cell transmission/resolution: inject +1 % into each coarse cell at 3 representative days (2013/2018/2023) and
     recover through the production day_solve (site terms included) -> diagonal transmission and crosstalk matrix;
 (2) posterior sigma per cell (median over frames) and the noise of 180-d means;
 (3) latitude-time figures of the DEEP and UP-DIP bands (90-d mean) and a per-cell table.
Outputs (coarse05/): coarse05_transmission.npz, coarse05_cells_summary.csv, coarse05_lat_time.png
"""
import os, sys, numpy as np, pandas as pd
os.environ["RESDIR"] = os.environ.get("RESDIR", "fault_tomography/inversion/coarse05")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "1")
sys.path.insert(0, "fault_tomography/inversion"); import invert_dvv_4d_daily as inv
RD = inv.OUT; z = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); days = pd.to_datetime(z["days"]); MODEL = np.asarray(z["MODEL"], float); SIGM = np.asarray(z["SIGM"], float)
cells = list(z["cell"]); M = len(cells); LAMF, LAMS = float(z["lam_f"]), float(z["lam_s"]); dep = np.asarray(z["depth"], float); clat = np.asarray(z["clat"], float)
rep = [int(days.get_indexer([pd.Timestamp(x)], method="nearest")[0]) for x in ("2013-06-15", "2018-06-15", "2023-06-15")]
Tq = np.zeros((len(rep), M, M))
for k, t in enumerate(rep):
    idx = inv.BY_DAY[t]
    for j in range(M):
        m_inj = np.zeros(M); m_inj[j] = 1.0; d = inv.Gc[inv.ROW[idx]] @ m_inj
        m, _, _, _ = inv.day_solve(t, np.zeros(inv.P), LAMF, LAMS, d_override=d); Tq[k, :, j] = m
diag = np.array([np.diag(Tq[k]) for k in range(len(rep))]); Tmean = diag.mean(0)
leak = np.array([np.abs(Tq[k] - np.diag(np.diag(Tq[k]))).sum(0) for k in range(len(rep))]).mean(0)      # sum of |off-diagonal| per injected cell
sig = np.nanmedian(SIGM, 0); n180 = 180 / 75.0
tab = pd.DataFrame(dict(cell=cells, lat=clat, depth_km=dep, T_diag=Tmean, T_2013=diag[0], T_2023=diag[2], crosstalk_sum=leak, sigma_frame=sig, sigma_180d=sig / np.sqrt(n180),
                        frames_valid=np.isfinite(MODEL).sum(0), model_sd=np.nanstd(MODEL, 0)))
tag = os.path.basename(RD); tab.to_csv(f"{RD}/{tag}_cells_summary.csv", index=False)
np.savez_compressed(f"{RD}/{tag}_transmission.npz", Tq=Tq.astype(np.float32), rep_days=days[rep].astype(str).values, cell=np.array(cells, object))
print(f"coarse cells {M} | lam_f {LAMF:.3f} | transmission diag: deep median {np.median(Tmean[dep>=32]):.2f} (2013 {np.median(diag[0][dep>=32]):.2f}, 2023 {np.median(diag[2][dep>=32]):.2f}), up-dip median {np.median(Tmean[dep<32]):.2f} | crosstalk sum median {np.median(leak):.2f}")
print(f"posterior sigma per frame: deep median {np.median(sig[dep>=32]):.3f}%, up-dip {np.median(sig[dep<32]):.3f}% | 180-d mean sigma: deep {np.median(sig[dep>=32])/np.sqrt(n180):.3f}%")
print(tab.round(3).sort_values(["depth_km", "lat"]).to_string(index=False))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
clon = np.asarray(z["clon"], float)
fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
for a, (nm, sel) in zip(ax, (("DEEP cells (>=32 km)", dep >= 32), ("UP-DIP cells (<32 km)", dep < 32))):
    idx = np.where(sel)[0]; order = idx[np.lexsort((clon[idx], clat[idx]))]; Z = pd.DataFrame(MODEL[:, order], index=days).rolling("90D", min_periods=20).mean().values
    im = a.pcolormesh(days, np.arange(len(order)), Z.T, cmap="RdBu_r", vmin=-0.6, vmax=0.6, shading="nearest"); a.set_yticks(np.arange(len(order))); a.set_yticklabels([cells[k] for k in order], fontsize=6); a.set_title(f"{tag} (corrected kernel): {nm}, 90-d mean δβ/β (%)", fontsize=10); plt.colorbar(im, ax=a, pad=.01)
ax[-1].set_xlim(pd.Timestamp("2010-01-01"), pd.Timestamp("2026-09-01")); fig.savefig(f"{RD}/{tag}_lat_time.png", dpi=110, bbox_inches="tight")
# annual-mean maps
yrs = list(range(2011, 2027)); fig, ax = plt.subplots(2, 8, figsize=(22, 9))
for k, y in enumerate(yrs):
    a = ax[k // 8, k % 8]; sel = (days.year == y); v = np.nanmean(MODEL[sel], 0); sc = a.scatter(clon, clat, c=v, cmap="RdBu_r", vmin=-0.6, vmax=0.6, s=180, marker="s", edgecolor="k", linewidth=.3); a.set_title(str(y), fontsize=9); a.set_aspect(1 / np.cos(np.radians(46))); a.set_xlim(-127, -121); a.set_ylim(39.5, 50.5); a.tick_params(labelsize=6)
fig.colorbar(sc, ax=ax, shrink=.5, label="annual mean δβ/β (%)"); fig.suptitle(f"{tag}: annual means per cell", fontsize=11); fig.savefig(f"{RD}/{tag}_annual_maps.png", dpi=100, bbox_inches="tight")
print(f"wrote {RD}/{tag}_lat_time.png, {tag}_annual_maps.png, {tag}_cells_summary.csv, {tag}_transmission.npz")
