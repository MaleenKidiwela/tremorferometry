"""Family 10 at PGC: dv/v by MWCS (moving-window cross-spectral) on SVD-Wiener denoised 3-component day stacks,
coda tS+2..tS+5 s (user 2026-10-12). Day stacks built exactly as s10_fam_3c_26.py (P>=0.5, per-event 3-C realign,
P-weighted day stacks, SVD-Wiener per component per epoch; reference = realigned epoch grand stack).
MWCS: 1.0-s sub-windows, 0.25-s step, Hann taper, zero-padded FFT, cross-spectrum smoothed over 5 bins,
coherence-weighted phase slope over 2-8 Hz -> delay dt per sub-window; dv/v = -(slope of dt vs lapse time after
picked S), one fit per day over all 3 components (weights = mean coherence / delay misfit). Quality = mean
coherence. Plot: no error bars, colour = coherence."""
import sys, numpy as np, pandas as pd, matplotlib
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import s10_measure as M, s10_svd_measure as V
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = 10; NS = "CN.PGC"; W0, W1 = 2.0, 5.0
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet"); rs = rs[rs.netsta == NS][["row_id", "epoch"]]
ser = pd.read_csv(f"{PIL}/tables/series.csv"); ser = ser[(ser.fam == FAM) & (ser.netsta == NS)].set_index("epoch")
z = np.load(f"{PIL}/win/{NS}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]; fs = 1 / (tau[1] - tau[0])
WL, STEP, NFFT, FMIN, FMAX, SM = 1.0, 0.25, 256, 2.0, 8.0, 5
f = np.fft.rfftfreq(NFFT, 1 / fs); band = (f >= FMIN) & (f <= FMAX)
def smooth(x): return np.convolve(x, np.ones(SM) / SM, mode="same")
def mwcs(cur, ref, tS):
    T, D, Wt, C = [], [], [], []
    nw = int(round(WL * fs)); hann = np.hanning(nw)
    for c in range(3):
        t0 = tS + W0
        while t0 + WL <= tS + W1 + 1e-9:
            i0 = int(np.searchsorted(tau, t0)); a = cur[c, i0:i0 + nw]; b = ref[c, i0:i0 + nw]
            if len(a) < nw: break
            A = np.fft.rfft((a - a.mean()) * hann, NFFT); B = np.fft.rfft((b - b.mean()) * hann, NFFT)
            X = A * np.conj(B); Xs = smooth(X.real) + 1j * smooth(X.imag)
            coh = np.abs(Xs) / np.sqrt(smooth(np.abs(A) ** 2) * smooth(np.abs(B) ** 2) + 1e-30)
            ph = np.unwrap(np.angle(Xs[band])); w = np.clip(coh[band], 1e-3, 0.999); w = w / (1 - w ** 2)  # coherence weights
            ff = f[band]; slope = np.sum(w * ff * ph) / np.sum(w * ff * ff)          # phase = -2 pi f dt (fit through origin)
            dt = -slope / (2 * np.pi)
            T.append(t0 + WL / 2 - tS); D.append(dt); Wt.append(coh[band].mean()); C.append(coh[band].mean())
            t0 += STEP
    T, D, Wt = map(np.array, (T, D, Wt))
    A = np.c_[np.ones_like(T), T]; Wm = np.diag(Wt ** 2)
    p = np.linalg.lstsq(A.T @ Wm @ A, A.T @ Wm @ D, rcond=None)[0]                # dt = a + b * lapse
    return -p[1] * 100, float(np.mean(C))
out = []
for ep_id, rr in rs.groupby("epoch"):
    g = rows[rows.row_id.isin(set(rr.row_id)) & rows.row_id.isin(pos)].reset_index(drop=True)
    if len(g) < 30: continue
    g["date"] = pd.to_datetime(g.date); W = z["W"][[pos[r] for r in g.row_id]].astype(np.float64)
    tS = float(ser.loc[ep_id, "tS_pick"]); w = g.P.to_numpy()
    G = np.array([np.nansum(np.nan_to_num(W[:, c]) * w[:, None], 0) / w.sum() for c in range(3)])
    Wa = M.realign(W, G, tS)[0]
    Ga = np.array([np.nansum(np.nan_to_num(Wa[:, c]) * w[:, None], 0) / w.sum() for c in range(3)])
    days = sorted(g.date.unique()); D = np.zeros((len(days), 3, len(tau))); keep = []
    for k, d in enumerate(days):
        gg = g.index[g.date == d]; ww = w[gg]
        if ww.sum() < 3: continue
        D[k] = [np.nansum(np.nan_to_num(Wa[gg, c]) * ww[:, None], 0) / ww.sum() for c in range(3)]; keep.append(k)
    D = D[keep]; days = [days[k] for k in keep]
    Df = np.stack([V.svdw(D[:, c], np.flatnonzero(tau < tS - 8)) for c in range(3)], 1)
    np.savez_compressed(f"{PIL}/tables/fam10_PGC_svd3c_{ep_id}.npz", Df=Df.astype(np.float32), D=D.astype(np.float32), days=np.array(days, dtype="datetime64[D]"), tau=tau, tS=tS, Ga=Ga)
    for k, d in enumerate(days):
        dv, coh = mwcs(Df[k], Ga, tS); out.append(dict(epoch=ep_id, date=d, dvv_mwcs=dv, coh=coh))
T = pd.DataFrame(out).sort_values("date"); T.to_csv(f"{PIL}/tables/fam10_PGC_mwcs_svd3c25_daily.csv", index=False)
K = T[(T.coh > 0.8) & (T.dvv_mwcs.abs() < 1)].copy()
dd = np.sort(K.date.unique()); epn = np.r_[0, np.cumsum(np.diff(dd).astype("timedelta64[D]").astype(int) > 20)]; K["ep"] = K.date.map(dict(zip(dd, epn)))
E = K.groupby("ep").apply(lambda q: pd.Series(dict(date=q.date.min(), dvv=np.average(q.dvv_mwcs, weights=q.coh), coh=q.coh.mean(), nd=len(q))), include_groups=False)
fig, ax = plt.subplots(figsize=(13, 4.5))
ax.scatter(K.date, K.dvv_mwcs, c=K.coh, cmap="viridis", vmin=0.8, vmax=1, s=12, alpha=0.6, lw=0)
ax.plot(E.date, E.dvv, "-", color="k", lw=1.2, zorder=3)
sc = ax.scatter(E.date, E.dvv, c=E.coh, cmap="viridis", vmin=0.8, vmax=1, s=25 + 10 * E.nd, edgecolor="k", lw=0.6, zorder=4)
ax.axhline(0, color="grey", lw=0.6, ls=":"); ax.set_ylabel("dv/v (%)"); fig.colorbar(sc, ax=ax, label="mean coherence (2–8 Hz)")
ax.set_title(f"Family 10 at PGC: dv/v by MWCS, SVD-Wiener 3-component day stacks, 2–5 s after picked S; small = days (coh > 0.8, {len(K)}), big + line = per-episode coherence-weighted mean\n"
             f"day sd {K.dvv_mwcs.std():.3f} %, episode sd {E.dvv.std():.3f} %", fontsize=9)
fig.tight_layout(); outp = f"{PIL}/figures/s10_fam10_PGC_mwcs_svd3c25.png"; fig.savefig(outp, dpi=120); print(outp)
print(f"days {len(T)}, coh>0.8 {len(K)}, day sd {K.dvv_mwcs.std():.3f} %, episode sd {E.dvv.std():.3f} %, median coh {T.coh.median():.3f}")
