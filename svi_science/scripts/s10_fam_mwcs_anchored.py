"""MWCS, S-anchored fit (dt = b * lapse after picked S, no intercept), on the saved SVD-Wiener 3-C day stacks of
family 10 at PGC (tables/fam10_PGC_svd3c_<epoch>.npz, built by s10_fam_mwcs.py). 2-5 s, 1-s windows, 0.25-s step,
2-8 Hz coherence-weighted phase slope. User 2026-10-12. Plot: no error bars, colour = coherence."""
import glob, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; W0, W1 = 2.0, 5.0
WL, STEP, NFFT, FMIN, FMAX, SM = 1.0, 0.25, 256, 2.0, 8.0, 5
out = []
for fn in sorted(glob.glob(f"{PIL}/tables/fam10_PGC_svd3c_*.npz")):
    Z = np.load(fn); tau = Z["tau"]; fs = 1 / (tau[1] - tau[0]); tS = float(Z["tS"]); Ga = Z["Ga"]
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
        b = np.sum(Wt ** 2 * T * D) / np.sum(Wt ** 2 * T * T)                       # dt = b * lapse (S-anchored)
        return -b * 100, float(np.mean(C))
    for k, d in enumerate(pd.to_datetime(Z["days"])):
        dv, coh = mwcs(Z["Df"][k].astype(np.float64), Ga, tS); out.append(dict(epoch=fn.split("_")[-1][:-4], date=d, dvv_mwcs=dv, coh=coh))
T = pd.DataFrame(out).sort_values("date"); T.to_csv(f"{PIL}/tables/fam10_PGC_mwcs_anchored_svd3c25_daily.csv", index=False)
K = T[(T.coh > 0.8) & (T.dvv_mwcs.abs() < 1)].copy()
dd = np.sort(K.date.unique()); epn = np.r_[0, np.cumsum(np.diff(dd).astype("timedelta64[D]").astype(int) > 20)]; K["ep"] = K.date.map(dict(zip(dd, epn)))
E = K.groupby("ep").apply(lambda q: pd.Series(dict(date=q.date.min(), dvv=np.average(q.dvv_mwcs, weights=q.coh), coh=q.coh.mean(), nd=len(q))), include_groups=False)
fig, ax = plt.subplots(figsize=(13, 4.5))
ax.scatter(K.date, K.dvv_mwcs, c=K.coh, cmap="viridis", vmin=0.8, vmax=1, s=12, alpha=0.6, lw=0)
ax.plot(E.date, E.dvv, "-", color="k", lw=1.2, zorder=3)
sc = ax.scatter(E.date, E.dvv, c=E.coh, cmap="viridis", vmin=0.8, vmax=1, s=25 + 10 * E.nd, edgecolor="k", lw=0.6, zorder=4)
ax.axhline(0, color="grey", lw=0.6, ls=":"); ax.set_ylabel("dv/v (%)"); fig.colorbar(sc, ax=ax, label="mean coherence (2–8 Hz)")
ax.set_title(f"Family 10 at PGC: dv/v by MWCS, S-anchored fit (zero delay at S), SVD-Wiener 3-component day stacks, 2–5 s; small = days (coh > 0.8, {len(K)}), big + line = per-episode mean\n"
             f"day sd {K.dvv_mwcs.std():.3f} %, episode sd {E.dvv.std():.3f} %", fontsize=9)
fig.tight_layout(); outp = f"{PIL}/figures/s10_fam10_PGC_mwcs_anchored_svd3c25.png"; fig.savefig(outp, dpi=120); print(outp)
print(f"days {len(T)}, coh>0.8 {len(K)}, day sd {K.dvv_mwcs.std():.3f} %, episode sd {E.dvv.std():.3f} %")
