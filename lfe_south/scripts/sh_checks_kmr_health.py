#!/usr/bin/env python
"""SH_CHECKS CHECK 1c: NC.KMR raw-data health per month (from the cut windows) + rate effect of KMR after the swaps."""
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_checks_kmr import load_month_kmr, CH, OUT, EXT, RES
from scipy.signal import butter, sosfiltfilt, welch
from sh_v3_common import apply_rules
SEN = lambda m: 629145000.0 if m < "2025-09-18" else (504678000.0 if m < "2026-03-13" else 316573000.0)
months = pd.period_range("2023-12", "2026-09", freq="M").astype(str)
rows = []; psd = {}
sb = butter(4, [2 / 50, 8 / 50], "band", output="sos"); sl = butter(2, [0.1 / 50, 0.5 / 50], "band", output="sos")
for m in months:
    r = load_month_kmr(m)
    if r is None: continue
    segs, dat = r
    for c in CH:
        G = segs[c]; G = G[G.n >= 4500]
        G = G.sample(n=min(150, len(G)), random_state=0)
        st = []; P = []
        for g in G.itertuples():
            x = dat[c][int(g.h5_start):int(g.h5_start) + 4500].astype(float)
            if (x == 0).mean() > 0.01: continue
            st.append((x.mean(), x.std(), len(np.unique(x)), sosfiltfilt(sb, x - x.mean()).std(), sosfiltfilt(sl, x - x.mean()).std()))
            f, p = welch(x - x.mean(), fs=100, nperseg=1024); P.append(p)
        if not st: continue
        a = np.median(np.array(st), 0); s = SEN(m)
        rows.append(dict(month=m, cha=c, n=len(st), mean_counts=a[0], std_counts=a[1], n_unique=a[2], rms_2_8Hz_nm_s=a[3] / s * 1e9,
                         rms_0p1_0p5Hz_nm_s=a[4] / s * 1e9, frac_flat=float(np.mean([q[1] == 0 for q in st]))))
        psd[(m, c)] = np.median(np.array(P), 0) / s ** 2
H = pd.DataFrame(rows); H.to_csv(f"{OUT}/kmr_health_monthly.csv", index=False)
pd.set_option("display.width", 250)
print(H.pivot(index="month", columns="cha", values="rms_0p1_0p5Hz_nm_s").round(1).join(H.pivot(index="month", columns="cha", values="rms_2_8Hz_nm_s").round(1), rsuffix="_2-8Hz").join(
      H.pivot(index="month", columns="cha", values="mean_counts").round(-3), rsuffix="_mean").to_string())
# PSD peak in 2-8 Hz, Horizon epoch
f = np.fft.rfftfreq(1024, 0.01)
for m in ("2024-06", "2025-10", "2026-01", "2026-05"):
    for c in CH:
        if (m, c) in psd:
            p = psd[(m, c)]; b = (f >= 2) & (f <= 8); j = np.argmax(p[b])
            print(m, c, "2-8 Hz PSD peak at %.2f Hz, peak/median %.1f" % (f[b][j], p[b][j] / np.median(p[b])))
np.savez_compressed(f"{OUT}/kmr_psd_monthly.npz", f=f, **{f"{k[0]}|{k[1]}": v for k, v in psd.items()})
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fg, ax = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
for a, c in zip(ax, CH):
    for m, col in (("2024-06", "k"), ("2025-06", "0.5"), ("2025-08", "#9467bd"), ("2025-12", "#d62728"), ("2026-06", "#2ca02c")):
        if (m, c) in psd: a.loglog(f[1:], np.sqrt(psd[(m, c)][1:]) * 1e9, color=col, lw=0.9, label=m)
    a.set_title(f"NC.KMR {c} (median of cut windows, StationXML gain)", fontsize=9); a.set_xlabel("Hz"); a.axvspan(2, 8, color="0.92", zorder=0)
ax[0].set_ylabel("amplitude spectral density, nm/s/√Hz"); ax[0].legend(fontsize=7)
fg.tight_layout(); fg.savefig(f"{OUT}/kmr_spectra.png", dpi=120)
# ---- rate effect: S1 9x counts with vs without KMR on the same days (B and C epochs)
rows = []
for ep, seg in (("B", "B"), ("C", "C")):
    NK = pd.read_parquet(f"{EXT}/cat_NK_{ep}.parquet"); CAT = pd.read_parquet(f"{EXT}/cat_S1_9x_internal.parquet")
    K = CAT[CAT.seg == seg]; days = sorted(set(K.t.dt.floor("D")) | set(NK.t.dt.floor("D")))
    rows.append(dict(epoch=ep, days=len(days), S1_without_KMR=len(NK), S1_with_new_KMR=len(K), ratio_with_over_without=len(K) / len(NK),
                     meancc_med_without=NK.meancc.median(), meancc_med_with=K.meancc.median(), nsta_med_without=NK.nsta.median(), nsta_med_with=K.nsta.median()))
    NK["fam"] = NK.fam.astype(str); K = K.assign(fam=K.fam.astype(str))
    pf = (K.groupby("fam").size() / NK.groupby("fam").size()); rows[-1].update(per_family_ratio_median=pf.median(), per_family_ratio_min=pf.min(), per_family_ratio_max=pf.max())
    # event overlap within 0.5 s
    a = np.sort(NK.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9); b = K.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    j = np.clip(np.searchsorted(a, b), 1, len(a) - 1); d = np.minimum(np.abs(a[j] - b), np.abs(a[j - 1] - b))
    rows[-1].update(frac_with_also_in_without_0p5s=float((d <= 0.5).mean()))
RR = pd.DataFrame(rows); RR.to_csv(f"{OUT}/kmr_rate_effect.csv", index=False); print(RR.round(3).T.to_string())
