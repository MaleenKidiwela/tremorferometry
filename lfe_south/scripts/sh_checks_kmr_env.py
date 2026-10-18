#!/usr/bin/env python
"""SH_CHECKS CHECK 1b: energy-envelope stacks of NC.KMR over the whole -15..+30 s cut window (polarity-, rotation- and
waveform-free).  Is there LFE energy at (or near) the expected S time after each sensor swap?"""
import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_checks_kmr import OUT, FS, IA, NW
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
W = np.load(f"{OUT}/kmr_windows.npz")["W"]; E = pd.read_parquet(f"{OUT}/kmr_windows_events.parquet")
sets = ["OLD_all", "P_OLD", "NK_B", "P_B", "NK_C", "P_C"]
k = int(0.5 * FS); ker = np.ones(k) / k
t = (np.arange(NW) - IA) / FS
out = {}; rows = []
for s in sets:
    envs = []
    for f in [str(i) for i in range(1, 28)]:
        g = E[(E.set == s) & (E.fam == f) & E.ok].sort_values("ccsum", ascending=False).index.values[:100]
        if len(g) < 10: continue
        X = W[g].astype(np.float64)
        e = (X ** 2).sum(1)                                           # (n, 900) 3-component energy
        nz = np.median(e[:, 20:200], 1, keepdims=True) + 1e-30          # each event normalised by its own -14..-5 s level
        e = np.apply_along_axis(lambda v: np.convolve(v, ker, "same"), 1, e / nz)
        envs.append(np.median(e, 0))
    A = np.array(envs); out[s] = A
    m = np.median(A, 0)
    pk = int(np.argmax(m[40:-40])) + 40
    rows.append(dict(set=s, n_fam=len(A), env_at_S_0_3s=float(m[IA:IA + 60].mean()), env_pre_m14_m5=float(m[20:200].mean()),
                     ratio_S_over_pre=float(m[IA:IA + 60].mean() / m[20:200].mean()), peak_time_s=float(t[pk]), peak_over_pre=float(m[pk] / m[20:200].mean()),
                     env_late_p20_p28=float(m[IA + 400:IA + 560].mean() / m[20:200].mean())))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/kmr_envelope_summary.csv", index=False); print(R.round(3).to_string())
f, ax = plt.subplots(1, 1, figsize=(10, 4))
for s, c in zip(sets, ["k", "0.5", "#d62728", "#ff9896", "#2ca02c", "#98df8a"]):
    ax.plot(t, np.median(out[s], 0), color=c, lw=1.1, label=s)
ax.axvline(0, color="0.7", ls=":"); ax.set_xlabel("time after expected S at NC.KMR (s)"); ax.set_ylabel("median energy / own pre-noise (2-8 Hz, 3-C)")
ax.set_title("NC.KMR energy stacks at the 100 strongest detections per family (median over families)", fontsize=10); ax.legend(fontsize=8)
f.tight_layout(); f.savefig(f"{OUT}/kmr_envelope_stacks.png", dpi=130)
