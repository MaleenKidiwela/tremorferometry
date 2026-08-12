#!/usr/bin/env python
"""Illustrate why reliable families are separable WITHOUT the reversed template:
the forward grand stack itself shows real causal coda (reliable) vs a symmetric blob (ringing).
Pick a high-causality (reliable) and a ~1 causality (ringing) family from B926 and plot both,
with the MIRROR (-2..0s, pre-arrival = noise floor) and CODA (2-4s) windows shaded."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

d = np.load("data/long_window_daily_B926p90f40_Z.npz", allow_pickle=True)
t, p, nd, S = d["t"], d["patches"], d["n_det"].astype(float), d["stacks"]
coda = (t >= 2) & (t <= 4); mir = (t >= -2) & (t < 0)
rms = lambda x: float(np.sqrt(np.mean(x**2)))


def grand(fam):
    m = p == fam; w = nd[m]
    return (S[m] * w[:, None]).sum(0) / w.sum()


rows = []
for fam in pd.unique(p):
    g = grand(fam); c, mm = rms(g[coda]), rms(g[mir])
    rows.append((fam, c / mm if mm else np.nan))
R = pd.DataFrame(rows, columns=["fam", "caus"]).dropna().sort_values("caus")
ringing = R.iloc[(R.caus - 1.0).abs().argmin()].fam          # causality ~ 1
reliable = R.iloc[-1].fam                                     # highest causality

fig, axes = plt.subplots(1, 2, figsize=(14, 4.2), sharey=False)
for ax, fam, kind in [(axes[0], reliable, "RELIABLE"), (axes[1], ringing, "RINGING")]:
    g = grand(fam); g = g / np.abs(g).max()
    c, mm = rms(g[coda]), rms(g[mir])
    ax.axvspan(-2, 0, color="#3b6ea5", alpha=0.16, label="mirror -2..0 s (pre-arrival = noise floor)")
    ax.axvspan(2, 4, color="#c1440e", alpha=0.16, label="coda 2..4 s")
    ax.axvline(0, color="0.4", lw=0.8, ls="--")
    ax.plot(t, g, color="#1a1a2e", lw=0.8)
    ax.set_title(f"{kind}: {fam}\ncoda/mirror = {c/mm:.2f}   (>1.5 = real causal coda)")
    ax.set_xlabel("lapse time (s), 0 = detection"); ax.set_xlim(-3, 8)
    ax.legend(fontsize=7, loc="upper right", frameon=False)
axes[0].set_ylabel("normalized amplitude")
fig.suptitle("Why no reversed template is needed: the forward stack shows causal coda (left) vs a symmetric blob (right)", fontsize=11)
fig.tight_layout()
fig.savefig("figures/causality_illustration_b926.png", dpi=130)
print(f"reliable {reliable} caus {R.iloc[-1].caus:.2f} | ringing {ringing} caus {R[R.fam==ringing].caus.iloc[0]:.2f}")
print("-> figures/causality_illustration_b926.png")
