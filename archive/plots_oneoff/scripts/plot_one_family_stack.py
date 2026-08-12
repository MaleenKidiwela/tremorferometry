#!/usr/bin/env python
"""Pick one LFE family in a P(LFE) band from the full 2010-2026 B011 run, stack its
3-comp member waveforms, and plot Z / EH1 / EH2 as SEPARATE panels sharing one time
axis (2-8 Hz LFE band, common amplitude scale), plus a Z-stack spectrum.

Usage: python scripts/plot_one_family_stack.py <LO> <HI>
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from scipy import signal

NET, STA, FS = "PB", "B011", 100.0
PRE, POST = 10.0, 30.0
LO, HI = float(sys.argv[1]), float(sys.argv[2])
PICK = pd.read_csv("data/family_picker_p70_2010_2026_m3.csv")
MEM = pd.read_parquet("data/b011_disc_p70_2010_2026_m3.members.parquet")

band = PICK[(PICK.p_lfe > LO) & (PICK.p_lfe < HI) & (PICK.pred == "LFE")]
fam = band.sample(1, random_state=7).iloc[0]
fid = fam.fam
sub = MEM[MEM.family_id == fid].copy()
sub["t"] = pd.to_datetime(sub.time, utc=True)
exp = int(round((PRE + POST) * FS))
print(f"family {fid}: P(LFE)={fam.p_lfe:.3f}, n_members={len(sub)}, loc {fam.lat:.3f},{fam.lon:.3f}")


def get(st, *ch):
    for c in ch:
        s = st.select(channel=f"*{c}")
        if len(s):
            try: s = s.merge(method=1, fill_value=0)
            except Exception: pass
            return s[0]
    return None


acc = np.zeros((3, exp)); n = 0
for _, r in sub.iterrows():
    d = r.t
    p = f"data/waveforms/{NET}.{STA}/{d.year}/{d.dayofyear:03d}.mseed"
    if not os.path.exists(p):
        continue
    try:
        st = read(p)
    except Exception:
        continue
    z, h1, h2 = get(st, "Z"), get(st, "1", "N"), get(st, "2", "E")
    if z is None or h1 is None or h2 is None:
        continue
    a = UTCDateTime(r.t.value / 1e9)
    try:
        Z = z.slice(a - PRE, a + POST).data.astype(float)[:exp]
        H1 = h1.slice(a - PRE, a + POST).data.astype(float)[:exp]
        H2 = h2.slice(a - PRE, a + POST).data.astype(float)[:exp]
    except Exception:
        continue
    if min(len(Z), len(H1), len(H2)) < exp or np.std(Z) == 0:
        continue
    m = max(np.abs(Z).max(), np.abs(H1).max(), np.abs(H2).max()) + 1e-20
    acc += np.vstack([Z, H1, H2]) / m; n += 1

sos = signal.butter(4, [2.0, 8.0], btype="band", fs=FS, output="sos")
chans = [signal.sosfiltfilt(sos, signal.detrend(w)) for w in (acc / max(n, 1))]
gmax = max(np.abs(c).max() for c in chans) + 1e-20   # common scale -> preserves H/V balance
chans = [c / gmax for c in chans]
tt = np.arange(exp) / FS - PRE

fig = plt.figure(figsize=(11, 9))
gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1.15], hspace=0.12)
axes = [fig.add_subplot(gs[i]) for i in range(3)]
labels = [("Z  (EHZ)", "#1a1a1a"), ("EH1", "#2c7d3f"), ("EH2", "#c0392b")]
for ax, w, (lab, col) in zip(axes, chans, labels):
    ax.plot(tt, w, color=col, lw=0.9)
    ax.axvline(0, color="#888", ls="--", lw=1)
    ax.axvspan(2, 4, color="#f0c000", alpha=0.16)
    ax.set_xlim(-3, 12); ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel(lab, fontsize=10)
    ax.axhline(0, color="#ccc", lw=0.5, zorder=0)
    if ax is not axes[-1]:
        ax.set_xticklabels([])
axes[0].set_title(f"B011 family {fid} — {n}-detection 3-comp stack (2–8 Hz, common scale)\n"
                  f"family-stack P(LFE)={fam.p_lfe:.3f} | loc {fam.lat:.2f},{fam.lon:.2f} | "
                  "dashed=anchor, yellow=2–4 s coda", fontsize=11)
axes[2].set_xlabel("time from detection (s)")

# Z-stack spectrum (broadband, S window) so the band content is visible
axs = fig.add_subplot(gs[3])
w = signal.detrend((acc / max(n, 1))[0][int((PRE - 1) * FS):int((PRE + 2) * FS)])
f, P = signal.welch(w, fs=FS, nperseg=min(256, len(w)))
axs.semilogy(f, P + 1e-20, color="#1a1a1a", lw=1.3)
axs.axvspan(2, 8, color="#2c7d3f", alpha=0.14, label="2–8 Hz LFE band")
axs.set_xlim(0, 25); axs.set_xlabel("frequency (Hz)"); axs.set_ylabel("Z PSD")
axs.legend(fontsize=8.5, loc="upper right")
out = f"figures/b011_family_{fid}_3ch.png"
plt.savefig(out, dpi=160, bbox_inches="tight")
print("saved", out)
