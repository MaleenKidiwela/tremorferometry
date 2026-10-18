#!/usr/bin/env python
"""q7_figfix (2026-10-07 PT) -- redraw panel 3 of q7_strain_stack.png (the first version used a signed-Mw axis that
turns pre-window noise into +/-6 spikes) and the map limits.  Display only: the ETS series are projected with the
P1 control covariances re-estimated here from 500 random 20-day windows per station (same estimator as q7_strain.py);
all numbers in RESULTS come from q7_strain.py, not from this script.

    Q7_DEBUG=1 python -P svi_science/scripts/q7_figfix.py      (Q7_DEBUG only keeps q7_strain's import-time log in /tmp)
"""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q7_strain as Q  # noqa: E402
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q7_strainmeter"
FIG = f"{ROOT}/svi_science/figures"
S = Q.load_stations()
NH = S[0]["Y"].shape[0]
Q.fault_setup(S)
E = Q.load_episodes()
Ex = pd.read_csv(f"{OUT}/positive_control_episodes.csv", parse_dates=["start", "end"])
P1 = Ex[Ex.P == "P1 big ETS"].reset_index(drop=True)
tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["h"] = Q.hour_index(tr.time)
boxes = {}
for reg, f in [("SVI", f"{ROOT}/svi_science/results/s0/families.csv"),
               ("WA", f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv")]:
    F = pd.read_csv(f)
    boxes[reg] = (F.lat.min() - 0.3, F.lat.max() + 0.3, F.lon.min() - 0.4, F.lon.max() + 0.4)
# covariances from random 20-day windows outside major episodes
majors = E[E.major.astype(bool)]
forb = Q.forbidden_mask(NH, [(a - 48, b + 48) for a, b in majors[["h0", "h1"]].values], 0)
D0 = 24 * 20
allowed = Q.window_allowed(forb, Q.W_ETS, D0)
cand = np.flatnonzero(allowed); cand = cand[cand % 24 == 0]
a_r = Q.RNG.choice(cand, 3000)
yrs = (Q.T0 + (a_r * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[Y]").astype(int) + 1970
st = np.full((a_r.size, len(S), 4), np.nan)
for j, s in enumerate(S):
    st[:, j], _ = Q.steps(s["Y"], a_r, D0, Q.W_ETS)
C = Q.cov_by_year(S, [st], [yrs])
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for i, r in P1.iterrows():
    b = boxes[r.region]
    t = tr[(tr.h >= r.a) & (tr.h < r.a + r.D) & tr.lat.between(b[0], b[1]) & tr.lon.between(b[2], b[3])]
    g = Q.template_tremor(t.lat.values, t.lon.values)
    a = np.array([int(r.a)])
    R = dict(st=np.full((1, len(S), 4), np.nan), ok=np.zeros((1, len(S), 4), bool),
             yrs=(Q.T0 + (a * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[Y]").astype(int) + 1970)
    for j, s in enumerate(S):
        R["st"][:, j], R["ok"][:, j] = Q.steps(s["Y"], a, int(r.D), Q.W_ETS)
    rel, y = Q.projected_series(S, R, 0, g, C, int(r.a), int(r.D), Q.W_ETS, np.ones(len(S), bool))
    if np.isfinite(y).all():
        ax[0].plot(rel / 24, y / 1e18, color="C1" if r.region == "WA" else "C0", alpha=0.45, lw=0.8)
ax[0].plot([], [], color="C1", label="WA big ETS"); ax[0].plot([], [], color="C0", label="SVI big ETS")
ax[0].axvline(0, color="k", lw=0.5); ax[0].axhline(0, color="k", lw=0.3)
for mwv in (6.0, 6.5):
    ax[0].axhline(Q.m0_of(mwv) / 1e18, color="0.5", ls=":", lw=0.8)
    ax[0].text(10.2, Q.m0_of(mwv) / 1e18, f"Mw {mwv}", fontsize=7, va="center")
ax[0].set_xlabel("days from episode (episode itself collapsed to 0)")
ax[0].set_ylabel("strain projected on tremor-weighted thrust pattern\n(moment, 10^18 N·m)")
ax[0].set_title("Positive control: each big ETS (n=30) is a strain step")
ax[0].legend(fontsize=8)
# single-episode Mw from strain, from q7_strain.py
for P, c in [("P1 big ETS", "C3"), ("P2 minor", "C2")]:
    d = Ex[Ex.P == P]
    ax[1].scatter(d.start, d.z_single, s=14, color=c, label=P)
ax[1].axhline(2, color="k", lw=0.6, ls=":"); ax[1].axhline(0, color="k", lw=0.3)
ax[1].set_ylabel("single-episode z vs random windows"); ax[1].set_title("Positive control per episode (q7_strain.py)")
ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/q7_positive_control.png", dpi=140); plt.close(fig)
print("saved")
