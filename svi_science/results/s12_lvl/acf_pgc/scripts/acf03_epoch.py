#!/usr/bin/env python
"""S12c §4c: superposed epoch of the differential LVL delay on the ALBH ETS midpoints, with scrambled control.
Usage: acf03_epoch.py --station CN.PGC --years 2018-2021 --band b05_2 --mask tremorfree --tag pilot [--ndraw 2000]
Out:  tables/epoch_<tag>_<STA>_<band>_<mask>.csv     (bin metrics + scrambled pct / sigma)
      tables/curve_<tag>_<STA>_<band>_<mask>.csv     (day-by-day curve, +-0/3/7 d boxcar)
      figures/epoch_<tag>_<STA>_<band>_<mask>.png
Sign: Vs increase in the LVL -> LVL echo earlier -> ddt < 0. pct_neg = fraction of scrambled draws with ddt <= observed.
"""
import sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from acf_common import *

ap = argparse.ArgumentParser()
ap.add_argument("--station", default="CN.PGC"); ap.add_argument("--years", default="2018-2021")
ap.add_argument("--band", default="b05_2"); ap.add_argument("--mask", default="tremorfree")
ap.add_argument("--tag", default="pilot"); ap.add_argument("--ndraw", type=int, default=2000)
ap.add_argument("--ncurve_draw", type=int, default=200); ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
y0, y1 = map(int, a.years.split("-"))
ns = a.station; W = windows_for(ns)
PAIRS = ["ZZ", "RR", "TT", "ZR+", "ZR-", "ZT+", "ZT-"]
t0 = time.time()
D, H, R = load(ns, range(y0, y1 + 1), a.band, pairs=["ZZ", "RR", "TT", "ZR", "ZT"])
ok = hour_mask(ns, D, R, a.mask)
# continuous calendar + cumulative sums for O(1) range stacks
cal = np.arange(D.min(), D.max() + np.timedelta64(1, "D"))
pos = ((D - cal[0]).astype(int))
CS, CN = {}, None
for pr in PAIRS:
    S, n = day_stack(H[pr.rstrip("+-")], ok)
    S = split_pair(S, pr)
    full = np.zeros((len(cal), S.shape[1])); cnt = np.zeros(len(cal))
    g = np.isfinite(S[:, 0]) & (n >= 4)
    full[pos[g]] = S[g]; cnt[pos[g]] = 1
    CS[pr] = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(full, 0)])
    if CN is None:
        CN = np.concatenate([[0], np.cumsum(cnt)])
del H
print(f"loaded {ns} {a.band} {a.mask}: {len(D)} days, {time.time() - t0:.0f} s", flush=True)


def stack(pr, anchors, k0, k1):
    """Mean of daily stacks over days anchor+k0 .. anchor+k1 (inclusive) for all anchors. Returns (trace, ndays)."""
    tot = 0; nd = 0
    for an in anchors:
        i0 = int((an - cal[0]).astype(int)) + k0; i1 = int((an - cal[0]).astype(int)) + k1 + 1
        i0 = max(i0, 0); i1 = min(i1, len(cal))
        if i1 <= i0: continue
        tot = tot + (CS[pr][i1] - CS[pr][i0]); nd += CN[i1] - CN[i0]
    return (tot / nd if nd > 0 else None), nd


PHASE = {"ZZ": "PP", "RR": "SS", "TT": "SS", "ZR+": "PS", "ZR-": "PS", "ZT+": "PS", "ZT-": "PS"}
BINS = {"post0_15": (0, 15), "post15_30": (15, 30), "post30_60": (30, 60), "post60_120": (60, 120)}
REF = (-30, -1)


def metrics(anchors):
    out = {}
    for pr in PAIRS:
        ref, nr = stack(pr, anchors, *REF)
        for bn, (k0, k1) in BINS.items():
            cur, nc = stack(pr, anchors, k0, k1)
            if ref is None or cur is None or nr < 5 or nc < 5:
                continue
            m = differential(ref, cur, W[PHASE[pr]])
            wl = W[PHASE[pr]]["lvl"]; mm = (LAG_AUTO >= wl[0]) & (LAG_AUTO < wl[1])
            m["amp_lvl"] = float(np.sqrt(np.mean(cur[mm] ** 2)) / np.sqrt(np.mean(ref[mm] ** 2)))
            for k, v in m.items():
                out[(pr, bn, k)] = v
            out[(pr, bn, "nref")] = nr; out[(pr, bn, "ncur")] = nc
    for bn in BINS:  # SS = mean of RR, TT (PLAN §4c); PS = mean of ZR+, ZR-
        for grp, mem in [("SS", ["RR", "TT"]), ("PSr", ["ZR+", "ZR-"])]:
            for k in ["ddt", "ddt_str", "ddt_ccpk", "dcc", "ddt_ctl", "amp_lvl", "cc_lvl"]:
                v = [out.get((p, bn, k), np.nan) for p in mem]
                out[(grp, bn, k)] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
    return out


mids = ets_mid()
mids = mids[(mids >= cal[0] + np.timedelta64(30, "D")) & (mids <= cal[-1] - np.timedelta64(15, "D"))]
print("episodes:", [str(m) for m in mids], flush=True)
obs = metrics(mids)

# scrambled control: same day of year +-30 d, random year with data, >60 d from every real midpoint
rng = np.random.default_rng(a.seed)
allm = ets_mid()
yrs = np.unique(D.astype("datetime64[Y]"))


def scramble():
    out = []
    for m in mids:
        for _ in range(200):
            y = rng.choice(yrs)
            doy = (m - m.astype("datetime64[Y]")).astype(int) + rng.integers(-30, 31)
            c = y.astype("datetime64[D]") + np.timedelta64(int(doy), "D")
            if np.abs((allm - c).astype(int)).min() > 60 and cal[0] <= c <= cal[-1]:
                out.append(c); break
    return np.array(out)


keys = [k for k in obs if k[2] in ("ddt", "ddt_str", "ddt_ccpk", "dcc", "ddt_ctl", "amp_lvl", "cc_lvl")]
SC = {k: [] for k in keys}
for i in range(a.ndraw):
    mm = metrics(scramble())
    for k in keys:
        SC[k].append(mm.get(k, np.nan))
rows = []
for k in keys:
    v = np.array(SC[k], float); v = v[np.isfinite(v)]
    o = obs[k]
    rows.append(dict(pair=k[0], bin=k[1], metric=k[2], obs=o, nref=obs.get((k[0], k[1], "nref")),
                     ncur=obs.get((k[0], k[1], "ncur")), scr_mean=v.mean() if len(v) else np.nan,
                     scr_sd=v.std() if len(v) else np.nan,
                     scr_rsd=1.4826 * np.median(np.abs(v - np.median(v))) if len(v) else np.nan, scr_p05=np.percentile(v, 5) if len(v) else np.nan,
                     scr_p95=np.percentile(v, 95) if len(v) else np.nan,
                     pct_neg=float((v <= o).mean()) if len(v) else np.nan, n_draw=len(v)))
T = pd.DataFrame(rows)
np.savez(f"{WS}/tables/epoch_draws_{a.tag}_{ns}_{a.band}_{a.mask}.npz", keys=np.array(["|".join(k) for k in keys]),
         draws=np.array([SC[k] for k in keys], float), obs=np.array([obs[k] for k in keys], float))
fo = f"{WS}/tables/epoch_{a.tag}_{ns}_{a.band}_{a.mask}.csv"; T.to_csv(fo, index=False)
print(f"scrambled done {time.time() - t0:.0f} s", flush=True)
pd.set_option("display.width", 220)
print(T[(T.metric.isin(["ddt", "ddt_str", "ddt_ctl", "dcc"])) & (T.bin == "post0_15")].round(4).to_string(index=False))

# day-by-day curve (reference = days -60..-10), boxcar half-widths 0, 3, 7
crow = []
for hw in [0, 3, 7]:
    for pr in ["RR", "TT", "ZZ", "ZR+", "ZR-"]:
        ref, _ = stack(pr, mids, -60, -10)
        for k in range(-60, 121):
            cur, nc = stack(pr, mids, k - hw, k + hw)
            if cur is None or nc < max(2, len(mids) // 2):
                continue
            m = differential(ref, cur, W[PHASE[pr]])
            crow.append(dict(hw=hw, pair=pr, k=k, ncur=nc, **m))
C = pd.DataFrame(crow)
C.to_csv(f"{WS}/tables/curve_{a.tag}_{ns}_{a.band}_{a.mask}.csv", index=False)
# scrambled curve band (SS mean, hw 3)
scr_curves = []
for i in range(a.ncurve_draw):
    an = scramble(); row = []
    refs = {pr: stack(pr, an, -60, -10)[0] for pr in ["RR", "TT"]}
    for k in range(-60, 121, 3):
        v = []
        for pr in ["RR", "TT"]:
            cur, nc = stack(pr, an, k - 3, k + 3)
            if cur is not None and refs[pr] is not None and nc >= 2:
                v.append(differential(refs[pr], cur, W["SS"])["ddt"])
        row.append(np.nanmean(v) if v else np.nan)
    scr_curves.append(row)
scr_curves = np.array(scr_curves); kk = np.arange(-60, 121, 3)
np.save(f"{WS}/tables/curve_scr_{a.tag}_{ns}_{a.band}_{a.mask}.npy", scr_curves)

# figure: SS (RR,TT mean) ddt per day, hw=3, coloured by CC_lvl; scrambled 5/95 % thin grey lines
fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
for ax, (lab, mem, ph) in zip(axs, [("SS (RR,TT)", ["RR", "TT"], "SS"), ("PS (ZR±)", ["ZR+", "ZR-"], "PS"), ("PP (ZZ)", ["ZZ"], "PP")]):
    c = C[(C.hw == 3) & (C.pair.isin(mem))].groupby("k")[["ddt", "cc_lvl", "ddt_ctl"]].mean().reset_index()
    ax.plot(c.k, 1e3 * c.ddt, color="0.6", lw=0.8, zorder=1)
    sc = ax.scatter(c.k, 1e3 * c.ddt, c=c.cc_lvl, cmap="viridis", s=14, zorder=2, vmin=0.3, vmax=1)
    ax.plot(c.k, 1e3 * c.ddt_ctl, color="tab:orange", lw=0.8, label="control (two windows above the LVL)")
    if ph == "SS":
        ax.plot(kk, 1e3 * np.nanpercentile(scr_curves, 5, 0), color="0.3", lw=0.6, ls="--", label="scrambled 5/95 %")
        ax.plot(kk, 1e3 * np.nanpercentile(scr_curves, 95, 0), color="0.3", lw=0.6, ls="--")
    tg = {"SS": -58, "PS": -34, "PP": -12}[ph]
    ax.axhline(tg, color="tab:red", lw=0.8, ls=":", label=f"target {tg} ms (Gosselin size)")
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5); ax.axvspan(0, 15, color="tab:red", alpha=0.06)
    ax.set_ylabel(f"{lab}\nΔt LVL − pre (ms)"); plt.colorbar(sc, ax=ax, label="CC (LVL window)")
    ax.legend(fontsize=7, loc="lower right")
axs[-1].set_xlabel("days after ETS midpoint (ALBH GPS)")
fig.suptitle(f"{ns} {a.band} {a.mask} {a.years}: differential LVL delay, ±3-day boxcar, {len(mids)} ETS; ref = days −60…−10", fontsize=9)
fig.tight_layout(); fig.savefig(f"{WS}/figures/epoch_{a.tag}_{ns}_{a.band}_{a.mask}.png", dpi=110)
print(f"done {time.time() - t0:.0f} s", flush=True)
