#!/usr/bin/env python
"""S12c PLAN amendment B: LVL echo AMPLITUDE (relative reflection strength) and waveform similarity, superposed epoch on
the ALBH ETS midpoints with the scrambled control.
Usage: acf05_amp.py --station CN.PGC --years 2003-2026 --band b05_2 --mask tremorfree --tag full [--ndraw 2000]
Out:  tables/amp_<tag>_<STA>_<band>_<mask>.csv (bins), tables/ampcurve_<tag>_..csv (±3-d curve),
      tables/ampcurve_scr_<tag>_..npy, figures/amp_<tag>_<STA>_<band>_<mask>.png
"""
import sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from acf_common import *

ap = argparse.ArgumentParser()
ap.add_argument("--station", default="CN.PGC"); ap.add_argument("--years", default="2003-2026")
ap.add_argument("--band", default="b05_2"); ap.add_argument("--mask", default="tremorfree")
ap.add_argument("--tag", default="full"); ap.add_argument("--ndraw", type=int, default=2000)
ap.add_argument("--ncurve_draw", type=int, default=200); ap.add_argument("--seed", type=int, default=2)
a = ap.parse_args()
y0, y1 = map(int, a.years.split("-")); ns = a.station; W = windows_for(ns)
B = pd.read_csv(f"{WS}/tables/predicted_lag_bands.csv")
PAIRS = ["ZZ", "RR", "TT", "ZR+", "ZR-"]
PHASE = {"ZZ": "PP", "RR": "SS", "TT": "SS", "ZR+": "PS", "ZR-": "PS"}
GROUP = {"SS": ["RR", "TT"], "PS": ["ZR+", "ZR-"], "PP": ["ZZ"]}
t0 = time.time()
D, H, R = load(ns, range(y0, y1 + 1), a.band, pairs=["ZZ", "RR", "TT", "ZR"])
ok = hour_mask(ns, D, R, a.mask)
cal = np.arange(D.min(), D.max() + np.timedelta64(1, "D")); pos = (D - cal[0]).astype(int)
CS, CN, LT = {}, None, {}
for pr in PAIRS:
    S, n = day_stack(H[pr.rstrip("+-")], ok); S = split_pair(S, pr)
    full = np.zeros((len(cal), S.shape[1])); cnt = np.zeros(len(cal))
    g = np.isfinite(S[:, 0]) & (n >= 4); full[pos[g]] = S[g]; cnt[pos[g]] = 1
    CS[pr] = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(full, 0)])
    LT[pr] = S[g].mean(0)                                   # long-term echo = all days
    if CN is None: CN = np.concatenate([[0], np.cumsum(cnt)])
del H
print(f"loaded {ns} {a.band} {a.mask}: {len(D)} days, {time.time() - t0:.0f} s", flush=True)

t = LAG_AUTO
WIN = {}
for ph in ["SS", "PS", "PP"]:
    for iface in ["top", "bottom"]:
        r = B[(B.station == ns) & (B.interface == iface)].iloc[0]
        WIN[(ph, iface)] = (t >= r[f"{ph}_lo"] - 0.3) & (t <= r[f"{ph}_hi"] + 0.3)
    WIN[(ph, "ref")] = (t >= W[ph]["pre"][0]) & (t < W[ph]["pre"][1])


def stack(pr, anchors, k0, k1):
    tot = 0; nd = 0
    for an in anchors:
        b = int((an - cal[0]).astype(int)); i0 = max(b + k0, 0); i1 = min(b + k1 + 1, len(cal))
        if i1 <= i0: continue
        tot = tot + (CS[pr][i1] - CS[pr][i0]); nd += CN[i1] - CN[i0]
    return (tot / nd if nd > 0 else None), nd


def feats(s, pr):
    """Per echo window: relative projection amplitude, relative RMS, CC to long-term echo, CC difference to ref window."""
    L = LT[pr]; ph = PHASE[pr]; o = {}
    mr = WIN[(ph, "ref")]
    a_ref = np.dot(s[mr], L[mr]) / np.dot(L[mr], L[mr]); r_ref = np.sqrt(np.mean(s[mr] ** 2))
    cc_ref = np.corrcoef(s[mr], L[mr])[0, 1]
    for iface in ["top", "bottom"]:
        m = WIN[(ph, iface)]
        o[f"A_{iface}"] = (np.dot(s[m], L[m]) / np.dot(L[m], L[m])) / a_ref
        o[f"rms_{iface}"] = np.sqrt(np.mean(s[m] ** 2)) / r_ref
        o[f"cc_{iface}"] = np.corrcoef(s[m], L[m])[0, 1]
        o[f"dcc_{iface}"] = o[f"cc_{iface}"] - cc_ref
    return o


BINS = {"post0_15": (0, 15), "post15_30": (15, 30), "post30_60": (30, 60), "post60_120": (60, 120)}
REF = (-30, -1)
MET = [f"{q}_{i}" for i in ["top", "bottom"] for q in ["A", "rms", "cc", "dcc"]]


def metrics(anchors):
    out = {}
    for pr in PAIRS:
        ref, nr = stack(pr, anchors, *REF)
        if ref is None or nr < 5: continue
        fr = feats(ref, pr)
        for bn, (k0, k1) in BINS.items():
            cur, nc = stack(pr, anchors, k0, k1)
            if cur is None or nc < 5: continue
            fc = feats(cur, pr)
            for q in MET:
                out[(pr, bn, q)] = fc[q] / fr[q] - 1 if q.startswith(("A_", "rms_")) else fc[q] - fr[q]
    for ph, mem in GROUP.items():
        for bn in BINS:
            for q in MET:
                v = [out.get((p, bn, q), np.nan) for p in mem]
                out[(ph, bn, q)] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
    return out


mids = ets_mid()
mids = mids[(mids >= cal[0] + np.timedelta64(30, "D")) & (mids <= cal[-1] - np.timedelta64(15, "D"))]
print("episodes:", len(mids), flush=True)
obs = metrics(mids)
rng = np.random.default_rng(a.seed); allm = ets_mid(); yrs = np.unique(D.astype("datetime64[Y]"))


def scramble():
    out = []
    for m in mids:
        for _ in range(200):
            y = rng.choice(yrs); doy = (m - m.astype("datetime64[Y]")).astype(int) + rng.integers(-30, 31)
            c = y.astype("datetime64[D]") + np.timedelta64(int(doy), "D")
            if np.abs((allm - c).astype(int)).min() > 60 and cal[0] <= c <= cal[-1]:
                out.append(c); break
    return np.array(out)


keys = list(obs)
SC = {k: [] for k in keys}
for i in range(a.ndraw):
    mm = metrics(scramble())
    for k in keys: SC[k].append(mm.get(k, np.nan))
rows = []
for k in keys:
    v = np.array(SC[k], float); v = v[np.isfinite(v)]; o = obs[k]
    med = np.median(v) if len(v) else np.nan
    rows.append(dict(pair=k[0], bin=k[1], metric=k[2], obs=o, scr_median=med,
                     scr_rsd=1.4826 * np.median(np.abs(v - med)) if len(v) else np.nan,
                     scr_p025=np.percentile(v, 2.5) if len(v) else np.nan, scr_p975=np.percentile(v, 97.5) if len(v) else np.nan,
                     pct_low=float((v <= o).mean()) if len(v) else np.nan,
                     pct_two=float(min(1, 2 * min((v <= o).mean(), (v >= o).mean()))) if len(v) else np.nan, n_draw=len(v)))
T = pd.DataFrame(rows); T.to_csv(f"{WS}/tables/amp_{a.tag}_{ns}_{a.band}_{a.mask}.csv", index=False)
pd.set_option("display.width", 220)
print(T[(T.pair.isin(["SS", "PS", "PP"])) & (T.bin == "post0_15") & T.metric.isin(["A_top", "A_bottom", "cc_top", "cc_bottom"])]
      .round(4).to_string(index=False), flush=True)

# ±3-d curve vs days -60..-10 (SS, PS, PP group means) + scrambled band for SS
def curve(anchors, ks, hw=3):
    res = {}
    for ph, mem in GROUP.items():
        refs = {pr: stack(pr, anchors, -60, -10)[0] for pr in mem}
        for k in ks:
            vals = {q: [] for q in MET}
            for pr in mem:
                cur, nc = stack(pr, anchors, k - hw, k + hw)
                if cur is None or refs[pr] is None or nc < 2: continue
                fc, fr = feats(cur, pr), feats(refs[pr], pr)
                for q in MET:
                    vals[q].append(fc[q] / fr[q] - 1 if q.startswith(("A_", "rms_")) else fc[q] - fr[q])
            for q in MET:
                res[(ph, k, q)] = np.mean(vals[q]) if vals[q] else np.nan
    return res
ks = np.arange(-60, 121)
cv = curve(mids, ks)
C = pd.DataFrame([dict(phase=p, k=k, metric=q, value=v) for (p, k, q), v in cv.items()])
C.to_csv(f"{WS}/tables/ampcurve_{a.tag}_{ns}_{a.band}_{a.mask}.csv", index=False)
kk = np.arange(-60, 121, 3); SCC = []
for i in range(a.ncurve_draw):
    r = curve(scramble(), kk)
    SCC.append([[r[("SS", k, q)] for k in kk] for q in ["A_top", "A_bottom", "cc_top", "cc_bottom"]])
SCC = np.array(SCC); np.save(f"{WS}/tables/ampcurve_scr_{a.tag}_{ns}_{a.band}_{a.mask}.npy", SCC)

fig, axs = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
for j, (ax, q, lab) in enumerate(zip(axs, ["A_top", "A_bottom", "cc_top", "cc_bottom"],
                                      ["S–S top echo strength change", "S–S bottom echo strength change",
                                       "S–S top echo CC change", "S–S bottom echo CC change"])):
    c = C[(C.phase == "SS") & (C.metric == q)]
    ax.plot(c.k, c.value, color="0.3", lw=1.0, label="S–S (RR, TT)")
    for ph, col in [("PS", "tab:blue"), ("PP", "tab:green")]:
        cc_ = C[(C.phase == ph) & (C.metric == q)]
        ax.plot(cc_.k, cc_.value, color=col, lw=0.7, alpha=0.7, label=ph)
    ax.plot(kk, np.nanpercentile(SCC[:, j], 2.5, 0), color="0.3", lw=0.5, ls="--", label="S–S scrambled 2.5/97.5 %")
    ax.plot(kk, np.nanpercentile(SCC[:, j], 97.5, 0), color="0.3", lw=0.5, ls="--")
    ax.axhline(0, color="k", lw=0.4); ax.axvline(0, color="k", lw=0.5); ax.axvspan(0, 15, color="tab:red", alpha=0.05)
    ax.set_ylabel(lab, fontsize=8); ax.legend(fontsize=6, loc="lower right")
axs[-1].set_xlabel("days after ETS midpoint (ALBH GPS)")
fig.suptitle(f"{ns} {a.band} {a.mask} {a.years}: LVL echo amplitude / similarity, ±3-d boxcar, {len(mids)} ETS; ref = days −60…−10", fontsize=9)
fig.tight_layout(); fig.savefig(f"{WS}/figures/amp_{a.tag}_{ns}_{a.band}_{a.mask}.png", dpi=110)
print(f"done {time.time() - t0:.0f} s", flush=True)
