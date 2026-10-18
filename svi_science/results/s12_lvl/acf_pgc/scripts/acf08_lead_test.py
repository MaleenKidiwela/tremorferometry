#!/usr/bin/env python
"""S12c acf08_lead_test (written Wed 2026-10-14 ~3:00 pm PT; PLAN §10 Amendment C, C1, written before running).
Confirming test of the post-hoc "weeks before ETS" lead on data not used to find it (VGZ; SNB descriptive). Copy of
acf07 with two changes: (1) bins = the lead contrast (days -21..-1 vs -60..-25) + the pre-set 0..15 vs -30..-1 as a
reproduction check; (2) SAME-MONTH scrambled control: each midpoint moved to another year with data at the same calendar
date +-U{-15..15} d, rejected if < 60 d from any real midpoint, < 30 d from a date already drawn, or outside the record.
Original acf07 docstring follows.
---
S12c POST HOC (written Wed 2026-10-14 ~1:40 pm PT, after the pre-set full-record results were known).
Compares the ACF with the RF recovery window from rf_gosselin RESULTS §6/§11 and target.json: at PGC the RF change lasts
~0-24 d after the ALBH midpoint, lowest at 6-15 d, and ~half of the pre/post contrast is a HIGH 2-3 weeks before.
Same data, hour mask, windows, delay estimator (acf_common.differential), echo-strength features (as acf05_amp.feats)
and scrambled-control rule as acf03/acf05; only the day bins differ. Not part of the pre-set test.
Usage: acf07_posthoc_rfwindow.py --station CN.PGC --band b05_2 --mask tremorfree [--ndraw 2000]
Out:  tables/posthoc_rfwindow_<STA>_<band>_<mask>.csv   [acf08: tables/leadtest_<STA>_<band>_<mask>.csv]
"""
import sys, os, argparse, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from acf_common import *

ap = argparse.ArgumentParser()
ap.add_argument("--station", default="CN.PGC"); ap.add_argument("--years", default="2003-2026")
ap.add_argument("--band", default="b05_2"); ap.add_argument("--mask", default="tremorfree")
ap.add_argument("--ndraw", type=int, default=2000); ap.add_argument("--seed", type=int, default=7)
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
    LT[pr] = S[g].mean(0)
    if CN is None: CN = np.concatenate([[0], np.cumsum(cnt)])
del H
print(f"loaded {ns} {a.band} {a.mask}: {len(D)} days, {time.time() - t0:.0f} s", flush=True)

t = LAG_AUTO; WIN = {}
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
    L = LT[pr]; ph = PHASE[pr]; o = {}; mr = WIN[(ph, "ref")]
    a_ref = np.dot(s[mr], L[mr]) / np.dot(L[mr], L[mr])
    for iface in ["top", "bottom"]:
        m = WIN[(ph, iface)]
        o[f"A_{iface}"] = (np.dot(s[m], L[m]) / np.dot(L[m], L[m])) / a_ref
        o[f"cc_{iface}"] = np.corrcoef(s[m], L[m])[0, 1]
    return o


# (current bin, reference bin)
CMP = {"m21_m1_vs_m60_m25 (LEAD TEST)": ((-21, -1), (-60, -25)),
       "0_15_vs_m30_m1 (pre-set, repro check)": ((0, 15), (-30, -1))}
TM = ["ddt", "dt_pre", "dt_lvl", "ddt_ctl"]; AM = ["A_top", "A_bottom", "cc_top", "cc_bottom"]


def metrics(anchors):
    out = {}
    for pr in PAIRS:
        for cn, (cb, rb) in CMP.items():
            ref, nr = stack(pr, anchors, *rb); cur, nc = stack(pr, anchors, *cb)
            if ref is None or cur is None or nr < 5 or nc < 5: continue
            m = differential(ref, cur, W[PHASE[pr]])
            for q in TM: out[(pr, cn, q)] = m[q]
            fr, fc = feats(ref, pr), feats(cur, pr)
            for q in AM: out[(pr, cn, q)] = fc[q] / fr[q] - 1 if q.startswith("A_") else fc[q] - fr[q]
    for ph, mem in GROUP.items():
        for cn in CMP:
            for q in TM + AM:
                v = [out.get((p, cn, q), np.nan) for p in mem]
                out[(ph, cn, q)] = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
    return out


mids = ets_mid()
mids = mids[(mids >= cal[0] + np.timedelta64(30, "D")) & (mids <= cal[-1] - np.timedelta64(15, "D"))]
print("episodes:", len(mids), flush=True)
obs = metrics(mids)
rng = np.random.default_rng(a.seed); allm = ets_mid(); yrs = np.unique(D.astype("datetime64[Y]"))


def scramble():
    out = []
    lo, hi = cal[0] + np.timedelta64(30, "D"), cal[-1] - np.timedelta64(15, "D")
    for m in mids:
        for _ in range(200):
            y = rng.choice(yrs); doy = (m - m.astype("datetime64[Y]")).astype(int) + rng.integers(-15, 16)
            c = y.astype("datetime64[D]") + np.timedelta64(int(doy), "D")
            if not (lo <= c <= hi) or np.abs((allm - c).astype(int)).min() < 60:
                continue
            if out and np.abs((np.array(out) - c).astype(int)).min() < 30:
                continue
            out.append(c); break
    return np.array(out, dtype="datetime64[D]")


keys = [k for k in obs if k[0] in GROUP]
SC = {k: [] for k in keys}
for i in range(a.ndraw):
    mm = metrics(scramble())
    for k in keys: SC[k].append(mm.get(k, np.nan))
rows = []
for k in keys:
    v = np.array(SC[k], float); v = v[np.isfinite(v)]; o = obs[k]; med = np.median(v)
    sc = 1e3 if k[2] in TM else 1.0
    rows.append(dict(phase=k[0], compare=k[1], metric=k[2], obs=o * sc, scr_median=med * sc,
                     scr_rsd=1.4826 * np.median(np.abs(v - med)) * sc, pct_low=float((v <= o).mean()),
                     pct_two=float(min(1, 2 * min((v <= o).mean(), (v >= o).mean()))), n_draw=len(v),
                     n_ets=len(mids), units="ms" if k[2] in TM else ("fraction" if k[2].startswith("A_") else "CC")))
T = pd.DataFrame(rows)
T.to_csv(f"{WS}/tables/leadtest_{ns}_{a.band}_{a.mask}.csv", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_rows", 500)
print(T.round(4).to_string(index=False))
print(f"done {time.time() - t0:.0f} s", flush=True)
