#!/usr/bin/env python
"""Merlin-prescribed shape-program kill-pilots (2026-09-01 ruling), run on the stored daily stacks.

Metrics per (station, family), on FULL / ODD-day / EVEN-day grand stacks (plain mean of L2-normalized
daily stacks; NO Wiener; horizontal energy H1^2+H2^2 primary, Z fallback):
  A1  early-coda decay rate: OLS slope of ln[E(t) t^2] over t = +4..+9 s (S at +1; starts 3 s after S,
      ends 1 s before edge), fit truncated where E < 3x mirror floor (E over t in [-3,-1]).
  A3  coda/direct energy ratio: E[+4..+9] / E[0..+2].
  A2r spectral centroid 2-8 Hz of the direct window (for co-located source-difference pilot only).
PILOT 1 reliability: corr(odd, even) across pairs per metric (kill < 0.5).
PILOT 2 variance decomposition: metric ~ station FE + distance + source depth; report R2 shares and
        whether residual spread exceeds the split-half noise.
PILOT 3 co-located pairs (same station, same 0.1-deg cell): |metric difference| vs split-half floor.
Sample: every multi-station family's (station,family) pairs + all their stations' co-located sets,
capped per station. Writes data/shape_pilot_metrics.csv + prints verdicts.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import glob, re
import numpy as np, pandas as pd

bl = pd.read_csv("data/family_user_blocklist.csv")
BLOCK = {}
for _, r in bl.iterrows(): BLOCK.setdefault(r.stem, set()).add(r.fam)
def cert(stem):
    p = f"data/{stem}_causality_cert.csv"
    if os.path.exists(p):
        c = pd.read_csv(p); s = set(c[c.reliable].fam)
    else:
        p = f"data/{stem}_fwd_vs_rev_coda.csv"
        if not os.path.exists(p): return set()
        c = pd.read_csv(p); s = set(c[c.ratio > 1.5].fam)
    return s - BLOCK.get(stem, set())

# station coords
COORD = {r.sta: (r.lat, r.lon) for _, r in pd.read_csv("data/candidate_stations_post2020.csv").iterrows()}
COORD.update({"PGC": (48.6498, -123.4521), "SHB": (49.5990, -123.8770), "CLRS": (48.8203, -124.1309)})
g = np.load("fault_tomography/inversion/res_catalog/G.npz", allow_pickle=True)
CELL_DEP = {(round(float(la), 2), round(float(lo), 2)): float(dp)
            for la, lo, dp in zip(g["cell_lat"], g["cell_lon"], g["depth_km"])}

# ---- choose the sample: multi-station families + co-located sets ----
zfiles = sorted(glob.glob("data/long_window_daily_*_Z.npz"))
sta_of = {}
for f in zfiles:
    tag = os.path.basename(f).split("long_window_daily_")[1][:-6]
    sta_of[tag] = f
fam_stations = {}
for tag in sta_of:
    stem = tag.lower().replace("p90f40", "")
    for fam in cert(stem):
        fam_stations.setdefault(fam, []).append(tag)
multi = {f: ts for f, ts in fam_stations.items() if len(ts) >= 2}
FULL = os.environ.get("FULL", "0") == "1"
target_st = sorted(sta_of) if FULL else sorted({t for ts in multi.values() for t in ts})
print(f"multi-station families: {len(multi)} across {len(target_st)} stations", flush=True)

def envelope2(x, fs, smooth_s=0.25):
    e = x.astype(np.float64) ** 2
    w = max(3, int(smooth_s * fs))
    k = np.ones(w) / w
    return np.convolve(e, k, mode="same")

def metrics(stack_mean, t, fs):
    E = envelope2(stack_mean, fs)
    mir = E[(t >= -3) & (t <= -1)].mean()
    dfit = (t >= 4) & (t <= 9) & (E > 3 * mir)
    out = {}
    if dfit.sum() >= 12:
        y = np.log(E[dfit] * t[dfit] ** 2)
        out["decay"] = float(np.polyfit(t[dfit], y, 1)[0])
    else:
        out["decay"] = np.nan
    Ec = E[(t >= 4) & (t <= 9)].mean(); Ed = E[(t >= 0) & (t <= 2)].mean()
    out["cdr"] = float(np.log10(Ec / Ed)) if Ed > 0 else np.nan
    # selection-immune variant (Merlin step-2): BOTH windows post-template, gap between
    E1 = E[(t >= 2.5) & (t <= 5)].mean(); E2 = E[(t >= 6.5) & (t <= 9)].mean()
    out["cor"] = float(np.log10(E2 / E1)) if E1 > 0 else np.nan
    seg = stack_mean[(t >= 0) & (t <= 2)]
    sp = np.abs(np.fft.rfft(seg * np.hanning(len(seg)))) ** 2
    fr = np.fft.rfftfreq(len(seg), 1 / fs)
    b = (fr >= 2) & (fr <= 8)
    out["cent"] = float((fr[b] * sp[b]).sum() / sp[b].sum()) if sp[b].sum() > 0 else np.nan
    return out

rows = []
MAXFAM_PER_ST = 10**9 if os.environ.get("FULL","0")=="1" else 40
for n_st, tag in enumerate(target_st):
    stem = tag.lower().replace("p90f40", "")
    fams = cert(stem)
    want = set(fams) if os.environ.get("FULL","0")=="1" else {f for f in fams if f in multi}
    # add co-located extras: any cell with >=2 certified fams
    cells = {}
    for f in fams:
        la, lo = f.split("_")[0], f.split("_")[1]
        cells.setdefault((la, lo), []).append(f)
    for c, fl in cells.items():
        if len(fl) >= 2: want.update(fl[:3])
    want = list(want)[:MAXFAM_PER_ST]
    if not want: continue
    try:
        comp = {}
        for cch in ("Z", "H1", "H2"):
            p = f"data/long_window_daily_{tag}_{cch}.npz"
            if os.path.exists(p):
                z = np.load(p, allow_pickle=True)
                comp[cch] = z
        z0 = comp.get("Z") or list(comp.values())[0]
        t = z0["t"]; fs = float(z0["fs"])
        pat = z0["patches"].astype(str)
        stacks = {c: comp[c]["stacks"] for c in comp}
        pats = {c: comp[c]["patches"].astype(str) for c in comp}
    except Exception as e:
        print("skip", tag, str(e)[:60], flush=True); continue
    sla, slo = COORD.get(stem.upper(), COORD.get(tag.upper(), (np.nan, np.nan)))
    for fam in want:
        try:
            idx = {c: np.where(pats[c] == fam)[0] for c in stacks}
            n = len(idx.get("Z", idx[list(idx)[0]]))
            if n < 200: continue
            def grand(sel):
                acc = None
                for c in ("H1", "H2"):
                    if c in stacks and len(idx[c]) == n:
                        m = stacks[c][idx[c][sel]].astype(np.float64).mean(0)
                        acc = m ** 2 if acc is None else acc + m ** 2
                if acc is not None:
                    return np.sqrt(acc)
                return stacks["Z"][idx["Z"][sel]].astype(np.float64).mean(0)
            full = metrics(grand(slice(None)), t, fs)
            odd = metrics(grand(slice(0, None, 2)), t, fs)
            even = metrics(grand(slice(1, None, 2)), t, fs)
            la, lo = float(fam.split("_")[0]), float(fam.split("_")[1])
            dep = CELL_DEP.get((round(round(la / 0.1) * 0.1, 2), round(round(lo / 0.1) * 0.1, 2)), np.nan)
            dist = np.hypot((la - sla) * 111, (lo - slo) * 111 * np.cos(np.radians(la))) if np.isfinite(sla) else np.nan
            rows.append(dict(tag=tag, fam=fam, n=n, lat=la, lon=lo, dep=dep, dist=dist,
                             **{k: full[k] for k in full},
                             **{k + "_o": odd[k] for k in odd}, **{k + "_e": even[k] for k in even}))
        except Exception:
            continue
    if n_st % 20 == 0: print(f"station {n_st}/{len(target_st)}", flush=True)

D = pd.DataFrame(rows)
D.to_csv("data/shape_metrics_all.csv" if os.environ.get("FULL","0")=="1" else "data/shape_pilot_metrics.csv", index=False)
print(f"\n{len(D)} (station,family) pairs measured", flush=True)
print("\nPILOT 1 — split-half reliability (kill < 0.5):")
for m in ("decay", "cdr", "cent"):
    ok = D[[m + "_o", m + "_e"]].dropna()
    r = ok[m + "_o"].corr(ok[m + "_e"])
    print(f"  {m:6s} n={len(ok):4d}  odd-even corr = {r:+.2f}  -> {'PASS' if r > 0.5 else 'KILL'}")
print("\nPILOT 2 — variance decomposition (per metric): R2 station-FE only / +geometry; residual vs split-half noise:")
for m in ("decay", "cdr"):
    ok = D.dropna(subset=[m, "dist", "dep"])
    if len(ok) < 100: print(f"  {m}: too few"); continue
    Xs = pd.get_dummies(ok.tag)
    import numpy.linalg as la_
    def r2(X, y):
        X = np.column_stack([np.ones(len(y)), X])
        b, *_ = la_.lstsq(X, y, rcond=None)
        res = y - X @ b
        return 1 - res.var() / y.var(), res
    y = ok[m].values.astype(float)
    r2s, _ = r2(Xs.values.astype(float), y)
    r2g, res = r2(np.column_stack([Xs.values.astype(float), ok.dist.values, ok.dep.values]), y)
    noise = np.nanstd(ok[m + "_o"] - ok[m + "_e"]) / np.sqrt(2)
    print(f"  {m:6s} R2(stationFE) = {r2s:.2f}  R2(+geometry) = {r2g:.2f}  residual sd = {res.std():.3g} "
          f"vs split-half noise = {noise:.3g}  -> {'medium signal ABOVE noise' if res.std() > 1.5*noise else 'no mappable signal'}")
print("\nPILOT 3 — co-located source differences (same station, same cell):")
diffs, floors = [], []
for (tag, la, lo), grp in D.groupby(["tag", "lat", "lon"]):
    if len(grp) >= 2:
        c = grp.cent.dropna().values
        if len(c) >= 2:
            diffs.append(np.abs(c[0] - c[1]))
            f0 = np.abs(grp.cent_o - grp.cent_e).mean() / np.sqrt(2)
            floors.append(f0)
if diffs:
    print(f"  co-located pairs n={len(diffs)}: median |centroid diff| = {np.median(diffs):.3f} Hz "
          f"vs split-half floor {np.nanmedian(floors):.3f} Hz -> "
          f"{'source signal present' if np.median(diffs) > 1.5*np.nanmedian(floors) else 'below floor — defer A2'}")
else:
    print("  no co-located pairs measured")
print("\nDONE")
