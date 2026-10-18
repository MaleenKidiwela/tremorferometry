#!/usr/bin/env python
"""Fleet 3-5 s measurability census (Merlin 2026-09-02 ruling, Step 2). READ-ONLY on the daily-stack npz.

Per station (data/long_window_daily_<TAG>_Z.npz), per certified family (data/<stem>_causality_cert.csv reliable,
minus data/family_user_blocklist.csv): build the 30-calendar-day trailing stacks EXACTLY as dvv_roll30cal.py does
(cumsum over date-sorted daily stacks, >= 5 stacks in window), then per family the MEDIAN over rows of the
amplitude ratio RMS(window)/RMS(8-10 s) for windows 2-4, 3-5, 4-5 s (stack time; S at +1.0 s), plus the mirror
(-2..0 s) ratio for the record. Also the grand-stack (mean of all certified daily stacks) RMS-by-0.5-s-bin
profile relative to the 8-10 s floor.

Writes NEW files only:  data/window_snr_census.csv  (tag, fam, nrows, snr24, snr35, snr45, snr_mirror)
                        data/window_snr_profile.csv (tag, bin_t0, rms_rel)
"""
import os, sys, glob, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

DAYS, MINST = 30, 5
WINS = {"snr24": (2.0, 4.0), "snr35": (3.0, 5.0), "snr45": (4.0, 5.0), "snr_mirror": (-2.0, 0.0)}
FLOOR = (8.0, 10.0)


def cert_fams(stem):
    f1 = f"data/{stem}_causality_cert.csv"
    if not os.path.exists(f1):
        return set()
    c = pd.read_csv(f1); out = set(c[c.reliable].fam)
    bp = "data/family_user_blocklist.csv"
    if os.path.exists(bp):
        b = pd.read_csv(bp); out -= set(b[b.stem == stem].fam)
    return out


def rms(S, t, w):
    m = (t >= w[0]) & (t < w[1])
    return np.sqrt(np.mean(S[:, m] ** 2, axis=1))


def one(tag):
    t0 = time.time()
    stem = tag.lower().replace("p90f40", "")
    fams = cert_fams(stem)
    fn = f"data/long_window_daily_{tag}_Z.npz"
    if not fams or not os.path.exists(fn):
        return tag, [], [], f"skip (fams {len(fams)}, npz {os.path.exists(fn)})"
    z = np.load(fn, allow_pickle=True)
    t = np.asarray(z["t"], float); patches = np.asarray(z["patches"]).astype(str)
    dates = np.asarray(z["dates"]).astype(str); stacks = z["stacks"]
    sel = np.isin(patches, list(fams))
    idx = np.where(sel)[0]
    if not len(idx):
        return tag, [], [], "no certified rows"
    S_all = np.asarray(stacks[idx], np.float32); P = patches[idx]; D = dates[idx]
    del stacks, z
    # grand-stack profile (mean of L2-normalised daily stacks of certified families)
    G = S_all.astype(np.float64).mean(0)
    fl = np.sqrt(np.mean(G[(t >= FLOOR[0]) & (t < FLOOR[1])] ** 2))
    prof = []
    for b0 in np.arange(-3.0, 10.0, 0.5):
        m = (t >= b0) & (t < b0 + 0.5)
        if m.any():
            prof.append(dict(tag=tag, bin_t0=b0, rms_rel=float(np.sqrt(np.mean(G[m] ** 2)) / max(fl, 1e-12))))
    rows = []
    for fam in np.unique(P):
        k = np.where(P == fam)[0]
        if len(k) < 8:
            continue
        o = np.argsort(D[k]); S = S_all[k[o]].astype(np.float64); Dn = pd.to_datetime(D[k[o]]).values.astype("datetime64[D]")
        n = len(S); ii = np.arange(n)
        cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
        los = np.searchsorted(Dn, Dn - np.timedelta64(DAYS, "D"), side="right")
        counts = ii + 1 - los; keep = counts >= MINST
        if int(keep.sum()) < 8:
            continue
        R = ((cs[ii + 1] - cs[los]) / counts[:, None])[keep]
        base = rms(R, t, FLOOR); base[base == 0] = np.nan
        rec = dict(tag=tag, fam=fam, nrows=int(keep.sum()))
        for kname, w in WINS.items():
            rec[kname] = float(np.nanmedian(rms(R, t, w) / base))
        rows.append(rec)
    return tag, rows, prof, f"{len(rows)} fams, {len(idx):,} rows, {time.time() - t0:.0f}s"


if __name__ == "__main__":
    fams_daily = pd.read_csv("fault_tomography/inversion/res_catalog/families_daily.csv")
    tags = sorted(fams_daily.tag.unique())
    workers = int(os.environ.get("WORKERS", "6"))
    allrows, allprof = [], []
    with ProcessPoolExecutor(workers) as ex:
        for tag, rows, prof, msg in ex.map(one, tags):
            print(f"{tag}: {msg}", flush=True)
            allrows += rows; allprof += prof
    pd.DataFrame(allrows).to_csv("data/window_snr_census.csv", index=False)
    pd.DataFrame(allprof).to_csv("data/window_snr_profile.csv", index=False)
    c = pd.DataFrame(allrows)
    print(f"\nFLEET: {c.tag.nunique()} stations, {len(c)} certified families | snr35 median {c.snr35.median():.2f} "
          f"| frac >=1.41: {(c.snr35 >= 1.41).mean():.2f} | frac >=2: {(c.snr35 >= 2).mean():.2f} "
          f"| snr24 median {c.snr24.median():.2f}")
    print("wrote data/window_snr_census.csv, data/window_snr_profile.csv")
