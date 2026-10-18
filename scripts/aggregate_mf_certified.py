#!/usr/bin/env python
"""One pass over data/mf_*.csv (1.5 TB): keep CERTIFIED families only (minus user blocklist,
station date vetoes applied), aggregate to (a) time x along-strike histograms at 3 cc thresholds
(1-day x 10-km bins) for the Fig-2A-style figure, and (b) a per-(station|family) DAILY COUNT matrix
at cc>=0.8 for duty/rate science. Along-strike axis = the tile-axis convention (G.npz cxy PCA).
Writes data/mf_certified_agg.npz. Parallel over files; each worker streams with usecols.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import glob, re
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

T0 = pd.Timestamp("2009-01-01"); T1 = pd.Timestamp("2026-08-01")
NDAY = (T1 - T0).days
SBINS = np.arange(-650, 651, 10.0)
CCS = [0.0, 0.8, 0.85]

g = np.load("fault_tomography/inversion/res_catalog/G.npz", allow_pickle=True)
cxy = g["cxy"]; lat0, lon0 = float(g["lat0"]), float(g["lon0"])
a = cxy - cxy.mean(0); u = np.linalg.eigh(a.T @ a)[1][:, -1]
if (a @ u)[np.argmax(g["cell_lat"])] < 0: u = -u
OFF = cxy.mean(0)

bl = pd.read_csv("data/family_user_blocklist.csv")
BLOCK = {}
for _, r in bl.iterrows(): BLOCK.setdefault(r.stem, set()).add(r.fam)
DBL = pd.read_csv("data/station_date_blocklist.csv")

def cert(stem):
    p = f"data/{stem}_causality_cert.csv"
    if os.path.exists(p):
        c = pd.read_csv(p); s = set(c[c.reliable].fam)
    else:
        p = f"data/{stem}_fwd_vs_rev_coda.csv"
        if not os.path.exists(p): return set()
        c = pd.read_csv(p); s = set(c[c.ratio > 1.5].fam)
    return s - BLOCK.get(stem, set())

def s_of(lat, lon):
    x = (lon - lon0) * 111 * np.cos(np.radians(lat0))
    y = (lat - lat0) * 111
    return (np.c_[x, y] - OFF) @ u

def work(fp):
    stem = re.sub(r"_(all|combined76|old32_capped|\d{4})\.csv$", "", os.path.basename(fp)[3:])
    stem = stem.replace("p90f40", "")
    fams = cert(stem)
    if not fams: return None
    H = [np.zeros((NDAY, len(SBINS)-1), np.int32) for _ in CCS]
    daily = {}
    try:
        for ch in pd.read_csv(fp, usecols=["template", "time", "cc"], chunksize=2_000_000):
            ch = ch[ch.template.isin(fams)]
            if not len(ch): continue
            t = pd.to_datetime(ch.time.values, errors="coerce", format="mixed")
            ok = np.asarray(pd.notna(t))
            ch = ch[ok]; t = t[ok]
            m = np.asarray((t >= T0) & (t < T1))
            ch = ch[m]; t = t[m]
            if not len(ch): continue
            for _, r in DBL[DBL.stem == stem].iterrows():
                k = np.asarray((t >= pd.Timestamp(r.t0)) & (t <= pd.Timestamp(r.t1)))
                ch = ch[~k]; t = t[~k]
            if not len(ch): continue
            sp = ch.template.str.split("_")
            la = sp.str[0].astype(float).values; lo = sp.str[1].astype(float).values
            s = s_of(la, lo)
            di = np.asarray((t - T0).days)
            si = np.clip(np.digitize(s, SBINS) - 1, 0, len(SBINS)-2)
            cc = ch.cc.values
            for j, th in enumerate(CCS):
                k = cc >= th
                np.add.at(H[j], (di[k], si[k]), 1)
            k = cc >= 0.8
            fam_keys = (stem + "|" + ch.template.values[k])
            for fk, d0 in zip(fam_keys, di[k]):
                arr = daily.setdefault(fk, np.zeros(NDAY, np.int32))
                arr[d0] += 1
    except Exception as e:
        return ("ERR", fp, str(e))
    return (stem, H, daily)

if __name__ == "__main__":
    files = sorted(glob.glob("data/mf_*.csv"))
    Hg = [np.zeros((NDAY, len(SBINS)-1), np.int64) for _ in CCS]
    DAILY = {}
    nerr = 0
    with ProcessPoolExecutor(20) as ex:
        for i, res in enumerate(ex.map(work, files, chunksize=1)):
            if res is None: continue
            if res[0] == "ERR":
                nerr += 1; print("ERR", res[1], res[2][:80], flush=True); continue
            stem, H, daily = res
            for j in range(len(CCS)): Hg[j] += H[j]
            for k, arr in daily.items():
                if k in DAILY: DAILY[k] += arr
                else: DAILY[k] = arr
            if i % 100 == 0: print(f"{i}/{len(files)} files", flush=True)
    keys = sorted(DAILY)
    M = np.stack([np.minimum(DAILY[k], 32767).astype(np.int16) for k in keys]) if keys else np.zeros((0, NDAY), np.int16)
    np.savez_compressed("data/mf_certified_agg.npz",
                        H_all=Hg[0], H_cc80=Hg[1], H_cc85=Hg[2], sbins=SBINS,
                        t0=str(T0.date()), nday=NDAY, fam_keys=np.array(keys, object),
                        daily_cc80=M, ccs=CCS)
    print(f"DONE: {len(keys):,} certified (station|family) series | errors {nerr} | "
          f"total detections cc>=0.8: {Hg[1].sum():,}")
