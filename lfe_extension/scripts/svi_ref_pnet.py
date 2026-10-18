#!/usr/bin/env python
"""SVI_REF PNET features (GPU) -- lfe_picker's SVI-trained PhaseNet-style S picker (results/p1/models/pnet_w05.pt; trained on
2008 + 2009 + 2010 + 2012 at S15 minus B926 / NLLB / LRIV, P1_RESULTS.md l.41) read at every candidate's expected S time.

Input contract (lfe_picker/scripts/p_scan.cf_day): 60-s windows, 20 Hz, the 1-8 Hz Butterworth-4 zero-phase chain, Z / N / E,
per-window peak normalisation -- the SVI 20 Hz cache, with PB borehole EH1 / EH2 ROTATED to N / E from the local StationXML
azimuths (lfe_picker/results/p4/stationxml/PB.xml; B926 from results/a4_stationxml_audit.csv).  The frozen detector chain does
not rotate; only this feature does.  Expected S at a station = the candidate's origin time + the frozen placement + 0.5 s; a
scrambled-control candidate is read at its site's ROLLED position (the MF's per-site shift).  Station value = max CF within
+-0.5 s; site value = max over the site's stations.  Features: pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05, pnet_mean,
pnet_rel_max (site value minus that site's median at 300 random times of the day), per station pn_<STA>.
The gate (Merlin) uses 2011 (fit) and 2013 (one look) only; production features only if the gate passes.

  python scripts/svi_ref_pnet.py --runs t2011_S15_nc,t2011_S14_nc,... --shard 0 --nshard 2
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings, traceback
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

PICKER = f"{R.PICK}/scripts"
PNET_W05 = f"{R.PICK}/results/p1/models/pnet_w05.pt"
FS = R.FS
OP = f"{R.OUT}/pnet"


def az_table():
    rows = []
    from obspy import read_inventory
    inv = read_inventory(f"{R.PICK}/results/p4/stationxml/PB.xml")
    for net in inv:
        for st in net:
            for ch in st:
                if ch.code[:2] == "EH" and ch.code[-1] in "12":
                    rows.append(dict(sta=st.code, cha=ch.code, start=str(ch.start_date)[:10],
                                     end=str(ch.end_date)[:10] if ch.end_date else "2100-01-01", az=float(ch.azimuth)))
    A = pd.read_csv(f"{R.RES}/a4_stationxml_audit.csv")
    for r in A[A.sta == "B926"].itertuples():
        if r.chan[-1] in "12":
            rows.append(dict(sta="B926", cha=r.chan, start=str(r.start)[:10], end="2100-01-01" if r.end == "open" else str(r.end)[:10], az=float(r.az)))
    return pd.DataFrame(rows)


def azimuth(AZ, sta, cha, date):
    g = AZ[(AZ.sta == sta) & (AZ.cha == cha) & (AZ.start <= date) & (AZ.end > date)]
    return float(g.az.iloc[0]) if len(g) else np.nan


def load_zne(ns, sta, date, caches, AZ):
    d = pd.Timestamp(date)
    p = next((f"{c}/{ns}/{d.year}/{d.dayofyear:03d}.npz" for c in caches if os.path.exists(f"{c}/{ns}/{d.year}/{d.dayofyear:03d}.npz")), None)
    if p is None:
        return None, "nodata"
    try:
        z = np.load(p); x = z["x"].astype(np.float32); t0 = float(z["t0"]); order = [str(c) for c in z["order"]]
    except Exception:
        return None, "error"
    if x.shape[0] != 3:
        return None, "short"
    if order[1][-1] == "1" and order[2][-1] == "2":
        a1, a2 = azimuth(AZ, sta, order[1], date), azimuth(AZ, sta, order[2], date)
        if not (np.isfinite(a1) and np.isfinite(a2)):
            return None, "error"
        r1, r2 = np.radians(a1), np.radians(a2)
        Ai = np.linalg.inv(np.array([[np.cos(r1), np.sin(r1)], [np.cos(r2), np.sin(r2)]])).astype(np.float32)
        x[1:3] = Ai @ x[1:3]
    return dict(x=x, t0=t0), "ok"


def maxfilt(v, half=10):
    from scipy.ndimage import maximum_filter1d
    return maximum_filter1d(v, size=2 * half + 1, mode="constant", cval=0.0)


def placement():
    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    return C.placement(tpl, sp, F, mv, C.station_table(), L)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--runs", required=True)
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--tag", default="gate")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    import torch
    sys.path.insert(0, PICKER); import p_scan
    model, mt = p_scan.load_model(PNET_W05)
    AZ = az_table(); OFF = placement(); nsof = R.netsta_of()
    SITE = C.site_map(R.S15)
    runs = [r for r in a.runs.split(",") if r]
    jobs = {}
    for run in runs:
        for p in sorted(glob.glob(f"{R.OUT}/feat/{run}/*.json")):
            ds = os.path.basename(p)[:10]
            jobs.setdefault(ds, []).append(run)
    dates = [d for j, d in enumerate(sorted(jobs)) if j % a.nshard == a.shard]
    R.ensure(f"{R.OUT}/status")
    status = []; t00 = time.time()
    for n, ds in enumerate(dates):
        todo = [r for r in jobs[ds] if not os.path.exists(f"{OP}/feat/{r}/{ds}.parquet")]
        if not todo: continue
        try:
            metas = {r: json.load(open(f"{R.OUT}/feat/{r}/{ds}.json")) for r in todo}
            y = int(ds[:4]); caches = R.TRUTH[y]["caches"] if y in R.TRUTH else [R.CACHE_R1, C.CACHE]
            stas = sorted({s for m in metas.values() for s in m["stations"]})
            CF = {}
            for s in stas:
                z, st = load_zne(nsof[s], s, ds, caches, AZ)
                status.append((ds, nsof[s], st, ""))
                if z is None: continue
                CF[s] = (maxfilt(p_scan.cf_day(model, mt, z["x"]), 10), z["t0"])
            rng = np.random.default_rng(int(ds.replace("-", "")))
            day0 = pd.Timestamp(ds).timestamp()
            for run in todo:
                meta = metas[run]
                fp = f"{R.OUT}/feat/{run}/{ds}.parquet"
                F = pd.read_parquet(fp, columns=["series", "k", "fam", "t", "i"]) if os.path.getsize(fp) > 0 else pd.DataFrame()
                out = F[["series", "k"]].copy() if len(F) else pd.DataFrame(columns=["series", "k"])
                if len(F):
                    FD = pd.read_parquet(f"{R.OUT}/feat/{run}/{ds}_fd.parquet")
                    FD = FD[FD.observable.astype(bool)].set_index(["fam", "series"]) if len(FD) and "tref" in FD else None
                    SH = {q: int(v) for q, v in meta.get("shifts", {}).items()}
                    sts = [s for s in meta["stations"] if s in CF]
                    cols = {f"pn_{s}": np.full(len(F), np.nan, np.float32) for s in sts}
                    for f, gi in F.groupby("fam").indices.items():
                        g = F.iloc[gi]
                        for s in sts:
                            if (float(f), s) not in OFF: continue
                            o = OFF[(float(f), s)] + C.PRE
                            cf, t0 = CF[s]
                            tS = g.t.to_numpy(float) + o
                            c1 = (g.series == "c1").to_numpy()
                            if c1.any() and FD is not None and (f, "c1") in FD.index:
                                tref = float(FD.loc[(f, "c1"), "tref"]); nvo = int(FD.loc[(f, "c1"), "NVo"])
                                iq = (g.i.to_numpy()[c1] + SH.get(SITE[s], 0) % nvo) % nvo
                                tS[c1] = tref + iq / FS + o
                            j = np.rint((tS - t0) * FS).astype(np.int64)
                            ok = (j >= 0) & (j < len(cf))
                            v = np.full(len(g), np.nan); v[ok] = cf[j[ok]]
                            cols[f"pn_{s}"][gi] = v
                    sites = sorted({SITE[s] for s in sts})
                    V = np.full((len(F), len(sites)), np.nan); MED = np.full(len(sites), np.nan)
                    for jq, q in enumerate(sites):
                        mem = [s for s in sts if SITE[s] == q]
                        V[:, jq] = np.nanmax(np.vstack([cols[f"pn_{s}"] for s in mem]), 0) if mem else np.nan
                        rt = rng.uniform(600, 86400 - 600, 300)
                        mv = [np.median(CF[s][0][np.clip(np.rint((day0 + rt - CF[s][1]) * FS).astype(int), 0, len(CF[s][0]) - 1)]) for s in mem]
                        MED[jq] = float(np.max(mv)) if mv else np.nan
                    with np.errstate(invalid="ignore", all="ignore"):
                        Vs = np.sort(np.nan_to_num(V, nan=-1.0), 1)[:, ::-1]
                        out["pnet_max"] = np.where(Vs[:, 0] >= 0, Vs[:, 0], np.nan)
                        out["pnet_2nd"] = np.where(Vs[:, 1] >= 0, Vs[:, 1], np.nan) if V.shape[1] > 1 else np.nan
                        out["n_sites_pnet03"] = (np.nan_to_num(V, nan=0) > 0.3).sum(1).astype(float)
                        out["n_sites_pnet05"] = (np.nan_to_num(V, nan=0) > 0.5).sum(1).astype(float)
                        out["pnet_mean"] = np.nanmean(V, 1)
                        out["pnet_rel_max"] = np.nanmax(V - MED[None], 1)
                    for c_, v_ in cols.items(): out[c_] = v_
                od = f"{OP}/feat/{run}"; R.ensure(od)
                out.to_parquet(f"{od}/{ds}.parquet.tmp", index=False); os.replace(f"{od}/{ds}.parquet.tmp", f"{od}/{ds}.parquet")
            del CF; torch.cuda.empty_cache()
            if n < 3 or n % 20 == 0:
                R.log(f"  [pnet {a.shard}] {ds}: runs {len(todo)}; {n + 1}/{len(dates)} {time.time() - t00:.0f}s")
        except Exception as e:
            R.log(f"  ERROR {ds}: {e}\n{traceback.format_exc()[-1200:]}")
            status += [(ds, "ALL", "error", f"day job: {str(e)[:100]}")]
        if n % 20 == 0: R.write_status(f"pnet_{a.tag}", f"s{a.shard}", status)
    R.write_status(f"pnet_{a.tag}", f"s{a.shard}", status)
    R.log(f"[pnet {a.shard}] done {len(dates)} dates {time.time() - t00:.0f}s")


if __name__ == "__main__":
    main()
