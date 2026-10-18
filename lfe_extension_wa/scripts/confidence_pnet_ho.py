#!/usr/bin/env python
"""CONFIDENCE (2026-09-30) helper: PNET features read ONLY at the 6 stations never used for PNET fine-tuning
(pnet_common.HELD_OUT = UW.HDW, UW.OSD, PB.B018, PB.B017, US.NLWA, UW.GNW).  Same rule as pnet_feat.py
(max CF within +-0.5 s of the template's S time; scrambled control read at the site's rolled position), but the
family-day station list is restricted to the held-out stations.  Columns get the suffix _ho.

  WA_SETUP=results/merlin_tests/t1_setup python scripts/confidence_pnet_ho.py --cf ft1 --workers 12
  -> results/confidence/pnet_ho_<cf>/feat_<date>.parquet  (fid, t, series, pnet_*_ho)
"""
from __future__ import annotations
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC
import pnet_feat as PF

HO = set(PC.HELD_OUT)
COLS = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_rel_max"]


def day_job(args):
    ds, tag, od = args
    t0 = time.time()
    TOFF, SITE = PF._G["T"], PF._G["SITE"]
    NV = W.NDAY - W.NT + 1
    D0 = pd.Timestamp(ds).timestamp()
    F = pd.read_parquet(f"{PF.FEAT0}/feat_{ds}.parquet", columns=["fid", "t", "series"])
    J1 = json.load(open(f"{W.RES}/mf_{PF.MF1}/run_{ds}.json"))
    SH = {"data": {}, "c1": J1["shift_site"], "c2": J1["shift_site2"]}
    p2 = f"{W.RES}/mf_{PF.MF2}/run_{ds}.json"
    if os.path.exists(p2):
        J2 = json.load(open(p2)); SH["c3"] = J2["shift_site"]; SH["c4"] = J2["shift_site2"]
    nodes = J1["nodes"]
    FD = pd.read_parquet(f"{W.RES}/mf_{PF.MF1}/famday_{ds}.parquet")
    FD = FD[FD.observable]
    fam_sta = {int(f): [x for x in s.split(",") if x in HO] for f, s in zip(FD.fid, FD.stations)}
    MFc = {}
    for ns in sorted({s for v in fam_sta.values() for s in v}):
        if nodes.get(ns) is None:
            continue
        cf = PC.load_cf(tag, ns, ds)
        if cf is not None:
            MFc[ns] = PC.maxfilt(cf, 10)
    rng = np.random.default_rng(1000 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    R = rng.integers(0, NV, PF.NRAND)
    for c in COLS:
        F[c + "_ho"] = np.nan
    ii = np.round((F.t.to_numpy() - D0) * W.FS).astype(np.int64)
    for f, gi in F.groupby("fid").groups.items():
        sts = [ns for ns in fam_sta.get(int(f), []) if ns in MFc and (int(f), ns, nodes[ns]) in TOFF]
        if not sts:
            continue
        g = F.loc[gi]
        sites = sorted({SITE[ns] for ns in sts})
        qi = {q: k for k, q in enumerate(sites)}
        V = np.full((len(g), len(sites)), -1.0); VR = np.full((PF.NRAND, len(sites)), -1.0)
        i = ii[np.asarray(gi)]
        for ns in sts:
            d = int(round(TOFF[(int(f), ns, nodes[ns])] * W.FS))
            q = SITE[ns]
            sh = np.array([int(SH[s].get(q, 0)) % NV if s != "data" else 0 for s in g.series], np.int64)
            j = ((i - sh) % NV) + d
            jr = R + d
            cf = MFc[ns]
            v = np.where((j >= 0) & (j < W.NDAY), cf[np.clip(j, 0, W.NDAY - 1)], np.nan)
            vr = np.where((jr >= 0) & (jr < W.NDAY), cf[np.clip(jr, 0, W.NDAY - 1)], np.nan)
            k = qi[q]
            V[:, k] = np.fmax(V[:, k], np.nan_to_num(v, nan=-1.0))
            VR[:, k] = np.fmax(VR[:, k], np.nan_to_num(vr, nan=-1.0))
        V[V < 0] = np.nan; VR[VR < 0] = np.nan
        Vs = np.sort(np.nan_to_num(V, nan=-1.0), 1)[:, ::-1]
        med = np.nanmedian(VR, 0) if np.isfinite(VR).any() else np.full(len(sites), np.nan)
        with np.errstate(invalid="ignore"):
            rel = V - med[None]
        F.loc[gi, "pnet_max_ho"] = np.where(Vs[:, 0] >= 0, Vs[:, 0], np.nan)
        F.loc[gi, "pnet_2nd_ho"] = np.where(Vs[:, 1] >= 0, Vs[:, 1], np.nan) if len(sites) > 1 else np.nan
        F.loc[gi, "n_sites_pnet03_ho"] = (np.nan_to_num(V, nan=0) > 0.3).sum(1)
        F.loc[gi, "n_sites_pnet05_ho"] = (np.nan_to_num(V, nan=0) > 0.5).sum(1)
        F.loc[gi, "pnet_rel_max_ho"] = np.nanmax(np.where(np.isfinite(rel), rel, -9), 1)
    F.loc[F.pnet_rel_max_ho <= -8, "pnet_rel_max_ho"] = np.nan
    os.makedirs(od, exist_ok=True)
    F.to_parquet(f"{od}/feat_{ds}.parquet.tmp", index=False)
    os.replace(f"{od}/feat_{ds}.parquet.tmp", f"{od}/feat_{ds}.parquet")
    return ds, len(F), len(MFc), float(F.pnet_max_ho.notna().mean()), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True)
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    od = f"{W.RES}/confidence/pnet_ho_{a.cf}"
    days = [d for w in ["V2008", "V2010", "A2011"] for d, _ in W.window_days(w)]
    days = [d for d in days if os.path.exists(f"{PF.FEAT0}/feat_{d}.parquet") and not os.path.exists(f"{od}/feat_{d}.parquet")]
    print(f"[pnet_ho:{a.cf}] {len(days)} days", flush=True)
    with ProcessPoolExecutor(a.workers, initializer=PF._init) as ex:
        fut = {ex.submit(day_job, (d, a.cf, od)): d for d in days}
        for fu in as_completed(fut):
            ds, nr, ncf, fr, secs = fu.result()
            print(f"[pnet_ho:{a.cf}] {ds}: {nr} cand, {ncf} held-out CFs, defined {fr:.3f}, {secs:.0f}s", flush=True)


if __name__ == "__main__":
    main()
