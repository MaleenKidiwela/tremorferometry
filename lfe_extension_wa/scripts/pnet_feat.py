#!/usr/bin/env python
"""Stage A/C: PNET features for every STEP-1 candidate (data and the 4 site-rolled scrambled controls).

For candidate origin t of family f on day d, at every station ns of the family-day network (the stations the
all-templates MF used: famday `stations`, node = the day's node, a >= 30-event build template exists):
  expected S sample j = ((i - shift_site(series)) mod NV) + round(offset_s * FS)
     i = candidate origin sample; offset_s = the family's template start at that station (= S envelope peak - 0.5 s,
     i.e. the onset; PNET's pick precedes the envelope peak, SVI time-base check); a scrambled-control candidate is
     read at its site's ROLLED position, exactly as the MF built that control (step1_features.py rule).
  station value = max of the PNET S probability over [j - 0.5 s, j + 0.5 s]
  site value    = max over the site's stations
Features: pnet_max (largest site value), pnet_2nd (2nd largest), n_sites_pnet03 / n_sites_pnet05 (sites > 0.3 / 0.5),
pnet_mean (mean over sites), pnet_nsites (sites with a CF), pnet_rel_max = max over sites of (site value - that
site's median over 1000 random valid origins of the same family-day)  [a noise-level-relative version, a column].

  python scripts/pnet_feat.py --cf w05 --workers 12
  -> results/pnet/feat_<cf>/feat_<date>.parquet = results/step1_prob/feat/feat_<date>.parquet + PNET columns
"""
from __future__ import annotations
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

FEAT0 = f"{W.RES}/step1_prob/feat"
MF1, MF2 = "mt1_all", "step1c34"
NRAND = 1000
PN = ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_mean", "pnet_nsites", "pnet_rel_max"]
_G = {}


def _init():
    import wa_mf
    T, _ = wa_mf.load_tpl(["build"], all_templates=True)
    S = W.load_stations()
    _G.update(T={k: v["off"] for k, v in T.items()}, SITE=dict(zip(S.netsta, S.site)))


def day_job(args):
    ds, tag = args
    t0 = time.time()
    TOFF, SITE = _G["T"], _G["SITE"]
    n = W.NT
    NV = W.NDAY - n + 1
    D0 = pd.Timestamp(ds).timestamp()
    F = pd.read_parquet(f"{FEAT0}/feat_{ds}.parquet")
    J1 = json.load(open(f"{W.RES}/mf_{MF1}/run_{ds}.json"))
    SH = {"data": {}, "c1": J1["shift_site"], "c2": J1["shift_site2"]}
    p2 = f"{W.RES}/mf_{MF2}/run_{ds}.json"
    if os.path.exists(p2):
        J2 = json.load(open(p2)); SH["c3"] = J2["shift_site"]; SH["c4"] = J2["shift_site2"]
    nodes = J1["nodes"]
    FD = pd.read_parquet(f"{W.RES}/mf_{MF1}/famday_{ds}.parquet")
    FD = FD[FD.observable]
    fam_sta = {int(f): s.split(",") for f, s in zip(FD.fid, FD.stations)}
    MFc = {}
    for ns in sorted({s for v in fam_sta.values() for s in v}):
        if nodes.get(ns) is None:
            continue
        cf = PC.load_cf(tag, ns, ds)
        if cf is not None:
            MFc[ns] = PC.maxfilt(cf, 10)
    rng = np.random.default_rng(1000 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    R = rng.integers(0, NV, NRAND)
    for c in PN:
        F[c] = np.nan
    ii = np.round((F.t.to_numpy() - D0) * W.FS).astype(np.int64)
    for f, gi in F.groupby("fid").groups.items():
        sts = [ns for ns in fam_sta.get(int(f), []) if ns in MFc and (int(f), ns, nodes[ns]) in TOFF]
        if not sts:
            continue
        g = F.loc[gi]
        sites = sorted({SITE[ns] for ns in sts})
        qi = {q: k for k, q in enumerate(sites)}
        V = np.full((len(g), len(sites)), -1.0); VR = np.full((NRAND, len(sites)), -1.0)
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
        F.loc[gi, "pnet_max"] = np.where(Vs[:, 0] >= 0, Vs[:, 0], np.nan)
        F.loc[gi, "pnet_2nd"] = np.where(Vs[:, 1] >= 0, Vs[:, 1], np.nan) if len(sites) > 1 else np.nan
        F.loc[gi, "n_sites_pnet03"] = (np.nan_to_num(V, nan=0) > 0.3).sum(1)
        F.loc[gi, "n_sites_pnet05"] = (np.nan_to_num(V, nan=0) > 0.5).sum(1)
        F.loc[gi, "pnet_mean"] = np.nanmean(V, 1) if np.isfinite(V).any() else np.nan
        F.loc[gi, "pnet_nsites"] = np.isfinite(V).sum(1)
        F.loc[gi, "pnet_rel_max"] = np.nanmax(np.where(np.isfinite(rel), rel, -9), 1)
    F.loc[F.pnet_rel_max <= -8, "pnet_rel_max"] = np.nan
    od = f"{PC.OUT}/feat_{tag}"
    os.makedirs(od, exist_ok=True)
    F.to_parquet(f"{od}/feat_{ds}.parquet.tmp", index=False)
    os.replace(f"{od}/feat_{ds}.parquet.tmp", f"{od}/feat_{ds}.parquet")
    return ds, len(F), len(MFc), float(F.pnet_max.notna().mean()), time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--days", default="")
    a = ap.parse_args()
    days = [d for w in ["V2008", "V2010", "A2011"] for d, _ in W.window_days(w)]
    if a.days:
        days = a.days.split(",")
    od = f"{PC.OUT}/feat_{a.cf}"
    days = [d for d in days if os.path.exists(f"{FEAT0}/feat_{d}.parquet") and not os.path.exists(f"{od}/feat_{d}.parquet")]
    W.log(PC.LOG, f"[feat:{a.cf}] {len(days)} days")
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = {ex.submit(day_job, (d, a.cf)): d for d in days}
        for fu in as_completed(fut):
            ds, nr, ncf, fr, secs = fu.result()
            W.log(PC.LOG, f"[feat:{a.cf}] {ds}: {nr} candidates, {ncf} station CFs, pnet_max defined {fr:.3f}, {secs:.0f}s")


if __name__ == "__main__":
    main()
