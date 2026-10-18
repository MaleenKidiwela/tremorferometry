#!/usr/bin/env python
"""WA PRODUCTION W8: PNET block for the small gate (truth days only; production only if the gate passes).

Model: the FROZEN WA fine-tune results/pnet/models/ft1.pt (no new training).  CF = lfe_picker p_scan.cf_day on the WA cache
rotated to ZNE (pnet_common.to_zne), computed on the fly per station-day (not stored).  For every candidate of a run day:
expected S at a station = the candidate's aligned position (scrambled control: the site's ROLLED position) + the template
offset; station value = max CF within +-0.5 s; site value = max over its stations (the network's own sites that day).
Block: pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05, pnet_mean, pnet_rel_max (site value minus that site's median at
1000 random origins), and the same six from the stations the fine-tune never trained on (B017, B018, GNW, HDW) only, '_ho'.

  python scripts/wa_prod_pnet.py --run truth [--days-file ...] [--shard 0 --nshard 1]
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.patch_wa()
import pnet_common as PC
sys.path.insert(0, PC.PICKER)
from wa_prod_mf import load_library

HO = {"PB.B017", "PB.B018", "UW.GNW", "UW.HDW"}
MODEL = f"{P.RES}/pnet/models/ft1.pt"
NV = W.NDAY - W.NT + 1


def block(D, sitev, ok_sites, rmed, prefix=""):
    """sitev: (rows, sites) station-max per site; rmed: (sites,) random-origin median."""
    S = np.where(ok_sites[None, :], sitev, np.nan)
    srt = -np.sort(-np.nan_to_num(S, nan=-1), 1)
    nval = np.isfinite(S).sum(1)
    out = {f"pnet_max{prefix}": np.where(nval > 0, srt[:, 0], np.nan),
           f"pnet_2nd{prefix}": np.where(nval > 1, srt[:, 1] if S.shape[1] > 1 else np.nan, np.nan),
           f"n_sites_pnet03{prefix}": (np.nan_to_num(S, nan=0) > 0.3).sum(1).astype(float),
           f"n_sites_pnet05{prefix}": (np.nan_to_num(S, nan=0) > 0.5).sum(1).astype(float),
           f"pnet_mean{prefix}": np.where(nval > 0, np.nansum(S, 1) / np.maximum(nval, 1), np.nan),
           f"pnet_rel_max{prefix}": np.nanmax(np.where(ok_sites[None, :], S - rmed[None, :], np.nan), 1) if ok_sites.any() else np.full(len(S), np.nan)}
    if not ok_sites.any():
        for k in out:
            out[k] = np.full(len(S), np.nan)
    return out


def main():
    import torch, p_scan
    from scipy.ndimage import maximum_filter1d
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--days-file", default="")
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    model, meta = p_scan.load_model(MODEL)
    TPL = load_library(rich=(a.run == "rich"))
    od = f"{P.OUT}/pnet/{a.run}"; P.ensure(od)
    dets = sorted(glob.glob(f"{P.OUT}/mf/{a.run}/det_*.parquet"))
    days = [os.path.basename(p)[4:14] for p in dets]
    if a.days_file:
        keep = set(json.load(open(a.days_file))); days = [d for d in days if d in keep]
    days = [d for j, d in enumerate(days) if j % a.nshard == a.shard]
    status = []
    t0 = time.time()
    for nd_, ds in enumerate(days):
        if os.path.exists(f"{od}/pnet_{ds}.parquet"):
            continue
        D = pd.read_parquet(f"{P.OUT}/mf/{a.run}/det_{ds}.parquet", columns=["k", "net", "series", "fid", "i"])
        J = json.load(open(f"{P.OUT}/mf/{a.run}/run_{ds}.json"))
        FD = pd.read_parquet(f"{P.OUT}/mf/{a.run}/famday_{ds}.parquet")
        tk = J["tkeys"]; sites = J.get("sites", {}); sh = J.get("shifts", {})
        CF = {}
        for ns in tk:
            z = PC.load_zne(ns, ds)
            if z is None:
                status.append((ds, ns, "missing", "cache")); continue
            try:
                cf = p_scan.cf_day(model, meta, z["x"])
                CF[ns] = maximum_filter1d(cf, size=21, mode="constant", cval=0.0)
                status.append((ds, ns, "ok", ""))
            except Exception as e:
                status.append((ds, ns, "error", str(e)[:80]))
        rng = np.random.default_rng(1000 + int(ds.replace("-", "")))
        RI = rng.integers(0, NV, P.NRAND)
        cols = {}
        for (net, ser, f), g in D.groupby(["net", "series", "fid"]).groups.items():
            r = FD[(FD.net == net) & (FD.fid == f)]
            if not len(r):
                continue
            stas = [s for s in str(r.stations.iloc[0]).split(",") if s in CF]
            if not stas:
                continue
            smap = sites.get(net, {})
            sl = sorted({smap.get(s, s) for s in stas})
            ii = D.loc[g, "i"].to_numpy(np.int64)
            SV = np.full((len(g), len(sl)), np.nan); RV = np.full((P.NRAND, len(sl)), np.nan)
            hoS = np.zeros(len(sl), bool)
            for s in stas:
                q = sl.index(smap.get(s, s))
                d = int(round(TPL[(int(f), s, tk[s])]["off"] * W.FS))
                if ser in ("c1", "c2"):
                    j = ((ii - int(sh.get(f"{net}|{ser}|{smap.get(s, s)}", 0))) % NV) + d
                else:
                    j = ii + d
                v = CF[s][np.clip(j, 0, W.NDAY - 1)]
                SV[:, q] = np.fmax(SV[:, q], v)
                RV[:, q] = np.fmax(RV[:, q], CF[s][np.clip(RI + d, 0, W.NDAY - 1)])
                if s in HO:
                    hoS[q] = True
            rmed = np.nanmedian(RV, 0)
            b = block(D.loc[g], SV, np.ones(len(sl), bool), rmed)
            bho = block(D.loc[g], SV, hoS, rmed, "_ho")
            for k, v in {**b, **bho}.items():
                cols.setdefault(k, np.full(len(D), np.nan))[np.asarray(g)] = v
        out = pd.DataFrame(dict(k=D.k.to_numpy(), **cols))
        out.to_parquet(f"{od}/pnet_{ds}.parquet", index=False)
        del CF
        torch.cuda.empty_cache()
        if nd_ % 10 == 0:
            P.write_status(f"pnet_{a.run}", f"s{a.shard}", status)
            P.log(f"[pnet {a.run} {a.shard}] {ds} {nd_ + 1}/{len(days)} {time.time() - t0:.0f}s")
    P.write_status(f"pnet_{a.run}", f"s{a.shard}", status)
    P.log(f"[pnet {a.run} {a.shard}] done {len(days)} days {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
