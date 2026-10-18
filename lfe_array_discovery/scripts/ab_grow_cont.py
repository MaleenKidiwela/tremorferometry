#!/usr/bin/env python
"""STAGE 5 CONTINUATION -- carry ab_grow.py's iteration loop past its --max-iters 8 cap, from the
saved iteration-K state, until the Bostock (2012) / Royer & Bostock (2014) sec.3.3 stopping rule:
"iterate until detection count and template SNR stop improving".

Method = ab_grow.py's, unchanged: the templates of iteration K+1 are restacked from det_it<K>
(ab_grow.restack's arithmetic: every pool station, 44 s long window at LONG_PRE, unit-max-3C,
gap-guard dead channels zeroed, cnt >= 2) then cut by ab_grow.cut_templates (top-10 SNR 3C
stations, common 15 s window), scanned by ab_grow.matched_filter (8 x MAD / 24 h), merged by
ab_grow.merge_duplicates.  The ONLY change is speed: the restack is run over date-chunks in
parallel worker processes (sums + counts per chunk, then added) -- identical sums up to float
summation order.

Stopping rule (both must hold at iteration K vs K-1):
    |N_K - N_{K-1}| / N_{K-1} < --tol (0.01)   AND   (SNR_K - SNR_{K-1}) / SNR_{K-1} < --snr-tol (0.01)
with N = merged detection count and SNR = median template SNR (iters.csv tpl_snr_median).
Stops also at --max-more extra iterations (default 6).

Outputs continue the numbering in the SAME folder (no existing file is modified):
  tpl_it<K>.npz, tpl_it<K>_stations.csv, detraw_it<K>.parquet, det_it<K>.parquet, families_it<K>.csv
  iters_cont.csv  (rows 0..K_start copied from iters.csv + the new rows, + converged flags)
  summary_cont.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_grow as G                                                         # noqa: E402


def _restack_chunk(args):
    """ab_grow.restack's loop body on a subset of dates -> raw sums and counts."""
    stas, det, long_pre, long_len = args
    dc = A.DayCache(stas)
    nL = int(round(long_len * A.FS))
    acc, cnt = {}, {}
    rc = A.RefCache(dc, int(round(G.TPL_LEN * A.FS)))
    for date, dd in det.sort_values("t").groupby("date"):
        refs = {s: rc.get(s, date) for s in stas}
        for r in dd.itertuples():
            g = int(r.grp)
            for s in stas:
                seg = dc.cut(s, r.t - long_pre, long_len, date)
                if seg is None:
                    continue
                okc = refs[s] > 0
                if not okc.any():
                    continue
                if not okc.all():
                    seg = np.where(okc[:, None], seg, 0.0)
                key = (g, s)
                if key not in acc:
                    acc[key] = np.zeros((3, nL))
                    cnt[key] = 0
                acc[key] += G.unit_max3(seg[:, :nL])
                cnt[key] += 1
    return acc, cnt


def restack_parallel(stas, det, workers=16, long_pre=G.LONG_PRE, long_len=G.LONG_LEN,
                     return_counts=False):
    dates = sorted(det.date.unique())
    chunks = [dates[i::workers] for i in range(workers)]
    jobs = [(stas, det[det.date.isin(c)], long_pre, long_len) for c in chunks if c]
    acc, cnt = {}, {}
    with ProcessPoolExecutor(len(jobs)) as ex:
        for a, c in ex.map(_restack_chunk, jobs):
            for k, v in a.items():
                if k in acc:
                    acc[k] += v
                    cnt[k] += c[k]
                else:
                    acc[k], cnt[k] = v, c[k]
    ST = {}
    for (g, s), v in acc.items():
        if cnt[(g, s)] >= 2:
            ST.setdefault(g, {})[s] = v / cnt[(g, s)]
    return (ST, cnt) if return_counts else ST


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default="NW")
    ap.add_argument("--dir", default="", help="grow folder holding det_it<K>/iters.csv")
    ap.add_argument("--from-iter", type=int, default=8)
    ap.add_argument("--max-more", type=int, default=6)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--snr-tol", type=float, default=0.01)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--dates", default="")
    a = ap.parse_args()
    out = a.dir or f"{A.RES}/{a.sub}/grow"
    S0 = json.load(open(f"{out}/summary.json"))
    stas = list(S0["pool"])
    dates = a.dates.split(",") if a.dates else A.scan_days(a.sub)
    it0 = a.from_iter
    rows = pd.read_csv(f"{out}/iters.csv")
    rows = rows[rows.iter <= it0].to_dict("records")
    det = pd.read_parquet(f"{out}/det_it{it0}.parquet")
    prevN = len(det)
    prevS = float(rows[-1]["tpl_snr_median"])
    dc = A.DayCache(stas)
    A.log(f"[grow-cont {a.sub}] from it{it0}: {det.grp.nunique()} families, {prevN} detections, "
          f"SNR {prevS:.3f}; pool {len(stas)}; {len(dates)} dates")
    converged = False
    for it in range(it0 + 1, it0 + 1 + a.max_more):
        t0 = time.time()
        ST = restack_parallel(stas, det, workers=a.workers)
        A.log(f"[grow-cont it{it}] restack {len(ST)} families ({time.time() - t0:.0f} s)")
        TPL, info = G.cut_templates(ST)
        info.to_csv(f"{out}/tpl_it{it}_stations.csv", index=False)
        if not TPL:
            A.log("[grow-cont] no scannable template")
            break
        np.savez_compressed(f"{out}/tpl_it{it}.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(json.dumps(
                                {f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                 for g in TPL for s in TPL[g]})))
        t0 = time.time()
        raw = G.matched_filter(dc, TPL, dates)
        raw.to_parquet(f"{out}/detraw_it{it}.parquet", index=False)
        A.log(f"[grow-cont it{it}] matched filter {len(raw)} raw ({time.time() - t0:.0f} s)")
        det, fam, mst = G.merge_duplicates(raw)
        fam.to_csv(f"{out}/families_it{it}.csv", index=False)
        det.to_parquet(f"{out}/det_it{it}.parquet", index=False)
        per = det.groupby("grp").size()
        snr = float(info.snr_med.median())
        dN = abs(len(det) - prevN) / prevN
        dS = (snr - prevS) / prevS
        converged = bool(dN < a.tol and dS < a.snr_tol)
        rows.append(dict(iter=it, n_templates=len(TPL), n_det_raw=int(len(raw)),
                         n_dup_edges=int(mst.get("n_dup_edges", 0)),
                         n_merge_groups=int(mst.get("n_components_merged", 0)),
                         n_templates_merged_away=int(mst.get("n_templates_merged_away", 0)),
                         n_families=int(per.size), n_detections=int(len(det)),
                         n_templates_with_det=int(per.size),
                         per_tpl_median=float(per.median()), per_tpl_max=int(per.max()),
                         tpl_snr_median=snr, rel_dN=dN, rel_dSNR=dS, stop_rule_met=converged))
        pd.DataFrame(rows).to_csv(f"{out}/iters_cont.csv", index=False)
        A.log(f"[grow-cont it{it}] {len(TPL)} templates -> {per.size} families, {len(det)} det; "
              f"|dN|/N {dN:.4f}, dSNR/SNR {dS:+.4f}, SNR {snr:.3f} -> stop={converged}")
        prevN, prevS = len(det), snr
        if converged or len(det) == 0:
            break
    last = int(rows[-1]["iter"])
    json.dump(dict(subarray=a.sub, pool=stas, dates=len(dates), from_iter=it0, final_iter=last,
                   converged=converged, tol=a.tol, snr_tol=a.snr_tol, max_more=a.max_more,
                   n_final_families=int(rows[-1]["n_families"]),
                   n_final_detections=int(rows[-1]["n_detections"]),
                   n_published=A.SUBARRAY[a.sub]["n_tpl_pub"], iterations=rows),
              open(f"{out}/summary_cont.json", "w"), indent=1, default=float)
    open(f"{out}/FINAL_ITER", "w").write(f"{last}\n")
    A.log(f"[grow-cont] final it{last}, converged={converged}")


if __name__ == "__main__":
    main()
