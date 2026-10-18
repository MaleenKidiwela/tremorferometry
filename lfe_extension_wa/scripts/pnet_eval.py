#!/usr/bin/env python
"""Stage B evaluation: per-station window AUC of PNET on HELD-OUT days (V2008 / V2010 TEST half, A2011 all), for
every CF tag (off-the-shelf SVI weights and the WA fine-tune), with the held-out stations as a separate stratum.

  score           max of the PNET S probability within +-0.5 s of a time
  positives       labelled arrivals (pnet_common.label_pairs rule, same as the training labels) of the held-out-day
                  array detections: t + offset_s at the station's node
  negatives (q)   random times on QUIET days (V2010 'quiet', V2008 'post'), test half, >= 3 s from any known arrival
  negatives (e)   random times on EPISODE days, same halves as the positives, >= 3 s from any known arrival (a harder
                  contrast: it holds tremor and uncatalogued LFEs; it shows whether the CF only detects 'tremor on')
  FP rate         declustered (+-2 s) CF peaks >= 0.3 / 0.5 per hour on quiet days (test half)

  python scripts/pnet_eval.py --tags w05,ft1   -> results/pnet/stageB_eval_by_station.csv, stageB_eval_summary.csv,
                                                  stageB_verdict.json
"""
from __future__ import annotations
import os, sys, json, time, argparse, zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

NR = 1500
_G = {}


def _init():
    import wa_mf
    T, _ = wa_mf.load_tpl(["build"], all_templates=True)
    L = PC.label_pairs(pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv"))
    TR = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    _G.update(TOFF={k: v["off"] for k, v in T.items()}, LAB={(int(r.fid), r.netsta, r.node) for r in L.itertuples()},
              TR={d: g for d, g in TR.groupby("date")})


def job(args):
    ns, ds, w, kind, tags = args
    TOFF, LAB = _G["TOFF"], _G["LAB"]
    p = W.cache_path(ns, ds)
    if not os.path.exists(p):
        return None
    z = np.load(p)
    node = str(z["node"]); m = z["m"]
    D0 = pd.Timestamp(ds).timestamp()
    tt = D0 + np.arange(W.NDAY) / W.FS
    hv = W.half_of(tt)
    use = (hv == 1) if w in ("V2008", "V2010") else (hv >= 0)          # scored samples (test half / A2011 all)
    tg = _G["TR"].get(ds)
    lab, msk = [], []
    if tg is not None:
        g = tg if w == "A2011" else tg[tg.half == 1]
        for f, gg in tg.groupby("fid"):
            k = (int(f), ns, node)
            if k in TOFF:
                msk.append(np.round((gg.t.to_numpy() - D0 + TOFF[k]) * W.FS).astype(np.int64))
        for f, gg in g.groupby("fid"):
            k = (int(f), ns, node)
            if k in LAB:
                lab.append(np.round((gg.t.to_numpy() - D0 + TOFF[k]) * W.FS).astype(np.int64))
    lab = np.concatenate(lab) if lab else np.zeros(0, np.int64)
    msk = np.sort(np.concatenate(msk)) if msk else np.zeros(0, np.int64)
    lab = lab[(lab >= 0) & (lab < W.NDAY)]
    lab = lab[use[lab] & m[lab]]
    rng = np.random.default_rng(int(pd.Timestamp(ds).strftime("%Y%m%d")) + zlib.crc32(ns.encode()) % 1000)
    r = rng.integers(0, W.NDAY, 4 * NR)
    r = r[use[r] & m[r]]
    if len(msk):
        j = np.clip(np.searchsorted(msk, r), 1, len(msk) - 1)
        dmin = np.minimum(np.abs(msk[j] - r), np.abs(msk[j - 1] - r)) if len(msk) > 1 else np.abs(msk[0] - r)
        r = r[dmin > 60]
    r = r[:NR]
    out = dict(netsta=ns, date=ds, window=w, kind=kind, node=node, n_lab=len(lab), n_rand=len(r))
    sc = {}
    for tag in tags:
        cf = PC.load_cf(tag, ns, ds)
        if cf is None:
            continue
        mf = PC.maxfilt(cf, 10)
        sc[tag] = (mf[lab], mf[r])
        if kind in ("quiet", "post"):
            from scipy.ndimage import maximum_filter1d
            pk = (cf >= maximum_filter1d(cf, 81)) & use & m
            hrs = (use & m).sum() / W.FS / 3600
            out[f"fp03_{tag}"] = float((pk & (cf >= 0.3)).sum() / hrs) if hrs > 0 else np.nan
            out[f"fp05_{tag}"] = float((pk & (cf >= 0.5)).sum() / hrs) if hrs > 0 else np.nan
            out[f"hrs_{tag}"] = hrs
    return out, sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", required=True)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    tags = a.tags.split(",")
    S = W.load_stations()
    jobs = [(ns, d, w, k, tags) for w in ["V2008", "V2010", "A2011"] for d, k in W.window_days(w) for ns in S.netsta]
    rows, POS, NEGQ, NEGE = [], {}, {}, {}
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = [ex.submit(job, j) for j in jobs]
        for fu in as_completed(fut):
            r = fu.result()
            if r is None:
                continue
            o, sc = r
            rows.append(o)
            q = o["kind"] in ("quiet", "post")
            for tag, (sp, sr) in sc.items():
                key = (o["netsta"], tag)
                POS.setdefault(key, []).append(sp)
                (NEGQ if q else NEGE).setdefault(key, []).append(sr)
                POS.setdefault((o["netsta"], tag, o["window"]), []).append(sp)
                if not q:
                    NEGE.setdefault((o["netsta"], tag, o["window"]), []).append(sr)
    D = pd.DataFrame(rows)
    W.log(PC.LOG, f"[eval] {len(D)} station-days {time.time()-t0:.0f}s")
    cat = lambda v: np.concatenate(v) if v else np.zeros(0)
    out = []
    for ns in S.netsta:
        for tag in tags:
            p = cat(POS.get((ns, tag), [])); nq = cat(NEGQ.get((ns, tag), [])); ne = cat(NEGE.get((ns, tag), []))
            dq = D[(D.netsta == ns) & D.kind.isin(["quiet", "post"])]
            r = dict(netsta=ns, tag=tag, held_out=ns in PC.HELD_OUT, n_pos=len(p), n_neg_quiet=len(nq), n_neg_episode=len(ne),
                     auc_vs_quiet=PC.auc(p, nq), auc_vs_episode=PC.auc(p, ne),
                     pos_med=float(np.median(p)) if len(p) else np.nan, pos_frac03=float((p > 0.3).mean()) if len(p) else np.nan,
                     quiet_frac03=float((nq > 0.3).mean()) if len(nq) else np.nan,
                     fp03_per_h=float(np.average(dq[f"fp03_{tag}"], weights=dq[f"hrs_{tag}"])) if f"fp03_{tag}" in dq and dq[f"hrs_{tag}"].sum() > 0 else np.nan,
                     fp05_per_h=float(np.average(dq[f"fp05_{tag}"], weights=dq[f"hrs_{tag}"])) if f"fp05_{tag}" in dq and dq[f"hrs_{tag}"].sum() > 0 else np.nan)
            for w in ["V2008", "V2010", "A2011"]:
                pw = cat(POS.get((ns, tag, w), [])); ew = cat(NEGE.get((ns, tag, w), []))
                r[f"auc_vs_quiet_{w}"] = PC.auc(pw, nq); r[f"auc_vs_episode_{w}"] = PC.auc(pw, ew); r[f"n_pos_{w}"] = len(pw)
            out.append(r)
    O = pd.DataFrame(out)
    nh = D.groupby("netsta").node.first()
    O.to_csv(f"{PC.OUT}/stageB_eval_by_station.csv", index=False)
    ok = O[O.n_pos >= 50]
    summ = ok.groupby(["held_out", "tag"]).agg(n_stations=("netsta", "size"), med_auc_vs_quiet=("auc_vs_quiet", "median"),
                                               med_auc_vs_episode=("auc_vs_episode", "median"),
                                               med_fp03_per_h=("fp03_per_h", "median"), med_fp05_per_h=("fp05_per_h", "median"),
                                               med_pos_frac03=("pos_frac03", "median")).reset_index()
    summ.to_csv(f"{PC.OUT}/stageB_eval_summary.csv", index=False)
    v = {}
    if len(tags) >= 2:
        h = ok[ok.held_out]
        a0 = h[h.tag == tags[0]].set_index("netsta"); a1 = h[h.tag == tags[1]].set_index("netsta")
        common = a0.index.intersection(a1.index)
        d_q = (a1.auc_vs_quiet[common] - a0.auc_vs_quiet[common]).median() if len(common) else np.nan
        d_e = (a1.auc_vs_episode[common] - a0.auc_vs_episode[common]).median() if len(common) else np.nan
        v = dict(base=tags[0], ft=tags[1], held_out_stations_scored=list(common), median_delta_auc_vs_quiet=d_q,
                 median_delta_auc_vs_episode=d_e,
                 works=bool(np.isfinite(d_q) and np.isfinite(d_e) and d_q > 0.02 and d_e > 0.0))
        json.dump(v, open(f"{PC.OUT}/stageB_verdict.json", "w"), indent=1, default=float)
    W.log(PC.LOG, f"[eval]\n{summ.round(3).to_string()}\nverdict {v}")


if __name__ == "__main__":
    main()
