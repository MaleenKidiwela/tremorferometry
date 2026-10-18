#!/usr/bin/env python
"""STAGE 3 -- Royer & Bostock (2014) sec.3.1 NETWORK AUTOCORRELATION, literally.

  "Each hour-long segment was divided into 15 s windows lagged by 0.5 s to produce a total of 7170
   individual windows per hour.  Each window was correlated with all other windows in the same
   1-hour record for an individual (station-component) channel.  The resulting time-series from all
   stations and all components corresponding to the same hour were stacked to create 'network
   autocorrelation' records.  When the network autocorrelation exceeded 8 times the median absolute
   deviation we registered the corresponding window pairs as LFE detections.  This catalogue was
   then culled to exclude overlapping detections and retain only those with high SNR."

Published yields to check against: 4915 initial detections (NW), 3306 (SW).

Differences from every earlier Washington run in this repo, all of them the paper's own values:
  15 s windows (not 6 s), 40 sps (not 20), the 7 named stations of the subarray (not 51 probe-
  admitted ones), 8 x MAD over the HOUR.

Two readings of "each window was correlated with all other windows" are implemented and BOTH are
reported, because the sentence constrains the row axis (7170 windows) and not the lag axis:
  --cols window   the 0.5 s x 0.5 s window-pair matrix, 7171 x 7171  (the literal reading)
  --cols sample   every sample lag, 7171 x 143401                    (what k_auto.py did at 20 Hz)

MASK DISCIPLINE (the Stage-0 fix).  A window enters the statistic only if EVERY one of its samples
is valid in the cache mask; the network function at (i, j) is the mean over the channels valid at
that pair, with a per-pair channel count, not a constant divisor; and the median/MAD pool contains
only valid pairs.  The relative-energy guard `ew > 1e-4 * median(ew)` of k_auto.py:70 is NOT used --
it degenerates to `ew > 0` whenever a channel is more than half dead.

  python scripts/ab_auto.py --sub NW --dates 2008-05-01 --hours 0-23 --out results/ab_royer/NW/auto
-> <out>/<date>_<hh>.parquet   i, j, t_i, t_j, cc, nmad, n_ch, snr_i, snr_j
-> <out>/<date>_<hh>_stat.json med, mad, thr, n_pair_valid, n_det, n_det_predecl, channels
-> <out>/_stats.csv, <out>/_params.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

MIN_LAG_S = A.WIN_S     # two windows whose starts are closer than one window length overlap
DECL_S = A.WIN_S        # "cull to exclude overlapping detections" -- the same length, on both axes


def _chan_tensors(H, dev, n, hop, cols):
    """Per-channel row and column matrices, mask-aware.

    Returns rows[(sta,c)] (NR, n) unit-norm float32 with invalid rows zeroed,
            colsM[(sta,c)] (NC, n) float32 raw windows,
            cden[(sta,c)] (NC,) float32 column norms, +inf where the column is invalid,
            rvalid[(sta,c)] (NR,) bool, cvalid[(sta,c)] (NC,) bool."""
    import torch
    rows, colsM, cden, rvalid, cvalid = {}, {}, {}, {}, {}
    for s, (X, m) in H.items():
        nh = X.shape[1]
        NR = (nh - n) // hop + 1
        ew, vw = A.window_stats(X, m, n, hop)                    # (3, NR), (NR,)
        Xt = torch.from_numpy(np.ascontiguousarray(X, np.float32)).to(dev)
        rw = Xt.unfold(-1, n, hop)[:, :NR].contiguous()          # (3, NR, n)
        rn = rw.pow(2).sum(-1).sqrt()                            # (3, NR)
        vwt = torch.from_numpy(vw).to(dev)
        if cols == "window":
            cw, cn, cv = rw, rn, vwt
        else:
            NC = nh - n + 1
            ewc, vwc = A.window_stats(X, m, n, 1)
            cw = Xt.unfold(-1, n, 1)[:, :NC].contiguous()
            cn = torch.from_numpy(np.sqrt(ewc).astype(np.float32)).to(dev)
            cv = torch.from_numpy(vwc).to(dev)
        for c in range(3):
            ok_r = vwt & (rn[c] > 0)
            rows[(s, c)] = torch.where(ok_r[:, None], rw[c] / rn[c].clamp_min(1e-30)[:, None],
                                       torch.zeros_like(rw[c]))
            colsM[(s, c)] = cw[c]
            ok_c = cv & (cn[c] > 0)
            cden[(s, c)] = torch.where(ok_c, cn[c], torch.full_like(cn[c], float("inf")))
            rvalid[(s, c)] = ok_r
            cvalid[(s, c)] = ok_c
        del Xt, rw, rn
    return rows, colsM, cden, rvalid, cvalid


def autocorr_hour(H, cols="window", nmad=A.NMAD, win_s=A.WIN_S, hop_s=A.HOP_S,
                  min_lag_s=MIN_LAG_S, decl_s=DECL_S, row_block=1024, dev=None, mad_pool="row"):
    """H = {sta: (X (3,nh) float32, m (nh,) bool)} for ONE hour.  -> (DataFrame, stats)."""
    import torch
    dev = dev or A.torch_dev()
    n = int(round(win_s * A.FS))
    hop = int(round(hop_s * A.FS))
    nh = min(X.shape[1] for X, _ in H.values())
    H = {s: (X[:, :nh], m[:nh]) for s, (X, m) in H.items()}
    NR = (nh - n) // hop + 1
    NC = NR if cols == "window" else nh - n + 1
    rows, colsM, cden, rvalid, cvalid = _chan_tensors(H, dev, n, hop, cols)
    if not rows:
        return None, None
    keys = sorted(rows)
    # column start sample of each column index, for the overlap exclusion
    cstart = (torch.arange(NC, device=dev) * (hop if cols == "window" else 1))
    minlag = int(round(min_lag_s * A.FS))
    # ---- pass 1: median / MAD of the network function over the WHOLE hour, valid pairs only
    sub = 1 if cols == "window" else 8
    pool, blocks = [], []
    for b0 in range(0, NR, row_block):
        b1 = min(b0 + row_block, NR)
        acc = torch.zeros((b1 - b0, NC), device=dev)
        cnt = torch.zeros((b1 - b0, NC), device=dev)
        for k in keys:
            acc += (rows[k][b0:b1] @ colsM[k].T) / cden[k][None]
            cnt += (rvalid[k][b0:b1, None] & cvalid[k][None, :]).float()
        net = torch.where(cnt > 0, acc / cnt.clamp_min(1.0),
                          torch.full_like(acc, float("nan")))
        rstart = torch.arange(b0, b1, device=dev)[:, None] * hop
        far = (cstart[None, :] - rstart).abs() >= minlag
        net = torch.where(far, net, torch.full_like(net, float("nan")))
        blocks.append((b0, b1, net, cnt))
        v = net[:, ::sub].reshape(-1)
        pool.append(v[torch.isfinite(v)].clone())
        if sum(x.numel() for x in pool) > 4e8:
            pool = [torch.cat(pool)[::4]]
    V = torch.cat(pool)
    if V.numel() < 1000:
        return None, None
    med = float(V.median())
    madv = float((V - med).abs().median())
    thr = med + nmad * madv
    stats = dict(med=med, mad=madv, thr=thr, sd=float(V.std()), n_pair_valid=int(V.numel()),
                 NR=int(NR), NC=int(NC), n_ch=len(keys), n_sta=len(H), cols=cols,
                 nmad_in_sd=float(nmad * madv / float(V.std())) if float(V.std()) > 0 else np.nan)
    del V, pool
    # mad_pool="row": the MAD of EACH window's own network-autocorrelation record (row of the
    # matrix: "the resulting time-series ... were stacked ... When the network autocorrelation
    # exceeded 8 x MAD"); a pair (i, j) registers if it exceeds the bar of row i OR of row j.
    rthr = None
    if mad_pool == "row" and cols == "window":
        full = torch.cat([b[2] for b in blocks], 0)                       # (NR, NC)
        rmed = torch.nanmedian(full, dim=1).values
        rmad = torch.nanmedian((full - rmed[:, None]).abs(), dim=1).values
        rthr = rmed + nmad * rmad
        del full
    # ---- pass 2: pick pairs above the bar, then cull overlapping ones (2-D max on the i,j plane)
    dec = max(int(round(decl_s * A.FS / (hop if cols == "window" else 1))), 1)
    out, n_pre = [], 0
    for (b0, b1, net, cnt) in blocks:
        if rthr is None:
            hit = torch.nan_to_num(net, nan=-1e9) >= thr
        else:
            bar = torch.minimum(rthr[b0:b1, None], rthr[None, :NC])
            hit = torch.nan_to_num(net, nan=-1e9) >= bar
        # AUDIT FIX A1: the window-pair matrix is symmetric -- keep each pair once (column after row)
        rstart = torch.arange(b0, b1, device=dev)[:, None] * hop
        hit = hit & (cstart[None, :] > rstart)
        n_pre += int(hit.sum())
        vals = torch.where(hit, net, torch.full_like(net, -1e9))
        # separable max-pool over the column axis then the row axis = 2-D neighbourhood max
        p = torch.nn.functional.max_pool1d(vals[None], 2 * dec + 1, stride=1, padding=dec)[0]
        p = torch.nn.functional.max_pool1d(p.T[None], 2 * max(int(round(decl_s / hop_s)), 1) + 1,
                                           stride=1,
                                           padding=max(int(round(decl_s / hop_s)), 1))[0].T
        keep = hit & (vals >= p - 1e-9)
        ij = torch.nonzero(keep, as_tuple=False)
        if ij.numel():
            v = net[ij[:, 0], ij[:, 1]]
            out.append(pd.DataFrame(dict(
                i=(ij[:, 0] + b0).cpu().numpy().astype(np.int32),
                j=ij[:, 1].cpu().numpy().astype(np.int32),
                cc=v.cpu().numpy().astype(np.float32),
                nmad=((v - med) / max(madv, 1e-12)).cpu().numpy().astype(np.float32),
                n_ch=cnt[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16))))
        del net, cnt, hit, vals, p, keep, ij
    for k in keys:
        del rows[k], colsM[k], cden[k]
    if dev == "cuda":
        torch.cuda.empty_cache()
    R = (pd.concat(out, ignore_index=True) if out
         else pd.DataFrame(dict(i=np.zeros(0, np.int32), j=np.zeros(0, np.int32),
                                cc=np.zeros(0, np.float32), nmad=np.zeros(0, np.float32),
                                n_ch=np.zeros(0, np.int16))))
    stats.update(n_det=int(len(R)), n_det_predecl=int(n_pre))
    return R, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI", "MENDO"])
    ap.add_argument("--dates", default="", help="default: every episode day of the subarray")
    ap.add_argument("--hours", default="all")
    ap.add_argument("--cols", default="window", choices=["window", "sample"])
    ap.add_argument("--nmad", type=float, default=A.NMAD)
    ap.add_argument("--win", type=float, default=A.WIN_S)
    ap.add_argument("--hop", type=float, default=A.HOP_S)
    ap.add_argument("--decl", type=float, default=DECL_S)
    ap.add_argument("--min-sta", type=int, default=4,
                    help="an hour is processed only if this many of the 7 are present (reported)")
    ap.add_argument("--out", default="")
    ap.add_argument("--mad-pool", default="row", choices=["hour", "row"])   # v3: row MAD (see AB_ROYER_LOG item 4)
    a = ap.parse_args()
    out = a.out or f"{A.RES}/{a.sub}/auto_{a.cols}"
    os.makedirs(out, exist_ok=True)
    stas = A.SUBARRAY[a.sub]["stas"]
    dates = a.dates.split(",") if a.dates else A.days(a.sub)
    A.log(f"[auto {a.sub}] {len(stas)} stations {stas}; {len(dates)} dates; cols={a.cols}")
    dc = A.DayCache(stas)
    for date in dates:
        for hh in (range(24) if a.hours == "all" else
                   [int(x) for x in a.hours.split(",")] if "," in a.hours else
                   range(int(a.hours.split("-")[0]), int(a.hours.split("-")[1]) + 1)
                   if "-" in a.hours else [int(a.hours)]):
            p = f"{out}/{date}_{hh:02d}.parquet"
            if os.path.exists(p):
                continue
            H = dc.hour(date, hh)
            H = {s: v for s, v in H.items() if v[1].mean() > 0}
            if len(stas) > A.N_AUTO:
                # NVI (sec.3.1): "different combinations of 7 of the 10 three-component stations"
                # -> per hour, the first 7 of the paper's listed 10 that recorded this hour
                # (three-component, i.e. non-zero horizontals)
                ok3 = [s for s in stas if s in H and (H[s][0][1:] ** 2).sum() > 0]
                H = {s: H[s] for s in ok3[:A.N_AUTO]}
            if len(H) < a.min_sta:
                json.dump(dict(date=date, hour=hh, n_sta=len(H), status="too few stations"),
                          open(f"{out}/{date}_{hh:02d}_stat.json", "w"))
                continue
            R, st = autocorr_hour(H, cols=a.cols, nmad=a.nmad, win_s=a.win, hop_s=a.hop,
                                  decl_s=a.decl, mad_pool=a.mad_pool)
            if R is None:
                json.dump(dict(date=date, hour=hh, n_sta=len(H), status="no valid pairs"),
                          open(f"{out}/{date}_{hh:02d}_stat.json", "w"))
                continue
            t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
            hopn = int(round(a.hop * A.FS))
            R["t_i"] = t_h + R.i.to_numpy(np.float64) * hopn / A.FS
            R["t_j"] = t_h + R.j.to_numpy(np.float64) * (hopn if a.cols == "window" else 1) / A.FS
            # the SNR half of Royer's cull, measured per detection and reported, never thresholded
            cache = {}

            def net_snr(t):
                key = round(float(t) * A.FS)
                if key not in cache:
                    vals = []
                    for s in H:
                        seg = dc.cut(s, t, a.win, date)
                        if seg is not None:
                            vals.append(A.snr_window(seg))
                    cache[key] = float(np.nanmedian(vals)) if vals else np.nan
                return cache[key]

            R["snr_i"] = [net_snr(t) for t in R.t_i.to_numpy()]
            R["snr_j"] = [net_snr(t) for t in R.t_j.to_numpy()]
            R["date"], R["hour"] = date, hh
            R.to_parquet(p, index=False)
            st.update(date=date, hour=hh, stas=",".join(sorted(H)),
                      frac_valid=float(np.mean([v[1].mean() for v in H.values()])), status="ok")
            json.dump(st, open(f"{out}/{date}_{hh:02d}_stat.json", "w"))
            A.log(f"[auto] {date} {hh:02d}h {len(H)} sta thr={st['thr']:.4f} "
                  f"({st['nmad_in_sd']:.2f} sd) -> {st['n_det_predecl']} pairs, "
                  f"{len(R)} after the overlap cull")
    S = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(f"{out}/*_stat.json"))])
    S.to_csv(f"{out}/_stats.csv", index=False)
    ok = S[S.get("status", "ok") == "ok"] if "status" in S else S
    json.dump(dict(subarray=a.sub, stations=stas, fs=A.FS, band=list(A.BAND), win_s=a.win,
                   hop_s=a.hop, nmad=a.nmad, cols=a.cols, min_lag_s=MIN_LAG_S, decl_s=a.decl,
                   n_windows_per_hour=int(ok.NR.median()) if len(ok) and "NR" in ok else None,
                   n_hours=int(len(ok)), n_det_total=int(ok.n_det.sum()) if "n_det" in ok else 0,
                   published_n_det=A.SUBARRAY[a.sub]["n_det_pub"]),
              open(f"{out}/_params.json", "w"), indent=1)
    if len(ok) and "n_det" in ok:
        A.log(f"[auto {a.sub}] {len(ok)} hours -> {int(ok.n_det.sum())} detections "
              f"(Royer: {A.SUBARRAY[a.sub]['n_det_pub']})")


if __name__ == "__main__":
    main()
