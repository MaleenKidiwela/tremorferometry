#!/usr/bin/env python
"""LITERAL Brown et al. (2008) / Bostock et al. (2012) network autocorrelation.  Prefix `k_`.

Implemented from the published descriptions, with every parameter our engine invented removed:

  Bostock et al. 2012 (G-cubed 13, 2012GC004391), as described there and in Royer & Bostock 2014:
    "computing autocorrelations of 6-second long windows for each component of 7 stations"
    on 2004 and 2005 tremor -> 5,775 LFEs -> 140 families by waveform cluster analysis.
  Bostock et al. 2015 (JGR 120, 2015JB012195) sec.2: templates filtered 1-8 Hz, 6 s long,
    P on Z, S on N+E; alignment by iterated correlation against the running stack.
  Detection rule (Shelly et al. 2007, restated in Lin et al. 2024 Seismica):
    "when [the summed cross-correlation function] exceeds a threshold (e.g. EIGHT TIMES THE
     MEDIAN ABSOLUTE DEVIATION) the window is considered a detection".

So, literally:
  * 6-s windows, 1-8 Hz;
  * per STATION and per COMPONENT, the normalised cc of every window against every other window
    at full sample lag resolution;
  * the network function = the PLAIN MEAN over all channels (no per-station z, no lower-half MAD
    renormalisation, no clipping, no site merging, no >=2-sites-at-z>=3 clause, no >=4-site
    clause);
  * threshold = median + 8 * MAD of that function over the whole block (the "day's distribution");
  * no quiet-day FA calibration anywhere.

  python scripts/k_auto.py --date 2005-09-13 --hours 0-23 --out results/repro/auto
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
import k_common as K                                                        # noqa: E402

WIN_S = 6.0            # Bostock 2012 / Royer & Bostock 2014: 6-second windows
HOP_S = 0.5            # row step (Brown: "running" autocorrelation; columns are 1 sample)
MIN_LAG_S = 6.0        # exclude the trivial self-overlap: one window length
NMAD = 8.0             # Shelly 2007 / Bostock: 8 x MAD
DECL_S = 4.0           # Bostock 2015: strongest detection in a 4-s interval is retained


def autocorr_block(D, fs, nmad=NMAD, win_s=WIN_S, hop_s=HOP_S, min_lag_s=MIN_LAG_S,
                   dev="cuda", row_block=512, max_out=400000, decl_s=DECL_S):
    """D = {sta: (3, N)} one block, same absolute time base.  Returns DataFrame(i, j, cc)."""
    import torch
    n = int(round(win_s * fs))
    hop = int(round(hop_s * fs))
    minlag = int(round(min_lag_s * fs))
    NH = min(x.shape[1] for x in D.values())
    NV = NH - n + 1
    NW = (NH - n) // hop + 1
    A, H, E = {}, {}, {}
    nch = 0
    for s, x in D.items():
        X = torch.from_numpy(np.ascontiguousarray(x[:, :NH], dtype=np.float32)).to(dev)
        xd = X.double()
        cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                     xd * xd], 1), 1)
        ew = (cs[:, n:NH + 1] - cs[:, :NV])[:, :NV]                 # (3, NV) window energies
        if float(ew.median()) <= 0:
            del X, xd, cs, ew
            continue
        # DATA QC (not a threshold invention): a window that is a gap / dead channel has
        # (near-)zero energy; normalising it gives a unit-norm noise vector whose cc with another
        # such window is ~1, and 7 of the 576 ETS hours blew the whole block's sd up to 470 and
        # produced 1.16 M of 1.32 M "detections".  Degenerate windows are excluded on BOTH axes.
        emed = ew.median(dim=1, keepdim=True).values
        good = ew > 1e-4 * emed                                     # (3, NV)
        a = X.unfold(-1, n, hop)[:, :NW].contiguous()               # (3, NW, n)
        r = a.pow(2).sum(-1, keepdim=True).sqrt()
        gr = good[:, ::hop][:, :NW].unsqueeze(-1)
        A[s] = torch.where((r > 0) & gr, a / r.clamp_min(1e-30), torch.zeros_like(a))
        H[s] = [X[c].unfold(0, n, 1)[:NV].contiguous() for c in range(3)]
        E[s] = torch.where(good, ew.clamp_min(1e-30).sqrt(),
                           torch.full_like(ew, float("inf"))).float()
        nch += 3
        del X, xd, cs, ew, a, r
    if not A:
        return None, None
    NW = min(NW, min(A[s].shape[1] for s in A))
    # ---- pass 1: median and MAD of the network function over the WHOLE block (subsampled cols)
    sub = 8
    samp = []
    rows_out = []
    acc_stats = []
    for b0 in range(0, NW, row_block):
        b1 = min(b0 + row_block, NW)
        nb = b1 - b0
        acc = torch.zeros((nb, NV), device=dev)
        for s in A:
            for c in range(3):
                acc += (A[s][c][b0:b1] @ H[s][c].T) / E[s][c][None]
        acc /= float(nch)
        acc_stats.append(acc[:, ::sub].reshape(-1).clone())
        rows_out.append((b0, b1, acc))
        if sum(x.numel() for x in acc_stats) > 6e8:
            acc_stats = [torch.cat(acc_stats)[::4]]
    V = torch.cat(acc_stats)
    med = float(V.median())
    mad = float((V - med).abs().median())
    thr = med + nmad * mad
    sd = float(V.std())
    del V, acc_stats
    # ---- pass 2: pick detections
    out = []
    jj = torch.arange(NV, device=dev)
    dec = int(round(decl_s * fs))
    for (b0, b1, acc) in rows_out:
        nb = b1 - b0
        ti = torch.arange(b0, b1, device=dev)[:, None] * hop
        m = ((jj[None, :] - ti).abs() >= minlag) & (acc >= thr)
        vals = torch.where(m, acc, torch.full((), -1e9, device=dev))
        pooled = torch.nn.functional.max_pool1d(vals[None], kernel_size=2 * dec + 1,
                                                stride=1, padding=dec)[0]
        keep = m & (vals >= pooled - 1e-9)
        ij = torch.nonzero(keep, as_tuple=False)
        if ij.numel():
            v = acc[ij[:, 0], ij[:, 1]]
            out.append(pd.DataFrame(dict(i=(ij[:, 0] + b0).cpu().numpy().astype(np.int32),
                                         j=ij[:, 1].cpu().numpy().astype(np.int32),
                                         cc=v.cpu().numpy().astype(np.float32),
                                         nmad=((v - med) / mad).cpu().numpy().astype(np.float32))))
        del acc, m, vals, pooled, keep, ij
    for s in list(A):
        del A[s], H[s], E[s]
    torch.cuda.empty_cache()
    R = pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["i", "j", "cc", "nmad"])
    if len(R) > max_out:
        R = R.nlargest(max_out, "nmad")
    stats = dict(n_ch=nch, n_sta=len(D), med=med, mad=mad, sd=sd, thr=thr,
                 nmad_in_sd=float(nmad * mad / sd) if sd > 0 else np.nan,
                 NW=int(NW), NV=int(NV), n_det=int(len(R)))
    return R, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dates", default="2005-09-13")
    ap.add_argument("--hours", default="all")
    ap.add_argument("--stas", default="")
    ap.add_argument("--nsta", type=int, default=0, help="keep the N nearest to the family field")
    ap.add_argument("--block", default="hour", choices=["hour", "day"])
    ap.add_argument("--nmad", type=float, default=NMAD)
    ap.add_argument("--win", type=float, default=WIN_S)
    ap.add_argument("--hop", type=float, default=HOP_S)
    ap.add_argument("--out", default=f"{K.OUT}/auto")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stas = a.stas.split(",") if a.stas else list(K.CLEAN11)
    if a.nsta:
        F = K.families()
        la, lo = float(F.lat.median()), float(F.lon.median())
        import s_common as S
        ST = K.station_table()
        stas = sorted(stas, key=lambda s: S.epi_km(la, lo, ST.lat[s], ST.lon[s]))[:a.nsta]
    dates = a.dates.split(",")
    K.log(f"stations ({len(stas)}): {stas}")
    dcache = K.DayCache(stas)
    allst = []
    for date in dates:
        D0 = dcache.day(date)
        have = [s for s in stas if s in D0]
        if len(have) < 4:
            K.log(f"[skip] {date}: {len(have)} stations")
            continue
        hours = (range(24) if a.hours == "all" else
                 [int(x) for x in a.hours.split(",")] if "," in a.hours else
                 range(int(a.hours.split("-")[0]), int(a.hours.split("-")[1]) + 1)
                 if "-" in a.hours else [int(a.hours)])
        for hh in hours:
            p = f"{a.out}/{date}_{hh:02d}.parquet"
            if os.path.exists(p):
                continue
            t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
            Db = {}
            for s in have:
                t0, X = D0[s]
                i0 = int(round((t_h - t0) * K.FS))
                seg = X[:, i0:i0 + int(round(3600 * K.FS))]
                if seg.shape[1] < int(3000 * K.FS) or not np.isfinite(seg).all():
                    continue
                Db[s] = seg
            if len(Db) < 4:
                continue
            R, st = autocorr_block(Db, K.FS, nmad=a.nmad, win_s=a.win, hop_s=a.hop)
            if R is None:
                continue
            hopn = int(round(a.hop * K.FS))
            R["t_i"] = t_h + R.i.to_numpy(np.float64) * hopn / K.FS
            R["t_j"] = t_h + R.j.to_numpy(np.float64) / K.FS
            R["date"], R["hour"] = date, hh
            R.to_parquet(p, index=False)
            st.update(date=date, hour=hh, stas=",".join(sorted(Db)))
            json.dump(st, open(f"{a.out}/{date}_{hh:02d}_stat.json", "w"))
            allst.append(st)
            K.log(f"[auto] {date} {hh:02d}h  {len(Db)} sta  thr={st['thr']:.4f} "
                  f"(={st['nmad_in_sd']:.2f} sd)  -> {len(R)} detections")
    sf = sorted(glob.glob(f"{a.out}/*_stat.json"))
    if sf:
        allst = [json.load(open(f)) for f in sf]
    if allst:
        S = pd.DataFrame(allst)
        S.to_csv(f"{a.out}/_stats.csv", index=False)
        K.log("\n" + str(S[["date", "hour", "n_sta", "thr", "nmad_in_sd", "n_det"]].describe()))
        json.dump(dict(win_s=a.win, hop_s=a.hop, nmad=a.nmad, stations=stas,
                       min_lag_s=MIN_LAG_S, decluster_s=DECL_S,
                       nmad_in_sd_median=float(S.nmad_in_sd.median())),
                  open(f"{a.out}/_params.json", "w"), indent=1)


if __name__ == "__main__":
    main()
