#!/usr/bin/env python
"""ARRAY-AWARE pair stage.  Prefix `o_`.

Three SITE TREATMENTS of the same literal Bostock/Brown network autocorrelation (6-s windows,
1-8 Hz, full-sample lags, |lag| >= 6 s, median + 8 x MAD of the block's own distribution, 4-s
decluster -- every constant comes from `k_auto` and none is restated here):

  member   every admitted station is one term of the network mean.  This is what the 2010
           N-Olympic run did, and it gives a 13-element cluster 13 votes: whatever those 13
           seismometers SHARE -- cultural noise, one glitch, wind on one ridge -- adds
           coherently in the network function.
  beam     one BEAM TRACE per array (o_beam.py), so each array is one term.  This is what
           Chestler & Creager (2017) / Ghosh et al. (2012) / Sweet et al. (2014) did.
  siteavg  every member still contributes its own correlation, but the network function is the
           plain mean over SITES of the per-site mean over that site's channels.  Each site has
           exactly one vote, as in `beam`, but a member's INDEPENDENT view of the signal is kept
           -- which matters because the arrays are only ~15-25 % coherent at 1-8 Hz even for a
           catalogued deep earthquake (results/WA_2010_ARRAY_RESULTS.md sec.0c).

`member` and `beam` are k_auto.autocorr_block verbatim on different station sets.  `siteavg`
is the same computation with per-channel weights w = 1/(3 * n_members(site) * n_sites) instead
of 1/n_channels; nothing else differs, and with one member per site it reduces to `member`.

  O_BED=2010 python scripts/o_auto.py --mode siteavg --hours tremor --out <dir>
  O_BED=2010 python scripts/o_auto.py --mode member  --hours tremor --only 2010-08-14:05
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
import o_common as O                                                        # noqa: E402

MIN_STA = 4


def autocorr_weighted(D, W, fs, nmad=8.0, win_s=6.0, hop_s=0.5, min_lag_s=6.0,
                      dev="cuda", row_block=512, max_out=400000, decl_s=4.0):
    """k_auto.autocorr_block with a per-STATION weight W[sta] on the network mean.  Every other
    line -- degenerate-window QC, the MAD-of-the-block threshold, the |lag| bar, the 4-s
    declustering -- is k_auto's, copied without change."""
    import torch
    import k_auto as KA                                                     # noqa: F401
    n = int(round(win_s * fs))
    hop = int(round(hop_s * fs))
    minlag = int(round(min_lag_s * fs))
    NH = min(x.shape[1] for x in D.values())
    NV = NH - n + 1
    NW = (NH - n) // hop + 1
    A, H, E, w = {}, {}, {}, {}
    wsum = 0.0
    for s, x in D.items():
        X = torch.from_numpy(np.ascontiguousarray(x[:, :NH], dtype=np.float32)).to(dev)
        xd = X.double()
        cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                     xd * xd], 1), 1)
        ew = (cs[:, n:NH + 1] - cs[:, :NV])[:, :NV]
        if float(ew.median()) <= 0:
            del X, xd, cs, ew
            continue
        emed = ew.median(dim=1, keepdim=True).values
        good = ew > 1e-4 * emed
        a = X.unfold(-1, n, hop)[:, :NW].contiguous()
        r = a.pow(2).sum(-1, keepdim=True).sqrt()
        gr = good[:, ::hop][:, :NW].unsqueeze(-1)
        A[s] = torch.where((r > 0) & gr, a / r.clamp_min(1e-30), torch.zeros_like(a))
        H[s] = [X[c].unfold(0, n, 1)[:NV].contiguous() for c in range(3)]
        E[s] = torch.where(good, ew.clamp_min(1e-30).sqrt(),
                           torch.full_like(ew, float("inf"))).float()
        w[s] = float(W[s])
        wsum += 3.0 * w[s]
        del X, xd, cs, ew, a, r
    if not A:
        return None, None
    NW = min(NW, min(A[s].shape[1] for s in A))
    sub = 8
    rows_out, acc_stats = [], []
    for b0 in range(0, NW, row_block):
        b1 = min(b0 + row_block, NW)
        acc = torch.zeros((b1 - b0, NV), device=dev)
        for s in A:
            for c in range(3):
                acc += w[s] * ((A[s][c][b0:b1] @ H[s][c].T) / E[s][c][None])
        acc /= float(wsum)
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
    out = []
    jj = torch.arange(NV, device=dev)
    dec = int(round(decl_s * fs))
    for (b0, b1, acc) in rows_out:
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
    stats = dict(n_ch=int(3 * len(A)), n_sta=len(D), med=med, mad=mad, sd=sd, thr=thr,
                 nmad_in_sd=float(nmad * mad / sd) if sd > 0 else np.nan,
                 NW=int(NW), NV=int(NV), n_det=int(len(R)))
    return R, stats


def bed_for(mode):
    """(list of trace keys, weight per key, loader)."""
    S = O.sites()
    SM = O.site_members()
    if mode == "beam":
        keys = list(S.site)
        return keys, {k: 1.0 for k in keys}, None
    keys, W = [], {}
    for r in S.itertuples():
        mem = SM[r.site]
        for m in mem:
            keys.append(m)
            W[m] = (1.0 / len(mem)) if mode == "siteavg" else 1.0
    return keys, W, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="siteavg", choices=["member", "beam", "siteavg"])
    ap.add_argument("--hours", default="tremor", choices=["tremor", "quiet"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--keep", default="", help="restrict to these trace codes (mode comparison)")
    a = ap.parse_args()
    K = O.activate()
    import k_auto as KA
    os.makedirs(a.out, exist_ok=True)
    keys, W, _ = bed_for(a.mode)
    if a.keep:
        kk = set(x for x in a.keep.split(",") if x)
        keys = [k for k in keys if k in kk]
    if a.mode == "beam":
        load = K.load_day
    else:
        import o_beam as OB
        load = OB.load_member
    H = O.candidate_hours() if a.hours == "tremor" else O.quiet_hours()
    HL = [(str(r.date), int(r.hour)) for r in H.itertuples()]
    if a.only:
        d, h = a.only.split(":")
        HL = [(d, int(h))]
    if a.limit:
        HL = HL[:a.limit]
    veto = O.health_veto() if a.mode == "beam" else O._member_veto()
    O.log(f"[{a.mode}] {len(keys)} traces, {len(HL)} {a.hours} hours -> {a.out}")
    cache = {}
    allst = []
    for (date, hh) in HL:
        p = f"{a.out}/{date}_{hh:02d}.parquet"
        if os.path.exists(p):
            continue
        if cache.get("date") != date:
            cache = {"date": date}
            for s in keys:
                v = load(s, date)
                if v is not None:
                    cache[s] = v
        bad = veto.get((date, hh), set())
        t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
        Db = {}
        for s in keys:
            if s in bad or s not in cache:
                continue
            t0, X = cache[s]
            i0 = int(round((t_h - t0) * K.FS))
            seg = X[:, i0:i0 + int(round(3600 * K.FS))]
            if seg.shape[1] < int(3000 * K.FS) or not np.isfinite(seg).all():
                continue
            Db[s] = seg
        if len(Db) < MIN_STA:
            O.log(f"[skip] {date} {hh:02d}h {len(Db)} traces")
            continue
        try:
            if a.mode == "siteavg":
                R, st = autocorr_weighted(Db, W, K.FS)
            else:
                R, st = KA.autocorr_block(Db, K.FS)
        except RuntimeError as e:
            O.log(f"[retry] {date} {hh:02d}h {type(e).__name__}: {str(e)[:120]}")
            import time
            import torch
            torch.cuda.empty_cache()
            time.sleep(90)
            R, st = (autocorr_weighted(Db, W, K.FS) if a.mode == "siteavg"
                     else KA.autocorr_block(Db, K.FS))
        if R is None:
            continue
        hopn = int(round(KA.HOP_S * K.FS))
        R["t_i"] = t_h + R.i.to_numpy(np.float64) * hopn / K.FS
        R["t_j"] = t_h + R.j.to_numpy(np.float64) / K.FS
        R["date"], R["hour"] = date, hh
        R.to_parquet(p, index=False)
        st.update(date=date, hour=hh, n_veto=len(bad), mode=a.mode,
                  stas=",".join(sorted(Db)))
        json.dump(st, open(f"{a.out}/{date}_{hh:02d}_stat.json", "w"))
        allst.append(st)
        O.log(f"[auto:{a.mode}] {date} {hh:02d}h  {len(Db)} traces (veto {len(bad)})  "
              f"thr={st['thr']:.5f} (={st['nmad_in_sd']:.2f} sd) -> {len(R)} detections")
    sf = sorted(glob.glob(f"{a.out}/*_stat.json"))
    if sf:
        S = pd.DataFrame([json.load(open(f)) for f in sf])
        S.to_csv(f"{a.out}/_stats.csv", index=False)
        json.dump(dict(win_s=KA.WIN_S, hop_s=KA.HOP_S, nmad=KA.NMAD, min_lag_s=KA.MIN_LAG_S,
                       decluster_s=KA.DECL_S, mode=a.mode, n_traces=len(keys),
                       hours=a.hours, n_blocks=int(len(S)),
                       nmad_in_sd_median=float(S.nmad_in_sd.median())),
                  open(f"{a.out}/_params.json", "w"), indent=1)
        O.log(f"{len(S)} blocks; n_det median {S.n_det.median():.0f}; "
              f"8xMAD = {S.nmad_in_sd.median():.2f} sd")


if __name__ == "__main__":
    main()
