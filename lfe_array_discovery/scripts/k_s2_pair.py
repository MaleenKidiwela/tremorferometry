#!/usr/bin/env python
"""Step 2 of the reproduce-Bostock loop: can the PAIR statistic see the loud tail?

For same-family LOUD event pairs inside one hour, at the Step-1 CORRECTED times, measure the
per-station 3-C cc and the network mean against the noise distribution of random window pairs in
the SAME hour (identical cut, identical lengths), in two chains and two cut rules:

  chain   20hz    the frozen cache, 20 Hz, 1-8 Hz   |  native  raw 100 Hz miniSEED, 2-8 Hz
  cut     moveout each station cut at its own predicted S   (oracle)
          absolute every station cut at the SAME absolute time (what w_autocorr does)

  python scripts/k_s2_pair.py --nfam 8 --chain 20hz --out results/repro/s2_20hz
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

PRE_W = 1.0          # the S arrival sits this far into the window
MIN_LAG_S = 10.0
N_NOISE = 4000


def norm_rows(Wm):
    r = np.sqrt((Wm ** 2).sum(1, keepdims=True))
    return np.where(r > 0, Wm / np.maximum(r, 1e-30), 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfam", type=int, default=8)
    ap.add_argument("--chain", default="20hz", choices=["20hz", "native"])
    ap.add_argument("--wins", default="6,8,10")
    ap.add_argument("--ndays", type=int, default=4)
    ap.add_argument("--maxpairs", type=int, default=400, help="per family per hour")
    ap.add_argument("--times", default=f"{K.OUT}/s1_times.csv")
    ap.add_argument("--loud", default=f"{K.OUT}/s1_loud.csv")
    ap.add_argument("--out", default=f"{K.OUT}/s2")
    a = ap.parse_args()
    wins = [float(x) for x in a.wins.split(",")]
    rng = np.random.default_rng(20260920)

    C = pd.read_csv(a.times)
    L = pd.read_csv(a.loud)
    C = C.merge(L[["bkey", "t_origin", "excess"]], on=["bkey", "t_origin"], how="left")
    C["t_corr"] = C.t_origin + C.dt
    C = C[C.loud].copy()
    fams = [int(b) for b in C.groupby("bkey").size().sort_values(ascending=False).index[:a.nfam]]
    C = C[C.bkey.isin(fams)].copy()
    C["hh"] = pd.to_datetime(C.t_corr, unit="s").dt.hour
    stas = K.CLEAN11
    TT = K.fam_tt(stas)
    fs = 20.0 if a.chain == "20hz" else 100.0
    dc = (K.DayCache(stas) if a.chain == "20hz" else K.NativeDayCache(stas))
    K.log(f"chain={a.chain} fs={fs} fams={fams} loud events={len(C)}")

    # busiest (family, date, hour) cells
    cells = (C.groupby(["bkey", "date", "hh"]).size().rename("n").reset_index())
    cells = cells[cells.n >= 4]
    cells = (cells.sort_values("n", ascending=False).groupby("bkey").head(a.ndays)
             .sort_values(["date", "hh"]).reset_index(drop=True))
    K.log(f"{len(cells)} family-hours, {int(cells.n.sum())} loud events")

    rows, netrows = [], []
    for (date, hh), gg in cells.groupby(["date", "hh"]):
        D = dc.day(date)
        have = [s for s in stas if s in D]
        if len(have) < 6:
            continue
        t_h = pd.Timestamp(f"{date} {int(hh):02d}:00:00").timestamp()
        # noise window bank per station per window length, on a 0.5-s grid inside the hour
        bank = {}
        for W in wins:
            n = int(round(W * fs))
            hop = int(round(0.5 * fs))
            for s in have:
                t0, X = D[s][0], D[s][1]
                i0 = int(round((t_h - t0) * fs))
                seg = X[:, i0:i0 + int(round(3600 * fs))]
                if seg.shape[1] < n + 10 * hop or not np.isfinite(seg).all():
                    continue
                nw = (seg.shape[1] - n) // hop + 1
                idx = np.arange(nw)[:, None] * hop + np.arange(n)[None, :]
                Wm = seg[:, idx].transpose(1, 0, 2).reshape(nw, 3 * n)
                bank[(s, W)] = (norm_rows(Wm.astype(np.float64)), nw, hop)
        if not bank:
            continue
        # --- noise: random (a, b) pairs with |lag| >= 10 s, SHARED across stations
        nwmin = min(v[1] for v in bank.values())
        ia = rng.integers(0, nwmin, N_NOISE)
        ib = rng.integers(0, nwmin, N_NOISE)
        ok = np.abs(ia - ib) >= int(round(MIN_LAG_S / 0.5))
        ia, ib = ia[ok], ib[ok]
        noise = {}
        for (s, W), (Wm, nw, hop) in bank.items():
            noise[(s, W)] = (Wm[ia] * Wm[ib]).sum(1)
        for b, g in gg.groupby("bkey"):
            E = C[(C.bkey == b) & (C.date == date) & (C.hh == hh)].sort_values("t_corr")
            t = E.t_corr.to_numpy()
            pi, pj = np.triu_indices(len(t), 1)
            keep = np.abs(t[pj] - t[pi]) >= MIN_LAG_S
            pi, pj = pi[keep], pj[keep]
            if len(pi) > a.maxpairs:
                sel = rng.choice(len(pi), a.maxpairs, replace=False)
                pi, pj = pi[sel], pj[sel]
            if not len(pi):
                continue
            ttmed = float(np.median([TT[(b, s)] for s in have]))
            for W in wins:
                n = int(round(W * fs))
                for cut in ("moveout", "absolute"):
                    per = {}
                    for s in have:
                        if (s, W) not in bank:
                            continue
                        off = (TT[(b, s)] if cut == "moveout" else ttmed) - PRE_W
                        segs = []
                        for k in range(len(t)):
                            c = dc.cut(s, t[k] + off, W, date)
                            segs.append(None if c is None or not np.isfinite(c).all()
                                        else c.reshape(-1))
                        okv = np.array([x is not None for x in segs])
                        if okv.sum() < 2:
                            continue
                        M = np.zeros((len(t), 3 * n))
                        M[okv] = np.stack([x for x in segs if x is not None])
                        M = norm_rows(M)
                        cc = (M[pi] * M[pj]).sum(1)
                        cc[~(okv[pi] & okv[pj])] = np.nan
                        per[s] = cc
                        nz = noise[(s, W)]
                        p999 = float(np.quantile(nz, 0.999))
                        v = cc[np.isfinite(cc)]
                        rows.append(dict(bkey=b, date=date, hh=int(hh), sta=s, win=W, cut=cut,
                                         chain=a.chain, n_pair=int(v.size),
                                         cc_med=float(np.median(v)) if v.size else np.nan,
                                         cc_p90=float(np.quantile(v, 0.9)) if v.size else np.nan,
                                         noise_med=float(np.median(nz)),
                                         noise_p999=p999,
                                         frac_gt_p999=float((v > p999).mean()) if v.size else
                                         np.nan,
                                         dist_km=K.epi_km(b, s)))
                    if len(per) >= 6:
                        Pm = np.vstack([per[s] for s in sorted(per)])
                        net = np.nanmean(Pm, 0)
                        Nz = np.vstack([noise[(s, W)] for s in sorted(per)]).mean(0)
                        p999 = float(np.quantile(Nz, 0.999))
                        v = net[np.isfinite(net)]
                        netrows.append(dict(bkey=b, date=date, hh=int(hh), win=W, cut=cut,
                                            chain=a.chain, n_sta=len(per), n_pair=int(v.size),
                                            net_med=float(np.median(v)),
                                            net_p90=float(np.quantile(v, 0.9)),
                                            noise_med=float(np.median(Nz)),
                                            noise_p999=p999,
                                            noise_sd=float(Nz.std()),
                                            z_med=float((np.median(v) - Nz.mean()) / Nz.std()),
                                            frac_gt_p999=float((v > p999).mean())))
        K.log(f"[pair] {date} {hh:02d}h  stations={len(have)} rows={len(rows)}")
    R = pd.DataFrame(rows)
    N = pd.DataFrame(netrows)
    R.to_csv(f"{a.out}_bysta.csv", index=False)
    N.to_csv(f"{a.out}_net.csv", index=False)
    if len(N):
        piv = N.groupby(["cut", "win"]).agg(
            n=("frac_gt_p999", "size"), frac=("frac_gt_p999", "median"),
            frac_mean=("frac_gt_p999", "mean"), z=("z_med", "median"),
            net_med=("net_med", "median"), noise_p999=("noise_p999", "median"))
        K.log("\nNETWORK (11 clean sites)\n" + str(piv))
        ps = R.groupby(["cut", "win"]).agg(frac=("frac_gt_p999", "median"),
                                           cc=("cc_med", "median"))
        K.log("\nPER STATION\n" + str(ps))
        json.dump(dict(chain=a.chain, n_cells=len(cells),
                       net=piv.reset_index().to_dict("records")),
                  open(f"{a.out}_summary.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
