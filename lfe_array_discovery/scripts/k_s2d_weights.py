#!/usr/bin/env python
"""Step 2d: how much network z do we throw away?  Station set and weighting.

Dumps, for loud same-family pairs in busy family-hours, the PER-STATION 3-C cc and the matched
per-station noise distribution, over EVERY 2005 station that has a coherent corrected stack
(including the four the blind design withholds and SHVB with its 4.76 Hz tone notched out), then
evaluates the network statistic under several station sets and weightings against the 8xMAD bar.

  python scripts/k_s2d_weights.py --chain 20hz --out results/repro/s2d_20hz
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

MIN_LAG_S = 10.0
N_NOISE = 20000
HOP_S = 0.5


def norm_rows(M):
    r = np.sqrt((M ** 2).sum(1, keepdims=True))
    return np.where(r > 0, M / np.maximum(r, 1e-30), 0.0)


def notch(X, fs, f0=4.7619, bw=0.25):
    from scipy.signal import iirnotch, filtfilt
    b, a = iirnotch(f0 / (fs / 2), f0 / bw)
    return filtfilt(b, a, X, axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="20hz", choices=["20hz", "native"])
    ap.add_argument("--win", type=float, default=6.0)
    ap.add_argument("--pre", type=float, default=0.5)
    ap.add_argument("--nfam", type=int, default=8)
    ap.add_argument("--ndays", type=int, default=3)
    ap.add_argument("--maxpairs", type=int, default=800)
    ap.add_argument("--times", default=f"{K.OUT}/s1_times.csv")
    ap.add_argument("--offsets", default=f"{K.OUT}/s2a_offsets_all.csv")
    ap.add_argument("--out", default=f"{K.OUT}/s2d")
    a = ap.parse_args()
    rng = np.random.default_rng(11)
    fs = 20.0 if a.chain == "20hz" else 100.0
    W, n = a.win, int(round(a.win * (20.0 if a.chain == "20hz" else 100.0)))

    C = pd.read_csv(a.times)
    C["t_corr"] = C.t_origin + C.dt
    C = C[C.loud].copy()
    C["hh"] = pd.to_datetime(C.t_corr, unit="s").dt.hour
    OFF = pd.read_csv(a.offsets)
    OFF = OFF[OFF.peak_mad >= 4.0]
    SOFF = {(int(r.bkey), r.sta): float(r.tt_pred + r.off_s) for r in OFF.itertuples()}
    fams = [int(b) for b in C.groupby("bkey").size().sort_values(ascending=False).index[:a.nfam]]
    stas = sorted({s for (b, s) in SOFF if b in fams})
    dc = (K.DayCache(stas) if a.chain == "20hz" else K.NativeDayCache(stas))
    K.log(f"chain={a.chain} win={W} stas={stas}")

    cells = (C[C.bkey.isin(fams)].groupby(["bkey", "date", "hh"]).size().rename("n")
             .reset_index())
    cells = cells[cells.n >= 5].sort_values("n", ascending=False).groupby("bkey").head(a.ndays)
    cells = cells.sort_values(["date", "hh"]).reset_index(drop=True)
    K.log(f"{len(cells)} family-hours")

    recs = []
    for (date, hh), gg in cells.groupby(["date", "hh"]):
        D = dc.day(date)
        t_h = pd.Timestamp(f"{date} {int(hh):02d}:00:00").timestamp()
        hours = {}
        for s in stas:
            if s not in D:
                continue
            t0, X = D[s][0], D[s][1]
            i0 = int(round((t_h - t0) * fs))
            seg = np.asarray(X[:, i0:i0 + int(round(3600 * fs))], float)
            if seg.shape[1] < int(3000 * fs) or not np.isfinite(seg).all():
                continue
            if s == "SHVB":
                seg = notch(seg, fs)
            hours[s] = seg
        if len(hours) < 6:
            continue
        hop = int(round(HOP_S * fs))
        noise, nwmin = {}, None
        grids = {}
        for s in hours:
            seg = hours[s]
            nw = (seg.shape[1] - n) // hop + 1
            idx = np.arange(nw)[:, None] * hop + np.arange(n)[None, :]
            grids[s] = norm_rows(seg[:, idx].transpose(1, 0, 2).reshape(nw, 3 * n))
            nwmin = nw if nwmin is None else min(nwmin, nw)
        ia = rng.integers(0, nwmin, N_NOISE)
        ib = rng.integers(0, nwmin, N_NOISE)
        ok = np.abs(ia - ib) >= int(round(MIN_LAG_S / HOP_S))
        ia, ib = ia[ok], ib[ok]
        for s in grids:
            noise[s] = (grids[s][ia] * grids[s][ib]).sum(1)
        del grids
        for b, g in gg.groupby("bkey"):
            E = C[(C.bkey == b) & (C.date == date) & (C.hh == hh)].sort_values("t_corr")
            t = E.t_corr.to_numpy()
            sset = [s for s in hours if (b, s) in SOFF]
            if len(sset) < 6:
                continue
            pi, pj = np.triu_indices(len(t), 1)
            keep = np.abs(t[pj] - t[pi]) >= MIN_LAG_S
            pi, pj = pi[keep], pj[keep]
            if len(pi) > a.maxpairs:
                sel = rng.choice(len(pi), a.maxpairs, replace=False)
                pi, pj = pi[sel], pj[sel]
            if len(pi) < 3:
                continue
            for s in sset:
                M = np.zeros((len(t), 3 * n))
                okv = np.zeros(len(t), bool)
                for k in range(len(t)):
                    c = dc.cut(s, t[k] + SOFF[(b, s)] - a.pre, W, date)
                    if c is not None and np.isfinite(c).all() and (c ** 2).sum() > 0:
                        cc = notch(np.asarray(c, float), fs) if s == "SHVB" else c
                        M[k] = np.asarray(cc).reshape(-1)
                        okv[k] = True
                M = norm_rows(M)
                p = (M[pi] * M[pj]).sum(1)
                p[~(okv[pi] & okv[pj])] = np.nan
                nz = noise[s]
                recs.append(dict(bkey=b, date=date, hh=int(hh), sta=s, dist_km=K.epi_km(b, s),
                                 n_pair=int(np.isfinite(p).sum()),
                                 cc_med=float(np.nanmedian(p)), cc_p90=float(np.nanquantile(p, .9)),
                                 noise_mu=float(nz.mean()), noise_sd=float(nz.std()),
                                 z_med=float((np.nanmedian(p) - nz.mean()) / nz.std()),
                                 cc=p.astype(np.float32).tolist()))
        K.log(f"[w] {date} {hh:02d}h stations={len(hours)} recs={len(recs)}")
    R = pd.DataFrame(recs)
    R.drop(columns=["cc"]).to_csv(f"{a.out}_bysta.csv", index=False)
    np.save(f"{a.out}_cc.npy", np.array(R.cc.tolist(), dtype=object), allow_pickle=True)
    K.log("\nPER STATION z of the median loud pair\n" +
          str(R.groupby("sta").agg(n=("z_med", "size"), dist=("dist_km", "median"),
                                   cc=("cc_med", "median"), z=("z_med", "median")).round(3)))

    # ------------------------------------------------------------ network schemes
    from scipy.stats import norm as _nrm
    z_fa = float(_nrm.isf(1.0 / (7190 * 7189 / 2)))
    SETS = dict(clean11=K.CLEAN11,
                surface9=["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC"],
                all_avail=sorted(R.sta.unique()))
    out = []
    for (b, date, hh), g in R.groupby(["bkey", "date", "hh"]):
        cc = {r.sta: np.asarray(r.cc, float) for r in g.itertuples()}
        sd = {r.sta: r.noise_sd for r in g.itertuples()}
        mu = {r.sta: r.noise_mu for r in g.itertuples()}
        zs = {r.sta: r.z_med for r in g.itertuples()}
        nmin = min(len(v) for v in cc.values())
        for name, ss in SETS.items():
            use = [s for s in ss if s in cc]
            if len(use) < 6:
                continue
            for wname in ("equal", "snr"):
                if wname == "equal":
                    w = np.ones(len(use))
                else:
                    w = np.array([max(zs[s], 0.0) for s in use])
                    if w.sum() <= 0:
                        continue
                w = w / w.sum()
                P = np.nansum([w[i] * (cc[s][:nmin] - mu[s]) for i, s in enumerate(use)], 0)
                sdn = float(np.sqrt(np.sum([(w[i] * sd[s]) ** 2 for i, s in enumerate(use)])))
                z = P / sdn
                out.append(dict(bkey=b, date=date, hh=int(hh), set=name, weight=wname,
                                n_sta=len(use), n_pair=int(np.isfinite(z).sum()),
                                z_med=float(np.nanmedian(z)), z_p90=float(np.nanquantile(z, .9)),
                                z_p99=float(np.nanquantile(z, .99)), z_max=float(np.nanmax(z)),
                                n_ge_fa=int(np.nansum(z >= z_fa)),
                                frac_ge_fa=float(np.nanmean(z >= z_fa))))
    N = pd.DataFrame(out)
    N.to_csv(f"{a.out}_net.csv", index=False)
    K.log(f"\nNETWORK vs 8xMAD bar z_fa={z_fa:.2f}\n" +
          str(N.groupby(["set", "weight"]).agg(
              n_cells=("z_med", "size"), n_sta=("n_sta", "median"), z_med=("z_med", "median"),
              z_p99=("z_p99", "median"), z_max=("z_max", "median"),
              pairs_ge_fa=("n_ge_fa", "sum"), cells_with_hit=("n_ge_fa", lambda x: (x > 0).sum()),
              ).round(3)))


if __name__ == "__main__":
    main()
