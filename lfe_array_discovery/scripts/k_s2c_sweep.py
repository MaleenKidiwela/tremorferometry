#!/usr/bin/env python
"""Step 2c: parameter sweep of the PAIR statistic and its matched-filter ceiling.

Everything is cut at the EMPIRICAL S offset read off the Step-2a corrected stack (not fam_tt),
and every event time is the Step-1 corrected time.  For each (window length, component set,
chain) we measure, on the same loud events and in the same hours:

  tpl   cc of the family's own stack template against each single event  (matched-filter ceiling)
  pair  cc of two single events of the family against each other        (what autocorrelation has)
  noise cc of random window pairs in the same hour                      (the null both are judged against)

and report the network mean over the clean sites as a z against the network noise, plus the
threshold an FA rule of 1 false pair per hour would demand.

  python scripts/k_s2c_sweep.py --chain 20hz --out results/repro/s2c_20hz
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="20hz", choices=["20hz", "native"])
    ap.add_argument("--wins", default="1.5,2,3,4,6,10")
    ap.add_argument("--comps", default="ZNE,NE")
    ap.add_argument("--nfam", type=int, default=8)
    ap.add_argument("--ndays", type=int, default=3)
    ap.add_argument("--pre", type=float, default=0.5, help="S peak this far into the window")
    ap.add_argument("--maxpairs", type=int, default=600)
    ap.add_argument("--times", default=f"{K.OUT}/s1_times.csv")
    ap.add_argument("--offsets", default=f"{K.OUT}/s2a_offsets.csv")
    ap.add_argument("--out", default=f"{K.OUT}/s2c")
    a = ap.parse_args()
    wins = [float(x) for x in a.wins.split(",")]
    comps = a.comps.split(",")
    CIDX = dict(ZNE=[0, 1, 2], NE=[1, 2], Z=[0])
    rng = np.random.default_rng(7)

    fs = 20.0 if a.chain == "20hz" else 100.0
    C = pd.read_csv(a.times)
    C["t_corr"] = C.t_origin + C.dt
    C = C[C.loud].copy()
    C["hh"] = pd.to_datetime(C.t_corr, unit="s").dt.hour
    OFF = pd.read_csv(a.offsets)
    # only stations whose corrected stack actually shows an S (peak >= 5 MAD)
    OFF = OFF[OFF.peak_mad >= 5.0]
    SOFF = {(int(r.bkey), r.sta): float(r.tt_pred + r.off_s) for r in OFF.itertuples()}
    fams = [int(b) for b in C.groupby("bkey").size().sort_values(ascending=False).index[:a.nfam]]
    stas = sorted({s for (b, s) in SOFF if b in fams})
    dc = (K.DayCache(stas) if a.chain == "20hz" else K.NativeDayCache(stas))
    K.log(f"chain={a.chain} fams={fams} stas={stas}")

    # templates from ALL corrected-time events (a single event contributes 1/N)
    TPL = {}
    Wmax = max(wins)
    ntp = int(round((Wmax + 4) * fs))
    accT = {}
    for date, g in C.groupby("date"):
        dc.day(date)
        for r in g.itertuples():
            b = int(r.bkey)
            for s in stas:
                if (b, s) not in SOFF:
                    continue
                c = dc.cut(s, r.t_corr + SOFF[(b, s)] - a.pre - 2.0, Wmax + 4, date)
                if c is None or (c ** 2).sum() <= 0:
                    continue
                A = accT.setdefault((b, s), [np.zeros((3, ntp)), 0])
                A[0] += K.norm3(c)
                A[1] += 1
    for k, (Sm, N) in accT.items():
        if N >= 30:
            TPL[k] = Sm / N
    K.log(f"{len(TPL)} templates")

    cells = (C.groupby(["bkey", "date", "hh"]).size().rename("n").reset_index())
    cells = cells[cells.n >= 4].sort_values("n", ascending=False).groupby("bkey").head(a.ndays)
    cells = cells.sort_values(["date", "hh"]).reset_index(drop=True)

    rows = []
    for (date, hh), gg in cells.groupby(["date", "hh"]):
        D = dc.day(date)
        t_h = pd.Timestamp(f"{date} {int(hh):02d}:00:00").timestamp()
        have = [s for s in stas if s in D]
        if len(have) < 6:
            continue
        hours = {}
        for s in have:
            t0, X = D[s][0], D[s][1]
            i0 = int(round((t_h - t0) * fs))
            seg = X[:, i0:i0 + int(round(3600 * fs))]
            if seg.shape[1] < int(3000 * fs) or not np.isfinite(seg).all():
                continue
            hours[s] = seg
        if len(hours) < 6:
            continue
        for W in wins:
            n = int(round(W * fs))
            hop = int(round(HOP_S * fs))
            grid = {}
            for s in hours:
                seg = hours[s]
                nw = (seg.shape[1] - n) // hop + 1
                idx = np.arange(nw)[:, None] * hop + np.arange(n)[None, :]
                grid[s] = seg[:, idx].transpose(1, 0, 2)          # (nw, 3, n)
            nwmin = min(v.shape[0] for v in grid.values())
            ia = rng.integers(0, nwmin, N_NOISE)
            ib = rng.integers(0, nwmin, N_NOISE)
            ok = np.abs(ia - ib) >= int(round(MIN_LAG_S / HOP_S))
            ia, ib = ia[ok], ib[ok]
            for cname in comps:
                ci = CIDX[cname]
                noise, noiseT = {}, {}
                for s in grid:
                    G = norm_rows(grid[s][:, ci, :].reshape(grid[s].shape[0], -1))
                    noise[s] = (G[ia] * G[ib]).sum(1)
                    noiseT[s] = G
                for b, g in gg.groupby("bkey"):
                    E = C[(C.bkey == b) & (C.date == date) & (C.hh == hh)].sort_values("t_corr")
                    t = E.t_corr.to_numpy()
                    sset = [s for s in grid if (b, s) in TPL]
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
                    ccp, cct, nzs = [], [], []
                    for s in sset:
                        M = np.zeros((len(t), len(ci) * n))
                        okv = np.zeros(len(t), bool)
                        for k in range(len(t)):
                            c = dc.cut(s, t[k] + SOFF[(b, s)] - a.pre, W, date)
                            if c is not None and np.isfinite(c).all() and (c ** 2).sum() > 0:
                                M[k] = c[ci].reshape(-1)
                                okv[k] = True
                        M = norm_rows(M)
                        tpl = TPL[(b, s)][ci, int(round(2.0 * fs)):int(round(2.0 * fs)) + n]
                        tpl = tpl.reshape(1, -1)
                        tpl = norm_rows(tpl)[0]
                        p = (M[pi] * M[pj]).sum(1)
                        p[~(okv[pi] & okv[pj])] = np.nan
                        tt = M @ tpl
                        tt[~okv] = np.nan
                        ccp.append(p)
                        cct.append(tt)
                        nzs.append(noise[s])
                    P = np.nanmean(np.vstack(ccp), 0)
                    T = np.nanmean(np.vstack(cct), 0)
                    NZ = np.vstack(nzs).mean(0)
                    mu, sd = float(NZ.mean()), float(NZ.std())
                    # FA rule: <= 1 false pair per hour over all window pairs in the hour
                    npair_hour = nwmin * (nwmin - 1) / 2
                    from scipy.stats import norm as _nrm
                    z_fa = float(_nrm.isf(1.0 / npair_hour))
                    rows.append(dict(chain=a.chain, win=W, comp=cname, bkey=b, date=date,
                                     hh=int(hh), n_sta=len(sset), n_ev=len(t), n_pair=len(pi),
                                     pair_med=float(np.nanmedian(P)),
                                     pair_p90=float(np.nanquantile(P, 0.9)),
                                     pair_max=float(np.nanmax(P)),
                                     tpl_med=float(np.nanmedian(T)),
                                     tpl_p90=float(np.nanquantile(T, 0.9)),
                                     noise_mu=mu, noise_sd=sd, z_fa=z_fa,
                                     pair_z_med=(float(np.nanmedian(P)) - mu) / sd,
                                     pair_z_p90=(float(np.nanquantile(P, 0.9)) - mu) / sd,
                                     pair_z_max=(float(np.nanmax(P)) - mu) / sd,
                                     frac_pair_ge_fa=float(np.nanmean((P - mu) / sd >= z_fa)),
                                     tpl_z_med=(float(np.nanmedian(T)) - mu) / sd))
        K.log(f"[sweep] {date} {hh:02d}h -> {len(rows)} rows")
    R = pd.DataFrame(rows)
    R.to_csv(f"{a.out}.csv", index=False)
    if len(R):
        g = R.groupby(["comp", "win"]).agg(
            n=("bkey", "size"), pair=("pair_med", "median"), pair_z=("pair_z_med", "median"),
            pair_z90=("pair_z_p90", "median"), pair_zmax=("pair_z_max", "median"),
            tpl=("tpl_med", "median"), tpl_z=("tpl_z_med", "median"),
            z_fa=("z_fa", "median"), frac_fa=("frac_pair_ge_fa", "mean"))
        K.log("\n" + str(g.round(3)))


if __name__ == "__main__":
    main()
