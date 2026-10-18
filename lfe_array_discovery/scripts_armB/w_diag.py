#!/usr/bin/env python
"""Stage 9 -- ENGINE SNR-MARGIN DIAGNOSTICS (why Step 1 came out the way it did).

Three measurements, all on KEPT (non-hidden) Bostock families only -- the hidden families are
never touched here:

  D1  DYNAMIC RANGE.  The per-station z of a PERFECT repeat (a window correlated with itself) is
      the ceiling of the whole autocorrelation statistic: z_self = 1 / (1.4826 x (median - Q25) of
      the hour's cc distribution).  If z_self is ~14, a real pair needs cc >= 0.35 to reach z = 5.
  D2  TRUE-PAIR z.  For the largest kept families, the per-station z at the TRUE lag of pairs of
      catalogued events in the same hour (moveout-aware cut at the family's own predicted S, i.e.
      the most favourable cut that exists), against a same-width jitter null.  This is the
      single-pair SNR the autocorrelation stage has to work with.
  D3  PAIR INVARIANT.  Of the pairs the engine actually reports, what fraction has a lag equal to
      an inter-event interval of one kept family at a common offset -- against a lag-jitter null.

Together these say whether a null Step-1 result is an ENGINE defect or the bed's physics.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import s_common as S                                                        # noqa: E402


def _prep(x, n, dev="cuda"):
    import torch
    NH = x.shape[1]
    NV = NH - n + 1
    X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
    HT = [X[c].unfold(0, n, 1)[:NV].contiguous() for c in range(3)]
    xd = X.double()
    cs = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd],
                                1), 1)
    ew = (cs[:, n:NH + 1] - cs[:, :NV])[:, :NV]
    E = ew.clamp_min(1e-10 * ew.median()).sqrt().float()
    return HT, E, NV


def _z(HT, E, NV, seg, dev="cuda"):
    import torch
    A = torch.from_numpy(np.ascontiguousarray(seg)).to(dev)
    r = A.pow(2).sum(-1, keepdim=True).sqrt()
    A = torch.where(r > 0, A / r.clamp_min(1e-30), torch.zeros_like(A))
    cc = torch.zeros(NV, device=dev)
    for c in range(3):
        cc += (A[c][None] @ HT[c].T)[0] / E[c]
    cc /= 3.0
    m = cc.median()
    q = torch.kthvalue(cc, max(1, NV // 4)).values
    sc = (1.4826 * (m - q)).clamp_min(1e-9)
    return ((cc - m) / sc).cpu().numpy(), float(1.0 / sc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--pairs", default="")
    ap.add_argument("--n-families", type=int, default=5)
    ap.add_argument("--n-pairs", type=int, default=25)
    ap.add_argument("--year", type=int, default=0)
    a = ap.parse_args()
    year = a.year or W.STEP1_YEAR
    sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
    import p_common as P
    mem, of, order, det = W.sites()
    # FIX 2026-09-20 (agent G).  `S.fam_tt(det)` defaults to the lfe_picker P3D roster station
    # table, which does not contain the Step-1b/1c bed stations (BPCB, GLBC, SHVB, TWKB, ...), so
    # the travel-time build raised KeyError and killed the diagnostic stage at the very end of both
    # the 1b and the 1c queue.  `w_common.station_table()` is the P3D roster EXTENDED with those
    # stations' FDSN coordinates -- the same table every other engine stage already uses.
    FT = S.fam_tt(det, ST=W.station_table())
    _ST = W.station_table()
    _FAM = P.families()
    hid, iso, _ = S.holdouts()

    def _dist(bk, sta):
        return float(S.epi_km(_FAM.lat[bk], _FAM.lon[bk], _ST.lat[sta], _ST.lon[sta]))

    B = P.bostock()
    B = B[(B.year == year) & (~B.bkey.isin(hid))]
    # FIX 2026-09-20 (agent G).  Restrict to the days the BED actually holds: without this the
    # "busiest hour" of a family can fall outside the bed window, every station load returns None
    # and the family is silently skipped (on the 2005 bed that emptied D1/D2 entirely).
    _beddays = {d.strftime("%Y-%m-%d") for d in W.bed_days(year)}
    B = B[B.time.dt.strftime("%Y-%m-%d").isin(_beddays)]
    top = B.groupby("bkey").size().sort_values(ascending=False).head(a.n_families).index.tolist()
    loader = W.DayLoader(maxn=len(det) + 4)
    n = int(round(W.AC_WIN_S * W.FS))
    rng = np.random.default_rng(20260920)
    d1, d2 = [], []
    for bk in top:
        bk = int(bk)
        g = B[B.bkey == bk].sort_values("t_origin")
        t = g.t_origin.to_numpy(float)
        tt = {s: FT[(bk, s)] for s in det}
        hh = (t // 3600).astype(int)
        byh = collections.defaultdict(list)
        for i in range(len(t) - 1):
            for j in range(i + 1, min(i + 10, len(t))):
                if hh[i] == hh[j] and 15 <= t[j] - t[i] <= 1200:
                    byh[hh[i]].append((t[i], t[j]))
        if not byh:
            continue
        k = max(byh, key=lambda z: len(byh[z]))
        ts = pd.Timestamp(k * 3600, unit="s")
        date, hour = ts.strftime("%Y-%m-%d"), int(ts.hour)
        D = {s: loader.hour(s, date, hour) for s in det}
        D = {s: v for s, v in D.items() if v is not None}
        if len(D) < 4:
            continue
        t0 = list(D.values())[0][0]
        NH = min(v[1].shape[1] for v in D.values())
        prep = {s: _prep(D[s][1][:, :NH], n) for s in D}
        # FIX 2026-09-20 (agent G).  D1 used the family's MEDIAN event over the whole bed, which
        # generally lies in a different hour from the loaded one -> every cut fell outside the hour
        # and D1 came out empty on any bed with more than a handful of events per family.  Use an
        # event that is actually inside the loaded hour.
        t_d1 = byh[k][0][0]
        for s in D:
            HT, E, NV = prep[s]
            ia = int(round((t_d1 + tt[s] - 1.0 - t0) * W.FS))
            if ia < 0 or ia + n > NH:
                continue
            z, inv = _z(HT, E, NV, D[s][1][:, ia:ia + n])
            d1.append(dict(bkey=bk, sta=s, site=of[s], dist_km=_dist(bk, s),
                           z_self=float(z[ia]), z_self_ceiling=float(inv),
                           hour_max_z=float(z.max())))
        pl = byh[k][: a.n_pairs]
        for (ta, tb) in pl:
            for s in D:
                HT, E, NV = prep[s]
                ia = int(round((ta + tt[s] - 1.0 - t0) * W.FS))
                jb = int(round((tb + tt[s] - 1.0 - t0) * W.FS))
                if ia < 0 or ia + n > NH or jb < 0 or jb >= NV:
                    continue
                z, _ = _z(HT, E, NV, D[s][1][:, ia:ia + n])
                w = 20
                zt = float(z[max(jb - w, 0):jb + w + 1].max())
                zn = []
                for _ in range(10):
                    jr = int(rng.integers(w + 1, NV - w - 1))
                    zn.append(float(z[jr - w:jr + w + 1].max()))
                d2.append(dict(bkey=bk, sta=s, site=of[s], dist_km=_dist(bk, s),
                               dt=float(tb - ta), z_true=zt,
                               z_null_median=float(np.median(zn)), z_null_max=float(max(zn))))
        for s in list(prep):
            prep[s] = None
        import torch
        torch.cuda.empty_cache()
    D1 = pd.DataFrame(d1)
    D2 = pd.DataFrame(d2)
    D1.to_csv(a.out + "_D1_dynamic_range.csv", index=False)
    D2.to_csv(a.out + "_D2_true_pair_z.csv", index=False)
    summ = dict(
        year=int(year), bed_mode=W.BED_MODE, sites=order,
        n_families=len(top), families=[int(x) for x in top],
        D1_median_z_self=float(D1.z_self.median()) if len(D1) else None,
        D1_median_hour_max_z=float(D1.hour_max_z.median()) if len(D1) else None,
        D2_median_z_true=float(D2.z_true.median()) if len(D2) else None,
        D2_median_z_null=float(D2.z_null_median.median()) if len(D2) else None,
        D2_frac_true_gt_null_max=float((D2.z_true > D2.z_null_max).mean()) if len(D2) else None,
    )
    # ---- D3 pair invariant on the engine's own output
    if a.pairs and os.path.exists(a.pairs):
        Pp = pd.read_parquet(a.pairs)
        if len(Pp):
            Pp = Pp.sort_values("net_mad", ascending=False).head(2000)
            ev = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}

            def match(ti, tj, tol=0.5):
                hits = 0
                for aa, bb in zip(ti, tj):
                    d = bb - aa
                    ok = False
                    for e in ev.values():
                        lo = np.searchsorted(e, aa - 32.0)
                        hi = np.searchsorted(e, aa + 2.0)
                        for x in e[lo:hi]:
                            y = x + d
                            kk = np.searchsorted(e, y - tol)
                            if kk < len(e) and abs(e[kk] - y) <= tol:
                                ok = True
                                break
                        if ok:
                            break
                    hits += ok
                return hits

            ti = Pp.t_i.to_numpy(float)
            tj = Pp.t_j.to_numpy(float)
            tr = match(ti, tj)
            nl = [match(ti, tj + rng.uniform(-60, 60, len(tj))) for _ in range(3)]
            summ["D3_n_pairs_tested"] = int(len(Pp))
            summ["D3_true_matches"] = int(tr)
            summ["D3_null_matches"] = [int(x) for x in nl]
    json.dump(summ, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
