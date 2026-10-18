#!/usr/bin/env python
"""F1184 Gate 1 at the withheld / never-summed sites -- DESCRIPTIVE, decision-irrelevant.

F1184 is the second template-free survivor of P4-0b arm (ii): 24 independent hits on 18 days across
3 episodes, +7.7..+19.1 dB alignment-free S excess at 9 of 10 DETECTION stations, N-moveout ratio
4.0, and NO match to any eligible kept family (best 1 match, null max 1).  Detection-site excess is
selection-biased (Merlin F3 of the P4 review), so the only non-circular read is at sites that were
never in any network sum.

  stack    alignment-free mean of the 24 independent hits at the bank's region-centre moveout,
           carried to CLRS, GOBB, LOPEZ, TURTL, LRIV, NLLB, B926
  statistic s_common.s_excess_db, maximised over a +-1.5 s constant offset
  null     200 same-day random-time draws, the SAME offset search charged to each
  PASS     >= 2 stations above their own null p99

Also reported: tremor co-timing of the hits against a same-day random base rate, so that a PASS
reads "coherent repeating signal of unknown class", not "LFE family".
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

FRZ = V.OUT
OUT = f"{V.RES}/p4/p4_0b_rescore"
WITHHELD = ["CLRS", "GOBB", "LOPEZ", "TURTL", "LRIV", "NLLB", "B926"]
HALF_W = 20.0
NW = int(round(2 * HALF_W * V.FS))
IDX_S = int(round(HALF_W * V.FS))
LAG = int(round(1.5 * V.FS))
N_NULL = 200
TID = "F1184"


def cut(dc, sta, t_abs):
    ts = pd.to_datetime(float(t_abs), unit="s")
    v = dc.get(sta, int(ts.year), int(ts.dayofyear))
    if v is None:
        return None
    t0, X = v
    i0 = int(round((t_abs - t0) * V.FS))
    if i0 < 0 or i0 + NW > X.shape[1]:
        return None
    seg = X[:, i0:i0 + NW].astype(np.float64)
    r = np.sqrt((seg ** 2).mean())
    if not np.isfinite(r) or r <= 0:
        return None
    return seg / r


def excess_with_offset(st):
    return max(S.s_excess_db(st, IDX_S + k) for k in range(-LAG, LAG + 1, 2))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(20260920)
    H = pd.read_parquet(f"{FRZ}/run_free_it2_indep.parquet")
    g = H[H.tid == TID].sort_values("t")
    times = g.t.to_numpy(float)
    lat, lon, dep, _ = V.region_centre()
    tt_reg_det = V.cell_tt(lat, lon, dep)
    tids, stas, W, off, meta = V.load_bank(f"{FRZ}/bank_free_it2.npz")
    bi = list(tids).index(TID)
    dts = [off[bi, si] + V.TPL1_PRE_S - tt_reg_det[s]
           for si, s in enumerate(stas) if np.isfinite(off[bi, si])]
    dt_cum = float(np.median(dts))
    tt_w = V.cell_tt(lat, lon, dep, stas=WITHHELD)
    day0 = np.array([pd.Timestamp(d).timestamp() for d in g.date])
    rows = []
    # One pass per (station, day): the day is loaded ONCE and every true and null cut for that day
    # is taken from it.  (A per-cut DayCache thrashes: 200 nulls x 18 days x 7 stations = 25k loads.)
    for w in WITHHELD:
        shift = tt_w[w] + dt_cum
        true_segs = []
        null_segs = [[] for _ in range(N_NULL)]
        for d in sorted(set(g.date)):
            m = (g.date.to_numpy() == d)
            ts = times[m]
            ddt = pd.Timestamp(d)
            d0 = ddt.timestamp()
            v = S.load_day_zne(w, int(ddt.year), int(ddt.dayofyear))
            if v is None:
                continue
            t0, X = v

            def take(t_abs, X=X, t0=t0):
                i0 = int(round((t_abs + shift - t0) * V.FS))
                if i0 < 0 or i0 + NW > X.shape[1]:
                    return None
                seg = X[:, i0:i0 + NW].astype(np.float64)
                r = np.sqrt((seg ** 2).mean())
                return None if (not np.isfinite(r) or r <= 0) else seg / r

            for t in ts:
                c = take(t)
                if c is not None:
                    true_segs.append(c)
            rt = d0 + rng.uniform(600.0, 86400.0 - 600.0, (N_NULL, len(ts)))
            for k in range(N_NULL):
                for t in rt[k]:
                    c = take(t)
                    if c is not None:
                        null_segs[k].append(c)
            del X
        if len(true_segs) < 5:
            rows.append(dict(sta=w, n_stacked=len(true_segs), excess_db=np.nan, null_n=0,
                             null_med=np.nan, null_p99=np.nan, above_p99=False))
            continue
        ex = float(excess_with_offset(np.stack(true_segs).mean(0)))
        nulls = np.asarray([excess_with_offset(np.stack(vv).mean(0))
                            for vv in null_segs if len(vv) >= 5])
        p99 = float(np.quantile(nulls, 0.99)) if len(nulls) else np.nan
        rows.append(dict(sta=w, n_stacked=len(true_segs), excess_db=ex, null_n=int(len(nulls)),
                         null_med=float(np.median(nulls)) if len(nulls) else np.nan,
                         null_p99=p99, above_p99=bool(np.isfinite(p99) and ex > p99)))
        V.log(f"  [{w}] n={len(true_segs)} excess {ex:+.2f} dB  null med "
              f"{rows[-1]['null_med']:+.2f} p99 {p99:+.2f} -> "
              f"{'ABOVE' if rows[-1]['above_p99'] else 'below'}")
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/F1184_gate1.csv", index=False)

    # tremor co-timing vs same-day random base rate
    tr = pd.read_csv(V.PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"])
    tr = tr[S.epi_km(tr.lat.values, tr.lon.values, lat, lon) <= 60.0]
    tt = np.sort(tr.time.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9)

    def frac_near(ts, win=1800.0):
        j = np.searchsorted(tt, ts)
        out = []
        for k, t in zip(j, ts):
            lo, hi = np.searchsorted(tt, t - win), np.searchsorted(tt, t + win)
            out.append(hi - lo > 0)
        return float(np.mean(out))

    obs = frac_near(times)
    base = [frac_near(day0 + rng.uniform(0, 86400, len(times))) for _ in range(200)]
    npass = int(R.above_p99.sum())
    summ = dict(label="DESCRIPTIVE -- decision-irrelevant; F1184 matched no eligible kept family",
                tid=TID, n_indep_hits=int(len(times)), n_days=int(g.date.nunique()),
                dt_cum_s=dt_cum, withheld_sites=WITHHELD,
                n_above_null_p99=npass, gate1_pass=bool(npass >= 2),
                per_station={r["sta"]: dict(n=r["n_stacked"], excess_db=r["excess_db"],
                                            null_p99=r["null_p99"], above=r["above_p99"])
                             for r in rows},
                tremor_cotiming=dict(observed=obs, base_rate_mean=float(np.mean(base)),
                                     base_rate_p05=float(np.quantile(base, 0.05)),
                                     base_rate_p95=float(np.quantile(base, 0.95)),
                                     window_s=1800.0, radius_km=60.0),
                reading="a PASS here means a coherent repeating signal of unknown class registering "
                        "at sites that were never in the detection sum; it is NOT evidence of an "
                        "LFE family, and the 30 hidden families were never matched against it")
    json.dump(summ, open(f"{OUT}/F1184_gate1.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
