#!/usr/bin/env python
"""ALL-STATION FINAL TEMPLATES -- Royer & Bostock (2014) sec.3.2 ("waveforms for all channels
available"; Fig.3 shows stations 1-71 km): for every final family, stack its detections on EVERY
station of the literal NW 25-station pool (Fig.1D) that has data, not only the 7-10 growth stations.
No extra public stations (coordinator 2026-09-26: keep it literal).

Stacking = ab_grow.restack's arithmetic (unit maximum amplitude across 3 components, gap-guard dead
channels zeroed, >= 2 detections), same 44 s long window [t - 12 s, t + 32 s] about each detection
time t.

MOVEOUT.  A detection time t is the network alignment of the family template: every station's
window is cut from the SAME absolute time t - 12 s, so each station's moveout relative to the
template (its predicted S/P delay from the current location) is already embedded in the common
time base -- the "stack at predicted moveout" step is exactly the common-time stack (pass 0).  The
44 s window holds S at the farthest pool station (72 km; S-P and inter-station S delays <= ~12 s
for 30-45 km deep sources, and the growth stations' arrivals lie in [t, t + 15 s]).

REFINEMENT (non-growth stations only; growth stations keep the network alignment that defines
the family):  pass p = 1..--passes: for each detection and each non-growth station, cross-correlate
the detection's 3C segment with the LEAVE-ONE-OUT station stack of pass p-1 (the stack minus this
detection's own contribution, so a detection is never aligned onto itself) over a +-4 s window
about that stack's 3C envelope peak, lag search +- --maxlag s (0.5 s, the same lag bound ab_hypo
uses for its phase-correlation differential times); the detection is re-stacked at the best lag.
No cc threshold (no invented gate): every detection enters at its best lag.  Per-station SNR
(peak / median 3C envelope, as ab_hypo.make_picks) and the shift statistics are written as columns
for common vs refined.

  python scripts/ab_allsta.py --det .../grow/det_it<K>.parquet --tpl-sta .../grow/tpl_it<K>_stations.csv \
      --out results/ab_royer/NW/hypo_allsta [--fams 5 --dates-max 20  (smoke)]
-> <out>/stacks_common.npz   pass-0 stacks (all pool stations, common time base)
-> <out>/stacks_loc.npz      refined stacks (read by ab_hypo.py as its stacks_loc.npz)
-> <out>/allsta_stack_stats.csv  per (family, station): growth?, n, snr_common, snr_refined,
                                 shift median/abs-median/std (s)
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.signal import hilbert

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_grow as G                                                         # noqa: E402
import ab_grow_cont as GC                                                   # noqa: E402

WIN_HALF = 4.0


def env3(st):
    return np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))


def snr_of(st):
    e = env3(st)
    m = np.median(e)
    return float(e.max() / m) if m > 0 else np.nan


def _refine_chunk(args):
    (stas, det, SUM, CNT, WIN, grow_sta, L, long_pre, long_len) = args
    dc = A.DayCache(stas)
    nL = int(round(long_len * A.FS))
    pad = L / A.FS
    rc = A.RefCache(dc, int(round(G.TPL_LEN * A.FS)))
    acc, cnt, shifts = {}, {}, []
    for date, dd in det.sort_values("t").groupby("date"):
        refs = {s: rc.get(s, date) for s in stas}
        for r in dd.itertuples():
            g = int(r.grp)
            for s in stas:
                if s in grow_sta.get(g, ()) or (g, s) not in SUM:
                    continue
                seg = dc.cut(s, r.t - long_pre - pad, long_len + 2 * pad, date)
                if seg is None:
                    continue
                okc = refs[s] > 0
                if not okc.any():
                    continue
                if not okc.all():
                    seg = np.where(okc[:, None], seg, 0.0)
                seg = seg[:, :nL + 2 * L]
                # LOO reference: the pass-(p-1) stack minus this detection's own contribution
                # (its unit_max3 window at its pass-(p-1) lag)
                l_prev = int(getattr(r, f"L_{s}"))
                own = G.unit_max3(seg[:, L + l_prev:L + l_prev + nL])
                n = CNT[(g, s)]
                ref = (SUM[(g, s)] - own) / max(n - 1, 1)
                i0, i1 = WIN[(g, s)]
                aw = ref[:, i0:i1]
                na = np.sqrt((aw ** 2).sum())
                bl, best = 0, np.nan
                if na > 0:
                    B = np.lib.stride_tricks.sliding_window_view(seg[:, i0:i1 + 2 * L], i1 - i0,
                                                                 axis=1)       # (3, 2L+1, W)
                    num = np.einsum("cw,clw->l", aw, B)
                    nb = np.sqrt(np.einsum("clw,clw->l", B, B))
                    c = np.where(nb > 0, num / (na * np.where(nb > 0, nb, 1.0)), -np.inf)
                    j = int(np.argmax(c))
                    bl, best = j - L, float(c[j])
                new = G.unit_max3(seg[:, L + bl:L + bl + nL])
                key = (g, s)
                if key not in acc:
                    acc[key] = np.zeros((3, nL))
                    cnt[key] = 0
                acc[key] += new
                cnt[key] += 1
                shifts.append((int(r.Index), s, bl, best))
    return acc, cnt, shifts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default="NW")
    ap.add_argument("--det", required=True)
    ap.add_argument("--tpl-sta", required=True, help="tpl_it<K>_stations.csv (growth stations)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--passes", type=int, default=2)
    ap.add_argument("--maxlag", type=float, default=0.5)
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--fams", type=int, default=0, help="smoke: first N families")
    ap.add_argument("--dates-max", type=int, default=0, help="smoke: first N dates")
    ap.add_argument("--null-shift", type=float, default=0.0,
                    help="NULL control: move every detection time by this many s (same day)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stas = list(A.POOL[a.sub])
    det = pd.read_parquet(a.det)
    if a.fams:
        keep = det.groupby("grp").size().sort_values(ascending=False).index[:a.fams]
        det = det[det.grp.isin(keep)]
    if a.dates_max:
        det = det[det.date.isin(sorted(det.date.unique())[:a.dates_max])]
    det = det.reset_index(drop=True)
    if a.null_shift:
        det["t"] = det.t + a.null_shift
    TS = pd.read_csv(a.tpl_sta)
    grow_sta = {int(r.cid): set(str(r.stations).split(",")) for r in TS.itertuples()
                if isinstance(r.stations, str)}
    A.log(f"[allsta] {det.grp.nunique()} families, {len(det)} detections, {len(stas)} pool stations")
    t0 = time.time()
    ST0, CNT0 = GC.restack_parallel(stas, det, workers=a.workers, return_counts=True)
    np.savez_compressed(f"{a.out}/stacks_common.npz", **{f"{g}_{s}": ST0[g][s].astype(np.float32)
                                                          for g in ST0 for s in ST0[g]})
    A.log(f"[allsta] pass 0 common-time stacks: {sum(len(v) for v in ST0.values())} "
          f"(family, station) ({time.time() - t0:.0f} s)")
    L = int(round(a.maxlag * A.FS))
    nL = int(round(G.LONG_LEN * A.FS))
    hw = int(round(WIN_HALF * A.FS))
    ST = ST0
    lagcol = {s: np.zeros(len(det), np.int16) for s in stas}
    SH = None
    PREVCNT = dict(CNT0)
    for p in range(1, a.passes + 1):
        t0 = time.time()
        # SUM/CNT of the previous pass's stack (exact counts) for the leave-one-out reference
        n_fam = det.groupby("grp").size().to_dict()
        SUM, CNT, WIN = {}, {}, {}
        for g, T in ST.items():
            for s, st in T.items():
                if s in grow_sta.get(g, ()):
                    continue
                n = PREVCNT.get((g, s), n_fam[g])
                SUM[(g, s)] = (st * n).astype(np.float64)
                CNT[(g, s)] = n
                k = int(np.argmax(env3(st)))
                i0, i1 = max(k - hw, 0), min(k + hw, nL)
                WIN[(g, s)] = (i0, i1)
        dd = det.copy()
        for s in stas:
            dd[f"L_{s}"] = lagcol[s]
        dates = sorted(dd.date.unique())
        chunks = [dates[i::a.workers] for i in range(a.workers)]
        jobs = [(stas, dd[dd.date.isin(c)], SUM, CNT, WIN, grow_sta, L, G.LONG_PRE, G.LONG_LEN)
                for c in chunks if c]
        acc, cnt, shifts = {}, {}, []
        with ProcessPoolExecutor(len(jobs)) as ex:
            for aa, cc, sh in ex.map(_refine_chunk, jobs):
                shifts += sh
                for k, v in aa.items():
                    if k in acc:
                        acc[k] += v
                        cnt[k] += cc[k]
                    else:
                        acc[k], cnt[k] = v, cc[k]
        PREVCNT = dict(cnt)
        SH = pd.DataFrame(shifts, columns=["idx", "sta", "lag", "cc"])
        for s, gsh in SH.groupby("sta"):
            lagcol[s][gsh.idx.to_numpy()] = gsh.lag.to_numpy()
        new = {}
        for g, T in ST0.items():
            for s, st in T.items():
                if s in grow_sta.get(g, ()):
                    new.setdefault(g, {})[s] = st           # growth station: network alignment
                elif cnt.get((g, s), 0) >= 2:
                    new.setdefault(g, {})[s] = acc[(g, s)] / cnt[(g, s)]
        ST = new
        A.log(f"[allsta] refine pass {p}: {len(SH)} (detection, station) shifts, median |lag| "
              f"{SH.lag.abs().median() / A.FS:.3f} s, median cc {SH.cc.median():.3f} "
              f"({time.time() - t0:.0f} s)")
    np.savez_compressed(f"{a.out}/stacks_loc.npz", **{f"{g}_{s}": ST[g][s].astype(np.float32)
                                                       for g in ST for s in ST[g]})
    rows = []
    SH = SH if SH is not None else pd.DataFrame(columns=["idx", "sta", "lag", "cc"])
    SH["grp"] = det.grp.to_numpy()[SH.idx.to_numpy()] if len(SH) else []
    agg = (SH.groupby(["grp", "sta"]).agg(n_ref=("lag", "size"),
                                          lag_med_s=("lag", lambda v: np.median(v) / A.FS),
                                          lag_absmed_s=("lag", lambda v: np.median(np.abs(v)) / A.FS),
                                          lag_std_s=("lag", lambda v: np.std(v) / A.FS),
                                          cc_med=("cc", "median"))
           if len(SH) else pd.DataFrame())
    for g, T in ST0.items():
        for s, st in T.items():
            r = dict(grp=g, sta=s, growth=s in grow_sta.get(g, ()), snr_common=snr_of(st),
                     snr_refined=snr_of(ST[g][s]) if s in ST.get(g, {}) else np.nan)
            if len(agg) and (g, s) in agg.index:
                r.update(agg.loc[(g, s)].to_dict())
            rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(f"{a.out}/allsta_stack_stats.csv", index=False)
    ng = R[~R.growth]
    A.log(f"[allsta] {R.grp.nunique()} families; stations/family median "
          f"{R.groupby('grp').size().median():.0f}; non-growth SNR common median "
          f"{ng.snr_common.median():.2f} -> refined {ng.snr_refined.median():.2f}")


if __name__ == "__main__":
    main()
