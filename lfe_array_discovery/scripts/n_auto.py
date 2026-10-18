#!/usr/bin/env python
"""Puget/S-Olympic pair stage: agent K's LITERAL autocorrelation engine (`k_auto.autocorr_block`,
unchanged) driven over an EXPLICIT (date, hour) list with the station-hour health veto applied.
Prefix `n_`.  Byte-for-byte the 2010 run's `l_auto.py` with `l_common` -> `n_common`; the hour
list is the AUDITED tremor-hour set of `n_hours.py` instead of the PNSN candidate hours (there is
no tremor catalogue before 2010-01-06).

The only thing this adds to `k_auto.main` is what the WA bed needs and the 2005 bed did not have:
  * an explicit list of hours (PNSN tremor hours, or the 40 hour-matched quiet hours), rather than
    "all hours of these dates";
  * the per-station-hour health mask of results/wa_station_hour_veto_<year>.csv (WA_V2_INPUTS sec.3:
    "adopt the station-hour veto table as a per-hour mask, not a station-level exclusion -- the WA
    problem is hour-shaped").
Every recipe constant (6 s window, 0.5 s row step, |lag| >= 6 s, plain per-component mean over all
channels, median + 8 x MAD of the block's own distribution, 4 s decluster) comes from k_auto and is
not restated here.

  python scripts/n_auto.py --hours tremor --out results/n_puget2008/2008/auto
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
import n_common as L                                                        # noqa: E402

MIN_STA = 4                                                                 # k_auto's own bar


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", default="tremor", choices=["tremor", "quiet"])
    ap.add_argument("--out", default=f"{L.OUT}/auto")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="", help="date:hour, for the smoke test")
    a = ap.parse_args()
    K = L.activate()
    import k_auto as KA
    os.makedirs(a.out, exist_ok=True)
    stas = L.stations()
    H = L.candidate_hours() if a.hours == "tremor" else L.quiet_hours()
    HL = [(str(r.date), int(r.hour)) for r in H.itertuples()]
    if a.only:
        d, h = a.only.split(":")
        HL = [(d, int(h))]
    if a.limit:
        HL = HL[:a.limit]
    veto = L.health_veto()
    L.log(f"{len(stas)} stations, {len(HL)} {a.hours} hours -> {a.out}")
    dc = K.DayCache(stas)
    n_done = 0
    for (date, hh) in HL:
        p = f"{a.out}/{date}_{hh:02d}.parquet"
        if os.path.exists(p):
            continue
        bad = veto.get((date, hh), set())
        D0 = dc.day(date)
        t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
        Db = {}
        for s in stas:
            if s in bad or s not in D0:
                continue
            t0, X = D0[s]
            i0 = int(round((t_h - t0) * K.FS))
            seg = X[:, i0:i0 + int(round(3600 * K.FS))]
            if seg.shape[1] < int(3000 * K.FS) or not np.isfinite(seg).all():
                continue
            Db[s] = seg
        if len(Db) < MIN_STA:
            L.log(f"[skip] {date} {hh:02d}h {len(Db)} stations")
            continue
        try:
            R, st = KA.autocorr_block(Db, K.FS)
        except RuntimeError as e:                                           # GPU shared with K
            L.log(f"[retry] {date} {hh:02d}h {type(e).__name__}: {str(e)[:120]}")
            import torch
            import time
            torch.cuda.empty_cache()
            time.sleep(60)
            R, st = KA.autocorr_block(Db, K.FS)
        if R is None:
            continue
        hopn = int(round(KA.HOP_S * K.FS))
        R["t_i"] = t_h + R.i.to_numpy(np.float64) * hopn / K.FS
        R["t_j"] = t_h + R.j.to_numpy(np.float64) / K.FS
        R["date"], R["hour"] = date, hh
        R.to_parquet(p, index=False)
        st.update(date=date, hour=hh, n_veto=len(bad), stas=",".join(sorted(Db)))
        json.dump(st, open(f"{a.out}/{date}_{hh:02d}_stat.json", "w"))
        n_done += 1
        L.log(f"[auto] {date} {hh:02d}h  {len(Db)} sta (veto {len(bad)})  thr={st['thr']:.5f} "
              f"(={st['nmad_in_sd']:.2f} sd) -> {len(R)} detections")
    sf = sorted(glob.glob(f"{a.out}/*_stat.json"))
    if sf:
        S = pd.DataFrame([json.load(open(f)) for f in sf])
        S.to_csv(f"{a.out}/_stats.csv", index=False)
        L.log(f"\n{len(S)} blocks; n_det median {S.n_det.median():.0f}; "
              f"8xMAD = {S.nmad_in_sd.median():.2f} sd (median)")
        json.dump(dict(win_s=KA.WIN_S, hop_s=KA.HOP_S, nmad=KA.NMAD, min_lag_s=KA.MIN_LAG_S,
                       decluster_s=KA.DECL_S, stations=stas, n_stations=len(stas),
                       hours=a.hours, n_blocks=int(len(S)),
                       nmad_in_sd_median=float(S.nmad_in_sd.median())),
                  open(f"{a.out}/_params.json", "w"), indent=1)


if __name__ == "__main__":
    main()
