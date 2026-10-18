#!/usr/bin/env python
"""Verification WITHOUT catalogue truth for the WA 2010 grown families.  Prefix `l_`.

There is no Washington LFE catalogue, so the north's Bostock >= 30 % one-to-one bar cannot be
applied.  Replaced by, per grown family:

  (a) LOCATION -- `l_locate.py` (agent K's k_locate, WA grid + intraslab S station terms):
      >= 5 picked stations, shuffled-station control, |dz| to Slab2.   [read from --loc]
  (b) REPEAT STATISTICS -- detections, distinct days, and TREMOR CO-TIMING: the fraction of the
      family's detections within +- 30 min of a PNSN epicentre <= 60 km of the family location,
      against the SAME-DAY RANDOM base rate (the same detection count redrawn uniformly on the
      same days, `--n-base` realisations).  A family passes only if the observed fraction exceeds
      the base rate.
  (c) the identical screen applied to the whole-chain SITE-SHIFT NULL, counted beside the real
      number (run separately with --loc/--det pointing at the null chain).
  (d) DESCRIPTIVE only: distance from each candidate to the nearest Chestler-derived coordinate
      (`derived_coordinate_rows_not_a_catalog.csv` -- NOT a catalogue, no truth is claimed of it)
      against 1,000 random draws from the sub-region's active cells.

CANDIDATE FAMILY = located on >= 5 stations AND shuffled control passes AND |dz| <= 8 km AND
>= 3 distinct days AND tremor co-timing above the same-day base rate.
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
import l_common as L                                                        # noqa: E402

CO_MIN = 30.0 * 60.0        # +- 30 min
CO_KM = 60.0
DZ_MAX = 8.0
MIN_DAYS = 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loc", required=True)
    ap.add_argument("--det", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default="real")
    ap.add_argument("--n-base", type=int, default=200)
    a = ap.parse_args()
    L.activate()
    import s_common as S
    R = pd.read_csv(a.loc)
    D = pd.read_parquet(a.det)
    T = L.pnsn()
    tla, tlo, tt = T.lat.to_numpy(), T.lon.to_numpy(), T.t.to_numpy()
    THR = L.tremor_hour_set()
    rng = np.random.default_rng(17)

    # ---- fraction of ALL grown detections inside a PNSN tremor-active hour, vs the base rate
    key = pd.to_datetime(D.t, unit="s")
    in_tr = np.array([(d.strftime("%Y-%m-%d"), int(d.hour)) in THR for d in key])
    days_run = sorted(set(D.date))
    base_hours = sum(1 for d in days_run for h in range(24) if (d, h) in THR) / max(
        24 * len(days_run), 1)
    L.log(f"[{a.tag}] {len(D)} grown detections; {in_tr.mean():.3f} in a PNSN tremor-active hour "
          f"(hour base rate {base_hours:.3f})")

    rows = []
    for r in R.itertuples():
        d = D[D.grp == r.grp]
        rec = dict(grp=int(r.grp), n_det=len(d), n_days=int(d.date.nunique()),
                   located=bool(r.located), n_sta=int(r.n_sta))
        if len(d):
            k = pd.to_datetime(d.t, unit="s")
            rec["frac_tremor_hour"] = float(np.mean(
                [(x.strftime("%Y-%m-%d"), int(x.hour)) in THR for x in k]))
        if not bool(r.located) or not len(d):
            rec.update(co_timing=np.nan, co_base=np.nan, candidate=False)
            rows.append(rec)
            continue
        rec.update(lat=float(r.lat), lon=float(r.lon), dep=float(r.dep),
                   slab2=float(r.slab2), dz=float(r.dz), rms=float(r.rms),
                   ctrl_pass=bool(r.ctrl_pass))
        dkm = S.epi_km(float(r.lat), float(r.lon), tla, tlo)
        te = np.sort(tt[dkm <= CO_KM])
        rec["n_pnsn_near"] = int(len(te))
        td = np.sort(d.t.to_numpy())
        if len(te) == 0:
            rec.update(co_timing=0.0, co_base=0.0, co_p=1.0)
        else:
            def frac(x):
                i = np.searchsorted(te, x)
                dd = np.full(len(x), np.inf)
                for kk in (-1, 0):
                    j = np.clip(i + kk, 0, len(te) - 1)
                    dd = np.minimum(dd, np.abs(te[j] - x))
                return float((dd <= CO_MIN).mean())
            obs = frac(td)
            dsec = np.array([pd.Timestamp(x).timestamp() for x in d.date])
            bs = np.array([frac(dsec + rng.random(len(dsec)) * 86400.0)
                           for _ in range(a.n_base)])
            rec.update(co_timing=obs, co_base=float(bs.mean()),
                       co_p=float((bs >= obs).mean()))
        rec["candidate"] = bool(r.located and r.ctrl_pass and np.isfinite(r.dz)
                                and abs(r.dz) <= DZ_MAX and rec["n_days"] >= MIN_DAYS
                                and rec["co_timing"] > rec["co_base"])
        rows.append(rec)
    V = pd.DataFrame(rows)
    V.to_csv(a.out, index=False)
    cand = V[V.candidate]
    L.log(f"[{a.tag}] {len(V)} templates; located {int(V.located.sum())}; "
          f"CANDIDATE FAMILIES {len(cand)}")

    rep = dict(tag=a.tag, n_templates=int(len(V)), n_located=int(V.located.sum()),
               n_ctrl_pass=int(V.get("ctrl_pass", pd.Series(dtype=bool)).fillna(False).sum()),
               n_candidate=int(len(cand)), n_det_total=int(len(D)),
               frac_det_in_tremor_hour=float(in_tr.mean()),
               tremor_hour_base_rate=float(base_hours))

    # ---- (d) descriptive Chestler comparison
    if len(cand) and os.path.exists(L.CHESTLER):
        C = pd.read_csv(L.CHESTLER)
        cla, clo = C.latitude.to_numpy(), C.longitude.to_numpy()
        dmin = np.array([np.min(S.epi_km(float(x), float(y), cla, clo))
                         for x, y in zip(cand.lat, cand.lon)])
        V.loc[cand.index, "km_to_chestler"] = dmin
        AC = L.active_cells()
        rr = np.random.default_rng(5)
        draws = []
        for _ in range(1000):
            i = rr.integers(0, len(AC), len(cand))
            draws.append(np.median([np.min(S.epi_km(float(AC.lat.iloc[j]), float(AC.lon.iloc[j]),
                                                    cla, clo)) for j in i]))
        draws = np.array(draws)
        rep["chestler"] = dict(n_ref=int(len(C)), median_km=float(np.median(dmin)),
                               min_km=float(dmin.min()),
                               random_cell_median_km=float(np.median(draws)),
                               random_cell_p05_km=float(np.quantile(draws, .05)),
                               p_better_than_random=float((draws <= np.median(dmin)).mean()))
        V.to_csv(a.out, index=False)
        L.log(f"[{a.tag}] Chestler: candidates median {np.median(dmin):.1f} km, "
              f"random active cells {np.median(draws):.1f} km "
              f"(p={rep['chestler']['p_better_than_random']:.3f})")
    json.dump(rep, open(a.out.replace(".csv", "_summary.json"), "w"), indent=1, default=float)
    L.log(json.dumps(rep, indent=1, default=float))


if __name__ == "__main__":
    main()
