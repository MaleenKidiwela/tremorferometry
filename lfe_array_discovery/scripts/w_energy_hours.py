#!/usr/bin/env python
"""Step 1c -- the ARRAY-ENERGY HOUR FLAG (hour selection only).

Merlin's Q3 ruling, verbatim: for a pre-2010 window with no tremor catalogue, "the array's own
hourly network-median 2-8 Hz energy relative to the window's quiet-day distribution -- a coarse
tremor-presence indicator, legitimate for HOUR selection because P4-0b's negative was about energy
as an EVENT SEED, not as an hour flag."

Statistic, per hour of the bed:
  per station   median over the hour of the 3-C summed squared amplitude of the FROZEN 1-8 Hz
                20 Hz cache (robust to transients), divided by that station's OWN median over the
                frozen quiet days -- so station gain and site noise cancel;
  per site      median over member stations;
  network       median over the observable sites (>= MIN_SITES required).
An hour is a CANDIDATE when the network ratio is at or above the ENERGY_HOUR_Q quantile of the
same statistic computed on the quiet days.

No Bostock event time and no family coordinate enters this file.
"""
from __future__ import annotations
import argparse, json, os, sys

import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402

ENERGY_HOUR_Q = 0.95
BAND_NOTE = "the cache is already Butterworth-4 zero-phase 1-8 Hz, so the 2-8 Hz of Merlin's rule " \
            "is the cache band itself; recorded as a documented deviation"


def hour_energy(loader, sta, date, hour):
    v = loader.hour(sta, date, hour)
    if v is None:
        return np.nan
    x = v[1].astype(np.float64)
    e = (x ** 2).sum(0)
    e = e[np.isfinite(e)]
    if e.size < 1000:
        return np.nan
    return float(np.median(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/step1c_energy_hours.csv")
    ap.add_argument("--q", type=float, default=ENERGY_HOUR_Q)
    a = ap.parse_args()
    mem, of, order, det = W.sites()
    loader = W.DayLoader(maxn=len(det) + 4)
    A = W.availability()
    Q = W.quiet_days()
    qdays = set(Q[Q.quiet].date)
    rows = []
    for r in A[A.usable].itertuples():
        for h in range(24):
            e = {s: hour_energy(loader, s, r.date, h) for s in det}
            sz = {}
            for q in order:
                v = [e[s] for s in mem[q] if np.isfinite(e.get(s, np.nan))]
                if v:
                    sz[q] = float(np.median(v))
            if len(sz) < W.MIN_SITES:
                continue
            rows.append(dict(date=r.date, hour=h, n_sites=len(sz),
                             **{f"e_{q}": sz.get(q, np.nan) for q in order}))
        W.log(f"[energy] {r.date} done ({len(rows)} hours)")
    E = pd.DataFrame(rows)
    if not len(E):
        raise SystemExit("no hours")
    cols = [f"e_{q}" for q in order]
    # per-station (here per-site) normalisation by its OWN quiet-day median
    qm = E[E.date.isin(qdays)][cols].median()
    R = E[cols] / qm
    E["net_ratio"] = R.median(axis=1, skipna=True)
    thr = float(np.nanquantile(E[E.date.isin(qdays)].net_ratio, a.q))
    E["candidate"] = (~E.date.isin(qdays)) & (E.net_ratio >= thr)
    E.to_csv(a.out, index=False)
    json.dump(dict(quantile=a.q, threshold=float(thr), n_hours=int(len(E)),
                   n_quiet_hours=int(E.date.isin(qdays).sum()),
                   n_candidate=int(E.candidate.sum()),
                   n_candidate_days=int(E[E.candidate].date.nunique()),
                   site_quiet_median={c: float(qm[c]) for c in cols},
                   band_note=BAND_NOTE,
                   rule="Merlin Q3: array energy is legitimate for HOUR selection only"),
              open(a.out.replace(".csv", ".json"), "w"), indent=1)
    W.log(f"[energy] threshold {thr:.3f} (q={a.q}); {int(E.candidate.sum())} candidate hours on "
          f"{int(E[E.candidate].date.nunique())} days of {len(E)} hours")


if __name__ == "__main__":
    main()
