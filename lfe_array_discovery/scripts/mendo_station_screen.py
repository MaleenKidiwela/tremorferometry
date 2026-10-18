#!/usr/bin/env python
"""MENDO (2026-09-29) single-station SELF-REPEAT screen for the autocorrelation station preference.

Why (smoke test 2024-07-29 00-05 h): network autocorrelation gave ~200 detections / hour (NW: ~1-5),
also in hours where Shelly's catalogue has none.  At the 60 strongest pairs PB.B933 alone carried
per-channel cc 0.67-0.73 (median 0.80) against <= 0.15 on every other station: a station whose own
record repeats (periodic noise, 1.9-2.7 Hz) makes the network autocorrelation fire on its own.  The
burst mute cannot see it (it is continuous, not a burst).  Same failure mode as W060 in NW
(memory noisy-station-hijacks-growth).

STATISTIC (per 3-C S1 station, per sampled hour, exactly the ab_auto window grid: 15 s windows, 0.5 s
hop, lags >= 15 s, mask-valid windows only): the station's own network function = mean over its live
channels of the window-pair cc; rate = fraction of pairs where it exceeds 0.6.  Station value = median
over the sampled hours.  FLAG: rate > max(10 x the median rate over the stations, 1e-5).  (Smoke test,
2 days x 8 hours: B933 2.9e-4, every other station 0; at cc > 0.3 / 0.5 the 1 Hz band-edge energy of
most stations dominates and nothing separates.)
A flagged station is NOT discarded: it moves to the END of the autocorrelation preference list (after
BRIC / RVIT), stays in the growth / location pool, and the table is kept as a column.

  python scripts/mendo_station_screen.py [--every 8] [--hours 0,6,12,18]
-> results/ab_royer/tremor/auto_screen_MENDO.csv  (sta, n_hours, rate_med, rate_p90, ratio, flag)
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_common as A  # noqa: E402

THR = 0.6
FLAG_X = 10.0
FLAG_FLOOR = 1e-5     # > ~250 self-repeating window pairs at cc > 0.6 per hour, as the median hour


def self_rate(X, m, dev):
    import torch
    n = int(round(A.WIN_S * A.FS)); hop = int(round(A.HOP_S * A.FS))
    ew, vw = A.window_stats(X, m, n, hop)
    Xt = torch.from_numpy(np.ascontiguousarray(X, np.float32)).to(dev)
    NR = (X.shape[1] - n) // hop + 1
    w = Xt.unfold(-1, n, hop)[:, :NR]
    nr = w.pow(2).sum(-1).sqrt()
    acc = torch.zeros((NR, NR), device=dev); k = 0
    for c in range(3):
        ok = torch.from_numpy(vw).to(dev) & (nr[c] > 0)
        if int(ok.sum()) < 100:
            continue
        u = torch.where(ok[:, None], w[c] / nr[c].clamp_min(1e-30)[:, None], torch.zeros_like(w[c]))
        acc += u @ u.T; k += 1
    if k == 0:
        return np.nan
    acc /= k
    i = torch.arange(NR, device=dev)
    far = (i[None, :] - i[:, None]).abs() * hop >= n
    okw = torch.from_numpy(vw).to(dev)
    sel = far & okw[:, None] & okw[None, :]
    if int(sel.sum()) == 0:
        return np.nan
    return float(((acc > THR) & sel).sum() / sel.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=8)
    ap.add_argument("--hours", default="0,6,12,18")
    ap.add_argument("--out", default=f"{A.RES}/tremor/auto_screen_MENDO.csv")
    a = ap.parse_args()
    import torch
    dev = A.torch_dev()
    days = [x for x, _ in A.SUBARRAY["MENDO"]["episodes"]][::a.every]
    hours = [int(h) for h in a.hours.split(",")]
    stas = A.SUBARRAY["MENDO"]["stas"]
    rows = []
    for d in days:
        dc = A.DayCache(stas)
        for hh in hours:
            H = dc.hour(d, hh)
            for s, (X, m) in H.items():
                if m.mean() < 0.5 or not (X[1:] ** 2).sum() > 0:
                    continue
                rows.append(dict(date=d, hour=hh, sta=s, rate=self_rate(X, m, dev)))
        A.log(f"[screen] {d} done ({len(rows)} station-hours)")
    R = pd.DataFrame(rows)
    R.to_csv(a.out.replace(".csv", "_hours.csv"), index=False)
    G = R.groupby("sta").rate.agg(n_hours="count", rate_med="median",
                                  rate_p90=lambda v: float(np.nanquantile(v, 0.9))).reset_index()
    ref = float(G.rate_med.median())
    G["ratio"] = G.rate_med / max(ref, 1e-12)
    G["flag"] = G.rate_med > max(FLAG_X * ref, FLAG_FLOOR)
    G.sort_values("ratio", ascending=False).to_csv(a.out, index=False)
    A.log("\n" + G.sort_values("ratio", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
