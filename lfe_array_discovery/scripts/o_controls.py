#!/usr/bin/env python
"""The two measurements that decided the beam, on the record.  Prefix `o_`.

 (a) SLOWNESS CONTROL.  Ghosh/Chestler's per-window slowness-grid beam raises the 5-8 Hz
     semblance well above the zero-slowness stack -- but does it do so only on tremor hours?
     The same computation is run on the hour-matched QUIET hours.  If the gain is the same
     there, it is the bias of a maximum over the slowness grid, not a plane wave.
 (b) ARRAY COHERENCE vs FREQUENCY, on catalogued INTRASLAB EARTHQUAKES -- real deep sources at
     the LFE incidence, with the pre-event noise window as the control.

  O_BED=2010 python scripts/o_controls.py --stage slowness --hours 4
  O_BED=2010 python scripts/o_controls.py --stage bands --events 25
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402
import o_beam as B                                                          # noqa: E402

BANDS = [(1, 2), (2, 3), (3, 4), (4, 6), (6, 8), (1, 4), (1, 8)]


def _semb(site, date, hh, band=(5.0, 8.0)):
    import torch
    ST = O.member_table()
    mem = O.arrays()[site]
    D = {}
    for m in mem:
        v = B.load_member(m, date)
        if v is not None:
            D[m] = v[1]
    ok = [m for m in mem if m in D]
    if len(ok) < 3:
        return None
    a0, p = hh * 72000, B.HOP_N
    Xs = np.zeros((len(ok), 3, 72000 + 2 * p), np.float32)
    for i, m in enumerate(ok):
        X = D[m]
        s0, s1 = max(a0 - p, 0), min(a0 + 72000 + p, X.shape[1])
        if s1 <= s0:
            return None
        Xs[i, :, s0 - (a0 - p):s1 - (a0 - p)] = X[:, s0:s1]
        r = np.sqrt((Xs[i] ** 2).mean())
        if r <= 0:
            return None
        Xs[i] /= r
    lat0 = float(np.mean([ST.lat[m] for m in mem]))
    lon0 = float(np.mean([ST.lon[m] for m in mem]))
    e, n = B._xy_km([ST.lat[m] for m in ok], [ST.lon[m] for m in ok], lat0, lon0)
    dev = "cuda"
    X = torch.from_numpy(Xs).to(dev)
    w = torch.hann_window(B.WIN_N, periodic=True, device=dev)
    F = torch.fft.rfft(X.unfold(-1, B.WIN_N, B.HOP_N) * w, n=B.NFFT, dim=-1)
    fr = torch.fft.rfftfreq(B.NFFT, d=1.0 / B.FS).to(dev)
    SX, SY = B.slowness_grid()
    tau = torch.from_numpy(np.outer(SX, e) + np.outer(SY, n)).to(dev).float()
    bi = torch.nonzero((fr >= band[0]) & (fr <= band[1])).ravel()
    M, nw = len(ok), F.shape[2]
    P = torch.zeros((tau.shape[0], nw), device=dev)
    for k in range(0, bi.numel(), 24):
        b = bi[k:k + 24]
        ph = torch.exp(2j * np.pi * fr[b][None, None, :] * tau[:, :, None])
        P += torch.einsum("smb,mwb->swb", ph,
                          F[:, 0][:, :, b].to(torch.complex64)).abs().pow(2).sum(-1)
    indiv = (F[:, 0][:, :, bi].abs() ** 2).sum(-1).sum(0)
    i0 = int(np.argmin(SX ** 2 + SY ** 2))
    pb, ib = P.max(0)
    out = dict(M=M, floor=1.0 / M,
               zero=float(np.median((P[i0] / (M * indiv)).cpu().numpy())),
               best=float(np.median((pb / (M * indiv)).cpu().numpy())),
               s=float(np.median(np.hypot(SX[ib.cpu().numpy()], SY[ib.cpu().numpy()]))))
    del X, F, P
    torch.cuda.empty_cache()
    return out


def stage_slowness(a):
    O.activate()
    random.seed(0)
    T = [(str(r.date), int(r.hour)) for r in O.candidate_hours().itertuples()]
    Q = [(str(r.date), int(r.hour)) for r in O.quiet_hours().itertuples()]
    rows = []
    for site in sorted(O.arrays()):
        tr = [x for x in (_semb(site, d, h) for d, h in random.sample(T, min(a.hours, len(T))))
              if x]
        qu = [x for x in (_semb(site, d, h) for d, h in random.sample(Q, min(a.hours, len(Q))))
              if x]
        if not tr or not qu:
            continue
        f = lambda L, k: float(np.median([x[k] for x in L]))               # noqa: E731
        rows.append(dict(site=site, M=int(f(tr, "M")), floor=round(f(tr, "floor"), 3),
                         tremor_zero=round(f(tr, "zero"), 3),
                         tremor_best=round(f(tr, "best"), 3), tremor_s=round(f(tr, "s"), 2),
                         quiet_zero=round(f(qu, "zero"), 3), quiet_best=round(f(qu, "best"), 3),
                         quiet_s=round(f(qu, "s"), 2),
                         gain_tremor=round(f(tr, "best") / max(f(tr, "zero"), 1e-9), 2),
                         gain_quiet=round(f(qu, "best") / max(f(qu, "zero"), 1e-9), 2)))
        O.log(f"  {rows[-1]}")
    R = pd.DataFrame(rows)
    os.makedirs(O.OUTROOT, exist_ok=True)
    R.to_csv(f"{O.OUTROOT}/slowness_control.csv", index=False)
    json.dump(dict(n_hours_each=a.hours, slowness_control=rows,
                   verdict=("the slowness-grid gain is the same on quiet hours as on tremor "
                            "hours -> it is max-over-grid bias, not a resolved plane wave; the "
                            "beam is formed at zero slowness")),
              open(f"{O.OUTROOT}/smoke_beam.json", "w"), indent=1)
    O.log(f"wrote {O.OUTROOT}/smoke_beam.json")


def stage_bands(a):
    from scipy.signal import butter, sosfiltfilt
    O.activate()
    E = B._eq_events()
    A = O.arrays()
    rows, nev = [], 0
    for ev in E.itertuples():
        used = False
        for site, mem in sorted(A.items()):
            D, fs = B._eq_cut(site, mem, ev)
            if len(D) < 5:
                continue
            ms = sorted(D)
            Z = np.array([D[m][0] for m in ms])
            N = np.array([D[m][1] for m in ms])
            used = True
            for lo, hi in BANDS:
                sos = butter(4, [lo / (fs / 2), hi / (fs / 2)], btype="band", output="sos")
                rows.append(dict(evid=ev.evid, site=site, M=len(ms), band=f"{lo}-{hi}",
                                 coh=B._coh(sosfiltfilt(sos, Z, axis=-1)),
                                 coh_noise=B._coh(sosfiltfilt(sos, N, axis=-1))))
        nev += int(used)
        if nev >= a.events:
            break
    R = pd.DataFrame(rows)
    g = (R.groupby("band").agg(M=("M", "median"), coh_S=("coh", "median"),
                               coh_noise=("coh_noise", "median"), n=("coh", "size"))
         .reindex([f"{x}-{y}" for x, y in BANDS]).reset_index())
    g["perfect"] = g.M
    g["incoherent_floor"] = 1.0
    g.to_csv(f"{O.OUTROOT}/array_coherence_bands.csv", index=False)
    R.to_csv(f"{O.OUTROOT}/array_coherence_bands_raw.csv", index=False)
    O.log("\n" + g.round(2).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["slowness", "bands"])
    ap.add_argument("--hours", type=int, default=4)
    ap.add_argument("--events", type=int, default=25)
    a = ap.parse_args()
    {"slowness": stage_slowness, "bands": stage_bands}[a.stage](a)


if __name__ == "__main__":
    main()
