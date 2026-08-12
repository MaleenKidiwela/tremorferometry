"""Convert per-day stacks into BACKWARD/trailing N-day moving-average stacks.

For each family and each anchor date t, combine all daily stacks in the window (t - N days, t]
(causal — only the PAST N days, never the future), count-weighted, then L2-normalize:
    trailing(t) = normalize( sum_{t-N < d <= t} daily_stack_d * n_det_d )
This is a higher-SNR, smoother stack series with no future leakage (good for precursor work).
Operates on a daily-stack npz (from build_long_window_*), writes the same schema.
"""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in-npz", default="data/long_window_daily_GNW.npz")
    p.add_argument("--out-npz", default="data/long_window_trail30_GNW.npz")
    p.add_argument("--window-days", type=int, default=30)
    p.add_argument("--min-det", type=int, default=20, help="min TOTAL detections in trailing window")
    p.add_argument("--epoch-bounds", default="2010-09-10,2019-05-07",
                   help="instrument-change dates; a trailing window never pools across one "
                        "(GNW: 50Hz<2010-09-10, 40Hz, 100Hz>=2019-05-07). '' to disable.")
    args = p.parse_args()

    d = np.load(args.in_npz, allow_pickle=True)
    S = d["stacks"]; patches = d["patches"]; dates = pd.to_datetime(d["dates"]).values.astype("datetime64[D]")
    ndet = d["n_det"].astype(np.float64); t = d["t"]; fs = float(d["fs"])
    N = np.timedelta64(args.window_days, "D")
    bounds = [np.datetime64(b) for b in args.epoch_bounds.split(",") if b.strip()]
    def _epoch(dd):  # epoch index = how many transition dates the day is on/after
        e = np.zeros(len(dd), int)
        for b in bounds:
            e += (dd >= b).astype(int)
        return e
    print(f"[in] {len(S)} daily stacks, {len(np.unique(patches))} patches; trailing {args.window_days} d; "
          f"epoch bounds {[str(b) for b in bounds]}")

    out_S=[]; out_p=[]; out_d=[]; out_n=[]
    for fam in np.unique(patches):
        m = patches == fam
        s = S[m]; dt = dates[m]; nd = ndet[m]
        order = np.argsort(dt); s=s[order]; dt=dt[order]; nd=nd[order]
        ep = _epoch(dt)
        # run_start[i] = first index of the contiguous same-epoch run containing i
        run_start = np.zeros(len(dt), int); start = 0
        for i in range(len(dt)):
            if i > 0 and ep[i] != ep[i-1]:
                start = i
            run_start[i] = start
        lo = 0
        for i in range(len(dt)):
            # trailing window (dt[i]-N, dt[i]]  -> advance lo
            while dt[lo] <= dt[i] - N:
                lo += 1
            lo_eff = max(lo, run_start[i])   # never cross an instrument-change boundary
            w = nd[lo_eff:i+1]
            tot = float(w.sum())
            if tot < args.min_det:
                continue
            acc = (s[lo_eff:i+1] * w[:, None]).sum(axis=0)
            nrm = np.linalg.norm(acc)
            if nrm <= 0:
                continue
            out_S.append((acc/nrm).astype(np.float32)); out_p.append(fam)
            out_d.append(str(pd.Timestamp(dt[i]).date())); out_n.append(int(tot))
    out_S=np.array(out_S,np.float32); out_p=np.array(out_p); out_d=np.array(out_d); out_n=np.array(out_n,np.int32)
    ordr = np.lexsort((out_d, out_p))
    np.savez(args.out_npz, stacks=out_S[ordr], patches=out_p[ordr], dates=out_d[ordr],
             n_det=out_n[ordr], t=t, fs=fs)
    print(f"[out] {len(out_S)} trailing stacks -> {args.out_npz}")
    print(f"  median detections/trailing-stack: {int(np.median(out_n))}, vs daily median {int(np.median(ndet))}")

if __name__ == "__main__":
    main()
