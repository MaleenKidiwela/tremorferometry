#!/usr/bin/env python
"""P3D DISCOVERY ARM C -- frozen-CF trigger: full-day CF, then the top peaks/station-day.

The CF is the frozen P2 model (`results/p2/models/p2nnpu.pt`) through `p2_cf.cf_day_sid`, every
station on the UNKNOWN-station embedding, days read from the 20 Hz caches with the component dip
honoured.  Candidates are the top `--peaks-per-day` CF peaks per station-day at the <=60 km sites,
declustered by `--min-sep-s`, pooled over sites and merged at `--merge-s`, which is what Merlin's
ARM C specifies as the trigger for the ARM-A chain.

  python scripts/s_discover_cf.py --net B --peaks-per-day 720
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_discover_pairs as DP                                              # noqa: E402

MODEL = f"{S.PICK}/results/p2/models/p2nnpu.pt"
OUT = S.P3D
# one station per site (the site's best-registering member), <=60 km sites of NET-B
SITE_STA = {"PGCq": "PGC", "VGZ": "VGZ", "SYMB": "SYMB", "BUCKS": "BUCKS",
            "LHLYB": "LHLYB", "ALBHB": "ALBHB", "DAOB": "DAOB"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="B")
    ap.add_argument("--peaks-per-day", type=int, default=720)
    ap.add_argument("--min-sep-s", type=float, default=6.0)
    ap.add_argument("--merge-s", type=float, default=3.0)
    ap.add_argument("--out", default=f"{OUT}/disc_C_candidates.parquet")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    import p_models as M
    import p2_cf as CF
    t0 = time.time()
    stas = [SITE_STA[q] for q in SITE_STA if SITE_STA[q] in DP.NETS[a.net]]
    C, idx, days, cat_end = DP.candidates(a.net, DP.EPISODES)
    daykeys = sorted({(int(y), int(d)) for y, d in zip(C.year, C.doy)})
    hours = set(pd.to_datetime(C.t.values, unit="s").floor("h"))
    z = torch.load(MODEL, map_location="cuda", weights_only=False)
    model = M.PhaseNetS(nsta=z.get("nsta", 0))
    model.load_state_dict(z["state"])
    model.cuda().eval()
    unk = int(z.get("nsta", 15))
    sep = int(round(a.min_sep_s * S.FS))
    idxall = S.cache_index(years=list(range(2019, 2027)), stas=stas)
    rows = []
    for k, (y, d) in enumerate(daykeys):
        for s in stas:
            v = S.load_day_zne(s, y, d, idxall)
            if v is None:
                continue
            t0d, x = v
            cf = CF.cf_day_sid(model, z, np.asarray(x, np.float32), unk)
            # keep only samples inside tremor hours, then take the top peaks with a min separation
            n = len(cf)
            tt = t0d + np.arange(n) / S.FS
            inh = pd.Series(pd.to_datetime(tt, unit="s")).dt.floor("h").isin(hours).to_numpy()
            c = np.where(inh, cf, 0.0)
            order = np.argsort(c)[::-1]
            taken, out = np.zeros(n, bool), []
            for i in order:
                if c[i] <= 0 or len(out) >= a.peaks_per_day:
                    break
                lo, hi = max(i - sep, 0), min(i + sep + 1, n)
                if taken[lo:hi].any():
                    continue
                taken[i] = True
                out.append(i)
            for i in out:
                rows.append((float(tt[i]), s, float(c[i])))
        if (k + 1) % 10 == 0:
            print(f"  [cf] {k+1}/{len(daykeys)} days  {len(rows)} peaks  {time.time()-t0:.0f}s",
                  flush=True)
    P = pd.DataFrame(rows, columns=["t", "sta", "cf"]).sort_values("t").reset_index(drop=True)
    # merge peaks that coincide across sites into single candidate times
    tv = P.t.to_numpy()
    grp = np.zeros(len(tv), np.int64)
    g = 0
    for i in range(1, len(tv)):
        if tv[i] - tv[i - 1] > a.merge_s:
            g += 1
        grp[i] = g
    P["grp"] = grp
    Cand = P.groupby("grp").agg(t=("t", "median"), n_sites=("sta", "nunique"),
                                cf_max=("cf", "max")).reset_index(drop=True)
    Cand = Cand.sort_values("t").reset_index(drop=True)
    Cand[["t", "n_sites", "cf_max"]].to_parquet(a.out, index=False)
    json.dump(dict(net=a.net, stations=stas, n_days=len(daykeys), n_peaks=int(len(P)),
                   n_candidates=int(len(Cand)), peaks_per_day=a.peaks_per_day,
                   min_sep_s=a.min_sep_s, merge_s=a.merge_s, runtime_s=time.time() - t0),
              open(f"{OUT}/disc_C_candidates.json", "w"), indent=1, default=str)
    print(f"[cf] {len(P)} peaks -> {len(Cand)} merged candidates -> {a.out}  "
          f"({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
