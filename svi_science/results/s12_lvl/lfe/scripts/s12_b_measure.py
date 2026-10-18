#!/usr/bin/env python
"""s12_b_measure -- PLAN §7: per-day measurement of EVERY store_ccw series of a region for the recharge confirmation.
Families without a round-7 front take the major's median front (column has_front). Per day: dv/v in 1-2, 2-3, 3-4, 4-6,
2-5 s (joint 3C S-anchored stretch vs the series grand stack, bank +-3 %) and cc; S-P change (ms; P-window lag after S
alignment, Z); S-window spectral ratio; events per day; time since the previous / to the next front; lag to the nearest
front; between-ETS flag (> 60 d after previous, > 30 d before next).
    S11_REGION=SVI|WA_NW python s12_b_measure.py [--workers 12]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, time, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s12_epoch import C, OUT, SUB, TB, log, load_store, front_tables, fronts_for, classify, StackObs
from s12_novel_measure import dvv_only, src_scalars

_G = {}


def job(key):
    st = _G["store"][key]; fam = int(key.split("|")[0])
    days = pd.to_datetime(pd.Series(st["days"])).dt.normalize()
    dn = (days.astype("datetime64[ns]").astype(np.int64) // 86400_000_000_000).to_numpy()
    M = _G["majors"]
    fr = fronts_for(fam, M, _G["start"], _G["medrel"], _G["rel"])
    midn = np.array([_G["mid"][m].normalize().value // 86400_000_000_000 for m in M], dtype=np.int64)
    cf, jf, cg, jg, betw, ib = classify(dn, fr, midn)
    prev = classify.prev.astype(float)
    lag = dn[:, None] - fr[None, :]
    nxt = np.where(lag < 0, -lag, 10 ** 6).min(1).astype(float)
    jj = np.argmin(np.where((lag >= -10) & (lag <= 60), np.abs(lag), 10 ** 6), 1)
    lj = lag[np.arange(len(dn)), jj].astype(float); lj[(lj < -10) | (lj > 60)] = np.nan
    pm = np.where(lag >= 0, lag, 10 ** 6).argmin(1)          # index of the previous major
    St = np.nan_to_num(st["St"].astype(np.float64)); Ws = st["Ws"]
    Wt = Ws.sum(0); ref = np.full((3, St.shape[-1]), np.nan)
    for c in range(3):
        if Wt[c] > 0:
            ref[c] = (Ws[:, c] @ St[:, c]) / Wt[c]
    SO = StackObs(ref, st["tau"], st["tS"], st["sp"], np.nan)
    o = dvv_only(SO, St)
    dsp = -1000 * SO.lag(St, (0,), SO.wP)
    sc = src_scalars(st["St"], st["tau"], st["tS"], st["sp"])
    out = dict(key=np.full(len(dn), key), fam=np.full(len(dn), fam), netsta=np.full(len(dn), key.split("|")[1]),
               date=np.asarray(st["days"]), prev=prev, nxt=nxt, lag=lj, prev_major=np.array(M)[pm], betw=betw,
               has_front=np.full(len(dn), fam in _G["fams"]), ws=Ws.mean(1),
               nd=st["Nd"] if st.get("Nd") is not None else np.full(len(dn), -1), dSP_ms=dsp, sr=sc["sr"],
               **{f"dvv_{k}": v for k, v in o.items()})
    return pd.DataFrame(out)


def init(g):
    _G.update(g)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    t0 = time.time()
    store = load_store()
    FT, majors, start, medrel, mid, rel = front_tables()
    g = dict(store=store, majors=majors, start=start, medrel=medrel, mid=mid, rel=rel, fams=set(FT.fam.astype(int)))
    keys = sorted(store)
    log(f"[b_measure] {C.REGION}: {len(keys)} series ({len({k.split('|')[0] for k in keys})} families, "
        f"{len(g['fams'])} with fronts)")
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.workers, initializer=init, initargs=(g,)) as p:
        F = pd.concat(list(p.imap_unordered(job, keys, chunksize=4)), ignore_index=True)
    F.to_parquet(f"{TB}/b_days.parquet", index=False)
    open(f"{OUT}/BMEAS_{SUB.upper()}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"[b_measure] {len(F)} days, between {int(F.betw.sum())} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
