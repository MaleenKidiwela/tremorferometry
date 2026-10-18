#!/usr/bin/env python
"""SH_V2 template (re)stacking (2026-09-28) -- new file, modelled on t9_templates_aa.py / aa_restack.py.

Every window is read through sh_v2_common.DayPrep (2-9 Hz, burst-muted, line-notched) and enters a
channel's stack only if it passes the same gap guard as the detector (every sample valid, e > 1e-3 x
the channel's reference energy, channel/station not dead).  aa_restack.py had no guard: it
rms-normalised dead-channel residue and glitch windows to unit amplitude and stacked them, and 100 %
of its seeds came from the broken family-days (SH_V2_LOG step 1).

Offsets are the stage-1 ones carried by the base table (tt_s + best_tau, cut at -pre) and are not
re-measured; 6 s windows with 2 s pre-arrival (row C of AA).  Each channel-window is normalised by its
own rms, then averaged.  A family is stacked ONLY on its OWN build days (stage2_days_build.json
family key), so its held-out test days are never seen.

  --times shelly : iteration 0 -- his published times on the family's build days
  --det X.parquet: iteration k -- our detections on the family's build days, the --nstack highest
                   correlation sums at >= --thr x MAD (Shelly sec. 2: "stacking the detections with
                   the highest correlation sum"; N is SILENT in the main text -> 200, as agent AA used)
"""
import os, sys, argparse, warnings, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v2_common import *

_G = {}
def _init(G):
    _G.update(G)


def _day(d):
    warnings.filterwarnings("ignore")
    E, off, NC, stas, PRE, L = (_G[k] for k in ["E", "off", "NC", "stas", "PRE", "L"])
    acc = {}; cnt = {}
    P = DayPrep(d)
    if P.z is None: return acc, cnt
    Ed = E[E.day == d]
    for ns in stas:
        fk = [f for f in Ed.fam.unique() if (f, ns) in off]
        if not fk: continue
        v = P.load(ns)
        if v is None: continue
        t0, X, mc, ref = v
        for f in fk:
            if X.shape[0] != NC[(f, ns)]: continue
            tt = Ed[Ed.fam == f].t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
            idx = np.round((tt + off[(f, ns)] - PRE - t0) * FS).astype(int)
            idx = idx[(idx >= 0) & (idx + L <= X.shape[1])]
            if not len(idx): continue
            ok = window_ok(X, mc, ref, L, idx)                # (C, n)
            W = X[:, idx[:, None] + np.arange(L)]              # (C, n, L)
            rms = np.sqrt((W ** 2).mean(-1)) + 1e-30
            W = np.where(ok[:, :, None], W / rms[:, :, None], 0.0)
            k = (f, ns)
            acc[k] = acc.get(k, 0) + W.sum(1); cnt[k] = cnt.get(k, 0) + ok.sum(1)
    return acc, cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="template table (csv/npz stem) giving offsets / channel counts")
    ap.add_argument("--times", default="", choices=["", "shelly"]); ap.add_argument("--det", default="")
    ap.add_argument("--days-file", required=True); ap.add_argument("--nstack", type=int, default=200)
    ap.add_argument("--thr", type=float, default=9.0); ap.add_argument("--out", required=True)
    ap.add_argument("--source", default="Shelly2025"); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--pre", type=float, default=2.0); ap.add_argument("--len", type=float, default=6.0)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    M = pd.read_csv(a.base + ".csv"); M = M[M.source == a.source].copy(); M["fam"] = M.fam.astype(str)
    Zb = np.load(a.base + ".npz") if os.path.exists(a.base + ".npz") else None
    dd = json.load(open(a.days_file))
    FD = {k.split("|")[1]: {pd.Timestamp(x) for x in v} for k, v in dd.items() if k.split("|")[0] == a.source}
    if a.times == "shelly":
        S = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
        S.columns = [c.strip('%') for c in S.columns]
        S["t"] = (pd.to_datetime(dict(year=S.year, month=S.month, day=S.day, hour=S.hour, minute=S.minute,
                  second=S.second.astype(int))).astype("datetime64[ns]") + pd.to_timedelta(S.second % 1, unit="s"))
        E = pd.DataFrame(dict(fam=S.idnum.astype(str), t=S.t, score=S.meancc))
    else:
        D = pd.read_parquet(a.det); D["fam"] = D.fam.astype(str)
        E = D[D.net_mad >= a.thr][["fam", "t", "net_mad"]].rename(columns={"net_mad": "score"})
    E["day"] = E.t.dt.floor("D")
    E = E[[f in FD and d in FD[f] for f, d in zip(E.fam, E.day)]]
    if a.det:
        E = E.sort_values("score", ascending=False).groupby("fam").head(a.nstack)
    log(f"[restack-v2] {len(E)} seed times over {E.fam.nunique()} families (own build days only); "
        f"per family {E.groupby('fam').size().describe()[['min', '50%', 'max']].to_dict()}")
    PRE, LEN = a.pre, a.len; L = int(round(LEN * FS))
    off = {(r.fam, r.netsta): float(r.tt_s) + float(r.best_tau) for r in M.itertuples()}
    NC = {(r.fam, r.netsta): int(r.ncomp) for r in M.itertuples()}
    ACC = {}; CNT = {}
    stas = sorted(M.netsta.unique())
    days = sorted(E.day.unique())
    G = dict(E=E, off=off, NC=NC, stas=stas, PRE=PRE, L=L)
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(G,)) as ex:
        for i, (acc, cnt) in enumerate(ex.map(_day, days)):
            for k in acc:
                ACC[k] = ACC.get(k, 0) + acc[k]; CNT[k] = CNT.get(k, 0) + cnt[k]
            if i % 50 == 0: log(f"  {i + 1}/{len(days)} stacks={len(ACC)}")
    TPL = {}; rows = []
    for r in M.itertuples():
        k = (r.fam, r.netsta)
        if k in ACC and np.any(CNT[k] > 0):
            T = ACC[k] / np.maximum(CNT[k], 1)[:, None]
            T[CNT[k] == 0] = 0.0
            src = "stacked"; n = CNT[k]
        elif a.det and Zb is not None and r.key in Zb:
            T = Zb[r.key].astype(np.float64); src = "kept_previous"; n = np.zeros(T.shape[0], int)
        else:
            continue
        TPL[r.key] = T.astype(np.float32)
        d_ = {c: getattr(r, c) for c in ["source", "fam", "net", "sta", "netsta", "key", "ncomp", "best_tau", "hypo_km", "tt_s"]}
        d_.update(pre_s=PRE, len_s=LEN, n_stack=int(np.max(n)) if len(n) else 0,
                  n_stack_ch=",".join(str(int(x)) for x in n), src=src)
        rows.append(d_)
    T = pd.DataFrame(rows); T.to_csv(a.out.replace(".npz", ".csv"), index=False)
    np.savez_compressed(a.out, **TPL)
    log(f"[restack-v2] {len(TPL)} templates ({T.src.value_counts().to_dict()}) over {T.fam.nunique()} families; "
        f"median n_stack {T.n_stack.median():.0f} -> {a.out}")
if __name__ == "__main__": main()
