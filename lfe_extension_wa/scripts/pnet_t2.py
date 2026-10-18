#!/usr/bin/env python
"""Stage D: T2 withheld-station template recovery with PNET-derived candidate times (instead of T2's permanent MF).

Same pairs, truth templates, ceiling and scorer as merlin_t2.py (results/merlin_tests/t2; linked into OD).
Candidates per (family, day): a PNET "beam" on the origin axis over the family's NEAREST permanent sites,
    B_f(i) = sum over the n nearest sites of max over the site's stations of maxfilt_0.5s(CF_sta)(i + offset_s(f, sta)),
sites = those with a >= 30-event build template of f on the day's node and a CF (distance family -> station, lat_use/lon_use);
n = 3 (series 'data'), n = 1 (series 'data_n1').  Picks = 4-s block maxima >= both neighbours, top K per day with
K = 2 per scored hour (T2's rule), origin in the TEST half only.  Scrambled control ('control', 'control_n1'): each
site's shifted CF rolled by its own shift in [300, 3000] s per day before the sum.  The withheld array station is stacked
at the picks exactly as in T2 (1-8 Hz cache, [t-20, t+32] s, unit max per window, valid windows only) and the round-0
template is cut with merlin_t2.tpl_from_stack.  Scoring: `MT2_OD=<OD> python scripts/merlin_t2.py --stage score`.

  python scripts/pnet_t2.py --cf w05 --od results/pnet/t2_w05
"""
from __future__ import annotations
import os, sys, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

T2 = f"{W.RES}/merlin_tests/t2"
NL = int(W.LONG_LEN * W.FS)
IPRE = int(W.LONG_PRE * W.FS)
DECL = int(4.0 * W.FS)
SER = ["data", "control", "data_n1", "control_n1"]
_G = {}


def _init(cf, od):
    import wa_mf
    T, _ = wa_mf.load_tpl(["build"], all_templates=True)
    S = W.load_stations()
    F = W.load_families()
    TR = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    _G.update(TOFF={k: v["off"] for k, v in T.items()}, SITE=dict(zip(S.netsta, S.site)),
              SLL={r.netsta: (r.lat, r.lon) for r in S.itertuples()}, FLL={r.fid: (r.lat_use, r.lon_use) for r in F.itertuples()},
              TR={d: g for d, g in TR.groupby("date")}, cf=cf, P=pd.read_csv(f"{od}/pairs.csv"))


def pick(v, ok, k):
    NB = (len(v) + DECL - 1) // DECL
    x = np.full(NB * DECL, -np.inf); x[:len(v)] = np.where(ok, v, -np.inf)
    xb = x.reshape(NB, DECL)
    ba = xb.argmax(1); bm = xb[np.arange(NB), ba]
    left = np.r_[-np.inf, bm[:-1]]; right = np.r_[bm[1:], -np.inf]
    bm2 = np.where((bm >= left) & (bm > right), bm, -np.inf)
    o = np.argsort(bm2)[::-1][:k]
    o = o[np.isfinite(bm2[o])]
    return np.sort(o * DECL + ba[o])


def day_job(args):
    w, ds = args
    t0 = time.time()
    TOFF, SITE, SLL, FLL = _G["TOFF"], _G["SITE"], _G["SLL"], _G["FLL"]
    P = _G["P"]; P = P[P.window == w]
    D0 = pd.Timestamp(ds).timestamp()
    NV = W.NDAY - W.NT + 1
    tt = D0 + np.arange(NV) / W.FS
    scored = W.half_of(tt) == 1
    K = int(round(2.0 * scored.sum() / (3600 * W.FS)))
    S = W.load_stations()
    MF, NODE = {}, {}
    for ns in S.netsta:
        cf = PC.load_cf(_G["cf"], ns, ds)
        p = W.cache_path(ns, ds)
        if cf is None or not os.path.exists(p):
            continue
        MF[ns] = PC.maxfilt(cf, 10); NODE[ns] = str(np.load(p)["node"])
    rng = np.random.default_rng(20260930 + 11 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    sites_all = sorted({SITE[ns] for ns in MF})
    shift = {q: int(rng.uniform(*W.ROLL_S) * W.FS) for q in sites_all}
    tg = _G["TR"].get(ds, None)
    t_any = np.sort(tg.t.to_numpy()) if tg is not None else np.zeros(0)
    ARR = {}
    acc, cands = {}, []
    for f, gp in P.groupby("fid"):
        sts = [ns for ns in MF if (int(f), ns, NODE[ns]) in TOFF]
        if not sts:
            continue
        la, lo = FLL.get(f, (np.nan, np.nan))
        dist = {ns: float(W.km(la, lo, *SLL[ns])) for ns in sts}
        bysite = {}
        for ns in sts:
            bysite.setdefault(SITE[ns], []).append(ns)
        order = sorted(bysite, key=lambda q: min(dist[s] for s in bysite[q]) if np.isfinite(la) else 0.0)
        own_all = np.sort(tg[tg.fid == f].t.to_numpy()) if tg is not None else np.zeros(0)
        beams = {}
        for nsite, suf in [(3, ""), (1, "_n1")]:
            use = order[:nsite]
            if len(use) < min(nsite, 1):
                continue
            bd = np.zeros(NV); bc = np.zeros(NV)
            for q in use:
                v = np.zeros(NV)
                for ns in bysite[q]:
                    d = int(round(TOFF[(int(f), ns, NODE[ns])] * W.FS))
                    a0, b0 = max(0, -d), min(NV, W.NDAY - d)
                    tmp = np.zeros(NV); tmp[a0:b0] = MF[ns][a0 + d:b0 + d]
                    v = np.maximum(v, tmp)
                bd += v; bc += np.roll(v, shift[q])
            beams["data" + suf] = bd; beams["control" + suf] = bc
        for sr, b in beams.items():
            ix = pick(b, scored, K)
            tsel = D0 + ix / W.FS
            po = int(np.sum(np.min(np.abs(own_all[None] - tsel[:, None]), 1) <= W.MATCH_S)) if len(own_all) and len(tsel) else 0
            pa = int(np.sum(np.min(np.abs(t_any[None] - tsel[:, None]), 1) <= W.MATCH_S)) if len(t_any) and len(tsel) else 0
            for p in gp.itertuples():
                if p.netsta not in ARR:
                    pp = W.cache_path(p.netsta, ds)
                    if os.path.exists(pp):
                        z = np.load(pp); ARR[p.netsta] = (z["x"], np.concatenate([[0], np.cumsum(~z["m"])]), z["live"])
                    else:
                        ARR[p.netsta] = None
                A = ARR[p.netsta]
                e = acc.setdefault((int(p.pair), sr), [np.zeros((3, NL)), 0, 0, 0, 0, np.zeros(3, bool)])
                e[2] += po; e[3] += pa; e[4] += len(tsel)
                if A is None:
                    continue
                x, cm, live = A
                st = ix - IPRE
                st = st[(st >= 0) & (st + NL <= W.NDAY)]
                st = st[(cm[st + NL] - cm[st]) == 0]
                if len(st):
                    seg = x[:, st[:, None] + np.arange(NL)[None]].transpose(1, 0, 2).astype(np.float64)
                    mx = np.abs(seg).max(axis=(1, 2)); g = mx > 0
                    e[0] += (seg[g] / mx[g, None, None]).sum(0); e[1] += int(g.sum())
                e[5] |= live
                cands.append((int(p.pair), sr, tsel))
    return ds, acc, cands, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True)
    ap.add_argument("--od", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--days", default="")
    a = ap.parse_args()
    os.makedirs(a.od, exist_ok=True)
    for f in ["pairs.csv", "truth_qc.csv", "truth_stacks.npz", "ceiling_qc.csv", "ceiling.npz"]:
        dst = f"{a.od}/{f}"
        if not os.path.exists(dst):
            os.symlink(f"{T2}/{f}", dst)
    import merlin_t2 as MT
    days = [(w, d) for w in MT.WINS for d, _ in W.window_days(w)]
    if a.days:
        days = [x for x in days if x[1] in set(a.days.split(","))]
    W.log(PC.LOG, f"[t2:{a.cf}] {len(days)} days -> {a.od}")
    acc, cands = {}, []
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(a.cf, a.od)) as ex:
        fut = [ex.submit(day_job, x) for x in days]
        for fu in as_completed(fut):
            ds, ac, cd, secs = fu.result()
            for k, v in ac.items():
                e = acc.setdefault(k, [np.zeros((3, NL)), 0, 0, 0, 0, np.zeros(3, bool)])
                e[0] += v[0]; e[1] += v[1]; e[2] += v[2]; e[3] += v[3]; e[4] += v[4]; e[5] |= v[5]
            cands += cd
            W.log(PC.LOG, f"[t2:{a.cf}] {ds} {len(ac)} (pair, series) {secs:.0f}s")
    ref = MT.wref()
    rows, tpls, stks = [], [], []
    for (pr, sr), (s_, c_, po, pa, nc, lv) in sorted(acc.items()):
        row = dict(pair=pr, series=sr, round=0, n_stacked=c_, n_cand=nc, pur_own=po / nc if nc else np.nan,
                   pur_any=pa / nc if nc else np.nan)
        st = s_ / c_ if c_ else np.zeros((3, NL))
        t = MT.tpl_from_stack(st, lv, ref) if c_ >= 1 else None
        if t is not None:
            row.update({k: v for k, v in t.items() if k != "tpl"}); tp = t["tpl"]
        else:
            tp = np.zeros((3, W.NT), np.float32)
        row["row"] = len(tpls)
        rows.append(row); tpls.append(tp); stks.append(st.astype(np.float32))
    pd.DataFrame(rows).to_csv(f"{a.od}/round0_qc.csv", index=False)
    np.savez(f"{a.od}/round0.npz", tpl=np.stack(tpls), stack=np.stack(stks))
    if cands:
        pd.DataFrame(dict(pair=np.concatenate([[p] * len(t) for p, _, t in cands]).astype(int),
                          series=np.concatenate([[s] * len(t) for _, s, t in cands]),
                          t=np.concatenate([t for _, _, t in cands]))).to_parquet(f"{a.od}/round0_cand.parquet", index=False)
    R = pd.DataFrame(rows)
    W.log(PC.LOG, f"[t2:{a.cf}] wrote {len(R)} (pair, series); median n_stacked {R.groupby('series').n_stacked.median().to_dict()}; "
          f"median purity own {R.groupby('series').pur_own.median().to_dict()}")


if __name__ == "__main__":
    main()
