#!/usr/bin/env python
"""T2B — T2 (withheld-array-station template recovery) with candidate times chosen by the CONFIDENCE model
instead of raw matched-filter peaks.  2026-09-30.

For every T2 pair (family f, window w, withheld XU/XG station s; results/merlin_tests/t2/pairs.csv):
  pool      the family's step-1 candidates (T1 all-templates PERMANENT-network MF, >= 8 MAD) in the TEST half of w,
            each with its out-of-fold confidence P (D2 own-family model, 5 folds grouped by day:
            results/confidence/candidates_confidence.parquet column D2_hgb_FULL_own).
  leakage   the MF, the step-1 features, PNET (off-the-shelf + fine-tuned) and the step-2 attribution are all
            computed on the 42-station permanent roster (UW/PB/CN/US; results/setup/stations.csv); no XU/XG station
            enters any of them, so s contributes nothing to the candidate times or to P -> no refit needed.  The
            PNET fine-tune used build-half windows of 36 permanent stations only.  D2 P of a test-half row comes
            from a model that never saw that day.  (The model's training LABELS on other days are array-catalogue
            matches; the catalogue was built with the array, so s helped define labels on OTHER days only.)
  variants  topN    top-N by P, N = the number of candidates T2 round 0 used for the pair (min(N, pool))
            P05     every pool candidate with P >= 0.5
            P03     every pool candidate with P >= 0.3
            mad_P03 / mad_P05   the same COUNT as P03 / P05, chosen by net_mad (raw MF strength) — isolates the
                    ranking gain of P from the pool itself
  control   the same selection on the scrambled control c1 (permanent sites shifted to other times of the day)
            with the same model's out-of-fold P; s is stacked unshifted at those times.
  stack     s's 20-Hz 1-8 Hz cache, [t-20, t+32] s windows, unit-max over 3 components, windows touching invalid
            data skipped (merlin_t2.stack_at rule).
  score     merlin_t2 compare rule: 2-8 Hz Butterworth-4 designed at fs = 20 Hz (the cache rate); cc_al = best lag
            +-2 s vs the own truth template (build half); rank among all families' truth templates at s.

  python scripts/t2b.py --workers 12      -> results/t2b/
"""
from __future__ import annotations
import os, sys, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
from merlin_t2 import bp28, tpl28, slide_cc, NL, IPRE

T2D = f"{W.RES}/merlin_tests/t2"
OD = f"{W.RES}/t2b"
LOG = f"{W.STUDY}/logs/t2b.log"
PCOL = "D2_hgb_FULL_own"
VARIANTS = ["topN", "P05", "P03", "mad_P05", "mad_P03"]
RS = {"Royer_NW_WA": "NW", "Royer_SW_WA": "SW", "Olympic_2010": "Olympic", "Olympic_2010_armB": "Arm B"}


def log(s):
    W.log(LOG, s)


def select(P, T2Q):
    """(pair, series, variant) -> sorted candidate times"""
    C = pd.read_parquet(f"{W.RES}/confidence/candidates_confidence.parquet",
                        columns=["fid", "t", "series", "window", "half", PCOL, "net_mad"])
    C = C[C.window.isin(["V2008", "V2010"]) & (C.half == 1) & C.series.isin(["data", "c1"])]
    C = C.rename(columns={PCOL: "P"})
    G = {k: g for k, g in C.groupby(["fid", "window", "series"])}
    ncand = {(int(x.pair), x.series): int(x.n_cand) for x in T2Q.itertuples()}
    sel, rows = {}, []
    for p in P.itertuples():
        for sr, cs in [("data", "data"), ("control", "c1")]:
            g = G.get((p.fid, p.window, cs))
            N = ncand.get((int(p.pair), sr), 0)
            pool = 0 if g is None else len(g)
            out = {v: np.zeros(0) for v in VARIANTS}
            if g is not None and len(g):
                gp = g.sort_values("P", ascending=False)
                gm = g.sort_values("net_mad", ascending=False)
                out["topN"] = gp.t.to_numpy()[:N]
                n05 = int((g.P >= 0.5).sum()); n03 = int((g.P >= 0.3).sum())
                out["P05"] = gp.t.to_numpy()[:n05]; out["P03"] = gp.t.to_numpy()[:n03]
                out["mad_P05"] = gm.t.to_numpy()[:n05]; out["mad_P03"] = gm.t.to_numpy()[:n03]
            for v in VARIANTS:
                sel[(int(p.pair), sr, v)] = np.sort(out[v])
                rows.append(dict(pair=int(p.pair), series=sr, variant=v, pool=pool, N_t2=N, n_cand=len(out[v]),
                                 meanP=float(g.set_index("t").P.reindex(out[v]).mean()) if len(out[v]) else np.nan))
    return sel, pd.DataFrame(rows)


def stack_job(job):
    """one (window, station): stack every (pair, series, variant) selection over its days"""
    w, ns, items = job      # items: list of (key, times)
    acc = {k: [np.zeros((3, NL)), 0, np.zeros(3, bool)] for k, _ in items}
    days = {}
    for k, tt in items:
        if not len(tt):
            continue
        d = pd.to_datetime(tt, unit="s").strftime("%Y-%m-%d").to_numpy()
        for ds in np.unique(d):
            days.setdefault(ds, []).append((k, tt[d == ds]))
    for ds, lst in sorted(days.items()):
        pth = W.cache_path(ns, ds)
        if not os.path.exists(pth):
            continue
        try:
            z = np.load(pth)
            X, m, live = z["x"], z["m"], z["live"]
        except Exception:
            continue
        D0 = pd.Timestamp(ds).timestamp()
        cm = np.concatenate([[0], np.cumsum(~m)])
        for k, tt in lst:
            i0 = np.round((tt - D0) * W.FS).astype(int) - IPRE
            ok = (i0 >= 0) & (i0 + NL <= W.NDAY)
            ok[ok] = (cm[i0[ok] + NL] - cm[i0[ok]]) == 0
            if not ok.any():
                continue
            seg = X[:, i0[ok][:, None] + np.arange(NL)[None]].transpose(1, 0, 2).astype(np.float64)
            mx = np.abs(seg).max(axis=(1, 2))
            g = mx > 0
            if not g.any():
                continue
            e = acc[k]
            e[0] += (seg[g] / mx[g, None, None]).sum(0); e[1] += int(g.sum()); e[2] |= live.astype(bool)
    return {k: (v[0] / v[1] if v[1] else None, v[1], v[2]) for k, v in acc.items()}


def build_truth(P):
    TQ = pd.read_csv(f"{T2D}/truth_qc.csv", dtype={"live": str})
    need = set(zip(P.window, P.netsta))
    B = TQ[(TQ.half == 0) & (TQ.n_events >= W.MIN_EVENTS) & TQ.offset_s.notna()].copy()
    B = B[[(w, s) in need for w, s in zip(B.window, B.netsta)]]
    B["live"] = B.live.astype(str).str.zfill(3)
    zt = np.load(f"{T2D}/truth_stacks.npz")
    TK = {k: i for i, k in enumerate(zt["keys"])}
    TS = zt["stacks"]
    TT = {}
    for (w, ns), g in B.groupby(["window", "netsta"]):
        lst = []
        for x in g.itertuples():
            k = f"{x.fid}|{w}|{ns}|0"
            if k not in TK:
                continue
            lv = np.array([c == "1" for c in x.live])
            t, ok = tpl28(TS[TK[k]], x.offset_s, lv)
            lst.append((int(x.fid), t, ok, int(round(x.offset_s * W.FS)) + IPRE))
        TT[(w, ns)] = lst
    return TT


def compare(p, x_long, live_x, TT, dk, fi):
    lag2 = int(W.MATCH_S * W.FS)
    j_lo, j_hi = IPRE + int(-5.5 * W.FS), IPRE + int(20 * W.FS)
    x28 = bp28(x_long)
    own, res = None, []
    for fid, t, ok, j0 in TT.get((p.window, p.netsta), []):
        c = slide_cc(t, ok, x28, live_x)
        if c is None:
            continue
        cany = float(c[j_lo:j_hi + 1].max())
        res.append((fid, cany))
        if fid == p.fid:
            a0, b0 = max(0, j0 - lag2), min(len(c) - 1, j0 + lag2)
            own = (float(c[a0:b0 + 1].max()), cany)
    if own is None:
        return {}
    res.sort(key=lambda v: -v[1])
    rank = 1 + sum(1 for f_, c_ in res if f_ != p.fid and c_ > own[1])
    best = res[0][0]
    bd = dk[fi[p.fid], fi[best]] if best in fi else np.nan
    return dict(cc_al=own[0], cc_any=own[1], rank=rank, n_competitors=len(res), best_fid=best, best_is_own=best == p.fid,
                best_is_nbr3=bool(best == p.fid or (np.isfinite(bd) and bd <= 3)))


def purity(sel, P):
    T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    own = {f: np.sort(g.t.to_numpy()) for f, g in T.groupby("fid")}
    fid = dict(zip(P.pair, P.fid))
    out = {}
    for (pr, sr, v), tt in sel.items():
        r = own.get(fid[pr], np.zeros(0))
        if not len(tt) or not len(r):
            out[(pr, sr, v)] = np.nan if not len(tt) else 0.0
            continue
        j = np.clip(np.searchsorted(r, tt), 1, max(len(r) - 1, 1))
        d = np.minimum(np.abs(r[j] - tt), np.abs(r[j - 1] - tt)) if len(r) > 1 else np.abs(r[0] - tt)
        out[(pr, sr, v)] = float((d <= W.MATCH_S).mean())
    return out


def summarise(O, T2):
    """per region x series x variant; fractions over ALL pairs (no stack = not recovered); median cc over stacked pairs;
    paired T2 round-0 and ceiling on the same stacked pairs"""
    t2d = T2[(T2.series == "data") & (T2["round"] == 0)].set_index("pair")
    t2c = T2[(T2.series == "control") & (T2["round"] == 0)].set_index("pair")
    ceil = T2[T2.series == "ceiling"].set_index("pair")
    rows = []
    O = O.assign(reg=O.region.map(RS))
    for rg, g0 in [("ALL", O)] + [(r, O[O.reg == r]) for r in ["NW", "SW", "Olympic", "Arm B"]]:
        for (sr, v), h in g0.groupby(["series", "variant"]):
            st = h[h.n_stacked > 0]
            ref = t2d if sr == "data" else t2c
            rp = ref.reindex(st.pair); cp = ceil.reindex(st.pair)
            rows.append(dict(region=rg, series=sr, variant=v, n_pairs=len(h), n_families=h.fid.nunique(),
                             frac_stacked=len(st) / len(h), med_windows=st.n_stacked.median() if len(st) else 0,
                             med_meanP=h.meanP.median(), med_purity_own=st.pur_own.median() if len(st) else np.nan,
                             med_cc=st.cc_al.median() if len(st) else np.nan,
                             frac_cc07=float((h.cc_al > 0.7).sum() / len(h)), own_best=float(h.best_is_own.fillna(False).astype(bool).sum() / len(h)),
                             top3=float((h["rank"] <= 3).sum() / len(h)), best_nbr3=float(h.best_is_nbr3.fillna(False).astype(bool).sum() / len(h)),
                             T2_med_cc_same=rp.cc_al.median(), T2_own_best_same=rp.best_is_own.astype(float).mean(),
                             T2_top3_same=(rp["rank"] <= 3).mean(), ceil_med_cc_same=cp.cc_al.median(),
                             frac_cc07_stacked=(st.cc_al > 0.7).mean() if len(st) else np.nan,
                             own_best_stacked=st.best_is_own.astype(float).mean() if len(st) else np.nan,
                             top3_stacked=(st["rank"] <= 3).mean() if len(st) else np.nan))
    # T2 references over all pairs
    for rg in ["ALL", "NW", "SW", "Olympic", "Arm B"]:
        for sr, lab in [("data", "T2_round0"), ("control", "T2_round0"), ("ceiling", "T2_ceiling")]:
            h = T2[(T2.series == sr) & (T2["round"].isin([0, -1]))].assign(reg=lambda d: d.region.map(RS))
            h = h if rg == "ALL" else h[h.reg == rg]
            rows.append(dict(region=rg, series=sr, variant=lab, n_pairs=len(h), n_families=h.fid.nunique(), frac_stacked=1.0,
                             med_windows=h.n_stacked.median(), med_purity_own=h.pur_own.median() if sr != "ceiling" else np.nan,
                             med_cc=h.cc_al.median(), frac_cc07=(h.cc_al > 0.7).mean(), own_best=h.best_is_own.astype(float).mean(),
                             top3=(h["rank"] <= 3).mean(), best_nbr3=h.best_is_nbr3.astype(float).mean()))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--max-pairs", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(OD, exist_ok=True)
    t0 = time.time()
    P = pd.read_csv(f"{T2D}/pairs.csv")
    if a.max_pairs:
        P = P.head(a.max_pairs)
    T2Q = pd.read_csv(f"{T2D}/round0_qc.csv")
    sel, SEL = select(P, T2Q)
    SEL.to_csv(f"{OD}/selection_counts.csv", index=False)
    pd.DataFrame([dict(pair=k[0], series=k[1], variant=k[2], t=t) for k, tt in sel.items() for t in tt]).to_parquet(
        f"{OD}/candidates.parquet", index=False)
    log(f"[t2b] selections {len(sel)} ({time.time()-t0:.0f}s); median n_cand by variant/series "
        f"{SEL.groupby(['series','variant']).n_cand.median().to_dict()}")
    jobs = []
    for (w, ns), g in P.groupby(["window", "netsta"]):
        items = [((int(pr), sr, v), sel[(int(pr), sr, v)]) for pr in g.pair for sr in ("data", "control") for v in VARIANTS]
        jobs.append((w, ns, items))
    STK = {}
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(stack_job, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            STK.update(f.result())
            if (i + 1) % 20 == 0:
                log(f"[t2b] stacked {i+1}/{len(jobs)} stations {time.time()-t0:.0f}s")
                open(f"{W.RES}/T2B_HEARTBEAT", "w").write(f"{time.time():.0f} stacked {i+1}/{len(jobs)}\n")
    log(f"[t2b] stacking done {time.time()-t0:.0f}s")
    TT = build_truth(P)
    F = W.load_families()
    la, lo = F.lat_use.to_numpy(float), F.lon_use.to_numpy(float)
    dk = W.km(la[:, None], lo[:, None], la[None], lo[None])
    fi = {f: i for i, f in enumerate(F.fid)}
    PU = purity(sel, P)
    SC = SEL.set_index(["pair", "series", "variant"])
    rows, stacks, keys = [], [], []
    for p in P.itertuples():
        for sr in ("data", "control"):
            for v in VARIANTS:
                k = (int(p.pair), sr, v)
                st, n, lv = STK.get(k, (None, 0, np.zeros(3, bool)))
                r = dict(pair=p.pair, fid=p.fid, family_uid=p.family_uid, region=p.region, window=p.window, netsta=p.netsta,
                         series=sr, variant=v, pool=int(SC.loc[k, "pool"]), N_t2=int(SC.loc[k, "N_t2"]),
                         n_cand=int(SC.loc[k, "n_cand"]), meanP=SC.loc[k, "meanP"], n_stacked=n, pur_own=PU[k])
                if st is not None:
                    r.update(compare(p, st, lv, TT, dk, fi))
                    keys.append(f"{p.pair}|{sr}|{v}"); stacks.append(st.astype(np.float32))
                rows.append(r)
    O = pd.DataFrame(rows)
    O.to_csv(f"{OD}/pair_results.csv", index=False)
    np.savez(f"{OD}/stacks.npz", keys=np.array(keys), stacks=np.stack(stacks) if stacks else np.zeros((0, 3, NL), np.float32))
    T2 = pd.read_csv(f"{T2D}/pair_results.csv")
    S = summarise(O, T2)
    S.to_csv(f"{OD}/summary.csv", index=False)
    log("[t2b] summary\n" + S.to_string())
    log(f"[t2b] done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
