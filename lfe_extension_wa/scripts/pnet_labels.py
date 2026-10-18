#!/usr/bin/env python
"""Stage B dataset: 60-s ZNE windows on the WA permanent stations for fine-tuning PNET (BUILD half only).

Scope: V2008 + V2010 days, windows lying entirely inside BUILD-half 2-h blocks (wa_common.half_of == 0 at both ends,
same block, outside the 60-s guards); the TEST half, A2011 and the held-out stations (pnet_common.HELD_OUT) are
never touched.
Label rule (pnet_common.label_pairs): (family, station node) pairs whose build-half stack (>= 30 events) has a
clear-ish arrival (imp >= 3 and dur < 3 s); weight 1.0 if clear (imp > 5), 0.5 otherwise.  S time = event time
t + offset_s (the template start = stack envelope peak - 0.5 s, the onset).  Events = every build-half array
detection of the family (truth.parquet), whether or not it is individually visible (the SVI recipe).
  positive windows  one per labelled arrival (cap per station-day), S placed uniformly in [5, 55] s; target =
                    Gaussian sigma 0.5 s at EVERY labelled arrival inside the window; loss weight = the label weight
                    within +-3 s of a labelled arrival, 0 elsewhere (the rest of an episode window is unlabelled).
  quiet negatives   random windows on quiet days (V2010 'quiet', V2008 'post'; tremor-free by the window design),
                    target 0, weight 1, except weight 0 within +-3 s of ANY known arrival (all >= 30-event templates
                    of any family at that station node, any clarity).
  control negatives the station's window at the position where a SITE-ROLLED scrambled-control candidate (step-1
                    series c1..c4, build half) read it, target 0, weight 1, same +-3 s mask around known arrivals.
Stored peak-normalised float16.  -> results/pnet/train/{X.npy, Y.npy, Wt.npy, meta.parquet}
"""
from __future__ import annotations
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

NW = 1200
SIG = 10.0          # 0.5 s
MASK = 60           # +-3 s
OD = f"{PC.OUT}/train"
_G = {}


def _init():
    import wa_mf
    T, Q = wa_mf.load_tpl(["build"], all_templates=True)
    L = PC.label_pairs(pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv"))
    S = W.load_stations()
    TR = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    _G.update(TOFF={k: v["off"] for k, v in T.items()},
              LAB={(int(r.fid), r.netsta, r.node): float(r.lab_w) for r in L.itertuples()},
              SITE=dict(zip(S.netsta, S.site)), TR={d: g for d, g in TR.groupby("date")})


def build_ok(s_abs):
    """window [s, s+60] entirely in one build-half block"""
    s_abs = np.atleast_1d(np.asarray(s_abs, float))
    a = W.half_of(s_abs); b = W.half_of(s_abs + 60.0)
    same = np.floor(s_abs / W.BLOCK_S) == np.floor((s_abs + 60.0) / W.BLOCK_S)
    return (a == 0) & (b == 0) & same


def job(args):
    ns, ds, kind, cap_pos, n_quiet, cap_ctl, seed = args
    TOFF, LAB, SITE = _G["TOFF"], _G["LAB"], _G["SITE"]
    z = PC.load_zne(ns, ds)
    if z is None:
        return None
    x, node = z["x"], z["node"]
    D0 = pd.Timestamp(ds).timestamp()
    rng = np.random.default_rng(seed)
    tg = _G["TR"].get(ds)
    lab_t, lab_w, msk_t = [], [], []
    if tg is not None:
        for f, g in tg.groupby("fid"):
            k = (int(f), ns, node)
            if k not in TOFF:
                continue
            s = (g.t.to_numpy() - D0) * W.FS + TOFF[k] * W.FS
            msk_t.append(s)
            if k in LAB:
                lab_t.append(s); lab_w.append(np.full(len(s), LAB[k]))
    lab_t = np.concatenate(lab_t) if lab_t else np.zeros(0)
    lab_w = np.concatenate(lab_w) if lab_w else np.zeros(0)
    msk_t = np.sort(np.concatenate(msk_t)) if msk_t else np.zeros(0)
    o = np.argsort(lab_t); lab_t, lab_w = lab_t[o], lab_w[o]
    valid = z["m"] & (np.abs(x).sum(0) > 0)
    cv = np.concatenate([[0], np.cumsum(~valid)])
    out_x, out_y, out_w, out_k = [], [], [], []

    def add(s0, kindc):
        s0 = int(s0)
        if s0 < 0 or s0 + NW > W.NDAY or not build_ok(D0 + s0 / W.FS)[0]:
            return False
        if cv[s0 + NW] - cv[s0] > NW // 10:            # > 10 % invalid
            return False
        seg = x[:, s0:s0 + NW]
        mx = np.abs(seg).max()
        if not mx > 0:
            return False
        y = np.zeros(NW, np.float32); w = np.zeros(NW, np.float32)
        tt = np.arange(NW)
        lo, hi = np.searchsorted(lab_t, s0 - MASK), np.searchsorted(lab_t, s0 + NW + MASK)
        if kindc == 0:
            for tl, wl in zip(lab_t[lo:hi], lab_w[lo:hi]):
                r = tl - s0
                y = np.maximum(y, np.exp(-0.5 * ((tt - r) / SIG) ** 2))
                m = np.abs(tt - r) <= MASK
                w[m] = np.maximum(w[m], wl)
        else:
            w[:] = 1.0
            ml, mh = np.searchsorted(msk_t, s0 - MASK), np.searchsorted(msk_t, s0 + NW + MASK)
            for tm in msk_t[ml:mh]:
                w[np.abs(tt - (tm - s0)) <= MASK] = 0.0
        w[~valid[s0:s0 + NW]] = 0.0
        out_x.append((seg / mx).astype(np.float16)); out_y.append((y * 255).round().astype(np.uint8))
        out_w.append((w * 255).round().astype(np.uint8)); out_k.append(kindc)
        return True

    # positives
    if len(lab_t):
        ip = rng.permutation(len(lab_t))[: 4 * cap_pos]
        npos = 0
        for i in ip:
            if npos >= cap_pos:
                break
            if add(lab_t[i] - rng.uniform(5, 55) * W.FS, 0):
                npos += 1
    # quiet negatives
    if kind in ("quiet", "post"):
        k = 0
        for s0 in rng.integers(0, W.NDAY - NW, 6 * n_quiet):
            if k >= n_quiet:
                break
            if add(s0, 1):
                k += 1
    # scrambled-control negatives
    fp = f"{W.RES}/step1_prob/feat/feat_{ds}.parquet"
    if cap_ctl and os.path.exists(fp):
        F = pd.read_parquet(fp, columns=["fid", "t", "series"])
        F = F[F.series != "data"]
        J1 = json.load(open(f"{W.RES}/mf_mt1_all/run_{ds}.json"))
        SH = {"c1": J1["shift_site"], "c2": J1["shift_site2"]}
        p2 = f"{W.RES}/mf_step1c34/run_{ds}.json"
        if os.path.exists(p2):
            J2 = json.load(open(p2)); SH["c3"] = J2["shift_site"]; SH["c4"] = J2["shift_site2"]
        F = F[F.series.isin(SH)]
        F = F[[(int(f), ns, node) in TOFF for f in F.fid]]
        if len(F):
            NV = W.NDAY - W.NT + 1
            F = F.iloc[rng.permutation(len(F))[: 4 * cap_ctl]]
            i = np.round((F.t.to_numpy() - D0) * W.FS).astype(np.int64)
            sh = np.array([int(SH[s].get(SITE[ns], 0)) % NV for s in F.series], np.int64)
            j = ((i - sh) % NV) + np.array([int(round(TOFF[(int(f), ns, node)] * W.FS)) for f in F.fid])
            k = 0
            for jj in j:
                if k >= cap_ctl:
                    break
                if add(jj - rng.uniform(5, 55) * W.FS, 2):
                    k += 1
    if not out_x:
        return ns, ds, None
    return ns, ds, (np.stack(out_x), np.stack(out_y), np.stack(out_w), np.array(out_k, np.int8), z["n_horiz"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--cap-pos", type=int, default=500)
    ap.add_argument("--n-quiet", type=int, default=250)
    ap.add_argument("--cap-ctl", type=int, default=100)
    ap.add_argument("--days", default="")
    ap.add_argument("--out", default=OD)
    a = ap.parse_args()
    S = W.load_stations()
    stas = [s for s in S.netsta if s not in PC.HELD_OUT]
    days = [(d, k) for w in ["V2008", "V2010"] for d, k in W.window_days(w)]
    if a.days:
        days = [x for x in days if x[0] in set(a.days.split(","))]
    jobs = [(ns, d, k, a.cap_pos, a.n_quiet, a.cap_ctl, 7 + i) for i, (d, k) in enumerate(days) for ns in stas
            if os.path.exists(W.cache_path(ns, d))]
    W.log(PC.LOG, f"[labels] {len(jobs)} station-days on {len(stas)} training stations; held out {PC.HELD_OUT}")
    X, Y, Wt, meta = [], [], [], []
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = [ex.submit(job, j) for j in jobs]
        for i, fu in enumerate(as_completed(fut)):
            r = fu.result()
            if r is None or r[2] is None:
                continue
            ns, ds, (x, y, w, k, nh) = r
            X.append(x); Y.append(y); Wt.append(w)
            meta.append(pd.DataFrame(dict(netsta=ns, date=ds, kind=k, n_horiz=nh)))
            if (i + 1) % 300 == 0:
                W.log(PC.LOG, f"[labels] {i+1}/{len(jobs)} {time.time()-t0:.0f}s")
    os.makedirs(a.out, exist_ok=True)
    M = pd.concat(meta, ignore_index=True)
    np.save(f"{a.out}/X.npy", np.concatenate(X)); np.save(f"{a.out}/Y.npy", np.concatenate(Y))
    np.save(f"{a.out}/Wt.npy", np.concatenate(Wt))
    M.to_parquet(f"{a.out}/meta.parquet", index=False)
    s = M.groupby(["kind", "n_horiz"]).size().to_dict()
    W.log(PC.LOG, f"[labels] {len(M)} windows (kind 0 pos / 1 quiet / 2 ctl, n_horiz): {s}; "
          f"per station positives {M[M.kind == 0].netsta.value_counts().head(12).to_dict()}; {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
