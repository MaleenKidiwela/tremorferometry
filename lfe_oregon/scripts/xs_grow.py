#!/usr/bin/env python
"""Candidate templates -> iterative network cross-correlation + stacking (Thomas & Bostock 2015 sec.2).

Paper (T&B 2015 sec.2), step by step, and where each lives here:
  "Once candidate templates are identified, we use network cross-correlation on 3 days of data to
   register additional detections defined as times when the network cross-correlation exceeds
   eight times the median absolute deviation."               -> ab_grow.matched_filter, nmad 8
  "Detections are then stacked to create a new template"      -> restack() below
  "phase-weighted stacking on the intermediate iterations (Schimmel and Paulssen, 1997)"
                                                              -> PWS, nu = --pws-nu (S&P use 2)
  "iteratively apply cross-correlation and stacking until no additional detections are registered
   (generally 10-15 iterations)"                              -> per-template convergence
  "If additional stations begin to show coherent energy, those stations are included"
                                                              -> admit_stations() (rule SILENT)
  "We require a minimum of 100 detections in the 3-day window to stack an LFE template."
                                                              -> column pass100 (no other gate)
Savard & Bostock 2015 Discussion: "Using our cross-station detection time windows as initial
templates" -> each cross-station detection is one candidate; its template is its own 18 s detection
window (xs_detect TLEN) on the detection's stations, one common absolute window (moveout kept).

SILENT choices (defaults here, all CLI columns; see METHOD_FROM_PAPER.md sec.5):
  --tpl-len 18 (the detection window) ; --comps ZNE ; --net mean (G&R channel-mean) ;
  --min-sta-valid 3 ; --decl-s = tpl-len ; --pws-nu 2 ; --max-iters 15 ;
  final template = LINEAR stack of the converged detections (PWS only "intermediate") ;
  station admission: odd/even half-stack cc in the template window above the station's own
  random-time null (--admit-q = 99.9th pct of --null-draws half-stack pairs; for incoherent
  noise the half-stack cc null does not depend on the stack size, so --null-half 10 suffices) -- a data-driven reading of
  "coherent energy"; the statistic is written as a column for every pool station every iteration.

The 10 %/5 s duplicate rule is applied AFTER the 46-day scan (xs_scan46.py), as the paper orders it;
--merge-each-iter applies it after every iteration instead (ablation; ab_grow.merge_duplicates).

  python scripts/xs_grow.py --out results/grow [--dates 2009-09-03,...] [--max-cand 0] [--iters 15]
-> <out>/det_it<k>.parquet, <out>/tpl_it<k>.npz, <out>/state_it<k>.csv, <out>/stations_it<k>.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time

import numpy as np
import pandas as pd
from scipy.signal import hilbert

import xs_common as C

A = C.patch_ab_common()
import ab_grow as G                                                          # noqa: E402

LONG_PRE, LONG_LEN = 12.0, 44.0          # stack window around the template start (ab_grow's)
POOL = [s for s in C.STATIONS.netsta if os.path.isdir(f"{C.CACHE}/{s}")]


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------ candidates ----------------
def load_candidates(xs_dir, dates, max_cand=0):
    D = pd.concat([pd.read_parquet(f"{xs_dir}/{d}_det.parquet").assign(date=d) for d in dates
                   if os.path.exists(f"{xs_dir}/{d}_det.parquet")], ignore_index=True)
    D = D.sort_values("t_win").reset_index(drop=True)
    D["grp"] = np.arange(len(D))
    if max_cand:
        D = D.sort_values("maxwin", ascending=False).head(max_cand).sort_values("t_win")
    return D


def initial_templates(D, dc, tpl_len, comps):
    """{g: {sta: (tpl (3,nt), off 0.0, snr nan)}}: the detection's own window on its stations."""
    TPL, rows = {}, []
    tcols = [c for c in D.columns if c.startswith("tS_")]
    for i, r in zip(D.index, D.itertuples()):
        stas = [c[3:] for c in tcols if np.isfinite(D.at[i, c])]
        T = {}
        for s in stas:
            seg = dc.cut(s, r.t_win, tpl_len, r.date)
            if seg is None:
                continue
            seg = seg * comps[:, None]
            if (seg ** 2).sum() > 0:
                T[s] = (seg, 0.0, np.nan)
        if len(T) >= 1:
            TPL[int(r.grp)] = T
        rows.append(dict(grp=int(r.grp), n_sta0=len(T), stations0=",".join(sorted(T))))
    return TPL, pd.DataFrame(rows)


# ------------------------------------------------------------------ stacking ------------------
def unit_max3(seg):
    m = float(np.abs(seg).max())
    return seg / m if m > 0 else seg


def restack(dc, det, stas, dates):
    """Stream the detections once: per (grp, sta) linear sum of unit-max 3-C segments, the phasor
    sum (for PWS) and odd/even half sums (for the admission statistic)."""
    nL = int(round(LONG_LEN * C.FS))
    acc = {}
    rc = A.RefCache(dc, nL)
    for date, dd in det.sort_values("t").groupby("date"):
        if date not in dates:
            continue
        refs = {s: rc.get(s, date) for s in stas}
        for r in dd.itertuples():
            g = int(r.grp)
            for s in stas:
                seg = dc.cut(s, r.t - LONG_PRE, LONG_LEN, date)
                if seg is None:
                    continue
                okc = refs[s] > 0
                if not okc.any():
                    continue
                seg = np.where(okc[:, None], seg, 0.0)[:, :nL]
                u = unit_max3(seg)
                an = hilbert(u, axis=1)
                ph = an / np.where(np.abs(an) > 0, np.abs(an), 1.0)
                k = (g, s)
                if k not in acc:
                    acc[k] = dict(lin=np.zeros((3, nL)), ph=np.zeros((3, nL), complex),
                                  h=[np.zeros((3, nL)), np.zeros((3, nL))], n=0)
                a = acc[k]
                a["lin"] += u
                a["ph"] += ph
                a["h"][a["n"] % 2] += u
                a["n"] += 1
    return acc


def pws(a, nu):
    lin = a["lin"] / a["n"]
    c = np.abs(a["ph"] / a["n"]) ** nu
    return lin * c


def half_cc(a, i0, nt):
    x, y = a["h"][0][:, i0:i0 + nt].ravel(), a["h"][1][:, i0:i0 + nt].ravel()
    d = np.sqrt((x ** 2).sum() * (y ** 2).sum())
    return float((x * y).sum() / d) if d > 0 else np.nan


def env_snr(st, i0, nt):
    e = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    med = float(np.median(e))
    return float(e[i0:i0 + nt].max() / med) if med > 0 else np.nan


def station_null(dc, stas, dates, tpl_len, draws, nhalf, q=99.9, seed=0):
    """Per station: `draws` values of the odd/even half-stack cc for random start times
    (nhalf windows per half) -> 99th percentile.  Same statistic as admit_stations.  Date-outer
    (the DayCache holds one date); windows touching an invalid sample are skipped."""
    rng = np.random.default_rng(seed)
    nt = int(round(tpl_len * C.FS))
    H = {s: np.zeros((draws, 2, 3, nt)) for s in stas}
    kd = rng.integers(len(dates), size=(draws, 2 * nhalf))
    for di, d in enumerate(dates):
        Dd = dc.day(d)
        dr, kk = np.nonzero(kd == di)
        for s in stas:
            if s not in Dd:
                continue
            t0, X, m = Dd[s]
            i = rng.integers(int(60 * C.FS), X.shape[1] - nt - int(60 * C.FS), size=dr.size)
            idx = i[:, None] + np.arange(nt)[None, :]
            ok = m[idx].all(1)
            seg = X[:, idx].transpose(1, 0, 2).astype(np.float64)             # (k,3,nt)
            mx = np.abs(seg).reshape(len(seg), -1).max(1)
            ok &= mx > 0
            seg = seg[ok] / mx[ok, None, None]
            np.add.at(H[s], (dr[ok], kk[ok] % 2), seg)
    out = {}
    for s in stas:
        x = H[s][:, 0].reshape(draws, -1)
        y = H[s][:, 1].reshape(draws, -1)
        dd = np.sqrt((x ** 2).sum(1) * (y ** 2).sum(1))
        v = np.where(dd > 0, (x * y).sum(1) / np.where(dd > 0, dd, 1), np.nan)
        out[s] = float(np.nanpercentile(v, q)) if np.isfinite(v).any() else np.inf
    return out


# ------------------------------------------------------------------ main loop -----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xs", default=f"{C.RES}/xs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dates", default=",".join(C.DETECT_DAYS))
    ap.add_argument("--max-cand", type=int, default=0, help="0 = every cross-station detection")
    ap.add_argument("--iters", type=int, default=15)
    ap.add_argument("--tpl-len", type=float, default=18.0)
    ap.add_argument("--comps", default="ZNE")
    ap.add_argument("--net", default="mean", choices=["mean", "sum"])
    ap.add_argument("--min-sta-valid", type=int, default=3)
    ap.add_argument("--nmad", type=float, default=8.0)
    ap.add_argument("--pws-nu", type=float, default=2.0)
    ap.add_argument("--min-det", type=int, default=100)
    ap.add_argument("--null-draws", type=int, default=2000)
    ap.add_argument("--null-half", type=int, default=10)
    ap.add_argument("--admit-q", type=float, default=99.9, help="null percentile a station must exceed")
    ap.add_argument("--no-admit", action="store_true", help="ablation: never add stations")
    ap.add_argument("--merge-each-iter", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    json.dump(vars(a), open(f"{a.out}/args.json", "w"), indent=1)
    dates = a.dates.split(",")
    comps = np.array([c in a.comps for c in "ZNE"], float)
    nt = int(round(a.tpl_len * C.FS))
    i0 = int(round(LONG_PRE * C.FS))
    dc = A.DayCache(POOL)
    D = load_candidates(a.xs, dates, a.max_cand)
    TPL, info = initial_templates(D, dc, a.tpl_len, comps)
    info.to_csv(f"{a.out}/candidates.csv", index=False)
    log(f"{len(D)} cross-station detections -> {len(TPL)} candidate templates, pool {len(POOL)}")
    nullp = station_null(dc, POOL, dates, a.tpl_len, a.null_draws, a.null_half, a.admit_q)
    pd.Series(nullp, name=f"half_cc_null{a.admit_q:g}").to_csv(f"{a.out}/station_null.csv")
    state = pd.DataFrame(dict(grp=sorted(TPL), n_prev=0, converged=False, it_conv=-1))
    state = state.set_index("grp")
    for it in range(1, a.iters + 1):
        active = [g for g in TPL if not state.at[g, "converged"]]
        if not active:
            log(f"it{it}: every template converged")
            break
        t0 = time.time()
        det = G.matched_filter(dc, {g: TPL[g] for g in active}, dates, nmad=a.nmad, net=a.net,
                               min_sta_valid=a.min_sta_valid, decl_s=a.tpl_len, tpl_len=a.tpl_len)
        if a.merge_each_iter and len(det):
            det, _, mst = G.merge_duplicates(det)
            log(f"  merge: {mst}")
        det.to_parquet(f"{a.out}/det_it{it}.parquet", index=False)
        n = det.groupby("grp").size().reindex(active, fill_value=0)
        acc = restack(dc, det, POOL, dates)
        rows, newT = [], {}
        for g in active:
            T = TPL[g]
            members = set(T)
            final = False
            conv = n[g] <= state.at[g, "n_prev"]                     # "no additional detections"
            last = conv or it == a.iters
            for s in POOL:
                k = (g, s)
                if k not in acc or acc[k]["n"] < 2:
                    continue
                hc = half_cc(acc[k], i0, nt)
                lin = acc[k]["lin"] / acc[k]["n"]
                snr = env_snr(lin, i0, nt)
                admit = (s not in members) and (not a.no_admit) and np.isfinite(hc) and hc > nullp[s]
                rows.append(dict(it=it, grp=g, sta=s, n=acc[k]["n"], half_cc=hc, null_q=nullp[s],
                                 snr=snr, member=s in members, admitted=admit))
                if s in members or admit:
                    st_ = lin if last else pws(acc[k], a.pws_nu)          # PWS on intermediate its
                    seg = st_[:, i0:i0 + nt] * comps[:, None]
                    if (seg ** 2).sum() > 0:
                        newT.setdefault(g, {})[s] = (seg, 0.0, snr)
            if conv and n[g] > 0:
                state.at[g, "converged"] = True
                state.at[g, "it_conv"] = it
            if n[g] == 0:
                state.at[g, "converged"] = True
                state.at[g, "it_conv"] = it
            state.at[g, "n_prev"] = int(n[g])
            if g in newT and not state.at[g, "converged"]:
                TPL[g] = newT[g]
            elif g in newT and state.at[g, "converged"]:
                TPL[g] = newT[g]                                        # final linear template
        S = pd.DataFrame(rows)
        S.to_csv(f"{a.out}/stations_it{it}.csv", index=False)
        st = state.copy()
        st["n_sta"] = [len(TPL[g]) for g in st.index]
        st["stations"] = [",".join(sorted(TPL[g])) for g in st.index]
        st["pass100"] = st.n_prev >= a.min_det
        st.to_csv(f"{a.out}/state_it{it}.csv")
        np.savez_compressed(f"{a.out}/tpl_it{it}.npz", **{f"{g}|{s}": v[0].astype(np.float32)
                                                          for g in TPL for s, v in TPL[g].items()})
        log(f"it{it}: {len(active)} active, {len(det)} det, median n {float(np.median(n)):.0f}, "
            f"converged {int(state.converged.sum())}/{len(state)}, pass100 {int(st.pass100.sum())}, "
            f"admitted {int(S.admitted.sum()) if len(S) else 0}, {time.time()-t0:.0f}s")
    open(f"{a.out}/GROW_DONE", "w").write(time.strftime("%FT%T\n"))


if __name__ == "__main__":
    main()
