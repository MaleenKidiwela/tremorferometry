#!/usr/bin/env python
"""Step 1 of the reproduce-Bostock loop: LOUD events + sub-second time correction.

For the N most productive KEPT families in the 2005 bed:
  a. loudness = network-median 3-s S-energy excess over the local tremor background at the 4
     nearest detection sites (computed, never magnitude);
  b. sub-second correction of the whole-second 2005 catalogue times by matched filtering each
     event against an iterated stack of the family's loudest events (Bostock's own bootstrap),
     with the ALIGNMENT stations and the VERIFICATION stations disjoint so that the coherence
     check is not the selection statistic;
  c. stack coherence, uncorrected vs corrected, at the verification sites, against a +300 s
     shifted-catalogue null.

  python scripts/k_s1_loud.py --nfam 20 --out results/repro/s1
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

FS = K.FS
PRE = 2.0          # template starts this far before the predicted S
LEN = 8.0          # template length (s)
LAGMAX = 1.5       # search half-width for the time correction (s)
BG_HALF = 60.0     # background half-window for the energy excess (s)
S_WIN = 3.0        # S energy window (s)


def energy_excess(seg_s, seg_bg):
    """E(3-s S window) / median E(3-s windows of the background), both 3-C mean square."""
    ns = int(round(S_WIN * FS))
    es = float((seg_s ** 2).mean())
    x = (seg_bg ** 2).sum(0)
    c = np.concatenate([[0.0], np.cumsum(x)])
    w = (c[ns:] - c[:-ns]) / (3.0 * ns)
    w = w[w > 0]
    if w.size < 10 or es <= 0:
        return np.nan
    return es / float(np.median(w))


def xcorr_lags(ref, seg, lagmax):
    """Normalised 3-C cc of `ref` (3,n) against `seg` (3, n+2*lagmax) at every integer lag."""
    n = ref.shape[1]
    L = 2 * lagmax + 1
    rn = np.sqrt((ref ** 2).sum())
    if rn <= 0:
        return None
    num = np.zeros(L)
    den = np.zeros(L)
    for k in range(L):
        w = seg[:, k:k + n]
        num[k] = float((ref * w).sum())
        den[k] = float(np.sqrt((w ** 2).sum()))
    ok = den > 0
    cc = np.zeros(L)
    cc[ok] = num[ok] / (rn * den[ok])
    return cc


def cut_event(dc, sta, t, pad_s=0.0):
    return dc.cut(sta, t - PRE - pad_s, LEN + 2 * pad_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfam", type=int, default=20)
    ap.add_argument("--out", default=f"{K.OUT}/s1")
    ap.add_argument("--nalign", type=int, default=120, help="loudest events used to build the ref")
    ap.add_argument("--iters", type=int, default=3)
    ap.add_argument("--shift", type=float, default=0.0, help="catalogue shift (s) for the null")
    ap.add_argument("--stas", default="")
    ap.add_argument("--hidden", action="store_true",
                    help="THE SEALED STRATUM -- only for the single blind check")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    B = K.bed_catalogue(kept_only=(False if a.hidden else True))
    if a.hidden:
        import s_common as _S
        hid, _, _ = _S.holdouts()
        B = B[B.bkey.isin(hid)].reset_index(drop=True)
    if a.shift:
        B = B.copy()
        B["t_origin"] = B.t_origin + a.shift
        B["date"] = pd.to_datetime(B.t_origin, unit="s").dt.strftime("%Y-%m-%d")
        B = B[B.date.isin(set(K.bed_days()))].reset_index(drop=True)
    cnt = B.groupby("bkey").size().sort_values(ascending=False)
    fams = [int(b) for b in cnt.index[:a.nfam]]
    DET = a.stas.split(",") if a.stas else K.DET_SITES
    TT = K.fam_tt(DET)
    near = {b: K.nearest_stas(b, DET, 8) for b in fams}
    ALIGN = {b: [near[b][i] for i in (0, 2, 4, 6)] for b in fams}
    VERIF = {b: [near[b][i] for i in (1, 3, 5, 7)] for b in fams}
    LOUD4 = {b: near[b][:4] for b in fams}
    need = sorted({s for b in fams for s in near[b]})
    K.log(f"{len(fams)} families, stations {need}")

    B = B[B.bkey.isin(fams)].reset_index(drop=True)
    dc = K.DayCache(need)

    # ---------------------------------------------------------------- a. loudness
    rows = []
    for date, g in B.groupby("date"):
        dc.day(date)
        for r in g.itertuples():
            b = int(r.bkey)
            ex = []
            for s in LOUD4[b]:
                tS = r.t_origin + TT[(b, s)]
                ss = dc.cut(s, tS - 0.5, S_WIN, date)
                bg = dc.cut(s, tS - BG_HALF, 2 * BG_HALF, date)
                if ss is None or bg is None:
                    continue
                v = energy_excess(ss, bg)
                if np.isfinite(v):
                    ex.append(v)
            rows.append(dict(bkey=b, t_origin=float(r.t_origin), date=date,
                             n_sta=len(ex), excess=(float(np.median(ex)) if ex else np.nan)))
        K.log(f"[loud] {date} {len(g)} events")
    L = pd.DataFrame(rows)
    L["loud"] = False
    for b in fams:
        m = (L.bkey == b) & np.isfinite(L.excess)
        if m.sum() >= 10:
            thr = float(np.nanquantile(L.excess[m], 0.90))
            L.loc[m & (L.excess >= thr), "loud"] = True
    L.to_csv(f"{a.out}_loud.csv", index=False)
    K.log(f"[loud] {int(L.loud.sum())} loud of {len(L)} events")

    # ---------------------------------------------------------------- b. reference bootstrap
    lm = int(round(LAGMAX * FS))
    n = int(round(LEN * FS))
    refs, refinfo = {}, []
    for b in fams:
        sel = L[(L.bkey == b) & L.loud].sort_values("excess", ascending=False).head(a.nalign)
        if len(sel) < 10:
            continue
        stas = ALIGN[b]
        segs = {s: [] for s in stas}          # padded segments, index-aligned with `keep`
        keep = []
        for date, g in sel.groupby("date"):
            dc.day(date)
            for r in g.itertuples():
                cut = {}
                for s in stas:
                    c = dc.cut(s, r.t_origin + TT[(b, s)] - PRE - LAGMAX,
                               LEN + 2 * LAGMAX, date)
                    if c is not None and (c ** 2).sum() > 0:
                        cut[s] = K.norm3(c)
                if len(cut) >= 2:
                    keep.append((float(r.t_origin), float(r.excess)))
                    for s in stas:
                        segs[s].append(cut.get(s))
        if len(keep) < 10:
            continue
        N = len(keep)
        order = int(np.argmax([k[1] for k in keep]))
        ref = {s: (segs[s][order][:, lm:lm + n].copy() if segs[s][order] is not None else None)
               for s in stas}
        dts = np.zeros(N)
        for it in range(a.iters):
            cc_net = np.full(N, np.nan)
            best = np.zeros(N, int)
            for i in range(N):
                acc, cnt = np.zeros(2 * lm + 1), 0
                for s in stas:
                    if segs[s][i] is None or ref[s] is None:
                        continue
                    c = xcorr_lags(ref[s], segs[s][i], lm)
                    if c is not None:
                        acc += c
                        cnt += 1
                if cnt:
                    acc /= cnt
                    best[i] = int(np.argmax(acc)) - lm
                    cc_net[i] = float(acc.max())
            dts = best / FS
            thr = np.nanquantile(cc_net, 0.5) if it < a.iters - 1 else -np.inf
            use = np.isfinite(cc_net) & (cc_net >= thr)
            newref = {}
            for s in stas:
                acc, cnt = np.zeros((3, n)), 0
                for i in np.nonzero(use)[0]:
                    if segs[s][i] is None:
                        continue
                    k0 = lm + best[i]
                    acc += segs[s][i][:, k0:k0 + n]
                    cnt += 1
                newref[s] = K.norm3(acc / cnt) if cnt >= 5 else None
            ref = newref
        refs[b] = dict(ref=ref, stas=stas, t=[k[0] for k in keep], dt=dts.tolist(),
                       cc=cc_net.tolist())
        refinfo.append(dict(bkey=b, n_ref=N, cc_med=float(np.nanmedian(cc_net)),
                            dt_std=float(np.nanstd(dts)),
                            dt_absmean=float(np.nanmean(np.abs(dts)))))
        K.log(f"[ref] fam {b}: N={N} cc_med={np.nanmedian(cc_net):.3f} "
              f"dt_std={np.nanstd(dts):.3f}s")
    pd.DataFrame(refinfo).to_csv(f"{a.out}_ref.csv", index=False)

    # ---------------------------------------------------------------- c. correct every event
    corr, acc = [], {}
    for b in refs:
        for s in near[b]:
            acc[(b, s, "raw")] = [np.zeros((3, n)), 0]
            acc[(b, s, "cor")] = [np.zeros((3, n)), 0]
    for date, g in B[B.bkey.isin(list(refs))].groupby("date"):
        dc.day(date)
        for r in g.itertuples():
            b = int(r.bkey)
            R = refs[b]
            segs, cutraw = {}, {}
            for s in near[b]:
                c = dc.cut(s, r.t_origin + TT[(b, s)] - PRE - LAGMAX, LEN + 2 * LAGMAX, date)
                if c is None or (c ** 2).sum() <= 0:
                    continue
                c = K.norm3(c)
                segs[s] = c
                cutraw[s] = c[:, lm:lm + n]
            accv, cnt = np.zeros(2 * lm + 1), 0
            for s in R["stas"]:
                if s in segs and R["ref"][s] is not None:
                    cc = xcorr_lags(R["ref"][s], segs[s], lm)
                    if cc is not None:
                        accv += cc
                        cnt += 1
            if not cnt:
                continue
            accv /= cnt
            k = int(np.argmax(accv)) - lm
            corr.append(dict(bkey=b, t_origin=float(r.t_origin), date=date,
                             dt=k / FS, cc=float(accv.max()), n_align=cnt,
                             loud=bool(L.loud[(L.bkey == b) &
                                              (L.t_origin == float(r.t_origin))].any())))
            for s in near[b]:
                if s not in segs:
                    continue
                A = acc[(b, s, "raw")]
                A[0] += K.norm3(cutraw[s])
                A[1] += 1
                A = acc[(b, s, "cor")]
                A[0] += K.norm3(segs[s][:, lm + k:lm + k + n])
                A[1] += 1
        K.log(f"[corr] {date} {len(g)}")
    C = pd.DataFrame(corr)
    C.to_csv(f"{a.out}_times.csv", index=False)

    out = []
    for b in refs:
        for s in near[b]:
            for tag in ("raw", "cor"):
                Sm, N = acc[(b, s, tag)]
                if N < 5:
                    continue
                out.append(dict(bkey=b, sta=s, kind=tag, n=N,
                                coh=float(N * ((Sm / N) ** 2).mean()),
                                role=("align" if s in ALIGN[b] else "verify"),
                                dist_km=K.epi_km(b, s)))
    O = pd.DataFrame(out)
    O.to_csv(f"{a.out}_coherence.csv", index=False)
    piv = O.pivot_table(index=["role", "kind"], values="coh", aggfunc=["median", "mean", "max"])
    K.log("\n" + str(piv))
    js = dict(n_fam=len(refs), n_events=len(C),
              dt_std=float(C.dt.std()), dt_absmed=float(C.dt.abs().median()),
              cc_med=float(C.cc.median()),
              coh=dict((f"{r}_{k}", float(O[(O.role == r) & (O.kind == k)].coh.median()))
                       for r in ("align", "verify") for k in ("raw", "cor")),
              shift=a.shift)
    json.dump(js, open(f"{a.out}_summary.json", "w"), indent=1)
    K.log(json.dumps(js, indent=1))


if __name__ == "__main__":
    main()
