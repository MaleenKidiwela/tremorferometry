#!/usr/bin/env python
"""IS THE HANDOFF'S ENERGY-ONSET PICKER THE BOTTLENECK?  (oracle post-mortem, 2026-09-20)

`w_grow.handoff` cuts the pass-0 template at the argmax of the 3-C stack envelope inside a 20-s
search window (HANDOFF_ONSET_SEARCH = -4 .. +16 s around the cluster time).  On this bed a 30-40
event stack shows its arrival at only ~3.4 background MAD, so that argmax can land on a tremor bump:
the oracle's own handoff diagnostic has onset_s spread over -4 .. +10 s (std 3.8 s).

This script builds the SAME oracle templates under four cut rules and measures, for each,
self-detection (does the template detect the events it was built from?) and one-to-one recall /
purity against the family's catalogued events over three ETS days:

  blind      the frozen rule: argmax of the envelope in -4 .. +16 s        (what 1c' ran)
  narrow     argmax of the envelope within +- 4 s of the PREDICTED S       (oracle-only)
  ttcut      no search at all: cut at predicted S - TPL0_PRE_S             (oracle-only)
  peakcut    argmax within +- 10 s of predicted S, on a per-station stack  (oracle-only)

If `blind` is far below the others, the onset picker -- not the bed -- is what returned C_30 = 0.
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
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import z_common as Z                                                        # noqa: E402
import z_diag_oracle as ZD                                                  # noqa: E402
import s_common as S                                                        # noqa: E402

FS = W.FS
RULES = ["blind", "narrow", "ttcut", "peakcut"]


def build(rule, bkeys, FT, det, of, out_dir, align_iters=2, align_lag=30):
    """Oracle handoff with one cut rule -> (bank path, seeds df, members df)."""
    import v_common as V
    import w_grow as G
    V_, E = W.patch_v()
    loader = E._W_LOADER
    n_cut = int(round(W.HANDOFF_CUT_LEN_S * FS))
    nom = int(round(W.HANDOFF_CUT_PRE_S * FS))
    n = int(round(W.TPL0_LEN_S * FS))
    FH = ZD.fam_hours(bkeys)
    rows, mems, diag = [], [], []
    cids, Wb, off = [], [], []
    for bk, (d, h, ev) in FH.items():
        tt = {s: FT[(bk, s)] for s in det if (bk, s) in FT}
        med = float(np.median(list(tt.values())))
        t = ev + med - 1.0
        cid = f"H{bk:04d}"
        rows.append(dict(cid=cid, date=d, n_members=len(t), n_in_hour=len(t), n_in_3h=len(t),
                         n_bins=1, t0=float(t.min()), t1=float(t.max()),
                         t_seed=float(np.median(t)), bkey=int(bk), hour=int(h),
                         axis_shift_s=float(med - 1.0)))
        mems.append(pd.DataFrame(dict(cid=cid, t=t)))
        Wrow = np.full((len(det), 3, n), np.nan, np.float32)
        orow = np.full(len(det), np.nan)
        for si, s in enumerate(det):
            if s not in tt:
                continue
            segs = [w for w in (G._cut(loader, s, x - W.HANDOFF_CUT_PRE_S, n_cut) for x in t)
                    if w is not None]
            if len(segs) < W.HANDOFF_MIN_SEGMENTS:
                continue
            Sg = np.stack(segs)
            st = (Sg.mean(0) if align_iters == 0
                  else S.align_stack(Sg, iters=align_iters, lag=align_lag)[0])
            # index of the PREDICTED S inside the cut
            i_pred = int(round((W.HANDOFF_CUT_PRE_S + tt[s] - med + 1.0) * FS))
            env = S.envelope(st)
            if rule == "blind":
                a0 = max(nom + int(round(W.HANDOFF_ONSET_SEARCH[0] * FS)), 0)
                b0 = min(nom + int(round(W.HANDOFF_ONSET_SEARCH[1] * FS)), len(env))
            elif rule == "narrow":
                a0, b0 = max(i_pred - int(4 * FS), 0), min(i_pred + int(4 * FS), len(env))
            elif rule == "peakcut":
                a0, b0 = max(i_pred - int(10 * FS), 0), min(i_pred + int(10 * FS), len(env))
            else:
                a0 = b0 = None
            i_on = i_pred if a0 is None else int(a0 + np.argmax(env[a0:b0]))
            a = i_on - int(round(W.TPL0_PRE_S * FS))
            if a < 0 or a + n > st.shape[1]:
                continue
            Wrow[si] = V_.norm_channels(st[:, a:a + n])
            orow[si] = -W.HANDOFF_CUT_PRE_S + a / FS
            diag.append(dict(rule=rule, cid=cid, sta=s, n_stacked=len(segs),
                             onset_minus_pred_s=float((i_on - i_pred) / FS),
                             s_snr_db=float(S.s_snr_db(st, i_on)),
                             s_snr_at_pred_db=float(S.s_snr_db(st, i_pred))))
        cids.append(cid)
        Wb.append(Wrow)
        off.append(orow)
    Sd = pd.DataFrame(rows)
    M = pd.concat(mems, ignore_index=True)
    Wb = np.stack(Wb)
    off = np.stack(off)
    ok = np.isfinite(Wb).all(axis=(2, 3))
    nsite = np.array([len({of[det[i]] for i in np.nonzero(ok[b])[0]}) for b in range(len(cids))])
    sel = np.nonzero(nsite >= W.MIN_SITES)[0]
    p = f"{out_dir}/bank_{rule}.npz"
    V_.save_bank(p, [cids[i] for i in sel], det, Wb[sel], off[sel],
                 dict(kind="oracle_" + rule, len_s=W.TPL0_LEN_S, pre_s=W.TPL0_PRE_S))
    pd.DataFrame(diag).to_csv(f"{out_dir}/diag_{rule}.csv", index=False)
    W.log(f"[{rule}] {len(sel)}/{len(cids)} templates with >= {W.MIN_SITES} sites")
    return p, Sd, M


def evaluate(bank, Sd, M, days, out_dir, rule, floor=5.0):
    import v_common as V
    import g_rediscover as GR
    import p_common as P
    V_, E = W.patch_v()
    tids, stas, Wt, off, meta = V_.load_bank(bank)
    H = E.grow(bank, days, f"{out_dir}/hits_{rule}.parquet", thr=None, floor=floor, block=32,
               tag=f"/{rule}")
    B = P.bostock()
    B = B[B.year == W.STEP1_YEAR]
    ev_all = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    lo = min(d.timestamp() for d in days)
    hi = max(d.timestamp() for d in days) + 86400
    rows = []
    for bi, tid in enumerate(tids):
        r0 = Sd[Sd.cid == tid].iloc[0]
        bk, shift = int(r0.bkey), float(r0.axis_shift_s)
        seed_t = np.sort(M[M.cid == tid].t.to_numpy(float))
        g = H[H.tpl == bi] if len(H) else H
        th = np.sort(g.t.to_numpy(float)) if len(g) else np.zeros(0)
        if len(th):
            d0 = np.array([np.min(np.abs(th - x)) for x in seed_t])
            self_frac = float((d0 <= 2.0).mean())
        else:
            self_frac = 0.0
        e = ev_all.get(bk, np.zeros(0))
        e = e[(e >= lo) & (e <= hi)]
        tor = th - shift
        if len(tor) and len(e):
            mr, _ = GR.one2one(e, tor, 2.0)
            mp, _ = GR.one2one(tor, e, 2.0)
            rec, pur = float(mr.mean()), float(mp.mean())
        else:
            rec, pur = 0.0, np.nan
        rows.append(dict(rule=rule, tid=tid, bkey=bk, n_hits=int(len(th)),
                         self_detect_frac=self_frac, n_catalog=int(len(e)),
                         recall=rec, purity=pur,
                         median_net_mad=float(np.median(g.net_mad)) if len(g) else np.nan))
        W.log(f"[{rule}] {tid}: hits {len(th)}, self {self_frac:.2f}, recall {rec:.3f}, "
              f"purity {pur if not np.isnan(pur) else float('nan'):.2f}")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/onset")
    ap.add_argument("--n-fam", type=int, default=8)
    ap.add_argument("--n-days", type=int, default=3)
    ap.add_argument("--rules", default=",".join(RULES))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    mem, of, order, det = W.sites()
    FT = Z.fam_tt_map(det, W.station_table())
    hid, iso, _ = S.holdouts()
    bkeys = sorted(int(x) for x in iso)[:a.n_fam]
    FH = ZD.fam_hours(bkeys)
    seed_days = sorted({v[0] for v in FH.values()})
    days = [d for d in W.bed_days() if d.strftime("%Y-%m-%d") in seed_days][:a.n_days]
    W.log(f"[onset] families {bkeys}; days {[d.strftime('%Y-%m-%d') for d in days]}")
    allr = []
    for rule in a.rules.split(","):
        bank, Sd, M = build(rule, bkeys, FT, det, of, a.out)
        allr.append(evaluate(bank, Sd, M, days, a.out, rule))
    R = pd.concat(allr, ignore_index=True)
    R.to_csv(f"{a.out}/onset_compare.csv", index=False)
    summ = {r: dict(n_templates=int(len(g)), median_hits=float(g.n_hits.median()),
                    median_self_detect=float(g.self_detect_frac.median()),
                    mean_self_detect=float(g.self_detect_frac.mean()),
                    median_recall=float(g.recall.median()),
                    mean_recall=float(g.recall.mean()),
                    median_purity=float(np.nanmedian(g.purity)) if g.purity.notna().any() else None)
            for r, g in R.groupby("rule")}
    json.dump(summ, open(f"{a.out}/onset_compare.json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))
    print(R.to_string())


if __name__ == "__main__":
    main()
