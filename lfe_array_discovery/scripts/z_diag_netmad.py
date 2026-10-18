#!/usr/bin/env python
"""DEFECT OR PHYSICS?  Two measurements the hit counts cannot distinguish.

  N  NET_MAD AT THE CATALOGUED EVENT TIMES.  The oracle template's detection statistic sampled at
     every catalogued event of its own family on the seed day, against the same statistic at random
     times.  A distribution centred just below the threshold means the bar is the problem; one that
     sits on top of the random distribution means there is no detectable signal at all.
  S  STACK COHERENCE vs N.  The per-station stack's S-SNR and envelope-peak height as a function of
     the number of stacked events (5 ... 160).  Coherent stacking grows as sqrt(N); a flat curve
     means the stack never becomes coherent -- catalogue-time jitter or the wrong time base.
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
import z_diag_onset as ZO                                                   # noqa: E402
import s_common as S                                                        # noqa: E402

FS = W.FS


def nm_series(bank, date, bi):
    """(tg0, net_mad (NG,), per-site z) for ONE template over one day -- run_day's own math."""
    import torch
    import v_common as V
    import v_engine as E
    V_, E_ = W.patch_v()
    tids, stas, Wt, off, meta = V_.load_bank(bank)
    n = Wt.shape[3]
    dd = E_.DayData(date, stas, n, dev="cuda", idx=V_.cache_idx())
    if len(dd.sites) < W.MIN_SITES:
        return None
    have = {s: i for i, s in enumerate(stas)}
    A = [s for s in dd.stas if s in have]
    t0max = max(dd.t0[s] for s in A)
    oo = off[[bi]][:, [have[s] for s in A]]
    offmin = float(np.nanmin(oo))
    tg0 = t0max + max(0.0, -offmin) + 1.0
    K, NG = {}, 10 ** 9
    for s in A:
        o = off[bi, have[s]]
        k = int(np.rint((tg0 + (o if np.isfinite(o) else 0.0) - dd.t0[s]) * FS))
        k = int(np.clip(k, 0, max(dd.NV[s] - 2, 0)))
        K[s] = k
        NG = min(NG, int(dd.NV[s] - k - 1))
    NG = int(min(NG, 86400 * FS))
    acc = torch.zeros((1, NG), device="cuda")
    nz3 = torch.zeros((1, NG), device="cuda")
    nsv = torch.zeros((1, NG), device="cuda")
    SZ = {}
    for q in dd.sites:
        mem = [s for s in V_.SITE_MEMBERS[q] if s in A]
        sbuf = torch.zeros((1, NG), device="cuda")
        cnt = torch.zeros((1, NG), device="cuda")
        for s in mem:
            w = Wt[bi, have[s]]
            if not np.isfinite(w).all():
                continue
            Zz = dd.station_z(s, torch.from_numpy(np.ascontiguousarray(w[None])).to("cuda"))
            vm = dd.valid[s].float()[K[s]:K[s] + NG]
            sbuf += Zz[:, K[s]:K[s] + NG] * vm[None]
            cnt += vm[None]
            del Zz
        okc = cnt > 0
        sz = torch.where(okc, sbuf / cnt.clamp_min(1.0), torch.zeros((), device="cuda"))
        sz = sz.clamp(max=W.AC_SITE_Z_CLIP)
        acc += torch.where(okc, sz, torch.zeros((), device="cuda"))
        nz3 += ((sz >= W.AC_SITE_Z_COINC) & okc).float()
        nsv += okc.float()
        SZ[q] = sz[0].detach().cpu().numpy()
        del sbuf, cnt, okc, sz
    m = acc.median(dim=-1, keepdim=True).values
    qq = torch.kthvalue(acc, max(1, NG // 4), dim=-1, keepdim=True).values
    sc = (1.4826 * (m - qq)).clamp_min(1e-9)
    nm = ((acc - m) / sc)[0].detach().cpu().numpy()
    z3 = nz3[0].detach().cpu().numpy()
    ns = nsv[0].detach().cpu().numpy()
    dd.free()
    del acc, nz3, nsv
    return tg0, nm, z3, ns, SZ


def part_N(out, bkeys, FT, det, of, work):
    import p_common as P
    bank, Sd, M = ZO.build("ttcut", bkeys, FT, det, of, work)
    import v_common as V
    V_, _ = W.patch_v()
    tids, *_ = V_.load_bank(bank)
    B = P.bostock()
    B = B[B.year == W.STEP1_YEAR]
    rows = []
    rng = np.random.default_rng(20260920)
    for bi, tid in enumerate(tids):
        r0 = Sd[Sd.cid == tid].iloc[0]
        bk, shift, d = int(r0.bkey), float(r0.axis_shift_s), str(r0.date)
        r = nm_series(bank, pd.Timestamp(d), bi)
        if r is None:
            continue
        tg0, nm, z3, ns, SZ = r
        ev = np.sort(B[(B.bkey == bk)].t_origin.to_numpy(float))
        ts = pd.Timestamp(d).timestamp()
        ev = ev[(ev >= ts) & (ev < ts + 86400)] + shift
        k = np.rint((ev - tg0) * FS).astype(int)
        k = k[(k >= 0) & (k < len(nm))]
        if len(k) < 5:
            continue
        # +- 1 s tolerance: take the max of the statistic in the window
        wj = int(FS)
        v = np.array([nm[max(x - wj, 0):x + wj + 1].max() for x in k])
        kr = rng.integers(wj + 1, len(nm) - wj - 1, 4000)
        vr = np.array([nm[x - wj:x + wj + 1].max() for x in kr])
        sites = {q: float(np.median([SZ[q][max(x - wj, 0):x + wj + 1].max() for x in k]))
                 for q in SZ}
        rows.append(dict(bkey=bk, date=d, n_events=len(k),
                         nm_med=float(np.median(v)), nm_p90=float(np.percentile(v, 90)),
                         nm_p99=float(np.percentile(v, 99)), nm_max=float(v.max()),
                         frac_ge5=float((v >= 5).mean()), frac_ge8=float((v >= 8).mean()),
                         null_med=float(np.median(vr)), null_p99=float(np.percentile(vr, 99)),
                         null_max=float(vr.max()),
                         frac_gt_null_p99=float((v > np.percentile(vr, 99)).mean()),
                         site_z=json.dumps({k2: round(v2, 2) for k2, v2 in sites.items()})))
        W.log(f"[N] {tid} bk{bk}: n {len(k)}, nm med {rows[-1]['nm_med']:.2f} "
              f"(null {rows[-1]['null_med']:.2f}), frac>=5 {rows[-1]['frac_ge5']:.3f}, "
              f"frac>null_p99 {rows[-1]['frac_gt_null_p99']:.3f}")
    R = pd.DataFrame(rows)
    R.to_csv(out + "_N_netmad_at_events.csv", index=False)
    return R


def part_S(out, bkeys, FT, det, of, n_list=(5, 10, 20, 40, 80, 160)):
    """Stack S-SNR and envelope-peak MAD vs the number of stacked events, at the near sites."""
    import p_common as P
    cut = Z.Cutter(maxn=40)
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = P.bostock()
    B = B[(B.year == W.STEP1_YEAR) & (B.bkey.isin(bkeys))].copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B = B[B.date.isin(days)]
    near = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB"]
    rows = []
    PRE, LEN = 20.0, 60.0
    nom = int(round(PRE * FS))
    for bk in bkeys:
        ev = np.sort(B[B.bkey == bk].t_origin.to_numpy(float))
        for s in near:
            tt = FT.get((int(bk), s))
            if tt is None:
                continue
            segs = []
            for t in ev:
                if len(segs) >= max(n_list):
                    break
                g = cut.cut(s, t + tt - PRE, LEN)
                if g is not None:
                    segs.append(g)
            if len(segs) < min(n_list):
                continue
            Sg = Z.norm_rows(np.stack(segs))
            for N in n_list:
                if len(Sg) < N:
                    continue
                st = Sg[:N].mean(0)
                env = Z.envelope_1d(st)
                a, b = max(nom - int(6 * FS), 0), min(nom + int(6 * FS), len(env))
                i = int(a + np.argmax(env[a:b]))
                rows.append(dict(bkey=int(bk), sta=s, N=N,
                                 snr_at_pred_db=float(S.s_snr_db(st, nom)),
                                 snr_at_peak_db=float(S.s_snr_db(st, i)),
                                 peak_mad=Z.mad_height(env, i),
                                 peak_offset_s=float((i - nom) / FS)))
        W.log(f"[S] family {bk}: {len(ev)} events in the window")
    R = pd.DataFrame(rows)
    R.to_csv(out + "_S_stack_vs_N.csv", index=False)
    G = R.groupby("N").agg(n=("peak_mad", "size"), med_peak_mad=("peak_mad", "median"),
                           med_snr_peak=("snr_at_peak_db", "median"),
                           med_snr_pred=("snr_at_pred_db", "median")).reset_index()
    print(G.to_string())
    return R, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/nmdiag/diag")
    ap.add_argument("--work", default="/tmp/nmdiag/work")
    ap.add_argument("--n-fam", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    os.makedirs(a.work, exist_ok=True)
    mem, of, order, det = W.sites()
    FT = Z.fam_tt_map(det, W.station_table())
    hid, iso, _ = S.holdouts()
    bkeys = sorted(int(x) for x in iso)[:a.n_fam]
    RN = part_N(a.out, bkeys, FT, det, of, a.work)
    RS, G = part_S(a.out, bkeys, FT, det, of)
    summ = dict(
        families=bkeys,
        netmad_at_events=dict(
            n=int(len(RN)),
            median_nm_at_events=float(RN.nm_med.median()) if len(RN) else None,
            median_null=float(RN.null_med.median()) if len(RN) else None,
            median_frac_ge5=float(RN.frac_ge5.median()) if len(RN) else None,
            median_frac_gt_null_p99=float(RN.frac_gt_null_p99.median()) if len(RN) else None,
            per_family=RN.to_dict("records") if len(RN) else []),
        stack_vs_N=G.to_dict("records"))
    json.dump(summ, open(a.out + "_netmad.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in summ.items() if k != "netmad_at_events"}, indent=1,
                     default=str))
    print(json.dumps({k: v for k, v in summ["netmad_at_events"].items() if k != "per_family"},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
