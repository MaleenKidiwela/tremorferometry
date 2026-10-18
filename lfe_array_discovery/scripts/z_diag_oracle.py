#!/usr/bin/env python
"""WHY DID THE ORACLE RECOVER NOTHING?  Diagnostic a/b/c (coordinator, 2026-09-20).

C_30 = 0 contradicts two measurements (Merlin's 4-15 MAD envelope stacks at predicted S on this
bed; P4-0b truth-stack recall ~0.45 on a worse bed), so it is treated as a suspected DEFECT, not a
bed verdict, until these three tests say otherwise.

  a  TIME BASE.  Per family: envelope-stack peak offset at each site against the predicted S, and
     the stack's S-SNR under three alignment settings -- NO align, align +- 1.5 s (frozen, lag 30)
     and align +- 6 s (lag 120).  Reports the fraction of families whose peak offset exceeds 1.5 s
     at half or more of their sites.
  b  TEMPLATE SANITY.  For 5 families, build the oracle template under each alignment setting and
     ask whether it detects the very events it was built from (net_mad at the seed times).
  c  RECALL.  The same templates over 3 ETS days at T = 5: one-to-one recall and purity (+- 2 s)
     against the family's catalogued events.

Nothing here writes into the frozen 1c' outputs; everything goes to --out.
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
import s_common as S                                                        # noqa: E402

FS = W.FS
VARIANTS = {"noalign": (0, 30), "align1.5s": (2, 30), "align6s": (2, 120)}
CUT_PRE, CUT_LEN = 20.0, 60.0          # per-station cut around the PREDICTED S
PEAK_SEARCH = 10.0


def fam_hours(bkeys, max_events=40):
    """{bkey: (date, hour, event origin times)} -- the family's busiest hour, as the oracle does."""
    import p_common as P
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = P.bostock()
    B = B[(B.year == W.STEP1_YEAR) & (B.bkey.isin(bkeys))].copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B["hh"] = B.time.dt.hour.astype(int)
    B = B[B.date.isin(days)]
    out = {}
    for bk in bkeys:
        g = B[B.bkey == bk]
        if not len(g):
            continue
        k = g.groupby(["date", "hh"]).size().sort_values(ascending=False)
        d, h = k.index[0]
        ev = np.sort(g[(g.date == d) & (g.hh == h)].t_origin.to_numpy(float))
        if len(ev) > max_events:
            ev = ev[np.linspace(0, len(ev) - 1, max_events).round().astype(int)]
        out[int(bk)] = (str(d), int(h), ev)
    return out


def part_a(out, bkeys, FT, det, of):
    cut = Z.Cutter(maxn=40)
    FH = fam_hours(bkeys)
    rows = []
    for bk, (d, h, ev) in FH.items():
        for s in det:
            tt = FT.get((bk, s))
            if tt is None:
                continue
            segs = [g for g in (cut.cut(s, t + tt - CUT_PRE, CUT_LEN) for t in ev)
                    if g is not None]
            if len(segs) < 5:
                continue
            Sg = Z.norm_rows(np.stack(segs))
            if len(Sg) < 5:
                continue
            nom = int(round(CUT_PRE * FS))
            r = dict(bkey=bk, date=d, hour=h, sta=s, site=of.get(s, s), n=len(Sg))
            for tag, (it, lag) in VARIANTS.items():
                st = (Sg.mean(0) if it == 0
                      else S.align_stack(Sg.copy(), iters=it, lag=lag)[0])
                env = Z.envelope_1d(st)
                a = max(nom - int(PEAK_SEARCH * FS), 0)
                b = min(nom + int(PEAK_SEARCH * FS), len(env))
                i = int(a + np.argmax(env[a:b]))
                r[f"off_{tag}"] = float((i - nom) / FS)
                r[f"mad_{tag}"] = Z.mad_height(env, i)
                r[f"snr_{tag}"] = float(S.s_snr_db(st, nom))       # at the PREDICTED S
                r[f"snrpk_{tag}"] = float(S.s_snr_db(st, i))       # at the measured peak
            rows.append(r)
        W.log(f"[a] family {bk} {d} {h:02d}h: {len(ev)} events")
    A = pd.DataFrame(rows)
    A.to_csv(out + "_a_offsets.csv", index=False)
    per = []
    for bk, g in A.groupby("bkey"):
        o = g.off_noalign.abs()
        per.append(dict(bkey=int(bk), n_sites=int(len(g)),
                        median_abs_offset_s=float(o.median()),
                        frac_sites_gt_1p5s=float((o > 1.5).mean()),
                        **{f"med_snr_{t}": float(g[f"snr_{t}"].median()) for t in VARIANTS},
                        **{f"med_mad_{t}": float(g[f"mad_{t}"].median()) for t in VARIANTS}))
    P_ = pd.DataFrame(per)
    P_.to_csv(out + "_a_per_family.csv", index=False)
    summ = dict(n_families=int(len(P_)),
                frac_families_offset_gt_1p5s_at_half_sites=float(
                    (P_.frac_sites_gt_1p5s >= 0.5).mean()),
                median_of_median_abs_offset_s=float(P_.median_abs_offset_s.median()),
                p90_median_abs_offset_s=float(P_.median_abs_offset_s.quantile(0.9)),
                median_mad_noalign=float(P_.med_mad_noalign.median()),
                median_stack_snr_db={t: float(P_[f"med_snr_{t}"].median()) for t in VARIANTS},
                median_peak_mad={t: float(P_[f"med_mad_{t}"].median()) for t in VARIANTS})
    print("[a]", json.dumps(summ, indent=1))
    return A, P_, summ


def build_banks(out, bkeys, FT, det, tag_dir):
    """The oracle handoff for `bkeys`, once per alignment variant."""
    import w_grow as G
    import v_common as V
    FH = fam_hours(bkeys)
    rows, mems = [], []
    for bk, (d, h, ev) in FH.items():
        tt = [FT[(bk, s)] for s in det if (bk, s) in FT]
        shift = float(np.median(tt)) - 1.0
        t = ev + shift
        cid = f"H{bk:04d}"
        rows.append(dict(cid=cid, date=d, n_members=len(t), n_in_hour=len(t), n_in_3h=len(t),
                         n_bins=1, t0=float(t.min()), t1=float(t.max()),
                         t_seed=float(np.median(t)), bkey=bk, hour=h, axis_shift_s=shift))
        mems.append(pd.DataFrame(dict(cid=cid, t=t)))
    Sd = pd.DataFrame(rows)
    M = pd.concat(mems, ignore_index=True)
    os.makedirs(tag_dir, exist_ok=True)
    Sd.to_csv(f"{tag_dir}/seeds.csv", index=False)
    M.to_parquet(f"{tag_dir}/clusters.parquet", index=False)
    banks = {}
    for tag, (it, lag) in VARIANTS.items():
        p = f"{tag_dir}/bank_{tag}.npz"
        W.ALIGN_ITERS, W.ALIGN_LAG = it, lag
        G.W.ALIGN_ITERS, G.W.ALIGN_LAG = it, lag
        if os.path.exists(p):
            os.remove(p)
        G.handoff(M, Sd, p)
        banks[tag] = p
    W.ALIGN_ITERS, W.ALIGN_LAG = 2, 30
    return Sd, M, banks


def part_bc(out, bkeys, FT, det, tag_dir, n_days=3):
    import v_common as V
    import g_rediscover as GR
    import p_common as P
    Sd, M, banks = build_banks(out, bkeys, FT, det, tag_dir)
    V_, E = W.patch_v()
    B = P.bostock()
    B = B[B.year == W.STEP1_YEAR]
    ev_all = {int(b): np.sort(g.t_origin.to_numpy(float)) for b, g in B.groupby("bkey")}
    seed_days = sorted({r.date for r in Sd.itertuples()})
    ets = [d for d in W.bed_days() if d.strftime("%Y-%m-%d") in seed_days][:n_days]
    days = ets if ets else W.bed_days()[-n_days:]
    W.log(f"[bc] growing on {[d.strftime('%Y-%m-%d') for d in days]}")
    rows = []
    for tag, bank in banks.items():
        tids, stas, Wt, off, meta = V_.load_bank(bank)
        H = E.grow(bank, days, f"{tag_dir}/hits_{tag}.parquet", thr=None, floor=5.0,
                   block=32, tag=f"/{tag}")
        for bi, tid in enumerate(tids):
            r0 = Sd[Sd.cid == tid].iloc[0]
            bk, shift = int(r0.bkey), float(r0.axis_shift_s)
            seed_t = np.sort(M[M.cid == tid].t.to_numpy(float))
            g = H[H.tpl == bi] if len(H) else H
            t_hit = np.sort(g.t.to_numpy(float)) if len(g) else np.zeros(0)
            # self-detection: is there a hit within +- 2 s of each seed time?
            if len(t_hit):
                k = np.searchsorted(t_hit, seed_t)
                d0 = np.full(len(seed_t), np.inf)
                for j, kk in enumerate(k):
                    for c in (kk - 1, kk):
                        if 0 <= c < len(t_hit):
                            d0[j] = min(d0[j], abs(t_hit[c] - seed_t[j]))
                self_frac = float((d0 <= 2.0).mean())
                nm = g.net_mad.to_numpy()
                med_nm = float(np.median(nm))
            else:
                self_frac, med_nm = 0.0, np.nan
            # recall / purity against the catalogued events on those days, on the ORIGIN axis
            t_org = t_hit - shift
            lo = min(d.timestamp() for d in days)
            hi = max(d.timestamp() for d in days) + 86400
            e = ev_all.get(bk, np.zeros(0))
            e = e[(e >= lo) & (e <= hi)]
            if len(t_org) and len(e):
                mr, _ = GR.one2one(e, t_org, 2.0)
                mp, _ = GR.one2one(t_org, e, 2.0)
                rec, pur = float(mr.mean()), float(mp.mean())
            else:
                rec, pur = 0.0, np.nan
            rows.append(dict(variant=tag, tid=tid, bkey=bk, n_seed=len(seed_t),
                             n_hits=int(len(t_hit)), median_net_mad=med_nm,
                             self_detect_frac=self_frac, n_catalog_on_days=int(len(e)),
                             recall=rec, purity=pur))
            W.log(f"[bc] {tag} {tid}: hits {len(t_hit)}, self {self_frac:.2f}, recall {rec:.3f}")
    R = pd.DataFrame(rows)
    R.to_csv(out + "_bc_templates.csv", index=False)
    summ = {t: dict(median_hits=float(g.n_hits.median()),
                    median_self_detect=float(g.self_detect_frac.median()),
                    median_recall=float(g.recall.median()),
                    median_purity=float(np.nanmedian(g.purity)) if g.purity.notna().any()
                    else None)
            for t, g in R.groupby("variant")}
    print("[bc]", json.dumps(summ, indent=1))
    return R, summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{W.RES}/step1c_prime/diag_oracle")
    ap.add_argument("--n-bc", type=int, default=5)
    ap.add_argument("--n-days", type=int, default=3)
    ap.add_argument("--skip-a", action="store_true")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    mem, of, order, det = W.sites()
    ST = W.station_table()
    FT = Z.fam_tt_map(det, ST)
    hid, iso, _ = S.holdouts()
    import p_common as P
    F = P.families()
    kept = [int(b) for b in F.index if int(b) not in hid][:10]
    bkeys = sorted(int(x) for x in hid) + kept
    summ = {}
    if not a.skip_a:
        A, PF, sa = part_a(a.out, bkeys, FT, det, of)
        summ["a"] = sa
    bc = sorted(int(x) for x in iso)[:a.n_bc]
    R, sbc = part_bc(a.out, bc, FT, det, f"{a.out}_bc", n_days=a.n_days)
    summ["bc"] = sbc
    summ["bc_families"] = bc
    json.dump(summ, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
