#!/usr/bin/env python
"""E6 step 3 -- template LEVERS (a) and (c) (Merlin E5r-audit prescription 3; results/merlin_audit_e5r_2026-09-11.md).
Reuses the E5r stacking routines (scripts/e5r_templates.py: cut_windows, align_stack, pick_s, s_snr_db) unchanged; writes NEW
e6_-prefixed template / weight / QC files; never touches results/e5r_*.  2012 BUILD half only (the test half is never read).

(a) --lever a : for UW.LRIV and PB.B001-B007 (F9: their S picks at +1.7..+3.3 s sit at the +-4 s search edge) rebuild the 2012
    templates with the S-pick search widened to +-8 s and centred on the station's MEDIAN residual (median s_pick_s of the
    station's E5r usable templates, e5r_template_qc_2012.csv); everything else (stack, alignment, null, usability rule, 8-s
    template from S-1) as E5r.  The other 8 stations keep their E5r templates.
(c) --lever c : for the 7 persistent stations (CN.PGC, CN.VGZ, CN.SNB, PB.B011, PB.B926, UW.LRIV, CN.NLLB) combine the 2012 BUILD
    stack with a 2010 stack of the SAME family (SeisBench id; Bostock 2010 events with an accepted E3 offset; cache20hz/*/2010):
    the 2010 stack is aligned to the 2012 build stack by cc over the 8-s S window (+-3 s), then N-weighted mean; S re-picked
    within +-4 s of the 2012 pick; null = the same procedure on random-time 2010 + 2012 stacks (same N each); usable = S-SNR >=
    null + 3 dB; where no 2010 events exist or the combination is not usable the 2012-only template is KEPT (reported).
Outputs results/e6_templates_2012_<lever>.npz (tpl, sp -- same layout as e5r_templates_2012.npz), e6_weights_snr_2012_<lever>.csv
        (recomputed over the merged usable set: w = 10^(S-SNR/20), cap 3x median, median-normalised), e6_template_qc_2012_<lever>.csv."""
from __future__ import annotations
import argparse, os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e5r_templates as T  # noqa: E402  (FS, HALF_W, PRE_TPL, PRE_TT, load_days, cut_windows, align_stack, pick_s, s_snr_db, epi_km)
import d3_vel1d as v1d  # noqa: E402
RES = T.RES
FS = T.FS
LENS = (8.0,)          # E6: 8-s templates only (GPU time)
STA_A = ["LRIV", "B001", "B003", "B004", "B005", "B006", "B007"]
STA_C = ["PGC", "VGZ", "SNB", "B011", "B926", "LRIV", "NLLB"]


def make_template(st, sp_idx):
    d, ok = {}, True
    for L in LENS:
        nL = int(round(L * FS)); i0 = sp_idx - int(round(T.PRE_TPL * FS))
        if i0 < 0 or i0 + nL > st.shape[1]:
            return None
        seg = st[:, i0:i0 + nL].astype(np.float32); nn = np.linalg.norm(seg, axis=1)
        if np.any(nn <= 0):
            return None
        d[L] = (seg / nn[:, None]).astype(np.float32)
    return d


def station_a(args):
    """lever (a): E5r do_station with the pick search +-8 s centred on the station's median residual."""
    sta, netsta, lat_s, lon_s, fam_tbl, ev, resid_s, seed = args
    t_st = time.time()
    D = T.load_days(netsta, 2012)
    rng = np.random.default_rng(seed)
    n = int(round(2 * T.HALF_W * FS)); idx_nom = int(round(T.HALF_W * FS)); idx_c = idx_nom + int(round(resid_s * FS))
    tmin, tmax = float(ev.t_origin.min()), float(ev.t_origin.max())
    rows, TPL, SP, STK = [], {}, {}, {}
    for fr in fam_tbl:
        g = ev[ev.bkey == fr["bkey"]]
        if len(g) < 30:
            continue
        dx = T.epi_km(fr["lat"], fr["lon"], lat_s, lon_s)
        ttP = float(v1d.tt(dx, fr["dep"], "P")); ttS = float(v1d.tt(dx, fr["dep"], "S"))
        S, keep = T.cut_windows(D, g.t_origin.to_numpy() + ttS - T.HALF_W, n)
        if S is None or len(S) < 20:
            continue
        st, Sal, sh = T.align_stack(S, iters=2)
        sp_idx = T.pick_s(st, idx_c, search_s=8.0); snr = T.s_snr_db(st, sp_idx)          # LEVER (a): +-8 s about nominal + residual
        Sn, _ = T.cut_windows(D, rng.uniform(tmin, tmax, len(S)) + ttS - T.HALF_W, n)
        snr_n = np.nan
        if Sn is not None and len(Sn) >= 10:
            stn, _, _ = T.align_stack(Sn, iters=2); snr_n = T.s_snr_db(stn, T.pick_s(stn, idx_c, search_s=8.0))
        usable = bool(np.isfinite(snr) and np.isfinite(snr_n) and snr >= snr_n + 3.0)
        d = make_template(st, sp_idx)
        if d is not None:
            key = (float(fr["bkey"]), sta)
            TPL[key] = d
            SP[key] = (ttS - ttP) + (sp_idx - idx_nom) / FS - (T.PRE_TPL - T.PRE_TT)     # same placement bookkeeping as E5r
            STK[f"{sta}|{int(fr['bkey'])}"] = st.astype(np.float32)
        rows.append(dict(sta=sta, bkey=int(fr["bkey"]), n_build=len(g), n_events=len(S), epi_km=dx, ttP=ttP, ttS=ttS,
                         resid_centre_s=resid_s, s_pick_s=(sp_idx - idx_nom) / FS, snr_db=snr, snr_null_db=snr_n,
                         usable=usable, template=d is not None, lever="a"))
    print(f"  [a {sta}] {len(rows)} stacks, {sum(r['usable'] and r['template'] for r in rows)} usable, {time.time()-t_st:.0f} s", flush=True)
    return sta, rows, TPL, SP, STK


def cc_align(ref, x, i0, nL, lag):
    """best shift k (samples) of x against ref over ref[:, i0:i0+nL]: max normalised cc of roll(x, k); returns (k, cc)."""
    A = ref[:, i0:i0 + nL]; A = A - A.mean(1, keepdims=True); na = np.sqrt((A ** 2).sum())
    best, kb = -9.0, 0
    for k in range(-lag, lag + 1):
        B = np.roll(x, k, axis=1)[:, i0:i0 + nL]; B = B - B.mean(1, keepdims=True)
        c = float((A * B).sum() / max(na * np.sqrt((B ** 2).sum()), 1e-20))
        if c > best:
            best, kb = c, k
    return kb, best


def station_c(args):
    """lever (c): 2012 build stack (E5r npz) + aligned 2010 stack of the same family -> combined template."""
    sta, netsta, lat_s, lon_s, fam_tbl, ev10, qc12, stk12, ev, seed = args     # ev = the 2012 BUILD events (random-null time span only)
    t_st = time.time()
    D10 = T.load_days(netsta, 2010)
    rng = np.random.default_rng(seed)
    n = int(round(2 * T.HALF_W * FS)); idx_nom = int(round(T.HALF_W * FS)); n8 = int(8 * FS)
    tmin10, tmax10 = float(ev10.t_origin.min()), float(ev10.t_origin.max())
    D12 = T.load_days(netsta, 2012)
    rows, TPL, SP, STK = [], {}, {}, {}
    for fr in fam_tbl:
        q = qc12[(qc12.sta == sta) & (qc12.bkey == fr["bkey"])]
        if not len(q) or not bool(q.usable.iloc[0]) or not bool(q.template.iloc[0]):
            continue                                  # no 2012 template at this station -> nothing to combine
        q = q.iloc[0]
        key12 = f"{sta}|{int(fr['bkey'])}"
        if key12 not in stk12:
            continue
        st12 = stk12[key12].astype(np.float64); n12 = int(q.n_events); sp12 = idx_nom + int(round(q.s_pick_s * FS))
        row = dict(sta=sta, bkey=int(fr["bkey"]), sb_fam_id=fr["sb_fam_id"], n12=n12, snr12_db=float(q.snr_db), s_pick12_s=float(q.s_pick_s),
                   n10=0, cc_align=np.nan, shift10_s=np.nan, snr_comb_db=np.nan, snr_comb_null_db=np.nan, usable_comb=False, used="2012-only", lever="c")
        g10 = ev10[ev10.fam_id == fr["sb_fam_id"]] if np.isfinite(fr["sb_fam_id"]) else ev10.iloc[0:0]
        if len(g10) >= 20:
            dx = T.epi_km(fr["lat"], fr["lon"], lat_s, lon_s)
            ttP = float(v1d.tt(dx, fr["dep"], "P")); ttS = float(v1d.tt(dx, fr["dep"], "S"))
            S10, _ = T.cut_windows(D10, g10.t_origin.to_numpy() + ttS - T.HALF_W, n)
            if S10 is not None and len(S10) >= 20:
                st10, _, _ = T.align_stack(S10, iters=2); n10 = len(S10)
                i0 = sp12 - int(1.0 * FS)
                if i0 >= 0 and i0 + n8 <= n:
                    k, cc = cc_align(st12, st10, i0, n8, int(3 * FS))       # align the 2010 stack to the 2012 BUILD stack (+-3 s)
                    st = (n12 * st12 + n10 * np.roll(st10, k, axis=1)) / (n12 + n10)
                    sp_idx = T.pick_s(st, sp12, search_s=4.0); snr = T.s_snr_db(st, sp_idx)
                    # null: random-time 2010 + 2012 stacks (same N each) through the SAME alignment + combination + pick
                    Sn10, _ = T.cut_windows(D10, rng.uniform(tmin10, tmax10, n10) + ttS - T.HALF_W, n)
                    t12 = ev.t_origin.to_numpy()
                    Sn12, _ = T.cut_windows(D12, rng.uniform(t12.min(), t12.max(), n12) + ttS - T.HALF_W, n)
                    snr_n = np.nan
                    if Sn10 is not None and Sn12 is not None and len(Sn10) >= 10 and len(Sn12) >= 10:
                        a10, _, _ = T.align_stack(Sn10, iters=2); a12, _, _ = T.align_stack(Sn12, iters=2)
                        kn, _ = cc_align(a12, a10, i0, n8, int(3 * FS))
                        stn = (len(Sn12) * a12 + len(Sn10) * np.roll(a10, kn, axis=1)) / (len(Sn12) + len(Sn10))
                        snr_n = T.s_snr_db(stn, T.pick_s(stn, sp12, search_s=4.0))
                    usable = bool(np.isfinite(snr) and np.isfinite(snr_n) and snr >= snr_n + 3.0)
                    row.update(n10=n10, cc_align=cc, shift10_s=k / FS, snr_comb_db=snr, snr_comb_null_db=snr_n, usable_comb=usable,
                               s_pick_comb_s=(sp_idx - idx_nom) / FS)
                    d = make_template(st, sp_idx)
                    if usable and d is not None:
                        key = (float(fr["bkey"]), sta)
                        TPL[key] = d
                        SP[key] = (ttS - ttP) + (sp_idx - idx_nom) / FS - (T.PRE_TPL - T.PRE_TT)
                        STK[f"{sta}|{int(fr['bkey'])}"] = st.astype(np.float32)
                        row["used"] = "combined"
        rows.append(row)
    print(f"  [c {sta}] {len(rows)} 2012 templates, {sum(r['used']=='combined' for r in rows)} combined, {time.time()-t_st:.0f} s", flush=True)
    return sta, rows, TPL, SP, STK


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lever", required=True, choices=["a", "c"])
    ap.add_argument("--workers", type=int, default=7)
    a = ap.parse_args()
    S = pd.read_csv(f"{RES}/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv")
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet")
    ev = gt[(gt.half == "build") & (~gt.guard)][["bkey", "t_origin"]].reset_index(drop=True)      # BUILD half only
    fam_tbl = F[F.n_build >= 30][["bkey", "sb_fam_id", "lat", "lon", "dep", "n_build"]].to_dict("records")
    qc12 = pd.read_csv(f"{RES}/e5r_template_qc_2012.csv")
    z = np.load(f"{RES}/e5r_templates_2012.npz", allow_pickle=True)
    tpl = z["tpl"].item(); sp = z["sp"].item()
    tpl = {k: {L: v[L] for L in LENS} for k, v in tpl.items()}          # 8 s only
    if a.lever == "a":
        resid = qc12[qc12.usable & qc12.template & qc12.sta.isin(STA_A)].groupby("sta").s_pick_s.median()
        print("[a] per-station median residual (s), the new search centre: " + ", ".join(f"{s} {resid.get(s, 0):+.2f}" for s in STA_A), flush=True)
        jobs = [(s, S.loc[s, "netsta"], float(S.loc[s, "lat"]), float(S.loc[s, "lon"]), fam_tbl, ev, float(resid.get(s, 0.0)), 300 + i) for i, s in enumerate(STA_A)]
        fn = station_a; rebuilt = set(STA_A)
    else:
        gt10 = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
        off10 = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv"); off10 = off10[off10.accept.astype(bool)]
        omap = dict(zip(off10.fam_id.astype(float), off10.offset_s.astype(float)))
        gt10 = gt10[gt10.fam_id.isin(omap)].copy(); gt10["time"] = gt10.time.astype("datetime64[ns]")
        gt10["t_origin"] = gt10.time.astype("int64") / 1e9 + gt10.fam_id.map(omap).to_numpy()
        ev10 = gt10[["fam_id", "t_origin"]].reset_index(drop=True)
        print(f"[c] 2010 events with an accepted E3 offset: {len(ev10)} in {ev10.fam_id.nunique()} SeisBench families "
              f"({int(F.sb_fam_id.isin(ev10.fam_id.unique()).sum())} of the 130 2012 families have 2010 events)", flush=True)
        stk12 = {k: z[k] for k in z.files if "|" in k and k.split("|")[0] in STA_C}
        jobs = [(s, S.loc[s, "netsta"], float(S.loc[s, "lat"]), float(S.loc[s, "lon"]), fam_tbl, ev10, qc12, stk12, ev, 400 + i) for i, s in enumerate(STA_C)]
        fn = station_c; rebuilt = set()                                  # (c) replaces pair by pair, keeps 2012-only otherwise
    with Pool(min(a.workers, len(jobs))) as pool:
        out = pool.map(fn, jobs)
    rows, TPLn, SPn, STK = [], {}, {}, {}
    for sta, rw, Tn, Pn, Kn in out:
        rows += rw; TPLn.update(Tn); SPn.update(Pn); STK.update(Kn)
    Q = pd.DataFrame(rows); Q.to_csv(f"{RES}/e6_template_qc_2012_{a.lever}.csv", index=False)
    if a.lever == "a":
        use_new = Q[Q.usable & Q.template]
        newkeys = {(float(r.bkey), r.sta) for r in use_new.itertuples()}
        merged = {k: v for k, v in tpl.items() if k[1] not in rebuilt}
        merged.update({k: TPLn[k] for k in newkeys})
        sp_m = {k: (SPn[k] if k in newkeys else sp[k]) for k in merged}
        old_n = sum(1 for k in tpl if k[1] in rebuilt)
        print(f"[a] rebuilt stations: {old_n} E5r usable templates -> {len(newkeys)} usable with the +-8 s / residual-centred pick; "
              f"merged set {len(merged)} (E5r 1,144)", flush=True)
        snr_of = {(float(r.bkey), r.sta): float(r.snr_db) for r in use_new.itertuples()}
    else:
        merged = dict(tpl); merged.update(TPLn); sp_m = dict(sp); sp_m.update(SPn)
        comb = Q[Q.used == "combined"]
        print(f"[c] {len(comb)} of {len(Q)} persistent-station templates replaced by 2010+2012 combined stacks "
              f"(median n10 {comb.n10.median() if len(comb) else float('nan'):.0f}, cc_align median {comb.cc_align.median() if len(comb) else float('nan'):.2f}, "
              f"S-SNR 2012-only -> combined median {comb.snr12_db.median() if len(comb) else float('nan'):.1f} -> {comb.snr_comb_db.median() if len(comb) else float('nan'):.1f} dB)", flush=True)
        snr_of = {(float(r.bkey), r.sta): float(r.snr_comb_db) for r in comb.itertuples()}
    # weights over the merged usable set (E5r rule): w = 10^(S-SNR/20), cap 3x median, median-normalised
    W0 = pd.read_csv(f"{RES}/e5r_weights_snr_2012.csv")
    snr_all = {(float(r.fam), r.sta): float(r.snr_db) for r in W0.itertuples()}
    snr_all = {k: v for k, v in snr_all.items() if k in merged and not (a.lever == "a" and k[1] in rebuilt)}
    snr_all.update({k: v for k, v in snr_of.items() if k in merged})
    keys = sorted(merged, key=lambda k: (k[1], k[0]))
    snr_v = np.array([snr_all[k] for k in keys]); w = 10 ** (snr_v / 20.0); w = np.minimum(w, 3 * np.median(w)); w = w / np.median(w)
    Wn = pd.DataFrame(dict(fam=[k[0] for k in keys], sta=[k[1] for k in keys], w=w, snr_db=snr_v))
    Wn.to_csv(f"{RES}/e6_weights_snr_2012_{a.lever}.csv", index=False)
    np.savez_compressed(f"{RES}/e6_templates_2012_{a.lever}.npz", tpl=np.array(merged, dtype=object), sp=np.array(sp_m, dtype=object), **STK)
    print(f"[weights] per-station median w: {Wn.groupby('sta').w.median().round(2).to_dict()}")
    print(f"[write] results/e6_templates_2012_{a.lever}.npz ({len(merged)} templates), e6_weights_snr_2012_{a.lever}.csv, e6_template_qc_2012_{a.lever}.csv")


if __name__ == "__main__":
    main()
