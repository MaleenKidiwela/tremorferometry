#!/usr/bin/env python
"""FINAL SPEC, STEP 2 (GPU) -- THE CORRECTED ORACLE CEILING and the 2005 reference v3.

Merlin's ruling (PLAN sec.12 step 2).  For each of the 30 HIDDEN families:

  * seed with up to MAX_STACK (400) catalogued events at the STEP-1-CORRECTED times, spread over
    the whole 24-day ETS window (the engine's own cap, w_grow.py:88-89, does the spreading);
  * t_seed / t0 / t1 are the median / min / max of the spread seed, so the near-seed burst guard
    (<= 20 % of hits within +- 10 min of EITHER end) cannot bind -- it was a burst guard applied
    to a seed chosen as a burst in 1c';
  * then the FROZEN chain unchanged: per-template quiet-day T_seed (FA <= 0.1 per template-day,
    engineering floor 5.0 reported), growth, the frozen survival gate with EVERY clause but
    counting ETS-WINDOW days ONLY (F5), 2 iterations, dedupe, locate/Gates 0-2, and the scorer.

Recall and purity are written against THREE denominators (Merlin Q4):
  (i)   all catalogued events in the window                              -- REPORTED ONLY
  (ii)  the 60-s-independent events the near-seed clause does not remove -- WHAT THE GATE COUNTS
  (iii) the LOUD subset: events whose network-median 3-s S-energy excess on the frozen cache at
        the family's four nearest admitted sites is >= 2.0 (the D8 statistic; the per-event
        numbers are written to s2_event_excess.parquet)           -- the FAIR detector ceiling

ACCEPT   C_30 >= 3.  If C_iso < 4 the pre-registered fallback bar governs: >= max(5, 0.5 C_30)
         of the 30 hidden families.
ABANDON  C_30 < 3 with Step-1 coherence verified -> the frozen chain cannot pass with the best
         possible templates -> STOP the study, and the write-up says "THIS CHAIN", not "the
         archive".
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import w_common as W                                                        # noqa: E402
import q_common as Q                                                        # noqa: E402
import q_chain as QC                                                        # noqa: E402

ABANDON_RC = 20


def build_seeds(out_dir, times_parquet):
    """Hidden-family seeds at the Step-1-corrected times, spread over the whole ETS window."""
    hidden, iso = Q.hidden_sets()
    T = pd.read_parquet(times_parquet)
    ets = set(Q.ets_dates())
    rows, mems = [], []
    for bk in sorted(hidden):
        g = T[T.bkey == int(bk)]
        if not len(g):
            continue
        t = np.sort(g.t_corrected.to_numpy(float))
        # keep only the ETS window (the quiet days are the FA calibration, never a seed)
        d = pd.to_datetime(t, unit="s").strftime("%Y-%m-%d")
        t = t[np.isin(d, list(ets))]
        if len(t) < 3:
            continue
        shift = Q.fam_shift(bk)
        tn = t + shift                       # the pair-seed NODE axis
        cid = f"H{int(bk):04d}"
        rows.append(dict(cid=cid, date=pd.Timestamp(tn[0], unit="s").strftime("%Y-%m-%d"),
                         n_members=int(len(tn)), n_in_hour=int(len(tn)), n_in_3h=int(len(tn)),
                         n_bins=int(len(tn)), t0=float(tn.min()), t1=float(tn.max()),
                         t_seed=float(np.median(tn)), max_net_mad=np.nan, n_pairs=0,
                         bkey=int(bk), hour=int(pd.Timestamp(tn[0], unit="s").hour),
                         n_events_window=int(len(tn)), axis_shift_s=float(shift),
                         n_corrected=int(g.corrected.sum()),
                         span_days=float((tn.max() - tn.min()) / 86400.0)))
        mems.append(pd.DataFrame(dict(cid=cid, t=tn)))
    Sd = pd.DataFrame(rows)
    M = pd.concat(mems, ignore_index=True) if mems else pd.DataFrame(columns=["cid", "t"])
    Sd.to_csv(f"{out_dir}/cl_seeds.csv", index=False)
    M.to_parquet(f"{out_dir}/cl_clusters.parquet", index=False)
    Q.log(f"[s2] {len(Sd)} hidden-family seeds, {int(Sd.n_members.sum()) if len(Sd) else 0} "
          f"corrected-time events, median span "
          f"{float(Sd.span_days.median()) if len(Sd) else 0:.1f} d "
          f"(the engine caps the stack at MAX_STACK = {W.MAX_STACK})")
    return Sd


def loud_events(out_parquet, Sd):
    """Denominator (iii): the per-event S-energy excess statistic, written once."""
    if os.path.exists(out_parquet):
        return pd.read_parquet(out_parquet)
    ev = Q.events_2005([int(b) for b in Sd.bkey])
    mem, of, order, det = W.sites()
    cut = Q.TwoYearCutter([2005], det, maxn=len(det) + 8)
    tabs = []
    for bk in sorted(ev):
        D = Q.event_s_excess(cut, bk, ev[bk])
        tabs.append(D)
        n_loud = int((D.excess_median >= Q.LOUD_EXCESS).sum()) if len(D) else 0
        Q.log(f"[s2/loud] family {bk}: {n_loud}/{len(D)} events with network-median 3-s S excess "
              f">= {Q.LOUD_EXCESS}")
    X = pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame(
        columns=["bkey", "t_origin", "n_sites_used", "excess_median"])
    X.to_parquet(out_parquet, index=False)
    return X


def score(d, Sd, X, out_csv, final_tag):
    import g_rediscover as GR
    import w_score_step1 as SC
    import s_common as S

    fam = Q.family_ll()
    hidden, iso = Q.hidden_sets()
    ev = Q.events_2005(hidden)
    evk = Q.events_2005([f for f in fam if f not in set(hidden)])
    CT = SC.cluster_truth(fam, {**ev, **evk}, hidden)
    kept = [f for f in fam if f not in set(hidden)]

    C = (pd.read_parquet(f"{d}/lg_clusters.parquet")
         if os.path.exists(f"{d}/lg_clusters.parquet")
         else pd.DataFrame(columns=["cid", "t", "lat", "lon"]))
    L = W.read_csv_safe(f"{d}/lg_locations.csv")
    SV = W.read_csv_safe(f"{d}/{final_tag}_survival.csv", W.SURV_COLS)
    surv = set(SV[SV.survives.astype(bool)].tid.astype(str)) if len(SV) else set()
    loc = {str(r.tid) for r in L.itertuples() if bool(getattr(r, "located", False))} if len(L) \
        else set()
    rows = []
    for r in Sd.itertuples():
        bk, cid = int(r.bkey), str(r.cid)
        g = C[C.cid == cid] if len(C) else C
        t = np.sort(g.t.to_numpy(float)) if len(g) else np.zeros(0)
        e_all = ev.get(bk, np.zeros(0))
        e_ind = Q.independent_outside(e_all, float(r.t0), float(r.t1))
        xb = X[X.bkey == bk] if len(X) else X
        e_loud = (np.sort(xb[xb.excess_median >= Q.LOUD_EXCESS].t_origin.to_numpy(float))
                  if len(xb) else np.zeros(0))
        rec = recall_ind = recall_loud = pur = pc3 = kb = np.nan
        dist = np.nan
        if len(t):
            rec = float(GR.one2one(e_all, t, W.SCORE_TOL_S)[0].mean()) if len(e_all) else 0.0
            recall_ind = (float(GR.one2one(e_ind, t, W.SCORE_TOL_S)[0].mean())
                          if len(e_ind) else np.nan)
            recall_loud = (float(GR.one2one(e_loud, t, W.SCORE_TOL_S)[0].mean())
                           if len(e_loud) else np.nan)
            pur = float(GR.one2one(t, e_all, W.SCORE_TOL_S)[0].mean())
            pc3 = float(GR.one2one(t, CT.get(bk, np.zeros(0)), W.SCORE_TOL_S)[0].mean())
            kb = max((GR.one2one(t, evk.get(k, np.zeros(0)), W.SCORE_TOL_S)[0].mean()
                      for k in kept), default=0.0)
            dist = float(S.epi_km(float(g.lat.iloc[0]), float(g.lon.iloc[0]), *fam[bk]))
        counts = bool(len(t) and np.isfinite(dist) and dist <= W.SCORE_R_KM and rec > 0
                      and (not np.isfinite(kb) or kb <= W.SCORE_MAX_EXPLAINED))
        sv = SV[SV.tid.astype(str) == cid] if len(SV) else SV
        rows.append(dict(
            bkey=bk, stratum=("isolated" if bk in set(iso) else "other"),
            n_seed_events=int(r.n_members), seed_span_days=float(r.span_days),
            n_events_all=int(len(e_all)), n_events_indep=int(len(e_ind)),
            n_events_loud=int(len(e_loud)),
            survived=bool(cid in surv), located=bool(cid in loc), n_cluster=int(len(t)),
            dist_km=dist,
            recall_ref=rec, recall_indep_ref=recall_ind, recall_loud_ref=recall_loud,
            purity_ref=pur, purity_cluster3km_ref=pc3, kept_explained=kb,
            n_indep=int(sv.n_indep.iloc[0]) if len(sv) else 0,
            n_days=int(sv.n_days.iloc[0]) if len(sv) else 0,
            n_indep_quiet=int(sv.n_indep_quiet.iloc[0])
            if (len(sv) and "n_indep_quiet" in sv.columns
                and pd.notna(sv.n_indep_quiet.iloc[0])) else -1,
            frac_near_seed=float(sv.frac_near_seed.iloc[0]) if len(sv) else np.nan,
            scramble_ratio=float(sv.scramble_ratio.iloc[0]) if len(sv) else np.nan,
            counts=counts))
    R = pd.DataFrame(rows)
    R.to_csv(out_csv, index=False)
    c_iso = int(R[(R.stratum == "isolated") & R.counts].shape[0]) if len(R) else 0
    c_30 = int(R.counts.sum()) if len(R) else 0
    fb = c_iso < 4
    bar = max(int(np.ceil(0.5 * c_30)), 5) if fb else W.ACCEPT_REDISCOVERED
    out = dict(
        n_seeded=int(len(R)),
        n_isolated_seeded=int((R.stratum == "isolated").sum()) if len(R) else 0,
        C_iso=c_iso, C_30=c_30,
        n_survived=int(R.survived.sum()) if len(R) else 0,
        n_located=int(R.located.sum()) if len(R) else 0,
        median_recall_all=float(np.nanmedian(R.recall_ref)) if len(R) else None,
        median_recall_indep=float(np.nanmedian(R.recall_indep_ref)) if len(R) else None,
        median_recall_loud=float(np.nanmedian(R.recall_loud_ref)) if len(R) else None,
        median_purity=float(np.nanmedian(R.purity_ref)) if len(R) else None,
        median_dist_km=float(np.nanmedian(R.dist_km)) if len(R) else None,
        median_quiet_day_indep=float(np.nanmedian(R.n_indep_quiet[R.n_indep_quiet >= 0]))
        if len(R) and (R.n_indep_quiet >= 0).any() else None,
        fallback_active=bool(fb),
        pass_bar=dict(mode=("all30_half_C30" if fb else "isolated_ge3"), value=int(bar),
                      note=("C_iso < 4: the isolated stratum is uninformative on this bed, so the "
                            "final bar becomes >= 0.5 x C_30 (minimum 5) of all 30 hidden families"
                            if fb else "C_iso >= 4: the isolated stratum is informative")),
        denominators=dict(
            i="all catalogued events in the 24-day ETS window -- REPORTED ONLY",
            ii=f"60-s-independent catalogued events excluding those within "
               f"+-{W.NEAR_SEED_S:.0f} s of either end of the seed window -- WHAT THE GATE COUNTS",
            iii=f"LOUD subset: network-median {Q.EX_WIN_S:.0f}-s S-energy excess >= "
                f"{Q.LOUD_EXCESS} at the family's {Q.LOUD_N_SITES} nearest admitted sites"),
        protocol=dict(seed="every catalogued event of the family in the ETS window, at the "
                           "STEP-1-CORRECTED time, spread over the whole 24-day window",
                      cap=f"the engine's own MAX_STACK = {W.MAX_STACK} (w_grow.py:88-89)",
                      near_seed="t_seed/t0/t1 = median/min/max of the spread seed, so the "
                                "near-seed burst guard cannot bind",
                      survival="the FROZEN gate with every clause, ETS-window days ONLY (F5); "
                               "quiet-day independent hits reported as the empirical FA",
                      chain=f"handoff -> quiet-day T_seed -> growth -> scramble twin -> survival "
                            f"-> {W.N_ITERATIONS} iterations -> dedupe -> locate/Gates 0-2"),
        accept=dict(C_30_ge=Q.S2_ACC_C30,
                    abandon="C_30 < 3 with Step-1 coherence verified -> STOP the study; the "
                            "write-up says 'this chain', not 'the archive'"),
        verdict="ACCEPT" if c_30 >= Q.S2_ACC_C30 else "ABANDON")
    Q.jdump(out, out_csv.replace(".csv", ".json"))
    Q.log(f"[s2] C_30 {c_30}/30, C_iso {c_iso}/10, fallback {fb}, bar {bar}; median recall "
          f"all/indep/loud "
          f"{out['median_recall_all']}/{out['median_recall_indep']}/{out['median_recall_loud']}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=f"{Q.OUT}/s2_oracle")
    ap.add_argument("--times", default=f"{Q.OUT}/corrected_times_2005.parquet")
    ap.add_argument("--out", default=f"{Q.OUT}/oracle_reference_2005_v3.csv")
    a = ap.parse_args()
    os.makedirs(a.dir, exist_ok=True)
    Sd = (pd.read_csv(f"{a.dir}/cl_seeds.csv") if os.path.exists(f"{a.dir}/cl_seeds.csv")
          else build_seeds(a.dir, a.times))
    if not len(Sd):
        raise SystemExit("no hidden family has >= 3 corrected events in the ETS window")
    X = loud_events(f"{Q.OUT}/s2_event_excess.parquet", Sd)
    r = QC.run(a.dir, f"{a.dir}/cl_seeds.csv", f"{a.dir}/cl_clusters.parquet", order="frozen")
    out = score(a.dir, Sd, X, a.out, r["final_tag"])
    if out["verdict"] != "ACCEPT":
        Q.marker(f"{Q.OUT}/S2_ABANDON.txt",
                 f"C_30 = {out['C_30']} (< {Q.S2_ACC_C30}) with Step-1 coherence verified: the "
                 f"frozen chain cannot pass with the best possible templates. STOP the study -- "
                 f"the write-up says 'this chain', not 'the archive'.\n"
                 f"C_iso {out['C_iso']}, survived {out['n_survived']}, located "
                 f"{out['n_located']}, median recall (all/indep/loud) "
                 f"{out['median_recall_all']}/{out['median_recall_indep']}/"
                 f"{out['median_recall_loud']}")
        return ABANDON_RC
    Q.marker(f"{Q.OUT}/S2_OK",
             f"ACCEPT: C_30 {out['C_30']}/30, C_iso {out['C_iso']}/10, fallback "
             f"{out['fallback_active']}, bar {out['pass_bar']['value']}; median recall "
             f"(all/indep/loud) {out['median_recall_all']}/{out['median_recall_indep']}/"
             f"{out['median_recall_loud']}; median purity {out['median_purity']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
