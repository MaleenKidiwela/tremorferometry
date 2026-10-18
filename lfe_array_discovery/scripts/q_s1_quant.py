#!/usr/bin/env python
"""FINAL SPEC, STEP 1 (CPU) -- QUANTISATION A/B and RE-ALIGNMENT.

Merlin's ruling (PLAN sec.12): the 2004-2005 rows of the Bostock export carry WHOLE-SECOND times
(100 % integer `sec`, 2.2-2.6 % in every other year) and `align_stack` cannot recover sub-second
alignment at per-event SNR < 1, so the 1c' oracle stacked NOISE.  This stage proves or kills that
before anything is rebuilt, and produces the corrected 2005 times the corrected oracle needs.

  (a) For the 5 hidden families with the most catalogued events in the 2005 bed (>= 300, and
      >= 100 catalogued 2012 events), at the shared permanents CN.PGC and PB.B003:
        * a PLAIN LINEAR 3-C stack of the 2012 events at t_origin + fam_tt  (sub-second times)
        * the same for the 2005 events                                      (quantised times)
      both through `z_common.stack_coherence` (1.0 = incoherent).
  (b) The 2012 stack is correlated against each 2005 event window within +- 1.0 s at those
      stations; the cc-peak lag is the CORRECTED 2005 time.
  (c) The 2005 events are restacked at the corrected times at the two nearest POLARIS sites.

ACCEPT   2012 coherence >= 2 at >= 1 shared station for >= 3 of the 5 families
         AND median 2005 quantised coherence <= 1.2
         AND corrected-lag spread (p90 - p10) >= 0.5 s
         AND corrected-time POLARIS coherence >= 2 at >= 2 distinct sites.
ABANDON  the best 2012 coherence anywhere is < 1.5 -> the stacking machinery or the time base is
         broken independent of quantisation (and the 2012 reference was itself obtained with
         incoherent templates): audit s_common.align_stack and fam_tt BEFORE any further step.
Anything else -> BELOW_ACCEPT, which the queue treats as a STOP ("report before Step 2", Merlin's
"Experiments I could not run": both coherent => F3 is wrong and the oracle failure is unexplained).

Also writes the corrected 2005 times for ALL 30 hidden families (needed by Step 2).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import z_common as Z                                                        # noqa: E402
import q_common as Q                                                        # noqa: E402

ABANDON_RC = 20
STOP_RC = 21


def pick_families(ev05, ev12, n=Q.S1_N_FAMILIES):
    """The n hidden families with the most 2005 in-window events among those with >= 300 in 2005
    AND >= 100 catalogued 2012 events.  Pre-registered, no discretion."""
    cand = [(len(v), int(k)) for k, v in ev05.items()
            if len(v) >= Q.S1_MIN_EVENTS_2005 and len(ev12.get(int(k), ())) >= Q.S1_MIN_EVENTS_2012]
    cand.sort(reverse=True)
    return [k for _n, k in cand[:n]]


def build_reference(cut12, bk, stas, ev12, max_events=400):
    """{sta: (stack, n, coherence)} of the 2012 sub-second-time linear stack."""
    FT = Q.fam_tt()
    out = {}
    for s in stas:
        tt = FT.get((int(bk), s))
        if tt is None:
            continue
        st, n, coh = Q.stack_at(cut12, s, ev12, tt, max_events=max_events)
        out[s] = (st, n, coh)
    return out


def corrected_times(cut05, bk, ref, ev05):
    """(corrected t_origin, per-event diagnostics) from the cc peak of the 2012 stack against each
    2005 event window within +- CC_LAG_S, summed over the usable shared stations."""
    FT = Q.fam_tt()
    lag = int(round(Q.CC_LAG_S * Q.FS))
    n = int(round(Q.QWIN_LEN_S * Q.FS))
    usable = [s for s, (st, nn, coh) in ref.items()
              if st is not None and np.isfinite(coh) and coh >= Q.S1_MIN_REF_COHERENCE]
    rows = []
    for tv in np.asarray(ev05, float):
        best_k, best_c, got = 0, np.nan, 0
        cc_sum = np.zeros(2 * lag + 1)
        for s in usable:
            tt = FT.get((int(bk), s))
            if tt is None:
                continue
            seg = cut05.cut(s, tv + tt - Q.QWIN_PRE_S - Q.CC_LAG_S,
                            Q.QWIN_LEN_S + 2 * Q.CC_LAG_S)
            if seg is None or seg.shape[1] < n + 2 * lag:
                continue
            r = ref[s][0]
            r = r - r.mean(1, keepdims=True)
            rn = np.sqrt((r ** 2).sum())
            if rn <= 0:
                continue
            for i, k in enumerate(range(-lag, lag + 1)):
                a = lag + k
                w = seg[:, a:a + n]
                w = w - w.mean(1, keepdims=True)
                wn = np.sqrt((w ** 2).sum())
                if wn > 0:
                    cc_sum[i] += float((r * w).sum() / (rn * wn))
            got += 1
        if got:
            i = int(np.argmax(cc_sum))
            best_k = i - lag
            best_c = float(cc_sum[i] / got)
        rows.append(dict(bkey=int(bk), t_origin=float(tv), n_sta=int(got),
                         lag_s=float(best_k) / Q.FS if got else np.nan,
                         cc_peak=best_c,
                         t_corrected=float(tv) + (float(best_k) / Q.FS if got else 0.0),
                         corrected=bool(got)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{Q.OUT}/s1_quantisation.csv")
    ap.add_argument("--times-out", default=f"{Q.OUT}/corrected_times_2005.parquet")
    ap.add_argument("--smoke", type=int, default=0, help="SMOKE: run only the first N families "
                                                         "and skip the all-30 correction")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    mem, of, order, det = W.sites()
    hidden, iso = Q.hidden_sets()
    ev05 = Q.events_2005(hidden)
    ev12 = Q.events_2012(hidden)
    fams = pick_families(ev05, ev12)
    if a.smoke:
        fams = fams[:a.smoke]
    Q.log(f"[s1] A/B families {fams} "
          f"(2005 {[len(ev05[k]) for k in fams]}, 2012 {[len(ev12[k]) for k in fams]})")

    shared = [s for s in Q.S1_SHARED if s in det]
    pol_stas = sorted({s for bk in fams for _q, s, _d in Q.polaris_sites(bk)})
    cut12 = Q.TwoYearCutter([2012], shared, maxn=24)
    cut05 = Q.TwoYearCutter([2005], sorted(set(shared) | set(pol_stas)), maxn=64)

    rows, tables = [], []
    for bk in fams:
        ref = build_reference(cut12, bk, shared, ev12[bk])
        # also build the POLARIS 2012 reference is NOT possible (POLARIS is gone in 2012) --
        # the corrected times come from the shared permanents alone, by construction.
        C = corrected_times(cut05, bk, ref, ev05[bk])
        tables.append(C)
        lags = C[C.corrected].lag_s.to_numpy(float)
        spread = float(np.percentile(lags, 90) - np.percentile(lags, 10)) if len(lags) else np.nan
        # PAIRED A/B: the quantised and the corrected stack use the SAME events, so the only
        # difference between the two coherences is the +- 1 s re-alignment.
        sub = C[C.corrected]
        if len(sub) > 400:
            sub = sub.iloc[np.linspace(0, len(sub) - 1, 400).round().astype(int)]
        tq = sub.t_origin.to_numpy(float) if len(sub) else np.asarray(ev05[bk], float)
        tc = sub.t_corrected.to_numpy(float) if len(sub) else np.zeros(0)
        FT = Q.fam_tt()
        for s in shared:
            st12, n12, c12 = ref.get(s, (None, 0, np.nan))
            tt = FT.get((int(bk), s))
            _st, n05, c05 = (Q.stack_at(cut05, s, tq, tt, max_events=400)
                             if tt is not None else (None, 0, np.nan))
            _stc, n05c, c05c = (Q.stack_at(cut05, s, tc, tt, max_events=400)
                                if (tt is not None and len(tc)) else (None, 0, np.nan))
            rows.append(dict(bkey=int(bk), sta=s, site=of.get(s, s), role="shared_permanent",
                             n_2012=n12, coherence_2012=c12,
                             n_2005=n05, coherence_2005_quantised=c05,
                             n_2005_corrected=n05c, coherence_2005_corrected=c05c,
                             lag_p10_s=float(np.percentile(lags, 10)) if len(lags) else np.nan,
                             lag_p90_s=float(np.percentile(lags, 90)) if len(lags) else np.nan,
                             lag_spread_s=spread,
                             lag_std_s=float(np.std(lags)) if len(lags) else np.nan,
                             median_cc_peak=float(np.nanmedian(C.cc_peak)) if len(C) else np.nan))
        for q, s, d in Q.polaris_sites(bk):
            tt = FT.get((int(bk), s))
            if tt is None:
                continue
            _st, n05, c05 = Q.stack_at(cut05, s, tq, tt, max_events=400)
            _stc, n05c, c05c = (Q.stack_at(cut05, s, tc, tt, max_events=400)
                                if len(tc) else (None, 0, np.nan))
            rows.append(dict(bkey=int(bk), sta=s, site=q, role="polaris", dist_km=d,
                             n_2012=0, coherence_2012=np.nan,
                             n_2005=n05, coherence_2005_quantised=c05,
                             n_2005_corrected=n05c, coherence_2005_corrected=c05c,
                             lag_spread_s=spread))
        Q.log(f"[s1] family {bk}: 2012 coh "
              f"{ {s: round(float(ref[s][2]), 3) for s in ref if np.isfinite(ref[s][2])} }; "
              f"lag spread {spread:.3f} s; corrected {int(C.corrected.sum())}/{len(C)}")

    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)

    # ---------------------------------------------------------------- verdict (mechanical)
    SH = R[R.role == "shared_permanent"]
    PO = R[R.role == "polaris"]
    fam_2012_ok = int(SH.groupby("bkey").coherence_2012.max().ge(Q.S1_ACC_COH_2012).sum())
    best_2012 = float(SH.coherence_2012.max()) if len(SH) else np.nan
    med_2005_q = float(SH.coherence_2005_quantised.median()) if len(SH) else np.nan
    spread = float(SH.groupby("bkey").lag_spread_s.first().median()) if len(SH) else np.nan
    pol_ok_sites = (PO.groupby("site").coherence_2005_corrected.median()
                    .ge(Q.S1_ACC_COH_POLARIS)) if len(PO) else pd.Series(dtype=bool)
    n_pol_sites = int(pol_ok_sites.sum())

    acc = dict(
        fam_2012_coherent=fam_2012_ok >= Q.S1_ACC_FAM_2012,
        quantised_incoherent=bool(np.isfinite(med_2005_q) and med_2005_q <= Q.S1_ACC_COH_2005_QUANT),
        lag_spread=bool(np.isfinite(spread) and spread >= Q.S1_ACC_LAG_SPREAD_S),
        polaris_corrected=n_pol_sites >= Q.S1_ACC_POLARIS_SITES)
    abandon = bool(not np.isfinite(best_2012) or best_2012 < Q.S1_ABANDON_COH_2012)
    verdict = "ABANDON" if abandon else ("ACCEPT" if all(acc.values()) else "BELOW_ACCEPT")

    out = dict(
        families=fams, shared_stations=shared,
        n_families_2012_coherent=fam_2012_ok, need=Q.S1_ACC_FAM_2012,
        best_2012_coherence=best_2012,
        median_2005_quantised_coherence=med_2005_q,
        median_lag_spread_p90_p10_s=spread,
        polaris_sites_corrected_ge2=n_pol_sites,
        polaris_site_median_corrected_coherence={k: float(v) for k, v in
                                                 (PO.groupby("site")
                                                  .coherence_2005_corrected.median().items()
                                                  if len(PO) else [])},
        accept_clauses=acc,
        accept=dict(coh_2012_ge=Q.S1_ACC_COH_2012, families_ge=Q.S1_ACC_FAM_2012,
                    coh_2005_quantised_le=Q.S1_ACC_COH_2005_QUANT,
                    lag_spread_ge_s=Q.S1_ACC_LAG_SPREAD_S,
                    polaris_coh_ge=Q.S1_ACC_COH_POLARIS,
                    polaris_sites_ge=Q.S1_ACC_POLARIS_SITES,
                    abandon_if_best_2012_lt=Q.S1_ABANDON_COH_2012),
        window=dict(pre_s=Q.QWIN_PRE_S, len_s=Q.QWIN_LEN_S, cc_lag_s=Q.CC_LAG_S,
                    ref_min_events=Q.S1_MIN_2012_FOR_CORRECTION,
                    ref_min_coherence=Q.S1_MIN_REF_COHERENCE),
        verdict=verdict,
        reason=("the best 2012 stack coherence anywhere is %.3f (< %.2f): the stacking machinery "
                "or the time base is broken independent of quantisation -- audit "
                "s_common.align_stack (489-505) and s_common.fam_tt BEFORE any further step"
                % (best_2012 if np.isfinite(best_2012) else -1.0, Q.S1_ABANDON_COH_2012))
        if abandon else ("" if verdict == "ACCEPT" else
                         "accept clauses not met: %s -- Merlin's 'Experiments I could not run' "
                         "requires this to be REPORTED before Step 2 (both coherent => F3 is "
                         "wrong and the 1c' oracle failure is unexplained)"
                         % {k: bool(v) for k, v in acc.items()}))
    Q.jdump(out, a.out.replace(".csv", ".json"))
    Q.log(f"[s1] VERDICT {verdict}: 2012 coherent {fam_2012_ok}/{len(fams)} families "
          f"(best {best_2012:.3f}); 2005 quantised median {med_2005_q:.3f}; lag spread "
          f"{spread:.3f} s; POLARIS sites corrected >= {Q.S1_ACC_COH_POLARIS}: {n_pol_sites}")

    if a.smoke:
        Q.log("[s1] SMOKE: corrected-time table for all 30 families NOT written")
        return 0
    if verdict == "ABANDON":
        Q.marker(f"{Q.OUT}/S1_ABANDON.txt", out["reason"])
        return ABANDON_RC
    if verdict != "ACCEPT":
        Q.marker(f"{Q.OUT}/S1_BELOW_ACCEPT.txt", out["reason"])
        return STOP_RC

    # ---------------------------------------------------------------- all 30 hidden families
    all_stas = sorted(set(shared))
    cut12b = Q.TwoYearCutter([2012], all_stas, maxn=24)
    cut05b = Q.TwoYearCutter([2005], all_stas, maxn=48)
    tabs = []
    for bk in sorted(hidden):
        e5 = ev05.get(int(bk))
        if e5 is None or not len(e5):
            continue
        e12 = ev12.get(int(bk), np.zeros(0))
        if len(e12) < Q.S1_MIN_2012_FOR_CORRECTION:
            tabs.append(pd.DataFrame(dict(bkey=int(bk), t_origin=e5, n_sta=0, lag_s=np.nan,
                                          cc_peak=np.nan, t_corrected=e5, corrected=False)))
            Q.log(f"[s1/all] family {bk}: only {len(e12)} 2012 events -- times left uncorrected")
            continue
        ref = build_reference(cut12b, bk, all_stas, e12)
        C = corrected_times(cut05b, bk, ref, e5)
        tabs.append(C)
        Q.log(f"[s1/all] family {bk}: corrected {int(C.corrected.sum())}/{len(C)} "
              f"(ref coherence { {s: round(float(ref[s][2]), 2) for s in ref if np.isfinite(ref[s][2])} })")
    T = pd.concat(tabs, ignore_index=True) if tabs else pd.DataFrame(
        columns=["bkey", "t_origin", "n_sta", "lag_s", "cc_peak", "t_corrected", "corrected"])
    T.to_parquet(a.times_out, index=False)
    Q.log(f"[s1] {len(T)} corrected 2005 times over {T.bkey.nunique()} hidden families -> "
          f"{a.times_out} ({int(T.corrected.sum())} actually re-aligned)")
    Q.marker(f"{Q.OUT}/S1_OK",
             f"ACCEPT: 2012 coherent {fam_2012_ok}/{len(fams)} families (best {best_2012:.3f}); "
             f"2005 quantised median {med_2005_q:.3f}; lag spread {spread:.3f} s; POLARIS sites "
             f"corrected >= {Q.S1_ACC_COH_POLARIS}: {n_pol_sites}; corrected times "
             f"{int(T.corrected.sum())}/{len(T)} over {T.bkey.nunique()} hidden families")
    return 0


if __name__ == "__main__":
    sys.exit(main())
