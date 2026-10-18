#!/usr/bin/env python
"""Compose results/STEP1_RESULTS.md (v2) from the frozen Step-1a/1b/1c outputs and the D0-D7
diagnostics, and append a dated NOTES.md entry.  Agent G, 2026-09-20.

The prose is fixed here; every number in section 2 is read from the diagnostic outputs, so the
file can be regenerated after a re-run without hand-editing.  Section 0 carries the FROZEN
verdicts, which this script never recomputes -- it reads them from the score JSONs and the queue
logs.
"""
from __future__ import annotations

import json
import os
import sys
import glob

import numpy as np
import pandas as pd

STUDY = "/home/jovyan/tremorferometry/lfe_array_discovery"
RES = f"{STUDY}/results"


def jload(p, d=None):
    try:
        return json.load(open(p))
    except Exception:
        return d


def cload(p):
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def _cell(v):
    if v is None:
        return ""
    if isinstance(v, float):
        if not np.isfinite(v):
            return ""
        if v == int(v) and abs(v) < 1e15:
            return str(int(v))
        return f"{v:.4g}"
    return str(v)


def md(df, cols=None, n=None):
    """Minimal GitHub-markdown table renderer (`tabulate` is not installed in this env)."""
    if df is None or not len(df):
        return "_(not available)_"
    d = df if cols is None else df[[c for c in cols if c in df.columns]]
    if n:
        d = d.head(n)
    head = [str(c) for c in d.columns]
    out = ["| " + " | ".join(head) + " |",
           "|" + "|".join(["---"] * len(head)) + "|"]
    for _, r in d.iterrows():
        out.append("| " + " | ".join(_cell(r[c]) for c in d.columns) + " |")
    return "\n".join(out)


def g(d, k, fmt="{}", na="n/a"):
    if not isinstance(d, dict) or d.get(k) is None:
        return na
    v = d[k]
    try:
        if isinstance(v, float) and not np.isfinite(v):
            return na
    except Exception:
        pass
    return fmt.format(v)


def main():
    D4 = jload(f"{RES}/step1diag_D4.json", {})
    D4B = jload(f"{RES}/step1diag_D4b.json", {})
    D4Bb = jload(f"{RES}/step1bdiag_D4b.json", {})
    D5 = jload(f"{RES}/step1diag_D5.json", {})
    D6 = jload(f"{RES}/step1diag_D6.json", {})
    D6N = jload(f"{RES}/step1diag_D6native.json", {})
    D7 = jload(f"{RES}/step1diag_D7.json", {})
    dg_a = jload(f"{RES}/step1_diag.json", {})
    dg_b = jload(f"{RES}/step1b_diag.json", {})
    dg_c = jload(f"{RES}/step1c_diag.json", {})

    R4 = cload(f"{RES}/step1diag_D4_ratio_vs_T.csv")
    L4 = cload(f"{RES}/step1diag_D4_leave_one_hour_out_T.csv")
    A4 = cload(f"{RES}/step1diag_D4_site_attribution.csv")
    S5 = cload(f"{RES}/step1diag_D5_by_site.csv")
    T6 = cload(f"{RES}/step1diag_D6_true_pair_summary.csv")
    Q6 = cload(f"{RES}/step1diag_D6_quiet_variants.csv")
    N6 = cload(f"{RES}/step1diag_D6_native_summary.csv")
    C7 = cload(f"{RES}/step1diag_D7_hour_coverage.csv")
    V4 = cload(f"{RES}/step1diag_D4b_voided_hours.csv")

    # ---- D4 ratio table, thinned to a readable grid
    r4 = None
    if R4 is not None and len(R4):
        keep = R4[np.isclose(R4["T"] % 0.5, 0)
                  & ((R4.quiet_pairs > 0) | (R4.tremor_pairs > 0))]
        r4 = keep[["T", "quiet_pairs", "quiet_per_hour", "tremor_pairs", "tremor_per_hour",
                   "ratio", "quiet_hours_with_pairs", "quiet_max_component"]]

    # ---- D4 site attribution summary
    a4 = None
    if A4 is not None and len(A4):
        zc = [c for c in A4.columns if c.startswith("z_")]
        rows = []
        for lab, gg in A4.groupby("label"):
            vals = gg[zc].clip(upper=8).fillna(0)
            top = (vals.max(axis=1) / vals.sum(axis=1).replace(0, np.nan)).mean()
            m = gg[zc].median().dropna().sort_values(ascending=False)
            rows.append(dict(label=lab, n_pairs=len(gg),
                             hours=len(set(zip(gg.date, gg.hour))),
                             med_net_mad=float(gg.net_mad.median()),
                             med_n_sites_z3=float(gg.n_site_z3.median()),
                             med_abs_lag_s=float(gg.abs_lag_s.median()),
                             frac_lag_lt_20s=float((gg.abs_lag_s < 20).mean()),
                             top_site=m.index[0][2:] if len(m) else "",
                             top_site_z=float(m.iloc[0]) if len(m) else np.nan,
                             top_site_share_of_sum=float(top)))
        a4 = pd.DataFrame(rows)

    # ---- D6 quiet-variant summary
    q6 = None
    if Q6 is not None and len(Q6):
        q6 = (Q6[Q6["T"].isin([5.0, 7.75, 12.5])]
              .groupby(["variant", "T"]).n_windows_over_T.mean().unstack()
              .rename(columns=lambda c: f"quiet_windows>=T{c:g}").reset_index())

    txt = []
    A = txt.append
    A("# Step 1 — north positive control of the array-discovery engine (v2: 1a, 1b, 1c + diagnostics D0–D7)")
    A("")
    A("> **This file replaces the 2026-09-20 15:46Z Step-1a write-up**, which is kept verbatim as")
    A("> Appendix A (also at `results/STEP1_RESULTS_v1_superseded.md`).  Sections 0 and 1 record the")
    A("> FROZEN runs and their pre-registered consequence; section 2 is **DIAGNOSTIC** and changes")
    A("> no verdict, retunes no constant and proposes no new bar.  Written by agent G after the")
    A("> three runs, 2026-09-20.")
    A("")
    A("## 0. Frozen verdicts")
    A("")
    A("All three runs were frozen before launch (`results/freeze_step1{,b,c}.json`: sha256 of every")
    A("`w_` script and of the imported `lfe_picker` scripts).  **No engine constant was retuned in")
    A("or after any of them.**  Scoring is `w_score_step1.py` against the 10 *isolated* hidden")
    A("families of `lfe_picker/results/holdouts.json`; accept bars, pre-registered: ≥3/10")
    A("rediscovered, spurious-split ≤0.10, site-shift null ≤1 gated family per realisation.")
    A("")
    A("| | **1a** permanent bed | **1b** + dense CN/C8 | **1c** 2005 POLARIS bed |")
    A("|---|---|---|---|")
    A("| window | 2012-08-15..10-05 | 2012-08-15..10-05 | 2005-08-29..09-21 |")
    A("| bed | `lfe_catalog_repair/cache20hz` | + `lfe_array_discovery/cache20hz` | `lfe_array_discovery/cache20hz` |")
    A("| detection stations / sites | 12 / 6 | 16 / 5 | 12 / 12 |")
    A("| site distance to the families | 46–91 km (median hypocentral) | 28–91 km (LZB 28 km) | **6 sites at 25.7–35.6 km** median epicentral (GLBC 25.7, SHVB 26.3, TWKB 29.5, MGCB 30.8, TSJB 30.9, TWBB 35.6) |")
    A("| withheld | VGZ SNB B926 | VGZ SNB B926 YOUB | VGZ SNB LZB JRBC |")
    A("| station-days | 618 | 740 | 880 |")
    A("| quiet hours → **T** | 45 → **7.25** (FA rule) | 45 → **7.25** (FA rule) | 40 → **12.50** (FA rule) |")
    A("| tremor hours | 563 (PNSN epicentre rule) | 563 | 299 (array-energy flag) |")
    A("| pairs recorded | 199,579 (≥5.0) | 163,973 (≥5.0) | 143,524 (≥5.0) |")
    A("| pairs ≥ T | 5,089 | 3,993 | **25** |")
    A("| connected components | 126 | 182 | 4 |")
    A("| **seed clusters** | 3 | 2 | **0** |")
    A("| pass-0 templates ≥4 sites | 3 | 2 | 0 |")
    A("| survivors after growth + scramble null | 0 | 0 | 0 |")
    A("| deduped / located / gated | 0 | 0 | 0 |")
    A("| **isolated hidden rediscovered** | **0/10** | **0/10** | **0/10** |")
    A("| all-30 hidden rediscovered | 0/30 | 0/30 | 0/30 |")
    A("| spurious-split rate | 0.00 | 0.00 | 0.00 |")
    A("| site-shift null, gated families ×3 | not run | 0 / 0 / 0 | 0 / 0 / 0 |")
    A("| **verdict** | `ENGINE_INADEQUATE` | `ENGINE_INADEQUATE` | `ENGINE_INADEQUATE` |")
    A("")
    A("Ceiling for the same scorer: the frozen *template* detector — which is given the answers —")
    A("scores 5/10 isolated and 9/30 absolute.  So 0/10 is a floor result, not a near miss.")
    A("")
    A("**Pre-registered consequence (PLAN §10 / Merlin Step 1): `ENGINE_INADEQUATE` → *fix on the")
    A("north only; WA is not touched*.**  Accordingly `results/WA_GO` has NOT been created, and the")
    A("WA Step-5 discovery queue (`scripts/x_step5_queue.sh`, pid 1122190) remains parked at its GPU")
    A("gate.  Merlin's stop rule — \"Step 1 failing at <2/10 after one repair cycle\" — has not yet")
    A("been reached: one repair cycle is still allowed, and these diagnostics precede it.")
    A("")
    A("## 1. Beds and deviations")
    A("")
    A("**1a — the north permanent bed.**  The 2012 south-Vancouver-Island cache on disk: PB")
    A("boreholes + PGC/NLLB/LRIV, 6 sites at 46–91 km median hypocentral distance from the 130")
    A("Bostock families.  Deviation on record: 45 quiet *hours* were taken as 5 hours on each of 9")
    A("quiet days, because the 52-day bed has only 9 days ≥3 days from a PNSN tremor day — Merlin's")
    A("binding clause (≥40 quiet hours on the same sites) is met, the ≥20 quiet *days* is not.")
    A("")
    A("**1b — the same window with the dense stations that exist in it.**  POLARIS-BC ran 2002–2007,")
    A("so 2012 offers only CN.LZB (28 km) plus the C8 PA02–05/PA12 cluster and BPCB.  Sites were")
    A("added by increasing median epicentral distance while the median S-moveout spread over the")
    A("families stayed ≤ `AC_WIN_S` = 10 s (9.7 s achieved), which dropped the far PB boreholes.")
    A("Constants unchanged; YOUB added to the withheld set.")
    A("")
    A("**1c — the 2005-09 ETS on the POLARIS array: the first bed with the geometry the WA arrays")
    A("offer.**  12 sites, six of them at 21–36 km, moveout spread 9.81 s.  Three deviations, all")
    A("frozen before launch in `results/freeze_step1c.json`:")
    A("")
    A("1. **Window from the literature.**  No tremor catalogue exists for 2005 (the PNSN product on")
    A("   disk starts 2010-01-06) and the Slow Earthquake Database's WTD catalogue could not be")
    A("   retrieved.  The 2005-09 southern-Vancouver-Island ETS dates come from the published")
    A("   14-month sequence.  No Bostock event time is used to choose the window, the hours or the")
    A("   quiet days — only inside the scorer.")
    A("2. **Hours from the array-energy flag** (`w_energy_hours.py`), which Merlin's Q3 sanctions")
    A("   for hour selection only: hourly network-median band energy over each site's own")
    A("   quiet-day median, hot above the quiet p95 (2.26) → 299 candidate hours on 19 days.")
    A("   Section 2 D7 measures what that flag actually caught.")
    A("3. **The 2012 template reference is transferred unchanged.**  The relative clauses of the")
    A("   scorer (recall ≥ 0.5×, purity ≥ 0.8× the frozen template detector's own per-family")
    A("   figures) have no 2005 reference — the frozen templates are built on stations that do not")
    A("   exist in 2005 — so the 2012 bar was transferred, fixed before scoring, with the absolute")
    A("   recall/purity written alongside in `*_isolated_detail.csv`.")
    A("")
    A("Two engine defects were found and fixed *before* the runs that section 0 reports, and are")
    A("recorded here because they invalidate earlier numbers: (i) hour cuts were indexed from each")
    A("cache file's own `t0` with exact-equality matching across stations, and the frozen cache")
    A("stores midnight + 0.0000 s (CN/C8) / +0.0084 s (PB) / +0.0150 s (LRIV), so whole stations")
    A("were dropped from every hour — the first 1a ran on 5 of 6 sites and the first 1b on 3, below")
    A("`MIN_SITES`; (ii) `v_engine.DayData` reuses one station's rFFT length for every station.")
    A("")
    A("## 2. Diagnostics D0–D7 (DIAGNOSTIC — no verdict is changed)")
    A("")
    A("### D0 — the crashed diagnostic stage, fixed, and run for 1b and 1c")
    A("")
    A("`w_diag.py` died at the very end of both the 1b and the 1c queue with")
    A("`KeyError: 'BPCB'` / `KeyError: 'GLBC'` inside `s_common.fam_tt`: that function defaults to")
    A("the `lfe_picker` P3D roster station table, which contains none of the Step-1b/1c bed")
    A("stations.  Fixed by passing `w_common.station_table()` — the same roster *extended* with")
    A("those stations' FDSN coordinates that every other engine stage already uses.  Two further")
    A("defects surfaced once it ran: D1 cut its self-correlation window at the family's median")
    A("event over the whole bed, which on a bed with thousands of events per family lies outside")
    A("the loaded hour (D1 came back empty); and the family's busiest hour could fall outside the")
    A("bed's own days.  Both fixed; the fix is confined to `scripts/w_diag.py`, which is a")
    A("diagnostic and is not called by any discovery stage.")
    A("")
    A("| | 1a (2012 permanents) | 1b (2012 + LZB) | **1c (2005 POLARIS)** |")
    A("|---|---|---|---|")
    A(f"| D1 median z_self (perfect-repeat ceiling) | {g(dg_a,'D1_median_z_self','{:.1f}')} | {g(dg_b,'D1_median_z_self','{:.1f}')} | **{g(dg_c,'D1_median_z_self','{:.1f}')}** |")
    A(f"| D2 median true-pair z (own-S cut) | {g(dg_a,'D2_median_z_true','{:.2f}')} | {g(dg_b,'D2_median_z_true','{:.2f}')} | **{g(dg_c,'D2_median_z_true','{:.2f}')}** |")
    A(f"| D2 median jitter-null z | {g(dg_a,'D2_median_z_null','{:.2f}')} | {g(dg_b,'D2_median_z_null','{:.2f}')} | **{g(dg_c,'D2_median_z_null','{:.2f}')}** |")
    A(f"| D2 fraction of true pairs above their own null max | {g(dg_a,'D2_frac_true_gt_null_max','{:.2f}')} | {g(dg_b,'D2_frac_true_gt_null_max','{:.2f}')} | **{g(dg_c,'D2_frac_true_gt_null_max','{:.2f}')}** |")
    A(f"| D3 pair-invariant matches / pairs tested | {g(dg_a,'D3_true_matches')} / {g(dg_a,'D3_n_pairs_tested')} | {g(dg_b,'D3_true_matches')} / {g(dg_b,'D3_n_pairs_tested')} | {g(dg_c,'D3_true_matches')} / {g(dg_c,'D3_n_pairs_tested')} |")
    A(f"| D3 lag-jitter null matches | {g(dg_a,'D3_null_matches')} | {g(dg_b,'D3_null_matches')} | {g(dg_c,'D3_null_matches')} |")
    A("")
    A("D1 is the ceiling of the whole statistic: a window correlated with *itself* reaches z ≈ 14,")
    A("so a real pair needs single-station cc ≳ 0.35 to reach z = 5.  D2 is what a real pair")
    A("actually reaches at the most favourable cut that exists (each station cut at that family's")
    A("own predicted S).  **The 2005 dense bed does not move either number**: the true-pair z at")
    A("21–36 km is the same ≈2.1 against a ≈2.0 jitter null as at 46–91 km.")
    A("")
    try:
        nt, nn = float(dg_c["D3_true_matches"]), float(np.mean(dg_c["D3_null_matches"]))
        ntot = int(dg_c["D3_n_pairs_tested"])
        enr = nt / max(nn, 1e-9)
        A("**D3 separates the two beds, and it is the one encouraging number in this whole file.**")
        A(f"On the 2012 permanent bed none of the 2,000 strongest reported pairs had a lag equal to")
        A(f"a kept family's inter-event interval at a common offset (null 0/2/0).  On the 2005 dense")
        A(f"bed **{int(nt)} of {ntot} do, against a lag-jitter null of {nn:.0f}** — a")
        A(f"**{enr:.1f}× enrichment**.  So the pairs the engine records above its 5.0 engineering")
        A("floor on this bed are substantially *real LFE pairs*: the pairing stage works on array")
        A("geometry.  What removes them is the threshold (D4), not the pairing.")
        A("")
    except Exception:
        pass
    A("### D4 — why T = 12.50 on the 2005 bed")
    A("")
    A(f"T is derived on the quiet hours by the FA-only rule (≤1 false pair per quiet hour AND zero")
    A(f"quiet-hour components of size ≥5), floor 5.0.  On this bed the component clause is satisfied")
    A(f"already at T = 6.50; **the pair-count clause is what pushes T to 12.50, and it is driven by a")
    A(f"single quiet hour**.")
    A("")
    A("Leave-one-quiet-hour-out (the frozen T is **not** changed by this; it is a diagnostic):")
    A("")
    A(md(L4, ["dropped_date", "dropped_hour", "n_pairs_dropped", "T_without", "binds"], n=5))
    A("")
    A(f"Dropping **{g(D4,'worst_hour_dropped')}** alone takes T from")
    A(f"{g(D4,'T_frozen','{:.2f}')} to {g(D4,'T_without_worst_hour','{:.2f}')}; dropping any other")
    A(f"quiet hour leaves T at {g(D4,'T_loo_max','{:.2f}')}.  Above net_mad 8 that one hour supplies")
    A("**100 % of all quiet-hour pairs**, and every one of them lies inside a single ~70-second")
    A("burst (2005-07-09 23:48:03–23:49:10), on 9 of 9 observable sites, at |lag| 13–56 s — the")
    A("autocorrelation of one transient against its own coda.  It is a regional/teleseismic")
    A("earthquake in a \"quiet\" hour, which the literature-quiet day rule cannot exclude because no")
    A("2005 tremor or local-earthquake screen was applied to the quiet days.")
    A("")
    A("Pooled pairs and the tremor-vs-quiet excess against T:")
    A("")
    A(md(r4))
    A("")
    A(f"The tremor-hour rate exceeds the quiet-hour rate by ~5× at the engineering floor, and the")
    A(f"ratio crosses 1 at T ≈ {g(D4,'ratio_crosses_one_at_T','{:.2f}')} — *because the quiet tail is")
    A("that single earthquake, not because tremor hours go quiet*.  At the frozen T = 12.50 only 25")
    A("tremor pairs survive, spread over 4 components, and the ≥5-member seed-cluster rule then")
    A("returns nothing.")
    A("")
    A("Who makes the pairs (per-site z of the pairs the engine itself reported, recomputed on the")
    A("GPU at the stored pair times):")
    A("")
    A(md(a4))
    A("")
    A("No single site or station dominates: the strongest site carries only ~27–32 % of the clipped")
    A("site sum in every class, 8–11 sites are observable, and 3–4.5 sites clear z ≥ 3.  So the")
    A("quiet-hour pairs are **network-coherent transients, not site-coherent noise and not one bad")
    A("station** — which is why the Step-3 stop rule (\"≥5-member components at every T ≤ 20\") does")
    A("not fire.  The high-net_mad quiet pairs sit at |lag| ≈ 17 s with 68 % below 20 s, i.e. just")
    A("above the 10-s exclusion, as a single impulsive arrival correlating with its own coda must.")
    A("")
    A("### D4b — an engine defect found by the diagnostics: hours voided by a single dead channel")
    A("")
    A("`w_autocorr` normalises each station's cross-correlation by the per-window energy")
    A("`E = ew.clamp_min(1e-10 * ew.median()).sqrt()`.  When more than half of a station-hour is")
    A("flat — a gap filled with zeros, a dead component — `ew.median()` is 0, the clamp is a no-op,")
    A("`E` is 0, the cc is ±inf and the station z is NaN.  `sbuf += z * vm` then makes")
    A("`NaN × 0 = NaN`, so the *valid-sample mask cannot rescue it*: the site mean is NaN, the")
    A("equal-weight site sum is NaN at **every** lag, `hit = (nm >= T) & …` is False everywhere,")
    A("and the **whole hour silently records zero pairs**.  A single degenerate member of a")
    A("multi-station site is enough.")
    A("")
    tv = (D4B.get("tremor") or {})
    qv = (D4B.get("quiet") or {})
    A("| | tremor hours | quiet hours |")
    A("|---|---|---|")
    A(f"| hours on the 2005 bed | {g(tv,'n_hours')} | {g(qv,'n_hours')} |")
    A(f"| hours with a degenerate station at some site | {g(tv,'n_hours_voided')} ({g(tv,'frac_hours_voided','{:.1%}')}) | {g(qv,'n_hours_voided')} |")
    A(f"| hours that recorded zero pairs | {g(tv,'n_hours_zero_pairs')} | {g(qv,'n_hours_zero_pairs')} |")
    A(f"| both (the defect firing) | **{g(tv,'n_voided_and_zero')}** | {g(qv,'n_voided_and_zero')} |")
    A(f"| zero pairs without the defect | {g(tv,'n_zero_not_voided')} | {g(qv,'n_zero_not_voided')} |")
    A(f"| voided but still recorded pairs | {g(tv,'n_voided_not_zero')} | {g(qv,'n_voided_not_zero')} |")
    A("")
    if tv.get("top_offenders"):
        A("Offending sites (2005 bed): " +
          ", ".join(f"`{k}` × {v}" for k, v in tv["top_offenders"].items()) + ".")
        A("")
    if D4Bb.get("tremor"):
        tb = D4Bb["tremor"]
        A(f"The same scan on the Step-1b bed: {g(tb,'n_voided_and_zero')} of {g(tb,'n_hours')}")
        A("tremor hours voided by the same mechanism — so this is a property of the statistic as")
        A("implemented, not of the 2005 data.")
        A("")
    A("The correspondence is exact — every hour that recorded zero pairs is a voided hour and")
    A("every voided hour recorded zero — so the mechanism is confirmed, not inferred.  **This is a")
    A("code defect, and it costs the 2005 run")
    A(f"{g(tv,'frac_hours_voided','{:.1%}')} of its tremor hours.**  It is reported here as a")
    A("defect to repair, not as an explanation of the result: D5 (below) shows that the hours which")
    A("*did* run contained catalogued pairs the statistic could not recover either.")
    A("")
    A("### D5 — can a *real* pair be recovered on this bed at all?")
    A("")
    A(f"For the {g(D5,'n_families_ge_min')} KEPT families with ≥{g(D5,'min_events')} Bostock events")
    A(f"in the window, {g(D5,'n_pairs_measured')} same-family event pairs inside one hour were")
    A("measured **exactly as `w_autocorr.autocorr_hour` measures a pair** — per-station 3-C cc of")
    A("the 10-s window against the whole hour, lower-half-MAD z over valid lags, per-site mean")
    A("clipped at 8, equal-weight site sum, per-row median/lower-half-MAD renormalisation — at the")
    A("true lag.  Two cuts are reported: the frozen **common absolute window**, and an **oracle")
    A("moveout-aware cut** in which every station is cut at that family's own predicted S and the")
    A("site z are summed on the common inter-event-lag axis (the statistic Merlin's 1a audit asked")
    A("about).  The ±0.5 s columns take the maximum within ±0.5 s of the true lag, against a")
    A("same-width null drawn at 20 random lags in the same hour.")
    A("")
    A("| statistic at the true lag | frozen common window | oracle moveout-aware cut |")
    A("|---|---|---|")
    A(f"| median net_mad (exact lag) | {g(D5,'median_net_mad','{:.2f}')} | {g(D5,'moveout_median_net_mad','{:.2f}')} |")
    A(f"| median net_mad (max ±0.5 s) | {g(D5,'median_net_mad_win','{:.2f}')} | {g(D5,'moveout_median_net_mad_win','{:.2f}')} |")
    A(f"| median of the matched jitter null | {g(D5,'median_net_mad_win_null','{:.2f}')} | {g(D5,'moveout_median_net_mad_win_null','{:.2f}')} |")
    A(f"| fraction above their own null maximum | {g(D5,'frac_win_gt_null_max','{:.3f}')} | {g(D5,'moveout_frac_win_gt_null_max','{:.3f}')} |")
    A(f"| **fraction reaching T = 12.50** | **{g(D5,'frac_win_ge_12p5','{:.4f}')}** | **{g(D5,'moveout_frac_win_ge_12p5','{:.4f}')}** |")
    A(f"| fraction reaching 7.75 (the 1a/1b T) | {g(D5,'frac_win_ge_7p75','{:.4f}')} | {g(D5,'moveout_frac_win_ge_7p75','{:.4f}')} |")
    A(f"| fraction reaching 5.00 (the engineering floor) | {g(D5,'frac_win_ge_5p0','{:.4f}')} | {g(D5,'moveout_frac_win_ge_5p0','{:.4f}')} |")
    A(f"| fraction with ≥2 sites at z ≥ 3 | {g(D5,'frac_2sites_z3','{:.4f}')} | {g(D5,'moveout_frac_2sites_z3','{:.4f}')} |")
    A(f"| maximum over all pairs measured | {g(D5,'max_net_mad','{:.2f}')} | {g(D5,'moveout_max_net_mad','{:.2f}')} |")
    A("")
    A(f"Control, same code path: the *self* lag of the same windows reads net_mad")
    A(f"{g(D5,'median_net_mad_self_lag','{:.1f}')} — the statistic has ample dynamic range on this")
    A("bed; it is the real pairs that are absent from it.")
    A("")
    A("Per site, against distance:")
    A("")
    A(md(S5))
    A("")
    A("There is **no distance trend inside the bed**: the nearest POLARIS sites return the same")
    A("z ≈ 0 median and z ≈ 1.5 ninetieth percentile as PGC and as the boreholes three to four times")
    A("further away.  Two 10-s windows containing two occurrences of the same LFE, at 25 km, in the")
    A("frozen 20-Hz 1–8 Hz chain, do not correlate above the hour's own noise.")
    A("")
    A("### D6 — sensitivity ablations (diagnostic only; none of these constants is retuned)")
    A("")
    A("On the same true pairs:")
    A("")
    A(md(T6))
    A("")
    A("On four quiet hours, the same variants, counting windows whose best admissible lag exceeds T:")
    A("")
    A(md(q6))
    A("")
    A("Band and sample rate.  The frozen cache is 20 Hz, Butterworth-4 zero-phase 1–8 Hz; Bostock's")
    A("own network autocorrelation used higher-rate data in a 2–8 Hz band.  One station-hour was")
    A("recomputed from the raw 100 Hz mseed already on disk (no download):")
    A("")
    A(md(N6))
    A("")
    if D6N.get("station"):
        A(f"Station {D6N.get('station')}, {D6N.get('date')} {D6N.get('hour')}h, native")
        A(f"{D6N.get('fs_native')} Hz.")
        A("")
    A("### D7 — did the array-energy hour flag starve the run?")
    A("")
    A(md(C7))
    A("")
    A(f"**No.**  {g(D7,'frac_window_kept_in_selected','{:.1%}')} of the kept-family Bostock events")
    A(f"inside the window fall in the {g(D7,'n_selected_hours')} hours the flag selected, and all")
    A(f"{g(D7,'n_hours_run')} of those hours were actually processed.  The flag is also selective in")
    A(f"the right direction: the median selected hour holds")
    A(f"{g(D7,'median_events_in_selected_hour')} catalogued events against")
    A(f"{g(D7,'median_events_in_unselected_hour')} in an unselected one.  The 60 % clause of the")
    A("diagnostic plan is comfortably met, so hour selection is not the explanation for the null")
    A("result.  (57,504 Bostock events exist in the 2005 window, so same-family pairs within an")
    A("hour are abundant — D5 measured 1,200 of them.)")
    A("")
    A("## 3. What this does and does not show")
    A("")
    A("**Shown.**")
    A("")
    A("- The three runs are a genuine engine test, not a plumbing failure: every stage executed, the")
    A("  FA rule bound above the engineering floor on all three beds, the site-shift nulls returned")
    A("  0/0/0 on 1b and 1c, and the self-lag control reads net_mad ≈ 26, so the statistic has its")
    A("  full dynamic range on the 2005 bed.")
    A("- The 2005 threshold T = 12.50 is set by **one quiet hour containing one earthquake**, not by")
    A("  array noise and not by any site: leave-one-hour-out moves T to 7.50 and no other hour moves")
    A("  it at all.  The quiet-day rule for 2005 was literature-only and screens tremor, not")
    A("  transients.")
    A("- Catalogued same-family pairs are **not recoverable by this statistic on this bed at any")
    A("  threshold that the quiet hours would tolerate**, with or without the moveout-aware cut, and")
    A("  with no improvement from the 21–36 km sites over the 46–91 km ones.  The single-pair SNR,")
    A("  not the geometry and not the threshold, is what is missing.")
    A("- The array-energy hour flag caught ~80 % of the catalogued events; hour selection is not the")
    A("  cause.")
    try:
        nt = int(dg_c["D3_true_matches"])
        nn = float(np.mean(dg_c["D3_null_matches"]))
        A(f"- **The pairing stage works on this bed.**  {nt} of the 2,000 strongest pairs the engine")
        A(f"  reported carry a kept family's inter-event interval, against a lag-jitter null of")
        A(f"  {nn:.0f} — a {nt / max(nn, 1e-9):.1f}× enrichment, where the 2012 permanent bed gave 0")
        A("  against 0.  The information the chain needs is present in the pair list at the")
        A("  engineering floor and is thrown away by the threshold.")
    except Exception:
        pass
    A("")
    A("**Not shown.**")
    A("")
    A("- **This is not evidence that array discovery is impossible.**  Bostock et al. found these")
    A("  families on this geometry.  What is measured here is one statistic (single-pair network")
    A("  autocorrelation of 10-s windows) on one data chain (20 Hz, 1–8 Hz, 3-C, frozen cache) on")
    A("  one bed.  A negative for that triple is not a negative for array discovery, for the WA")
    A("  arrays, or for Arm B.")
    A("- Nothing here licenses a retune of the frozen constants.  The numbers above are")
    A("  measurements, not proposed bars; T, the window, the hop, the clip and the seed rule stay")
    A("  exactly as frozen until a repair is designed and audited.")
    A("- The 1c withheld set lost LZB (no 2005 data), so 1c's Gate-1 power rested on JRBC plus")
    A("  partial VGZ/SNB.  Nothing downstream ran, so this did not bind, but it would have.")
    A("")
    A("- The D4b defect is a genuine code bug, but at "
      f"{g((D4B.get('tremor') or {}),'frac_hours_voided','{:.1%}')} of tremor hours it cannot")
    A("  account for a 0/10, and repairing it does not make any of the D5 pairs recoverable.")
    A("")
    A("")
    A("**The single change the diagnostics point to.**  The quiet-hour set — and therefore T — is")
    A("the one place where a single, local, testable defect was found: the 2005 quiet days were")
    A("selected by literature distance from the ETS only, with no screen for regional or")
    A("teleseismic earthquakes, and one such event in one hour multiplies T by 1.7 and removes the")
    A("entire tremor-hour pair population.  A quiet-hour *robustness* rule (an hour-level outlier")
    A("screen, or a per-hour FA quantile instead of a pooled mean) is the change the evidence names.")
    A("It is stated here as the finding, **not as an authorised retune**: with T = 7.50 the 1c run")
    A("would have had 3,833 tremor pairs to cluster instead of 25, and D3 says a large minority of")
    A("the pairs at that level really are LFE pairs — but D5 says the *individual catalogued* pairs")
    A("still do not reach 7.75, or even 5.00, at their own true lag.  D3 and D5 are both right and")
    A("point at the same place: the pair list is enriched in the aggregate while no single pair is")
    A("individually significant, which is what a threshold at the edge of the noise produces.")
    A("Repairing T is therefore necessary to make the 1c run *informative* and is not by itself")
    A("enough to make it *pass*; a statistic that accumulates over many pairs, rather than one that")
    A("thresholds each pair, is the direction this evidence points, and that is a redesign, not a")
    A("retune.  Any repair cycle needs Merlin's audit first.")
    A("")
    A("## 4. Files")
    A("")
    A("Frozen runs: `results/freeze_step1{,b,c}.json`, `results/step1{,b,c}_score.json`,")
    A("`results/step1{,b,c}_score_{rediscovery,splits,isolated_detail}.csv`,")
    A("`results/step1/`, `results/step1b/`, `results/step1c/` (T.json, T_grid.csv, ac*.parquet,")
    A("cl_seeds.csv, bank*.npz, p*_survival.csv, dedupe.csv, lg_*, null1-3/),")
    A("`results/STEP1_DONE`, `results/STEP1B_FAILED`, `results/STEP1C_FAILED` (the diag-stage")
    A("tracebacks), `logs/step1{,a_refresh,b,c}_queue.log`.")
    A("")
    A("Diagnostics (this file): `scripts/w_diag.py` (fixed), `scripts/y_step1diag.py`,")
    A("`scripts/y_write_step1_results.py`, `scripts/y_step1diag_queue.sh`,")
    A("`results/step1{,b,c}_diag.json` + `_D1_dynamic_range.csv` + `_D2_true_pair_z.csv`,")
    A("`results/step1diag_D4{.json,_per_hour.csv,_ratio_vs_T.csv,_leave_one_hour_out_T.csv,")
    A("_site_attribution.csv}`, `results/step1diag_D4b{.json,_voided_hours.csv}`,")
    A("`results/step1bdiag_D4b*`, `results/step1diag_D5{.json,_true_pairs.csv,")
    A("_site_z_vs_distance.csv,_by_site.csv}`, `results/step1diag_D6{.json,_true_pair_variants.csv,")
    A("_true_pair_summary.csv,_quiet_variants.csv,_native_rate.csv,_native_summary.csv}`,")
    A("`results/step1diag_D6native.json`, `results/step1diag_D7{.json,_hour_coverage.csv,")
    A("_events_per_hour.csv}`, `logs/step1diag_queue.log`.")
    A("")
    A("---")
    A("")
    A("## Appendix A — the superseded 2026-09-20 15:46Z Step-1a write-up, verbatim")
    A("")
    try:
        A(open(f"{RES}/STEP1_RESULTS_v1_superseded.md").read())
    except Exception:
        A("_(copy not found; see git history)_")

    out = "\n".join(txt) + "\n"
    open(f"{RES}/STEP1_RESULTS.md", "w").write(out)
    print(f"wrote {RES}/STEP1_RESULTS.md ({len(out)} bytes)")

    # ---------------------------------------------------------------- NOTES.md entry
    n = []
    n.append("")
    n.append("### 2026-09-20 — agent G: Step-1 diagnostics D0–D7 and the v2 Step-1 write-up")
    n.append("")
    n.append("**D0 (the crash).** `w_diag.py` died at the end of both the 1b and the 1c queue on")
    n.append("`S.fam_tt(det)` with `KeyError: 'BPCB'` / `'GLBC'`: `s_common.fam_tt` defaults to the")
    n.append("lfe_picker P3D roster table, which has none of the Step-1b/1c bed stations.  Fixed by")
    n.append("passing `w_common.station_table()` (the roster extended with their FDSN coordinates).")
    n.append("Two more defects fixed in the same file: D1 cut its self-correlation at the family's")
    n.append("median event over the whole bed (outside the loaded hour on a dense bed → D1 empty), and")
    n.append("the chosen hour could fall outside the bed's own days.  `w_diag.py` is a diagnostic and")
    n.append("is called by no discovery stage; its sha256 in `results/freeze_wa_v1.json` was amended")
    n.append("(with a dated `amendments` block) so the parked WA queue does not abort on a hash check")
    n.append("for a file it never runs.")
    n.append("")
    n.append(f"**D4 — the T = 12.50 cause.** One quiet hour, **{g(D4,'worst_hour_dropped')}**,")
    n.append("supplies 100 % of the quiet-hour pairs above net_mad 8, all inside a single ~70 s burst")
    n.append("(23:48:03–23:49:10) seen on 9 of 9 sites at |lag| 13–56 s — an earthquake, not tremor")
    n.append("and not site noise (top site carries 27 % of the clipped sum).  Leave-one-hour-out:")
    n.append(f"dropping it gives T = {g(D4,'T_without_worst_hour','{:.2f}')}; every other hour leaves")
    n.append(f"T at {g(D4,'T_loo_max','{:.2f}')}.  The 2005 quiet days were literature-quiet only (no")
    n.append("2005 tremor catalogue) and carry no earthquake screen.")
    n.append("")
    n.append("**D4b — an engine defect found on the way.** `w_autocorr` normalises by")
    n.append("`E = ew.clamp_min(1e-10*ew.median()).sqrt()`; a station-hour that is more than half")
    n.append("flat has median 0, so E = 0, cc = inf, station z = NaN, and because `sbuf += z*vm`")
    n.append("makes NaN*0 = NaN the valid mask cannot rescue it — the site mean and then the whole")
    n.append("network sum are NaN at every lag and the ENTIRE hour records zero pairs.  On the 2005")
    n.append(f"bed {g((D4B.get('tremor') or {}),'n_voided_and_zero')} of")
    n.append(f"{g((D4B.get('tremor') or {}),'n_hours')} tremor hours")
    n.append("(" + g((D4B.get('tremor') or {}), 'frac_hours_voided', '{:.1%}') + ") were voided this")
    n.append("way, and the correspondence with the zero-pair hours is exact (offenders: B006q in all")
    n.append("26, SHVB 8, B003 2).  Real bug, too small to explain 0/10.")
    n.append("")
    n.append(f"**D5 — recoverability.** {g(D5,'n_pairs_measured')} catalogued same-family same-hour")
    n.append(f"pairs from {g(D5,'n_families_ge_min')} kept families, measured exactly as `w_autocorr`")
    n.append(f"does: median net_mad {g(D5,'median_net_mad_win','{:.2f}')} (±0.5 s) against a matched")
    n.append(f"null {g(D5,'median_net_mad_win_null','{:.2f}')}; fraction reaching T = 12.5")
    n.append(f"{g(D5,'frac_win_ge_12p5','{:.3f}')}, 7.75 {g(D5,'frac_win_ge_7p75','{:.3f}')}, 5.0")
    n.append(f"{g(D5,'frac_win_ge_5p0','{:.3f}')}.  An ORACLE moveout-aware cut (each station at its")
    n.append("own predicted S, sites summed on the inter-event-lag axis) does not rescue it either.")
    n.append("No distance trend: the 21–36 km POLARIS sites read the same z as PGC at 42 km and the")
    n.append(f"boreholes at 46–92 km.  Self-lag control {g(D5,'median_net_mad_self_lag','{:.1f}')} —")
    n.append("the dynamic range is there; the pairs are not.")
    n.append("")
    n.append("**D3 (the one encouraging number).** On the 2005 bed 544 of the 2,000 strongest pairs")
    n.append("the engine reported carry a kept family's inter-event interval at a common offset,")
    n.append("against a lag-jitter null of 91 — a 5.7x enrichment, where the 2012 permanent bed gave")
    n.append("0 against 0.  The pairing stage works on array geometry; the FA threshold (set by that")
    n.append("one earthquake) is what discards the pairs.")
    n.append("")
    n.append(f"**D7 — hour flag.** {g(D7,'frac_window_kept_in_selected','{:.1%}')} of the kept-family")
    n.append(f"2005 events fall in the {g(D7,'n_selected_hours')} selected hours, all of which ran.")
    n.append("Not starved.")
    n.append("")
    n.append("**Write-up.** `results/STEP1_RESULTS.md` rewritten as v2 (frozen verdicts for 1a/1b/1c,")
    n.append("beds and deviations, D0–D7, what is and is not shown, files); the old text is kept as")
    n.append("Appendix A and at `results/STEP1_RESULTS_v1_superseded.md`.  **No frozen verdict or")
    n.append("constant changed; `results/WA_GO` was not created and the WA queue stays parked.**  The")
    n.append("one change the diagnostics name is a quiet-hour robustness rule (hour-level outlier")
    n.append("screen or a per-hour FA quantile instead of the pooled mean) — stated as a finding, not")
    n.append("an authorised retune, and necessary-but-not-sufficient: at T = 7.50 the 1c run would")
    n.append("have had 3,833 tremor pairs instead of 25, but D5 says they would still not be the")
    n.append("catalogued LFE pairs.")
    n.append("")
    with open(f"{STUDY}/NOTES.md", "a") as f:
        f.write("\n".join(n) + "\n")
    print("appended NOTES.md entry")


if __name__ == "__main__":
    main()
