#!/usr/bin/env python
"""CONFIDENCE report: results/confidence/*.csv -> results/CONFIDENCE_RESULTS.md (tables; the plain-language summary
at the top is read from results/confidence/summary.md if present)."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OD = f"{W.RES}/confidence"
OUTF = f"{W.RES}/CONFIDENCE_RESULTS.md"
SUMMARY = r"""Every matched-filter repeat (step-1 candidate, ≥ 8.0 MAD, nothing cut) gets a confidence P = the probability that it matches the array catalogue, learned directly from the known answer rather than from the data-vs-scrambled-control contrast used in step 1. Inputs: the step-1 features, the step-1 P, PNET (off-the-shelf and fine-tuned on Washington) and the step-2 attribution. One model per label: own family (for firing rates), any family ≤ 3 km, any family (for "is it a real LFE").

**One design change from the spec, and why.** The spec asked for a fit on the build half. But the templates are stacks of build-half array detections, so on the build half a catalogued repeat partly correlates with itself: the raw feature n_z2 separates matched from unmatched with AUC 0.92 on the build half but only 0.67 on the test half, and the own-family match rate is 55–61 % there vs 27–29 % (NW). A model fitted there learns an over-optimistic mapping. It was run as specified (**D1**), and it comes out over-confident on independent data (A2011 calibration slope 0.64). The design to use is **D2**: fit on the TEST half with 5 folds grouped by day, so every test-half row is scored by a model that never saw its day. A2011, a different year, is scored by the model fitted on the whole test half and is the fully independent check. The PNET fine-tune used build-half windows at 36 stations only, so it has no test-half, A2011 or held-out-station leakage (verified from its training metadata).

- **Separation (D2, own family), vs step-1 P / fine-tuned PNET alone / net_mad alone**: A2011 (independent) AUC **0.87** vs 0.64 / 0.78 / 0.68, average precision 0.61 vs 0.38. V2010 test: NW 0.85 (step-1 P 0.71), Olympic 0.90 (0.52), Arm B 0.82 (0.79). V2008 test: NW 0.87 (0.71), SW 0.89 (0.82; only 296 candidates). The day-block bootstrap gain over step-1 P is +0.16 to +0.29 AUC, and every interval excludes 0. Label "any": A2011 0.80 vs 0.57; V2010 0.78 (NW 0.78, Olympic 0.76, Arm B 0.66).
- **Within one family-day** (does P rank the repeats of one family on one day?): 0.73–0.85 for the model in NW / Olympic / SW (Arm B 0.57) vs **0.50–0.63 for step-1 P** and 0.61–0.73 for PNET alone. Step-1 P mainly told busy family-days from quiet ones; the new P also ranks individual repeats.
- **Calibration (D2)**: A2011 slope own / any **0.93 / 0.98** (target 0.8–1.2), but calibration in the large (matched / ΣP) is 0.76 / 0.70. P runs about 25–30 % high on A2011 because its own-family match rate (14 %) is half the training rate, so a base-rate shift needs a per-epoch correction. On the test half (out-of-fold) the slope is 0.95–1.19, except Arm B (own 1.40, any 0.74; only 86 / 498 positives). Step-1 P for comparison: 0.21–0.68.
- **Scrambled control (c)**, own model: P > 0.5 for **0.7 % of control vs 2.2 % of data** (V2010 test, all regions); NW 3.7 % vs 13.6 %; A2011 NW 3.5 % vs 11.1 %; Olympic 0.3 % vs 0.4 %, i.e. no separation, as in step 1. **But the control's total ΣP is 50–65 % of the data's in NW and ~100 % in Olympic.** Much of the low-level P comes from context (a family's busy hours), not from evidence in the individual repeat. The "any family" model is worse: control P > 0.5 in 10–18 % vs data 12–24 %. It has largely learned "busy hour" and must not be used alone as a per-repeat "is it real" score. Adding the scrambled controls as known negatives (**D3**) lowers control P > 0.5 (NW own 3.6 → 2.3 %, any 13.6 → 11.1 %), costs 0.005 AUC and leaves ΣP unchanged.
- **Quiet vs episode days (d)**: on quiet days data and control look the same (NW P > 0.5 0.1 % vs 0.1 %; ΣP ratio 0.99), so there is no false excess on tremor-free days. All the separation comes from episode days (NW 18 % vs 5 %).
- **Expected count (e)**: ΣP per family-day vs the true own-family array count, Spearman V2010 0.52 / V2008 0.46 / A2011 0.47, vs step-1 ΣP 0.49 / 0.45 / 0.48 and the raw candidate count 0.46 / 0.43 / 0.40. Only a small gain, except Olympic (0.75 vs 0.37). ΣP / truth is 0.16–0.35, which reflects recall: most array detections never produce a permanent-station candidate. Subtracting the scrambled-control ΣP makes the ranking worse (ρ 0.34 at A2011), because the context part tracks activity too.
- **Does PNET contribute? (f)** Yes. The **fine-tuned PNET is the most important input**. Permuting it on A2011 drops AUC by 0.16 (own) and 0.20 (any), and removing all PNET features drops A2011 AUC from 0.870 to 0.781 (own) and from 0.804 to 0.699 (any). The off-the-shelf PNET adds nothing (≈ 0.001). **The gain is station-specific**: with PNET read only at the 6 stations the fine-tune never saw, A2011 own AUC is 0.778, no better than without PNET (0.781). Fine-tuned PNET alone at those stations gives AUC 0.58–0.69 vs 0.48–0.65 off-the-shelf. The step-2 attribution adds nothing (AUC drop 0.003); the step-1 P adds a little (0.02).
- **Caveats**: (1) the catalogue is incomplete, so P = P(matches the catalogue), a lower bound on P(real); (2) most of the P mass is family/hour context, so a per-repeat P > 0.5 is informative while a ΣP firing rate still carries a context floor, which the scrambled control measures; (3) the fit is array-era (2008–2011), and after 2012 the station roster, instrument responses, PNET on changed stations, base rate and density features all need re-checking (list at the end).
"""
EV = ["V2010_test", "V2008_test", "A2011"]
REG = ["ALL", "NW", "SW", "Olympic", "Arm B"]


def md(df, fl=3):
    return df.to_markdown(index=False, floatfmt=f".{fl}f")


def order(df):
    df = df.copy()
    df["_e"] = df.eval_set.map({e: i for i, e in enumerate(EV)}); df["_r"] = df.region.map({r: i for i, r in enumerate(REG)})
    return df.sort_values(["_e", "_r"]).drop(columns=["_e", "_r"])


def main():
    L = []
    s = SUMMARY
    L += ["# CONFIDENCE — a per-repeat confidence fitted against the KNOWN ANSWER (array-era truth)", "",
          "Written by `scripts/confidence_report.py` from `results/confidence/` (fit: `scripts/confidence_fit.py`, "
          "held-out-station PNET features: `scripts/confidence_pnet_ho.py`). Scrambled control = the same matched filter "
          "with each site's data shifted to a different time of the same day.", "", "## Summary (plain language)", "", s, ""]
    L += ["## Method", "",
          "- Candidates: every step-1 candidate (T1 all-templates MF, record tier ≥ 8.0 MAD, nothing cut) on V2008 / V2010 / A2011, "
          "`results/step1_prob/candidates.parquet`; 4 site-rolled scrambled controls c1–c4.",
          "- Features: the 21 step-1 features; the step-1 lfdr P (`P1`; on the build half cross-fitted by day, so no build row "
          "gets a P from a model fitted on its own day); PNET off-the-shelf (`_w`) and WA fine-tuned (`_f`): pnet_max, pnet_2nd, "
          "n_sites_pnet03, n_sites_pnet05, pnet_rel_max; step 2: group size, P(member) and n eligible members of the nearest "
          "group-MF candidate within ± 2 s (group frame t + δ_f; exists for data and c1 only; NaN = none / singleton). "
          "Tremor-day labels are NOT features.",
          "- Labels (data candidates on days the family's array run scanned): matched within ± 2 s (nearest) to an array detection "
          "of the own family / any family of the ≤ 3-km group / any family. One model per label.",
          "- Model: HistGradientBoosting (15 leaves, ≥ 300 per leaf, L2 = 2, early stopping) + isotonic calibration "
          "(3 folds grouped by day); logistic regression (C = 0.1, + isotonic) for comparison.",
          "- **D1 (as specified)**: fit on the BUILD half of V2008 + V2010; evaluated on the V2010 / V2008 TEST half and A2011.",
          "- **D2 (template-leak-free)**: the templates are stacks of BUILD-half array detections, so on the build half a "
          "truth-matched candidate is partly correlating with itself (table below). D2 fits on the TEST half with 5 outer folds "
          "grouped by day (every test-half row, data and scrambled control, is scored by a model that never saw its day) and scores "
          "A2011 with a model fitted on the whole test half. A2011 is independent in both designs.",
          "- PNET leakage: the fine-tune used BUILD-half windows of V2008 + V2010 only (`pnet_labels.build_ok`: window inside one "
          "build 2-h block), 36 training stations; the test half, A2011 and the 6 held-out stations (HDW, OSD, B018, B017, NLWA, "
          "GNW) never entered it (`results/pnet/train/meta.parquet`: 36 stations, none held out; no A2011 date). "
          "`PNET_HO` = the same model with PNET read only at those 6 stations.", ""]
    # leakage
    K = pd.read_csv(f"{OD}/leakage_build_vs_test.csv")
    K = K[K.region.isin(["ALL", "NW", "SW", "Olympic", "Arm B"])]
    L += ["### Why D2: build half vs test half (own-family label, raw features, AUC matched vs unmatched)", "",
          md(K[["window", "half", "region", "n", "frac_own", "frac_any", "auc_own_net_mad", "auc_own_n_z2", "auc_own_site_cc_max",
                "auc_own_P1_orig", "auc_own_P1", "auc_own_pnet_max_f", "auc_own_pnet_max_w"]]),
          "", "`P1_orig` = step-1 P as published (in-sample on the build half); `P1` = cross-fitted (used as the feature).", ""]
    # (a) AUC / AP
    A = order(pd.read_csv(f"{OD}/auc_ap.csv"))
    Aa = A[A.subset == "all"]
    for lab, title in [("own", "own family (firing-rate use)"), ("any", "any family ('is it a real LFE')"), ("nbr3", "≤ 3-km group")]:
        g = Aa[Aa.label == lab]
        t = g[["eval_set", "region", "n", "n_pos", "base_rate", "auc_D2_hgb_FULL", "auc_D1_hgb_FULL", "auc_P1", "auc_pnet_max_f",
               "auc_pnet_max_w", "auc_net_mad", "ap_D2_hgb_FULL", "ap_D1_hgb_FULL", "ap_P1", "ap_pnet_max_f", "ap_net_mad"]]
        L += [f"## (a) AUC and average precision on held-out rows — label: {title}", "", md(t), ""]
    B = pd.read_csv(f"{OD}/bootstrap_dauc.csv")
    L += ["Day-block bootstrap (200 resamples of days) of AUC(model) − AUC(step-1 P), all regions:", "", md(B), ""]
    # ablations
    cols = [c for c in Aa.columns if c.startswith("auc_D") ]
    ab = Aa[Aa.region.isin(["ALL", "NW"]) & Aa.label.isin(["own", "any"])][["eval_set", "region", "label"] + cols]
    ab.columns = [c.replace("auc_", "") for c in ab.columns]
    L += ["### Ablations (AUC; D = design, hgb / logit, feature set)", "",
          "FULL = all features; NO_PNET = without any PNET feature; NO_STEP2 = without step-2 attribution; NO_P1 = without the "
          "step-1 P; PNET_HO = PNET only from the 6 never-fine-tuned stations; S1_ONLY = the 21 step-1 features only.", "", md(ab), ""]
    ho = A[(A.subset == "near_ho") & A.label.isin(["own", "any"])][
        ["eval_set", "region", "label", "n", "n_pos", "auc_pnet_max_ho_f", "auc_pnet_max_ho_w", "auc_pnet_max_f", "auc_pnet_max_w",
         "auc_D2_hgb_PNET_HO", "auc_D2_hgb_NO_PNET", "auc_D2_hgb_FULL", "auc_P1"]]
    L += ["### Held-out PNET stations: families with HDW / OSD / B018 / B017 / NLWA / GNW among their 2 nearest stations", "",
          "`pnet_max_ho_*` = PNET read only at the 6 held-out stations (never used in the fine-tune).", "", md(ho), ""]
    # (b) calibration
    C = order(pd.read_csv(f"{OD}/calibration_slopes.csv"))
    for lab in ["own", "any", "nbr3"]:
        g = C[(C.label == lab) & C.model.isin(["D2_hgb_FULL", "D1_hgb_FULL", "D2_logit_FULL", "D1_logit_FULL", "step1_lfdr_P1",
                                                "D2_hgb_NO_PNET"])]
        t = g.pivot_table(index=["eval_set", "region"], columns="model", values=["slope", "obs_over_exp"], sort=False)
        t.columns = [f"{a} {b}" for a, b in t.columns]
        t = t.reset_index()
        L += [f"## (b) Calibration on held-out rows — label {lab} (slope target 0.8–1.2; obs/exp = matched / ΣP)", "", md(t, 2), ""]
    R = pd.read_csv(f"{OD}/reliability.csv")
    for es, rg, lab in [("V2010_test", "ALL", "own"), ("V2008_test", "ALL", "own"), ("A2011", "ALL", "own"), ("V2010_test", "ALL", "any"),
                        ("A2011", "ALL", "any")]:
        g = R[(R.eval_set == es) & (R.region == rg) & (R.label == lab)]
        t = g.pivot_table(index="bin", columns="model", values=["n", "mean_P", "frac"])
        t.columns = [f"{a} {b}" for a, b in t.columns]
        t = t.reset_index(); t["bin"] = t.bin.map(lambda b: f"{b/10:.1f}-{(b+1)/10:.1f}")
        L += [f"Reliability, {es}, {rg}, label {lab}:", "", md(t, 3), ""]
    # (c)(d) controls + quiet
    Q = order(pd.read_csv(f"{OD}/control_quiet.csv"))
    for key in ["D2_hgb_FULL_own", "D2_hgb_FULL_any", "D2_hgb_NO_STEP2_own", "D2_hgb_NO_STEP2_any", "D1_hgb_FULL_own", "D1_hgb_FULL_any", "P1"]:
        g = Q[Q.model == key][["eval_set", "region", "days", "controls", "n_data", "n_ctl_per_real", "data_fracP_gt05", "ctl_fracP_gt05",
                               "data_meanP", "ctl_meanP", "data_sumP", "ctl_sumP_per_real", "ctl_over_data_sumP"]]
        L += [f"## (c)/(d) Scrambled control and quiet vs episode days — model {key}", "",
              "FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. "
              "`ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.", "", md(g), ""]
    T = order(pd.read_csv(f"{OD}/tremor_label_split.csv"))
    L += ["### By PNSN regional tremor-day label (label only, never a feature)", "", md(T), ""]
    # (e)
    E = order(pd.read_csv(f"{OD}/expected_count.csv"))
    L += ["## (e) Expected count — ΣP per family-day vs the true array count (episode family-days the family's array run scanned; "
          "scored half)", "", "The truth counts every array detection, including the many the permanent-station MF never "
          "produced a candidate for, so ratio < 1 is recall, not mis-calibration; the rank agreement is the test.", "", md(E), ""]
    # (f)
    F = pd.read_csv(f"{OD}/feature_importance.csv")
    for key in F.model.unique():
        g = F[F.model == key].sort_values("auc_drop", ascending=False)
        gb = g[g.is_block]; gf = g[~g.is_block].head(15)
        L += [f"## (f) Permutation importance — {key} (rows: {g.eval_rows.iloc[0]}, n = {g.n.iloc[0]})", "",
              "Blocks (features permuted jointly):", "", md(gb[["feature", "auc_drop", "logloss_gain"]], 4), "",
              "Top 15 single features:", "", md(gf[["feature", "auc_drop", "logloss_gain"]], 4), ""]
    # D3 add-on
    if os.path.exists(f"{OD}/d3_auc_calib.csv"):
        A3 = order(pd.read_csv(f"{OD}/d3_auc_calib.csv"))
        A3 = A3[A3.label.isin(["own", "any"])]
        L += ["## D3 — scrambled-control candidates added to the fit as known negatives (`scripts/confidence_d3.py`)", "",
              "Same as D2 with the NO_STEP2 features, but c1–c3 candidates of the same family-days enter the fit with label 0 "
              "(weighted like the data rows); isotonic calibration on data rows only; c4 never enters any fit. "
              "`withinFD` = AUC computed inside each family-day and averaged (pairs-weighted): does P rank the repeats of one "
              "family on one day, beyond knowing which family-days are busy?", "",
              md(A3[["eval_set", "region", "label", "n", "n_pos", "auc_D3_hgb_NO_STEP2", "auc_D2_hgb_NO_STEP2", "ap_D3_hgb_NO_STEP2",
                     "ap_D2_hgb_NO_STEP2", "slope_D3_hgb_NO_STEP2", "slope_D2_hgb_NO_STEP2", "obs_exp_D3_hgb_NO_STEP2",
                     "obs_exp_D2_hgb_NO_STEP2", "auc_withinFD_D3_hgb_NO_STEP2", "auc_withinFD_D2_hgb_NO_STEP2",
                     "auc_withinFD_P1", "auc_withinFD_pnet_max_f"]]), ""]
        Q3 = order(pd.read_csv(f"{OD}/d3_control_quiet.csv"))
        Q3 = Q3[Q3.region.isin(["ALL", "NW", "Olympic"]) & Q3.days.isin(["all", "quiet"])]
        L += ["D3 vs D2 on the scrambled control (c1–c3: seen in training on OTHER days only; c4: never seen):", "",
              md(Q3[["eval_set", "region", "days", "model", "controls", "data_fracP_gt05", "ctl_fracP_gt05", "data_sumP",
                     "ctl_sumP_per_real", "ctl_over_data_sumP"]]), ""]
        E3 = order(pd.read_csv(f"{OD}/d3_expected_count.csv"))
        L += ["D3 expected count:", "", md(E3), ""]
    if os.path.exists(f"{OD}/excess_count.csv"):
        X = order(pd.read_csv(f"{OD}/excess_count.csv"))
        L += ["### Expected count with the scrambled control subtracted (ΣP data − mean ΣP of c1–c4, per family-day)", "",
              "`rho_*` = Spearman vs the true own-family array count.", "", md(X), ""]
    L += ["## Caveats", "",
          "- The array catalogue is incomplete: P is the probability that a candidate MATCHES THE CATALOGUE, a lower bound on P(real). "
          "A real LFE the array run missed counts as a negative.",
          "- 'any family' matches include chance coincidences on busy tremor days (the scrambled control matches 'any' at ~5–8 %), so "
          "P(any) mixes real-LFE evidence with 'busy hour' evidence; P(own) is the cleaner label for firing rates.",
          "- The model is fitted in the array era (2008–2011). Before using it after 2012, re-check: (1) station roster / epochs "
          "(the site-count and n_z features, dom_top and wsum depend on which stations are live; new, dead or replaced stations shift "
          "them); (2) instrument / response changes at PB boreholes and UW stations (PNET inputs, HF ratio); (3) the fine-tuned PNET "
          "on stations and years it never saw (use PNET_HO-type checks with a post-2012 known answer, e.g. PNSN tremor-day "
          "coincidence or a later temporary array); (4) the region mix and base rate (the isotonic map assumes the array-era "
          "matched fraction; a quieter or busier epoch shifts calibration in the large); (5) the density features "
          "(dens_own / dens_grp) whose meaning depends on the candidate rate per day; (6) templates — post-2012 data never "
          "contributed to them, like the test half here (D2 is the right analogue, D1 is not).", ""]
    open(OUTF, "w").write("\n".join(L))
    print(OUTF)


if __name__ == "__main__":
    main()
