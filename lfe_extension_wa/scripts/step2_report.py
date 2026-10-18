#!/usr/bin/env python
"""STEP 2 report -> results/STEP2_GROUPS_RESULTS.md + figures/step2_groups.png (numbers; interpretation after review)."""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

G = f"{W.RES}/step2_groups"
SC = f"{G}/score"
OUT = f"{W.RES}/STEP2_GROUPS_RESULTS.md"
FIG = f"{W.STUDY}/figures/step2_groups.png"
C_GRP, C_T1, C_GRAY = "#2a78d6", "#eb6834", "#8a8a85"


def f3(x):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.3f}"


def f2(x):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2f}"


def md(df):
    return df.to_markdown(index=False, floatfmt=".3f")


def main():
    L = []
    GR = pd.read_csv(f"{G}/groups.csv")
    M = pd.read_csv(f"{G}/members.csv")
    Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_step2grp.csv")
    sz = GR.n_members.value_counts().sort_index()
    S = pd.read_csv(f"{SC}/summary.csv") if os.path.exists(f"{SC}/summary.csv") else None
    SEP = pd.read_csv(f"{SC}/sep_curve.csv") if os.path.exists(f"{SC}/sep_curve.csv") else None
    ACC = pd.read_csv(f"{G}/attrib_accuracy.csv") if os.path.exists(f"{G}/attrib_accuracy.csv") else None
    CAL = pd.read_csv(f"{G}/attrib_calibration.csv") if os.path.exists(f"{G}/attrib_calibration.csv") else None
    EXP = pd.read_csv(f"{G}/attrib_expected_count.csv") if os.path.exists(f"{G}/attrib_expected_count.csv") else None
    PK = pd.read_csv(f"{G}/attrib_peakedness.csv") if os.path.exists(f"{G}/attrib_peakedness.csv") else None
    MT = pd.read_csv(f"{G}/attrib_multi_true.csv") if os.path.exists(f"{G}/attrib_multi_true.csv") else None
    BG = pd.read_csv(f"{G}/attrib_by_group.csv") if os.path.exists(f"{G}/attrib_by_group.csv") else None
    MJ = json.load(open(f"{G}/attrib_model.json")) if os.path.exists(f"{G}/attrib_model.json") else None
    SCK = pd.read_csv(f"{SC}/singleton_check.csv") if os.path.exists(f"{SC}/singleton_check.csv") else None
    TCK = json.load(open(f"{G}/singleton_template_check.json")) if os.path.exists(f"{G}/singleton_template_check.json") else None
    S1 = f"{G}/step1"
    nonref = M[~M.is_ref]
    alok = nonref.codet_prior.astype(bool) | (nonref.align_ratio >= 1.2)
    q = Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna()].merge(GR[["fid", "n_members"]], on="fid")
    qm = q[q.n_members > 1]
    L.append("# STEP 2 — detect ≤ 3-km GROUPS of families, then attribute which member fired\n")
    L.append(f"Written {time.strftime('%Y-%m-%d %H:%MZ', time.gmtime())} by `scripts/step2_report.py` (numbers; interpretation "
             "after review). Queue `scripts/step2_groups_queue.sh`; outputs `results/step2_groups/`; figure `figures/step2_groups.png`.\n")
    # ------------------------------------------------ summary
    L.append("## Summary (plain language)\n")
    L.append("Families closer than 3 km cannot be told apart reliably by location, and T1 showed that their permanent-station "
             "detections overlap. Step 2 treats each such cluster as one GROUP: one template per station built from all members' "
             "array-era events, one matched filter per group, scored against the union of the members' array detections. Then, "
             "for every group detection, it asks which member fired, using only the stations of the 2–3 nearest sites, as a "
             "probability. The scrambled control (each site shifted to a different time of the same day) is the same one T1 used.\n")
    L.append(f"- **Groups** (complete linkage, every pair ≤ 3 km): {len(GR)} groups from 329 families; sizes "
             + ", ".join(f"{k}: {v}" for k, v in sz.items()) + f". {int((GR.n_members > 1).sum())} multi-member groups hold "
             f"{int(GR.n_members[GR.n_members > 1].sum())} families. The 23-member group is 21 UNLOCATED Olympic Arm A families whose "
             "approximate position is the B013/B943 borehole pair, plus 2 NW families there — a group by position rule, not by evidence.")
    L.append(f"- **Group templates**: {len(q)} ≥ 30-event templates for {q.fid.nunique()} groups; in multi-member groups the median "
             f"template stacks {qm.n_events.median():.0f} events vs {qm.n_events_best_member.median():.0f} for the best single member "
             f"(median gain ×{(qm.n_events / qm.n_events_best_member.clip(lower=1)).median():.2f}); clear templates {int(Q.clear.sum())}. "
             f"Member alignment reliable for {int(alok.sum())} of {len(nonref)} non-reference members (co-detection prior or stack-cc peak "
             "≥ 1.2× the next peak); the others' origin shift is uncertain (their stacks are noise-like on permanent stations).")
    if S is not None:
        for rg in ["ALL", "Royer_NW_WA", "Royer_SW_WA", "Olympic_2010", "Olympic_2010_armB"]:
            a = S[(S.run == "grp") & (S.window == "ALL") & (S.region == rg) & (S.size_class == "ALL")]
            b = S[(S.run == "t1union") & (S.window == "ALL") & (S.region == rg) & (S.size_class == "ALL")]
            if not len(a) or not len(b):
                continue
            a, b = a.iloc[0], b.iloc[0]
            L.append(f"- **Detection {rg}** (core tier, raw): observable groups {int(b.n_groups_obs)} (T1 summed) → {int(a.n_groups_obs)} "
                     f"(group MF); recall over all scanned truth {f3(b.recall_all_raw_core)} → {f3(a.recall_all_raw_core)}; precision "
                     f"{f3(b.prec_raw_core)} → {f3(a.prec_raw_core)} (scrambled control's chance {f3(b.prec_ctl_core)} / {f3(a.prec_ctl_core)}); "
                     f"per group-day data/control {f2(b.data_per_gday_core)}/{f2(b.ctl_per_gday_core)} → {f2(a.data_per_gday_core)}/"
                     f"{f2(a.ctl_per_gday_core)}; group-windows ≥ 3× above the control {int(b.sep3_core_data)} (control {int(b.sep3_core_ctl)}) → "
                     f"{int(a.sep3_core_data)} (control {int(a.sep3_core_ctl)}).")
    stop_line = None
    if ACC is not None:
        base = ACC[(ACC.region == "ALL") & (ACC.size_class == "ALL") & (ACC.tier == "all")]
        for es in ["V2010_test", "V2008_test", "A2011"]:
            a = base[(base.eval_set == es) & (base.model == "cc_k3")]
            p = base[(base.eval_set == es) & (base.model == "cc_k3_prior")]
            if not len(a):
                continue
            a = a.iloc[0]; p = p.iloc[0] if len(p) else None
            L.append(f"- **Attribution {es}** ({int(a.n_cand)} group detections matched to one member's array detection): the "
                     f"top member is the true one {f3(a.top1)} [{f2(a.top1_lo)}–{f2(a.top1_hi)}] vs chance {f3(a.chance)} and vs 'always the "
                     f"most active member' {f3(a.most_active_top1)}; mean P(true) {f3(a.mean_P_true)} (uniform {f3(a.mean_P_unif)}); "
                     f"log-loss skill vs uniform {f3(a.logloss_skill_vs_unif)}" + (f"; with the activity prior added top-1 {f3(p.top1)}, "
                                                                                  f"skill {f3(p.logloss_skill_vs_unif)}" if p is not None else "") + ".")
        te = base[base.eval_set.isin(["V2010_test", "V2008_test", "A2011"]) & (base.model == "cc_k3")]
        if len(te):
            n = te.n_cand.sum()
            top = (te.top1 * te.n_cand).sum() / n; ch = (te.chance * te.n_cand).sum() / n
            ma = (te.most_active_top1 * te.n_cand).sum() / n
            sk = (te.logloss_skill_vs_unif * te.n_cand).sum() / n
            if sk < 0.05 and top < 1.2 * ch:
                stop_line = (f"**Merlin's stop signal FIRES**: at the nearest sites attribution is ≈ uniform (top-1 {top:.3f} vs chance {ch:.3f}, "
                             f"skill {sk:.3f}) → member identity is not carried past the array era on permanent stations.")
            elif top <= ma + 0.02:
                stop_line = (f"Stop-signal check: attribution beats chance (top-1 {top:.3f} vs {ch:.3f}, skill {sk:.3f}) but NOT the trivial "
                             f"'most active member' rule ({ma:.3f}) → the cc adds little identity beyond member activity.")
            else:
                stop_line = (f"Stop-signal check (pooled test sets, n = {int(n)}): attribution is NOT uniform — top-1 {top:.3f} vs chance "
                             f"{ch:.3f} and vs most-active {ma:.3f}, log-loss skill {sk:.3f} → member identity is carried by the nearest "
                             "permanent sites for these (array-era) events.")
            L.append(f"- {stop_line}")
    if EXP is not None and len(EXP):
        e = EXP[(EXP.window == "ALL") & (EXP.region == "ALL")]
        for r in e.itertuples():
            L.append(f"- **Family-level expected count** ({r.families}, {r.n_famdays} episode family-days): Spearman with the true array "
                     f"count Σ P(real)·P(member) {f2(r.spearman_E)} (uniform split {f2(r.spearman_E_unif)}; step-1 family ΣP {f2(r.spearman_P_T1)}); "
                     f"ΣE {r.sum_E:.0f} vs truth {r.truth} (step-1 ΣP {r.sum_P_T1:.0f}).")
    L.append("\nCaveats: the array catalogue is the truth but incomplete; merged co-detections (the same LFE in two members' catalogues) are "
             "counted once; attribution temperature is fitted on the BUILD half, whose events are inside the member templates "
             "(in-sample for the templates) — the TEST half and A2011 are the honest numbers.\n")
    # ------------------------------------------------ method
    L.append("## Method\n")
    L.append("- **Group rule**: complete-linkage clustering of the epicentral distances (lat_use/lon_use: location, else the catalogue's "
             "approximate position), cut at 3.0 km → every pair in a group ≤ 3 km, no chaining. Family without a position (OLY:526) = own group. "
             "Single-linkage (for reference) would chain into groups of 20, 24 and 33.")
    L.append("- **Alignment** (deviation, data-driven): each family's detection time is on its own origin convention (offsets between "
             "co-located families up to ±9 s). A member's shift δ to the reference member (most build events) = the lag maximising the summed "
             "3-C normalised cross-correlation of their build-half long stacks over all shared station-nodes (±10 s), searched within ±1 s of "
             "the co-detection mode when ≥ 20 of the member's detections have a reference detection within 8 s and the mode bin holds ≥ max(5, "
             "3× uniform); per-station residual τ within δ ± 1 s. The own S picks ('template offsets') were not used: on the ≥ 30-event stacks "
             "they are noise-dominated at most stations. Where δ is strong it equals the co-detection offset to ≤ 0.15 s.")
    L.append("- **Truth**: union of member array detections in the group frame (t + δ), merged ±2 s; a group-day is scanned when every "
             "member whose windows include it was scanned. Templates: all members' build-half events (V2008 + V2010 half 0, co-detections "
             "once) cut at t + τ, wa_templates rule, weights on the T1 build scale. Every ≥ 30-event template kept (clear = a column).")
    L.append("- **Group MF**: `wa_mf.py --all-templates --drop-dead 0.05 --n-control 2` (seed as T1 → identical site shifts) + a second run "
             "with seed 20261001 (2 more scrambled controls, as step 1). Scoring `scripts/step2_score.py` (see its docstring).")
    L.append("- **Attribution** `scripts/step2_attrib.py`: member templates cut from each member's own build stack at the group window "
             "(≥ 30 events at that node); nearest 3 sites with ≥ 2 member templates live that day; eligible members = template at every "
             "chosen station; cc per station at the aligned position and ±0.5 s (sub-sample lag by parabola); P = softmax(β Σ cc0), β by "
             "maximum likelihood on build-half candidates matched to one member.")
    L.append("- **P(real)**: step-1 lfdr (`step1_features.py` + `step1_fit.py`, same features, 4 scrambled controls, fit on build half) "
             "refitted on the group candidates → `results/step2_groups/step1/`.\n")
    if TCK or SCK is not None:
        L.append("Checks: " + (f"singleton group templates vs T1 build templates: {TCK['n_compared']} compared, max |diff| "
                               f"{TCK['max_abs_diff']}, offset mismatches {TCK['offset_mismatch']}. " if TCK else "") +
                 (f"Singleton group detections vs T1: {int(SCK.n_grp.iloc[0])} vs {int(SCK.n_t1.iloc[0])}, identical = {bool(SCK.identical.iloc[0])}." if SCK is not None else "") + "\n")
    L.append("### Group sizes\n")
    t = GR.groupby("n_members").agg(groups=("fid", "size"), families=("n_members", "sum"),
                                    max_pair_km=("max_pair_km", "max"), with_unlocated=("n_unlocated", lambda v: int((v > 0).sum()))).reset_index()
    L.append(md(t) + "\n")
    t = GR.groupby("region").agg(groups=("fid", "size"), multi=("n_members", lambda v: int((v > 1).sum())),
                                 families=("n_members", "sum")).reset_index()
    L.append(md(t) + "\n")
    t = nonref.assign(reliable=alok).groupby("region").agg(members=("fid", "size"), codet_prior=("codet_prior", "sum"),
                                                          reliable=("reliable", "sum"), unaligned=("aligned", lambda v: int((~v.astype(bool)).sum()))).reset_index()
    L.append("Member alignment (non-reference members):\n\n" + md(t) + "\n")
    # ------------------------------------------------ detection
    if S is not None:
        L.append("## Group detection vs T1 summed over members\n")
        L.append("`grp` = group matched filter; `t1union` = T1 family detections of the members in the group frame, merged ±2 s. Recall = union "
                 "truth matched one-to-one ±2 s on observable scanned days (`recall all` = over ALL scanned truth); `recall al` = truth of "
                 "well-aligned members only; precision own = matched / ours, any = within ±2 s of any group's truth; control = the scrambled "
                 "control's chance precision. sep = (data+1)/(control+1) on label days; control column = the same on control 2 vs control 1.\n")
        for tier, lab in [("core", "raw"), ("core", "cat"), ("core_ep", "raw"), ("ext", "raw")]:
            k = f"{lab}_{tier}"
            x = S[(S.window == "ALL") & (S.size_class == "ALL")].copy()
            x = x.sort_values(["region", "run"])
            cols = ["region", "run", "n_groups_obs", "n_families_obs", "n_truth_test", "n_truth_scored_all", f"recall_{k}", f"recall_all_{k}",
                    f"recall_alignok_{k}", f"prec_{k}", f"prec_any_{k}", f"prec_ctl_{tier}", f"data_per_gday_{tier}", f"ctl_per_gday_{tier}",
                    f"sep2_{tier}_data", f"sep2_{tier}_ctl", f"sep3_{tier}_data", f"sep3_{tier}_ctl", f"sep5_{tier}_data", f"sep5_{tier}_ctl"]
            cols = [c for c in cols if c in x]
            L.append(f"**Tier `{tier}`, {lab}** ({'tier only' if lab == 'raw' else 'tier & HF veto & cross-group rule'})\n")
            L.append(md(x[cols].rename(columns=lambda c: c.replace(f"_{k}", "").replace(f"_{tier}", ""))) + "\n")
        x = S[(S.window == "ALL") & (S.region == "ALL") & (S.size_class != "ALL")].sort_values(["size_class", "run"])
        L.append("By group size (core, raw):\n\n" + md(x[["size_class", "run", "n_groups_obs", "n_truth_test", "recall_raw_core", "recall_all_raw_core",
                                                           "prec_raw_core", "prec_ctl_core", "data_per_gday_core", "ctl_per_gday_core",
                                                           "sep3_core_data", "sep3_core_ctl"]]) + "\n")
        x = S[(S.window != "ALL") & (S.region != "ALL")].sort_values(["window", "region", "run"])
        L.append("By window and region (core, raw):\n\n" + md(x[["window", "region", "run", "n_groups_obs", "n_truth_test", "recall_raw_core",
                                                                  "recall_all_raw_core", "prec_raw_core", "prec_ctl_core", "data_per_gday_core",
                                                                  "ctl_per_gday_core", "sep3_core_data", "sep3_core_ctl"]]) + "\n")
        if os.path.exists(f"{SC}/tiers_by_epoch.csv"):
            L.append("Group-run per-epoch tiers:\n\n" + md(pd.read_csv(f"{SC}/tiers_by_epoch.csv")) + "\n")
    if SEP is not None:
        x = SEP[SEP.tier == "core"].pivot_table(index="sep", columns="run", values=["n_data", "n_ctl"]).reset_index()
        x.columns = ["sep ≥"] + [f"{a} {b}" for a, b in x.columns[1:]]
        L.append("Separation curve (core, raw; number of group-windows with label days whose data/control ratio ≥ x, and the same on the "
                 "control):\n\n" + md(x) + "\n")
    # ------------------------------------------------ attribution
    if ACC is not None:
        L.append("## Attribution — which member fired\n")
        if MJ:
            L.append(f"Fitted on {MJ['n_fit_candidates']} build-half candidates: " + "; ".join(f"{k} β = {', '.join(f'{v:.2f}' for v in p)}"
                                                                                           for k, p in MJ["params"].items()) + ".\n")
        x = ACC[(ACC.region == "ALL") & (ACC.size_class == "ALL") & (ACC.tier == "all") & (ACC.model != "unif")]
        L.append("All groups (candidates matched to exactly one eligible member's array detection). `top1` = the highest-P member is the "
                 "true one; `chance` = mean 1/n eligible; `most_active` = always pick the member with the most build events; skill = "
                 "1 − log-loss / uniform log-loss.\n")
        L.append(md(x[["eval_set", "model", "n_cand", "top1", "top1_lo", "top1_hi", "chance", "most_active_top1", "mean_P_true",
                       "mean_P_unif", "logloss_skill_vs_unif", "mean_n_elig"]]) + "\n")
        x = ACC[(ACC.model == "cc_k3") & ((ACC.region != "ALL") | (ACC.size_class != "ALL") | (ACC.tier != "all"))]
        L.append("Primary model (cc, 3 nearest sites) by region / group size / tier:\n\n" +
                 md(x[["eval_set", "region", "size_class", "tier", "n_cand", "top1", "chance", "most_active_top1", "mean_P_true",
                       "logloss_skill_vs_unif"]]) + "\n")
        if MT is not None and len(MT):
            L.append("Co-detected events (the same LFE in ≥ 2 members' catalogues): top member in the true set vs chance |set|/n:\n\n" + md(MT) + "\n")
        if PK is not None:
            L.append("How peaked P(member) is on data vs the scrambled control (all record-tier candidates; the control has no true member):\n\n" + md(PK) + "\n")
        if CAL is not None:
            x = CAL[CAL.eval_set != "fit"].pivot_table(index="pbin", columns="eval_set", values=["n", "mean_P", "frac_true"]).reset_index()
            x.columns = ["P bin"] + [f"{a} {b}" for a, b in x.columns[1:]]
            L.append("Reliability of P(member) (every candidate × member pair):\n\n" + md(x) + "\n")
        if BG is not None:
            L.append("Per group (test sets pooled, ≥ 20 candidates):\n\n" + md(BG[BG.n_cand >= 20].head(40)) + "\n")
    # ------------------------------------------------ P(real) refit + combination
    if os.path.exists(f"{S1}/calibration_slopes.csv"):
        L.append("## P(real) on group candidates (step-1 lfdr refitted)\n")
        cs = pd.read_csv(f"{S1}/calibration_slopes.csv")
        c1 = pd.read_csv(f"{W.RES}/step1_prob/calibration_slopes.csv")
        rows = []
        for es in ["V2010_test", "V2008_test", "A2011"]:
            for lab, cc in [("group", cs), ("family (step 1)", c1)]:
                r = cc[(cc.eval_set == es) & (cc.region == "ALL") & (cc.series == "data")]
                if len(r):
                    r = r.iloc[0]
                    rows.append(dict(eval_set=es, level=lab, n=int(r.n), sum_P=r.sum_P, slope_own=r.slope_own, matched_over_sumP=r.obs_over_exp_own))
        L.append(md(pd.DataFrame(rows)) + "\n")
        cq = pd.read_csv(f"{S1}/control_and_quiet.csv")
        cq1 = pd.read_csv(f"{W.RES}/step1_prob/control_and_quiet.csv")
        rows = []
        for es in ["V2010_test", "V2008_test", "A2011"]:
            for lab, cc in [("group", cq), ("family (step 1)", cq1)]:
                r = cc[(cc.eval_set == es) & (cc.region == "ALL") & (cc.days == "all")]
                if len(r):
                    r = r.iloc[0]
                    rows.append(dict(eval_set=es, level=lab, data_n=int(r.n_data), data_P_gt05=r.data_fracP_gt05, data_sumP=r.data_sumP,
                                     ctl_P_gt05=r.ctl_fracP_gt05, ctl_P_gt05_never_seen=r.ctl_loco_fracP_gt05, ctl_sumP_per_real=r.ctl_sumP_per_real))
        L.append("Scrambled control vs data:\n\n" + md(pd.DataFrame(rows)) + "\n")
        if os.path.exists(f"{S1}/expected_count.csv"):
            e = pd.read_csv(f"{S1}/expected_count.csv")
            L.append("Group-level expected count (ΣP per group-day vs the union truth):\n\n" +
                     md(e[["eval_set", "region", "n_famdays", "sumP", "sumP_ctl", "truth_own", "ratio_P_over_own", "spearman_own",
                           "spearman_ncand_own"]].rename(columns={"n_famdays": "group_days"})) + "\n")
    if EXP is not None:
        L.append("### Family-level expected count Σ P(real)·P(member)\n")
        L.append("Per member family-day (episode days the family's array run scanned; scored half): `E` = attributed, `E_unif` = P(real) "
                 "split equally over the group, `P_T1` = step 1's family-level ΣP; Spearman against the family's own array count.\n")
        L.append(md(EXP) + "\n")
    L.append("Files: `results/step2_groups/` — groups.csv, members.csv (δ, alignment quality), align_station.parquet (τ), setup/ (group "
             "families / truth / scan days), score/ (by_group, summary, sep_curve, tiers, HF cuts, singleton check), attrib/ (per-day cc), "
             "attrib_pairs.parquet (every candidate × member with P), attrib_*.csv, attrib_expected_*.parquet, step1/ (group P(real)). "
             "Templates `results/templates/*_step2grp.*`; MF `results/mf_step2grp`, `results/mf_step2grp_c34`.\n")
    # ------------------------------------------------ figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
        if SEP is not None:
            for run, col, lab in [("grp", C_GRP, "group MF"), ("t1union", C_T1, "T1 summed over members")]:
                x = SEP[(SEP.run == run) & (SEP.tier == "core")]
                ax[0].plot(x.sep, x.n_data, "-o", color=col, lw=2, ms=5, label=f"{lab}: data / control 1")
                ax[0].plot(x.sep, x.n_ctl, "--o", color=col, lw=2, ms=5, mfc="white", label=f"{lab}: control 2 / control 1")
            ax[0].set_xscale("log"); ax[0].set_xlabel("separation ratio ≥ x"); ax[0].set_ylabel("group-windows")
            ax[0].set_title("Separation from the scrambled control (core, raw)"); ax[0].legend(fontsize=8, frameon=False)
        if ACC is not None:
            x = ACC[(ACC.region == "ALL") & (ACC.size_class == "ALL") & (ACC.tier == "all") & (ACC.model == "cc_k3")]
            x = x.set_index("eval_set").reindex(["fit", "V2010_test", "V2008_test", "A2011"]).dropna(subset=["top1"])
            xi = np.arange(len(x))
            ax[1].bar(xi - 0.27, x.top1, 0.25, color=C_GRP, label="cc, 3 nearest sites")
            ax[1].bar(xi, x.most_active_top1, 0.25, color=C_T1, label="most active member")
            ax[1].bar(xi + 0.27, x.chance, 0.25, color=C_GRAY, label="chance 1/n")
            ax[1].set_xticks(xi); ax[1].set_xticklabels([s.replace("_", "\n") for s in x.index])
            ax[1].set_ylim(0, 1); ax[1].set_ylabel("top-1 = true member"); ax[1].set_title("Attribution (matched group detections)")
            ax[1].legend(fontsize=8, frameon=False)
        if CAL is not None:
            for es, col in [("V2010_test", C_GRP), ("V2008_test", C_T1), ("A2011", "#1baf7a")]:
                x = CAL[(CAL.eval_set == es) & (CAL.n >= 10)]
                ax[2].plot(x.mean_P, x.frac_true, "-o", color=col, lw=2, ms=5, label=es)
            ax[2].plot([0, 1], [0, 1], color=C_GRAY, lw=1)
            ax[2].set_xlabel("P(member)"); ax[2].set_ylabel("fraction true member"); ax[2].set_title("Reliability of P(member)")
            ax[2].legend(fontsize=8, frameon=False)
        for a in ax:
            a.spines[["top", "right"]].set_visible(False); a.grid(alpha=0.25, lw=0.5)
        fig.tight_layout(); os.makedirs(os.path.dirname(FIG), exist_ok=True); fig.savefig(FIG, dpi=130)
        L.append("![step2](../figures/step2_groups.png)\n")
    except Exception as e:  # figure is optional
        L.append(f"(figure failed: {e})\n")
    open(OUT, "w").write("\n".join(L))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
