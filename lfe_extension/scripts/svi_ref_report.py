#!/usr/bin/env python
"""SVI_REF report -> $SVI_REF_REPORT (default results/SVI_REF_RESULTS.md): plain summary first, then one section per stage,
every number with the file it comes from.  Robust to stages that were skipped (it says so)."""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R

O = R.OUT
REPORT = os.environ.get("SVI_REF_REPORT", f"{R.RES}/SVI_REF_RESULTS.md")


def rel(p):
    return os.path.relpath(p, R.EXT) if p.startswith(R.EXT) else p


def J(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def f3(x, n=3):
    try:
        return f"{float(x):.{n}f}"
    except Exception:
        return "n/a"


def main():
    L, S = [], []        # long sections, short summary bullets
    r0w = J(f"{O}/r0/watcher.json"); r0v = J(f"{O}/r0/v11_verify.json")
    r1 = J(f"{O}/r1/r1_summary.json"); pil = J(f"{O}/cut_pilot/PILOT.json")
    pc = J(f"{O}/r2/r2_positive_control.json"); r3 = J(f"{O}/r3/r3_decision.json")
    bars = J(f"{O}/conf/sealed_bars.json"); gate = J(f"{O}/conf/pnet_gate.json"); frz = J(f"{O}/conf/MODELS_FROZEN.json")
    opened = J(f"{O}/conf/SEALED_OPENED.json"); val = J(f"{O}/val/val_summary.json"); runs = J(f"{O}/conf/runs.json")
    # ---------------------------------------------------------------- summary
    S.append("The southern-Vancouver-Island refinement queue (results/SVI_REFINE_PLAN.md as amended by Merlin) ran on the data "
             "already on disk. Nothing was downloaded, nothing was deleted, and no v1.1 file was overwritten"
             + (f" (re-hash of {r0v['n_files']} v1.1 files: {r0v['n_changed_or_missing']} changed; `{rel(O)}/r0/v11_verify.json`)." if r0v else "."))
    if r1:
        S.append(f"**Completeness (R1).** {r1['n_expected']:,} expected station-days 2014-2026. v1.1 had silently lost "
                 f"{r1.get('cache_rebuilt_ok', 0) + r1.get('cache_rebuild_errors', 0)} station-days at the cache stage; they were rebuilt into "
                 f"`cache20hz_r1/` ({r1.get('cache_rebuilt_ok', 0)} built, {r1.get('cache_rebuild_errors', 0)} errors). After the rebuild "
                 f"{r1['months_fail_cache']} station-months lose > 3 % at the cache stage, {r1['months_fail_raw']} at the raw stage. "
                 f"{r1['station_months_masked_gt3']} station-months have > 3 % of windows masked by the energy floor (dead or gappy channels; handed to R2). "
                 f"(`{rel(O)}/r1/r1_summary.json`)")
    if pc:
        pa = pc.get("a2_reproduction", {})
        a2 = pa.get("a2_inputs_LRIV", {}); e2 = pa.get("r2_events_LRIV", {}); b2 = pa.get("r2_events_B926", {})
        S.append(f"**Station health (R2), LRIV positive control.** Track A-2's LRIV number is reproduced exactly from A-2's own inputs "
                 f"(cc median {f3(a2.get('cc_median'))} against the stored {f3(a2.get('stored_a2_median'))}). The same A-2 procedure on the R2 "
                 f"event set (v1.1 catalogue repeats inside tremor windows, 2021 doy 276-302) gives LRIV {f3(e2.get('cc_median'))} "
                 f"({e2.get('n_fam')} families) and B926 {f3(b2.get('cc_median'))}. R2's own test (cc at the detector's fixed placement, "
                 f"against the station read on other days) puts LRIV after 2019-05-30 at a median {f3(pc.get('cc_med_median'))} against a control "
                 f"95th percentile of {f3(pc.get('ctl_p95_median'))}; it flags LRIV as at control level in {pc.get('n_control_level')} of "
                 f"{pc.get('n_tested_months')} tested months. So LRIV is much weaker than B926 or PGC, but whether it is at control level depends "
                 f"on the event set. LRIV stays dropped by the v1.1 policy. (`{rel(O)}/r2/r2_positive_control.json`, `{rel(O)}/r2/pc_a2.json`)")
    if r3:
        S.append(f"**Cleaning A/B (R3).** Adopt the southern cleaning: **{r3['adopt_cleaning']}** (recall at matched control rate "
                 f"{r3['A_recall_at_matched_control']}, known-family veto {r3['V_known_family_veto']}, control falls more than data on the "
                 f"flagged months {r3['B_control_falls_more']}: data x{f3(r3['months_data_ratio'])}, control x{f3(r3['months_control_ratio'])}). "
                 + ("The production pass therefore applied the cleaning (this is the R4 re-run)." if r3["adopt_cleaning"] else
                    "R4 was skipped; the production pass kept the v1.1 data path and stores the cleaning flags as columns.")
                 + f" (`{rel(O)}/r3/r3_decision.json`)")
    if bars:
        rows = []
        for k, b in bars.items():
            rows.append(f"{k}: AUC {f3(b['auc_WAVE'])} vs net_mad {f3(b['auc_net_mad'])}, slope {f3(b['slope'], 2)}, matched/sumP "
                        f"{f3(b['matched_over_sumP'], 2)}, control P>0.5 {f3(100 * b['ctl_frac_P_gt05'], 2)} %")
        S.append(f"**Confidence (R5), the one-time {R.SEALED_YEAR} check** (opened {opened.get('opened_utc', '?') if opened else '?'}; nothing fitted on it). "
                 "Headline model WAVE: " + "; ".join(rows) + f". (`{rel(O)}/conf/sealed_bars.json`, `{rel(O)}/conf/eval_sealed.csv`)")
        weak = [k for k, b in bars.items() if not b.get("wave_beats_net_mad_by_0p05")]
        if weak:
            S.append(f"WAVE does not beat the plain network score by 0.05 AUC for {', '.join(weak)}: for those, the net_mad ranking with the "
                     "scrambled control is the honest product (Merlin's rule).")
    if gate:
        S.append(f"**PNET gate.** " + (f"dAUC(same) {f3(gate['dAUC'][0])} [{f3(gate['dAUC'][1])}, {f3(gate['dAUC'][2])}], slope {f3(gate['slope_with'], 2)} "
                                       f"-> {'ADOPTED: production PNET features computed' if gate.get('adopt') else 'not adopted: no production PNET'}"
                                       if gate.get("evaluated") else "not evaluated (no PNET features)") + f" (`{rel(O)}/conf/pnet_gate.json`)")
    if val:
        S.append(f"**Inter-ETS validation.** withheld S-energy {val['energy']['pass']}, control sum of P flat {val['control_sumP']['pass']}, "
                 f"station-set invariance {val['invariance'].get('pass')} -> inter-ETS rates reported as **{val['inter_ets_mode'].replace('_', ' ')}**. "
                 f"(`{rel(O)}/val/val_summary.json`)")
    if os.path.exists(f"{O}/rates/monthly_total.csv"):
        MT = pd.read_csv(f"{O}/rates/monthly_total.csv")
        SC = pd.read_csv(f"{O}/rates/step_checks.csv") if os.path.exists(f"{O}/rates/step_checks.csv") else pd.DataFrame()
        S.append(f"**Rates (R7).** Bostock-equivalent repeats per day, network total, median month {f3(MT.rate_per_day.median(), 1)} "
                 f"(range {f3(MT.rate_per_day.min(), 1)}-{f3(MT.rate_per_day.max(), 1)}); scrambled-control sum of P subtracted and shown. "
                 f"Step checks: {int(SC.flagged.sum()) if len(SC) else 0} of {len(SC)} flagged (|z| > 2 vs the same calendar windows of other years), "
                 f"{int((SC.flagged & SC.sse_adjacent).sum()) if len(SC) else 0} of them next to tremor days (not separable from activity). "
                 f"(`{rel(O)}/rates/monthly_total.csv`, `{rel(O)}/rates/step_checks.csv`)")
    cs = f"{R.CUT}/cut_summary.csv"
    if os.path.exists(cs):
        Cm = pd.read_csv(cs); vt = open(f"{R.CUT}/VERIFY.txt").read().strip() if os.path.exists(f"{R.CUT}/VERIFY.txt") else "no VERIFY.txt"
        st = J(f"{R.CUT}/STACKCHECK.json") or {}
        S.append(f"**Repeat windows (R8).** {int(Cm.n_windows.sum()):,} windows, {Cm.cut_GB.sum():.1f} GB from {Cm.raw_GB.sum():.0f} GB raw "
                 f"({len(Cm)} months). Verify: {vt[:120]}. Stack check: {st.get('n_pass')}/{st.get('n_families')} families pass. (`{rel(R.CUT)}/`)")
    elif pil:
        S.append(f"**Repeat windows (R8).** Pilot only: {pil['extrapolation']}")
    dl = f"{O}/r9/DELETE_CANDIDATES_svi.txt"
    if os.path.exists(dl):
        S.append("**Delete candidates (R9).** " + open(dl).read().splitlines()[1].strip() + ". Nothing was deleted; the list is input to a "
                 f"joint SVI + Washington decision. (`{rel(O)}/r9/`)")
    # ---------------------------------------------------------------- long sections
    L.append("## Stages and files\n")
    L.append(f"* R0 watcher: {r0w}  (`{rel(O)}/r0/watcher.json`); v1.1 hash snapshot `{rel(O)}/r0/v11_hashes.csv`")
    if r1:
        L.append(f"* R1: failing station-months (cache) {r1['fail_cache_list']}; (mf, re-run by the SVI_REF pass) {r1['fail_mf_list'][:20]}; "
                 f"energy-floor > 3 % {r1['masked_gt3_list'][:20]}  (`{rel(O)}/r1/`)")
    if pil:
        L.append(f"* R8 pilot: verify {pil['verify']}, stack check {pil['stackcheck']}, extrapolation {pil['extrapolation']}  (`{rel(O)}/cut_pilot/`)")
    if os.path.exists(f"{O}/r2/r2_health_monthly.csv"):
        M = pd.read_csv(f"{O}/r2/r2_health_monthly.csv")
        L.append(f"* R2: station-month status {M.status_eff.value_counts().to_dict()}; exclusions applied to the production pass: "
                 f"{len(pd.read_csv(f'{O}/r2/r2_exclusions.csv')) if os.path.exists(f'{O}/r2/r2_exclusions.csv') else 0}  "
                 f"(`{rel(O)}/r2/r2_health_monthly.csv`, `r2_effective_sets.csv`, `r2_swap_checks.csv`, `r2_bad_for_dvv.csv`)")
        if os.path.exists(f"{O}/r2/r2_swap_checks.csv"):
            Sw = pd.read_csv(f"{O}/r2/r2_swap_checks.csv")
            L.append("\n  CNSN swap checks (cc at placement, before -> after; control p95; best channel transform after):\n")
            for r in Sw.itertuples():
                L.append(f"  * {r.sta} {r.date}: {f3(r.cc_before)} -> {f3(r.cc_after)} (control {f3(r.ctl95_before)} / {f3(r.ctl95_after)}), "
                         f"lag {f3(r.lag_before, 2)} -> {f3(r.lag_after, 2)} s, identity transform best after: {r.perm_identity_after}, status after {r.status_after}")
    if r3:
        L.append(f"\n* R3 recall table: {r3['recall_table']}\n* R3 months: {r3['months']}")
    if frz:
        L.append(f"\n* Models frozen {frz['frozen_utc']} (sha256 {frz['sha256_models'][:16]}), fit rows {frz['n_fit']} on {frz['fit_days']} days; "
                 f"features `{rel(O)}/conf/MODELS_FROZEN.json`; never features: {frz['never_features']}")
    if os.path.exists(f"{O}/conf/eval_dev.csv"):
        E = pd.read_csv(f"{O}/conf/eval_dev.csv")
        E = E[E.model.isin(["WAVE", "FULL", "NO_FAM", "NO_DENSITY", "MF_ONLY", "NO_PNET", "WITH_PNET"])]
        L.append("\n### Development (out-of-fold, fit days) -- `results/refine/conf/eval_dev.csv`\n")
        L.append(E[["eval", "sset", "label", "model", "n", "n_pos", "auc", "auc_net_mad", "slope", "matched_over_sumP", "ctl_frac_P_gt05"]]
                 .round(3).to_markdown(index=False))
    if os.path.exists(f"{O}/conf/eval_sealed.csv"):
        E = pd.read_csv(f"{O}/conf/eval_sealed.csv")
        L.append(f"\n### Sealed {R.SEALED_YEAR} (one look) -- `results/refine/conf/eval_sealed.csv`\n")
        L.append(E[["sset", "label", "model", "n", "n_pos", "auc", "auc_net_mad", "slope", "matched_over_sumP", "ctl_frac_P_gt05"]]
                 .round(3).to_markdown(index=False))
    if val:
        L.append(f"\n### Inter-ETS validation\n\n```\n{json.dumps(val, indent=1, default=str)[:4000]}\n```")
    if os.path.exists(f"{O}/rates/k_factors.csv"):
        K = pd.read_csv(f"{O}/rates/k_factors.csv")
        L.append(f"\n* k = his / sumP per (family, set), shrunk: median by set {K.groupby('sset').k.median().round(3).to_dict()}; "
                 f"tau2 {K.groupby('sset').tau2.first().round(3).to_dict()}  (`{rel(O)}/rates/k_factors.csv`)")
    if os.path.exists(f"{O}/rates/step_checks.csv"):
        SC = pd.read_csv(f"{O}/rates/step_checks.csv")
        f = SC[SC.flagged & (SC.series == "rate")]
        L.append(f"\n### Flagged steps (network rate)\n\n" + (f[["date", "change", "window_d", "ratio_after_before", "z_vs_other_years",
                                                                   "tremor_days_in_window"]].round(3).to_markdown(index=False) if len(f) else "none"))
    L.append(f"\n* Runs used: {runs}")
    txt = ("# SVI_REF results (auto-written)\n\n## Plain summary\n\n" + "\n\n".join(S) + "\n\n" +
           "Caveat: P is the probability of matching Bostock's ETS-only catalogue. Between ETS windows it is a lower bound on 'real', "
           "not a measure of it; the inter-ETS checks above say how far the rates can be read there. 'Scrambled control' = the same "
           "detector run on the same day with each station site's data shifted to a different time of that day.\n\n" + "\n".join(L) + "\n")
    open(REPORT, "w").write(txt)
    R.log(f"[report] {REPORT} ({len(txt)} chars)")


if __name__ == "__main__":
    main()
