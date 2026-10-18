#!/usr/bin/env python
"""DUC_FIN report -> results/DUC_FIN_RESULTS.md (plain-language summary first; every number with its file)."""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
from t_common import RES, log

O = F.OUT; CO = F.CUT_OUT
REP = os.environ.get("DUC_FIN_REPORT", f"{RES}/DUC_FIN_RESULTS.md")


def rel(p): return os.path.relpath(p, os.path.dirname(REP))
def rd(p, **k): return pd.read_csv(p, **k) if os.path.exists(p) else None
def md(D, n=None): return "_(missing)_" if D is None else (D.head(n) if n else D).to_markdown(index=False, floatfmt=".3f")


def fig_monthly():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    J = rd(f"{O}/joins_by_month.csv"); M = rd(f"{O}/monthly.csv")
    if J is None: return None
    fig, ax = plt.subplots(figsize=(11, 3.8)); x = pd.to_datetime(J.month)
    ax.plot(x, J.his_rule, "k-", lw=2, label="his catalogue (his rule)")
    ax.plot(x, J.ours_rule, "C0-", label="ours, his rule (final scan)")
    ax.plot(x, J.ours_rule_mf1, "C0:", label="ours, his rule (pass 1, before re-templating)")
    if M is not None:
        m = M[M.fam == "all"].copy(); m["x"] = pd.to_datetime(m.month)
        days = m.days.to_numpy()
        ax.plot(m.x, m.sumP_same_rate * days, "C3-", label="sum P(same family) per month")
        ax.plot(m.x, m.ctl_sumP_same_rate * days, "C3--", label="scrambled control, sum P")
    ax.axvspan(pd.Timestamp("2011-05-01"), pd.Timestamp("2012-01-01"), color="0.9", zorder=0, label="his catalogue lost stations")
    ax.set_ylabel("events per month"); ax.legend(fontsize=7, ncol=2); ax.set_title("Ducellier families 2008-01..2012-07")
    fig.tight_layout(); os.makedirs(f"{O}/figures", exist_ok=True); p = f"{O}/figures/monthly_counts.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main():
    L = []; S = []
    # ---------------- numbers
    ES = rd(f"{O}/retpl/epoch_summary.csv"); TV = rd(f"{O}/retpl/template_vs_his.csv")
    JY = rd(f"{O}/joins_by_year.csv"); JM = rd(f"{O}/joins_by_month.csv")
    CS = rd(f"{CO}/cut_summary.csv")
    ver = open(f"{CO}/VERIFY.txt").read().strip() if os.path.exists(f"{CO}/VERIFY.txt") else "(missing)"
    stk = json.load(open(f"{CO}/STACKCHECK.json")) if os.path.exists(f"{CO}/STACKCHECK.json") else {}
    osc = json.load(open(f"{CO}/offset_source_counts.json")) if os.path.exists(f"{CO}/offset_source_counts.json") else {}
    AUC = rd(f"{O}/eval_auc.csv"); CAL = rd(f"{O}/eval_calibration.csv"); PER = rd(f"{O}/eval_permutation.csv")
    PST = rd(f"{O}/eval_pnet_station.csv"); BOO = rd(f"{O}/eval_bootstrap.csv"); WIT = rd(f"{O}/eval_within_famday.csv")
    TRN = rd(f"{O}/pnet/models/train_ft1.csv")
    dtxt = open(f"{F.OUT}/D2011.md").read() if os.path.exists(f"{F.OUT}/D2011.md") else ""
    v = dtxt.split("## Plain-language verdict (part A, catalogue files only)")[-1].split("## 1.")[0].strip() if dtxt else "(missing)"
    figm = fig_monthly()
    def auc(es, lab, m):
        if AUC is None: return np.nan
        r = AUC[(AUC.eval_set == es) & (AUC.label == lab)]
        return float(r[f"auc_{m}"].iloc[0]) if len(r) and f"auc_{m}" in r else np.nan
    def cal(es, lab, m, c):
        if CAL is None: return np.nan
        r = CAL[(CAL.eval_set == es) & (CAL.label == lab) & (CAL.model == m)]
        return float(r[c].iloc[0]) if len(r) else np.nan
    # ---------------- summary
    S.append("## Plain-language summary\n")
    S.append("**1. Why we find more events than Ducellier in May-Dec 2011 (decision point).** " +
             f"His catalogue lost stations; our scan did not (full diagnosis: [`D2011.md`]({rel(f'{O}/D2011.md')})).\n\n" + v + "\n")
    if ES is not None:
        od = ES[ES.on_disk]
        S.append(f"**2. Re-templating of the new 2011-10..12 station epochs** ([`retpl/epoch_summary.csv`]({rel(f'{O}/retpl/epoch_summary.csv')})): "
                 f"{len(od)} new epochs start inside the data on disk; {int((od.kept > 0).sum())} re-templated, {int((od.kept == 0).sum())} dropped "
                 f"(no family-channel beat its scrambled control); {int(od.kept.sum())} of {int(od.family_channels.sum())} family-channels kept. "
                 f"{int((~ES.on_disk).sum())} epochs start after 2012-07 and stay pending (no data on disk; no downloads in this run). "
                 + (f"Kept templates vs his template of the old channel: median cc {TV.cc_max.median():.2f}, median |cc| {TV.abscc_max.median():.2f}; "
                    f"{int((TV.sign_at_abs < 0).sum())} of {len(TV)} best matched with reversed polarity ([`retpl/template_vs_his.csv`]({rel(f'{O}/retpl/template_vs_his.csv')})). " if TV is not None and len(TV) else ""))
    if JM is not None:
        j12 = JM[JM.month >= "2012-01"]
        S.append(f"Templates per family-day in 2012: {j12.ntpl_famday_mf1.mean():.1f} before -> {j12.ntpl_famday.mean():.1f} after; "
                 f"events passing his rule 2012-01..07: {int(j12.ours_rule_mf1.sum())} before -> {int(j12.ours_rule.sum())} after "
                 f"([`joins_by_month.csv`]({rel(f'{O}/joins_by_month.csv')})).\n")
    if JY is not None:
        S.append("Ours / his under his rule (final scan): " + ", ".join(f"{int(r.year)} {r.ratio_rule:.3f}" for r in JY.itertuples() if np.isfinite(r.ratio_rule))
                 + f" ([`joins_by_year.csv`]({rel(f'{O}/joins_by_year.csv')})).\n")
    if CS is not None:
        S.append(f"**3. Windows for dv/v** ([`{rel(CO)}/`]({rel(CO)})): {int(CS.n_windows.sum()):,} windows (-15..+30 s, raw counts) around "
                 f"{int(pd.read_parquet(f'{CO}/detections.parquet', columns=['det_id']).shape[0]):,} candidates, {CS.cut_GB.sum():.1f} GB cut from "
                 f"{CS.raw_GB.sum():.0f} GB raw. Arrival offsets per window: {osc}. Check: {ver.split(':')[0]} ({ver}). "
                 f"Stack check: {stk.get('fam')} at {stk.get('netsta')} {stk.get('month')}, {stk.get('n_windows')} windows, stack vs his template "
                 f"cc {stk.get('cc_stack_vs_his_template')} (single windows {stk.get('cc_single_median')}, shifted control {stk.get('cc_control_shift_5_10s')}) "
                 f"-> {'PASS' if stk.get('pass') else 'FAIL'}.\n")
    if AUC is not None:
        S.append(f"**4. Per-repeat confidence P** (probability that a candidate is one of his catalogue events; fitted on 2008-2010 days, "
                 f"out-of-fold by day). Independent checks, AUC same family / any family: held-out months {auc('ho_month','same','FULL'):.3f} / "
                 f"{auc('ho_month','any','FULL'):.3f}, Jan-Apr 2011 {auc('ind2011a','same','FULL'):.3f} / {auc('ind2011a','any','FULL'):.3f}; "
                 f"his rule flag alone {auc('ho_month','same','pass_his'):.3f} / {auc('ho_month','any','pass_his'):.3f}, peak height (xmad) alone "
                 f"{auc('ho_month','same','xmad'):.3f}. Without PNET: {auc('ho_month','same','NO_PNET'):.3f} / {auc('ho_month','any','NO_PNET'):.3f}. "
                 f"Calibration slope (target 0.8-1.2) held-out months {cal('ho_month','same','FULL','slope'):.2f}, Jan-Apr 2011 "
                 f"{cal('ind2011a','same','FULL','slope'):.2f}; matched / sum P {cal('ho_month','same','FULL','obs_over_exp'):.2f} and "
                 f"{cal('ind2011a','same','FULL','obs_over_exp'):.2f}. Scrambled control: P > 0.5 for {cal('ho_month','same','FULL','ctl_fracP_gt05'):.1%} of "
                 f"control candidates vs {cal('ho_month','same','FULL','data_fracP_gt05'):.1%} of data candidates (held-out months). "
                 f"May-Dec 2011 (his catalogue incomplete, not fitted): AUC {auc('gap2011','same','FULL'):.3f}, matched / sum P "
                 f"{cal('gap2011','same','FULL','obs_over_exp'):.2f} (below 1 = we find more than his catalogue, as expected from step 1). "
                 f"Note: our detector reproduces his and his catalogue is his rule on the same statistic, so P is close to the "
                 f"his-rule flag; the gain is the ranking inside / near his threshold and a P that stays defined when channel counts change "
                 f"(2011-2012). Against his RAW detections (all 66 families) the FULL model reaches AUC {auc('ho_month','same_raw','FULL'):.3f} / "
                 f"{auc('ho_month','any_raw','FULL'):.3f} (same / any family). "
                 f"P is written for every candidate, 2012-01..07 included ([`P_by_candidate.parquet`]({rel(f'{O}/P_by_candidate.parquet')}), "
                 f"joinable [`conf_P_by_det_id.parquet`]({rel(f'{CO}/conf_P_by_det_id.parquet')})).\n")
    if TRN is not None:
        k = TRN.val_auc.idxmax()
        S.append(f"PNET fine-tune on his events: internal validation AUC {TRN.val_auc.iloc[0]:.3f} (off-the-shelf) -> {TRN.val_auc.max():.3f} "
                 f"(kept step {int(TRN.step.iloc[k])}); learning curve [`figures/pnet_learning_curve.png`]({rel(f'{O}/figures/pnet_learning_curve.png')}).\n")
    dl = f"{os.path.dirname(os.path.abspath(CO))}/DELETE_CANDIDATES_ducellier.txt"
    if os.path.exists(dl):
        t = open(dl).read().splitlines()
        S.append(f"**5. Delete candidates (nothing deleted)**: {t[1]}; {t[2]} ([`DELETE_CANDIDATES_ducellier.txt`]({rel(dl)})).\n")
    L += ["# DUC_FIN: Ducellier & Creager (2022) families finished on the data on disk (2008-01 .. 2012-07)", "",
          f"_Auto-written by `scripts/duc_fin_report.py`, {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%M} UTC; queue `scripts/duc_fin_queue.sh`; "
          f"outputs `{rel(O)}/` and `{rel(CO)}/`._", ""] + S
    if figm: L += [f"![monthly]({rel(figm)})", ""]
    # ---------------- sections
    L += ["## 1. 2011 diagnosis", "", f"See [`D2011.md`]({rel(f'{O}/D2011.md')}) (one-to-one matching by month, his channel counts, stations "
          "his catalogue lost, the extras' margin / clustering / co-occurrence, part B: extras re-summed on his station set and the scrambled control).", ""]
    L += ["## 2. Re-templating and the re-scan", "", md(ES), "", "Joins by year:", "", md(JY), "", "Joins by month (2011-2012):", "",
          md(JM[JM.month >= "2011-01"] if JM is not None else None), ""]
    L += ["## 3. Windows", "", md(CS), "", f"Station-day read guard: `{rel(f'{O}/status/GUARD_cut.txt')}`.", ""]
    L += ["## 4. Confidence", "", f"Day roles fixed before fitting ({rel(f'{O}/day_table.csv')}): held-out months {', '.join(F.HO_MONTHS)}; "
          f"PNET training days = every 10th 2008-2010 day (never fitted / evaluated); PNET held-out stations {', '.join(F.PNET_HO)} "
          "(NC.LBP is 1-component, fed as Z, 0, 0).", "", "AUC (selected columns; all in eval_auc.csv):", ""]
    if AUC is not None:
        cols = ["eval_set", "label", "n", "n_pos", "base_rate", "n_days", "auc_FULL", "ap_FULL", "auc_NO_PNET", "auc_PNET_W_ONLY", "auc_NO_FAM",
                "auc_MF_RULE", "auc_xmad", "auc_pass_his", "auc_margin", "auc_site_z_max", "auc_pnet_max_w", "auc_pnet_max_f", "auc_FULL_x"]
        L += [md(AUC[[c for c in cols if c in AUC]]), ""]
    if CAL is not None:
        L += ["Calibration and scrambled control (FULL, NO_PNET):", "", md(CAL[CAL.model.isin(["FULL", "NO_PNET"])]), ""]
    L += ["Within one family-day:", "", md(WIT), "", "Day-block bootstrap of AUC(FULL) - AUC(other):", "", md(BOO), ""]
    if PER is not None:
        L += ["Permutation importance (independent rows; blocks and top single features):", "",
              md(PER[PER.is_block].sort_values(["label", "auc_drop"], ascending=[True, False])), "",
              md(PER[~PER.is_block].sort_values("auc_drop", ascending=False), 15), ""]
    if PST is not None:
        P2 = PST.pivot_table(index=["netsta", "held_out_of_finetune"], columns="model", values=["auc_same", "auc_any"]).reset_index()
        P2.columns = ["_".join([c for c in col if c]) for col in P2.columns]
        L += ["PNET per station, before (off-the-shelf) / after (fine-tuned) -- AUC of the single-station value:", "", md(P2), ""]
    L += [f"![pnet learning curve]({rel(f'{O}/figures/pnet_learning_curve.png')})", ""]
    L += ["## 5. Delete candidates", ""]
    if os.path.exists(dl): L += ["```", open(dl).read().strip()[:6000], "```", ""]
    L += ["## Station-day read guards", ""] + [f"* {open(p).read().strip()}" for p in sorted(glob.glob(f"{O}/status/GUARD_*.txt"))] + \
         ([f"* {open(f'{O}/status/MFCHECK.txt').read().strip()}"] if os.path.exists(f"{O}/status/MFCHECK.txt") else []) + [""]
    open(REP, "w").write("\n".join(L) + "\n")
    log(f"[report] -> {REP}")


if __name__ == "__main__":
    main()
