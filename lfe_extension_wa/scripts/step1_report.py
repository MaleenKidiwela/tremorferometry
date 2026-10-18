#!/usr/bin/env python
"""STEP 1 report: results/STEP1_PROB_RESULTS.md + figures/step1_reliability.png from results/step1_prob/*.csv.
Numbers only in the tables; the plain-language summary is generated from the same numbers."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OD = os.environ.get("STEP1_OUT", f"{W.RES}/step1_prob")
FIG = f"{W.STUDY}/figures/step1_reliability.png" if "STEP1_OUT" not in os.environ else f"{OD}/step1_reliability.png"
MD = f"{W.RES}/STEP1_PROB_RESULTS.md" if "STEP1_OUT" not in os.environ else f"{OD}/STEP1_PROB_RESULTS.md"
REGN = {"ALL": "ALL", "Royer_NW_WA": "NW", "Royer_SW_WA": "SW", "Olympic_2010": "Olympic", "Olympic_2010_armB": "Arm B"}
SETS = ["V2010_test", "A2011", "V2008_test", "fit_in_sample"]
COL = {"own": "#2a78d6", "nbr3": "#eb6834", "any": "#1baf7a"}


def f3(x):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.3f}"


def f2(x):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2f}"


def figure(CAL):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sets = [s for s in SETS[:3] if s in set(CAL.eval_set)]
    fig, ax = plt.subplots(2, len(sets), figsize=(4.6 * len(sets), 8.2), squeeze=False)
    for j, es in enumerate(sets):
        a = ax[0, j]
        g = CAL[(CAL.eval_set == es) & (CAL.region == "ALL")]
        d = g[(g.series == "data") & (g.n >= 20)]
        c = g[(g.series == "control") & (g.n >= 20)]
        a.plot([0, 1], [0, 1], color="#999999", lw=1, ls="--", zorder=1)
        for lev, lab in [("own", "own family"), ("nbr3", "≤3-km group"), ("any", "any family")]:
            a.plot(d.mean_P, d[f"frac_{lev}"], "-o", color=COL[lev], lw=2, ms=6, label=f"data: {lab}", zorder=3)
        a.plot(c.mean_P, c.frac_any, "-s", color="#777777", lw=1.5, ms=5, label="scrambled control: any family", zorder=2)
        a.set_xlim(0, 1); a.set_ylim(0, 1)
        a.set_title(es.replace("_", " "), fontsize=11)
        a.set_xlabel("P(real) bin mean"); a.set_ylabel("fraction within ±2 s of an array detection")
        a.grid(color="#e6e6e6", lw=0.6); a.set_axisbelow(True)
        if j == 0:
            a.legend(fontsize=8, frameon=False, loc="upper left")
        b = ax[1, j]
        for s, colr, lab in [("data", "#2a78d6", "data"), ("control", "#777777", "scrambled control (4 realizations)")]:
            gg = g[g.series == s].sort_values("pbin")
            x = np.round(gg.pbin.str.split("-").str[0].astype(float) * 10).to_numpy() + (0.2 if s == "control" else -0.2)
            b.bar(x, gg.n / gg.n.sum(), width=0.38, color=colr, label=lab)
        b.set_yscale("log"); b.set_ylim(1e-5, 1.5)
        b.set_xticks(range(10)); b.set_xticklabels([f"{k/10:.1f}" for k in range(10)], fontsize=8)
        b.set_xlabel("P(real) bin (lower edge)"); b.set_ylabel("fraction of candidates")
        b.grid(axis="y", color="#e6e6e6", lw=0.6); b.set_axisbelow(True)
        if j == 0:
            b.legend(fontsize=8, frameon=False)
    fig.suptitle("STEP 1 — per-repeat P(real) by local FDR (scrambled control = noise model); truth-scanned family-days", fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(FIG), exist_ok=True)
    fig.savefig(FIG, dpi=130)


def main():
    info = json.load(open(f"{OD}/model_info.json"))
    CAL = pd.read_csv(f"{OD}/calibration.csv")
    SL = pd.read_csv(f"{OD}/calibration_slopes.csv")
    CQ = pd.read_csv(f"{OD}/control_and_quiet.csv")
    TL = pd.read_csv(f"{OD}/tremor_label_split.csv")
    EC = pd.read_csv(f"{OD}/expected_count.csv")
    FI = pd.read_csv(f"{OD}/feature_importance.csv")
    R = pd.read_parquet(f"{OD}/rate_family_day.parquet")
    figure(CAL)
    L = []
    L.append("# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)\n")
    L.append(f"Written {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%MZ} by `scripts/step1_report.py` (numbers only; interpretation after review). "
             "Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.\n")
    # ---------------- summary
    def sl(es, rg, s, col):
        r = SL[(SL.eval_set == es) & (SL.region == rg) & (SL.series == s)]
        return float(r[col].iloc[0]) if len(r) else np.nan

    def cq(es, rg, days, col):
        r = CQ[(CQ.eval_set == es) & (CQ.region == rg) & (CQ.days == days)]
        return float(r[col].iloc[0]) if len(r) else np.nan

    def ec(es, rg, col):
        r = EC[(EC.eval_set == es) & (EC.region == rg)]
        return float(r[col].iloc[0]) if len(r) else np.nan
    L.append("## Summary (plain language)\n")
    L.append("Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a "
             "real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each "
             "site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how "
             "data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, "
             "HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled "
             "by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half "
             "of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).\n")
    for es in ["V2010_test", "A2011", "V2008_test"]:
        L.append(f"- **{es.replace('_', ' ')}**: calibration slope own / ≤3 km / any = {f2(sl(es,'ALL','data','slope_own'))} / "
                 f"{f2(sl(es,'ALL','data','slope_nbr3'))} / {f2(sl(es,'ALL','data','slope_any'))} (target 0.8–1.2); matched / ΣP own "
                 f"{f2(sl(es,'ALL','data','obs_over_exp_own'))}, ≤3 km {f2(sl(es,'ALL','data','obs_over_exp_nbr3'))}. "
                 f"Scrambled control with P > 0.5: {f3(cq(es,'ALL','all','ctl_fracP_gt05'))} (never-seen control realization "
                 f"{f3(cq(es,'ALL','all','ctl_loco_fracP_gt05'))}); data with P > 0.5: {f3(cq(es,'ALL','all','data_fracP_gt05'))}. "
                 f"Quiet days: data P > 0.5 {f3(cq(es,'ALL','quiet','data_fracP_gt05'))} vs control {f3(cq(es,'ALL','quiet','ctl_fracP_gt05'))}. "
                 f"Expected count ΣP vs true array count per family-day: Spearman {f2(ec(es,'ALL','spearman_own'))} "
                 f"(raw candidate count {f2(ec(es,'ALL','spearman_ncand_own'))}), ratio ΣP / truth {f2(ec(es,'ALL','ratio_P_over_own'))}.")
    L.append("")
    L.append("Caveat for reading calibration: the array catalogue is the truth but is itself incomplete, and a real LFE of a "
             "neighbouring family counts as 'real' for P but matches only at the ≤ 3-km / any levels. Matched fraction below P "
             "can therefore mean real-but-uncatalogued; matched fraction ABOVE P means P is conservative.\n")
    # ---------------- method
    L.append("## Method\n")
    L.append(f"- Candidates: `results/mf_mt1_all` (T1 all ≥ 30-event build templates, w-weighted, 42 stations, `--drop-dead 0.05`), record "
             f"tier ≥ 8.0 MAD, series data + control + control2; plus `results/mf_step1c34` = the same run with seed 20261001 → "
             f"two more independent site-rolled scrambled controls (its data series is identical to mt1_all's: checked per day). "
             f"K = {info['K']} scrambled-control realizations.")
    L.append("- Features (all computable after the array era): " + ", ".join(f"`{f}`" for f in info["features"]) +
             ". Site values are the within-site w-weighted station cc at the candidate's aligned position (scrambled control: at the "
             "site's rolled position, exactly as the MF built it); the recomputed network cc equals the MF's to ≤ 2e-6 on every day. "
             "Site z = (site cc − median)/(1.4826 MAD) over 1000 random valid origins of the same family-day. `hfr_rel` = log10 "
             "E(10–19)/E(2–8) at the family-day's HF station minus the median of 300 random windows of that station-day. "
             "Densities count candidates of the same series.")
    L.append(f"- Density ratio: {info['method']} classifier data-vs-control on the fit set; f0/f = [(1−p)/p]·N_data/N_ctl.")
    L.append("- π0 rule: on quiet days (V2010 'quiet', V2008 'post'; tremor-free by the window design) every data candidate is "
             "noise → r_q = N_data(quiet) / (N_ctl(quiet)/K) per region (global if < 500 quiet data candidates), "
             "π0 = r_q·(N_ctl/K)/N_data, so P = 1 − r_q·(1−p)/(p·K), clipped to [0, 1]. Real LFEs on quiet days count as noise → conservative.")
    L.append("- Fit set: BUILD half of V2008 + V2010 (every series). Held out: V2010 / V2008 TEST half, A2011 (all). "
             "Leave-one-control-out refits (K = 3) give the model spread of the rate product and a P for a control realization "
             "the model never saw. Matching: within ±2 s (nearest, not one-to-one) of an array detection of the own family / "
             "any family ≤ 3 km (lat_use/lon_use) / any family, only on days the family's own array run scanned.\n")
    rq = info["r_q"]
    L.append("| region | quiet data candidates | quiet control candidates (all K) | r_q | source |")
    L.append("|---|---|---|---|---|")
    for rg in W.REGIONS + ["ALL"]:
        r = rq[rg]
        L.append(f"| {REGN.get(rg, rg)} | {r['n_data_quiet']} | {r['n_ctl_quiet']} | {r['r_q']:.3f} | {r['source']} |")
    L.append("")
    L.append("Feature ranking (fit set; permutation log-loss gain on 200k candidates; univariate AUC data vs control):\n")
    L.append("| feature | perm. log-loss gain | AUC data vs control | NaN frac |")
    L.append("|---|---|---|---|")
    for r in FI.itertuples():
        L.append(f"| {r.feature} | {r.perm_logloss_gain:.4f} | {f3(r.auc_data_vs_ctl)} | {r.frac_nan:.3f} |")
    L.append("")
    # ---------------- (a) calibration
    L.append("## (a) Calibration — fraction matched to the array catalogue by P bin\n")
    L.append("Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = "
             "calibration in the large. Control rows are the chance level.\n")
    L.append("| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |")
    L.append("|---|---|---|---|---|---|---|---|")
    for es in SETS:
        for rg in ["ALL"] + W.REGIONS:
            for s in ["data", "control"]:
                r = SL[(SL.eval_set == es) & (SL.region == rg) & (SL.series == s)]
                if not len(r):
                    continue
                r = r.iloc[0]
                L.append(f"| {es} | {REGN[rg]} | {s} | {int(r.n)} | {r.sum_P:.1f} | {f2(r.slope_own)} · {f2(r.slope_nbr3)} · {f2(r.slope_any)} | "
                         f"{f2(r.icpt_own)} · {f2(r.icpt_nbr3)} · {f2(r.icpt_any)} | "
                         f"{f2(r.obs_over_exp_own)} · {f2(r.obs_over_exp_nbr3)} · {f2(r.obs_over_exp_any)} |")
    L.append("")
    for es in ["V2010_test", "A2011", "V2008_test"]:
        g = CAL[(CAL.eval_set == es) & (CAL.region == "ALL")]
        if not len(g):
            continue
        L.append(f"Reliability table, {es}, ALL regions:\n")
        L.append("| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |")
        L.append("|---|---|---|---|---|---|---|---|")
        for b in sorted(set(g.pbin)):
            d = g[(g.series == "data") & (g.pbin == b)]; c = g[(g.series == "control") & (g.pbin == b)]
            dn = int(d.n.iloc[0]) if len(d) else 0; cn = int(c.n.iloc[0]) if len(c) else 0
            L.append(f"| {b} | {dn} | {f3(d.mean_P.iloc[0]) if dn else '–'} | {f3(d.frac_own.iloc[0]) if dn else '–'} | "
                     f"{f3(d.frac_nbr3.iloc[0]) if dn else '–'} | {f3(d.frac_any.iloc[0]) if dn else '–'} | {cn} | "
                     f"{f3(c.frac_any.iloc[0]) if cn else '–'} |")
        L.append("")
    L.append(f"![reliability](../figures/step1_reliability.png)\n")
    # ---------------- (b)+(c)+(d)
    L.append("## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region\n")
    L.append("`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = "
             "expected 'real' count the method would wrongly claim in one scrambled network.\n")
    L.append("| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for es in SETS + ["ALL"]:
        for rg in ["ALL"] + W.REGIONS:
            for dy in ["all", "quiet", "episode"]:
                r = CQ[(CQ.eval_set == es) & (CQ.region == rg) & (CQ.days == dy)]
                if not len(r):
                    continue
                r = r.iloc[0]
                L.append(f"| {es} | {REGN[rg]} | {dy} | {int(r.n_data)} | {f3(r.data_fracP_gt05)} | {f3(r.data_fracP_gt09)} | {r.data_sumP:.1f} | "
                         f"{r.n_ctl_per_real:.0f} | {f3(r.ctl_fracP_gt01)} | {f3(r.ctl_fracP_gt05)} | {f3(r.ctl_loco_fracP_gt05)} | "
                         f"{r.ctl_sumP_per_real:.1f} | {f2(r.ctl_loco_sumP_per_real)} |")
    L.append("")
    L.append("By PNSN regional tremor-day label (label only):\n")
    L.append("| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in TL.itertuples():
        L.append(f"| {r.eval_set} | {REGN[r.region]} | {r.tremor_day} | {r.n_famdays} | {f3(r.data_sumP_per_famday)} | {f3(r.ctl_sumP_per_famday)} | "
                 f"{f3(r.data_fracP_gt05)} | {f3(r.ctl_fracP_gt05)} |")
    L.append("")
    # ---------------- (e)
    L.append("## (e) Expected-count check — ΣP per family-day vs the true array count\n")
    L.append("Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array "
             "detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.\n")
    L.append("| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in EC.itertuples():
        L.append(f"| {r.eval_set} | {REGN[r.region]} | {r.n_famdays} | {r.sumP:.1f} | {r.sumP_ctl:.1f} | {r.truth_own} | {r.truth_nbr3} | "
                 f"{f3(r.ratio_P_over_own)} | {f3(r.ratio_P_over_nbr3)} | {f2(r.spearman_own)} | {f2(r.spearman_nbr3)} | {f2(r.spearman_ncand_own)} |")
    L.append("")
    # ---------------- rate product
    L.append("## Rate product — expected real repeats per family and day\n")
    L.append("`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): "
             "`sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model "
             "spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), "
             "truth counts where the family's array run scanned the day. Tremor-day label is a column only.\n")
    Rg = R.groupby(["region", "kind"]).agg(famdays=("fid", "size"), sumP_mean=("sumP", "mean"), sd_mean=("sd_total", "mean"),
                                         ctl_mean=("sumP_ctl_per_real", "mean"),
                                         fams_gt5=("sumP", lambda x: int((x >= 5).sum()))).reset_index()
    L.append("| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in Rg.itertuples():
        L.append(f"| {REGN[r.region]} | {r.kind} | {r.famdays} | {r.sumP_mean:.2f} | {r.sd_mean:.2f} | {r.ctl_mean:.2f} | {r.fams_gt5} |")
    L.append("")
    L.append("Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, "
             "`calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, "
             "`feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.")
    open(MD, "w").write("\n".join(L) + "\n")
    print("wrote", MD, FIG)


if __name__ == "__main__":
    main()
