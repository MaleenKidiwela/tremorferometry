#!/usr/bin/env python
"""Final stage: results/WA_TRIAL_RESULTS.md + figures/wa_trial_summary.png from the stage outputs."""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W


def f3(x):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.3f}"


def main():
    R = W.RES
    F = W.load_families()
    out = ["# WA permanent-station trial extension — results", "",
           f"Written {pd.Timestamp.utcnow():%Y-%m-%d %H:%MZ} by `scripts/wa_report.py`. Design: `PLAN.md`; log: `NOTES.md`.", ""]
    # ---------------- data
    fm = pd.read_csv(f"{R}/fetch_manifest.csv") if os.path.exists(f"{R}/fetch_manifest.csv") else pd.DataFrame()
    cm = pd.read_csv(f"{R}/cache_manifest.csv") if os.path.exists(f"{R}/cache_manifest.csv") else pd.DataFrame()
    out += ["## Data", ""]
    if len(fm):
        out.append(f"- Downloaded: {int((fm.status=='data').sum())} station-days, {fm.bytes.sum()/1e9:.1f} GB; "
                   f"nodata {int((fm.status=='nodata').sum())}, errors {int((fm.status=='error').sum())}.")
    if len(cm):
        out.append(f"- 20 Hz cache: {int(cm.status.isin(['ok','cached']).sum())} station-days "
                   f"({cm[cm.status.isin(['ok','cached'])].netsta.nunique()} stations).")
    out.append("")
    # ---------------- templates
    Q = pd.read_csv(f"{R}/templates/tpl_qc_build.csv")
    c = Q[Q.clear.fillna(False).astype(bool)]
    fs = c.groupby("fid").site.nunique()
    out += ["## Permanent templates (build half of the array-era detections)", "",
            f"- (family, station-epoch) stacks with >= 1 event: {len(Q)}; with >= 30 events: {int((Q.n_events>=30).sum())}; "
            f"clear arrival (imp > 5, dur < 3 s at 2–8 Hz, fs = 20 Hz): **{len(c)}**.",
            f"- Families with >= 1 clear permanent template: {c.fid.nunique()} of {len(F)}; with >= 2 clear SITES "
            f"(the observability rule): **{int((fs>=2).sum())}**.", ""]
    g = c.merge(F[["fid", "region"]]).groupby("region")
    t = pd.DataFrame(dict(n_fam=F.groupby("region").size(), fam_clear=g.fid.nunique(),
                          fam_2sites=c.merge(F[["fid", "region"]]).groupby(["region", "fid"]).site.nunique().ge(2).groupby("region").sum()))
    out += ["| region | families | ≥1 clear template | ≥2 clear sites |", "|---|---|---|---|"]
    for rg, r in t.fillna(0).iterrows():
        out.append(f"| {rg} | {int(r.n_fam)} | {int(r.fam_clear)} | {int(r.fam_2sites)} |")
    out += ["", "Clear templates per station: " + ", ".join(f"{k} {v}" for k, v in c.groupby("netsta").size().sort_values(ascending=False).items()), ""]
    # ---------------- HF cuts
    if os.path.exists(f"{R}/hf_cuts.csv"):
        H = pd.read_csv(f"{R}/hf_cuts.csv")
        out += ["## HF veto cuts (99th percentile of build-half own-matched ratios, per HF station node)", "",
                "| HF node | n calib | cut | median matched |", "|---|---|---|---|"]
        for r in H.sort_values("n_calib", ascending=False).itertuples():
            out.append(f"| {r.hf_node} | {r.n_calib} | {f3(r.cut)} | {f3(r.median_matched)} |")
        out.append("")
    # ---------------- validation
    p = f"{R}/score_val/validation_summary.csv"
    if os.path.exists(p):
        S = pd.read_csv(p)
        out += ["## Validation — the known answer (array-era detections of the same families)", "",
                "Test half only (2-h blocks, 60-s guards), observable family-days (≥ 2 permanent sites) the family's own "
                "array run scanned. One-to-one ±2 s, own family. `cat` = tier & HF pass & cross-family kept; `raw` = tier only. "
                "Scrambled control = site-rolled series (each site shifted by its own 300–3000 s), detections per family-day at "
                "the same tier (HF veto / cross-family rule not applicable to it).", "",
                "`control prec` = fraction of scrambled-control detections (test half) that land within ±2 s of an own-family array "
                "detection by chance (the precision baseline); `control recall` = array detections matched by the control.", "",
                "| window | region | fam rows | fam obs | truth (test) | recall core | prec own core | prec any core | control recall / prec core | recall ext | prec own ext | data/famday core | control/famday core | quiet data / control core |",
                "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in S.itertuples():
            out.append(f"| {r.window} | {r.region} | {r.n_fam} | {r.n_fam_observable} | {r.n_truth_test} | "
                       f"{f3(r.recall_cat_core)} | {f3(r.prec_own_cat_core)} | {f3(r.prec_any_cat_core)} | "
                       f"{f3(r.recall_ctl_core)} / {f3(r.prec_own_ctl_core)} | "
                       f"{f3(r.recall_cat_ext)} | {f3(r.prec_own_cat_ext)} | {f3(r.data_per_famday_core)} | "
                       f"{f3(r.ctl_per_famday_core)} | {f3(r.quiet_data_per_famday_core)} / {f3(r.quiet_ctl_per_famday_core)} |")
        TT = pd.read_csv(f"{R}/score_val/tiers_by_epoch.csv") if os.path.exists(f"{R}/score_val/tiers_by_epoch.csv") else pd.DataFrame()
        out += ["", "**Per-network-epoch tiers** (SVI `c_tiers_by_epoch.py` rule, FA only: T_core = smallest T ≥ 8.0 on a 0.25 grid with "
                "quiet-day FA ≤ 0.1 per template-day, T_ext ≤ 1.0; quiet = the window's non-episode observable family-days; "
                "V2008's are post-episode days with no tremor catalogue):", "",
                "| window | quiet fam-days | FA at 8.25 | FA at 8.75 | T_ext_ep | T_core_ep |", "|---|---|---|---|---|---|"]
        for r in TT.itertuples():
            out.append(f"| {r.window} | {r.n_quiet_famdays} | {f3(r.fa_at_8p25)} | {f3(r.fa_at_8p75)} | {r.T_ext_ep} | {r.T_core_ep} |")
        out += ["", "| window | region | recall core_ep | prec own core_ep | prec any core_ep | control recall / prec core_ep | recall ext_ep | prec own ext_ep | data/famday core_ep | control/famday core_ep |",
                "|---|---|---|---|---|---|---|---|---|---|"]
        for r in S.itertuples():
            out.append(f"| {r.window} | {r.region} | {f3(r.recall_cat_core_ep)} | {f3(r.prec_own_cat_core_ep)} | {f3(r.prec_any_cat_core_ep)} | "
                       f"{f3(r.recall_ctl_core_ep)} / {f3(r.prec_own_ctl_core_ep)} | {f3(r.recall_cat_ext_ep)} | {f3(r.prec_own_cat_ext_ep)} | "
                       f"{f3(r.data_per_famday_core_ep)} | {f3(r.ctl_per_famday_core_ep)} |")
        V = pd.read_csv(f"{R}/score_val/validation_by_family.csv")
        v = V[V.n_truth_test >= 20].sort_values("recall_cat_core", ascending=False)
        out += ["", f"Per family: `results/score_val/validation_by_family.csv`. Families with ≥ 20 test-half array detections "
                f"on observable days: {len(v)}; recall(core) median {f3(v.recall_cat_core.median())}, "
                f"≥ 0.2: {int((v.recall_cat_core>=0.2).sum())}, ≥ 0.5: {int((v.recall_cat_core>=0.5).sum())}.", "",
                "Top 15 by recall (core):", "",
                "| family | window | med sites | truth test | recall core | prec own core | data/famday | control/famday |", "|---|---|---|---|---|---|---|---|"]
        for r in v.head(15).itertuples():
            out.append(f"| {r.family_uid} | {r.window} | {r.med_sites:.0f} | {r.n_truth_test} | {f3(r.recall_cat_core)} | "
                       f"{f3(r.prec_own_cat_core)} | {f3(r.data_per_famday_core)} | {f3(r.ctl_per_famday_core)} |")
        out.append("")
    # ---------------- trial
    tabs = {}
    for lab, mf in [("pass 1 (validation-epoch templates only)", "t2012p1"), ("pass 2 (+ re-templated new epochs)", "t2012p2")]:
        p = f"{R}/score_{mf}/trial_by_family.csv"
        if os.path.exists(p):
            tabs[lab] = pd.read_csv(p)
    if tabs:
        nn = pd.read_csv(f"{R}/trial/new_nodes.csv") if os.path.exists(f"{R}/trial/new_nodes.csv") else pd.DataFrame()
        out += ["## Trial — 2012-09 ETS, after the arrays (XU ended 2012-05-31, XG 2011-09)", "",
                "Window 2012-08-20..24 (quiet), 2012-09-05..10-13 (episode), 2012-10-18..22 (quiet). No truth: tremor days "
                "(PNSN, the region's tremor-day table, ≥ 1 window) are a LABEL, never a gate.", ""]
        if len(nn):
            out.append("New station epochs re-templated from pass-1 detections: " + ", ".join(nn.node) + ".")
            if os.path.exists(f"{R}/templates/tpl_qc_t2012new.csv"):
                qn = pd.read_csv(f"{R}/templates/tpl_qc_t2012new.csv")
                out.append(f"  → {len(qn)} stacks, {int(qn.clear.fillna(False).astype(bool).sum())} clear templates.")
            out.append("")
        out += ["| pass | region | families | with ≥2 live permanent sites | cat core | on tremor days | other days | rate/famday tremor | rate/famday other | control rate tremor | control rate other |",
                "|---|---|---|---|---|---|---|---|---|---|---|"]
        for tier in ["core", "core_ep"]:
          for lab, Tt in tabs.items():
            for rg, g in [("ALL", Tt)] + list(Tt.groupby("region")):
                o = g[g.n_days_obs > 0]
                ntd, nod = o.n_tremor_days.sum(), o.n_other_days.sum()
                ctr = (o[f"ctl_rate_{tier}_tremor"].fillna(0) * o.n_tremor_days).sum() / ntd if ntd else np.nan
                cot = (o[f"ctl_rate_{tier}_other"].fillna(0) * o.n_other_days).sum() / nod if nod else np.nan
                out.append(f"| {lab} [{tier}] | {rg} | {len(g)} | {len(o)} | {int(g[f'n_cat_{tier}'].sum())} | {int(g[f'n_cat_{tier}_tremor'].sum())} | "
                           f"{int(g[f'n_cat_{tier}_other'].sum())} | {f3(g[f'n_cat_{tier}_tremor'].sum()/ntd if ntd else np.nan)} | "
                           f"{f3(g[f'n_cat_{tier}_other'].sum()/nod if nod else np.nan)} | {f3(ctr)} | {f3(cot)} |")
        for mf in ["t2012p1", "t2012p2"]:
            pt = f"{R}/score_{mf}/tiers_by_epoch.csv"
            if os.path.exists(pt):
                x = pd.read_csv(pt).iloc[0]
                out.append(f"")
                out.append(f"{mf} epoch tiers (trial quiet days, {x.n_quiet_famdays} family-days): FA at 8.75 = {f3(x.fa_at_8p75)}, "
                           f"T_ext_ep {x.T_ext_ep}, T_core_ep {x.T_core_ep}.")
        last = list(tabs.values())[-1]
        o = last[last.n_days_obs > 0].sort_values("n_cat_core", ascending=False)
        out += ["", f"Per family: `results/score_t2012p{len(tabs)}/trial_by_family.csv`. Top 20 by catalogue count (core):", "",
                "| family | region | obs days | med sites | cat core | tremor days / rate | other days / rate | control rate tremor / other |", "|---|---|---|---|---|---|---|---|"]
        for r in o.head(20).itertuples():
            out.append(f"| {r.family_uid} | {r.region} | {r.n_days_obs} | {r.med_sites_obs:.0f} | {r.n_cat_core} | "
                       f"{r.n_tremor_days} / {f3(r.rate_core_tremor)} | {r.n_other_days} / {f3(r.rate_core_other)} | "
                       f"{f3(r.ctl_rate_core_tremor)} / {f3(r.ctl_rate_core_other)} |")
        out.append("")
    open(f"{R}/WA_TRIAL_RESULTS.md", "w").write("\n".join(out) + "\n")
    # ---------------- figure
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        if os.path.exists(f"{R}/score_val/validation_by_family.csv"):
            V = pd.read_csv(f"{R}/score_val/validation_by_family.csv")
            V = V[V.n_truth_test >= 20]
            for rg, g in V.groupby("region"):
                ax[0].scatter(g.med_sites, g.recall_cat_core, s=14, label=rg, alpha=.8)
            ax[0].set_xlabel("median permanent sites (observable days)"); ax[0].set_ylabel("recall (T_core catalogue, test half)")
            ax[0].set_title("Validation: permanent-only vs array detections"); ax[0].legend(fontsize=7)
        if tabs:
            Tt = list(tabs.values())[-1]; o = Tt[Tt.n_days_obs > 0]
            ax[1].scatter(o.rate_core_other, o.rate_core_tremor, s=14, alpha=.8, label="data")
            ax[1].scatter(o.ctl_rate_core_other, o.ctl_rate_core_tremor, s=10, alpha=.6, marker="x", label="scrambled control")
            m = np.nanmax([o.rate_core_tremor.max(), o.rate_core_other.max(), 1])
            ax[1].plot([0, m], [0, m], "k:", lw=.8)
            ax[1].set_xlabel("detections / family-day, other days"); ax[1].set_ylabel("detections / family-day, tremor days")
            ax[1].set_title("Trial 2012-09 (labels only)"); ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f"{W.STUDY}/figures/wa_trial_summary.png", dpi=130)
    except Exception as e:
        print("figure failed", e)
    print("wrote", f"{R}/WA_TRIAL_RESULTS.md")


if __name__ == "__main__":
    main()
