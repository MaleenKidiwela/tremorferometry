#!/usr/bin/env python
"""Tables, figures and results/stack_activity/STACK_ACTIVITY_RESULTS.md from the scored stack_activity runs.

Validation (known answer = the families' own array detections): per run / window label / region / W
  AUC of 'active window' (n_truth >= N_ACTIVE = 5 test-half detections; N = 1 and 20 as sensitivity columns)
  vs the activity score A (Stouffer z over stations vs random times), for DATA and for the SCRAMBLED CONTROL
  (control1 scored exactly the same way), plus A_w (network-weighted) and A_clear (clear templates only);
  pooled Spearman(score, n_truth); median per-family Spearman; fraction of windows with p_ctl <= 0.05.
Trials (2012-2014, no truth): the same tables with the PNSN tremor-day LABEL (window has >= 1 tremor day of the
  family's region) in place of truth — a label, never a gate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
from stack_activity import OUT

N_ACTIVE = 5
RUNS = [("valA", "validation", "V2008 + V2010, 42 permanent stations, build-half templates, TEST half"),
        ("valB", "validation", "A2011 (independent episode: templates are V2008/V2010), 42-station build templates, all halves"),
        ("valC", "validation", "V2010 with the step-4b stations (build + step4b_a templates), TEST half"),
        ("valS", "validation", "V2008 + V2010, continuous-series configuration (selected families + stations), TEST half"),
        ("valSB", "validation", "A2011, continuous-series configuration"),
        ("tr2012", "trial", "T2012 (2012-08-20..10-22), build + step4b_a templates, 94 stations"),
        ("tr2013", "trial", "T2013 (2013-08-30..10-18)"),
        ("tr2014", "trial", "T2014 (2014-11-10..12-19)"),
        ("cont", "trial", "continuous 2012-06-01..2014-12-31, selected families + stations")]
C_DATA, C_CTL, C_INK, C_MUTED = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def auc(score, lab):
    from scipy.stats import rankdata
    m = np.isfinite(score) & np.isfinite(lab)
    s = score[m]; y = lab[m].astype(bool)
    n1, n0 = y.sum(), (~y).sum()
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(s)
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def spear(a, b):
    from scipy.stats import spearmanr
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 4 or np.nanstd(a[m]) == 0 or np.nanstd(b[m]) == 0:
        return np.nan
    return float(spearmanr(a[m], b[m])[0])


def load(run):
    p = f"{OUT}/score_{run}/family_window.parquet"
    return pd.read_parquet(p) if os.path.exists(p) else None


def metrics(g, target, active):
    row = dict(n_famwin=len(g), n_families=g.fid.nunique(), n_active=int(active.sum()))
    for sc in ["A", "A_w", "A_clear"]:
        for s, lab in [("data", "data"), ("control1", "ctl")]:
            row[f"auc_{sc}_{lab}"] = auc(g[f"{sc}_{s}"].to_numpy(), active.to_numpy().astype(float))
    for s, lab in [("data", "data"), ("control1", "ctl")]:
        row[f"spearman_A_{lab}"] = spear(g[f"A_{s}"].to_numpy(), target.to_numpy().astype(float))
        pf = [spear(x[f"A_{s}"].to_numpy(), target[x.index].to_numpy().astype(float)) for _, x in g.groupby("fid")]
        row[f"median_family_spearman_{lab}"] = float(np.nanmedian(pf)) if np.isfinite(pf).any() else np.nan
        pp = g[f"p_A_{s}"]
        row[f"frac_p05_active_{lab}"] = float((pp[active] <= 0.05).mean()) if active.any() else np.nan
        row[f"frac_p05_inactive_{lab}"] = float((pp[~active] <= 0.05).mean()) if (~active).any() else np.nan
    return row


def main():
    od = OUT
    figd = f"{W.STUDY}/figures"
    os.makedirs(figd, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    F = pd.read_csv(f"{W.RES}/setup/families.csv")
    UID = dict(zip(F.fid, F.family_uid))
    have = {r: load(r) for r, _, _ in RUNS}
    md = ["# Family activity by stacking — results", "",
          f"Written {pd.Timestamp.utcnow():%Y-%m-%d %H:%MZ} by `scripts/stack_activity_report.py`. Design: `PLAN.md` "
          "§Stack activity; engine `scripts/stack_activity.py`, scoring `scripts/stack_activity_score.py`.", "",
          "**Statistic.** Per family, station-epoch template (ALL >= 30-event array-era stacks; clear flag a column) and "
          "day: the top-K (2 per scored hour) declustered (+-4 s) peaks of the family's network matched filter computed "
          "WITHOUT the station's own site (leave-one-site-out); the station's 2-8 Hz waveforms (filter designed at fs = 20 Hz) "
          "are stacked at those times (each window unit-max normalised); rho = stack-template correlation. Windows of "
          "W = 1, 5, 10, 30 days = sums of daily stacks. z = (rho - mean of 6 random-times stacks) / pooled random-times SD; "
          "family activity **A** = sum_s z_s / sqrt(n_s) (Stouffer over stations; A_w network-weighted; A_clear clear "
          "templates only). The SCRAMBLED CONTROL (3 independent site-rolled series, 300-3000 s per site per day, candidates "
          "picked on the rolled network) gets exactly the same statistic; p_ctl = empirical p of A against the family's "
          "control pool (3 controls x all windows at that W; each control is scored against the other two).", "",
          f"**Active window** = >= {N_ACTIVE} of the family's own array detections in the window (test half where the run "
          "scored the test half; only days its array run scanned). N = 1 and 20 in the sensitivity columns. Nothing tuned.", ""]
    # ---------------- validation tables
    VT, FT = [], []
    for run, mode, desc in RUNS:
        FW = have[run]
        if FW is None or mode != "validation":
            continue
        FW = FW[FW.n_truth.notna() & (FW.n_scanned_days > 0)].copy()
        for (lab, w), g0 in FW.groupby(["label", "W"]):
            for rg in ["ALL"] + sorted(g0.region.dropna().unique()):
                g = g0 if rg == "ALL" else g0[g0.region == rg]
                if len(g) < 5:
                    continue
                row = dict(run=run, label=lab, region=rg, W=w)
                row.update(metrics(g, g.n_truth, g.n_truth >= N_ACTIVE))
                for nn in (1, 20):
                    act = g.n_truth >= nn
                    row[f"auc_A_data_N{nn}"] = auc(g.A_data.to_numpy(), act.to_numpy().astype(float))
                    row[f"auc_A_ctl_N{nn}"] = auc(g.A_control1.to_numpy(), act.to_numpy().astype(float))
                VT.append(row)
        for (w, f), g in FW.groupby(["W", "fid"]):
            act = g.n_truth >= N_ACTIVE
            FT.append(dict(run=run, W=w, fid=f, family_uid=UID.get(f), region=g.region.iloc[0], n_windows=len(g),
                           n_active=int(act.sum()), truth_total=int(g.n_truth.sum()),
                           med_n_sta=float(g.n_sta.median()), med_n_sta_clear=float(g.n_sta_clear.median()),
                           auc_A_data=auc(g.A_data.to_numpy(), act.to_numpy().astype(float)),
                           auc_A_ctl=auc(g.A_control1.to_numpy(), act.to_numpy().astype(float)),
                           auc_A_clear_data=auc(g.A_clear_data.to_numpy(), act.to_numpy().astype(float)),
                           spearman_A_data=spear(g.A_data.to_numpy(), g.n_truth.to_numpy()),
                           spearman_A_ctl=spear(g.A_control1.to_numpy(), g.n_truth.to_numpy()),
                           med_A_active=float(g.A_data[act].median()) if act.any() else np.nan,
                           med_A_inactive=float(g.A_data[~act].median()) if (~act).any() else np.nan,
                           cand_truth_frac_data=float(g.get("cand_truth_frac_data", pd.Series(dtype=float)).mean()),
                           cand_truth_frac_ctl=float(g.get("cand_truth_frac_control1", pd.Series(dtype=float)).mean())))
    VT = pd.DataFrame(VT); FT = pd.DataFrame(FT)
    VT.to_csv(f"{od}/validation_by_W.csv", index=False)
    FT.to_csv(f"{od}/validation_by_family.csv", index=False)
    md += ["## Validation — the known answer", ""]
    for run, mode, desc in RUNS:
        if mode == "validation" and have[run] is not None:
            md.append(f"- `{run}`: {desc}")
    md.append("")
    if len(VT):
        cols = ["run", "label", "region", "W", "n_famwin", "n_families", "n_active", "auc_A_data", "auc_A_ctl",
                "auc_A_w_data", "auc_A_w_ctl", "auc_A_clear_data", "auc_A_clear_ctl", "spearman_A_data", "spearman_A_ctl",
                "median_family_spearman_data", "median_family_spearman_ctl", "frac_p05_active_data", "frac_p05_inactive_data",
                "frac_p05_active_ctl", "frac_p05_inactive_ctl", "auc_A_data_N1", "auc_A_ctl_N1", "auc_A_data_N20", "auc_A_ctl_N20"]
        md += ["`ctl` = control1 scored identically (its AUC/correlation against the same truth is the chance level). "
               "`frac_p05_*` = fraction of windows with p_ctl <= 0.05 (the control's own rate is its false-alarm rate).", "",
               "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for r in VT.itertuples():
            md.append("| " + " | ".join(f"{getattr(r, c):.3f}" if isinstance(getattr(r, c), float) else str(getattr(r, c))
                                         for c in cols) + " |")
        md.append("")
        if len(FT):
            md += ["### Per family (W = 1, 5, 30), best by AUC of A (data) — full table `validation_by_family.csv`", ""]
            for w in (1, 5, 30):
                x = FT[(FT.W == w) & (FT.n_active >= 2) & (FT.n_windows - FT.n_active >= 2)]
                if not len(x):
                    continue
                md.append(f"W = {w}: {len(x)} family-runs with >= 2 active and >= 2 inactive windows; AUC(data) median "
                          f"{x.auc_A_data.median():.3f} (control {x.auc_A_ctl.median():.3f}); AUC >= 0.8: data "
                          f"{(x.auc_A_data >= 0.8).sum()}, control {(x.auc_A_ctl >= 0.8).sum()}; Spearman(A, n_truth) median "
                          f"data {x.spearman_A_data.median():.3f}, control {x.spearman_A_ctl.median():.3f}.")
                md.append("")
                c2 = ["run", "family_uid", "region", "n_windows", "n_active", "truth_total", "med_n_sta", "med_n_sta_clear",
                      "auc_A_data", "auc_A_ctl", "spearman_A_data", "spearman_A_ctl", "med_A_active", "med_A_inactive",
                      "cand_truth_frac_data", "cand_truth_frac_ctl"]
                md += ["| " + " | ".join(c2) + " |", "|" + "---|" * len(c2)]
                for r in x.sort_values("auc_A_data", ascending=False).head(15).itertuples():
                    md.append("| " + " | ".join(f"{getattr(r, c):.3f}" if isinstance(getattr(r, c), float) else str(getattr(r, c))
                                                 for c in c2) + " |")
                md.append("")
    # ---------------- trials (label only)
    TT = []
    for run, mode, desc in RUNS:
        FW = have[run]
        if FW is None or mode != "trial":
            continue
        for w, g0 in FW.groupby("W"):
            for rg in ["ALL"] + sorted(g0.region.dropna().unique()):
                g = g0 if rg == "ALL" else g0[g0.region == rg]
                if len(g) < 5:
                    continue
                lab = g.n_tremor_days > 0
                row = dict(run=run, region=rg, W=w, n_famwin=len(g), n_families=g.fid.nunique(),
                           n_tremor_label=int(lab.sum()),
                           auc_tremor_A_data=auc(g.A_data.to_numpy(), lab.to_numpy().astype(float)),
                           auc_tremor_A_ctl=auc(g.A_control1.to_numpy(), lab.to_numpy().astype(float)),
                           med_A_data=float(g.A_data.median()), med_A_ctl=float(g.A_control1.median()),
                           q95_A_ctl=float(np.nanpercentile(g[["A_control1", "A_control2", "A_control3"]].to_numpy(), 95)),
                           frac_A_data_gt_ctl_q95=float((g.A_data > np.nanpercentile(
                               g[["A_control1", "A_control2", "A_control3"]].to_numpy(), 95)).mean()),
                           frac_p05_data=float((g.p_A_data <= 0.05).mean()), frac_p05_ctl=float((g.p_A_control1 <= 0.05).mean()))
                TT.append(row)
    TT = pd.DataFrame(TT)
    TT.to_csv(f"{od}/trial_by_W.csv", index=False)
    if len(TT):
        md += ["## Trials 2012-2014 — no truth; tremor days are a LABEL", "",
               "`auc_tremor_*` = AUC of A for windows containing >= 1 PNSN tremor day of the family's region (a label; the "
               "control's value is its chance level). `frac_A_data_gt_ctl_q95` = fraction of data windows above the pooled "
               "control 95th percentile (5 % expected by chance).", ""]
        cols = list(TT.columns)
        md += ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for r in TT.itertuples(index=False):
            md.append("| " + " | ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in r) + " |")
        md.append("")
    # ---------------- figures
    figs = []
    if len(VT):
        fig, axs = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
        V = VT[VT.region == "ALL"]
        for (run, lab), g in V.groupby(["run", "label"]):
            g = g.sort_values("W")
            axs[0].plot(g.W, g.auc_A_data, "-o", color=C_DATA, lw=2, ms=6, alpha=0.9)
            axs[0].plot(g.W, g.auc_A_ctl, "--o", color=C_CTL, lw=2, ms=6, alpha=0.9)
            axs[0].annotate(f"{run}:{lab}", (g.W.iloc[-1], g.auc_A_data.iloc[-1]), fontsize=7, color=C_MUTED,
                            xytext=(3, 0), textcoords="offset points")
            axs[1].plot(g.W, g.median_family_spearman_data, "-o", color=C_DATA, lw=2, ms=6)
            axs[1].plot(g.W, g.median_family_spearman_ctl, "--o", color=C_CTL, lw=2, ms=6)
        for ax, t in zip(axs, [f"AUC, active window (>= {N_ACTIVE} array detections) vs A", "median per-family Spearman(A, n_truth)"]):
            ax.set_xscale("log"); ax.set_xticks([1, 5, 10, 30]); ax.set_xticklabels(["1", "5", "10", "30"])
            ax.set_xlabel("window W (days)"); ax.set_title(t, fontsize=10); ax.grid(alpha=0.2)
        axs[0].axhline(0.5, color=C_MUTED, lw=0.8)
        axs[1].axhline(0.0, color=C_MUTED, lw=0.8)
        axs[0].plot([], [], "-o", color=C_DATA, label="data"); axs[0].plot([], [], "--o", color=C_CTL, label="scrambled control")
        axs[0].legend(fontsize=8, frameon=False)
        p = f"{figd}/stack_activity_validation.png"; fig.savefig(p, dpi=130); plt.close(fig); figs.append(p)
    # control distributions
    fig, axs = plt.subplots(2, 4, figsize=(14, 6), constrained_layout=True)
    for i, mode in enumerate(["validation", "trial"]):
        X = pd.concat([have[r] for r, m, _ in RUNS if m == mode and have[r] is not None and r != "cont"] or [pd.DataFrame()])
        for j, w in enumerate([1, 5, 10, 30]):
            ax = axs[i, j]
            if not len(X) or "W" not in X:
                ax.axis("off"); continue
            g = X[X.W == w]
            d = g.A_data.dropna(); c = pd.concat([g.A_control1, g.A_control2, g.A_control3]).dropna()
            if not len(d):
                ax.axis("off"); continue
            bins = np.linspace(min(d.quantile(0.005), c.quantile(0.005)), max(d.quantile(0.995), c.quantile(0.995)), 50)
            ax.hist(c, bins, density=True, histtype="step", color=C_CTL, lw=2, label="scrambled control (3)")
            ax.hist(d, bins, density=True, histtype="step", color=C_DATA, lw=2, label="data")
            ax.set_title(f"{mode}, W = {w} d  (n = {len(d)})", fontsize=9); ax.set_xlabel("A"); ax.grid(alpha=0.2)
            if i == 0 and j == 0:
                ax.legend(fontsize=7, frameon=False)
    p = f"{figd}/stack_activity_control_dist.png"; fig.savefig(p, dpi=130); plt.close(fig); figs.append(p)
    # time series of the top families (validation W = 1 + trials W = 5)
    if len(FT):
        top = (FT[(FT.W == 1) & (FT.run == FT.run.iloc[0]) & (FT.n_active >= 3) & (FT.n_windows - FT.n_active >= 3)]
               .sort_values("auc_A_data", ascending=False).fid.head(6).tolist())
        if not top:
            top = FT[FT.W == 1].sort_values("auc_A_data", ascending=False).fid.head(6).tolist()
        runs_ts = [r for r in ["valA", "valB", "tr2012", "tr2013", "tr2014"] if have.get(r) is not None]
        segs = []
        for r in runs_ts:
            for lab in sorted(have[r].label.unique()):
                segs.append((r, lab))
        if top and segs:
            fig, axs = plt.subplots(2 * len(top), len(segs), figsize=(3.2 * len(segs), 2.3 * len(top)),
                                    squeeze=False, constrained_layout=True,
                                    gridspec_kw=dict(height_ratios=[2, 1] * len(top)))
            for i, f in enumerate(top):
                for j, (r, lab) in enumerate(segs):
                    ax, bx = axs[2 * i, j], axs[2 * i + 1, j]
                    g = have[r]
                    g = g[(g.fid == f) & (g.label == lab) & (g.W == 1)].sort_values("start")
                    if not len(g):
                        ax.axis("off"); bx.axis("off"); continue
                    t = pd.to_datetime(g.start)
                    for k in (1, 2, 3):
                        ax.plot(t, g[f"A_control{k}"], "-", color=C_CTL, lw=1, alpha=0.5, label="control" if k == 1 else None)
                    ax.plot(t, g.A_data, "-o", color=C_DATA, lw=2, ms=3, label="data")
                    ax.axhline(0, color=C_MUTED, lw=0.6)
                    ax.set_title(f"{UID.get(f)} — {r}:{lab}", fontsize=8)
                    if "n_truth" in g and g.n_truth.notna().any():
                        bx.bar(t, g.n_truth.fillna(0), color=C_MUTED, width=0.8)
                        bx.set_ylabel("array det.", fontsize=7)
                    else:
                        bx.bar(t, g.n_tremor_days, color=C_MUTED, width=0.8)
                        bx.set_ylabel("tremor day", fontsize=7)
                    for x_ in (ax, bx):
                        x_.tick_params(labelsize=6)
                    if i == 0 and j == 0:
                        ax.legend(fontsize=6, frameon=False)
                    if j == 0:
                        ax.set_ylabel("A (W = 1 d)", fontsize=7)
            p = f"{figd}/stack_activity_top_families.png"; fig.savefig(p, dpi=110); plt.close(fig); figs.append(p)
    if have.get("cont") is not None:
        C = have["cont"]
        fams = C.fid.unique()[:12]
        fig, axs = plt.subplots(len(fams), 1, figsize=(12, 1.6 * len(fams) + 1), squeeze=False, sharex=True,
                                constrained_layout=True)
        for i, f in enumerate(fams):
            ax = axs[i, 0]
            for w, lw, a_ in [(5, 1, 0.5), (30, 2.5, 1.0)]:
                g = C[(C.fid == f) & (C.W == w)].sort_values("start")
                t = pd.to_datetime(g.start) + pd.Timedelta(days=w / 2)
                ax.plot(t, g.A_data, "-", color=C_DATA, lw=lw, alpha=a_, label=f"data W={w}")
                ax.plot(t, g.A_control1, "-", color=C_CTL, lw=lw, alpha=a_ * 0.7, label=f"control W={w}")
            g1 = C[(C.fid == f) & (C.W == 1) & (C.n_tremor_days > 0)]
            for d in pd.to_datetime(g1.start):
                ax.axvspan(d, d + pd.Timedelta(days=1), color=C_MUTED, alpha=0.12, lw=0)
            ax.axhline(0, color=C_MUTED, lw=0.6)
            ax.set_ylabel(UID.get(f), fontsize=7, rotation=0, ha="right")
            if i == 0:
                ax.legend(fontsize=7, frameon=False, ncol=4)
        axs[0, 0].set_title("activity A, continuous 2012-06..2014-12 (grey = PNSN tremor days of the region, label only)", fontsize=9)
        p = f"{figd}/stack_activity_continuous.png"; fig.savefig(p, dpi=110); plt.close(fig); figs.append(p)
    sel = f"{od}/cont_selection.json"
    if os.path.exists(sel):
        s = json.load(open(sel))
        md += ["## Continuous series 2012-06 -> 2014-12", "", f"Selection rule and volume: `{sel}`.",
               f"- families ({len(s['fids'])}): {', '.join(s['family_uids'])}",
               f"- stations ({len(s['stations'])}): {', '.join(s['stations'])}",
               f"- station-days {s.get('n_station_days')} ({s.get('n_to_fetch')} to fetch, ~{s.get('est_GB')} GB)", ""]
    md += ["## Figures", ""] + [f"- `{os.path.relpath(p, W.STUDY)}`" for p in figs] + [""]
    md += ["## Files", "", "- `results/stack_activity/run_<tag>/` daily stacks (npz), family-days, candidates",
           "- `results/stack_activity/score_<tag>/station_window.parquet`, `family_window.parquet|csv`",
           "- `results/stack_activity/validation_by_W.csv`, `validation_by_family.csv`, `trial_by_W.csv`", ""]
    open(f"{od}/STACK_ACTIVITY_RESULTS.md", "w").write("\n".join(md))
    print("\n".join(md[:5]))


if __name__ == "__main__":
    main()
