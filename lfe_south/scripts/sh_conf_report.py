#!/usr/bin/env python
"""SH_CONF report: figures + results/SH_CONF_RESULTS.md from OUT (sh_conf_fit.py outputs).
  python scripts/sh_conf_report.py
"""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_conf_common as C
O = C.OUT; FIG = f"{O}/figures"
HIS, OURS, PCOL, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"


def md(df, fmt=".3f"):
    try:
        return df.to_markdown(index=False, floatfmt=fmt)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def style(ax):
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def fig_monthly(Mo, out, fams):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    Mo = Mo.copy(); Mo["tm"] = pd.to_datetime(Mo.month + "-15")
    n = len(fams); nc = 4 if n > 1 else 1; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(15 if n > 1 else 12, 2.4 * nr + 1.2 if n > 1 else 5.5), sharex=True, squeeze=False)
    for ax, f in zip(axs.flat, fams):
        g = Mo[Mo.fam == f].sort_values("tm")
        ax.plot(g.tm, g.his_rate, color=HIS, lw=1.4, label="Shelly 2025 catalogue (all days)")
        gg = g[g.days > 0]
        if "his_sameday_rate" in gg:
            ax.plot(gg.tm, gg.his_sameday_rate, color=HIS, lw=1.0, ls="--", label="Shelly 2025, on our processed days")
        ax.plot(gg.tm, gg.S1_rate, color=OURS, lw=1.2, label="ours: S1 9×MAD catalogue count")
        ax.plot(gg.tm, gg.sumP_same_rate, color=PCOL, lw=1.6, label="ours: ΣP(real, same family)")
        ax.fill_between(gg.tm, gg.sumP_same_rate - 2 * gg.sumP_same_rate_sd, gg.sumP_same_rate + 2 * gg.sumP_same_rate_sd,
                        color=PCOL, alpha=0.2, lw=0)
        ax.plot(gg.tm, gg.ctl_sumP_same_rate, color=GREY, lw=1.0, ls=":", label="scrambled control ΣP (context floor)")
        ax.axvline(C.T_HIS1, color="0.2", ls="--", lw=1.0)
        ax.axvline(pd.Timestamp("2023-12-01"), color="0.6", ls=":", lw=0.8)
        ax.set_title("all 27 families" if f == "all" else f"family {f}", fontsize=9, loc="left"); ax.tick_params(labelsize=7)
        style(ax)
    for ax in list(axs.flat)[n:]: ax.set_visible(False)
    axs.flat[0].legend(fontsize=7, frameon=False, loc="upper left")
    for ax in axs[:, 0]: ax.set_ylabel("events / day", fontsize=8)
    fig.suptitle("Mendocino LFE families, monthly rate. Before 2023-12 our values are on the v3b test days only; from 2023-12-01 (dotted line)\n"
                 "the continuous forward run. Dashed vertical = his catalogue end (2024-12-01). Band = ±2σ (Bernoulli + fold-model spread).",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97)); fig.savefig(out, dpi=110); plt.close(fig)


def fig_reliability(R, out):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4))
    cols = {"fit_oof": GREY, "ho2024": HIS, "fwd_ovl": OURS}
    for ax, lab in zip(axs, ["same", "any"]):
        ax.plot([0, 1], [0, 1], color="0.5", lw=0.8, ls="--")
        for es, c in cols.items():
            g = R[(R.label == lab) & (R.eval_set == es) & (R.n >= 20)]
            ax.plot(g.mean_P, g.frac, "-o", color=c, ms=4, lw=1.6, label={"fit_oof": "2018-23 fit days (out-of-fold)",
                                                                         "ho2024": "2024 held-out days", "fwd_ovl": "forward run, 2023-12..2024-11"}[es])
        ax.set_xlabel("predicted P"); ax.set_ylabel("matched fraction"); ax.set_title(f"label: {lab} family", loc="left", fontsize=10)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); style(ax)
    axs[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def main():
    os.makedirs(FIG, exist_ok=True)
    A = pd.read_parquet(f"{O}/candidates_conf.parquet", columns=["run", "series", "role", "date", "fam", "t_epoch", "P_same", "P_any",
                                                              "in_S1_9x", "m_same", "m_any", "net_mad"])
    E = pd.read_csv(f"{O}/eval_auc.csv"); K = pd.read_csv(f"{O}/eval_calibration.csv"); R = pd.read_csv(f"{O}/eval_reliability.csv")
    W = pd.read_csv(f"{O}/eval_within_famday.csv"); B = pd.read_csv(f"{O}/eval_bootstrap.csv"); Fm = pd.read_csv(f"{O}/eval_by_family.csv")
    Pm = pd.read_csv(f"{O}/eval_permutation.csv"); PS = pd.read_csv(f"{O}/eval_pnet_station.csv")
    Mo = pd.read_csv(f"{O}/monthly.csv", dtype={"fam": str}); G = pd.read_parquet(f"{O}/famday.parquet")
    MI = json.load(open(f"{O}/model_info.json")); DT = pd.read_csv(f"{O}/day_table.csv", dtype={"date": str})
    fig_monthly(Mo, f"{FIG}/monthly_rate_total.png", ["all"])
    fig_monthly(Mo, f"{FIG}/monthly_rate_per_family.png", C.FAMS)
    fig_reliability(R, f"{FIG}/reliability.png")
    # ---- checks
    metas = [json.load(open(p)) for p in glob.glob(f"{O}/feat/*/*.json")]
    err_mad = max(m["max_err_mad"] for m in metas); err_cs = max(m["max_err_ccsum"] for m in metas)
    n_stored = sum(m["n_cand_stored"] for m in metas); n_data = sum(m["n_data"] for m in metas); n_c1 = sum(m["n_c1"] for m in metas)
    guards = {os.path.basename(p)[6:-4]: open(p).read().strip() for p in sorted(glob.glob(f"{O}/status/GUARD_*.txt"))}
    ft = pd.read_csv(f"{O}/pnet/models/train_ft1.csv") if os.path.exists(f"{O}/pnet/models/train_ft1.csv") else None
    have_f = "PNET_W_ONLY" in MI["features"]
    # ---- headline numbers
    def g(es, lab, col):
        r = E[(E.eval_set == es) & (E.label == lab)]
        return float(r[col].iloc[0]) if len(r) and col in r else np.nan
    def k(es, lab, col, model="FULL"):
        r = K[(K.eval_set == es) & (K.label == lab) & (K.model == model)]
        return float(r[col].iloc[0]) if len(r) else np.nan
    Pb = Pm[Pm.is_block]
    def blk(lab, name):
        r = Pb[(Pb.label == lab) & (Pb.feature == name)]
        return float(r.auc_drop.iloc[0]) if len(r) else np.nan
    fwd = G[G.run == "fwd"]; ext = fwd[fwd.date >= "2024-12-01"]
    jn = fwd[(fwd.date < "2024-12-01")]
    S = C.his_catalogue(); his_full = S[S.t < C.T_HIS1].groupby("fam").size() / (C.T_HIS1 - C.T_HIS0).days
    nd_ext = ext.date.nunique(); nd_jn = jn.date.nunique()
    PF = pd.DataFrame(dict(
        his_rate_2018_2024=his_full,
        join_his_rate=jn.groupby("fam").his_n.sum() / nd_jn,
        join_sumP_same_rate=jn.groupby("fam").sumP_same.sum() / nd_jn,
        join_S1_rate=jn.groupby("fam").n_S1_9x.sum() / nd_jn,
        ext_S1_rate=ext.groupby("fam").n_S1_9x.sum() / nd_ext,
        ext_sumP_same_rate=ext.groupby("fam").sumP_same.sum() / nd_ext,
        ext_sumP_same_sd=np.sqrt((ext.groupby("fam").sdB_same.apply(lambda s: (s ** 2).sum())
                                  + ext.groupby("fam").sdM_same.apply(lambda s: (s ** 2).sum()))) / nd_ext,
        ext_ctl_sumP_same_rate=ext.groupby("fam").ctl_sumP_same.sum() / nd_ext,
        ext_sumP_any_rate=ext.groupby("fam").sumP_any.sum() / nd_ext)).reset_index().rename(columns={"index": "fam"})
    PF["join_obs_over_exp"] = PF.join_his_rate / PF.join_sumP_same_rate
    PF["ext_vs_his_2018_2024"] = PF.ext_sumP_same_rate * PF.join_obs_over_exp / PF.his_rate_2018_2024
    PF = PF.sort_values("fam", key=lambda s: s.astype(int))
    PF.to_csv(f"{O}/per_family_forward.csv", index=False)
    tot_ext_P = ext.sumP_same.sum(); tot_ext_S1 = int(ext.n_S1_9x.sum()); tot_ext_c = int(ext.n_cand.sum())
    roles = DT.groupby(["run", "role"]).size().reset_index(name="days")
    nrows = A[A.series == "data"].groupby(["run", "role"]).size().reset_index(name="data_rows")
    roles = roles.merge(nrows, on=["run", "role"], how="left")
    roles = roles.merge(A[A.series == "c1"].groupby(["run", "role"]).size().reset_index(name="control_rows"), on=["run", "role"], how="left")
    lab_rates = A[(A.series == "data") & A.m_same.notna()].groupby("role")[["m_same", "m_any"]].mean().reset_index()
    up = int((PF.ext_vs_his_2018_2024 > 1.25).sum()); dn = int((PF.ext_vs_his_2018_2024 < 0.8).sum())
    ftline = ""
    if ft is not None:
        kept = ft.loc[ft.val_auc.idxmax()]
        ftline = (f"PNET fine-tune: internal validation AUC {ft.val_auc.iloc[0]:.3f} (SVI weights) → {kept.val_auc:.3f} (kept step {int(kept.step)}). ")
    hs = PS.pivot_table(index=["held_out_of_finetune"], columns="model", values="auc_same", aggfunc="median")
    L = ["# SH_CONF: a per-repeat confidence for Shelly (2025)'s 27 Mendocino families, fitted against his catalogue", "",
         f"(auto-generated by `scripts/sh_conf_report.py`, {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%MZ}; queue `scripts/sh_conf_queue.sh`; "
         f"outputs `results/sh_conf/`)", "",
         "## Plain-language summary", "",
         "Every candidate our detector makes (every ≥ 8×MAD peak, nothing cut) now carries a probability that it is one of "
         "Shelly's catalogued events: **P(same)** = an event of the same family within ±2 s, **P(any)** = an event of any family within "
         "±2 s. The model was learned from his catalogue on 2018–2023 test days and then checked on days it never saw. "
         "The scrambled control is the same detector run with each station site's data shifted to a different time of the same day, so "
         "nothing lines up across sites except by chance.", "",
         f"- **Separation on independent days.** 2024 held-out days: AUC {g('ho2024', 'same', 'auc_FULL'):.2f} (same family) / "
         f"{g('ho2024', 'any', 'auc_FULL'):.2f} (any family); our raw peak height (net_mad) alone gives {g('ho2024', 'same', 'auc_net_mad'):.2f} / "
         f"{g('ho2024', 'any', 'auc_net_mad'):.2f} and the S1-catalogue flag {g('ho2024', 'same', 'auc_in_S1_9x'):.2f} / {g('ho2024', 'any', 'auc_in_S1_9x'):.2f}. "
         f"Forward run, Dec 2023–Nov 2024 (non-build days): {g('fwd_ovl', 'same', 'auc_FULL'):.2f} / {g('fwd_ovl', 'any', 'auc_FULL'):.2f} "
         f"(net_mad {g('fwd_ovl', 'same', 'auc_net_mad'):.2f} / {g('fwd_ovl', 'any', 'auc_net_mad'):.2f}). Average precision "
         f"{g('fwd_ovl', 'same', 'ap_FULL'):.2f} / {g('fwd_ovl', 'any', 'ap_FULL'):.2f} at base rates {g('fwd_ovl', 'same', 'base_rate'):.2f} / "
         f"{g('fwd_ovl', 'any', 'base_rate'):.2f}.",
         f"- **Calibration.** Slope (target 0.8–1.2): 2024 held-out {k('ho2024', 'same', 'slope'):.2f} / {k('ho2024', 'any', 'slope'):.2f}; "
         f"forward overlap {k('fwd_ovl', 'same', 'slope'):.2f} / {k('fwd_ovl', 'any', 'slope'):.2f}. Matched / ΣP: "
         f"{k('ho2024', 'same', 'obs_over_exp'):.2f} / {k('ho2024', 'any', 'obs_over_exp'):.2f} (2024) and "
         f"{k('fwd_ovl', 'same', 'obs_over_exp'):.2f} / {k('fwd_ovl', 'any', 'obs_over_exp'):.2f} (forward overlap).",
         f"- **Scrambled control.** P(same) > 0.5 for {k('fwd_ovl', 'same', 'ctl_fracP_gt05'):.1%} of control candidates vs "
         f"{k('fwd_ovl', 'same', 'data_fracP_gt05'):.1%} of data candidates (forward overlap; 2024 held-out {k('ho2024', 'same', 'ctl_fracP_gt05'):.1%} vs "
         f"{k('ho2024', 'same', 'data_fracP_gt05'):.1%}). The control's total ΣP is {k('fwd_ovl', 'same', 'ctl_over_data_sumP'):.0%} of the data's "
         f"(any family: P > 0.5 {k('fwd_ovl', 'any', 'ctl_fracP_gt05'):.1%} vs {k('fwd_ovl', 'any', 'data_fracP_gt05'):.1%}, ΣP ratio "
         f"{k('fwd_ovl', 'any', 'ctl_over_data_sumP'):.0%}). Control candidates are also much rarer than data candidates "
         f"({n_c1} vs {n_data} in all).",
         f"- **PNET.** {ftline}Permuting all PNET inputs on the independent rows lowers AUC by {blk('same', 'PNET (all)'):.3f} (same) / "
         f"{blk('any', 'PNET (all)'):.3f} (any); a model without PNET reaches {g('fwd_ovl', 'same', 'auc_NO_PNET'):.3f} / "
         f"{g('fwd_ovl', 'any', 'auc_NO_PNET'):.3f} vs {g('fwd_ovl', 'same', 'auc_FULL'):.3f} / {g('fwd_ovl', 'any', 'auc_FULL'):.3f} with it (forward overlap)"
         + (f"; off-the-shelf PNET only: {g('fwd_ovl', 'same', 'auc_PNET_W_ONLY'):.3f} / {g('fwd_ovl', 'any', 'auc_PNET_W_ONLY'):.3f}." if have_f else "."),
         f"- **Beyond his catalogue** (2024-12-01 .. {ext.date.max()}, {nd_ext} days): {tot_ext_c} candidates, ΣP(same) = {tot_ext_P:.0f} expected "
         f"same-family repeats vs {tot_ext_S1} S1 catalogue detections. Scaled by each family's his/ΣP ratio at the join year, {up} families run more "
         f"than 25 % above their 2018–2024 average and {dn} more than 20 % below it (table in section 5).",
         "- **Leakage handled.** The templates were stacked from his events on the v3b build days, so every build day of any family is "
         "flagged in-sample and kept out of fitting and of every evaluation (P is still written for those rows, flagged). Fit = 132 "
         "strict 2018–2023 test days (out-of-fold by day); the 2024 held-out days and the forward overlap are scored by a model that never "
         "saw them. The PNET fine-tune used only build days (and never the 3 held-out stations).",
         "- Caveats: P is the probability of matching HIS catalogue, a lower bound on 'real' (his catalogue is incomplete). Forward use assumes "
         "the 2018–2023 relation between features and his catalogue carries over (station roster, KMR re-templates, base rates); the "
         "join-year matched/ΣP above is the check, and per-family rates are scaled by it.", "",
         "![monthly total](sh_conf/figures/monthly_rate_total.png)", "",
         "![monthly per family](sh_conf/figures/monthly_rate_per_family.png)", "",
         "![reliability](sh_conf/figures/reliability.png)", "",
         "## 1. Design", "",
         "- Candidates: every ≥ 8×MAD 4-s-bin maximum of the v3b it0 run on the 385 Shelly-2025 test days 2018–2024 and of the forward run "
         "2023-12-01..2026-09-28 (all segments: cont/A it0, NK = no KMR on the two KMR switch days, B/C = KMR re-templated). Scrambled control "
         "'c1': the same detector (t11_mf_v2 functions, same prep, same templates) with each ≤ 3-km site's per-channel cc series rolled by its "
         "own random shift of the same day (10 min..1 day−10 min, ≥ 2 min apart; NC.KBU + PB.B049 are one site), ≥ 8×MAD, 4-s bins.",
         f"- Check: the re-run reproduces every stored candidate: max |net_mad error| {err_mad:.2e}, max |ccsum error| {err_cs:.2e} "
         f"({n_data} of {n_stored} stored candidates found).",
         "- Features (no truth-based feature): MF net_mad, mean cc, summed cc, live channels / fraction, nchan, nsta; per ≤ 3-km site "
         "(live-channel-weighted mean of the station cc at the aligned position): site_cc_max/2nd, site_z_max/2nd/min (z vs that site's cc at "
         "1000 random origins of the same family-day), n sites z ≥ 2/3/5, frac z ≥ 3, top-site share of the sum, median site noise; "
         "hfr_rel = log10 E(10–19 Hz)/E(2–8 Hz) over the 12-s window at the family's HF station (highest pick_q with raw data) minus the median of "
         "300 random windows of the same station-day; candidate density of the family and of its ≤ 3-km group within ±1 h and ±10 min (same "
         "day), other-family candidates within ±4 s; S1 flags (9×MAD, mean cc ≥ 0.12, summed cc ≥ 2.75, 4-s cross-family winner, in the S1 "
         "catalogue); PNET (off-the-shelf SVI `pnet_w05`, suffix _w; Mendocino fine-tune, suffix _f): max / 2nd site value, sites > 0.3 / 0.5, "
         "max site value above its random-origin median, at each station's expected S (template pick, ±0.5 s); family id (categorical). "
         "Per-station cc, z and PNET values are kept as columns.",
         "- Labels: nearest of his events within ±2 s, same family (`m_same`) and any family (`m_any`); one model per label.",
         f"- Model: {MI['model']}. Feature sets: " + "; ".join(f"{k_} ({len(v)})" for k_, v in MI["features"].items()) + ".",
         "- PNET input: the cache is the model's 20 Hz, 1–8 Hz chain; horizontals rotated to N/E with the FDSN azimuths of the day "
         "(PB EH1/EH2, BK.BRIC 7°/97° before 2024-05-21); NC.KPR / NC.KBU are vertical-only → (Z, 0, 0). Fine-tune (WA recipe): SVI start, "
         f"windows on v3b BUILD days only at 14 stations (held out {', '.join(C.PNET_HO)}), positives = his events at stations whose template "
         "pick_q ≥ 3 (weight 1 if ≥ 10 else 0.5), negatives = random windows with ±3 s masks around every catalogued S, best internal-AUC "
         "checkpoint (every 10th build date).", "",
         "### Day roles and rows", "", md(roles, ".0f"), "",
         "Matched fractions by role (data rows):", "", md(lab_rates), "",
         "### Data-read guard (stage fails if > 3 % of station-days error)", "",
         "\n".join(f"- {k_}: {v}" for k_, v in guards.items()), "",
         "## 2. Separation (AUC / average precision)", ""]
    cols = ["eval_set", "label", "n", "n_pos", "base_rate", "n_days", "auc_FULL", "ap_FULL", "auc_NO_PNET"] + \
           (["auc_PNET_W_ONLY"] if have_f else []) + ["auc_NO_FAM", "auc_MF_S1", "auc_FULL_ownbuild", "auc_net_mad", "auc_meancc",
                                                       "auc_site_z_max", "auc_pnet_max_w", "auc_pnet_max_f", "auc_in_S1_9x", "ap_net_mad"]
    L += [md(E[[c for c in cols if c in E]]), "",
          "`FULL_ownbuild` = sensitivity: fitted on 2018–2023 test days excluding only the family's OWN build days (other families' build days "
          "allowed); evaluated on the same independent rows.", "",
          "Within one family-day (does P rank the repeats of one family on one day?), weighted mean AUC:", "", md(W), "",
          "Day-block bootstrap (200 resamples of days) of AUC(FULL) − AUC(other):", "", md(B), "",
          "## 3. Calibration and the scrambled control", "",
          md(K[["eval_set", "label", "model", "n", "slope", "icpt", "matched", "sumP", "obs_over_exp", "brier", "data_fracP_gt05",
                "ctl_fracP_gt05", "data_meanP", "ctl_meanP", "n_ctl", "ctl_over_data_sumP"]]), "",
          "## 4. PNET contribution", "",
          "Permutation importance of the final FULL model on the independent rows (2024 held-out + forward overlap); blocks permuted jointly:", "",
          md(Pb[["label", "feature", "auc_drop", "logloss_gain", "n"]]), "",
          "Top single features:", "",
          md(Pm[~Pm.is_block].sort_values(["label", "auc_drop"], ascending=[True, False]).groupby("label").head(12)[["label", "feature", "auc_drop", "logloss_gain"]]), "",
          "Single-station PNET value (matched same family vs unmatched), independent rows; median over stations by fine-tune status:", "",
          md(hs.reset_index()), "", md(PS), ""]
    if ft is not None:
        L += ["Fine-tune internal validation (build days only):", "", md(ft), ""]
    L += ["## 5. Forward run: expected repeats per family", "",
          "join = forward days 2023-12..2024-11 (all, incl. build days); ext = 2024-12-01 onward. join_obs_over_exp = his events / ΣP(same) at the "
          "join; ext_vs_his_2018_2024 = ext ΣP rate × that ratio / his 2018–2024 rate. sd = Bernoulli + 5 fold-model spread.", "",
          md(PF), "",
          "### By family on the independent rows (label same; families with ≥ 10 positives and ≥ 10 negatives)", "",
          md(Fm[(Fm.label == "same") & (Fm.eval_set != "fit_oof")].drop(columns=["label"])), "",
          "## 6. Files (results/sh_conf/)", "",
          "- `candidates_conf.parquet`: every candidate row (run v3b / fwd, series data / c1, role, all features incl. per-station cc_/z_/pnw_/pnf_ "
          "columns, labels, P_same / P_any (= FULL) and every ablation P_<label>_<set>).",
          "- `forward_P_by_det_id.parquet` (copy `results/cut/sh_ext/conf_P_by_det_id.parquet`): one row per cut-window det_id (join to "
          "`results/cut/sh_ext/index_*.parquet` and `detections.parquet` on det_id) with P_same, P_any, the no-PNET versions, role and in_sample.",
          "- `famday.parquet`: run × day × family: n_cand, n_S1_9x, ΣP same/any with Bernoulli sd and (forward) fold-model sd, scrambled-control "
          "ΣP and count, his count, matched counts. `monthly.csv`, `per_family_forward.csv`.",
          "- `feat/<run>/<date>.{parquet,json}` (+ `_noise.parquet`: per family-day-station noise median/MAD), `hf/`, `pnet/feat_*`, "
          "`pnet/models/ft1.pt`, `eval_*.csv`, `models.joblib`, `model_info.json`, `day_table.csv`, `status/`.", "",
          "## 7. Deviations", "",
          "- One scrambled control series (c1) per day, as specified; WA used four.",
          "- Density features are counted within the same calendar day for both runs (the v3b test days are isolated days, the forward days are "
          "contiguous), so they mean the same thing in fit and application.",
          "- 'wsum' (WA: sum of template weights) is represented by the live channel count and live fraction at the peak (all Shelly channels carry "
          "equal weight in his sum).",
          "- Family groups use the families_south.csv locations, ≤ 3 km epicentral (as WA); most families have 14–18 others within 3 km.",
          "- HF stations come from the raw miniSEED day files (traces/), as WA; a family whose ranked stations lack raw data that day gets NaN.", ""]
    import fcntl
    with open(f"{C.WS}/locks/sh_results_md.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX); open(os.environ.get("SH_CONF_REPORT", f"{C.RES}/SH_CONF_RESULTS.md"), "w").write("\n".join(L))
    C.log("[report] written")


if __name__ == "__main__":
    main()
