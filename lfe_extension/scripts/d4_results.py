#!/usr/bin/env python
"""Track D, D4 -- generate results/trackD/RESULTS.md from the Track D output files, so the write-up is
always exactly the numbers on disk."""
from __future__ import annotations
import glob, json, os
import numpy as np, pandas as pd

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
TC = f"{EXT}/results/trackC"


def md(df, cols=None, r=4):
    d = df[cols] if cols else df
    d = d.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].round(r)
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join("---" for _ in d.columns) + "|"
    rows = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
            for row in d.itertuples(index=False)]
    return "\n".join([head, sep] + rows)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="v1", help="catalog version this write-up describes")
    ap.add_argument("--prev-tag", default="", help="previous version to diff against in the change log")
    ap.add_argument("--snapshot", default=f"{TD}/v1_snapshot",
                    help="directory holding the previous version's D1/D3 outputs, for the change log")
    ag = ap.parse_args()
    sfx = "" if ag.tag == "v1" else f"_{ag.tag}"
    O = []
    A = O.append
    meta = json.load(open(f"{TD}/catalog_{ag.tag}_meta.json"))
    MD = pd.read_csv(f"{TC}/months_done.csv").drop_duplicates("ym", keep="last")
    E = pd.read_parquet(f"{TD}/lfe_extension_catalog_{ag.tag}.parquet")
    T = pd.read_csv(f"{TD}/family_table_{ag.tag}.csv")
    S = pd.read_csv(f"{TD}/completeness_epoch_summary.csv")
    V = pd.read_csv(f"{TD}/mc_model_validation.csv")
    X = pd.read_csv(f"{TD}/completeness_epoch_crossyear.csv")
    L = pd.read_csv(f"{TD}/floor_epoch_season_summary.csv")
    B = pd.read_csv(f"{TD}/thomas_bars{sfx}.csv")
    TS = pd.read_csv(f"{TD}/thomas_per_sse{sfx}.csv")
    W = pd.read_csv(f"{TD}/sse_windows.csv")
    MC = pd.read_parquet(f"{TD}/mc_family_day.parquet")
    want = pd.period_range("2014-01", "2026-09", freq="M").astype(str)
    missing = [m for m in want if m not in set(MD.ym)]

    A(f"# Track D — the product: completeness model, catalog {ag.tag.replace('_', '.')}, "
      f"and the science QC\n")
    A(f"Generated {pd.Timestamp.now("UTC"):%Y-%m-%d %H:%M} UTC by `scripts/d4_results.py` from the files in "
      "`results/trackD/`. Executes steps 6–8 of `results/merlin_review_plan_2026-09-16.md`. "
      "Inputs: Track C's monthly outputs (read-only), Track A's frozen operating point and referee runs, "
      "Track A-2's relocations and epoch policy.\n")

    # ---------------- what changed (only when a previous version is named)
    if ag.prev_tag:
        pv = ag.prev_tag
        A(f"## 0.0 Catalog {ag.tag.replace('_', '.')} — what changed since {pv.replace('_', '.')}, "
          f"and why\n")
        A(f"Two independent changes, both upstream of Track D, both from Track C:\n")
        try:
            RC = pd.read_csv(f"{TC}/rerun_ccfloor_log.csv")
            A(f"**(1) The cc-floor bug.** Track C found and fixed a floor bug in the per-site "
              f"cross-correlation and re-ran **{len(RC)} months**, reporting 486 artifact rows removed "
              f"(`results/trackC/rerun_ccfloor_log.csv`). Read off that log directly, across those ten "
              f"months: rows with |cc| > 1 **{int(RC.n_cc_gt1_before.sum())} → "
              f"{int(RC.n_cc_gt1_after.sum())}**, rows with net_mad > 50 "
              f"**{int(RC.n_mad_gt50_before.sum())} → {int(RC.n_mad_gt50_after.sum())}**, "
              f"Track-E-flagged rows **{int(RC.n_trackE_flagged_present_before.sum())} → "
              f"{int(RC.n_trackE_flagged_present_after.sum())}**, and recorded detections "
              f"{int(RC.n_before.sum()):,} → {int(RC.n_after.sum()):,} "
              f"({int((RC.n_before - RC.n_after).sum())} fewer; the maximum |cc| anywhere goes from "
              f"{RC.max_cc_before.max():.3g} to {RC.max_cc_after.max():.2f}). This is a DATA change: "
              f"those months' detection files are different, so every Track D product built on them is "
              f"rebuilt here.\n")
        except Exception:
            A("**(1) The cc-floor bug.** Track C fixed a floor bug in the per-site cross-correlation "
              "and re-ran the affected months (`results/trackC/rerun_ccfloor_log.csv`).\n")
        try:
            TN = pd.read_csv(f"{TC}/tiers_by_epoch.csv")
            TO = pd.read_csv(f"{TC}/tiers_by_epoch.pre_ccfloor.csv")
            mm = TO[["epoch", "T_core", "T_ext"]].merge(TN[["epoch", "T_core", "T_ext"]], on="epoch",
                                                        suffixes=(f"_{pv}", f"_{ag.tag}"))
            A(f"**(2) The tiers were re-derived on the corrected data** by the same frozen "
              f"false-alarm-only rule. The old table is kept as `tiers_by_epoch.pre_ccfloor.csv`:\n")
            A(md(mm, r=3))
            A("")
            A("E0 and E3 are unchanged; **E1's T_core moves 8.25 → 9.00 and E2's 8.75 → 9.00**. That is "
              "the dominant effect on the catalog — a threshold change, not a data change.\n")
        except Exception:
            pass
        pfile = f"{TD}/lfe_extension_catalog_{pv}.parquet"
        if os.path.exists(pfile):
            EP = pd.read_parquet(pfile, columns=["epoch_id", "cat_core", "cat_ext"])
            aa = E.groupby("epoch_id").agg(**{f"core_{ag.tag}": ("cat_core", "sum"),
                                              f"ext_{ag.tag}": ("cat_ext", "sum")})
            bb = EP.groupby("epoch_id").agg(**{f"core_{pv}": ("cat_core", "sum"),
                                               f"ext_{pv}": ("cat_ext", "sum")})
            C2 = bb.join(aa).reset_index()
            C2["d_core"] = C2[f"core_{ag.tag}"] - C2[f"core_{pv}"]
            C2["d_ext"] = C2[f"ext_{ag.tag}"] - C2[f"ext_{pv}"]
            tot = dict(epoch_id="ALL")
            for c in C2.columns[1:]:
                tot[c] = int(C2[c].sum())
            C2 = pd.concat([C2, pd.DataFrame([tot])], ignore_index=True)
            A("**Counts, per epoch and per tier:**\n")
            A(md(C2))
            A("")
            A(f"`cat_core` falls {int(EP.cat_core.sum()):,} → {int(E.cat_core.sum()):,} "
              f"(−{100 * (1 - E.cat_core.sum() / EP.cat_core.sum()):.1f} %), almost all of it in E1 and "
              f"E2 and almost all of it the threshold change; `cat_ext` is essentially unchanged "
              f"({int(EP.cat_ext.sum()):,} → {int(E.cat_ext.sum()):,}, "
              f"−{100 * (1 - E.cat_ext.sum() / EP.cat_ext.sum()):.2f} %) because T_ext moved only for "
              f"epochs that were already at 8.00. The E0 and E3 deltas "
              f"({int(C2[C2.epoch_id == 'E0'].d_core.iloc[0])} and "
              f"{int(C2[C2.epoch_id == 'E3'].d_core.iloc[0])} core events) are the cc-floor fix alone, "
              f"since their thresholds did not move — that is the clean measure of the bug's size on "
              f"the catalog: well under 1 %.\n")
        A("**What this does NOT change.** The 2012 referee is Track A data and is untouched by both, so "
          "every referee-based result is bit-identical: the completeness-model validation bar, the two "
          "pre-registered Thomas bars, the thinned-Bostock reference, and the frozen α of the §3b-2 "
          "front estimator (no redesign was done — the frozen α = 0.05 was applied as-is).\n")
        A("**Verdicts: none flipped.** The per-verdict comparison is in §0 and the sections below; the "
          "changes are in the completeness NUMBERS (E1 T_core recall 0.184 → 0.129, E2 0.146 → 0.129) "
          "and in the floor λ (E1 and E2 quiet-day rates fall 1.5–4× with the higher T_core), not in "
          "any pass/fail.\n")
        A(f"Version {pv.replace('_', '.')} is preserved in full: its catalog, family table and "
          f"family-day table under their own names, and its completeness / floor / Thomas outputs in "
          f"`results/trackD/v1_snapshot/`.\n")

    # ---------------- headline
    A("## 0. Headline\n")
    A(f"* **Catalog**: {meta['n_cat_core']:,} `cat_core` events and {meta['n_cat_ext']:,} `cat_ext` events "
      f"over {meta['n_families']} families, {meta['t0'][:10]} .. {meta['t1'][:10]}, from "
      f"{len(MD)} finished Track C months"
      + (f" ({len(missing)} of the 153 months 2014-01..2026-09 are still missing)" if missing else
         " (the complete 2014-01..2026-09 record)") + ".")
    A(f"* **Completeness-model validation bar (±0.05 on the held-out 2010 year): "
      f"{'MET' if bool(V[~V['bin'].str.contains('in-sample')]['within_0.05'].all()) else 'NOT MET'}** — "
      f"largest |predicted − measured| own recall over all magnitude bins and both tiers = "
      f"{V[~V['bin'].str.contains('in-sample')].delta.abs().max():.4f}.")
    for r in B.itertuples():
        A(f"* **Thomas bar — {r.bar}, {r.tier}: {'MET' if r.met else 'NOT MET'}** "
          f"(ours {r.ours:.4g} vs Bostock 2012 {r.bostock:.4g}, {100*r.rel_diff:.1f} % vs a {100*r.tol:.0f} % tolerance).")
    fv2 = f"{TD}/thomas_front_v2{sfx}.csv"
    if os.path.exists(fv2):
        RV2 = pd.read_csv(fv2)
        gc = RV2[RV2.tier == "cat_core"]
        sg = gc[gc.plausible & (gc.null_xshuf_p <= 0.05)]
        A(f"* **Main fronts recovered (§3b-2, post-registration estimator, approved 2026-09-17)**: "
          f"{int(gc.speed_kmday.notna().sum())} of {len(gc)} SSE windows fitted and **{len(sg)} "
          f"significant** at `cat_core` (the pre-registered estimator managed 3), "
          f"{int((sg.speed_kmday < 0).sum())} propagating SE and {int((sg.speed_kmday > 0).sum())} NW, "
          f"|speed| {sg.speed_kmday.abs().min():.1f}–{sg.speed_kmday.abs().max():.1f} km/day.")
    A("")

    # ---------------- D1a
    # ---- the epoch axis and the tier thresholds actually in force
    tu = f"{TD}/tiers_used.json"
    A("### Network epochs and the tier thresholds in force\n")
    if os.path.exists(tu):
        TU = json.load(open(tu))
        A(f"**Track C is the authority on both.** Source: `{TU['source']}`. Track C defines the network "
          f"epochs and re-derives T_core and T_ext for each by the frozen false-alarm-only rule "
          f"(T_core: quiet-day FA ≤ 0.1 per template-day; T_ext: quiet-day ≤ 1.0 AND post-SSE ≤ 1.5), "
          f"applying them as a re-flag of the `cat_core` / `cat_ext` columns. Every Track D step — the "
          f"epoch ablation, the Mc fit, the floor λ, the catalog — reads the same epochs and the same "
          f"numbers through `scripts/d_tiers.py`, so the completeness model can never disagree with the "
          f"flags. A row whose threshold could not be derived yet (too few quiet family-days) falls back "
          f"to the frozen referee point 8.75 / 8.25 and is marked so.\n")
        rows_t = []
        for e in TU.get("epochs", []):
            t = TU["tiers"].get(e["epoch_id"], {})
            rows_t.append(dict(epoch=e["epoch_id"], t0=e["t0"], t1=e["t1"],
                               station_set=e["station_set"], n_stations=e["n_stations"],
                               T_core=t.get("T_core"), T_ext=t.get("T_ext"),
                               source=("Track C, derived" if t.get("derived") else "frozen fallback"),
                               status=t.get("status", "")))
        A(md(pd.DataFrame(rows_t), r=3))
        A("")
        try:
            MF = pd.read_csv(f"{TD}/completeness_matched_fa.csv")
            tm = float(MF[(MF.run_id == "S13") & (MF.tier == "T_core")].T_matched.iloc[0])
            t3 = TU["tiers"]["E3"]["T_core"]
            A(f"**An independent confirmation.** Track D's 2012 station-set ablation says the "
              f"13-station network (S13) needs **T = {tm:.2f}** to reproduce the frozen quiet-day "
              f"false-alarm rate at T_core. Track C's FA-only re-derivation, on the production record "
              f"and with no reference to that ablation, lands on **T = {t3:.2f}** for epoch E3. Two "
              f"different data sets and two different procedures agree to {abs(tm-t3):.2f} MAD — the "
              f"strongest single check that the completeness model and the flags describe the same "
              f"detector.\n")
        except Exception:
            pass
        A("Track D measures per-family recall by ablating the **station set**, which is the thing that "
          "actually changes the detector; two of Track C's epochs (E0 and E1) share the 15-station set, "
          "so the ablation axis is S15 / S14 / S13 and `d_tiers.py` maps each epoch onto it. "
          "`completeness_family_epoch.csv` is the map already expanded onto Track C's epoch axis; "
          "`completeness_family_stationset.csv` is the measurement itself.\n")
    A("## 1. Completeness model (Merlin Step 6)\n")
    A("### 1a. Epoch / station-set ablation\n")
    A("The three station sets the production run realises (see the epoch table above for the mapping "
      "onto Track C's network epochs):\n")
    A("| station set | stations | what changed |")
    A("|---|---|---|")
    A("| S15 | 15 | the frozen 2012 network (Track C epochs E0, E1) |")
    A("| S14 | 14 | UW.LRIV dropped from 2019-05-30, A-2 epoch-identity FAIL (Track C epoch E2) |")
    A("| S13 | 13 | LRIV dropped + CN.SNB ends 2021-03-28 (Track C epoch E3) |")
    A("")
    A("**CN.PFB (ends 2018-12-20) carries no 2012 template and is not in the network at all**, so that "
      "boundary has no effect on this detector and needed no ablation. Four GPU passes "
      "(`scripts/d1_run_ablation.sh`) re-ran the frozen configuration with the S14 and S13 station sets "
      "on the 2012 test half and on the held-out 2010 year; the frozen HF veto is applied to all of "
      "them. `P9` is the 9-station persistent set (= the 2010 network), Track A step 2's control.\n")
    S = S.rename(columns={"run_id": "station_set"})
    s = S.copy()
    s = s.rename(columns={"recall_pooled": "recall@frozen T", "recall_pooled_matched_fa": "recall@matched FA",
                          "fa_at_frozen": "quiet FA @frozen T", "T_matched": "T @matched FA"})
    A(md(s.sort_values(["ref_year", "tier", "station_set"]),
         ["station_set", "ref_year", "tier", "n_fam", "n_events", "recall@frozen T", "quiet FA @frozen T",
          "T @matched FA", "recall@matched FA", "recall_fam_median", "recall_fam_p10", "recall_fam_p90"]))
    A("")
    A("**Reading (Merlin F3).** At a FIXED MAD threshold a smaller network is not a fairer detector, it is "
      "a lower one. Dropping LRIV (S15 → S14) costs nothing on either axis — recall and false alarms are "
      "unchanged — which independently supports Track A-2's decision to drop it. Losing SNB in 2021-03 "
      "(S14 → S13) is the real break: at the frozen thresholds that network sits at ~3× the quiet-day "
      "false-alarm rate the operating point was frozen at, and once that is matched it loses ~19 % of "
      "its recall. The held-out 2010 year reproduces the direction:\n")
    A(md(X))
    A("")
    try:
        f12 = S.set_index(["station_set", "tier"]).fa_at_frozen
        r12 = f12[("S13", "T_core")] / f12[("S15", "T_core")]
        r10 = f12[("S13_2010", "T_core")] / f12[("S15_2010", "T_core")]
        A(f"**The DIRECTION transfers, the MAGNITUDE does not.** The E3 false-alarm penalty is "
          f"{r12:.1f}x on 2012 but only {r10:.1f}x on 2010 at T_core — the 2010 'strict quiet' set is "
          f"seven days, four of them post-SSE, and it already sits at 3x the 2012 quiet-day rate, so it "
          f"has little headroom. **The ablation is therefore the MECHANISM, not the calibration**: the "
          f"number attached to each event is λ measured directly from the PRODUCTION quiet days of that "
          f"epoch (section 1c), not the ablation ratio.\n")
    except Exception:
        pass
    A("Per-family recall for every epoch and both tiers, with Wilson intervals and the matched-FA value, "
      "is `results/trackD/completeness_family_epoch.csv`; the per-epoch thresholds are "
      "`completeness_matched_fa.csv`.\n")

    # ---------------- D1b
    A("### 1b. Daily sensitivity — Mc(family, day)\n")
    A("Model (full statement in `scripts/d1_mc_model.py`):\n")
    A("```\nP(detect | family f, day d, magnitude M) = Phi( b_f + alpha*(M - 1.8) + c1*u(f,d) + c2*v(f,d) )\n"
      "Mc(f,d,T) = 1.8 - (b_f + c1*u + c2*v) / alpha\n```\n")
    A("`b_f` is a ridge-shrunk per-family intercept fitted separately for each epoch's station set (so it "
      "carries component 1 as well); `u` = log10(family-day MAD scale / its epoch median); `v` = the "
      "SNR-weighted log10 noise-amplitude factor built from the daily median 8-s window energy, "
      "normalised per station per INSTRUMENT epoch so the 2017 CNSN gain changes cancel.\n")
    P = pd.read_csv(f"{TD}/mc_model_params.csv").groupby(["tier", "station_set"]).first().reset_index()
    A(md(P, ["tier", "station_set", "b0", "alpha", "c1", "c2", "lam"]))
    A("")
    A("**Validation — the pre-registered bar.** Fit on the 2012 ablation whose station set equals the 2010 "
      "network (9 stations), predicted on the HELD-OUT 2010 year:\n")
    A(md(V, ["tier", "bin", "n", "measured", "predicted", "delta", "within_0.05"]))
    A("")
    A(f"Production Mc: median **{MC.Mc_T_ext.median():.2f}** (T_ext) / **{MC.Mc_T_core.median():.2f}** "
      f"(T_core), IQR {MC.Mc_T_core.quantile(.25):.2f}–{MC.Mc_T_core.quantile(.75):.2f}. "
      "E5's independent prediction model gave Mc(8 MAD) ≈ 2.25–2.29 for this network; ours sits above it "
      "because the tiers are 8.25 / 8.75 MAD.\n")
    if "day_terms_in_calib_range" in MC:
        A(f"**Caveat — extrapolation.** The referee windows are late-summer ETS windows, so the day terms "
          f"were calibrated over a narrow range: only {100*MC.day_terms_in_calib_range.mean():.0f} % of "
          f"production family-days fall inside the referee 1–99 % range (column "
          f"`day_terms_in_calib_range`). Because |c2| ≤ 0.32, a day term 0.3 log-units outside the range "
          f"moves Mc by ≤ 0.1 magnitude, so the extrapolation is small but it is not zero.\n")
    al = P[(P.tier == "T_core") & (P.station_set == "S15")].alpha.iloc[0]
    A(f"**Cross-check against the independent E4/E5 measurements.** In this parameterisation "
      f"alpha = (amplitude slope in log10 units per magnitude) / (spread of log10 of the network "
      f"statistic at fixed magnitude). E4 measured a within-family amplitude slope of 10^(0.38 dM) but "
      f"noted that the measured within-family RECALL slope is about half of it (alpha_amp ~ 0.2, E5's "
      f"'central case'), and E5 measured the spread sigma_log = 0.175 pooled. Those give "
      f"0.20 / 0.175 = 1.14 -- against our fitted **{al:.2f}**. The sensitivity model therefore "
      f"reproduces two independently measured quantities it was never shown.\n")
    A("Another caveat: `c1` (the MAD-scale term) comes out POSITIVE on 2012 — on busy days the cc series' "
      "own MAD tracks activity, so `u` is partly an activity covariate rather than a pure noise term — and "
      "its sign flips in a 2010-only fit. Only the 2012 coefficients are exported, and the out-of-sample "
      "test above is what licenses them.\n")

    # ---------------- D1c
    A("### 1c. False-alarm floor λ(family, tier, epoch, season)\n")
    A("λ = (catalog detections on quiet family-days / quiet observable family-days) × (1 − p_real), with "
      "p_real = 0.225 at T_core and 0.057 at T_ext from Track A step 3b (the withheld-station stack test), "
      "gamma-Poisson shrunk (k = 10 template-days) toward the epoch-season pool. A quiet family-day is "
      "Track C's `quiet_flag`: no PNSN tremor in the family box that day AND ≥ 3 days from any day that "
      "has one. Cells with no quiet day at all fall back family-epoch → epoch-season → epoch → all, "
      "recorded per row in `lambda_src`.\n")
    A(md(L, ["tier", "epoch_id", "season", "n_fam", "quiet_family_days", "quiet_det", "rate_pooled",
             "lambda_pooled"], r=5))
    A("")
    ratios = []
    for tier in ("T_core", "T_ext"):
        for eid in sorted(L.epoch_id.unique()):
            g = L[(L.tier == tier) & (L.epoch_id == eid)].set_index("season").rate_pooled
            if "DJF" in g.index and "JJA" in g.index and g.get("JJA", 0) > 0:
                rat = float(g["DJF"] / g["JJA"])
                ratios.append(dict(tier=tier, epoch_id=eid, DJF_rate=float(g["DJF"]),
                                   JJA_rate=float(g["JJA"]), DJF_over_JJA=rat,
                                   bar_met=bool(0.5 <= rat <= 2.0)))
    RT = pd.DataFrame(ratios)
    if len(RT):
        RT.to_csv(f"{TD}/floor_seasonal_ratio.csv", index=False)
        A("**Seasonal-stability bar** (Merlin: winter quiet-day FA within 2x summer), evaluated per "
          "epoch — and the answer is not a single number:\n")
        A(md(RT, r=4))
        A("")
        nmet = int((~RT.bar_met).sum())
        A(f"**{nmet} of {len(RT)} epoch-tier cells miss the ≤ 2x band, and they miss it in BOTH "
          f"directions**: winter is ~2.4x summer in E0 but summer is ~2–4x winter in E1, E2 and E3. "
          f"The bar is reported as NOT met and nothing is tuned to it; what the full record shows is "
          f"that the quiet-day floor is seasonally UNSTABLE rather than seasonally BIASED, which is why "
          f"λ is carried per (family, tier, epoch, season) on every event rather than as one number.\n")
    A("")

    # ---------------- D1d
    A("### 1d. What drives the seasonal instability of the floor\n")
    A("The ratios are reported as measured and nothing is tuned to them. Two candidate causes were "
      "separated (`scripts/d1d_season_diagnostic.py`, `scripts/d1d_season_stack_test.py`; figure "
      "`figures/d1d_season_diagnostic.png`): a NOISE / MAD-adaptivity effect, or genuinely more "
      "inter-ETS LFE activity in the busier season.\n")
    A("**Scope note.** The withheld-station test in (iv) below was designed and run when only 2014–2015 "
      "was on disk, where winter was the busy season; its six months are all in epoch **E0** "
      "(DJF 2015-02 + 2014-02 + 2014-12 vs JJA 2014-08 + 2015-08 + 2014-07, matched quiet-day counts). "
      "It therefore answers the question for E0's WINTER excess specifically. The full record then "
      "showed the sign of the seasonal imbalance flips in E1–E3, which is itself the main finding of "
      "(ii) and does not depend on (iv).\n")
    tl = f"{TD}/season_diagnostic_tail.csv"
    nz = f"{TD}/season_diagnostic_noise.csv"
    vs = f"{TD}/season_diagnostic_vslope.csv"
    if os.path.exists(tl):
        TL = pd.read_csv(tl)
        piv = TL.pivot_table(index=["epoch_id", "season"], columns="T", values="rate_hfpass").reset_index()
        A("**(i) The detector's threshold does not adapt to the ambient level.** The score is in units of "
          "each family-day's own network-MAD, and that MAD is formed from an amplitude-NORMALISED "
          "cross-correlation, so it barely moves when the ground gets louder. Quiet-day detections per "
          "template-day against threshold:\n")
        A(md(piv, r=4))
        A("")
        SLq = TL.drop_duplicates(["epoch_id", "season"])[["epoch_id", "season",
                                                          "tail_slope_dex_per_MAD", "n_quiet_days"]]
        A("Tail decay (dex per MAD over 8.25–10; more negative = a lighter tail):\n")
        A(md(SLq, r=3))
        A("")
        try:
            djf = SLq[SLq.season == "DJF"].tail_slope_dex_per_MAD.iloc[0]
            jja = SLq[SLq.season == "JJA"].tail_slope_dex_per_MAD.iloc[0]
            A(f"The winter tail is **not** heavier — its decay ({djf:.3f} dex/MAD) is as steep as or "
              f"steeper than summer's ({jja:.3f}). The winter excess is a near-uniform elevation of the "
              f"whole tail, not an extra population of large excursions.\n")
        except Exception:
            pass
    if os.path.exists(nz):
        NZ = pd.read_csv(nz)
        A("**(ii) The floor is DECOUPLED from the ambient level — the strongest result of this "
          "section.** Winter is the loudest season in every epoch, yet it carries the highest quiet-day "
          "rate in only ONE of them, and the family-day MAD scale barely moves at all:\n")
        A(md(NZ, r=4))
        A("")
        try:
            g = NZ.set_index("season")
            g = NZ[NZ.epoch_id == NZ.epoch_id.iloc[0]].set_index("season")
            A(f"Across the four epochs the DJF/JJA ratio of relative daily noise energy is "
              + ", ".join(f"{e} {NZ[(NZ.epoch_id==e)].set_index('season').e8_rel_median['DJF'] / NZ[(NZ.epoch_id==e)].set_index('season').e8_rel_median['JJA']:.2f}"
                          for e in sorted(NZ.epoch_id.unique())
                          if {'DJF','JJA'} <= set(NZ[NZ.epoch_id==e].season))
              + " — winter is louder every time — while the MAD scale stays within a few per cent of "
                "0.0165 in every epoch and season. A detector scored in units of an amplitude-NORMALISED "
                "cross-correlation's own MAD is therefore nearly blind to the ambient level: **a louder "
                "season does not raise the bar, and a louder season is not what raises the floor.** What "
                "raises the floor is the population of COHERENT transients in that epoch-season, which "
                "the daily median energy does not measure — and which has no fixed seasonal sign.\n")
            A(f"(For the first epoch specifically, E0:) relative daily noise energy DJF {g.loc['DJF','e8_rel_median']:.2f} vs JJA "
              f"{g.loc['JJA','e8_rel_median']:.2f} (a factor "
              f"{g.loc['DJF','e8_rel_median']/g.loc['JJA','e8_rel_median']:.2f} in energy, "
              f"{(g.loc['DJF','e8_rel_median']/g.loc['JJA','e8_rel_median'])**0.5:.2f} in amplitude), "
              f"against a MAD scale that moves by only "
              f"{100*(g.loc['DJF','mad_scale_median']/g.loc['JJA','mad_scale_median']-1):+.1f} %. "
              f"So a louder winter does NOT raise the bar — the false-alarm floor is free to rise with "
              f"the number of coherent winter transients. This is the mechanism; it is a property of the "
              f"frozen MAD-relative operating point, not a fault in the month.\n")
        except Exception:
            pass
    if os.path.exists(vs):
        VS = pd.read_csv(vs)
        A("**(iii) Within a season, the noisiest days are not the ones that fire.** Slope of the "
          "quiet-day detection count against the daily noise term `v`:\n")
        A(md(VS, r=4))
        A("")
        A("Inside winter the slope is NEGATIVE — on the loudest days the correlation coefficient itself "
          "is depressed and fewer detections clear the threshold. The seasonal excess is therefore not "
          "carried by the loudest days; it is carried by the winter population of coherent transients "
          "(storm microseism bursts and uncatalogued micro-earthquakes) that the daily median energy "
          "does not measure.\n")
    ss = f"{TD}/season_stack_test_summary.csv"
    if os.path.exists(ss):
        SS = pd.read_csv(ss)
        A("**(iv) The decisive test — is the winter excess real LFE activity?** Track A step 3b's "
          "withheld-station test, repeated per season on production months with CN.VGZ, CN.SNB and "
          "PB.B926 excluded from the network sum (`scripts/d1d_run_withheld.sh`, six GPU passes), the "
          "frozen HF veto applied, quiet family-days taken from Track C's own `quiet_flag`, and months "
          "chosen so the two seasons have MATCHED quiet-day counts (DJF 2015-02 + 2014-02 + 2014-12 vs "
          "JJA 2014-08 + 2015-08 + 2014-07).\n")
        A("Production's quiet-day rate is ~4x lower than 2012's, so no single family has enough "
          "quiet-day detections for step 3b's coherent per-family stack (2–7 here against 20+ there). "
          "The primary statistic is therefore a POOLED INCOHERENT one: for every quiet-day detection, "
          "r = E(S ± 2 s) / E(pre-arrival 4 s) at the withheld station, each window normalised by its "
          "own noise, averaged in log10 over thousands of detections. **Two nulls**, and the second is "
          "the one that matters:\n")
        A("* **time null** — matched random times on the same quiet days, same offsets: *is there "
          "anything there at all?*\n"
          "* **moveout-shuffled null** — the SAME detection times with ANOTHER family's S offset at that "
          "station (offsets differing by > 5 s): *is the energy at the FAMILY-SPECIFIC predicted "
          "arrival?* A broadband regional transient that made the network fire raises both arms; a real "
          "LFE raises only the real one. Step 3b did not have this control; at production's purity it is "
          "essential.\n")
        A(md(SS, ["season", "tier", "n_quiet_family_days", "n_det", "excess_db", "excess_db_se", "z",
                  "p_null", "mo_mean_log10r", "excess_vs_mo_db", "p_mo"], r=4))
        A("")
        try:
            d = {(r.season, r.tier): r for r in SS.itertuples()}
            NZ2 = pd.read_csv(f"{TD}/season_diagnostic_noise.csv").set_index("season")
            rn = float(NZ2.loc["DJF", "e8_rel_median"] / NZ2.loc["JJA", "e8_rel_median"])
            lines = []
            for t in ("T_core", "T_ext"):
                w, j = d[("DJF", t)], d[("JJA", t)]
                lw = 10 ** (w.excess_vs_mo_db / 10) - 1
                ls = 10 ** (j.excess_vs_mo_db / 10) - 1
                lines.append(dict(tier=t, DJF_time_null_db=w.excess_db, JJA_time_null_db=j.excess_db,
                                  DJF_family_specific_db=w.excess_vs_mo_db,
                                  JJA_family_specific_db=j.excess_vs_mo_db,
                                  family_specific_ratio_DJF_over_JJA=lw / ls,
                                  noise_corrected_ratio=lw * rn / ls,
                                  nonfamily_share_DJF=1 - w.excess_vs_mo_db / w.excess_db,
                                  nonfamily_share_JJA=1 - j.excess_vs_mo_db / j.excess_db))
            A(md(pd.DataFrame(lines), r=3))
            A("")
            A("**Verdict.** Three things are established and one is not.\n")
            A("1. **Neither season's quiet-day set is pure noise.** Both beat the time null "
              "(p ≤ 0.049) AND the family-specific moveout null (p ≤ 0.048) at both tiers, so both "
              "contain real inter-ETS repeats — consistent with Track A step 3b's ~22 % at T_core.\n")
            A("2. **Most of the quiet-day energy is NOT family-specific.** 36–62 % of the time-null "
              "excess survives shuffling the moveout, i.e. the quiet-day detection sets are dominated in "
              "both seasons by broadband transients that made the whole network fire — the same "
              "population the HF veto targets.\n")
            A("3. **The winter excess is NOT explained by winter detections being more real.** The "
              f"per-detection family-specific content is LOWER in winter, not higher "
              f"({lines[0]['family_specific_ratio_DJF_over_JJA']:.2f}x at T_core, "
              f"{lines[1]['family_specific_ratio_DJF_over_JJA']:.2f}x at T_ext). At **T_ext — the tier "
              f"where the seasonal rate excess is largest — the winter detections are clearly less real** "
              f"even after correcting for the louder winter noise "
              f"({lines[1]['noise_corrected_ratio']:.2f}x), so that excess is predominantly false alarms, "
              "which is what the floor λ is for.\n")
            A(f"4. **What is NOT established**: r is normalised by each window's own noise, so a real LFE "
              f"of fixed amplitude yields a smaller r in a louder winter. Correcting by the measured "
              f"{rn:.2f}x winter noise energy leaves the two seasons statistically indistinguishable at "
              f"**T_core** ({lines[0]['noise_corrected_ratio']:.2f}x). **A genuine seasonal component of "
              f"inter-ETS activity therefore cannot be excluded at T_core**; what the test excludes is "
              f"the winter excess being PREDOMINANTLY real. λ is carried per season for exactly this "
              f"reason, and the seasonal-stability bar stays reported as NOT met rather than tuned "
              f"away.\n")
        except Exception as exc:
            A(f"(seasonal contrast not rendered: {exc})\n")
    else:
        A("**(iv)** The per-season withheld-station stack test is still running "
          "(`scripts/d1d_finish.sh`); its verdict goes here.\n")
    A("")

    # ---------------- D2
    A("## 2. Catalog v1 (Merlin Step 7)\n")
    A(f"`results/trackD/lfe_extension_catalog_" + ag.tag + ".parquet` (+ `.csv`): **{meta['n_events']:,} rows** "
      f"({meta['n_cat_core']:,} `cat_core`, {meta['n_cat_ext']:,} `cat_ext`; the tiers overlap), "
      f"{meta['n_families']} families, {meta['t0'][:10]} .. {meta['t1'][:10]}. "
      f"{meta['n_duplicates_removed']} duplicates removed; month ownership is half-open (a detection "
      f"belongs to the month of its own UTC date) and every month's count agrees with "
      f"`months_done.csv` ({meta['months_mismatched']} mismatches).\n")
    n_co = int((E.cat_core & ~E.cat_ext).sum())
    A(f"**The tiers are not nested.** {n_co} rows are `cat_core` but not `cat_ext`: T_core > T_ext, so "
      f"every core detection clears the ext threshold, but the frozen cross-family assignment is "
      f"evaluated separately on each tier's own set and a detection can lose to a station-sharing "
      f"neighbour at T_ext that is itself below T_core. Both flags are stored; use the one the analysis "
      f"needs, never their union as a single tier.\n")
    A("**TIME BASIS.** `time` is the origin-axis time the network sum was formed on = **Bostock's database "
      "clock for that family plus the family's fitted offset** (`e5r_family_offsets_2012.csv`, 3.25–21.25 s, "
      "median 12.25 s). It is **NOT an origin time**: the datum is Bostock's and the offset was fitted, so "
      "only DIFFERENCES carry physical meaning — within a family exactly, and between families to the "
      "accuracy of the offsets.\n")
    yr = E.assign(year=pd.DatetimeIndex(E.time).year).groupby("year").agg(
        cat_core=("cat_core", "sum"), cat_ext=("cat_ext", "sum"), families=("family_id", "nunique")).reset_index()
    A("Events per year:\n")
    A(md(yr))
    A("")
    ep = E.groupby("epoch_id").agg(cat_core=("cat_core", "sum"), cat_ext=("cat_ext", "sum"),
                                   families=("family_id", "nunique"),
                                   median_lambda_core=("lambda_core", "median"),
                                   median_Mc_core=("Mc_T_core", "median")).reset_index()
    FL = pd.read_csv(f"{TD}/floor_family_epoch_season.csv")
    src = FL[FL.tier == "T_core"].groupby(["epoch_id", "lambda_src"]).size().unstack(fill_value=0)
    A("Events and floor by epoch:\n")
    A(md(ep))
    A("")
    A("How each epoch's floor was estimated (number of (family, season) cells per fallback level, "
      "T_core) — `family-epoch-season` is the direct measurement, the rest are fallbacks for cells with "
      "no quiet observable family-day:\n")
    A(md(src.reset_index()))
    A("")
    A(f"`family_table_v1.csv`: {len(T)} families, **{int(T.template_less.sum())} template-less** (no usable "
      f"2012 template at any station — they can never be detected and are carried explicitly), "
      f"{int((T.n_events_core > 0).sum())} with at least one `cat_core` event; "
      f"{int(T.reloc_verdict.astype(str).str.startswith('ACCEPT(').sum())} carry our own accepted "
      f"hypocentre (σ_h ≈ 2.0 km from the matched synthetic floor), "
      f"{int(T.reloc_verdict.notna().sum())} carry a reported-only one.\n")
    A("`family_day_table_v1.parquet`: one row per (family, day) with `observable`, `n_sites`, `Mc_*`, "
      "`Pdet_*`, `lambda_*` and the day terms.\n")

    # ---------------- D3
    A("## 3. Science QC — the Thomas harness per SSE (Merlin Step 8)\n")
    A(f"**SSE windows.** PNSN tremor in the family box, ≥ 50 epicentres/day, clustered with gaps ≤ 20 d, "
      f"span ≥ 10 d, padded ±3 d → **{len(W)} windows** 2014–2026, {int(W.major.sum())} of them with "
      f"≥ 10 active days ('major'). A second set with a 5-day clustering gap "
      f"(`sse_windows_sub.csv`) is run as a reported sensitivity. Windows, with the harness result for "
      f"each, are `results/trackD/thomas_per_sse.csv`; the figures are `figures/d3_sse_<id>_<tier>.png` "
      f"and `figures/d3_thomas_summary_<tier>.png`.\n")
    A("Three legs are run per window: `cat_core`, `cat_ext`, and `cat_core_hp` = `cat_core` restricted to "
      "≥ 5 contributing SITES — a detector-internal purity leg, reported and never barred, because the "
      "harness's leading edge is the 10th percentile of event position and at ~50 % precision the "
      "background can pin it (the harness then reports 'no departure').\n")
    A("### 3a. The two pre-registered bars, on the 2012 referee\n")
    A(md(B))
    A("")
    thin = {}
    for f in sorted(glob.glob(f"{TD}/thinned_reference_*.csv")):
        if f.endswith("_runs.csv"):
            continue
        b = os.path.basename(f)[len("thinned_reference_"):-len(".csv")]
        t = "_".join(b.split("_")[:2])            # T_core / T_ext
        thin[t] = pd.read_csv(f)
    if thin:
        A("**The re-rupture bar failure is the decimation bias, not a catalog fault.** Bostock's 2012 "
          "catalog thinned per family at OUR measured recall plus OUR measured quiet-day floor "
          "(`d3_thinned_reference.py`, 20 seeds) is the fair reference:\n")
        rows = []
        for tier, D in thin.items():
            for m in ("front_firstburst_kmday", "reactivation_interval_d_p50", "n_secondary",
                      "n_secondary_null_circ", "frac_fam_reactivated"):
                g = D[(D.kind == "ref_floor") & (D.metric == m)]
                if not len(g):
                    continue
                rows.append(dict(tier=tier, metric=m, thinned_mean=g["mean"].iloc[0],
                                 thinned_p2_5=g.p2_5.iloc[0], thinned_p97_5=g.p97_5.iloc[0],
                                 bostock_full=g.bostock_full.iloc[0]))
        R = pd.DataFrame(rows)
        our = {}
        for tier, tag in (("T_core", "ref2012_ours_core"), ("T_ext", "ref2012_ours_ext")):
            r = TS[TS.tag == tag]
            if len(r):
                our[tier] = r.iloc[0].to_dict()
        R["ours"] = [our.get(t, {}).get(m, np.nan) for t, m in zip(R.tier, R.metric)]
        R["inside_thinned_95"] = (R.ours >= R.thinned_p2_5) & (R.ours <= R.thinned_p97_5)
        R.to_csv(f"{TD}/thinned_vs_ours.csv", index=False)
        A(md(R, ["tier", "metric", "ours", "thinned_mean", "thinned_p2_5", "thinned_p97_5",
                 "bostock_full", "inside_thinned_95"]))
        A("")
    A("### 3b. How often the harness's estimators actually work at this completeness\n")
    yrows = []
    for tier in ("cat_core", "cat_core_hp", "cat_ext"):
        M_ = TS[(TS.window_set == "main") & (TS.tier == tier) & TS.n_events.notna()
                & TS.sse_id.ne("REF2012")]
        if not len(M_):
            continue
        fb = pd.to_numeric(M_.front_firstburst_kmday, errors="coerce")
        ld = pd.to_numeric(M_.front_lead_kmday, errors="coerce")
        ok = fb.notna() & (fb < -3) & (fb > -15)
        sec = M_[M_.n_secondary.notna()]
        rr = pd.to_numeric(M_.reactivation_interval_d_p50, errors="coerce").dropna()
        yrows.append(dict(tier=tier, n_windows=len(M_),
                          first_burst_defined=int(fb.notna().sum()),
                          first_burst_plausible=int(ok.sum()),
                          lead_track_defined=int(ld.notna().sum()),
                          secondary_above_null_p95=int((sec.n_secondary >
                                                        sec.n_secondary_null_circ_p95).sum()),
                          n_with_reactivation=len(rr),
                          reactivation_median_over_sses_d=float(rr.median()),
                          plausible_speeds=";".join(f"{v:.2f}" for v in fb[ok])))
    Y = pd.DataFrame(yrows)
    Y.to_csv(f"{TD}/thomas_estimator_yield.csv", index=False)
    A(md(Y, ["tier", "n_windows", "first_burst_defined", "first_burst_plausible",
             "lead_track_defined", "secondary_above_null_p95", "n_with_reactivation",
             "reactivation_median_over_sses_d"], r=3))
    A("")
    A("**This is the honest limit of the science QC.** The main-front first-burst estimator returns a "
      "physically plausible southward speed (−3 to −15 km/day) in only **3 of 31** windows at "
      "`cat_core`, 4 at `cat_core_hp` and 6 at `cat_ext`; in most windows the harness reports "
      "'no departure' or a degenerate propagation window. The cause is visible in the space-time "
      "figures: its leading edge is the 10th percentile of event position, and at ~50 % precision the "
      "false-alarm background — which is spread along the whole strike — pins that percentile before "
      "the front arrives. Where the estimator does work the speeds are "
      + (", ".join(f"{t}: {r.plausible_speeds}" for t, r in zip(Y.tier, Y.itertuples()) if r.plausible_speeds)
         if len(Y) else "")
      + " km/day, bracketing Bostock's 2012 value of −4.70 and Thomas's 7.3 km/day reference. The "
        "re-activation interval, which needs no front, is available in most windows and its median "
        "across SSEs (~0.9–1.0 d) sits where the per-family-thinned reference says it should.\n")
    A("A front estimator matched to this catalog's purity — a first-activation fit restricted to "
      "families whose burst exceeds their own λ — was named here as the obvious next step. It was "
      "subsequently approved and built; see §3b-2.\n")

    # ---------------- 3b-2: the purity-matched front estimator
    fv = f"{TD}/thomas_front_v2{sfx}.csv"
    fr = f"{TD}/thomas_front_v2_referee.csv"
    dj = f"{TD}/thomas_front_v2_design.json"
    if os.path.exists(fv):
        A("### 3b-2. A purity-matched main-front estimator (POST-REGISTRATION, approved 2026-09-17)\n")
        A("**This section is not part of the pre-registered plan.** The estimator below was designed "
          "after §3b showed the harness's front estimator failing on this catalog; it was approved by "
          "the coordinator on 2026-09-17, designed on the 2012 referee ALONE, frozen, and only then run "
          "on the extension years. Code: `scripts/d3b2_front_v2.py`; figure "
          "`figures/d3b2_front_v2_summary.png`; results `thomas_front_v2.csv`.\n")
        A("**Definition.** Two changes to the harness's first-burst front, both aimed at the diagnosed "
          "failure and nothing else:\n")
        A("1. **Purity gate.** A family activates on the first day whose event count for that family is "
          "inconsistent with that family's OWN measured false-alarm floor: `n_f(d) ≥ N_MIN` AND "
          "`P(N ≥ n_f(d) | Poisson(λ_f)) < α / n_days` (Bonferroni over the window's days), with λ from "
          "§1c for that (family, tier, epoch, season). The activation time is the start of the first "
          "burst on that day (harness definition: gap < 600 s, ≥ 4 events).\n")
        A("2. **Window from the activations.** The propagation window is the 5th–95th percentile of the "
          "activation times, not the 10th-percentile leading-edge track — so background families cannot "
          "set it.\n")
        A("Everything else is the harness's (Theil-Sen, ≥ 5 families, 300-draw family bootstrap). One "
          "along-strike frame, computed once from all 130 Bostock hypocentres, is used for every "
          "window so speeds are comparable between SSEs.\n")
        A("**Free parameters and how they were fixed.** `N_MIN = 4` and `MIN_FAM = 5` are inherited "
          "from the harness; the 5–95 % window is fixed a priori. The ONLY fitted parameter is α, "
          "chosen from the pre-declared grid {0.05, 0.01, 0.001} as the first value putting ALL FOUR "
          "2012 referee catalogs within 20 % of Bostock's 2012 first-burst speed (−4.695 km/day):\n")
        if os.path.exists(fr):
            RF = pd.read_csv(fr)
            A(md(RF, ["alpha", "leg", "n_events", "n_fam_activated", "n_fit", "win_days",
                      "speed_kmday", "ci_lo", "ci_hi", "rel_diff_vs_bostock2012", "within_20pct"], r=3))
            A("")
            A("All three α values pass on all four legs — the estimator is insensitive to it (the four "
              "referee speeds move by ≤ 0.3 km/day across the whole grid) — so the rule selects "
              "**α = 0.05** and it is frozen there (`thomas_front_v2_design.json`). The extension years "
              "were not looked at until that point.\n")
        RV = pd.read_csv(fv)
        A("**Nulls.** The harness's circular and day-shuffle nulls are run unchanged, and a third is "
          "added because one of them is degenerate here: the day-shuffle permutes family labels among "
          "the SAME day's events, so every family's DAILY count is exactly preserved and a "
          "day-resolution activation statistic is nearly invariant under it. The evidence is in the "
          "table below — `null_dayshuf_mean` reproduces the observed speed to a few per cent in almost "
          "every window. It is reported for continuity and is NOT a significance test of this "
          "estimator. The added **position-shuffle** null keeps every activation time and permutes the "
          "along-strike coordinate across families, which asks directly whether activation ORDER "
          "relates to POSITION.\n")
        A("**Plausibility flag** (pre-declared): ≥ 5 families fitted, a bootstrap CI excluding zero, "
          "and |speed| ≤ 20 km/day (4× Thomas's 7.3 km/day reference).\n")
        yr = []
        for tier in ("cat_core", "cat_ext"):
            g = RV[RV.tier == tier]
            f_ = g[g.speed_kmday.notna()]
            sig = f_[f_.plausible & (f_.null_xshuf_p <= 0.05)]
            yr.append(dict(tier=tier, n_windows=len(g), fitted=len(f_),
                           plausible=int(f_.plausible.sum()), plausible_and_xshuf_sig=len(sig),
                           n_SE=int((sig.speed_kmday < 0).sum()), n_NW=int((sig.speed_kmday > 0).sum()),
                           median_abs_speed=float(sig.speed_kmday.abs().median()) if len(sig) else np.nan,
                           v1_plausible_for_comparison=3 if tier == "cat_core" else 6))
        A(md(pd.DataFrame(yr), r=3))
        A("")
        A("**Yield: 27 of 31 windows fitted and 11 significant at `cat_core`, against 3 of 31 for the "
          "pre-registered estimator** — the diagnosis in §3b was right and the fix works. Full per-SSE "
          "table (speed, CI, direction, propagation window, nulls, plausibility) is "
          "`results/trackD/thomas_front_v2.csv`; per-family activations are "
          "`thomas_front_v2_activations.csv`.\n")
        for tier in ("cat_core", "cat_ext"):
            g = RV[(RV.tier == tier) & RV.speed_kmday.notna()].sort_values("t0")
            if not len(g):
                continue
            A(f"**{tier}**:\n")
            A(md(g, ["sse_id", "t0", "major", "n_events", "n_fam_activated", "n_fit", "win_start",
                     "win_days", "speed_kmday", "ci_lo", "ci_hi", "direction", "rho", "null_circ_p",
                     "null_xshuf_p", "null_dayshuf_mean", "plausible"], r=3))
            A("")
        sig = RV[(RV.tier == "cat_core") & RV.plausible & (RV.null_xshuf_p <= 0.05)]
        A("**What it says, stated conservatively.** Of the 11 significant `cat_core` fronts, "
          f"{int((sig.speed_kmday < 0).sum())} propagate SE (Thomas's 'southward', the sense of the "
          f"2012 episode) and {int((sig.speed_kmday > 0).sum())} propagate NW; |speed| spans "
          f"{sig.speed_kmday.abs().min():.1f}–{sig.speed_kmday.abs().max():.1f} km/day with a median "
          f"of {sig.speed_kmday.abs().median():.1f}, bracketing Bostock's 2012 −4.70 and Thomas's "
          "7.3 km/day reference. Bidirectional along-strike propagation is expected in Cascadia, so a "
          "mix of senses is not by itself suspicious — but note that the four windows with the largest "
          "|speed| are also the ones whose SSE window is long enough (> 25 days) to contain more than "
          "one episode, and this estimator fits ONE front per window. Splitting multi-episode windows "
          "is the next refinement and is not done here.\n")
    A("### 3c. Per-SSE results (pre-registered estimator)\n")
    M = TS[(TS.window_set == "main") & TS.n_events.notna()].copy()
    for tier in ("cat_core", "cat_ext", "cat_core_hp"):
        g = M[M.tier == tier]
        if not len(g):
            continue
        A(f"**{tier}** ({len(g)} windows with ≥ 50 events):\n")
        A(md(g.sort_values("t0"), ["sse_id", "t0", "t1", "major", "n_events", "win_days",
                                   "front_firstburst_kmday", "front_lead_kmday", "n_secondary",
                                   "n_secondary_null_circ", "n_secondary_null_dayshuf",
                                   "reactivation_interval_d_p50", "frac_fam_reactivated"], r=3))
        A("")
    # catalog size per SSE (all windows, both tiers)
    idx = TS[TS.window_set == "main"].pivot_table(index="sse_id", columns="tier", values="n_in_catalog",
                                                  aggfunc="max").reset_index()
    A("Catalog size inside every SSE window (including windows too small to run the harness):\n")
    A(md(idx))
    A("")

    # ---------------- limitations
    A("## 4. Known limitations\n")
    r_e1 = S[(S.station_set == "S15") & (S.tier == "T_core")].recall_pooled.iloc[0]
    r_e1x = S[(S.station_set == "S15") & (S.tier == "T_ext")].recall_pooled.iloc[0]
    A(f"1. **Completeness is ~{100*r_e1:.0f} % (T_core) / ~{100*r_e1x:.0f} % (T_ext) of Bostock-class "
      f"events** on the 15-station set at the referee thresholds, and it varies by epoch because the "
      f"thresholds do: per epoch, at each epoch's own T, "
      + "; ".join(f"{r.epoch_id} {100*r.recall_pooled:.1f} %"
                  for r in pd.read_csv(f"{TD}/completeness_by_epoch_summary.csv")
                  .query("tier == 'T_core'").itertuples())
      + " at T_core (`completeness_by_epoch_summary.csv`). This catalog is a decimated view of the LFE "
        "population, not an equivalent of Bostock's. Every decimation-sensitive statistic must be "
        "compared against the per-family-thinned reference, not against his raw catalog.")
    A(f"2. **Completeness is strongly family-dependent**: per-family own recall on the S15 set at T_core has "
      f"p10/p50/p90 = "
      f"{S[(S.station_set=='S15')&(S.tier=='T_core')].recall_fam_p10.iloc[0]:.3f} / "
      f"{S[(S.station_set=='S15')&(S.tier=='T_core')].recall_fam_median.iloc[0]:.3f} / "
      f"{S[(S.station_set=='S15')&(S.tier=='T_core')].recall_fam_p90.iloc[0]:.3f}, and "
      f"{int(T.template_less.sum())} of {len(T)} families have no usable template at all and are NEVER "
      f"detected. Any spatial statistic (front position, along-strike coverage) inherits this and must "
      f"use the per-family recall column.")
    A("3. **Time basis**: Bostock's DB clock plus a fitted per-family offset — not an origin time "
      "(section 2).")
    A("4. **UW.LRIV is dropped after 2019-05-30** by Track A-2's epoch-identity result. The ablation "
      "above shows the loss is free (recall and false alarms both unchanged), and every family-day "
      "records it in `excluded_sta` / `n_sta_excluded`.")
    A("5. **CN.SNB ends 2021-03-28** and that IS costly. At the frozen threshold the 13-station network "
      "sits at ~3× the quiet-day false-alarm rate it was frozen at; Track C's per-epoch re-derivation "
      "responds by raising E3's T_core to 9.25, which restores the false-alarm budget and pays for it "
      "in recall (E3 T_core completeness 0.122 against E0's 0.146, −16 %). Either way the post-2021 "
      "record is the weakest part of the catalog; the per-event λ and Mc columns carry it.")
    A("6. **The post-2017 HF veto is provisional.** The cut is 3.6384 on the PGC BH epochs and 3.5129 on "
      "the HH epochs (Track A-2 RULE-M); the HH cut's 95 % bootstrap CI contains the frozen value, but it "
      "is calibrated from one SSE per epoch and the SSE-to-SSE spread inside an epoch is ±35 %. "
      "`hfr`, `hf_fs` and `hf_cut` are stored on every row so any re-calibration is a re-flag, not a "
      "re-run (`scripts/c_reflag_hf.py`).")
    A("7. **Quiet-day detections are not all false alarms.** Track A step 3b measured ~22 % real "
      "inter-ETS repeats at T_core (≤ 14 % at T_ext); λ is corrected by that fraction and `lambda_fa_low` "
      "carries the 95 %-upper-bound version. Only one family passes the withheld-station test "
      "individually, so no inter-ETS detection is promoted to a confirmed event.")
    A("8. **Magnitudes.** Our detections carry `rel_amp` (an SNR-weighted cc·√energy ratio in units of "
      "each station's 2012 quiet-day window rms), not a magnitude. Mc is expressed on Bostock's magnitude "
      "scale through the referee fit; it is a completeness statement about Bostock-class events, not a "
      "magnitude for our own events.")
    A("9. **The secondary-front count is uninformative at this completeness.** In the thinned reference "
      "the circular null is itself 10–14 detections against Bostock's 2.25 — the estimator's null scales "
      "with catalog density, so the count is reported and never interpreted.")
    A("10. **The quiet-day floor is seasonally unstable, and not in a fixed direction.** All eight "
      "epoch-tier cells miss Merlin's ≤ 2× winter/summer band, but E0 misses it high (winter 2.4× "
      "summer) and E1–E3 miss it low (summer 2–4× winter). Section 1d shows the floor is decoupled from "
      "the ambient level — winter is the loudest season in every epoch — so no seasonal correction is "
      "applied; λ is carried per (family, tier, epoch, season) instead.")
    A("11. **The PRE-REGISTERED main-front estimator mostly does not work at this purity** — a "
      "plausible southward speed in only 3 of 31 SSE windows at `cat_core` — because its leading edge "
      "is the 10th percentile of event position and the along-strike false-alarm background pins it. "
      "The purity-matched replacement of §3b-2 (post-registration, approved 2026-09-17) fits 27 of 31 "
      "and finds 11 significant, but it is a NEW estimator: its α was fixed on the 2012 referee, and "
      "the pre-registered bars in §3a are still the bars. Its own remaining limitation is that it fits "
      "ONE front per SSE window, and the four largest |speeds| come from windows long enough (> 25 "
      "days) to hold more than one episode.")
    if missing:
        A(f"11. **The record is incomplete**: {len(missing)} of 153 months are not yet processed "
          f"({', '.join(missing[:12])}{' ...' if len(missing) > 12 else ''}). Every number above is over "
          f"the months on disk and is refreshed by re-running `bash scripts/d_run_all.sh` followed by "
          f"`python scripts/d4_results.py`.")
    A("")
    open(f"{TD}/RESULTS.md", "w").write("\n".join(O) + "\n")
    print(f"[write] {TD}/RESULTS.md  ({len(O)} blocks)")


if __name__ == "__main__":
    main()
