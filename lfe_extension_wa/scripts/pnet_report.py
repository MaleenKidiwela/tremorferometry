#!/usr/bin/env python
"""PNET report -> results/PNET_RESULTS.md (numbers only; every section tolerates a missing stage)."""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
import pnet_common as PC

G = PC.OUT
MD = f"{W.RES}/PNET_RESULTS.md"
REG = ["ALL"] + W.REGIONS
SHORT = {"ALL": "ALL", "Royer_NW_WA": "NW", "Royer_SW_WA": "SW", "Olympic_2010": "Olympic", "Olympic_2010_armB": "Arm B"}
TAGS = [("w05", "off-the-shelf (SVI pnet_w05)"), ("ft1", "WA fine-tune (ft1)")]


def md(df, floatfmt=".3f"):
    try:
        return df.to_markdown(index=False, floatfmt=floatfmt)
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def rd(p, **k):
    return pd.read_csv(p, **k) if os.path.exists(p) else None


def main():
    L, S = [], {}          # lines, summary facts
    # ------------------------------------------------ A / features AUC
    secA = []
    for tag, name in TAGS:
        A = rd(f"{G}/stageA_{tag}_auc.csv")
        if A is None:
            secA.append(f"*{name}: not available.*\n"); continue
        A["region"] = A.region.map(SHORT)
        for con, lab in [("data_vs_control", "data vs scrambled control"), ("matched_own_vs_unmatched", "own-family truth-matched vs unmatched data"),
                         ("matched_3km_vs_unmatched", "≤ 3 km truth-matched vs unmatched data")]:
            t = A[(A.contrast == con) & A.eval_set.isin(["V2010_test", "V2008_test", "A2011", "fit_in_sample"])]
            t = t[["eval_set", "region", "n_pos", "n_neg", "pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05", "pnet_rel_max", "n_z2", "P"]]
            secA.append(f"**{name} — AUC, {lab}** (reference columns: step-1 `n_z2` and step-1 `P`)\n\n" + md(t.dropna(subset=["pnet_max"])) + "\n")
        h = A[(A.eval_set == "ALL_heldout")]
        for con in ["data_vs_control", "matched_own_vs_unmatched"]:
            for rg in ["ALL", "NW", "SW", "Olympic", "Arm B"]:
                x = h[(h.contrast == con) & (h.region == rg)]
                if len(x):
                    S[(tag, con, rg)] = (float(x.pnet_max.iloc[0]), float(x.n_z2.iloc[0]), float(x.P.iloc[0]))
    # ------------------------------------------------ B
    secB = []
    M = f"{G}/train/meta.parquet"
    if os.path.exists(M):
        m = pd.read_parquet(M)
        c = m.groupby(["kind", "n_horiz"]).size().rename("windows").reset_index()
        c["kind"] = c.kind.map({0: "positive", 1: "quiet-day negative", 2: "scrambled-control negative"})
        secB.append(f"Training windows (BUILD half of V2008 + V2010, {m.netsta.nunique()} training stations, {m.date.nunique()} days):\n\n" + md(c, ".0f") + "\n")
        S["n_train"] = len(m); S["n_pos"] = int((m.kind == 0).sum())
    H = rd(f"{G}/models/train_ft1.csv")
    if H is not None:
        secB.append("Fine-tune history (internal validation = every 10th date of the BUILD half; step 0 = SVI weights):\n\n" + md(H) + "\n")
        try:
            import torch
            z = torch.load(f"{G}/models/ft1.pt", weights_only=False, map_location="cpu")
            secB.append(f"Kept checkpoint: step {z['step']} (internal val AUC {z['val_auc']:.3f}).\n")
            S["ft_step"] = int(z["step"])
        except Exception:
            pass
    E = rd(f"{G}/stageB_eval_summary.csv")
    if E is not None:
        secB.append("Held-out days (V2008/V2010 TEST half, A2011): per-station AUC summary, stations with ≥ 50 labelled arrivals "
                    "(`held_out` = the 6 stations never used for fine-tuning):\n\n" + md(E) + "\n")
    B = rd(f"{G}/stageB_eval_by_station.csv")
    if B is not None:
        cols = ["netsta", "tag", "held_out", "n_pos", "auc_vs_quiet", "auc_vs_episode", "auc_vs_episode_V2008", "auc_vs_episode_V2010",
                "auc_vs_episode_A2011", "pos_frac03", "fp03_per_h", "fp05_per_h"]
        b = B[B.n_pos >= 50][cols].sort_values(["held_out", "netsta", "tag"], ascending=[False, True, True])
        secB.append("Per station (≥ 50 labelled arrivals on held-out days):\n\n" + md(b) + "\n")
        bq = B.groupby(["held_out", "tag"]).agg(fp03=("fp03_per_h", "median"), fp05=("fp05_per_h", "median")).reset_index()
        secB.append("False positives on QUIET days (test half): median over stations of declustered CF peaks per hour:\n\n" + md(bq) + "\n")
        for tag, _ in TAGS:
            x = B[(B.tag == tag) & (B.n_pos >= 50)]
            h = x[x.held_out]
            S[("B", tag)] = (h.auc_vs_quiet.median(), h.auc_vs_episode.median(), h.fp05_per_h.median(), len(h),
                             x[~x.held_out].auc_vs_episode.median())
    V = f"{G}/stageB_verdict.json"
    if os.path.exists(V):
        v = json.load(open(V))
        secB.append(f"Verdict rule (fixed before the run): fine-tune 'works' if, over held-out stations, the median AUC gain is > 0.02 vs quiet "
                    f"AND > 0 vs episode random times. Median ΔAUC vs quiet {v['median_delta_auc_vs_quiet']:.3f}, vs episode "
                    f"{v['median_delta_auc_vs_episode']:.3f} → **works = {v['works']}**.\n")
        S["works"] = v["works"]
    # ------------------------------------------------ C
    secC = []
    base = f"{W.RES}/step1_prob"
    runs = [("baseline (step 1)", base), ("+ PNET off-the-shelf", f"{G}/step1_w05"), ("+ PNET fine-tuned", f"{G}/step1_ft1")]
    rows_sl, rows_cq, rows_ec, rows_fi = [], [], [], []
    for name, d in runs:
        sl = rd(f"{d}/calibration_slopes.csv"); cq = rd(f"{d}/control_and_quiet.csv"); ec = rd(f"{d}/expected_count.csv")
        fi = rd(f"{d}/feature_importance.csv")
        if sl is None:
            continue
        for es in ["V2010_test", "V2008_test", "A2011"]:
            for rg in REG:
                x = sl[(sl.eval_set == es) & (sl.region == rg) & (sl.series == "data")]
                q = cq[(cq.eval_set == es) & (cq.region == rg) & (cq.days == "all")] if cq is not None else pd.DataFrame()
                e = ec[(ec.eval_set == es) & (ec.region == rg)] if ec is not None else pd.DataFrame()
                if not len(x) and not len(q):
                    continue
                rows_sl.append(dict(run=name, eval_set=es, region=SHORT[rg],
                                    slope_own=x.slope_own.iloc[0] if len(x) else np.nan, slope_3km=x.slope_nbr3.iloc[0] if len(x) else np.nan,
                                    matched_over_sumP_own=x.obs_over_exp_own.iloc[0] if len(x) else np.nan,
                                    data_P05=q.data_fracP_gt05.iloc[0] if len(q) else np.nan,
                                    ctl_P05=q.ctl_fracP_gt05.iloc[0] if len(q) else np.nan,
                                    ctl_P05_never_seen=q.ctl_loco_fracP_gt05.iloc[0] if len(q) else np.nan,
                                    data_sumP=q.data_sumP.iloc[0] if len(q) else np.nan,
                                    ctl_sumP_per_real=q.ctl_sumP_per_real.iloc[0] if len(q) else np.nan,
                                    spearman_own=e.spearman_own.iloc[0] if len(e) else np.nan,
                                    spearman_3km=e.spearman_nbr3.iloc[0] if len(e) else np.nan))
        if fi is not None:
            f = fi[fi.feature.str.startswith(("pnet", "n_sites_pnet"))]
            for r in f.itertuples():
                rows_fi.append(dict(run=name, feature=r.feature, perm_logloss_gain=r.perm_logloss_gain, auc_data_vs_ctl=r.auc_data_vs_ctl,
                                    rank=int(fi.perm_logloss_gain.rank(ascending=False)[fi.feature == r.feature].iloc[0]), of=len(fi)))
    if rows_sl:
        C = pd.DataFrame(rows_sl)
        secC.append("Same lfdr design (HistGradientBoosting data-vs-control + isotonic, 4 scrambled controls, π0 from quiet days, fit on the BUILD "
                    "half); only the feature list differs (+ `pnet_max`, `pnet_2nd`, `n_sites_pnet03`, `n_sites_pnet05`). Held-out sets only. "
                    "slope = calibration slope of the truth-matched fraction on P (target 0.8–1.2); P05 = fraction of candidates with P > 0.5; "
                    "Spearman = expected count ΣP per family-day vs the true array count.\n\n" + md(C) + "\n")
        for name, _ in runs:
            for es in ["V2010_test", "A2011", "V2008_test"]:
                for rg in ["ALL", "NW", "SW"]:
                    x = C[(C.run == name) & (C.eval_set == es) & (C.region == rg)]
                    if len(x):
                        S[("C", name, es, rg)] = x.iloc[0].to_dict()
    if rows_fi:
        secC.append("PNET features inside the refits (permutation log-loss gain, rank among all features):\n\n" + md(pd.DataFrame(rows_fi), ".4f") + "\n")
    # ------------------------------------------------ D
    secD = []
    T2 = rd(f"{W.RES}/merlin_tests/t2/summary.csv")
    rows = []
    if T2 is not None:
        t = T2[(T2.group == "region") & (T2["round"].isin([-1, 0]))].copy(); t["source"] = "T2 permanent MF (round 0)"
        rows.append(t)
    for tag, name in TAGS:
        x = rd(f"{G}/t2_{tag}/summary.csv")
        if x is not None:
            x = x[(x.group == "region") & (x["round"] == 0)].copy(); x["source"] = f"PNET beam, {name}"
            rows.append(x)
    if rows:
        D = pd.concat(rows, ignore_index=True)
        D["value"] = D.value.map(lambda v: SHORT.get(v, v))
        D = D[["source", "value", "series", "n_pairs", "med_cc_al", "frac_cc_gt07", "frac_cc_gt05", "own_best", "best_nbr3", "med_rank", "med_pur_own", "med_n_stacked"]]
        D = D.rename(columns={"value": "region"}).sort_values(["region", "source", "series"])
        secD.append("Withheld XU/XG station, template rebuilt from candidate times only; cc (2–8 Hz) vs the station's true template. "
                    "`data` = beam over the 3 nearest permanent sites, `data_n1` = nearest site only; `control*` = the same on the site-rolled "
                    "scrambled control; `ceiling` = stack at the true test-half detections.\n\n" + md(D) + "\n")
        for (src, rg, sr), g in D.groupby(["source", "region", "series"]):
            S[("D", src, rg, sr)] = float(g.med_cc_al.iloc[0])
    # ------------------------------------------------ summary
    def g(k, i=None, fmt="{:.2f}"):
        v = S.get(k)
        if v is None:
            return "–"
        v = v[i] if i is not None else v
        try:
            return fmt.format(v)
        except Exception:
            return str(v)
    summ = ["## Summary (plain language)", "",
            "PNET is a neural-network S-wave picker (PhaseNet-style) trained on Bostock's southern Vancouver Island LFEs. We asked whether its "
            "S-wave probability at the times a Washington family's template says the S should arrive gives independent evidence that a "
            "matched-filter repeat is real, and whether it can make candidate times on its own. The scrambled control is the same matched filter "
            "with each site's data shifted to a different time of the same day (so nothing lines up across sites except by chance).", ""]
    summ.append(f"- **Off-the-shelf PNET (A)** — AUC of `pnet_max`, held-out days pooled: data vs scrambled control {g(('w05','data_vs_control','ALL'),0)} "
                f"(NW {g(('w05','data_vs_control','NW'),0)}, SW {g(('w05','data_vs_control','SW'),0)}, Olympic {g(('w05','data_vs_control','Olympic'),0)}); "
                f"truth-matched vs unmatched {g(('w05','matched_own_vs_unmatched','ALL'),0)} (NW {g(('w05','matched_own_vs_unmatched','NW'),0)}). "
                f"For comparison the step-1 P gives {g(('w05','data_vs_control','ALL'),2)} / {g(('w05','matched_own_vs_unmatched','ALL'),2)}. "
                "On 1-component EHZ stations the SVI model returns ~0 (it was trained on 3-C only).")
    if ("B", "ft1") in S:
        b0, b1 = S.get(("B", "w05"), (np.nan,) * 5), S[("B", "ft1")]
        summ.append(f"- **Fine-tuned on WA (B)** — {S.get('n_pos','–')} labelled positive windows; held-out stations ({b1[3]}): median per-station AUC "
                    f"(labelled arrival vs quiet-day random time) {b0[0]:.2f} → {b1[0]:.2f}; vs episode-day random time {b0[1]:.2f} → {b1[1]:.2f}; "
                    f"training stations on held-out days vs episode {b0[4]:.2f} → {b1[4]:.2f}; quiet-day peaks > 0.5 per hour {b0[2]:.1f} → {b1[2]:.1f}. "
                    f"Pre-set rule 'works': **{S.get('works','–')}**.")
    if ("ft1", "data_vs_control", "ALL") in S:
        summ.append(f"- Fine-tuned `pnet_max` AUC data vs scrambled control {g(('ft1','data_vs_control','ALL'),0)} (NW {g(('ft1','data_vs_control','NW'),0)}, "
                    f"SW {g(('ft1','data_vs_control','SW'),0)}); truth-matched vs unmatched {g(('ft1','matched_own_vs_unmatched','ALL'),0)} "
                    f"(NW {g(('ft1','matched_own_vs_unmatched','NW'),0)}).")
    for name in ["baseline (step 1)", "+ PNET off-the-shelf", "+ PNET fine-tuned"]:
        x = S.get(("C", name, "V2010_test", "ALL")); y = S.get(("C", name, "A2011", "ALL")); z = S.get(("C", name, "V2010_test", "NW"))
        if x:
            summ.append(f"- **Step-1 refit (C) {name}**: calibration slope own V2010 test {x['slope_own']:.2f}, A2011 {y['slope_own'] if y else np.nan:.2f}; "
                        f"P > 0.5 data {x['data_P05']:.3f} vs scrambled control {x['ctl_P05']:.3f} (V2010 test); NW data {z['data_P05'] if z else np.nan:.3f} "
                        f"vs control {z['ctl_P05'] if z else np.nan:.3f}; expected-count Spearman V2010 {x['spearman_own']:.2f}, A2011 "
                        f"{y['spearman_own'] if y else np.nan:.2f}.")
    dl = []
    for rg in ["NW", "SW", "Olympic"]:
        a0 = S.get(("D", "T2 permanent MF (round 0)", rg, "data"))
        a1 = S.get(("D", "PNET beam, off-the-shelf (SVI pnet_w05)", rg, "data"))
        a2 = S.get(("D", "PNET beam, WA fine-tune (ft1)", rg, "data"))
        c1 = S.get(("D", "PNET beam, off-the-shelf (SVI pnet_w05)", rg, "control"))
        f = lambda v: "–" if v is None else f"{v:.2f}"
        dl.append(f"{rg} {f(a0)} (MF) / {f(a1)} (PNET off-the-shelf) / {f(a2)} (PNET fine-tuned), scrambled control {f(c1)}")
    summ.append("- **PNET candidate times for the withheld-station test (D)** — median recovered-template cc: " + "; ".join(dl) + ".")
    summ.append("")
    txt = ["# PNET as independent per-repeat evidence — WA extension (Merlin step 6)", "",
           f"Written by `scripts/pnet_report.py`; queue `scripts/pnet_queue.sh`; outputs `results/pnet/`. Numbers only; interpretation after review.", ""]
    txt += summ
    txt += ["## Input handling", "",
            "- Weights: `lfe_picker/results/p1/models/pnet_w05.pt` (P1b selected model: PhaseNet-style U-Net, 1.82 M parameters, 60-s windows at 20 Hz, "
            "3-C Z/N/E, predicted-label weight 0.5, trained 2008–2010 + 2012 on SVI stations).",
            "- Inference exactly as `lfe_picker/scripts/p_scan.cf_day` (step 20 s, central third kept, per-window 3-C peak normalisation, bf16).",
            "- Band/rate: the WA cache is the same frozen 20 Hz, 1–8 Hz Butterworth-4 zero-phase chain the model was trained on — no re-filtering.",
            "- Rotation to ZNE (SVI store was rotated with StationXML azimuths): PB EH1/EH2 rotated (B001 335.5/65.5, B003 239/329, B004 167.5/257.5, "
            "B005 90.5/180.5, B006 201/291, B007 36.5/126.5, B009 0/90, B010 228.5/318.5, B011 11/101, B013 31.5/301.5, B014 170/260, B017 206/296, "
            "B018 21.5/111.5, B019 305/35, B941 173.5/263.5, B943 301.5/31.5); US.NLWA 4/94; UW.HDW has one horizontal (EH2 at 45°) → Z only; "
            "1-C EHZ stations enter as (Z, 0, 0). Invalid samples (cache mask) set to 0.",
            "- Expected S at a station = candidate origin + the family's template offset at that station (template start = stack envelope peak − 0.5 s ≈ onset); "
            "value = max CF within ±0.5 s; scrambled-control candidates read at their site's rolled position (as the MF built them).", ""]
    txt += ["## A. Off-the-shelf check (and the same features from the fine-tune)", ""] + secA
    txt += ["## B. Fine-tune on Washington", "",
            f"Label rule: (family, station node) pairs whose build-half stack (≥ 30 events) has a clear-ish arrival (imp ≥ 3 and dur < 3 s; catalogue clear "
            f"= imp > 5); label weight 1.0 if clear, 0.5 otherwise; every build-half array detection of the family gives a Gaussian (σ 0.5 s) target at "
            f"t + template offset; loss only within ±3 s of a labelled arrival in positive windows. Negatives: quiet-day windows and windows at the "
            f"stations' scrambled-control candidate positions (build half), masked ±3 s around every known arrival. Held-out stations: "
            f"{', '.join(PC.HELD_OUT)} (each its own site); held-out days: V2008/V2010 test half and A2011.", ""] + secB
    txt += ["## C. Step-1 refit with PNET features", ""] + secC
    txt += ["## D. PNET-derived candidate times → T2 withheld-station template recovery", ""] + secD
    open(MD, "w").write("\n".join(txt) + "\n")
    print(f"wrote {MD}")


if __name__ == "__main__":
    main()
