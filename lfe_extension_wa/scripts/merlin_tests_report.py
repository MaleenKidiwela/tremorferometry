#!/usr/bin/env python
"""Merlin tests report -> results/MERLIN_TESTS_RESULTS.md (tolerates a missing T1 or T2: reports what exists).

  WA_SETUP=results/merlin_tests/t1_setup python scripts/merlin_tests_report.py
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

MT = f"{W.RES}/merlin_tests"
OUT = f"{W.RES}/MERLIN_TESTS_RESULTS.md"
RG = [("Royer_NW_WA", "NW"), ("Royer_SW_WA", "SW"), ("Olympic_2010", "Olympic"), ("Olympic_2010_armB", "Arm B"), ("ALL", "ALL")]
RUNS = [("val", "existing trial run (clear only, V2008+V2010)"), ("mt1_before", "before: clear only (rerun)"),
        ("mt1_all", "after: ALL templates")]


def f2(x, n=3):
    try:
        return "–" if x is None or not np.isfinite(float(x)) else f"{float(x):.{n}f}"
    except Exception:
        return "–"


def i0(x):
    try:
        return "–" if not np.isfinite(float(x)) else f"{int(x)}"
    except Exception:
        return "–"


def load_t1():
    S = {}
    for tag, _ in RUNS:
        p = f"{MT}/t1/score_{tag}/summary.csv"
        if os.path.exists(p):
            S[tag] = pd.read_csv(p)
    return S


def t1_tables(S, tier, kind="cat"):
    what = ("catalogue = tier & HF veto & cross-family rule" if kind == "cat" else
            "RAW = tier only, no HF veto / cross-family rule (with every template nearly all families share >= 50 % of "
            "their stations, so the frozen cross-family rule becomes a competition among all families)")
    L = [f"**Tier `{tier}`, {kind}** ({what}; the scrambled control = tier only). "
         "`obs fam` = families with ≥ 1 observable scanned episode day; `truth obs / all` = own array detections on "
         "observable days / on all scanned days (scored half); recall = on observable days, `recall all` = over ALL "
         "scanned truth (so a run that makes more families observable is not penalised). Levels own / ≤3 km / any.", "",
         "| region | run | obs fam | truth obs / all | recall own · 3km · any | recall all (own) | prec own · 3km · any | "
         "control prec own · 3km · any | data / control per fam-day | sep ≥1.5 · 2 · 3 · 5 (data / control) |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for rg, lab in RG:
        for tag, desc in RUNS:
            if tag not in S:
                continue
            x = S[tag]
            x = x[(x.window == "ALL") & (x.region == rg)]
            if not len(x):
                continue
            r = x.iloc[0]
            k = f"{kind}_{tier}"; c = f"ctl_{tier}"
            sep = " · ".join(f"{i0(r.get(f'sep{s:g}_{tier}_data'))}/{i0(r.get(f'sep{s:g}_{tier}_ctl'))}" for s in (1.5, 2, 3, 5))
            L.append(f"| {lab} | {desc} | {i0(r.n_fam_observable)} | {i0(r.n_truth_test)} / {i0(r.n_truth_scored_all)} | "
                     f"{f2(r.get(f'recall_own_{k}'))} · {f2(r.get(f'recall_nbr3_{k}'))} · {f2(r.get(f'recall_any_{k}'))} | "
                     f"{f2(r.get(f'recall_all_own_{k}'))} | "
                     f"{f2(r.get(f'prec_own_{k}'))} · {f2(r.get(f'prec_nbr3_{k}'))} · {f2(r.get(f'prec_any_{k}'))} | "
                     f"{f2(r.get(f'prec_own_{c}'))} · {f2(r.get(f'prec_nbr3_{c}'))} · {f2(r.get(f'prec_any_{c}'))} | "
                     f"{f2(r.get(f'data_per_famday_{tier}'), 2)} / {f2(r.get(f'ctl_per_famday_{tier}'), 2)} | {sep} |")
    return L


def t1_window_table(S, tier):
    L = ["| window | region | run | obs fam | recall own (obs) | recall all (own) | prec own · 3km · any | control prec own | data / ctl per fam-day | sep≥2 data / ctl |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for w in ["V2008", "V2010", "A2011"]:
        for rg, lab in RG:
            for tag, desc in RUNS:
                if tag not in S:
                    continue
                x = S[tag]; x = x[(x.window == w) & (x.region == rg)]
                if not len(x):
                    continue
                r = x.iloc[0]; k = f"cat_{tier}"
                L.append(f"| {w} | {lab} | {tag} | {i0(r.n_fam_observable)} | {f2(r.get(f'recall_own_{k}'))} | "
                         f"{f2(r.get(f'recall_all_own_{k}'))} | {f2(r.get(f'prec_own_{k}'))} · {f2(r.get(f'prec_nbr3_{k}'))} · "
                         f"{f2(r.get(f'prec_any_{k}'))} | {f2(r.get(f'prec_own_ctl_{tier}'))} | "
                         f"{f2(r.get(f'data_per_famday_{tier}'), 2)} / {f2(r.get(f'ctl_per_famday_{tier}'), 2)} | "
                         f"{i0(r.get(f'sep2_{tier}_data'))} / {i0(r.get(f'sep2_{tier}_ctl'))} |")
    return L


def t1_summary(S):
    if "mt1_all" not in S or "mt1_before" not in S:
        return ["- T1: not available (see `results/merlin_tests/T1_FAILED` or the queue log)."]
    out = []
    for rg, lab in RG[:4]:
        b = S["mt1_before"]; a = S["mt1_all"]
        b = b[(b.window == "ALL") & (b.region == rg)]; a = a[(a.window == "ALL") & (a.region == rg)]
        if not len(a) or not len(b):
            continue
        b, a = b.iloc[0], a.iloc[0]
        out.append(f"- T1 {lab}: observable families {i0(b.n_fam_observable)} → {i0(a.n_fam_observable)}; own-family recall over ALL "
                   f"scanned truth (core) {f2(b.recall_all_own_cat_core)} → {f2(a.recall_all_own_cat_core)} (≤3 km: "
                   f"{f2(b.recall_all_nbr3_cat_core)} → {f2(a.recall_all_nbr3_cat_core)}); own precision {f2(b.prec_own_cat_core)} → "
                   f"{f2(a.prec_own_cat_core)} (scrambled control's chance level {f2(a.prec_own_ctl_core)}); detections per family-day "
                   f"data / control {f2(b.data_per_famday_core, 1)} / {f2(b.ctl_per_famday_core, 1)} → "
                   f"{f2(a.data_per_famday_core, 1)} / {f2(a.ctl_per_famday_core, 1)}; families ≥ 2× above the control "
                   f"{i0(b.sep2_core_data)} (control {i0(b.sep2_core_ctl)}) → {i0(a.sep2_core_data)} (control {i0(a.sep2_core_ctl)}).")
    return out


def load_t2():
    p = f"{MT}/t2/summary.csv"
    return pd.read_csv(p) if os.path.exists(p) else None


def t2_tables(T):
    L = ["| group | value | series | round | pairs | families | median cc (±2 s) | cc > 0.7 | cc > 0.5 | own best | top 3 | best ≤3 km | median rank / competitors | candidate purity (own) | windows stacked |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in T.itertuples():
        L.append(f"| {r.group} | {r.value} | {r.series} | {r.round} | {r.n_pairs} | {r.n_families} | {f2(r.med_cc_al)} | "
                 f"{f2(r.frac_cc_gt07, 2)} | {f2(r.frac_cc_gt05, 2)} | {f2(r.own_best, 2)} | {f2(r.top3, 2)} | {f2(r.best_nbr3, 2)} | "
                 f"{f2(r.med_rank, 1)} / {f2(r.med_competitors, 0)} | {f2(r.med_pur_own, 3)} | {f2(r.med_n_stacked, 0)} |")
    return L


def t2_summary(T):
    if T is None:
        return ["- T2: not available (see `results/merlin_tests/T2_FAILED` or the queue log)."]
    out = []
    g = T[T.group == "region"]
    rmax = int(g[g.series == "data"]["round"].max()) if (g.series == "data").any() else 0
    for val, lab in RG:
        h = g[g.value == val]
        if not len(h):
            continue
        d0 = h[(h.series == "data") & (h["round"] == 0)]; dl = h[(h.series == "data") & (h["round"] == rmax)]
        c0 = h[(h.series == "control") & (h["round"] == rmax)]; ce = h[h.series == "ceiling"]
        if not len(d0):
            continue
        d0, dl = d0.iloc[0], dl.iloc[0] if len(dl) else d0.iloc[0]
        out.append(f"- T2 {lab} ({int(d0.n_pairs)} family–station pairs, {int(d0.n_families)} families): recovered-vs-true template "
                   f"median cc {f2(d0.med_cc_al, 2)} (round 0) → {f2(dl.med_cc_al, 2)} (round {rmax}); cc > 0.7 in "
                   f"{f2(dl.frac_cc_gt07, 2)} of pairs; own family best-matching in {f2(dl.own_best, 2)} (top 3: {f2(dl.top3, 2)}; "
                   f"best is a ≤3-km neighbour: {f2(dl.best_nbr3, 2)}). Scrambled control: median cc "
                   f"{f2(c0.iloc[0].med_cc_al, 2) if len(c0) else '–'}, own best {f2(c0.iloc[0].own_best, 2) if len(c0) else '–'}. "
                   f"Ceiling (stack at the true test-half detections): median cc {f2(ce.iloc[0].med_cc_al, 2) if len(ce) else '–'}.")
    return out


def main():
    S = load_t1()
    T = load_t2()
    L = ["# Merlin's two deciding tests — T1 all templates, T2 withheld-station template recovery", "",
         f"Written {time.strftime('%Y-%m-%d %H:%MZ', time.gmtime())} by `scripts/merlin_tests_report.py` (numbers only; "
         "interpretation to be added after review). Queue `scripts/merlin_tests_queue.sh`; outputs `results/merlin_tests/`; "
         "log `NOTES.md`.", "", "## Summary (plain language)", "",
         "T1 asks: if the network matched filter keeps EVERY station template that stacked ≥ 30 array-era events (weighted by its "
         "signal-to-noise weight) instead of only the 386 clear ones, do we find more of the known array detections while the "
         "scrambled control (stations shifted to different times) stays flat? T2 asks: if we pretend an array station is new, "
         "can we rebuild its template for a family only from times the permanent network picks, i.e. without the station?", ""]
    L += t1_summary(S) + t2_summary(T)
    L += ["", "## T1 — all templates vs clear only (V2008 + V2010 test half, A2011 all)", "",
          "Same chain as WA_TRIAL: frozen 1–8 Hz network matched filter, site merge B, ≥ 2 sites per family-day, "
          "record ≥ 8.0 MAD, T_ext 8.25 / T_core 8.75, plus per-epoch tiers (SVI FA-only rule on quiet days; A2011 has "
          "none → V2010's), HF veto calibrated on V2008/V2010 build-half own matches, frozen cross-family rule on the "
          "stations of the templates in use, one-to-one ±2 s own-family scoring. Two site-rolled scrambled controls "
          "(control 1 is the chance level; control 2 applies the separation rule to a control). One rule added for BOTH "
          "before and after: a station-day with > 5 % of its samples guarded is treated as absent (under the "
          "all-stations-valid network rule one half-dead station, e.g. SP2 on 2010-08-20, blanked half of the day for "
          "all 274 families once every template is used). The existing trial run (`val`) is scored by the same scorer "
          "for reference (1 control → no control separation column).", ""]
    if S:
        L += t1_tables(S, "core") + [""] + t1_tables(S, "core", "raw") + [""] + t1_tables(S, "core_ep") + [""] + \
            t1_tables(S, "ext") + ["", "Per window (tier core, catalogue):", ""] + t1_window_table(S, "core")
        for tag, _ in RUNS:
            p = f"{MT}/t1/score_{tag}/tiers_by_epoch.csv"
            if os.path.exists(p):
                tt = pd.read_csv(p)
                L += ["", f"Per-epoch tiers `{tag}`: " + "; ".join(f"{r.window} T_ext_ep {f2(r.T_ext_ep, 2)} / T_core_ep {f2(r.T_core_ep, 2)} "
                                                              f"(FA at 8.75 = {f2(r.fa_at_8p75, 2)})" for r in tt.itertuples())]
        L += ["", "Per family: `results/merlin_tests/t1/score_<run>/by_family.csv` (every family-window, every column)."]
    else:
        L += ["T1 results not available."]
    L += ["", "## T2 — withheld-station template recovery (array era: V2008 XU, V2010 XU + XG)", "",
          "Pairs = (family, window, XU/XG station with a clear truth template from the family's BUILD-half array detections, "
          "≥ 30 events), ≤ 4 stations per family-window drawn at random. Round 0: candidates = top-K (2 per scored hour, ±4 s "
          "declustered) peaks of the family's PERMANENT-network matched filter (all ≥ 30-event templates; identical to "
          "stack_activity's candidates), TEST-half origin times only; the withheld station is stacked there and a template is "
          "cut by the template rule. Rounds 1–2: the station joins the network with its recovered template (its own site, "
          "weight n_ch × SNR weight) and candidates are re-picked. cc = 2–8 Hz (designed at fs = 20 Hz) normalised cc between "
          "the recovered stack and the TRUE template (build half), best lag within ±2 s; rank among all families' truth "
          "templates at that station (any lag). Scrambled control = the same procedure on the site-rolled network. Ceiling = "
          "the station stacked at the family's own TEST-half array detections. `perm_clear_sites` = median number of permanent "
          "sites with a CLEAR template for the family (the carrying network).", ""]
    if T is not None:
        L += t2_tables(T[T.group == "region"]) + ["", "By number of permanent sites carrying the family (clear templates):", ""] + \
            t2_tables(T[T.group == "perm_clear_sites"])
        L += ["", "Per pair and round: `results/merlin_tests/t2/pair_results.csv`; candidates `round<r>_cand.parquet`; "
              "recovered stacks/templates `round<r>.npz` + `round<r>_qc.csv`; truth `truth_qc.csv`, `truth_stacks.npz`."]
    else:
        L += ["T2 results not available."]
    open(OUT, "w").write("\n".join(L) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
