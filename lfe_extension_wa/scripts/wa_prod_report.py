#!/usr/bin/env python
"""WA PRODUCTION report -> results/WA_PROD_RESULTS.md (plain summary first; every number with its file).
Robust to missing stages (a stage that did not run is reported as such).  "Scrambled control" = the same detector with each
station site's data shifted to a different time of the same day; "time-reversed control" = the same detector with every
template reversed in time (the main false-alarm reference)."""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P

O = P.OUT
REP = os.environ.get("WA_PROD_REPORT", f"{P.RES}/WA_PROD_RESULTS.md")


def rel(p):
    return p.replace(P.STUDY + "/", "")


def jl(p):
    try:
        return json.load(open(p))
    except Exception:
        return None


def rd(p, **k):
    try:
        return pd.read_csv(p, **k)
    except Exception:
        return None


def md(df, n=40):
    if df is None or not len(df):
        return "(none)\n"
    return df.head(n).to_markdown(index=False, floatfmt=".3f") + "\n"


def main():
    L = []
    L.append("# Washington production run (W0-W13) -- results\n")
    L.append(f"Auto-written by `scripts/wa_prod_report.py` at {pd.Timestamp.now('UTC').strftime('%Y-%m-%d %H:%MZ')}. "
             f"Outputs under `{rel(O)}/` and `{rel(P.CUT)}/`; queue `scripts/wa_prod_queue.sh` (log `logs/wa_prod_queue.log`). "
             "No downloads, no deletions, no git, nothing existing overwritten. 'Scrambled control' = the same detector with each "
             "station site's data shifted to a different time of the same day; 'time-reversed control' = every template reversed "
             "in time (the main false-alarm reference, user 2026-10-03).\n")
    # ---------------- gather
    CS = jl(f"{O}/setup/core_sets.json") or {}
    det = rd(f"{O}/cache/determinism.csv")
    guards = {os.path.basename(p)[6:-4]: open(p).read().strip() for p in glob.glob(f"{O}/status/GUARD_*.txt")}
    ev = rd(f"{O}/conf/eval_dev.csv"); es = rd(f"{O}/conf/eval_sealed.csv")
    bars = jl(f"{O}/conf/sealed_bars.json"); gate = jl(f"{O}/conf/pnet_gate.json"); rr = jl(f"{O}/conf/refit_risk.json")
    lro = rd(f"{O}/conf/olympic_lro.csv"); cut = rd(f"{O}/cutoff_check/false_share_by_cutoff.csv")
    K = rd(f"{O}/rates/k_factors.csv"); RM = rd(f"{O}/rates/rates_monthly.csv"); SC = rd(f"{O}/rates/step_checks.csv")
    w4 = rd(f"{O}/templates/tpl_qc_w4b.csv"); w4t = rd(f"{O}/templates/w4b_targets.csv")
    pil = jl(f"{O}/cut_pilot/PILOT.json"); stk = jl(f"{P.CUT}/STACKCHECK.json")
    ver = open(f"{P.CUT}/VERIFY.txt").read().strip() if os.path.exists(f"{P.CUT}/VERIFY.txt") else None
    csum = rd(f"{P.CUT}/cut_summary.csv")
    dele = open(f"{O}/delete/DELETE_CANDIDATES_wa.txt").read() if os.path.exists(f"{O}/delete/DELETE_CANDIDATES_wa.txt") else None
    mb = jl(f"{O}/post/member_beta.json"); b14 = jl(f"{O}/setup/b014_check.json")
    # ---------------- plain summary
    L.append("## In short\n")
    if CS:
        L.append(f"- **Fixed rate networks** (rule written before any run, `setup/core_sets.json`): headline E3 2015-2026 = "
                 f"{len(CS.get('E3', []))} stations ({', '.join(s.split('.')[1] for s in CS.get('E3', []))}); E2 2012-2014 = "
                 f"{len(CS.get('E2', []))}; E1 2007-2011 = {len(CS.get('E1', []))}.")
    if ev is not None and len(ev):
        h = ev[(ev.model == "WAVE") & (ev.label == "own") & (ev.net.isin(["E3", "all"]))]
        if len(h):
            L.append("- **Confidence (WAVE, own family, out-of-fold on fit days)**: " + "; ".join(
                f"{r.net} {r.reg} AUC {r.auc:.3f} (net_mad {r.auc_net_mad:.3f}), slope {r.slope:.2f}, matched/ΣP {r.matched_over_sumP:.2f}, "
                f"control P>0.5 scrambled {r.ctl_frac_P_gt05:.3f} / reversed {r.rev_frac_P_gt05:.3f}" for r in h.itertuples()) + " (`conf/eval_dev.csv`).")
    if bars:
        L.append("- **Sealed one-look check** (2007-01 headline, 2009-05 semi-independent; `conf/sealed_bars.json`): " + "; ".join(
            f"{k}: AUC {v.get('auc_WAVE_pooled', float('nan')):.3f} vs net_mad {v.get('auc_net_mad_pooled', float('nan')):.3f}"
            for k, v in bars.items() if k.endswith("|own")) + ".")
    else:
        L.append("- Sealed check: not opened (stage did not run).")
    if gate:
        L.append(f"- **PNET gate**: {'evaluated' if gate.get('evaluated') else 'not evaluated'}; adopt = {gate.get('adopt')} "
                 f"(dAUC {gate.get('dAUC')}, held-out-station block {gate.get('dAUC_heldout_station_block')}; `conf/pnet_gate.json`).")
    if cut is not None and len(cut):
        h = cut[(cut.run == "prod") & (cut.net == "core") & (cut.cut == P.P_CUT) & (cut.day_class.isin(["quiet", "tremor60"]))]
        h = h[h.label == h.region.map({"NW": "own", "SW": "own", "ArmB": "own", "OlyA": "nbr3"})]
        L.append("- **False share at P >= 0.75 (core network)**, time-reversed / scrambled control vs data: " + "; ".join(
            f"{r.region} {r.day_class}: {r.rev_control}/{r.scrambled_control} of {r.data} ({r.false_share_rev:.3f} / {r.false_share_scrambled:.3f})"
            for r in h.itertuples()) + " (`cutoff_check/false_share_by_cutoff.csv`).")
    if rr:
        L.append(f"- **Refit risk** (P_own_WAVE >= 0.75 but the older model < 0.2): {rr.get('refit_risk', float('nan')):.4f} -> window floor "
                 f"P_max >= {rr.get('floor')} (`conf/refit_risk.json`).")
    if csum is not None and len(csum):
        L.append(f"- **Windows for dv/v**: {int(csum.n_windows.sum()):,} windows, {csum.cut_GB.sum():.1f} GB from {csum.raw_GB.sum():.0f} GB raw "
                 f"({len(csum)} months); verify: {ver}; stack check pass = {stk.get('pass') if stk else None} (`{rel(P.CUT)}/`).")
    elif pil:
        L.append(f"- Windows: pilot only ({pil.get('extrapolation')}); full cut not run.")
    if dele:
        L.append("- **Delete-candidate list** (nothing deleted): " + dele.splitlines()[1].strip() + " (`delete/`).")
    L.append("- Caveat: P = probability of matching the array catalogue on its scanned days. Between ETS windows it is a lower "
             "bound on 'real', not a measure of it; inter-ETS rates are a control-subtracted bracket.\n")
    # ---------------- sections
    L.append("## W0-W3 freeze, completeness, health, caches\n")
    L.append(f"- Frozen inputs: `w0/frozen_inputs.csv`; re-hash at the end: {sorted(glob.glob(f'{O}/w0/verify_*.csv'))[-1:] }")
    if det is not None:
        L.append(f"- Cache determinism: {int(det.identical.sum())}/{len(det)} station-days rebuilt bit-identical (`cache/determinism.csv`).")
    for k, v in sorted(guards.items()):
        L.append(f"- guard {k}: {v}")
    if b14:
        L.append(f"- B014 2023 check: {b14}")
    L.append("")
    L.append("## W4 templates\n")
    lib = rd(f"{O}/setup/tpl_library.csv", usecols=["tag", "node", "fid"])
    if lib is not None:
        L.append(f"- W4a library: {len(lib)} (family, node) templates on {lib.node.nunique()} nodes; by tag {lib.tag.value_counts().to_dict()} "
                 "(`setup/tpl_library.csv`; 'v10b' = the step-4b templates rebuilt from the V2010 build half so A2011 stays leak-free).")
    nm = rd(f"{O}/setup/tpl_nodemap.csv")
    L.append("- Node map (W2 verdicts):\n\n" + md(nm))
    if w4 is not None and len(w4):
        u = w4[w4.use.astype(bool)]
        L.append(f"- W4b: {len(w4t) if w4t is not None else 0} target nodes; {len(w4)} (family, node) stacks, {len(u)} used on {u.node.nunique()} nodes; "
                 f"median cc of a used template vs the scrambled-control-time stack {u.cc_ctl_stack.median():.2f} (`templates/tpl_qc_w4b.csv`).")
        L.append("\n" + md(u.groupby(["kind", "netsta"]).agg(n_templates=("fid", "size"), med_cc_ctl=("cc_ctl_stack", "median")).reset_index()))
    L.append("## W5-W8 confidence\n")
    if ev is not None:
        L.append(md(ev[(ev.model.isin(["WAVE", "FULL", "MF_ONLY"]))][["net", "reg", "label", "model", "n", "n_pos", "auc", "auc_net_mad",
                                                                     "slope", "matched_over_sumP", "ctl_frac_P_gt05", "rev_frac_P_gt05"]], 80))
    af = rd(f"{O}/conf/eval_dev_by_truth_array_frac.csv")
    L.append("By truth_array_frac (dev, WAVE own):\n\n" + md(af, 60))
    L.append("Olympic leave-region-out (credible if within 0.03 AUC of out-of-fold):\n\n" + md(lro))
    if es is not None:
        L.append("Sealed (one look):\n\n" + md(es[es.model.isin(["WAVE", "FULL", "MF_ONLY"])][["eval", "net", "reg", "label", "model", "n", "n_pos",
                                                                                         "auc", "auc_net_mad", "slope", "matched_over_sumP",
                                                                                         "ctl_frac_P_gt05", "rev_frac_P_gt05"]], 80))
    if mb:
        L.append(f"- W9 member probability: beta {mb.get('beta')} from {mb.get('n_groups')} truth groups (`post/member_beta.json`).")
    L.append("\n## W10 false share by cut-off\n")
    if cut is not None:
        L.append(md(cut[(cut.run == "prod") & (cut.net == "core") & (cut.cut.isin([0.5, 0.75, 0.9]))]
                    [["region", "label", "day_class", "cut", "data", "rev_control", "scrambled_control", "false_share_rev",
                      "false_share_scrambled", "families_covered"]], 120))
    L.append("## W11 rates\n")
    if K is not None:
        L.append(f"- k factors: {len(K)} (unit, era set) rows; median k by era set and region:\n\n" +
                 md(K.groupby(["net", "unit_kind", "reg"]).agg(n=("unit", "size"), k_med=("k", "median"), k_raw_med=("k_raw", "median")).reset_index()))
    if RM is not None and len(RM):
        t = RM.groupby(["unit_kind", "label"]).agg(units=("unit", "nunique"), months=("month", "nunique"), net_count=("net_count", "sum"),
                                                     rate_scaled=("rate_scaled", "sum")).reset_index()
        L.append("- Rate tables `rates/rates_daily.parquet`, `rates/rates_monthly.csv`:\n\n" + md(t))
    if SC is not None and len(SC):
        L.append(f"- Step checks: {int(SC.flag.sum())} of {len(SC)} flagged (|z| > 2 vs the same calendar windows of other years); "
                 "a step next to an ETS cannot be told from activity (`rates/step_checks.csv`).\n\n" + md(SC[SC.flag], 30))
    L.append("## W12 windows\n")
    L.append(f"- pilot: {pil}\n- verify: {ver}\n- stack check: {json.dumps(stk)[:1500] if stk else None}\n")
    L.append("## W13 delete-candidate list (nothing deleted)\n")
    L.append("```\n" + (dele or "not produced") + "```\n")
    L.append("## Deviations from the plan (with reasons)\n")
    dv = f"{O}/DEVIATIONS.md"
    L.append(open(dv).read() if os.path.exists(dv) else "(see the launch note)\n")
    open(REP, "w").write("\n".join(L) + "\n")
    print(f"wrote {REP}")


if __name__ == "__main__":
    main()
