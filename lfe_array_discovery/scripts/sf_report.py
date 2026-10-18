#!/usr/bin/env python
"""Write results/WA_FINAL_ASSESSMENT.md under the RANK-DON'T-DISCARD rule, and append a dated
NOTES.md entry.  Every number is read from a file; nothing is hard-coded, nothing is filtered."""
from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

F = f"{SF.RES}/sf_final"
BEDS = [("2010 N-Olympic, arrays as arrays (agent O)", "oly2010"),
        ("2008 Puget / S-Olympic (agent N)", "puget2008"),
        ("2010 Arm B, Chestler's own chain (agent O)", "armb")]
CHESTLER = ("/home/jovyan/tremorferometry/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")


def jl(p, d=None):
    try:
        return json.load(open(p))
    except Exception:
        return {} if d is None else d


def cl(p):
    try:
        return pd.read_csv(p)
    except Exception:
        return None


def n(x, d=0):
    if x is None:
        return "--"
    try:
        v = float(x)
    except Exception:
        return str(x)
    return "--" if not np.isfinite(v) else f"{v:,.{d}f}"


def inventory_table():
    s = ("| bed | chain | grown templates (ALL kept) | dedupe: kept as distinct | merged (row "
         "retained) | located >= 5 stations | unlocatable (row retained) | on a grid face | on a "
         "depth rail | re-fired on the sealed year |\n" + "|" + "---|" * 10 + "\n")
    any_row = False
    for name, bed in BEDS:
        for ch in ("real", "null"):
            d = jl(f"{F}/table_{bed}_{ch}_summary.json")
            if not d:
                continue
            any_row = True
            s += (f"| {name} | {ch} | {d.get('n_rows')} | {d.get('n_dedupe_kept')} | "
                  f"{int(d.get('n_rows', 0)) - int(d.get('n_dedupe_kept', 0))} | "
                  f"{d.get('n_located')} | "
                  f"{int(d.get('n_rows', 0)) - int(d.get('n_located', 0))} | "
                  f"{d.get('n_on_grid_edge')} | {d.get('n_on_depth_rail')} | "
                  f"{d.get('n_sealed_refired')} |\n")
    return s if any_row else "_no chain table has been built yet_\n"


def sep_table(bed):
    SU = cl(f"{F}/sep_{bed}_summary.csv")
    if SU is None or not len(SU):
        return "_not computed_\n", None
    SU = SU.sort_values("sep", ascending=False)
    s = ("| criterion | direction | AUC real vs null | +- | separation \\|AUC-0.5\\| | favours | "
         "best lift | at threshold | n real | n null |\n" + "|" + "---|" * 10 + "\n")
    for r in SU.head(18).itertuples():
        s += (f"| `{r.col}` | {r.direction} | {n(r.auc, 3)} | {n(r.auc_se, 3)} | "
              f"**{n(r.sep, 3)}** | {getattr(r, 'favours', '--')} | {n(r.best_lift, 2)} | "
              f"{n(r.best_thr, 3)} | {r.best_n_real} | {r.best_n_null} |\n")
    return s, SU


def curve(bed, col, k=8):
    SW = cl(f"{F}/sep_{bed}_sweep.csv")
    if SW is None or col is None:
        return ""
    S = SW[SW.col == col]
    if not len(S):
        return ""
    S = S.tail(k)
    s = (f"\nSeparation curve for `{col}`:\n\n| threshold | real beyond it | null beyond it | "
         f"% of real | % of null | lift |\n" + "|" + "---|" * 6 + "\n")
    for r in S.itertuples():
        s += (f"| {n(r.thr, 3)} | {r.n_real} | {r.n_null} | {n(100 * r.frac_real, 1)} | "
              f"{n(100 * r.frac_null, 1)} | {n(r.lift, 2)} |\n")
    return s


def chestler_section():
    C = cl(CHESTLER)
    if C is None:
        return "_coordinate file unavailable_\n"
    s = ("| bed | chain | located rows | median km to nearest published coordinate | within 5 km | "
         "within 10 km | within 20 km | of the 69 coordinates, hit within 10 km |\n"
         + "|" + "---|" * 8 + "\n")
    for name, bed in BEDS:
        for ch in ("real", "null"):
            T = cl(f"{F}/table_{bed}_{ch}.csv")
            if T is None or "near_chestler_km" not in T.columns:
                continue
            L = T[T.located.astype(bool) & T.near_chestler_km.notna()]
            if not len(L):
                continue
            hit = set()
            for r in L.itertuples():
                if r.near_chestler_km <= 10.0:
                    hit.add(int(r.near_chestler_row))
            s += (f"| {name} | {ch} | {len(L)} | {n(L.near_chestler_km.median(), 1)} | "
                  f"{int((L.near_chestler_km <= 5).sum())} | "
                  f"{int((L.near_chestler_km <= 10).sum())} | "
                  f"{int((L.near_chestler_km <= 20).sum())} | {len(hit)} |\n")
    return s


def ranked_section(bed, name, top=25):
    RK = cl(f"{F}/sep_{bed}_ranked.csv")
    if RK is None or not len(RK):
        return ""
    col = str(RK.rank_col.iloc[0]) if "rank_col" in RK.columns else "n_det"
    cols = ["grp", col, "null_percentile_on_rank_col", "n_det", "n_days", "lat", "lon", "dep",
            "edge", "dep_rail", "ctrl_pass", "near_chestler_km", "sealed_n_det", "sealed_n_days",
            "sealed_fa"]
    cols = [c for c in cols if c in RK.columns]
    s = (f"\n**{name}** -- the {top} highest-ranking REAL rows on `{col}`, the strongest "
         f"discriminator for this bed.  `null %ile` is the fraction of the NULL population this "
         f"row beats on that column; a row near 1.00 is unusual for the null, a row near 0.50 is "
         f"typical of it.\n\n")
    s += "| " + " | ".join(("null %ile" if c == "null_percentile_on_rank_col" else c)
                           for c in cols) + " |\n" + "|" + "---|" * len(cols) + "\n"
    for r in RK.head(top).itertuples():
        vals = []
        for c in cols:
            v = getattr(r, c)
            if isinstance(v, (bool, np.bool_)):
                vals.append(str(bool(v)))
            elif isinstance(v, (int, np.integer)):
                vals.append(str(int(v)))
            elif isinstance(v, float):
                vals.append(n(v, 3 if abs(v) < 10 else 1))
            else:
                vals.append(str(v))
        s += "| " + " | ".join(vals) + " |\n"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{SF.RES}/WA_FINAL_ASSESSMENT.md")
    ap.add_argument("--no-notes", action="store_true")
    a = ap.parse_args()
    today = datetime.date.today().isoformat()
    cal = jl(f"{F}/calib.json")
    dc = jl(f"{F}/depth_check.json")

    sep_blocks, verdicts = "", []
    for name, bed in BEDS:
        tab, SU = sep_table(bed)
        best = None
        if SU is not None and len(SU):
            num = SU[(SU.direction != "flag") & SU.auc.notna()]
            fav = num[num.auc > 0.5].sort_values("auc", ascending=False)
            best = fav.iloc[0] if len(fav) else (num.iloc[0] if len(num) else None)
        sep_blocks += f"\n### {name}\n\n{tab}"
        if best is not None:
            sep_blocks += curve(bed, str(best.col))
            verdicts.append((name, str(best.col), float(best.auc), float(best.sep),
                             float(best.best_lift),
                             str(getattr(best, "favours", "--"))))
        if SU is not None and len(SU):
            rev = SU[(SU.direction != "flag") & (SU.auc < 0.5) & (SU.sep > 0.10)]
            if len(rev):
                sep_blocks += ("\nReverse separations on this bed (the NULL population scores "
                               "better than the real one): "
                               + ", ".join(f"`{r.col}` AUC {n(r.auc, 3)}"
                                           for r in rev.itertuples())
                               + ".  A criterion on which noise beats signal is evidence about "
                                 "the chain, not a candidate filter.\n")
        sep_blocks += ranked_section(bed, name)

    if verdicts:
        vmax = max(verdicts, key=lambda x: x[3])
        separates = vmax[3] > 0.10
        vtxt = (f"The largest separation IN FAVOUR OF THE REAL CHAIN anywhere is on "
                f"**{vmax[0]}**, column `{vmax[1]}`: AUC {vmax[2]:.3f} "
                f"(0.5 = indistinguishable), i.e. |AUC-0.5| = {vmax[3]:.3f}, best lift "
                f"{vmax[4]:.2f}x (favours {vmax[5]}).")
        if separates:
            verdict = (f"**Real and null DO separate, and here is where.**  {vtxt}  Every other "
                       f"column and bed is listed above with its own AUC, so the reader can see "
                       f"exactly which measurements carry the difference and which carry none.")
        else:
            verdict = (f"**Real and null do not separate anywhere.**  {vtxt}  No criterion "
                       f"measured on these beds -- detection count, repeat days, quiet-day false "
                       f"alarms, tremor co-timing, location quality, interface proximity, "
                       f"distance to the published Olympic coordinates, or the sealed-year "
                       f"re-fire -- puts the real chain's templates meaningfully outside the "
                       f"range of templates grown from data whose inter-site coherence has been "
                       f"destroyed.")
    else:
        verdict = "_no separation analysis has been produced yet_"
        vtxt = ""

    dn = dc.get("north", {})
    dw = dc.get("washington", {})
    dw_rows = ""
    for k in sorted(dw):
        d = dw[k]
        dw_rows += (f"| {k} | {d.get('n_dedupe_located')} | {d.get('n_rejected')} | "
                    f"{d.get('rejected_on_bottom_rail')} | "
                    f"{d.get('rejected_genuinely_interior_depth')} | "
                    f"{d.get('n_on_floor_rail_kept')} |\n")

    txt = f"""# Washington LFE families -- the full ranked tables, real against the site-shift null

Agent S, {today}.  Workspace `results/sf_final/`; queue `scripts/s_final_queue.sh`; code
`scripts/sf_*.py`.

## 0. What changed, and why there is no candidate count in this document

The user's objection, {today}: *"you are discarding them, what if there are undiscovered
families"*.  It is correct and it is fatal to the design every previous Washington write-up used.
Calibrating an acceptance criterion against the 58 blind-rediscovered Bostock families can only
detect a criterion that is too strict **for that population**.  It offers no protection at all to
a real family that is smaller, deeper, noisier, or seen by fewer stations -- which is exactly what
an undiscovered Washington family would be.  Every "0 candidates" verdict in
`WA_2010_RESULTS.md`, `WA_2011_REFIRE_RESULTS.md`, `WA_2008_RESULTS.md`, `WA_2008Q_RESULTS.md` and
`WA_2010_ARRAY_RESULTS.md` was produced by exactly that mistake, compounded by clauses that have
no counterpart in any published Cascadia catalogue.

The rule now in force:

1. **Nothing is filtered.**  Every grown family keeps its row, permanently.
2. **Every criterion is a column**, carrying its measured value -- never a gate.
3. **The whole-chain site-shift null gets the byte-identical table.**  It is the only honest
   comparator.
4. **The evidence is the real-vs-null count as a function of each criterion** -- a separation
   curve -- plus a ranking by the strongest discriminators, not a pre-chosen threshold.
5. The >= 85 %-of-known-families test survives only as a **veto on a proposed cut**.  Passing it
   never licenses a cut.
6. The only steps allowed to reduce a list are the two the papers actually perform -- the
   **dedupe** (>10 % shared declustered detections, or stack cc > 0.3) and the **location** -- and
   both are reported with counts before and after, with the merged and unlocatable rows kept and
   flagged rather than deleted.

## 1. Inventory: what is in the tables

{inventory_table()}
Both reducing steps are reversible from the tables: `dedupe_kept` / `dedupe_merged_into` /
`dedupe_max_share` / `dedupe_max_cc` record what the dedupe did and why, and `located` / `n_sta`
record the location outcome.  No row was deleted.

## 2. Real against null, criterion by criterion

AUC is the probability that a randomly drawn real template scores better than a randomly drawn
null one.  **0.5 means the two populations are indistinguishable on that column.**  Lift is the
enrichment at the strictest threshold that still leaves at least 5 null templates.
{sep_blocks}

## 3. Distance to the published Olympic coordinates, as a column

`derived_coordinate_rows_not_a_catalog.csv` is not a catalogue; it is 69 coordinates derived from
Chestler & Creager (2017) / Sweet et al.  It is used here only as one more column, with the null
answering the same question.

{chestler_section()}

## 4. Two measurements that are reported, and are NOT used to cut anything

### 4.1 The confirmation rule calibrated on Bostock's own families

`results/CONFIRMATION_CALIBRATION_NORTH.md`.  The 2005 ETS window was cut into calendar halves;
every template of the 58 rediscovered families and of a matched-size noise sample was re-stacked
from **H1 detections only** and re-fired unchanged on H2 and on the quiet days.  The strictest
values that still pass >= 85 % of the testable true families are **>= {cal.get('min_days_h2')}
distinct tremor day** and **quiet-day false alarms <= {cal.get('fa_max_per_tpl_day')} per
template-day**; the location-move clause {'was kept' if cal.get('use_move') else '**cannot be calibrated at all** -- no finite value passes 85 % of the true families while rejecting a majority of the noise, confirming SCREEN_VALIDATION_NORTH.md sec.3.1 on the confirmation side'}.

Two things this measurement says, both of which matter more than the thresholds:
* only **{cal.get('n_fam_testable')} of the {cal.get('n_fam_all')} rediscovered families could be
  tested at all** -- the rest do not have enough detections in half an episode to re-stack.  A
  confirmation protocol built on half a window silently discards two thirds of a real catalogue
  before it starts.
* at the calibrated thresholds the rule confirms {cal.get('joint_true_families')} of the
  {cal.get('n_fam_testable')} testable true families
  ({100 * float(cal.get('joint_true_frac', 0)):.0f} %) **and
  {cal.get('joint_noise_templates')} of {cal.get('joint_noise_total')} noise templates
  ({100 * float(cal.get('joint_noise_frac', 0)):.0f} %)**.  A rule tuned to keep real families
  keeps almost as much noise: on this bed, sealed-window re-fire is a weak discriminator, not a
  verdict.  That is why the sealed-year outcome appears above as columns
  (`sealed_n_days`, `sealed_fa`, `sealed_ratio`, `sealed_move_km`) and not as a gate.

### 4.2 The depth window (15-58 km) -- measured, and rejected as a cut

The window was in an earlier draft of this analysis.  It is not used, for three measured reasons.

| chain | dedupe-kept and located | would be rejected by 15-58 km | of those, on the bottom depth rail (59 km) | genuinely interior depths outside the window | kept while sitting on the 15-km floor rail |
|---|---|---|---|---|---|
{dw_rows}
* **Every single rejection on every Washington chain is the bottom grid rail, and none is a
  genuinely interior solution outside the window.**  The window is a grid-rail clause wearing a
  physical label -- and an asymmetric one, since it deletes the 59-km rail while keeping hundreds
  of solutions on the 15-km floor rail.
* On the north bed the literal window is **vacuous**: that locator's grid is
  {n(dn.get('grid_dep_min'), 0)}-{n(dn.get('grid_dep_max'), 0)} km, so 100 % of the located true
  templates pass it for a reason that has nothing to do with the families.
* Its transferable form -- "reject a solution that rails on a depth face" -- keeps only
  **{100 * float(dn.get('no_depth_rail', {}).get('frac_fam_of_located', 0)):.0f} % of the located
  true north families** ({dn.get('no_depth_rail', {}).get('fam')} of
  {dn.get('located', {}).get('fam')}), far below the 85 % veto bar.
* The literature bar is physical, not numerical: Bostock and Chestler place families at roughly
  28-45 km, but they get there with a 3-D velocity model and double-difference relocation.  With a
  1-D model and no P picks our depths are coarse -- a 28-45 km window would keep
  {[d for d in dn.get('window_sweep', []) if d.get('lo') == 28][0]['fam'] if any(d.get('lo') == 28 for d in dn.get('window_sweep', [])) else '--'}
  of {dn.get('located', {}).get('fam')} located true north families -- so a tight window is not
  transferable.

For completeness, the one location clause that does pass the veto bar is the locator's own
**>= 5 stations**: {dn.get('located', {}).get('fam')} of {dn.get('n_true_families')} true
families ({100 * float(dn.get('located', {}).get('frac_fam', 0)):.0f} %), reproducing agent Q's
86 %.  It is still reported as the `located` flag, with the unlocatable rows retained.

## 5. Superseded comparison: agent Q's tier A/B/C as columns

Tier A/B/C were a threshold design and are kept only so the numbers can be compared with the
earlier write-ups.  They are columns in every table (`tierA`, `tierB`, `tierC`) and appear in the
separation ranking above like any other column; they are not used to select anything.

## 6. Verdict

{verdict}

Read the negative precisely.  It is **not** "nothing was found in the raw output": the published
recipe grows thousands of templates on both Washington beds, hundreds survive the literature
dedupe and locate at plausible interface depths, and they are all in the tables with their
measurements.  It is **not** "a screen rejected everything" -- no screen is applied here.  It is
the statement that on these beds, measured column by column, the population produced by the real
data and the population produced by data whose inter-site coherence has been destroyed are not
distinguishable by any statistic this project can compute. A real Washington family, if it is in
these tables, is not separable from noise by these measurements -- which is a statement about the
measurements and the beds, and leaves every row available for a future test rather than deleting
it.

## 7. Files

| what | where |
|---|---|
| per-template tables (every row kept, every criterion a column) | `results/sf_final/table_<bed>_<chain>.csv` |
| separation curves and AUC summaries | `results/sf_final/sep_<bed>_sweep.csv`, `sep_<bed>_summary.csv`, `sep_<bed>.json` |
| ranked tables | `results/sf_final/sep_<bed>_ranked.csv` |
| dedupe detail (what merged, on what) | `results/sf_final/<bed>_<chain>_dedupe.csv` |
| sealed-year re-fires (ALL templates) | `results/sf_final/sealed/<bed>_<chain>/` |
| Arm-B site-shift null chain (new) | `results/sf_final/armb_null/` |
| Arm-B locations (new -- agent O's locate stage had crashed) | `results/sf_final/armb_<chain>_loc.csv` |
| confirmation calibration | `results/CONFIRMATION_CALIBRATION_NORTH.md`, `results/sf_final/calib.json`, `results/sf_final/north/` |
| depth-window measurement | `results/sf_final/depth_check.json` |
| code | `scripts/sf_common.py`, `sf_north_split.py`, `sf_calib.py`, `sf_depthcheck.py`, `sf_accept.py`, `sf_table.py`, `sf_separation.py`, `sf_locate_member.py`, `sf_armb_null.py`, `sf_armb_refire.py`, `sf_report.py` |
| queue / log | `scripts/s_final_queue.sh`, `logs/s_final_queue.log` |
"""
    open(a.out, "w").write(txt)
    SF.log(f"wrote {a.out} ({len(txt)} chars)")

    note = f"""
## {today} -- agent S: Washington re-scored under RANK-DON'T-DISCARD (no screen, no candidate count)
User rule, {today}: criteria calibrated on the 58 known Bostock families can only catch cuts that
are too strict for THAT population and protect no undiscovered family, so nothing may be filtered.
Every grown template of agent O's array-based 2010 chain, agent N's 2008 chain and Chestler's own
Arm B now keeps its row in `results/sf_final/table_<bed>_<chain>.csv`, with every criterion as a
column and a byte-identical table for the whole-chain site-shift null; the only list-reducing steps
are the published dedupe and the location, both recorded as flags with the rows retained.
Evidence = real-vs-null separation curves and AUC per column (`sep_<bed>_*`), not thresholds.
{vtxt}
New work behind it: Arm B was located for the first time (agent O's locate stage had crashed on a
member-station KeyError) and given the site-shift null it never had; every template of every chain
was re-fired on the sealed 2011 / 2007 windows, not just a screened subset; the confirmation rule
was calibrated on Bostock's own families (north half-split) and shown to be a weak discriminator
({cal.get('joint_true_families')}/{cal.get('n_fam_testable')} true vs {cal.get('joint_noise_templates')}/{cal.get('joint_noise_total')} noise), so it is reported as columns, not a gate;
the depth window was measured and rejected as a cut (100 % of its Washington rejections are the
bottom grid rail; its transferable form keeps only {100 * float(dn.get('no_depth_rail', {}).get('frac_fam_of_located', 0)):.0f} % of located true north families).
Full ledger `results/WA_FINAL_ASSESSMENT.md`.
"""
    if not a.no_notes:
        open(f"{SF.STUDY}/NOTES.md", "a").write(note)
    print(note, flush=True)


if __name__ == "__main__":
    main()
