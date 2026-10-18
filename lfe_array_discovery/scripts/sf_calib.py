#!/usr/bin/env python
"""STAGE 1 analysis -- set the CONFIRMATION rule from the north half-split, and write
results/CONFIRMATION_CALIBRATION_NORTH.md + results/sf_final/calib.json.

Rule fixed before the numbers (results/SCREEN_VALIDATION_NORTH.md sec.3's pre-stated rule): a
confirmation clause is set at the value that passes >= 85 % of the TRUE Bostock families that
could be tested at all; a clause with no discriminating power against the noise sample is
reported as such; a clause that cannot reach 85 % at any usable value is DROPPED.

Unit of account = the family (a family passes a clause if any of its templates does), the same
convention under which the 58 count as rediscovered.  Template-level fractions are reported beside
it.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

OUT = f"{SF.RES}/sf_final/north"
PASS_FRAC = 0.85
DAY_GRID = list(range(1, 13))
FA_GRID = [0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 1e9]
MOVE_GRID = [5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 120.0, 1e9]


def build():
    S = pd.read_csv(f"{OUT}/set.csv")
    H1 = pd.read_csv(f"{OUT}/h1_stackinfo.csv").set_index("grp")
    qs = json.load(open(f"{OUT}/det_quiet_summary.json"))
    hs = json.load(open(f"{OUT}/det_h2_summary.json"))
    DH = pd.read_parquet(f"{OUT}/det_h2.parquet")
    DQ = pd.read_parquet(f"{OUT}/det_quiet.parquet")
    L1 = pd.read_csv(f"{OUT}/loc_h1.csv").set_index("grp")
    L2 = pd.read_csv(f"{OUT}/loc_h2.csv").set_index("grp")
    tplh1 = set(np.load(f"{OUT}/tpl_h1.npz", allow_pickle=True).files)
    have_h1 = {int(k.split("_", 1)[0]) for k in tplh1 if k != "meta"}
    nq, nh = max(qs["n_dates"], 1), max(hs["n_dates"], 1)
    nd_h2 = DH.groupby("grp").size()
    ndays_h2 = DH.groupby("grp").date.nunique()
    nd_q = DQ.groupby("grp").size()
    rows = []
    for r in S.itertuples():
        g = int(r.grp)
        rec = dict(grp=g, pop=r.pop, bkey=int(r.bkey), n_det_disc=int(r.n_det),
                   n_det_h1=int(H1.n_det_h1.get(g, 0)),
                   refired=bool(g in have_h1),
                   n_det_h2=int(nd_h2.get(g, 0)), n_days_h2=int(ndays_h2.get(g, 0)),
                   n_det_quiet=int(nd_q.get(g, 0)),
                   fa_quiet=float(nd_q.get(g, 0)) / nq,
                   rate_h2=float(nd_h2.get(g, 0)) / nh)
        rec["ratio_h2_quiet"] = (rec["rate_h2"] / rec["fa_quiet"]) if rec["fa_quiet"] > 0 else (
            np.inf if rec["rate_h2"] > 0 else np.nan)
        for tag, L in (("h1", L1), ("h2", L2)):
            if g in L.index and bool(L.located.get(g, False)):
                rec[f"lat_{tag}"] = float(L.lat[g])
                rec[f"lon_{tag}"] = float(L.lon[g])
                rec[f"dep_{tag}"] = float(L.dep[g])
                rec[f"loc_{tag}"] = True
            else:
                rec[f"lat_{tag}"] = rec[f"lon_{tag}"] = rec[f"dep_{tag}"] = np.nan
                rec[f"loc_{tag}"] = False
        if rec["loc_h1"] and rec["loc_h2"]:
            rec["move_epi_km"] = float(SF.epi_km(rec["lat_h1"], rec["lon_h1"],
                                                 rec["lat_h2"], rec["lon_h2"]))
            rec["move_3d_km"] = float(SF.hyp_km(rec["lat_h1"], rec["lon_h1"], rec["dep_h1"],
                                                rec["lat_h2"], rec["lon_h2"], rec["dep_h2"]))
        else:
            rec["move_epi_km"] = rec["move_3d_km"] = np.nan
        rows.append(rec)
    T = pd.DataFrame(rows)
    T.to_csv(f"{OUT}/calib_templates.csv", index=False)
    return T, nq, nh


def fam_pass(T, mask):
    """fraction of TRUE families with >= 1 template passing `mask`, over the testable families."""
    t = T[T["pop"] == "true"]
    test = t[t.refired]
    fams = sorted(set(test.bkey))
    ok = sorted(set(t[mask.reindex(t.index, fill_value=False) & t.refired].bkey))
    return len(ok), len(fams), sorted(set(t.bkey))


def sweep(T, col, grid, direction):
    """direction 'ge' = keep value >= thr (a minimum clause); 'le' = keep value <= thr."""
    out = []
    t = T[(T["pop"] == "true") & T.refired]
    z = T[(T["pop"] == "noise") & T.refired]
    nf = len(set(t.bkey))
    for thr in grid:
        m = (T[col] >= thr) if direction == "ge" else (T[col] <= thr)
        m = m & T.refired
        ft = len(set(t[m.reindex(t.index, fill_value=False)].bkey))
        out.append(dict(thr=thr, fam_pass=ft, fam_total=nf,
                        frac_fam=ft / max(nf, 1),
                        tpl_pass=int(m[T["pop"] == "true"].sum()),
                        tpl_total=int(len(t)),
                        noise_pass=int(m[T["pop"] == "noise"].sum()),
                        noise_total=int(len(z)),
                        frac_noise=float(m[T["pop"] == "noise"].sum()) / max(len(z), 1)))
    return pd.DataFrame(out)


def pick(SW, direction):
    """the strictest threshold that still passes >= 85 % of testable true families."""
    ok = SW[SW.frac_fam >= PASS_FRAC]
    if not len(ok):
        return None
    return float(ok.thr.max() if direction == "ge" else ok.thr.min())


def md_table(df, cols, head):
    s = "| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n"
    for r in df.itertuples():
        s += "| " + " | ".join(str(getattr(r, c)) for c in cols) + " |\n"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{SF.RES}/CONFIRMATION_CALIBRATION_NORTH.md")
    a = ap.parse_args()
    T, nq, nh = build()
    t = T[T["pop"] == "true"]
    z = T[T["pop"] == "noise"]
    fam_all = sorted(set(t.bkey))
    fam_test = sorted(set(t[t.refired].bkey))

    SD = sweep(T, "n_days_h2", DAY_GRID, "ge")
    SF_ = sweep(T, "fa_quiet", FA_GRID, "le")
    SM3 = sweep(T, "move_3d_km", MOVE_GRID, "le")
    SME = sweep(T, "move_epi_km", MOVE_GRID, "le")
    SR = sweep(T, "ratio_h2_quiet", [0.5, 1.0, 1.5, 2.0, 3.0, 5.0], "ge")

    min_days = pick(SD, "ge")
    fa_max = pick(SF_, "le")
    move3 = pick(SM3, "le")
    movee = pick(SME, "le")
    # a location-move clause is kept only if some finite value both passes 85 % of the true
    # families AND rejects a real share of the noise sample
    def useful(SW, thr):
        if thr is None or not np.isfinite(thr) or thr >= 1e8:
            return False
        row = SW[SW.thr == thr].iloc[0]
        return bool(row.frac_noise <= 0.5)
    use_move = useful(SM3, move3) or useful(SME, movee)

    rule = dict(min_days_h2=(int(min_days) if min_days is not None else None),
                fa_max_per_tpl_day=(float(fa_max) if fa_max is not None else None),
                move_km_3d=(float(move3) if move3 is not None else None),
                move_km_epi=(float(movee) if movee is not None else None),
                use_move=bool(use_move),
                pass_frac_target=PASS_FRAC,
                n_quiet_days=nq, n_h2_days=nh,
                n_fam_all=len(fam_all), n_fam_testable=len(fam_test))

    m = pd.Series(True, index=T.index) & T.refired
    if min_days is not None:
        m &= T.n_days_h2 >= min_days
    if fa_max is not None and np.isfinite(fa_max) and fa_max < 1e8:
        m &= T.fa_quiet <= fa_max
    if use_move:
        thr = move3 if useful(SM3, move3) else movee
        col = "move_3d_km" if useful(SM3, move3) else "move_epi_km"
        m &= T[col] <= thr
        rule["move_clause"] = dict(col=col, thr=float(thr))
    joint_fam = len(set(t[m.reindex(t.index, fill_value=False)].bkey))
    joint_noise = int(m[T["pop"] == "noise"].sum())
    rule.update(joint_true_families=joint_fam, joint_true_frac=joint_fam / max(len(fam_test), 1),
                joint_noise_templates=joint_noise,
                joint_noise_total=int(len(z[z.refired])),
                joint_noise_frac=joint_noise / max(len(z[z.refired]), 1))
    os.makedirs(f"{SF.RES}/sf_final", exist_ok=True)
    json.dump(rule, open(f"{SF.RES}/sf_final/calib.json", "w"), indent=1, default=float)
    T["confirmed_by_rule"] = m
    T.to_csv(f"{OUT}/calib_templates.csv", index=False)

    q = lambda s, p: float(np.nanquantile(s.astype(float), p)) if len(s.dropna()) else float("nan")
    tt = t[t.refired]
    zz = z[z.refired]
    txt = f"""# Calibrating the CONFIRMATION rule on Bostock's own families (north half-split)

Agent S, 2026-09-23.  Workspace `results/sf_final/north/`; code `scripts/sf_north_split.py`,
`scripts/sf_calib.py`; queue `scripts/s_final_queue.sh`.

**Why.**  Every Washington "confirmation" so far (`WA_2011_REFIRE_RESULTS.md` sec.5,
`WA_2008_RESULTS.md` sec.5, `WA_2010_ARRAY_RESULTS.md` sec.6) used a rule nobody had ever measured
against a known catalogue: >= 3 tremor days in the sealed year, quiet-day false alarms
<= 0.1/template-day, and an INTERIOR location within 10 km (3-D) of the discovery location.
`SCREEN_VALIDATION_NORTH.md` sec.3.1 already showed the location leg is hopeless (a true family's
re-stacked location moves a median 35 km), but the rule as a whole was never calibrated.  This is
that measurement.

**The experiment.**  The 24-day 2005 POLARIS ETS window is cut into calendar halves --
**H1 = 2005-08-29 ... 09-09**, **H2 = 2005-09-10 ... 09-21**.  Every template in the two
populations is re-stacked (`k_grow.build_templates`, verbatim) **from its H1 detections only**,
then fired unchanged at **8 x MAD** (`k_grow.matched_filter`, verbatim) on H2 and on the
{nq} usable literature-quiet 2005 days.  The H1-only stack and a stack built from that template's
own H2 detections are each located on `k_locate`'s grid, and the distance between them is the
half-to-half location move.  This is a Washington confirmation, run where truth exists.

**Populations.**  TRUE = the {len(t)} templates of the **{len(fam_all)} blind-rediscovered kept
Bostock families** (`results/repro/classify15.csv`, cls = kept_family).  NOISE = a matched-size
random sample ({len(z)}) of the 1,201 unlocatable/unmatched templates, seed 20260923.

**Testability.**  {int(tt.refired.sum())} of the {len(t)} true templates and
{int(zz.refired.sum())} of the {len(z)} noise templates had >= 3 H1 detections and produced an
H1-only stack, covering **{len(fam_test)} of the {len(fam_all)} true families**.  A family with no
H1 stack cannot be tested here and is excluded from the denominators below; that exclusion is
itself a result -- {len(fam_all) - len(fam_test)} of Bostock's own families would be lost by a
confirmation protocol that needs half the episode to build the template.

---

## 1. What the re-fire looks like

| | TRUE templates (n = {int(tt.refired.sum())}) | NOISE templates (n = {int(zz.refired.sum())}) |
|---|---|---|
| H2 detections, median (q25-q75) | {q(tt.n_det_h2, .5):.0f} ({q(tt.n_det_h2, .25):.0f}-{q(tt.n_det_h2, .75):.0f}) | {q(zz.n_det_h2, .5):.0f} ({q(zz.n_det_h2, .25):.0f}-{q(zz.n_det_h2, .75):.0f}) |
| H2 distinct days, median | {q(tt.n_days_h2, .5):.0f} | {q(zz.n_days_h2, .5):.0f} |
| quiet-day FA per template-day, median | {q(tt.fa_quiet, .5):.2f} | {q(zz.fa_quiet, .5):.2f} |
| quiet-day FA, 85th percentile | {q(tt.fa_quiet, .85):.2f} | {q(zz.fa_quiet, .85):.2f} |
| H2-rate / quiet-rate ratio, median | {q(tt.ratio_h2_quiet.replace(np.inf, np.nan), .5):.2f} | {q(zz.ratio_h2_quiet.replace(np.inf, np.nan), .5):.2f} |
| H1-to-H2 move, 3-D, median km | {q(tt.move_3d_km, .5):.1f} | {q(zz.move_3d_km, .5):.1f} |
| H1-to-H2 move, epicentral, median km | {q(tt.move_epi_km, .5):.1f} | {q(zz.move_epi_km, .5):.1f} |

Network rates: the H1 templates fire {json.load(open(f'{OUT}/det_h2_summary.json'))['per_tpl_day']:.2f}
per template-day on the H2 tremor days against
{json.load(open(f'{OUT}/det_quiet_summary.json'))['per_tpl_day']:.2f} on the quiet days.

## 2. Clause (i) -- minimum distinct tremor days in the re-fire window

{md_table(SD.assign(frac_fam=SD.frac_fam.round(3), frac_noise=SD.frac_noise.round(3)),
          ["thr", "fam_pass", "fam_total", "frac_fam", "noise_pass", "frac_noise"],
          [">= days", "true families", "of", "fraction", "noise templates", "noise fraction"])}

**Chosen: >= {min_days} distinct tremor days** (the strictest value passing >= 85 % of the
testable true families).

## 3. Clause (ii) -- maximum quiet-day false-alarm rate

{md_table(SF_.assign(thr=[("no cap" if x >= 1e8 else x) for x in SF_.thr],
                     frac_fam=SF_.frac_fam.round(3), frac_noise=SF_.frac_noise.round(3)),
          ["thr", "fam_pass", "fam_total", "frac_fam", "noise_pass", "frac_noise"],
          ["FA <= /tpl-day", "true families", "of", "fraction", "noise templates",
           "noise fraction"])}

**Chosen: <= {fa_max} detections per template-day on the quiet days.**  For scale, the
0.1/template-day cap used in every Washington confirmation so far passes
{float(SF_[SF_.thr == 0.1].frac_fam.iloc[0]):.0%} of the true families, and Thomas & Bostock
(2015) accept 2-3 false detections per day for their own high-quality templates.

Diagnostic (not part of the rule -- no published counterpart): the H2-rate/quiet-rate ratio.

{md_table(SR.assign(frac_fam=SR.frac_fam.round(3), frac_noise=SR.frac_noise.round(3)),
          ["thr", "fam_pass", "fam_total", "frac_fam", "noise_pass", "frac_noise"],
          [">= ratio", "true families", "of", "fraction", "noise templates", "noise fraction"])}

## 4. Clause (iii) -- maximum H1-to-H2 location move

3-D:

{md_table(SM3.assign(thr=[("any" if x >= 1e8 else x) for x in SM3.thr],
                     frac_fam=SM3.frac_fam.round(3), frac_noise=SM3.frac_noise.round(3)),
          ["thr", "fam_pass", "fam_total", "frac_fam", "noise_pass", "frac_noise"],
          ["move <= km", "true families", "of", "fraction", "noise templates", "noise fraction"])}

epicentral:

{md_table(SME.assign(thr=[("any" if x >= 1e8 else x) for x in SME.thr],
                     frac_fam=SME.frac_fam.round(3), frac_noise=SME.frac_noise.round(3)),
          ["thr", "fam_pass", "fam_total", "frac_fam", "noise_pass", "frac_noise"],
          ["move <= km", "true families", "of", "fraction", "noise templates", "noise fraction"])}

**Verdict on clause (iii): {'KEPT at ' + str(rule.get('move_clause')) if use_move else 'DROPPED'}.**
{'' if use_move else "No finite move threshold both passes 85 % of the true families and rejects a majority of the noise sample: the value that reaches 85 % of true families also admits most of the noise, so the clause carries no information. This confirms SCREEN_VALIDATION_NORTH.md sec.3.1 on the confirmation side: location stability under re-stacking is not a property Bostock's families have on this bed, and the 10-km clause used in every Washington confirmation so far was rejecting real families for a property they do not possess."}

## 5. THE CALIBRATED CONFIRMATION RULE

> A template is **CONFIRMED** on a sealed window if it fires on **>= {min_days} distinct tremor
> days** of that window **and** its quiet-day false-alarm rate there is **<= {fa_max} detections
> per template-day**{'' if not use_move else ', and its re-located position moves by no more than ' + str(rule.get('move_clause'))}.

Measured on this bed: it confirms **{joint_fam} of the {len(fam_test)} testable true families
({joint_fam / max(len(fam_test), 1):.0%})** and **{joint_noise} of the {len(zz)} re-fired noise
templates ({joint_noise / max(len(zz), 1):.0%})**.

Caveats.  (1) The quiet days carry 7-12 cached stations against 9-15 on the ETS days, so the
quiet-day rate is, if anything, over-estimated on both populations.  (2) The noise population is
not certified noise (SCREEN_VALIDATION_NORTH.md sec.6), so its confirmation fraction is an upper
bound on the false-confirmation rate.  (3) The thresholds are rates per template-day and day
counts, so they transfer between beds with different day counts; the FA cap does NOT transfer to a
bed with a different detector threshold, and Arm B (45 x MAD on a 2nd-root stack) is reported
separately for that reason.

## 6. Files

| what | where |
|---|---|
| per-template measurements | `results/sf_final/north/calib_templates.csv` |
| H1-only templates / H2-only templates | `results/sf_final/north/tpl_h1.npz`, `tpl_h2.npz` |
| H2 and quiet-day re-fire detections | `results/sf_final/north/det_h2.parquet`, `det_quiet.parquet` |
| H1 / H2 locations | `results/sf_final/north/loc_h1.csv`, `loc_h2.csv` |
| the rule, machine-readable | `results/sf_final/calib.json` |
"""
    open(a.out, "w").write(txt)
    SF.log(f"wrote {a.out} ({len(txt)} chars); rule = {json.dumps(rule, default=float)}")


if __name__ == "__main__":
    main()
