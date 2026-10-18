# Calibrating the CONFIRMATION rule on Bostock's own families (north half-split)

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
20 usable literature-quiet 2005 days.  The H1-only stack and a stack built from that template's
own H2 detections are each located on `k_locate`'s grid, and the distance between them is the
half-to-half location move.  This is a Washington confirmation, run where truth exists.

**Populations.**  TRUE = the 177 templates of the **58 blind-rediscovered kept
Bostock families** (`results/repro/classify15.csv`, cls = kept_family).  NOISE = a matched-size
random sample (177) of the 1,201 unlocatable/unmatched templates, seed 20260923.

**Testability.**  37 of the 177 true templates and
84 of the 177 noise templates had >= 3 H1 detections and produced an
H1-only stack, covering **18 of the 58 true families**.  A family with no
H1 stack cannot be tested here and is excluded from the denominators below; that exclusion is
itself a result -- 40 of Bostock's own families would be lost by a
confirmation protocol that needs half the episode to build the template.

---

## 1. What the re-fire looks like

| | TRUE templates (n = 37) | NOISE templates (n = 84) |
|---|---|---|
| H2 detections, median (q25-q75) | 10 (2-21) | 44 (9-114) |
| H2 distinct days, median | 3 | 9 |
| quiet-day FA per template-day, median | 0.05 | 0.28 |
| quiet-day FA, 85th percentile | 0.49 | 4.70 |
| H2-rate / quiet-rate ratio, median | 8.33 | 8.92 |
| H1-to-H2 move, 3-D, median km | 26.3 | 42.7 |
| H1-to-H2 move, epicentral, median km | 12.9 | 41.9 |

Network rates: the H1 templates fire 5.67
per template-day on the H2 tremor days against
2.31 on the quiet days.

## 2. Clause (i) -- minimum distinct tremor days in the re-fire window

| >= days | true families | of | fraction | noise templates | noise fraction |
|---|---|---|---|---|---|
| 1 | 18 | 18 | 1.0 | 81 | 0.964 |
| 2 | 15 | 18 | 0.833 | 74 | 0.881 |
| 3 | 14 | 18 | 0.778 | 73 | 0.869 |
| 4 | 11 | 18 | 0.611 | 68 | 0.81 |
| 5 | 10 | 18 | 0.556 | 64 | 0.762 |
| 6 | 10 | 18 | 0.556 | 63 | 0.75 |
| 7 | 6 | 18 | 0.333 | 60 | 0.714 |
| 8 | 4 | 18 | 0.222 | 56 | 0.667 |
| 9 | 4 | 18 | 0.222 | 45 | 0.536 |
| 10 | 3 | 18 | 0.167 | 34 | 0.405 |
| 11 | 2 | 18 | 0.111 | 22 | 0.262 |
| 12 | 2 | 18 | 0.111 | 12 | 0.143 |


**Chosen: >= 1.0 distinct tremor days** (the strictest value passing >= 85 % of the
testable true families).

## 3. Clause (ii) -- maximum quiet-day false-alarm rate

| FA <= /tpl-day | true families | of | fraction | noise templates | noise fraction |
|---|---|---|---|---|---|
| 0.0 | 7 | 18 | 0.389 | 12 | 0.143 |
| 0.05 | 11 | 18 | 0.611 | 20 | 0.238 |
| 0.1 | 13 | 18 | 0.722 | 28 | 0.333 |
| 0.25 | 15 | 18 | 0.833 | 42 | 0.5 |
| 0.5 | 15 | 18 | 0.833 | 55 | 0.655 |
| 1.0 | 16 | 18 | 0.889 | 63 | 0.75 |
| 1.5 | 16 | 18 | 0.889 | 64 | 0.762 |
| 2.0 | 17 | 18 | 0.944 | 68 | 0.81 |
| 3.0 | 17 | 18 | 0.944 | 70 | 0.833 |
| 5.0 | 17 | 18 | 0.944 | 71 | 0.845 |
| 10.0 | 17 | 18 | 0.944 | 76 | 0.905 |
| no cap | 18 | 18 | 1.0 | 84 | 1.0 |


**Chosen: <= 1.0 detections per template-day on the quiet days.**  For scale, the
0.1/template-day cap used in every Washington confirmation so far passes
72% of the true families, and Thomas & Bostock
(2015) accept 2-3 false detections per day for their own high-quality templates.

Diagnostic (not part of the rule -- no published counterpart): the H2-rate/quiet-rate ratio.

| >= ratio | true families | of | fraction | noise templates | noise fraction |
|---|---|---|---|---|---|
| 0.5 | 16 | 18 | 0.889 | 78 | 0.929 |
| 1.0 | 15 | 18 | 0.833 | 72 | 0.857 |
| 1.5 | 15 | 18 | 0.833 | 67 | 0.798 |
| 2.0 | 15 | 18 | 0.833 | 63 | 0.75 |
| 3.0 | 15 | 18 | 0.833 | 60 | 0.714 |
| 5.0 | 14 | 18 | 0.778 | 53 | 0.631 |


## 4. Clause (iii) -- maximum H1-to-H2 location move

3-D:

| move <= km | true families | of | fraction | noise templates | noise fraction |
|---|---|---|---|---|---|
| 5.0 | 2 | 18 | 0.111 | 0 | 0.0 |
| 10.0 | 3 | 18 | 0.167 | 1 | 0.012 |
| 15.0 | 3 | 18 | 0.167 | 1 | 0.012 |
| 20.0 | 3 | 18 | 0.167 | 2 | 0.024 |
| 30.0 | 5 | 18 | 0.278 | 5 | 0.06 |
| 50.0 | 6 | 18 | 0.333 | 7 | 0.083 |
| 80.0 | 9 | 18 | 0.5 | 12 | 0.143 |
| 120.0 | 9 | 18 | 0.5 | 12 | 0.143 |
| any | 9 | 18 | 0.5 | 13 | 0.155 |


epicentral:

| move <= km | true families | of | fraction | noise templates | noise fraction |
|---|---|---|---|---|---|
| 5.0 | 2 | 18 | 0.111 | 0 | 0.0 |
| 10.0 | 4 | 18 | 0.222 | 1 | 0.012 |
| 15.0 | 4 | 18 | 0.222 | 1 | 0.012 |
| 20.0 | 4 | 18 | 0.222 | 2 | 0.024 |
| 30.0 | 5 | 18 | 0.278 | 5 | 0.06 |
| 50.0 | 6 | 18 | 0.333 | 8 | 0.095 |
| 80.0 | 9 | 18 | 0.5 | 12 | 0.143 |
| 120.0 | 9 | 18 | 0.5 | 12 | 0.143 |
| any | 9 | 18 | 0.5 | 13 | 0.155 |


**Verdict on clause (iii): DROPPED.**
No finite move threshold both passes 85 % of the true families and rejects a majority of the noise sample: the value that reaches 85 % of true families also admits most of the noise, so the clause carries no information. This confirms SCREEN_VALIDATION_NORTH.md sec.3.1 on the confirmation side: location stability under re-stacking is not a property Bostock's families have on this bed, and the 10-km clause used in every Washington confirmation so far was rejecting real families for a property they do not possess.

## 5. THE CALIBRATED CONFIRMATION RULE

> A template is **CONFIRMED** on a sealed window if it fires on **>= 1.0 distinct tremor
> days** of that window **and** its quiet-day false-alarm rate there is **<= 1.0 detections
> per template-day**.

Measured on this bed: it confirms **16 of the 18 testable true families
(89%)** and **60 of the 84 re-fired noise
templates (71%)**.

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
