# SVI backward extension 2003–2013 — results in plain words

Written Sun 2026-10-11, ~12:10 am PT, after P6 and P7. Details and every number: `p6/P6.md`, `p7/P7.md` (rules written
before computing, results below them). Earlier stages: `PLAN.md`, `p4b/P4B.md`, `../backward_gates/GATES.md`.
"Scrambled control" = the same detector on the same day with each station site's data shifted to another time of that day.

## What was made

- **Catalogue v1.4-backward**, `data/catalogue/svi_catalogue_v1_4_backward.parquet`: 886,249 candidates
  2003-01-01 → 2013-12-31 (614,591 data, 271,658 scrambled control), with the same columns as the forward SVI catalogue
  (`refine/conf/P_by_candidate.parquet`) plus `era`, `net_set`, `tpl_epoch`, `detector`, `P`, `burst_day`, `n_z2`,
  `multi_site_z2`, `counted_once`, `primary`, `ctl_heavy_set`.
- **Rate tables** in the forward format, `data/catalogue/svi_rates_v1_4_backward_*` (per family and day, monthly per family,
  monthly total), each with an uncorrected and a completeness-corrected version. The `primary == True` rows give one series
  per day that joins the forward 2014–2026 tables.
- Confidence P for every candidate. B0 (CN4), B1 and B2 were refitted per station set on the 2011–12 truth days. B3 and B4
  use the frozen forward model unchanged. The era detector's model was fitted on the 2004 ETS build blocks only.

## Passed

- **P exists everywhere, and for the frozen network it behaves like forward.** The scrambled control gets P > 0.5 in
  0.8–2.2 % of candidates for CN4, B1, B2, B3 and B4, about the forward level. The CN4/B1 worry from G1 (5–6 %) did not
  carry over to the real 2003–05 data, so no refit was triggered.
- **The confidence P ranks real events well for the frozen network.** On the new truth windows AUC is 0.84–0.87 for
  B1, B2 and B3, against 0.71–0.74 for the raw detection strength. CN4 in 2003 reaches 0.78 against 0.69.
- **The high-P candidates between episodes look like real LFEs.** On quiet days the two sites that hardly triggered a
  candidate still see S energy at its predicted time. The scrambled control shows no such energy. This holds for CN4,
  B2, B3, B4 and the era detector. It is the same test that passed forward.
- **The quiet-day flag used before 2010 agrees with PNSN.** On 2010–13, between-episode rates on env-quiet days divided
  by rates on PNSN-quiet days give 0.97–1.11 per year and 1.04 pooled.
- **No step at the join to the forward catalogue.** The daily excess 60 d before 2014-01-01 vs 60 d after gives
  ratio 0.62, z −1.6. Activity alone flags 35 % of month starts in this test (see "Failed").
- **The era detector finds much more of Bostock's 2003–2005 catalogue.** At P ≥ 0.75 it recovers 14 % (2003) and 21 %
  (2005) of his events, with 79–84 % precision. The frozen CN4/B1 network recovers 1 %.
- **Yearly Bostock-equivalent totals are the right size.** The total is ΣP excess × each set's Bostock k, which includes
  between-episode activity. It comes to 34k (2003), 60k (2004), 52k (2005), 24k (2008), 17k (2009) and 22k / 21k /
  37k / 20k (2010–13). His catalogue has 31k, 55k, 55k, 14k, 22k and 18k / 18k / 27k / 14k in the ETS windows.

## Failed

- **SCOPE §5 bars on the one-look windows.** No set passes all three bars:
  - The control bar (≤ 1 % P > 0.5) fails everywhere, at 2.3–17 %. The forward model also missed it in 2013.
  - Matched/ΣP fails for CN4 2003 (1.30), B2 2007 (0.54) and B3 2008 (0.77).
  - By G5, every era is therefore reported as a ΣP bracket only.
- **The era detector's P is not calibrated outside big ETS.**
  - Its scrambled control gets P > 0.5 in 13.5 % of candidates. It was fitted on the dense 2004 ETS, where 76 % of
    candidates are real.
  - On the 2003 and 2005 windows its P ranks events no better than the raw detection strength (AUC 0.74–0.77 against
    0.75–0.76).
  - Its rows are flagged `ctl_heavy_set`. For single repeats use P ≥ 0.90 (P\*), not 0.75. Its rates are ΣP excess only.
- **Step checks (b) fail, and the test cannot separate network steps from activity.**
  - Flagged boundaries:
    - CN4 → B1 (z −2.2);
    - B1 → B2 (−2.6, right after the 2005 ETS);
    - LRIV's first data on 2009-11-20 (−3.3);
    - the era → frozen switch on 2006-07-01 (−7.8);
    - 2009-08-27, where no station changed (+3.7).
  - With no change at all, the same test flags 35 % of month starts on between-episode days.
- **Control fairness (c) fails for every frozen set.** On quiet days the excess rises with the control (slope 0.4–1.5).
  The scrambled control removes only part of what is not LFE. The era detector passes on all candidates.
- **Recall vs ladder (d) fails, except B2 2007 and B3 2009.** The era/frozen recall ratio on the same windows is
  15–17×. The ETS ladder factor is 8–10×.
- **Station-set invariance fails** for CN4, B1 and B2 on the 2014–17 quiet days (median monthly ratio 1.7–1.9). B3
  passes (1.04). So, as forward, every era is "control-subtracted bracket", never "rate" (G6).

## What the catalogue can be used for

- **ETS episodes 2003–2013: yes.**
  - Use the `primary` detector:
    - era for 2003-02 → 2006-06, at P ≥ 0.90;
    - frozen otherwise, at P ≥ 0.75.
  - Count with `counted_once`. Per family, use `rate` (Bostock-equivalent) for episode sizes.
  - On the same 2003–2005 ETS days the era detector sees 8–10× more ΣP excess than the frozen network. Episode sizes
    across the 2006-07 switch are only comparable through `rate`, not through raw counts.
- **Family identity and timing through time: yes.**
  - Each candidate has its family, time and P, so families can be followed 2003 → 2026 together with the forward
    catalogue.
  - `tpl_epoch` and `net_set` say which network saw it.
- **Firing rates and stacks for dv/v: yes, for ETS repeats.** P ≥ 0.75 frozen and P ≥ 0.90 era give 452 (2007) to
  10,891 (2004) counted-once repeats per year, with 0–20 scrambled-control candidates per year at the same cuts.
- **Between-episode rates: only as a bracket.**
  - Use `excess` (data − control).
  - The "rate" reading is not earned (G6). The control is not fair on active quiet days (c). The network-size factors
    differ between quiet days and ETS days.
  - The 2006 inter-ETS bursts (burst2006 note) are in the catalogue and flagged by `burst_day`. They are not dropped.
- **Not for:**
  - long-term trends of between-episode rate across 2005-07, 2005-11, 2006-07 or 2009-11 (the step checks cannot rule
    out network steps);
  - era-detector single repeats at P 0.5–0.9 (control-heavy);
  - CN4/B1-frozen counts as completeness estimates (1 % of Bostock's events).

## Brackets that remain

| quantity | bracket | why |
|---|---|---|
| CN4/B1 T_core (headline counts only) | [8.79, 9.82] / [8.69, 9.56] | P4-2b TRIM rule failed (D-P4b-2 kept as bracket) |
| K (to the S15 level) CN4 | 1.77 (95 % 0.86–3.65), and quiet-day 2014–17 value 0.56 | ladder vs invariance disagree ×3 |
| K B1 | 2.01 (1.16–3.48), quiet-day 0.67 | same |
| K B2 | 0.68 (0.41–1.12) | ladder and invariance agree within ×2 |
| K era | [0.66, 1.59] (via CN4 / B2 / B1 overlaps) | the three routes disagree; ETS-day k is 3–10× the quiet-day k |
| era vs frozen on ETS days | ×8–10 (overlap k), ×15–17 (recall ratio) | step check (d) |
| between-episode rates, every era | `excess` (control-subtracted) up to `n_core_data − n_core_c1` | G6 = control_subtracted_bracket |
| B2/B3 between rates | [all rows, n_z2 ≥ 2 rows] (1–4 % apart) | one-site-signal control (P4-2b item 4) |
| LRIV effect (2009-11-20) | not measured by the ladder | LRIV had no data in the ladder window |

## Open decision points

- **D-P4b-1, per-site cap:** needs a capped-weight matched-filter GPU run on 2005-11..2006-06. Not run.
- **Era P recalibration:** the era model could be refitted with quiet-day scrambled-control rows as negatives, or
  calibrated on non-ETS days. That would be a new pre-registered step; nothing was tuned here.
- **Step-check statistic:** a test that accounts for week-to-month activity swings (e.g. the placebo distribution as the
  reference instead of the day bootstrap) is needed before any 2003–13 vs 2014–26 between-episode claim.
