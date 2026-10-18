# SVI backward extension — first look (sanity check, not science)

Written Sat 2026-10-10 (PT). **PROVISIONAL**: the confidence model P is not yet applied, the P4 decision and the
step checks are not done, and the era-detector threshold is a placeholder. Nothing here is scored against Bostock's
event catalogue; only his ETS window *dates* (`backward_scope/truth_episodes.csv`) are shaded as background.

## In plain words

- **Do the ETS stand out?** Yes, every Bostock window with output stands out clearly in both detectors. The frozen detector
  shows a median of 17–165 extra candidates per day inside the windows, against 0–4 between them. The era detector
  shows 683–996 extra per day in 2003–2005, about 23–166 times its scrambled control.
- **Excess between episodes?** Yes, but small. On days more than 30 d from any window, data beats the scrambled control
  every year. Over each year the frozen detector has 1.5–4.1 times as many data candidates as control candidates, and
  the era detector 5–9 times. The median daily excess is only 0–4. Before P is applied we cannot tell how much of this is
  minor tremor and how much is detector behaviour that the scrambled control does not reproduce.
- **Steps at station-set changes?** Data minus control shows no obvious jump. The *levels* do change, though. The frozen
  detector's control rate per observable family-day is about 0.01 in 2003–2005 (CN4/B1, T 10.75/10.5). It is about
  0.04 after B2 (2005-11-02, T 9.0) and 0.10–0.22 in 2010–2013 (S15, T 8.75). So the per-set thresholds do not yet
  give a constant control rate. P4-2 matched them to about 0.09 using station subsets of 2010–2013 data. On the real
  2003–2005 data the CN4/B1 thresholds come out about 10 times stricter than that target. The era detector's control
  rate falls from 0.054 to 0.021 when the POLARIS/CN array comes on (Aug 2003).

## What was counted

- Runs: frozen `p5f_2003`…`p5f_2013` (2007–2009 = ETS windows only: 2007-01-01→04-03, 2008-03-10→08-02,
  2009-03-10→08-07; 2010–2013 are the forward-detector years) and era `p5e_2003_kelb`, `p5e_2004`, `p5e_2005`,
  `p5e_2006` (2006 still running, to 2006-06-30). There is no era output after mid-2006. In total 4,538 run-days;
  a few days with unreadable or empty files were skipped.
- Per UTC day: candidates with `net_mad ≥ T`, series `data` and `c1` (scrambled control: each station site shifted to
  another time of the same day), and observable family-days per series (from `_fd.parquet`). There is no HF veto
  (same as the P4 control-rate rule) and no cross-family de-duplication, so one event seen by several families
  counts several times. The PT day starts 7–8 h after the UTC day.
- Frozen T_core per station set (`p4/tables/p4_2_decision.json`): CN4 10.75 (to 2005-06-30), B1 10.5 (from 2005-07-01),
  B2 9.0 (from 2005-11-02), B3 9.0 (from 2007-11-13), S15 8.75 (from 2009-08-27). Boundaries are the station on-dates in
  `bw_gate_common.py`.
- Era detector: `p4_3_decision.json` has no value (flag validation failed; `p4_3_era.csv` is empty), so the counts are
  shown at **8.75 and 9.0, both provisional**.

## Table: median daily excess, ETS-window days vs days > 30 d from any window

Excess = data − scrambled control per day. Ratio = sum data / sum control over those days. Full table with rates per
observable family-day: `table_ets_vs_between.csv`.

| detector | year | ETS days: n, median excess, ratio | between (>30 d): n, median excess, ratio |
|---|---|---|---|
| frozen | 2003 | 13, 32, 7.9 | 290, 0, 2.2 |
| frozen | 2004 | 26, 17, 20.4 | 279, 0, 1.8 |
| frozen | 2005 | 23, 25, 11.4 | 281, 0, 2.2 |
| frozen | 2006 | — (no window) | 357, 2, 1.9 |
| frozen | 2007 | 11, 53, 2.4 | 29, 0, 1.5 |
| frozen | 2008 | 24, 150, 3.3 | 62, 1, 1.8 |
| frozen | 2009 | 29, 42, 6.9 | 62, 1, 4.1 |
| frozen | 2010 | 24, 95, 3.2 | 281, 2, 1.7 |
| frozen | 2011 | 22, 165, 2.7 | 283, 2, 1.9 |
| frozen | 2012 | 27, 154, 2.9 | 278, 4, 1.7 |
| frozen | 2013 | 31, 81, 2.8 | 274, 3.5, 1.8 |
| era @ 8.75 (prov.) | 2003 | 13, 683, 22.7 | 264, 1, 7.0 |
| era @ 8.75 (prov.) | 2004 | 26, 996, 166 | 275, 1, 6.0 |
| era @ 8.75 (prov.) | 2005 | 23, 891, 93 | 279, 1, 6.4 |
| era @ 8.75 (prov.) | 2006 (to Jun) | — | 180, 2, 5.4 |
| era @ 9.0 (prov.) | 2003 | 13, 607, 32.7 | 264, 1, 8.9 |
| era @ 9.0 (prov.) | 2004 | 26, 918, 226 | 275, 1, 8.3 |
| era @ 9.0 (prov.) | 2005 | 23, 797, 132 | 279, 0, 8.0 |
| era @ 9.0 (prov.) | 2006 (to Jun) | — | 180, 1, 6.7 |

The 2007 frozen ETS days are only 11 because the frozen 2007 run ends 2007-04-03. The window itself is 2007-01-24→02-03.

## Things to look at (not claims)

- **2006 noise months (frozen, B2 set, T 9.0).** January, March and May 2006 have control rates of 0.19–0.27 per
  family-day, against about 0.04 in the neighbouring months. These are burst days on which the control tracks the data
  (e.g. 2006-01-15: data 454, control 378), so the scrambled control is doing its job. One exception: 2006-03-28/29
  (data 731/732, control 243/278) has a large excess and no Bostock window. Worth an eye review later.
- **Era detector, late Sep–Oct 2005.** Several days have control ≥ data (up to −33 per day), which looks like noise bursts.
- The frozen control rate in 2003–2005 is about 10 times below the P4-2 target, so the CN4/B1 thresholds look
  conservative on real data. This is an input to the P4 decision, not a decision.
- The era detector's between-episode data/control ratio (5–9) is much higher than the frozen detector's (1.5–4).
  In panel (c) the era data comes in blocks of weeks. That is either minor inter-ETS tremor that the denser POLARIS
  array sees, or something site-specific. P and the step checks should settle which.

## Files

- `fig_daily_excess_2003_2013.png`: (a) frozen daily data − control, (b) era daily data − control at 8.75 and 9.0,
  (c) 30-d mean rates per observable family-day for data and control, between-episode days only.
- `daily_counts.csv` (per run-day counts at every threshold used), `table_ets_vs_between.csv`, `table_steps.csv`
  (60 d before/after each station-set change, between-episode days).
- Scripts: `svi_science/scripts/bw_firstlook_counts.py` (4 CPU workers, read-only on `mf/feat`) and
  `svi_science/scripts/bw_firstlook_plot.py`.
