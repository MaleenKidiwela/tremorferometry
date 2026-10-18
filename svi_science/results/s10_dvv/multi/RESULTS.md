# S10 multi-family dv/v: do LFE families agree on which slow-slip episodes are fast or slow?

Written Sun 2026-10-11, ~10:45 am PT. Rules: `PLAN.md` (written first, ~10:10 am PT). Revision 1 was added after the
first numbers were read and is labelled post hoc wherever it is used. Scripts: `svi_science/scripts/s10_multi_*.py`
(queues `s10_multi_queue*.sh`). Tables are in `tables/`, figures in `figures/`. dv/v = −ε, in %. "Scrambled control"
means noise windows: the detector's c1 rows read at their frame time, passed through the same path. All numbers are
provisional.

## What the data show (plain words)
- **Yes, the families agree, and what they agree on is the time of year.** Across 76 families and 45 major episodes
  (2003–2026), dv/v is higher (faster) in episodes that happen in late summer and lower (slower) in episodes in late
  winter. Pooled over all families, this annual cycle has an amplitude of **0.07–0.08 %** (peak to peak ≈ 0.15 %).
  The peak is around 21 Aug – 3 Sep and the low around late Feb – early Mar
  (`figures/fig_seasonal.png`, `figures/fig_pooled_series_es.png`).
- It is seen at PGC (surface broadband) and at the boreholes. It is seen in 2003–2014 and again in 2015–2026 with the
  same size and timing. It is seen in both coda windows, smaller at 3–5 s (≈ 0.6 ×).
- **Once the annual cycle is removed, the families no longer agree on anything.** The cross-family correlation drops
  to the level of the shuffled episodes. The remaining episode-to-episode scatter is what the errors predict
  (χ²_red 1.10, p = 0.31). So we see no episode-specific dv/v, for example a big ETS that is slower than the season
  would predict, larger than ~0.05–0.08 % per episode.
- **Within an episode, nothing changes as the front passes.** Pooled over 75 families and 5,689 station-days, the
  slope is −0.0002 ± 0.0004 %/day (0–20 d after the front), and days 0–3 minus later days is +0.003 ± 0.008 %. That
  bounds any within-episode change to ≲ 0.025 % (3σ).
- **No long-term trend:** −0.0004 ± 0.002 %/yr over 2003–2026 (pooled majors, episode stacks).

## Do families agree on episode-to-episode changes, and how big?
- Each family's episode value is noisy (median σ 0.22–0.25 %), and the common signal is ~0.06–0.08 %. So a single pair of
  families barely agrees: the median pairwise r is 0.03–0.07, and 54–58 % of pairs are positive. Each family against the
  mean of all the others gives a median r of 0.16–0.19, positive for 64–72 % of families. The pooled series is where
  the common signal stands out.
- **Pre-registered headline (CC²-weighted means of daily values, core 5 stations, 2–4 s): narrowly NOT passed.**
  M = 0.034 against a 99th shuffle percentile of 0.036 (p = 0.014, 10,000 draws). The PGC-only subset gives M = −0.014
  (only 29 eligible pairs, because a PGC-only family-episode value needs ≥ 3 measurable days at one station). The
  borehole-only subset gives M = +0.009 (p = 0.31).
- **Per-episode stacks (secondary unit in the plan; full battery added post hoc): passes the same rule.** One stack per
  family, station and episode (median 21 repeats), against the leave-episode-out reference:
  - core 5: M = 0.073, above the 99th percentile 0.039, p = 0.0002;
  - PGC only: M = 0.044, p = 0.002;
  - boreholes only: M = 0.058, p = 0.001;
  - 3–5 s: M = 0.029, p = 0.023;
  - majors + bursts: M = 0.073, p = 0.0002.
  The two units agree on the same family-episodes (r = 0.82, slope 1.0). The stack is the less noisy unit: fewer days
  are dropped as "no peak", and one stretch fit is made on all repeats.
- The pooled PGC-only and borehole-only series correlate across major episodes (r = 0.34–0.38, p = 0.02–0.03, 40
  episodes), so the agreement does not live at one instrument type.
- Size: the pooled series' SD over major episodes is 0.10 % with errors of 0.08–0.09 %, so the common signal is
  ≈ 0.04–0.06 % RMS. Almost all of it is the annual cycle. The share of family-episode variance explained by the
  episode (ICC) is 0.03 (ANOVA F = 1.7, p = 0.003–0.006).
- Leave-one-family-out: the pooled series changes by a median 0.02–0.03 % per episode (max 0.21 %, in an episode with
  few families). The minimum r between a leave-one-out series and the full series is 0.94–0.96, so no single family
  drives it.

## What is new
- An **annual dv/v cycle in the coda of deep LFEs**, recovered from 23 years of SVI slow-slip episodes sampled at
  ~14-month intervals that drift through the calendar. It is coherent across 76 families, PGC and four boreholes, and
  both halves of the record:

  | unit, stations | amplitude (%) | peak (day of year) | permutation p |
  |---|---|---|---|
  | episode stacks, core 5 | 0.081 ± 0.011 | 233 (~21 Aug) | 0.0001 |
  | day means, core 5 | 0.066 ± 0.010 | 246 (~3 Sep) | 0.0002 |
  | 2003–2014 (stacks) | 0.083 ± 0.018 | 225 | — |
  | 2015–2026 (stacks) | 0.087 ± 0.018 | 236 | — |
  | PGC only (stacks) | 0.085 ± 0.020 | 217 | — |
  | B011 (stacks) | 0.096 ± 0.033 | 210 | — |
  | B009 (stacks) | 0.126 ± 0.034 | 183 | — |
  | B010 (stacks) | 0.162 ± 0.016 | 279 | — |
  | B926 (stacks) | 0.077 ± 0.030 | 240 | — |

  Per-station and split-half fits: `tables/posthoc_seasonal.csv`. Permutation = pooled values shuffled among
  episodes, 10,000 draws.
- Fast in late summer and slow in late winter is the pattern expected from groundwater / pore pressure in the shallow
  crust (wet winters). This would make the 2–4 s LFE coda a near-receiver monitor of the same shallow seasonal cycle
  that ambient-noise studies see. It is a sensitivity check of the method. It is not a deep-fault signal. This
  interpretation is not tested here.
- The deep-fault answer: **no episode-specific or within-episode dv/v** above ~0.05–0.08 % per episode or ~0.025 %
  within episodes. This is consistent with the earlier ETS bound (notes/2026-09-04_SSE_vs_dvv_RESULTS.md), now from
  P-clean LFE repeats on single response epochs.
- VGZ and SNB (added stations) do not show the cycle (amplitudes 0.04 ± 0.05 / 0.04 ± 0.06 %). They are noisier, with
  only 37–59 % of series passing causality. Adding them dilutes the agreement (core 5 + VGZ/SNB: M = 0.013 day means,
  0.030 stacks).

## Confirming tests and outcomes
| test | rule | outcome |
|---|---|---|
| cross-family agreement, headline (day means) | M > 99th pct shuffle AND M > 0 in PGC-only AND borehole-only | **not passed** (0.034 vs 0.036; PGC −0.014; boreholes +0.009) |
| same rule, episode stacks (post hoc) | same | **passed** (0.073 vs 0.039; PGC +0.044; boreholes +0.058) |
| scrambled control through the same path | control must not pass | **not measurable**: only 7–14 % of control day values and 49 % of control episode stacks have a cc peak inside ±1 %. 0 family pairs share ≥ 8 control episodes, so no pair test. Control pooled series: no annual cycle (p 0.35–0.90), but it is too noisy (σ ≈ 0.19 %) to exclude a 0.08 % cycle. |
| agreement after removing the annual cycle (post hoc) | — | gone: M = 0.008 (p 0.29) stacks, 0.013 (p 0.17) day means; boreholes −0.005 |
| lapse, pooled 3–5 s vs 2–4 s | slope ≥ 0.6 and sign ≥ 70 % | stacks **pass** (0.64, 86 %, r 0.78); day means **fail** (0.49, 74 %, r 0.71). Noise baseline from the shared 3–4 s: unit-level slope 0.49–0.52, sign 67–68 %. |
| annual cycle in subsets | same sign/phase | PGC, B011, B009, B010, B926 and both halves: all positive in summer, peaks Jul–early Oct |
| leave-one-family-out | stable | min r 0.94–0.96 |
| within-episode change | > 3σ, control within 2σ, same sign PGC/boreholes | **no change** (slope −0.0002 ± 0.0004 %/d; contrast +0.003 ± 0.008 %) |
| causality (coda/pre-P > 1.5) | per series | 716 of 1,077 series; boreholes 62–69 of 71–76, PGC 146/185, VGZ 118/273, SNB 190/320 |
| reuse check | pilot ε identical | family 10/87 at B011 reproduced the pilot ε24 exactly (max diff 0) |

## Caveats
- B009, B010 and B011 sit at the PGC site, so "PGC vs boreholes" separates instruments, not sites. B926 is the only
  second site, and it shows the cycle too (0.077 ± 0.030 %, p = 0.05).
- The scrambled control cannot test for a path artefact of this size (see table). The argument against an artefact is
  that the cycle appears on two instrument types, in two independent halves of the record, at both lapse windows,
  with the expected hydrological phase, and that noise windows show no cycle.
- Each (station, epoch) series is zeroed on its own epoch reference. PGC E2 (2015-01 → 2017-08) therefore has a short
  zero, which can offset its values. The borehole subset is free of this.
- The episode-stack verdict, the deseasoning and the annual-cycle tests were added after the headline result was
  read (PLAN.md Revision 1).

## Files
- Tables:
  - `families.csv`: 76 families;
  - `daily.parquet`: 20,165 station-day values;
  - `episode_stacks.parquet`;
  - `fe_*.csv`: family-episode values, σ in the tables;
  - `agreement.csv`: every set;
  - `pairs_*`, `fam_vs_rest_*`, `perm_*.npy`;
  - `pooled_series_{dm,es}.csv`, with family-bootstrap σ and leave-one-out;
  - `within_episode*.csv`, `posthoc_seasonal.csv`, `result.json`.
- Figures:
  - `fig_episode_panels_{a,b}.png`;
  - `fig_heatmap_{dm,es}.png`;
  - `fig_pooled_series_{dm,es}.png`;
  - `fig_corr_matrix_{dm,es}.png`;
  - `fig_permutation.png`;
  - `fig_seasonal.png`;
  - `fig_within_episode.png`;
  - `fig_lapse_{dm,es}.png`.
- Windows: `win/` (no new cut stage; nothing deleted).
