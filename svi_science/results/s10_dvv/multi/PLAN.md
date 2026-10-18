# S10 multi-family dv/v — do LFE families agree on which slow-slip episodes are fast or slow? (rules fixed first)

Written Sun 2026-10-11, ~10:10 am PT, before any window was read or any number computed. User-approved scale-up of the
pilot (`../pilot/RESULTS.md`, `../pilot/PLAN.md`). ε = −dv/v (positive ε = slower); every plot shows dv/v = −ε.
"Scrambled control" = the detector's c1 rows (each station site's data shifted to another time of the same day), read
here at their detector frame time t + offset (pilot Revision 1, `c2`): windows at same-day times unrelated to what the
station contributed, i.e. noise through the same path.

## 0. What is reused unchanged from the pilot
Rows (forward `s0/analysis_table.parquet` prod + keep_xfam8, P = P_same; backward `svi_catalogue_v1_4_backward`
primary + keep_xfam8, P = P), active day = Σ P of counted-once data rows with P ≥ 0.5 ≥ 3 (UTC date), P ≥ 0.5 rows
P-weighted; window reader (cut store where the row is in it, raw day files otherwise), 2–8 Hz 4-pole zero-phase, 40 Hz,
despike 8 × MAD, per-window RMS normalisation; grand stack per (family, station, epoch), S re-pick, cc-realign ±0.5 s,
S re-pick; S-anchored stretch ±1 % (0.01 % grid, parabolic), windows τ = 2–4 s (headline) and 3–5 s after S, Z
component; reference = leave-day-out epoch grand stack; 200-draw bootstrap σ; Revision 2 rail rule (|ε| ≥ 0.99 % = no
peak = missing). **Raw day stacks, no SVD** (the SVD step shrank steps to 0.44). Response epochs never crossed (each
(family, station, epoch) series is its own zero). Code: `s10_windows.task`, `s10_measure.{Ref, realign, pick_S, wmean,
boot_stacks}` imported, not copied. Not run here (characterised in the pilot): step injection, pairwise graph, precision
vs n, horizontals.

## 1. Families and stations
- Families: every SVI family with ≥ 10 active days spread over ≥ 10 distinct calendar years (2003–2026). Counted
  before this plan from the catalogue alone: **76 families, 3,834 family-days** (the user expected ~30–50; the rule as
  given yields 76 and all 76 are used).
- Stations (by response epoch, pilot table): PB.B011, PB.B009, PB.B010, PB.B926, CN.PGC (E1/E2/E3) = **core 5**; CN.VGZ
  (E1–E4) and CN.SNB (E1–E5) are added (same reader, cheap) and reported as a with/without check, not in the headline.
- Windows are read only for the 73 new families; the pilot windows of families 10, 87, 17 are reused by row key.

## 2. Per family
- Daily series: one value per (family, station, epoch, active day) as in the pilot (`tables/daily.parquet`).
- **Episode units** (shared by all families):
  - major episodes: 2003–2009 from `p3_bostock/episodes_bostock_final.csv` (6, all major); 2010+ from
    `b_ets_record/episodes_b.csv`, region SVI, major == True (46). A day belongs to an episode if it lies in
    [start − 2 d, end + 2 d] (UTC dates). Class `major`.
  - between-episode bursts: every active day of any selected family outside all major windows; their union is split
    into bursts wherever two consecutive dates are > 7 d apart (the pilot's cluster rule). Label `B<first date>`. Class
    `burst`.
- **Family-episode value** (headline): dv/v_fe = CC²-weighted mean of −ε24 over all non-railed (station, day) values
  of family f in episode e at the core 5 stations; weight = cc24² (cc of the day stack with its reference). Kept if
  ≥ 3 values. Error (tables/text only): σ_fe = √(Σ w²σ²)/Σ w with the bootstrap σ, and the weighted scatter / √n_eff
  as a second column; the larger is used. Same for 3–5 s (weight cc35²), for the PGC-only and borehole-only
  (B009/B010/B011/B926) subsets, for core 5 + VGZ/SNB, and for the scrambled control (−c2_eps24, weight c2_cc24²).
- Per-episode stacks (all P ≥ 0.5 windows of the family in the episode, leave-episode-out reference, pilot code path)
  are computed as a secondary column only.

## 3. Cross-family agreement (the key test)
- Headline set: major episodes, core 5, 2–4 s.
- Pairwise: for each family pair sharing ≥ 8 episodes (both values present), Pearson r of their episode values
  (Spearman as a column). Statistic **M = median r over eligible pairs**.
- Family vs rest: for each family, Pearson r of its episode values with the CC²-weighted mean of all *other*
  families' values in the same episodes (≥ 8 shared). Median reported.
- Permutation: 10,000 draws; in each draw every family's values are permuted among that family's own episodes
  (seed 20261011); M recomputed with the same eligible pairs. p = (1 + #{M* ≥ M}) / 10,001.
- Subsets: M recomputed on PGC-only and borehole-only family-episode values (each with its own 10,000-draw
  permutation, reported).
- **PASS** = M (core 5) > 99th percentile of its permutation distribution **AND** M > 0 in the PGC-only subset **AND**
  M > 0 in the borehole-only subset (the station-specific systematic trap of notes/2026-09-04_SSE_vs_dvv_RESULTS.md:
  the agreement must not live at one instrument type only). Reported alongside: r between the PGC-only and the
  borehole-only pooled episode series.
- **Scrambled control** through the same path: M_ctrl and its permutation percentile. If the control also passes, the
  agreement is called a path artefact, whatever the data show.
- **Lapse check**: (i) M in 3–5 s; (ii) pooled episode series 3–5 s against 2–4 s: Theil–Sen slope and sign agreement;
  pass = slope ≥ 0.6 and sign agreement ≥ 70 % (pilot rule d). Because 2–4 s and 3–5 s share 3–4 s, noise alone
  makes them correlate; the same slope and sign agreement for the daily values (the pilot's noise level) are reported
  as the baseline.
- Secondary sets (reported, not the verdict): majors + bursts; bursts only; core 5 + VGZ/SNB; P-weight unchanged
  but cc-unweighted mean (sensitivity to the weighting).
- Size: amplitude of the common episode signal = √max(0, var(pooled series) − mean σ²_pooled), and the fraction of
  family-episode variance explained by the episode (ICC from a one-way random-effects ANOVA on the family-episode
  values).

## 4. Pooled episode series
- Per episode: CC²-weighted mean of −ε24 over all families and core stations (all non-railed day values). Error
  (text/tables only) = SD over 1,000 bootstrap resamples of families. Leave-one-family-out: max |Δ| per episode and
  the minimum Pearson r between a LOFO series and the full series. Also for PGC-only and borehole-only.
- Long-term: weighted LS slope of the pooled major series vs time (per epoch-free borehole subset as the clean case).

## 5. Within episodes
- x = days since the front reached the family: S5/S3 front arrivals (`s5/step3_front_arrivals.csv`, 10 big episodes
  2015–2024) where available; else days since the episode start. Days are UTC dates; x = date − front date.
- y = day value minus the CC²-weighted mean of that (family, station, epoch, episode) (needs ≥ 2 days). Pooled
  CC²-weighted means in 1-day bins from −2 to 30 d; weighted LS slope over 0–20 d and the contrast (days 0–3) −
  (days ≥ 4); errors by family bootstrap (1,000). Same for the scrambled control and for the PGC-only / borehole-only
  subsets. "Consistent change" = |slope or contrast| > 3σ AND control within 2σ AND same sign in both subsets.

## 6. Figures (user plot rules)
dv/v plots: points coloured by CC with the reference, marker = station, **no error bars**, averages CC²-weighted;
uncertainties in text and tables only.
- `fig_episode_panels_*.png`: one subplot per major episode, all families and core stations, x = days since episode
  start, colour = CC, marker = station, red line = CC²-weighted mean.
- `fig_heatmap.png`: family × episode values (majors, then bursts).
- `fig_pooled_series.png`: pooled 2003–2026 episode series, colour = mean CC, marker size = number of families; second
  panel PGC-only vs borehole-only (marker = station group).
- `fig_corr_matrix.png` + `fig_permutation.png`: pairwise r matrix; permutation histograms (core 5, PGC-only,
  borehole-only, 3–5 s, scrambled control) with the observed M.
- `fig_within_episode.png`, `fig_lapse.png`.

## 7. Resources and conduct
GPU not needed (CPU stretch); ≤ 12 workers; RAM < 100 GB; `scripts/mem_watchdog.sh` kept alive; kill by PID only;
nothing deleted; no new cut stage (windows read into `multi/win/`). Long steps run as a nohup queue
(`scripts/s10_multi_queue.sh`) with DONE markers.

## Revision 1 (Sun 2026-10-11, ~10:30 am PT — AFTER the headline agreement numbers were read; post-hoc, labelled so)
What the first run showed that the plan did not foresee:
1. **The scrambled control cannot be run through the pair test.** Only 7–14 % of control day values have a correlation
   peak inside ±1 % (noise windows do not correlate with the coda reference; that is what noise should do), so 0
   family pairs share ≥ 8 control episodes. The pre-registered control result is therefore "no measurable values",
   not "pass" or "fail". Added: scrambled-control **episode stacks** (all c2 windows of the family's episode days,
   equal weights, at most the episode's number of data windows, same leave-episode-out reference), which carry more
   windows per value; measure re-run for this (daily values unchanged, same seeds).
2. **The per-episode stack variant (secondary column) agrees more strongly than the headline day-mean variant.**
   Added, post hoc, the same battery for it: PGC-only, borehole-only, 3–5 s, scrambled control, majors + bursts,
   pooled series, lapse. A family-episode value needs ≥ 2 station values (≥ 1 for the one-station PGC-only set). The
   pre-registered verdict stays the day-mean one; the episode-stack battery is reported as a post-hoc finding.
3. Bug fix: the family-bootstrap σ of a pooled episode with < 3 families is 0 or meaningless, which made the trend
   fit singular. Amplitude, χ² and trend use only episodes with ≥ 3 families.
4. Added check (post hoc): annual (seasonal) sinusoid fit to the pooled episode series, since a shallow seasonal
   velocity change would be common to all families at the same station.
5. Added (post hoc, after a test run showed a strong annual cycle in the pooled series): annual-sinusoid amplitude,
   peak day of year, F-test and a 10,000-draw permutation p (pooled values shuffled among episodes); and the
   cross-family agreement re-run after subtracting the pooled annual cycle from every family-episode value
   (`*_deseason` rows in `tables/agreement.csv`), to ask whether families agree on anything beyond the season.
