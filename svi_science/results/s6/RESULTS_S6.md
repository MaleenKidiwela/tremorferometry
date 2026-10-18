# S6 — LFE size per detection, southern Vancouver Island (written 2026-10-07 ~1:30 pm PT by the main session from the agent's outputs)

## Plain answer
- **LFE sizes did not change over 2014–2026.** Per family, the trend is −0.008 log units per decade (95 % CI −0.038 to
  +0.023, p 0.70): sizes stayed within about ±8 % per decade. This is the clean result of S6, and it bears on the
  "does the deep fault change over a decade?" question: by LFE size, no detectable change.
- **The two "seen" differences are not about the LFEs.** LFEs in big episodes read 1.5× bigger (ΔM +0.18) and
  between-episode LFEs read smaller (−0.15) — but the noise reference (6,000 low-confidence detections) shows the same
  shifts (+0.18 and −0.20). The measured amplitude rides on the tremor background: busy episodes are noisier. Test 2
  passes the pre-registered rule only by 0.002 (noise effect 0.175 vs 0.177) — treated as not interpretable.
- **No size change behind the front** (+0.03, CI −0.01 to +0.07): LFEs are not bigger or smaller once the front has passed.
- **Depth: unresolved.** Bostock's family magnitudes fall with depth (Spearman −0.62); our distance-corrected sizes show
  no depth trend (+0.01) and do not correlate with his (−0.02). Either our distance correction or his magnitude scale
  differs; size-vs-depth claims are on hold.
- Typical LFE magnitude ~1.7 (5–95 %: 1.25–2.16, approximate, anchored on Bostock's family magnitudes); within a
  family sizes vary ~6× (ΔM −0.40 to +0.39). Only 34 % of events are above the completeness amplitude (32 % in big
  episodes, 46 % outside them).

## Tests (pass rule in PLAN_S6.md: permutation p < 0.01, CI excludes 0, sign stable, noise effect smaller than the ΔM effect)
| test | effect (log10 amplitude) | 95 % CI | p | noise reference | verdict |
|---|---|---|---|---|---|
| 1 behind − at front | +0.031 | −0.009 … +0.070 | 0.062 | +0.111 | not seen |
| 2 big − small episodes | +0.177 | +0.118 … +0.244 | 0.0005 | +0.175 | passes by 0.002; read as tremor-background effect |
| 3 between-episode − ETS | −0.147 | −0.269 … −0.016 | 0.0005 | −0.196 | not seen (noise effect larger) |
| 4 trend per decade | −0.008 | −0.038 … +0.023 | 0.70 | −0.021 | not seen → bound ±~8 %/decade |
| 5 depth (ours) | Spearman +0.009 | −0.24 … +0.26 | 0.97 | — | not seen; Bostock's mags −0.62 (conflict) |

## Measurement notes
- Deviation 1 (before any test): plain RMS amplitude is buried in tremor (median S/N 1.17); headline = matched-filter
  amplitude ratio to each family × station stack (median S/N 1.69; Spearman 0.80 with RMS).
- Station steps (full data, matched filter): PGC 2015-01-12 −0.25 log units (PGC 1.7× high before; listed gain change 1.4 %)
  — absorbed by the per-station-epoch correction; NLLB 0.006, VGZ 0.025, SNB 0.008, PGC 2017 0.004, LRIV −0.046.

Files: tests_table.csv, tests.json, sizes_annotated.parquet, station_steps_mf.csv; figures svi_science/figures/s6_*.png.

## Addendum (2026-10-07 ~1:45 pm PT): one shared completeness threshold — post-hoc, not pre-registered
Why: the pre-registered "complete" flag uses each family-month's own threshold, which rises with the tremor background.
A noisier group (big ETS, mid-episode) therefore keeps only its bigger events by construction. Here each family gets ONE
threshold shared by the two groups being compared (the higher of the two groups' median thresholds). Only events above
it count, and each group needs ≥ 3 events (script `svi_science/scripts/s6_common_threshold.py`, 300 shuffles, so the
smallest possible p is 0.003; output `results/s6/common_threshold_check.json`).
- **Big vs small ETS:** +0.044 log units (≈ 11 % larger amplitude), CI +0.022…+0.067, p 0.003, 12 families. About three
  quarters of the original +0.177 was the threshold / tremor-background effect. A small real increase may remain, but it
  rests on 12 families and is provisional until Washington-NW repeats it.
- **Behind vs at the front:** +0.038 (≈ 9 %), CI +0.008…+0.071, p 0.013 (47 family × episode cells). A hint that later
  LFEs are slightly bigger; it misses the p < 0.01 bar.
- **Outside episodes vs ETS:** −0.039, CI −0.080…+0.002, p 0.027, 12 families. Not seen. With the shared threshold,
  between-episode LFEs are about the same size as ETS LFEs. The 76 S7 events are too few for this check (0 families with
  ≥ 3).
- **Plain words:** once both groups are cut at the same size floor, LFE size barely changes. At most it changes ~10 %
  with episode size and position in the episode, compared with the ~6× spread within a family.
