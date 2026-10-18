# S3b — do the deepest SVI families fire ahead of the slip front? (pre-registered 2026-10-07 ~12:40 pm PT, before computing)
Front arrival is taken from PNSN tremor (independent of our LFE detection): for each big SVI episode with a migrating front
(S3 'migrates' = True), each 5-km along-strike bin's arrival = time of the bin's 10th-percentile tremor epicentre (bins with
≥ 20 epicentres). For every headline LFE (data, P ≥ 0.75, counted once) in that episode: lag = time − arrival(bin of its
along-strike position). "Early" = lag < −1 day. Deep = families in the deepest third by down-dip coordinate (dd_km),
shallow = shallowest third (thirds over the 115 families).
Statistic: Δ = early fraction(deep) − early fraction(shallow), pooled over episodes (event-weighted) and per episode.
Control: permute the deep/middle/shallow label among families within the same 10-km along-strike bin (5,000 draws) → p.
Second control: same analysis with each episode's tremor arrival map taken from a different episode (mismatched fronts).
Claim "deep families fire ahead of the front" only if pooled Δ > 0 with permutation p < 0.01, the sign holds in ≥ 7 of the
migrating episodes, it holds with flagged episodes (2017-03, 2025-10) removed, and the mismatched-front control shows no such Δ.

## Result (2026-10-07 ~12:45 pm PT): NOT supported by the pre-registered rule
Pooled Δ = +0.017 (deep early 2.4 % vs shallow 0.7 %), permutation p = 0.004, but positive in only 5 of 10 migrating
episodes (rule ≥ 7) and driven by 2015-12 (deep 13.6 % early); the mismatched-front control gives a larger Δ (0.18), so
the measure is not specific to the real front. Median lag deep 2.3 d vs shallow 2.9 d (not a pre-registered test).
Verdict: the "deep families fire ahead of the front" impression from the S3 figures is one episode, not a general pattern.
Files: result.json, per_episode.csv; script scripts/s3b_deep_early.py.
