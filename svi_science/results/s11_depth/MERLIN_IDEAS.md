# MERLIN_IDEAS — frank assessment of four ideas (2026-10-12)
Model: Claude Fable 5.1, read-only (file saved by the main session). "Scrambled control" = station data shifted to another time of the same day.

## One fact that reorders the ideas
In `s10_dvv/pilot/tables/fam10_PS_timing_by_episode.csv`, PGC and B011 (same site, independent instruments) are each referenced to their own epoch grand stack. Their S-lag difference over the 8 shared episodes 2006–2014 is −10.3, −10.5, −8.4, −8.5, −12.0, −10.3, −9.5, −9.4 ms, a scatter of ≈ 1.2 ms.
- (a) cc timing on episode stacks reaches about 1–2 ms (a twentieth of a sample).
- (b) The 7 ms episode-to-episode S scatter at PGC is common to both instruments, so it is real (detector-frame drift or physical), not pick noise.
- (c) The 13–16 ms S−P scatter is limited by P: P cc 0.66–0.86, and B011 P cycle-skips in 5 of 19 episodes with the ±500 ms search. 0.2 % is a P-pick floor, not the method's floor.

## Idea 3 — ballistic P/S times. Plausible; highest upside; rank 1
- **Physics:** δ(tS−tP)/(tS−tP) = (−1.75 δβ/β + δα/α)/0.75 at Vp/Vs 1.75, so S−P amplifies a Vs change ×2.3.
- **Precision:**
  - 13 ms on 6.5 s already means path-averaged δβ/β ≈ 0.09 % per family–station–episode, and ≈ 0.25 % for the 25–45 km leg alone. That equals the coda cell, but along a straight ray.
  - It also sees the mid-crust (≈ 0.5 of the S time), which the coda cannot (a_m 0.05).
  - With S at 1–2 ms the deep leg reaches ~0.05 %, if P is brought to the same quality (full P wavetrain, 3C/polarised P, ±100 ms search, n ≥ 50).
- **Obstacles:**
  1. There is no LFE origin time, so only S−P and inter-station S differences are usable.
  2. Source-locus wander: 1 km shifts S−P by ~115 ms, so a 3-parameter locus nuisance per family–episode is needed, separated by its moveout pattern. That needs ≥ 4 well-spread stations, which SVI has only marginally.
  3. The same per-series offset problem: only episode contrasts are possible.

## Idea 2 — REDPy-style repeater dv/v. Adopt one piece; rank 2
- **What they do:** Hotovec-Ellis et al. measure coda delays between individual event pairs, then solve all pairs by least squares (Brenguier-type all-pairs graph), with closure residuals as the error.
- **What doesn't transfer:** single-event LFE coda is below noise, so individual repeats and cluster-centroid picking don't carry over.
- **What transfers:** the all-pairs graph, with episode stacks as nodes. S10 found it passes injection better (39/42 vs 28/42). The gain is calibration and an internal error; no new precision.

## Idea 1 — per-SSE patch maps. Not realistic at patch scale; rank 3
- **Why:** only 14 % of the deep coefficient sits within 10 km of the source; the rest is smeared along the path.
- **Precision:**
  - per-patch deep σ ≈ 0.7 % per episode, ≈ 0.2 % pooled over 15 majors;
  - coarse along-strike zones (3) may reach 0.2–0.3 %.
- **Confound:** the front. Each family's sample is a 2–5-day window at front arrival.

## Idea 4 — full-waveform/adjoint. No; the middle ground yes; rank 4
- **Why not:** 2–8 Hz over 50–80 km needs ~50 m grids (10⁹–10¹⁰ cells) per pair, and the coda scattering isn't in any model.
- **Middle ground:** linear joint inversion of ray-traced S and S−P rows plus the existing residence rows, with the same unknowns and the same Σ = 1 test. That is how idea 3 should be built.

## Ranked cheapest tests
1. **Ballistic:** extend s10_ps_timing.py to all core stations and families, with a ±100 ms search; form inter-station S differences and S−P per episode.
   - Common across families → path change.
   - Family-specific → locus wander.
   - Flat at ~1–2 ms → nothing there.
2. **Graph estimator:** Q3 with episode stacks as nodes. Pass = deep-pulse recovery 0.9–1.1 and misclosure within 20 % of bootstrap σ.
3. **Front confound:** superposed-epoch coda ε vs days-since-front, with the scrambled control beside it (~0.03 % per lag-day).
4. **Zone resolution:** design-only σ for 3 along-strike deep zones; go only if σ ≤ 0.2 %.

What would change my mind: S differences < 2 ms scatter with no family-common signal (ballistic has nothing to measure), or per-patch design σ < 0.3 % (idea 1 live).
