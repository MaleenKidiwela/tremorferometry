# S11-REDPy pilot — results in plain words

Written Mon 2026-10-12, ~3:40 pm PT. Rules: `PLAN.md` (written 2:55 pm PT, before any number); deviations and post-hoc
additions are listed at the end. Scripts `svi_science/scripts/s11_redpy_*.py`; tables `tables/`; figures `figures/`.
Families 10, 87, 17 at CN.PGC (E1/E2/E3) and PB.B009/B010/B011, P ≥ 0.5 repeats, raw stacks, 3C S-anchored stretch
2–5 s, 2–8 Hz. dv/v = −ε. "Scrambled control" = detector run with station data shifted to another time of the same day
(frame-time windows, i.e. noise). Everything is a pilot on 3 families and 4 co-located instruments; treat it as provisional.

## What the data show (short)
1. **Stacking does rescue the weak single LFEs.**
   - One LFE window correlates with another at a median cc of only 0.08–0.23.
   - A 50-event stack measures dv/v to ≈ 0.2 % per node.
   - Episode means at PGC for family 10 scatter by **0.08–0.10 %**. The old raw Z 2–4 s episode stacks gave 0.20 %, and raw 3C daily stacks gave 0.12 %.
   - The gain comes from 3 components + the 2–5 s window + per-event realignment. The reference-stack method on the same nodes gets the same 0.09–0.10 %.
2. **There are no distinct REDPy-type sub-families.**
   - The cc matrix is a smooth continuum of event strength (`figures/cluster_fam*_*.png`), not blocks.
   - The biggest cluster is simply the strongest events. Its stack has the same waveform as the rest (cc 0.86–0.99).
   - Its share changes from episode to episode, so the event mix (by strength) is episode-locked.
   - **Resolving dv/v by sub-family does not reduce the scatter.** Median SD ratio is 0.92–1.13 with the pre-registered labels and 0.90–0.99 with a finer waveform split (path B). Family 17 is mixed (0.5–1.7).
   - So composition, as clustering sees it, is not the main noise.
3. **One real thing clustering found: B009 in 2007 and 2009.**
   - Almost every 2007 and 2009 event at B009 falls into its own waveform group: family 10 61/65 and 47/49 events, family 87 36/42 and 31/32.
   - The same events at PGC, B010 and B011 do not.
   - That is a B009 station-era effect, not a source change. B009's 2007/2009 dv/v points (the −0.25 to −0.47 % values in `figures/fam10_graph_vs_reference.png`) are suspect.
   - Confirming test: B009 response/orientation metadata for 2007–2009.
4. **Selective stacking helps modestly.**
   - Weighting by cc with the sub-family stack lowers the half-split σ by a median 13 % (cc/P 0.87; cc² 0.89; better in 68 % of 38 cases).
   - It lowers the episode scatter by 0–20 %. Family 10 PGC, N = 50 nodes: P 0.083 %, cc 0.077 %, cc² 0.070 %.
   - Step injection still passes with cc weights.
5. **The all-pairs graph is calibrated only for nodes of ≥ 50 events. Its misclosure error is not honest.**
   - Recovery at N ≥ 50 and on episode nodes is 0.92–1.00. It shrinks changes on day and N = 20 nodes (0.64–0.91).
   - The misclosure-based σ is 2.5–4 × smaller than the event-bootstrap σ in every series.
   - With misclosure σ, events shuffled in time give χ²_red 10–25 (should be ~1).
   - Calibrated for its shrinkage, the graph is about as precise as the reference-stack method. This matches Merlin's forecast: calibration and consistency, no new precision.
6. **Nothing time-structured stands above the noise at the ~0.1 % (PGC) / ~0.15–0.2 % (borehole) episode level.**
   - Family 10 PGC episode means scatter as the bootstrap predicts (χ²_red 0.66–0.68).
   - Shuffling the events randomly in time gives the same episode scatter as the real order.
   - The boreholes show a small excess (χ² 1–3, ≈ 0.1 %). That is consistent with the known ~0.07–0.08 % annual cycle, but not tested here.
7. **Co-located instruments agree chunk by chunk on identical events, also within an episode.**
   - ρ = 0.28–0.43 for PGC vs boreholes in families 10/87, weaker in family 17.
   - Within one episode (hours apart) the chunk-to-chunk residuals still correlate (family 10, N = 20: ρ 0.24–0.61).
   - A velocity change is not expected to flip back and forth within hours. This common part is most likely source-side: which LFEs and how much overlapping tremor go into a chunk.
   - It is barely explained by S-arrival timing (site-mean realign lag), event strength or P (|ρ| ≤ 0.34 at N = 20, signs not consistent across families).
   - More events per stack averages it down. Averaging co-located instruments does not.
   - Confirming test below.

## 1. Sub-clustering (PLAN §1)
`tables/subfam_summary_fam{F}.csv`, `subfam_members_*`, `subfam_station_ari_*`, `subfam_stackcc_*`, `composition_sub1_share_*`,
`subfam_pathB_cv_*`; figures `figures/cluster_fam{F}_{NET.STA}_{E}.png` (cc matrix in dendrogram order, dendrogram,
data vs scrambled-control cc histograms, membership in time).

- **Single-event cc:**
  - data–data median: family 10 0.13–0.17, family 87 0.16–0.23, family 17 0.08–0.11;
  - scrambled control–data 99th percentile: 0.19–0.25 (the pre-registered threshold);
  - data pairs above it: 8–43 %.
- **Path A (average linkage at the threshold):**
  - family 10: 29–49 % of events in groups of ≥ 20; one large group (339–770 events) plus 0–7 groups of 20–49;
  - family 87: 44–64 %;
  - family 17: 0–19 %, so path B was used (as pre-registered).
  - Single linkage chains everything into one cluster (100 %); complete linkage leaves 0–15 % in groups.
- **What the large group is:**
  - median single-event cc with the grand stack 0.45–0.56 vs 0.25–0.40 for the rest;
  - its stack vs the rest's stack: cc 0.97–0.99 (families 10/87) and 0.86–0.93 (family 17, smaller groups).
  - So it is the strong events, not a different waveform. The small groups' stacks reach cc 0.72–0.82 with others, as expected for 20–50-event stacks.
- **Persistence:** the large group appears every year (family 10 PGC-E1 12/12 years; boreholes 16–19 of 17–20 years). Small groups appear in 1–9 years.
- **Station consistency (same events, co-located instruments):**
  - large-group membership agrees with ARI median 0.21 (family 10), 0.14 (87), 0.38 (17); 95 % of shuffles < 0.015; p = 0.002. It is set by event size, so it is shared.
  - Path B (3–5 balanced groups, CV gain ≥ 0.005): ARI 0.03–0.07. Mostly station-specific, only a trace shared.
- **Episode-locked membership:**
  - among clustered events, path A groups vs episodes: Cramér's V 0.13–0.18, p 0.13–0.70 at PGC, B010 and B011 for family 10. Not locked.
  - Exception B009: V 0.46 (family 10) and 0.72 (family 87, path B), p 0.0005, from the 2007/2009 group (above).
  - The strong-event share per episode is locked: V 0.11–0.39, p ≤ 0.007 in 11 of 18 series and p < 0.05 in 14. The share ranges e.g. 13–47 % across PGC-E1 episodes for family 10.
  - So the mix of strong vs weak events changes between episodes; the waveform type does not.

## 2. Selective stacking (PLAN §2)
`tables/halfsplit_fam{F}.csv`, `halfsplit_pivot_*`.
- **Half-split σ at n = 20/50/100 in the busiest episode:**
  - cc weights / P: median 0.87 (PGC 0.87, B009 0.95, B010 0.87, B011 0.84); cc²/P 0.89; cc better in 68 % of 38 cases.
  - Family 10 PGC-E1: P 0.34/0.19/0.16 % vs cc 0.25/0.14/0.13 %.
- **Episode SD, graph (N = 50 nodes), P → cc → cc²:**
  - family 10: PGC 0.083/0.077/0.070, B010 0.145/0.117/0.110, B011 flat;
  - family 87: −4 to −11 %;
  - family 17: mixed.
- Injection with cc/cc² weights: recovery 0.98–1.01, leakage ≤ 0.01 (graph, N50/episode). Calibrated.

## 3–4. Nodes and all-pairs graph vs reference stack (PLAN §3–4)
`tables/nodes_fam{F}.parquet` (every node: v, misclosure σ, bootstrap σ, reference ε, cc), `episode_means_fam{F}.parquet`;
figures `fam{F}_graph_vs_reference.png`, `fam{F}_subfamily_series.png`, `overview_episode_scatter.png`.
- **Pair cc** (median, family 10): N = 20 nodes 0.38–0.49; scrambled-control nodes 0.07–0.09.
- **Railed pairs** (no peak in ±1 %): few at N ≥ 20.
- **Triplet misclosure RMS** (family 10): 0.6–0.85 % at N = 20 and day nodes, 0.3–0.5 % at N = 50, 0.15–0.4 % for episodes.
- **Graph and reference series track each other closely** at every station (fam10 figure). The graph's episode means are a little smoother because of the shrinkage at small N.

## 5. Tests (pass rules from PLAN §5)
| Test | Result | Verdict |
|---|---|---|
| T1 injected step, graph | N50: recovery 0.955–0.996, episode: 0.916–1.001, leakage ≤ 0.04 → **pass** (all 3 families). N20: 0.71–0.91 (pass only family 87), day: 0.79–0.90 → **fail** (leakage day 0.05–0.10). Path-B sub-family day nodes 0.64–0.73. | graph usable at N ≥ 50 |
| T1 injected step, reference | recovery 1.00–1.08, leakage ≤ 0.06 → **pass** everywhere except family 17 sub-family N50 (1.35) | pass |
| T2 misclosure σ vs bootstrap σ | ratio 0.25–0.39 (median per family and node type), all series | **fail** |
| T3 scrambled-control nodes | 72–75 % railed (no peak = noise, as expected); non-railed within 2σ: 96 / 96 / 90 % (families 10/87/17); pair cc gap data − control 0.26–0.37; series passing both clauses 8/12, 11/12, 6/11 | consistent with zero |
| T3b time-shuffled data | χ²_red of episode means vs zero: misclosure σ 10–25 → **fail**; bootstrap σ: N50 1.4–1.5 (borderline pass), N20 1.7–2.8 (fail) | error model: use bootstrap σ; still ~1.2–1.6 × low at N20 |
| T4 scatter, family 10 PGC | episode SD: graph N50 0.083 %, episode nodes 0.094 %; reference 0.090 / 0.101 %; vs baselines 0.20 % (raw Z 2–4 s episode stacks) and 0.122 % (raw 3C 2–5 s day→episode, CC > 0.5). Within-episode day SD: graph 0.21 %, reference 0.34 % (baseline 0.28 % at CC > 0.5) | ×2 better than raw Z; graph ≈ reference once calibrated |
| T4 other stations/families | episode SD 0.14–0.25 % at the boreholes; PGC 0.15–0.17 % (family 87), 0.12–0.25 % (family 17) | — |
| T5 sub-family resolved | episode SD ratio sub/all median 0.92 (10), 0.92–1.13 (87), 0.80 (17); path B column 0.90–0.99; excess-scatter ratio < 0.8 in 16/32, 7/35, 19/31 comparisons, no consistent station or node type | **no** consistent reduction |
| T6 co-located | N = 20 chunks on identical events: ρ PGC–borehole 0.28–0.37 (family 10), 0.34–0.43 (87), −0.02–0.37 (17); borehole–borehole up to 0.60; within-episode residuals still correlated (family 10: 0.24–0.61) | common, source-side component |

- Calibrated precision (σ / recovery), family 10, node level:
  - N = 20: graph 0.227/0.78 = 0.29 % vs reference 0.31 %;
  - N = 50: 0.20 % vs 0.20 %.
  - Equal.
- Real vs time-shuffled episode SD, family 10 (per series, N20/N50):
  - PGC-E1 real 0.06/0.07 % vs shuffled 0.09/0.13 %;
  - B009 0.12/0.17 vs 0.12/0.13; B010 0.10/0.15 vs 0.09/0.11; B011 0.14/0.15 vs 0.12/0.13.
  - The real time order adds at most ~0.1 % (boreholes, N50) over random event order.

## Reading
- Stacking works: weak single LFEs (cc ≈ 0.1–0.2) become 0.2 % measurements at 50 events and ≈ 0.1 % episode values at PGC.
- The REDPy idea of separating sub-families does not add precision here, because LFE family 10/87/17 waveforms form one population that differs mainly in strength.
- Weighting by cc (selective stacking) is the useful part of it: ~10–15 %.
- The all-pairs graph is a good consistency check at N ≥ 50 with bootstrap errors. Its closure errors under-state the noise ×3, so they must not be used as the error bar.
- The scatter that remains is mostly noise plus a component shared by co-located instruments for the same events, which no clustering removed.
- What would change this reading: a sub-family split that is shared across co-located instruments (ARI ≫ 0.1) AND lowers the episode SD by > 20 % at all four instruments.

## Confirming tests (not run)
1. **Source-side floor:** draw 200 random 20/50-event subsets from one big episode (2004-07, 2005-09), measure each at PGC and B011.
   - The part of the subset-to-subset spread shared by the two instruments is the composition / overlapping-tremor floor.
   - If it equals the within-episode common component found here, only more events per stack lowers it.
2. **B009 2007/2009:** metadata (response, orientation, gain) for B009 2007–2009; rerun B009 without that group (column).
3. **Borehole excess (~0.1 %):** fit the annual cycle to the graph episode means with bootstrap σ; the excess should disappear if it is the known seasonal term.

## Deviations from PLAN and post-hoc additions (logged)
- The graph's a-priori pair σ uses the Weaver et al. (2011) constant times √(1−cc²)/(2cc). The scale is set by the plain-LS misclosure before the Huber passes. Only the relative weights matter, as planned; the constant makes the Huber cut meaningful.
- Graph zero: plain mean over connected nodes (S10 `graph_solve` convention), not weighted.
- Bootstrap σ: all four node types for all events (P), N50/episode for sub-family sets; none for injected or shuffled runs (they use the natural run's σ per node).
- Injection boundary: if the median episode is the series' first, the next episode start is used (fix made before the family 87/17 runs; family 10 PGC-E2 has only 2 episodes and no injection value).
- T3b was run on all-event series only (not per sub-family).
- Post hoc (after first numbers):
  - episode means weighted by bootstrap σ ("graph_boot");
  - T3b with bootstrap σ;
  - calibrated SD;
  - strong-event share per episode;
  - within-episode residual correlation and the composition-driver diagnostic (`s11_redpy_composition.py`, `tables/composition_drivers_*`);
  - path-B sub-families as a column for families 10 and 87 (`*_labB.*`).
- The first family-10 launch over-subscribed BLAS threads (load ~45); it was killed by PID after 1.5 min and rerun with 1 thread per worker. No numbers came from it.
