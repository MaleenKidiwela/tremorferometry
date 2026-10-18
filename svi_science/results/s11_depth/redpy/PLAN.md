# S11-REDPy pilot: sub-families, selective stacking and an all-pairs graph for LFE coda dv/v (rules fixed before any number)

Written Mon 2026-10-12, 2:55 pm PT, before any sub-family, node or graph value was computed. Background: `../MERLIN_IDEAS.md`
(Idea 2), `../REASONING.md`, the S10 pilot (`../../s10_dvv/pilot/`). User: "stop giving up"; stacking should rescue weak
single LFEs. ε = −dv/v (positive ε = slower); figures show dv/v = −ε. "Scrambled control" = the detector run with each
station site's data shifted to another time of the same day; at the PGC site its site-shifted windows are real LFEs (S10
Revision 1), so the scrambled-control windows used here are the S10 frame-time ones (`win/<NET.STA>__c1f.npz`), which hold
noise. Times in PT; dates are the catalogue's UTC dates.

## 0. Inputs (read only)
- Windows: `s10_dvv/pilot/win/<NET.STA>.npz` (2–8 Hz, 40 Hz, S−14 … S+10 s about the expected S, 3C, RMS-normalised).
- Rows: `s10_dvv/pilot/tables/rows.parquet`, data rows, family 10, **P ≥ 0.5**; epoch per row from `row_station.parquet`;
  picked S `tS_pick` per (family, station, epoch) from `series.csv`.
- Series (station-epoch): family 10 at CN.PGC E1, E2, E3; PB.B009, PB.B010, PB.B011 E1. Families 87 and 17 at the same
  stations if time allows (same rules, no change).
- Episode label per row: `s10_measure.clusters` with `r23/tables/episodes_2003_2026.csv` (the S10 labels).
- Realignment: `s10_measure.realign` (3C, ±0.5 s on S−1 … S+2 s) against the series' P-weighted grand stack, once
  (as `s10_fam_3c_26.py`). Raw stacks only (no SVD).
- Measurement: S-anchored stretch, 3C (each component demeaned and normalised, cc = mean of component cc's, missing
  components left out), lapse 2–5 s after picked S, ε grid ±1 % in 0.005 % steps + parabolic refinement. |ε| ≥ 0.99 % =
  no peak inside the search = non-measurement (S10 Revision 2), kept as a column.

## 1. Sub-clustering (REDPy core)
- Per series: zero-lag 3C cc between every pair of realigned single-event windows, S−1 … S+5 s (same per-component
  normalisation). Same cc between every scrambled-control window of that station-epoch (any of families 10/87/17,
  frame-time windows, realigned the same way against the same grand stack) and every data window.
- **Threshold (pre-registered):** thr = 99th percentile of the control–data cc of that series (if < 30 control windows in
  the series, pool the station's controls across epochs). A data pair above thr has ≤ 1 % chance of being noise-like.
- **Path A (headline, REDPy-like):** average-linkage hierarchical clustering on 1 − cc, cut at 1 − thr. A sub-family is a
  cluster with ≥ 20 members; other events are "unclustered" (kept as a pool, never dropped). Single and complete linkage
  at the same cut are reported as columns (counts/sizes only).
- **Path B (stack assignment, column; becomes the downstream sub-families only if Path A puts < 30 % of events into
  sub-families of ≥ 20):** k = 1 … 5 clusters, initialised from the average-linkage tree cut to k groups (groups < 20 merged
  into the nearest by stack cc), then 10 iterations of assigning each event to the leave-self-out cluster stack with the
  highest cc. k is chosen by 2-fold cross-validation (events split odd/even in time): stacks from one half, mean
  assigned-cc on the other half; k = the smallest k after which the CV gain is < 0.005. If k = 1, there are no sub-families.
- Reported per series: number of sub-families, sizes, fraction of events in them; per sub-family the years with ≥ 3
  members over years with any event (persistence); sub-family × episode contingency (episodes with ≥ 20 events):
  Cramér's V and permutation p (2,000 shuffles of labels across events of the series) = **episode-locked membership =
  composition effect**. Station consistency: adjusted Rand index between the PGC-E1 and B011/B009/B010 clusterings on
  the events they share, against 500 label shuffles (real source sub-families must agree across co-located instruments).
- Figures: cc matrix ordered by the dendrogram + dendrogram, per series.

## 2. Selective stacking (weights)
- w_P = P (baseline); w_cc = max(cc_i, 0) with the event's leave-self-out sub-family stack (the series stack for unclustered
  events / no sub-families), on S−1 … S+5 s; w_cc2 = w_cc².
- Compared by (i) half-split precision: 50 random disjoint half splits of the busiest episode into two stacks of
  n = 20/50/100, both stretched against the leave-episode-out series stack; σ = SD(ε_A − ε_B)/√2; (ii) the episode scatter
  of §5. Weight choice for the headline nodes is fixed now: **w_P** (the S10 baseline); w_cc and w_cc2 run as columns of
  the same nodes and graph. (Weights are a column, not a cut.)

## 3. Nodes (per series and per sub-family; "all" = every event of the series, ignoring sub-families)
- (a) day stacks: one per date with ≥ 3 events;
- (b) chunks of N = 20 and N = 50 events in time order **within an episode label**; a remainder < N/2 joins the previous
  chunk; an episode with fewer than N events gives one node if it has ≥ N/2 events (≥ 10 / 25), else none (its events
  appear in the day nodes).
- (c) episode stacks (all events of one episode label, ≥ 10 events) — used for the like-for-like comparison with S10.

## 4. All-pairs graph
- Every pair (i < j) of nodes of one node set (series × sub-family × node type × weight): ε_ij = stretch of node j
  against node i as reference, with cc_ij. Equation v_j − v_i = ε_ij. Railed pairs excluded (count kept).
- Weight 1/σ_ij², σ_ij = √(1 − cc_ij²)/cc_ij (relative weights; the scale comes from the misclosure). Solved with the
  S10 `graph_solve` (4 Huber IRLS passes, zero weighted mean per node set). **Misclosure error** σ_v = formal σ ×
  √max(χ²_red, 1). Triplet misclosure ε_ij + ε_jk − ε_ik reported (RMS over up to 20,000 random triplets).
- **Reference-stack method (comparison):** each node stretched against the leave-node-out P-weighted stack of the same
  node set; σ = 100-draw event bootstrap.
- Bootstrap: 50 draws resampling events with replacement inside every node, re-measuring all pairs and re-solving →
  σ_boot per node (node types (b) N = 50 and (c) only, and N = 20 / day for PGC E1 and B011 if time allows).

## 5. Tests (pass rules fixed now)
- **T1 injected step** (waveform level, every event dated on or after the boundary stretched by a = +0.1 % and +0.3 %
  about the picked S before node building; boundary = first date of the episode that holds the series' median event
  date). Recovery = Δ[wmean(post) − wmean(pre)]/a; leakage = Δ[mean(last pre-boundary episode) − mean(earlier pre
  episodes)]/a (S10 definitions; Δ = injected − natural). **Pass: median recovery 0.9–1.1 and |leakage| < 0.2 for both
  amplitudes**, per estimator (graph, reference) and node type.
- **T2 misclosure vs bootstrap:** pass if the median of σ_v(misclosure)/σ_boot over nodes is 0.8–1.25.
- **T3 scrambled-control nodes:** control windows of the station-epoch made into nodes the same way (day stacks ≥ 3 and
  chunks N = 20 in time order) and measured against the data series reference (leave-none-out stack of the node set's
  events). Pass if ≥ 90 % of non-railed control values satisfy |ε|/σ_boot < 2, AND the data-pair cc median exceeds the
  control-pair cc median by ≥ 0.2. The railed fraction is reported (it is the expected outcome for noise).
- **T3b time-shuffled data (not a scrambled control; a check of the error model):** events randomly permuted in time,
  then chunked (N = 20/50) and solved. Pass if episode means have χ²_red vs zero within 0.5–1.5 (the graph σ describes a
  series with no time structure).
- **T4 scatter vs the existing family-10 PGC numbers** (baselines registered now, computed from existing tables):
  - raw Z 2–4 s episode stacks (pilot `episodes.parquet`, 25 episodes): episode SD **0.20 %**;
  - raw 3C 2–5 s day stacks (`fam10_PGC_3c25_daily.csv`): day CC > 0.8 → 20 days, day SD 0.198 %, within-episode SD
    0.088 %, CC-weighted episode-mean SD 0.199 % (16 episodes); CC > 0.5 → 99 days, day SD 0.292 %, within-episode SD
    0.281 %, episode-mean SD 0.122 % (25 episodes); (SVD 3C 2–5 s, for context only: episode-mean SD 0.075 %).
  - Like-for-like here: reference-stack method on the same nodes.
  - Metric: episode-to-episode SD of the episode means (inverse-variance mean of the nodes in an episode; each series its
    own zero) and the within-episode day-to-day SD; plus the **excess scatter** √max(SD² − mean σ², 0).
  - "Better" = lower episode SD than the reference method on the same nodes AND T1 passed.
- **T5 sub-family-resolved dv/v:** sub-families solved separately (own zero each), combined per episode by inverse
  variance. Composition was part of the noise if its episode excess scatter is < 0.8 × the "all events" graph's excess
  scatter on the same episodes, and T3b passes for the sub-family series.
- **T6 co-located check (diagnostic):** PGC-E1 vs B011 node values for identical event sets (chunk nodes built on the
  shared events, 2005-11 → 2015-01): Pearson ρ, permutation p (2,000). ρ > 0 with p < 0.01 = a signal common to both
  instruments (medium change or source composition; not noise). Repeated after sub-family resolution.

## 6. Plots (user style)
dv/v points coloured by CC, no error bars, line graphs; cc matrix + dendrogram per series; sub-family time series.

## 7. Compute
CPU only, ≤ 8 workers, RAM < 40 GB; nothing written outside `redpy/` except NOTES lines; the WA-NW queue and the
ballistic test files are not touched.
