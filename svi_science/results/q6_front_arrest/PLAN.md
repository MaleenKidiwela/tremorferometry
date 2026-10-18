# Q6 — What stops an ETS front? Front arrest vs time since last slip, and the Strait of Juan de Fuca as a boundary

Pre-registered 2026-10-07 4:55 pm PT, before any Q6 number (episode, start, stop, regression or control) was computed.
Only the PNSN catalogue's latitude/longitude distribution was looked at (to place the margin axis). Changes after this
point go in §9 with a PT time and reason. Scripts `svi_science/scripts/q6_*.py`; outputs here; figures
`svi_science/figures/q6_*.png`. CPU only, < 30 GB RAM. Merlin big picture item 3 / prescription step 8.

## 1. Questions and hypotheses

- **Q6.1 front arrest.** H-refractory: a front fails to enter, or stops at, a stretch of fault that slipped recently
  (coefficient of log10 days-since-last-active on "front continues / bin active" > 0). H-reverse: fronts stop where the
  fault has been quiet longest (coefficient < 0). Both directions are tested two-sided; the sign says which.
- **Q6.2 Strait boundary.** (a) ETS starts cluster at the Strait of Juan de Fuca more than a uniform-along-margin control
  and a control matched to where tremor (or LFE families) exists. (b) Fronts stop at the Strait more than those controls;
  fronts that reach the Strait cross it less often than they cross other stretches of the same width.

## 2. Data

- **Primary (detection-independent of the LFE work): PNSN tremor** `catalogs/pnsn_tremor_cascadia_full.csv`
  (2010-01-06 → 2026-05-27 in the file; 748,002 epicentres; 39.5–50.9 N). All tremor, not only big episodes.
- **Margin axis.** Polyline through the tremor band: nodes = (median lat, median lon) of epicentres in 0.2° latitude bins
  with ≥ 2,000 epicentres, longitude smoothed by a 3-node running mean. Each epicentre → s = arc length of its projection
  (km), d = perpendicular distance; |d| > 60 km dropped. s = 0 at the projection of 48.30 N 123.40 W (mid-Strait between
  Port Angeles and Victoria, the S3 cross-region origin); s > 0 toward NW (Vancouver Island), s < 0 toward SE.
  **Strait zone** = |s| ≤ 15 km (columns 10 and 25 km). Domain = first to last node; **end zones** = 25 km at each end
  (starts/stops there are "at a domain end" and excluded from boundary tests; stops there are censored, not arrests).
  **Northern domain** (SVI + Washington, the S3 cross-region axis) = s ∈ [−200, +200] km.
- **Bins** 5 km in s; days in PT.
- **Episodes (components).** A bin-day is active if it has ≥ 3 epicentres. Active bin-days are linked if |Δbin| ≤ 3
  (15 km) and |Δday| ≤ 2; connected components (union-find) = episodes. Every component counts as an episode of the
  history ("any episode incl. minor"); components with < 10 epicentres are kept in the history only.
  **Big** = upper class of the optimal two-class split of log10(epicentres) over components with ≥ 30 epicentres
  (exhaustive, minimum within-class sum of squares; as B). Columns: big = n ≥ 300; n ≥ 1000.
- A bin is **active in an episode** if ≥ 3 of the episode's epicentres fall in it (column ≥ 1, ≥ 10).
  Median time of those epicentres = the bin's activity time in that episode.
- **Start** of an episode = median s of its first 10 % of epicentres (column: first 2 days). **Front run** = maximal
  set of active bins reached from the start bin walking NW and SE, allowing single inactive bins (a gap of one bin is
  crossed, two consecutive inactive bins end the run). **Stops** = last active bin of the run on each side (+ bin half
  width). Episodes whose start bin is inactive use the nearest active bin.
- **LFE family-resolved checks.**
  - SVI ours 2014+: `s0/analysis_table.parquet`, series data, run prod, P_same ≥ 0.75 & keep_xfam8 (counted once);
    5-km bins of family `along_strike_km` (S0 SVI frame); episodes = `episodes_svi_final.csv` (79, MAJOR + minor),
    padded as S3; big = S3 big (10). Bin active = ≥ 3 events. Start = S3 `lfe_start10`.
  - Bostock 2003–2013: `p3_bostock/bostock_events_mapped.parquet` (his catalogue covers only his 10 ETS windows, so
    history = those 10 only); 5-km bins of `along_km`; start = median of first 10 %.
  - Joint SVI 2003–2026: Bostock 10 + our big 10, history from big episodes only (like-for-like).
  - NW E3 (column only; Merlin pitfall: ~29 % of events from one source, lumpy positions): `s3/events_nw.parquet`
    headline (selected_v2 & counted_once), 5-km bins of the A/B source positions (NW frame); episodes
    `episodes_nw_final.csv` E3; big = S3 big (14).
  - For Q6.2 the LFE check puts SVI and NW events on the PNSN margin axis; overlapping big SVI/NW episodes
    (S3 `cross_region_overlaps` rule: padded windows overlap) are merged into one joint episode before the start is taken,
    so a front crossing the Strait is not split into two "starts" at the two box edges.

## 3. Q6.1 tests (front arrest)

Predictors for (bin b, big episode e): **T_any** = days from e's first epicentre back to b's activity time in the most
recent earlier episode (any size) where b was active; **T_big** = same, earlier big episodes only. Predictor =
z-scored log10 T. Rows with no earlier activity since 2010-01-06 are censored → dropped (column: kept at the lower
bound). Column: T_any7 = only earlier episodes that ended ≥ 7 days before e started (removes e's own precursors).
Covariates: dist = |s_b − s_start| (km); log10 size of e (epicentres); prop = log10 of b's all-time epicentre count
excluding e (detectability / usual activity).

- **Primary statistic (walk / stop hazard).** For each big episode and each direction, trials = every bin of the front
  run beyond the start bin (y = 1, front continued into it) plus the first bin after the run end (y = 0, front stopped),
  unless that bin is in an end zone or outside the domain (censored: no stop trial). Logistic regression
  y ~ z(log10 T) + dist + log10 size + prop; statistic β_T.
- **Secondary (enter / fail to enter).** All bins within ±150 km of the start inside the domain (excluding end zones):
  y = bin active in e; same model; statistic β_T.
- **Referee = episode-order shuffle**: component footprints (their bins and activity times, shifted with the slot) are
  permuted among the episodes' start dates — big among big slots, all others among other slots — and T is recomputed
  for every row; outcomes, distances, sizes and prop stay with the episode. 2,000 draws, seed 0. This keeps each bin's
  usual rate (so a "quiet bins are poorly detected" artefact appears in the shuffle too) and destroys memory.
- **Pass rules.** H-refractory SUPPORTED if β_T(T_any) > 0 with one-sided order-shuffle p < 0.01 in the PNSN primary
  (full margin, walk) AND β_T > 0 in the northern-domain subset AND β_T > 0 in the SVI LFE check (ours 2014+).
  H-reverse SUPPORTED with the same rule for β_T < 0. PARTIAL if only the PNSN primary passes. Otherwise NOT SUPPORTED
  (neither direction). The same rule is reported (not as a separate pass) for T_big, the secondary model, and columns.
  Bostock and joint SVI: reported with order-shuffle p (p < 0.05 noted).
- **Figure statistic**: stop fraction (y = 0 share among walk trials) in quintiles of T_any, real vs the 5–95 % band of
  the shuffle.

## 4. Q6.2 tests (Strait boundary)

Big PNSN episodes (primary); all episodes ≥ 30 epicentres as a column. Starts and stops in end zones excluded.

- **(a) Starts.** Statistic = number of starts in the Strait zone, for (i) the northern domain (primary: starts in
  s ∈ [−200, 200]) and (ii) the full margin. Controls (10,000 draws each): **uniform** = each start uniform over the
  domain (minus end zones); **matched** = each start drawn from the all-time PNSN epicentre density along s (same
  domain); **own-footprint** = each start replaced by the median of a random 10 % of the episode's own epicentres
  (S3 control). Starts cluster at the Strait if p < 0.01 against uniform AND matched in the northern domain;
  own-footprint reported (p < 0.05 noted).
- **(b) Stops.** Statistic = number of front stops (both ends) in the Strait zone. Controls: uniform, matched, and
  **rigid shift** = each episode's run moved by a random offset keeping its length inside the domain. Fronts stop at
  the Strait if p < 0.01 against matched AND rigid shift in the northern domain.
- **(c) Crossing.** For every 30-km zone of the northern domain (step 5 km), fronts that ARRIVED (run reaches the zone's
  near edge while moving toward it from outside) and CROSSED (run passes the far edge). Strait crossing rate vs the
  distribution over other zones (rank) and a Fisher test vs all other zones pooled; descriptive.
- **Detectability check (diagnostic).** All-time epicentre density in the Strait zone vs neighbouring 30-km zones; if
  the Strait is a density minimum, stops/starts there are flagged as possibly detection-driven.
- **LFE check.** Joint SVI+NW LFE big episodes (2015+), start on the margin axis; controls uniform over the
  family-covered span and matched to the family positions (SVI families + NW A/B sources weighted equally per family /
  source), own-events random 10 %. Reported with p; small n, no pass rule. S3's "5 of 24 starts near the Strait"
  restated on the margin axis with and without merging.

## 5. Pitfalls handled

NW LFE lumpiness → PNSN primary, NW LFE column only. Detectability at region ends → end zones excluded / censored;
prop covariate + order shuffle carry detection-rate effects. LFE box edges at the Strait (SVI SE end, NW N end) →
joint episodes. PNSN network change over time → order shuffle keeps each episode's size; reported per era
(2010–2014, 2015–2026) as a column.

## 6. Outputs

`q6_build.py` → `axis_nodes.csv`, `episodes_pnsn.csv` (all components ≥ 10 epicentres: start, run ends, size, big),
`bin_episode.parquet`; `q6_tests.py` → `walk_rows_*.csv`, `tests.json`, `tests_table.csv`; `q6_lfe.py` → LFE checks;
`q6_figures.py` → `figures/q6_fronts_margin.png` (front extents per big episode over time, start/stop marks),
`figures/q6_arrest_vs_tsince.png`, `figures/q6_strait.png`. RESULTS.md: plain answer first.

## 9. Changes after pre-registration
- 2026-10-07 ~4:58 pm PT (code, before the full run; no rule changed): the first margin-axis build used a flat x = lon·cos(lat)
  projection that stretched the axis (~1,900 km); replaced by an azimuthal-equidistant projection centred on the Strait
  reference point (axis now −961 … +317 km; Victoria +9 km, Port Angeles −17 km, Sequim −31 km, Olympia −147 km).
- ~5:00 pm PT (implementation details the plan left open): own-footprint start control computed as per-episode
  probabilities (300 random-10 % subsamples each) combined exactly (Poisson-binomial) instead of 10,000 joint draws;
  rigid-shift stop control uses offsets U(−100, +100) km per episode (an offset that would push the run outside the domain
  is set to 0); in the LFE walk, a stop trial is only made where the next bin holds families (else censored at the box
  edge). LFE Q6.2 p-values are binomial for uniform / matched.
- Facts found when running (no rule changed): SVI ours has events in 52 of the 79 listed episodes; Bostock's 10 ETS
  cover every 5-km bin of his families (no stop and no "not entered" bin) so the Bostock and the Bostock half of the joint
  test carry no arrest information; the NW events file covers major windows only, so the NW LFE history = majors.
- 2026-10-07 ~5:00 pm PT POST-HOC, after the main results: `q6_posthoc.py` re-ran the walk model on rows with
  T_any ≥ 14 d and ≥ 60 d (same shuffle, 1,000 draws), because the linking rule (2 d / 15 km) mechanically prevents a
  front entering a bin another component used within ~2 days.
