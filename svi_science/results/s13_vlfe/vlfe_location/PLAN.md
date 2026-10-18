# vlfe_location — PLAN (pre-registered 2026-10-16 ~10:02 am PT, before any offset number)

Trigger: user 2026-10-16 "these [VLFEs] are not located right". Question: how far are the catalogue VLFE locations from
where the concurrent LFEs are, is the difference systematic, how sharp is the location actually, and does a better
Green's function help?

Inputs (read only): `../catalogue/tables/` final event set = `prod<Y>_events_v2.csv` for 2010–2015 and
`prodr<Y>_events_v2.csv` for 2016–2018 (the set behind `prod_*_final.csv`, 242 tier events). Location there =
`lat_best/lon_best` = scan node of the highest stored scan z within ±30 s of the matched-filter time (`cat_e2.cheap`),
on a 10-km grid whose depth is fixed to Slab2 (so depth is not free). Grid templates = whole-space thrust (315/12/90),
flat ×2 free surface (`scripts/s13_green.py`, `scan_common.template`), per-station-day gains `g` from `daycfg`.
Nothing in `../catalogue/` is written or modified; every output goes to this folder.

## Test 1 — independent location check against concurrent LFEs
- Events: every tier event (ETS and between) with ≥ 1 located LFE within ±10 min.
- LFE set (primary, the same set the catalogue uses for "with LFEs"): SVI P ≥ 0.75 (`s13_common.svi_rows`, data series,
  counted once) + WA-NW `catalogue_v2.parquet` `selected_v2`, each LFE placed at its family location
  (`data/catalogue/families_v1_3.csv`; SVI by family_id in region S_Vancouver_Is, WA by family_uid). Families with no
  location are dropped (count kept).
- Centroid = mean of family lat/lon over the LFEs in the window (each LFE weight 1, so busy families weigh more).
  Columns: ±10 min (primary), ±2 min, SVI-only ±10 min, median-family location, n LFEs, n families.
- Offset = VLFE − LFE centroid: east km, north km, distance km. Reported per year and overall: mean vector, median
  distance, spread (SD east, SD north, distance IQR).
- Scrambled control: each VLFE location is paired with the LFE centroid of a different event of the same episode
  (1000 random within-episode re-pairings; between-ETS events are re-paired within the same calendar year).
  - Tracking: real median distance vs the scrambled distribution (one-sided p = fraction of re-pairings with median ≤ real).
  - Bias: mean offset vector, 95 % interval from a bootstrap over episodes (resampling whole episodes). A bias is called
    systematic if this interval excludes 0 in either component. The scrambled re-pairing shows what mean offset arises
    when VLFE and LFE locations are unrelated within an episode.
- Raw figure first: map of VLFE locations with lines to their concurrent LFE centroids, coloured by year (grid nodes,
  stations, family cloud shown); offset histograms (east, north, distance) real vs scrambled. No error bars.

## Test 2 — achieved location resolution
- For each tier event: recompute the production scan z over all grid nodes at the event's best time (`cat_inj.local_z`,
  production run `real_v2`, same grid, gains, MAD normalisation). Check: z at the best node equals the stored `z_best`.
- Width = the region with z ≥ z_best − 1 (Δz ≤ 1): number of nodes, area (nodes × 100 km²), equivalent radius
  √(area/π), and the largest distance from the best node; also Δz ≤ 2 as a column. Also: whether the LFE centroid's
  nearest node lies inside the Δz ≤ 1 region.
- Compare the median width with the LFE-offset spread from Test 1.

## Test 3 — template dependence
- Template location: H&G T1–T5 from `tables/hg2018_table1.csv`; C0–C5 from `../catalogue/tables/e3_bank.csv`
  (member mean); relay RC0–RC7 from the 2018 relay bank member mean. Mapped to the nearest grid node.
- Fraction of events whose best node is within one grid step (≤ 10 km, edge neighbour; ≤ 14.2 km incl. diagonal as a
  column) of the detecting template's node (`template` column) — and of either agreeing template (`agree`).
- Scrambled control: same fraction when each event is assigned a random template of the same bank instead.
- Also: is the best node nearer the template node than the LFE centroid (fraction)?

## Test 4 — quick improved-location test (2011–2013 tier events with LFE centroids; the 14 vetted/near-miss events
are a named subset)
- Same data windows (production band-passed displacement, cleaned windows, station weights 1/σ², on/off flags), same
  node grid, same Slab2 depths; statistic for all variants = weighted variance reduction of a thrust template with free
  M0 ≥ 0 and one common time shift (±10 s) — so the variants differ only in the Green's function.
  - (0) baseline: whole-space GF (production templates, production gains g).
  - (a) 1-D layered model (USGS CVM v1.7 / Crosson 1976): only if a layered GF code is available offline (pyrocko
    qseis/qssp backend binaries or stores, CPS). If none: not run, reported as such (no downloads).
  - (b) half-space Okada static × band-passed step (pyrocko.modelling.okada; point-like 1×1 km thrust patch at the
    node depth; Poisson solid; per-station delay r/β, β = 3.95 km/s), production gains g.
  - (c) (b) with empirical station-component amplitude terms calibrated on the S13 pooled ETS stacks: per
    station-component factor = observed stack amplitude / Okada prediction at the stack's LFE count-weighted centroid,
    pooled over the episodes whose stacks carry that station; replaces g. Also (c0): the same stack terms with the
    whole-space GF.
- Outcome: median and mean distance to the LFE centroid per variant, and the mean offset vector; a variant "reduces
  offsets" if the median distance drops by ≥ 3 km (≈ 1/3 grid step) AND the paired per-event change is negative in
  ≥ 60 % of events. The baseline (0) must reproduce the production best node for most events (report the fraction);
  if it does not, that is reported first.

## Rules
- ≤ 2 CPU workers, < 30 GB memory, no downloads; PT times; "scrambled control", never "null"; no error bars in figures.
- Every number goes to RESULTS.md with its table; failures and changes to this plan go to NOTES.md as dated amendments.
