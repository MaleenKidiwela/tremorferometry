# 2026-08-26..31 — User visual-QC pass: family/date vetoes, era ruling, tooling

## Standing veto mechanisms (NEW, both honored by assembler + map builder)
- `data/family_user_blocklist.csv` — per-family veto (stem, full fam id, map label, date, reason).
  Subtracted at the certification step in BOTH `assemble_res_catalog_daily.py::cert_fams` and
  `build_borehole_dvv_map.py::certfams`.
- `data/station_date_blocklist.csv` — per-station DATE-RANGE veto; rows dropped pre-gates in the
  assembler (logged as `date_veto` in drops) and pre-_blocks in the map builder.
- `scripts/family_from_map_index.py` — resolves the HTML map's "N/M" family label to the stored id,
  replicating the map EXACTLY (preferred CSV unif05>eps05>orig; cycle-skip prefilter on wide files;
  cc>0.6; certified minus blocklist; >=60 distinct days; sorted ids; **label is 1-BASED**:
  fam[fi-1]). The /M is a checksum — mismatch refuses to resolve.
  ⚠ HISTORY: an earlier 0-based reading mis-identified BUCK "55/87" as c96; correct is
  `44.500_-123.600__c23` (methods artifact Fig 1 corrected 2026-08-31).

## Vetoed (user visual QC of family dv/v series in the HTML map)
FAMILIES (229 total; see the CSV for full ids/labels):
- KBO 2 (c87, c286) | MOGU 4 (c135,c255,c32,c44) | SBAR 12 (c155 BOTH, c28 9/34 only, c129,c84,
  c109,c113,c69,c100,c22,c65,c7) | TRIN 4 (c67,c107,c86 BOTH) | KHBB 15 (13 + c28 BOTH) |
  JRO 7 | WLOO 5 (c6,c78,c38,c93,c66)
- WHOLE STATIONS: KMR (2 fams), SEAS2 (172 fams — big Oregon/SW-WA coverage change 45.1-46.9N),
  LEID (6 fams)
DATE WINDOWS (4): PETY 2022-08-30..2025-08-08 (middle instrument era + tail);
B031 after 2014-10-01 (the fleet's worst offset outlier — post-2014 era was driving |c|~1.5%);
SSO2 after 2025-12-07; BTB after 2018-01-01 (renders its 2020 era boundary moot).

## USER RULING — era splitting restricted (2026-08-31)
Only stations with a VISUALLY CONFIRMED cross-era dv/v change keep per-era offsets:
`data/era_split_whitelist.txt` = WISH, SCOT, KHBB, TRIN, LCCR, YACT.
The other 26 era_table stations are now SINGLE-era in the assembler (no straddle drops, one c per
family). Rationale: no per-era reference exists in the measurement; where no step is visible the
extra offsets only cost slow-signal capacity (c-channel).
RESOLVED (user, 2026-08-31, after one reversal): **PETY ADDED to the era whitelist** — its middle
era's data stays vetoed AND its remaining sensor-A (pre-2022) and sensor-C (post-2025) data get
SEPARATE offsets. Whitelist final: WISH, SCOT, KHBB, TRIN, LCCR, YACT, PETY (7). Chain restarted
2026-08-31 with this; map ERAS lines updated to the 7.

## Map + explorer changes (2026-08-28..31)
- borehole map: OSM tiles (CARTO began requiring API key); family id shown in title/slider/info;
  dotted red era lines (NOW filtered to the 6 whitelisted stations only); rebuilt 2026-08-31 with
  all vetoes → 183 stations (was 186), KBO 9, KHBB 12, TRIN 28, MOGU 17, SBAR 22, JRO 11, WLOO 21.
- `figures/interface_dvv_explorer.html` (NEW tool): map+series explorer of the corrected daily
  product; slider/play; interface (cells) vs site (stations) toggle; robust autoscale; color-scale
  selector with clipped labels; era lines per station (site mode) and per contributing station
  (cell mode). REBUILD PENDING after the rerun (splice new META/DATA, eras filtered to whitelist).

## POST-VETO RERUN RESULTS (chain completed 2026-08-31, lam_f 0.464 — adjacent grid step to old
## 0.215, stable across two oracle runs; some amplitude diffs are lambda-driven, remember in
## comparisons)
- assembler: 186 stns, 12,872 family-eras, 32.35M rows; date_veto 18,490; straddle 8,865 (7-stn split)
- solver: monotone 30 sweeps, honest dc metric ends 0.147 (no gauge-shuttle plateau — diagnostic fix
  working); products true-sign as before
- TILES old->new: sigma improved nearly everywhere (47.2: .0397->.0359; 48.2: .0364->.0326; 49.2:
  .1248->.1042); 49.9 worse (.376->.424). SAME survivors 47.2+48.2; merged 41.6+42.5 passes at
  .044 (was .047); **NEW: merged 48.2+49.2 now passes (.049)** — claimable deep coverage extends
  into S. Vancouver Island as a ~250-km merged tile. 180-d series corr(old,new) 0.93-0.99 on
  claimable tiles (cleanup didn't change the story); deep-index trend -0.034 -> -0.037 %/yr
  (the decline SURVIVES the QC cleanup — robustness point). Aggregate sigma .0168 -> .0153.
- **ALIASING PROBE now 0.75 (bar 0.8-1.2, was 0.83) — MARGINAL FAIL** -> per Merlin's
  pre-registered rule this is a SCOPE STATEMENT: tile amplitudes valid for structure at >=125-km
  scale aligned with the tiling; in-between structure under-recovered up to ~25% (one-sided LOW).
  Likely lambda-driven (0.464 smooths more). Consider Merlin ruling if it matters downstream;
  staggered second tiling is the remedy if ever needed.
- explorer rebuilt on new products (181 stations w/ coords, era lines at the 7 whitelisted).
  Movies + checkerboard_daily + methods-artifact figures still PRE-VETO (stale) — regenerate next.

## Rerun state
Full chain launched 2026-08-31 (assembler → solver → transmission_T → transmission_deconv),
log `logs/rerun_blocklist.log`. Map already rebuilt. Movies + explorer to regenerate after.
Products on disk BEFORE this chain completes still include vetoed families.
Expect: southern tiles (41.6/42.5) + mid-margin illumination (SEAS2 gone) to move most; the
before/after comparison on tiles + deep index is the post-chain deliverable.

## 2026-08-31 (late) — CATALOG-ALONE EXPLORATION (Thomas et al. 2025 prompt) [PROVISIONAL]
Discovery-member detections (33,960 after dedup/vetoes) tested for periodicity (Schuster):
- **lunar M2 (12.421 h): amp 26%, p~1e-248** — vs solar S2 (12.00 h) amp 9%. M2>>S2 is the classic
  genuine-tidal signature (cultural noise lives at 24 h [amp 38%!] and weekly [12%], neither of which
  can generate 12.42 h). Replicates known LFE tidal sensitivity IN OUR OWN CATALOG.
- **Along-margin gradient: N 34% / C 36% / S 13% M2 amplitude** — potentially interesting
  (along-margin variation in tidal sensitivity), NOT yet a claim.
- Inter-episode test: 100% of discovery detections lie within tremor days (by construction —
  PNSN-seeded discovery windows). Inter-ETS LFE activity question CANNOT be answered from the
  discovery members; would need the densify-era event lists (not retained) or a re-detection pass.
CONTROLS NEEDED before any claim (don't-over-conclude): (a) per-station detection-threshold vs
tidal-noise coupling (ocean noise is itself tidal); (b) seed-window inheritance (tremor is tidally
modulated; if discovery windows were sub-day, phase could be inherited — check window lengths);
(c) proper tidal-stress phase per location (not wall-clock folding); (d) Merlin audit.
Figure: figures/lfe_tidal_modulation.png; fronts figure: figures/lfe_fronts_allyears.png.
