# NW event sharing between families — results

Run 2026-10-07, ~12:05–12:20 pm PT. Rules: `PLAN.md` (frozen before any pair count; one revision, recorded there, before
the episode re-run). Script: `lfe_extension_wa/scripts/nw_sharing.py` (sub-commands pairs, table, graph, array, loc,
episodes, figures). Log: `nw_sharing.log`. Nothing was deleted; no production file was edited or re-run
(`wa_post_dedup.py`, `wa_prod_post_v2.py`, `episode_rates.py` imported read-only). All results provisional.

"Scrambled control" = each station site's data shifted to a different time of the same day; "time-reversed control" = the
same detector with every template reversed in time. "Same-event pair" = the variant D pair condition: within ±2 s, the two
families share ≥ 50 % of their era-core template stations, and both detections predict the same arrival (≤ 0.5 s) at their
best shared station.

## The answer

1. **On the permanent network the 132 NW families behave as 93 distinct sources (83 among the 122 families with any
   2015–26 repeat).** 58 families fall into 19 multi-family "source groups" (12 pairs, 3 triples, one of 4, two of 5, one
   chain of 11); 74 are on their own. Linked pairs share almost everything: median sharing fraction 0.92 of the smaller
   family's repeats (minimum 0.34), against a time-shift chance level of at most 0.018. The sharing is bimodal (fig. 2):
   linked pairs sit at s ≈ 0.34–1.0, every pair with ≥ 10 shared events and no link sits below 0.3. 45 further pairs reach
   s ≥ 0.3 but share fewer than 10 events, so they are not linked.
2. **Most of these are the same source found twice by different discovery runs and located 10–140 km apart.**
   - 53 of the 58 links (91 %) join families from **different** union runs (base rate over all pairs: 72 %).
   - In the dense-array era (2007–11, 25-station 2-8 Hz stacks), 38 of 58 linked pairs are **array-confirmed copies**:
     the array stacks match on 5–25 stations with lag spread ≤ 0.15 s, and/or the two families' original discovery
     detections coincide (co-detection excess 0.1–0.94, median 0.28; matched-distance unlinked control: median 0.002,
     4 of 174 ≥ 0.1).
   - The NW union removed copies only for pairs **≤ 10 km apart** (`lfe_array_discovery/scripts/nw_union.py`). These copies
     sit a median 21.5 km apart (90 % > 10 km), so the union test never compared them.
   - By the array test, the 132 contain **100 distinct sources**: 22 sources have 2–4 copies (54 families), i.e. 32
     families are extra copies.
3. **The other 20 links are sources that are distinct on the array but not separable on the permanent network.** The
   array test calls them "neighbour" or "unclear" (lag spread 0.4–3.5 s) with no array co-detection (excess < 0.1). Most
   hang off the chain G001 (members 19, 46, 49, 67, 101, 122) and the pairs 11/69, 11/113, 113/127, 94/95, 95/115,
   66/115. Their shared events go through B013/B943.
4. **Why the permanent network merges them: one site dominates.** For 97 % of linked pairs the best shared station is
   B013 or B943 (one site; base 86 %). Weighted template-station Jaccard is 0.87 for linked pairs against 0.76 for all
   pairs. 95 % of grouped families carry B013/B943, against 73 % of singletons. The arrival check at one shared station
   proves the same event, not the same place.
5. **It is location error, often well beyond the quoted erh.**
   - Only 53 % of the array-confirmed copies lie within 2σ = 2·√(erh_i² + erh_j²) of each other. 29 % include a quality-C
     family; group C001 (fids 46/80/87/122) has erh 50/5/34/53 km and spans 136 km.
   - But 26 % are quality A on both sides. Examples: 66–95 are 43 km apart with erh 2.8/3.6 km and a 35 km depth
     difference; 75–126 are 39 km apart with erh 4.2/4.5 km and a 30 km depth difference; 25–95 are 51 km apart with
     erh 2.5/3.6 km. All three are waveform duplicates on 9–25 stations.
   - So these are one source with two Hyp2000 solutions, not neighbouring sources. erh is not a safe guide for this.
6. **Grouping follows neither location nor the catalogue flags.**
   - Linked pairs within 10 km: 7 % (base 4.5 %); within 20 km: 33 % (base 17 %).
   - Catalogue v1.2 `dup_group` (only cross-catalogue tests) touches 5 % of links (base 11 %).
   - `union_flag_lt2km` touches 29 % of links (base 24 %).
7. **The science result does not change; the station-coupling part shrinks.** NW E3 per-episode consistency at source-group
   level: 36 of 54 tested groups are over-dispersed (0.67), against 58/89 (0.65) at family level and 50/73 (0.68)
   de-duplicated. Median χ²/df is 4.3, against 4.03 and 3.67. The distance term gets stronger (−0.16, p = 0.001, against
   −0.11). The carrier-overlap term weakens: +0.17 (p = 0.018) against +0.23 (p = 0.001), and it is not significant at link
   s ≥ 0.5 (+0.10, p = 0.11). Part of the "partial station coupling" in the family-level result was copies of one source
   counted as two families.
8. **Rates.** Counting each event once inside a group removes 21 % of the NW headline repeats in every era:
   - E3 18,525 → 14,564; E2 6,311 → 5,006; E1 12,225 → 9,403;
   - array-confirmed groups: the same to within 0.1 %.
   - 93 % of all NW E3 variant-D flags (4,563 of 4,894) fall on the 58 linked pairs. The counted-once D column of today's
     dedup run is therefore essentially this grouping.
   - Per-family rates of a member are not wrong, but a group's members are one source (or one inseparable cluster).
     Firing rates, dv/v stacks and family counts should use the group: one representative, or the union of the members'
     events counted once.
   - Royer-truth matches: 46 families match a Royer NW template within 5 km. In 5 array-confirmed groups, two copies each
     match a different Royer template, so 41 distinct sources match.

## Numbers

### Source groups (E3, pre-registered link: s ≥ 0.30, s ≥ 5× chance, ≥ 10 shared events)
- Files: `source_groups_E3.csv` (per family), `source_groups_E3_summary.csv` (per group), `links_E3.csv`,
  `groups_summary.json`, `groups_sensitivity.csv`.

| era | link s | links | groups (all 132) | groups among active families | largest | families in multi-groups |
|---|---|---|---|---|---|---|
| E3 | 0.2 | 66 | 88 | 78 / 122 | 11 | 63 |
| **E3** | **0.3** | **58** | **93** | **83 / 122** | **11** | **58** |
| E3 | 0.5 | 47 | 97 | 87 / 122 | 6 | 56 |
| E2 | 0.3 | 37 | 100 | 90 / 122 | 4 | 55 |
| E1 | 0.3 | 51 | 92 | 91 / 131 | 5 | 66 |

- **Persistence across eras.** Of the 56 E3 links testable in E1 (≥ 10 shared), 40 also link in E1; in E2, 34 (51 have
  s ≥ 0.3).
- **The 11-family chain G001.** It is a chain, not a clique: 20 of its 55 member pairs are linked, and it spans 143 km.
  At s ≥ 0.5 it breaks up (largest group 6). By the array test it is C001 (46, 80, 87, 122) + C008 (19, 67) + C018
  (52, 120) + 3 singletons.

### Array-era tests (`array_era_tests.csv`, `array_era_tests_summary.csv`, `links_E3_classified.csv`)

| set | n | waveform duplicate | neighbour | unclear | median stack cc | median co-detection excess | excess > 0.2 | median km |
|---|---|---|---|---|---|---|---|---|
| linked | 58 | 33 | 20 | 5 | 0.40 | 0.28 | 36 | 27.6 |
| unlinked control (distance-matched) | 174 | 9 (all on 2–3 stations) | 43 | 122 | 0.17 | 0.002 | 4 | 27.7 |

Array-confirmed rule (PLAN Revision 1, set after seeing these numbers, diagnostic): co-detection excess ≥ 0.1 or waveform
duplicate on ≥ 5 stations → 38 same-source links, giving 100 groups (`source_groups_E3_array_confirmed.csv`,
`array_confirmed_groups_summary.csv`). Copies lie 11.6 / 25.5 / 50.7 km apart (10 / 50 / 90 % of the groups' maximum
separation).

**Surprise.** Four unlinked control pairs (7–61, 13–14, 11–112, 95–116) co-detect strongly in the array era (excess
0.44–0.62) but have no waveform match and share no events on the permanent network. These are not investigated here.

### Location vs sharing (`location_vs_sharing.json`, fig. 3)

| pairs | n | median km | within 2σ_h | median σ_h | median max erh | median min picks | any quality C | both A |
|---|---|---|---|---|---|---|---|---|
| linked, array-confirmed copies | 38 | 21.5 | 0.53 | 9.7 | 7.5 | 9 | 0.29 | 0.26 |
| linked, distinct on array | 20 | 46.2 | 0.40 | 14.7 | 14.2 | 6 | 0.45 | 0.05 |
| unlinked (N_small ≥ 10) | 5,093 | 39.9 | 0.27 | 8.8 | 7.7 | 8 | 0.29 | 0.30 |

- Pick counts are not low for the copies: median of the smaller count 9; 71 % have ≥ 8 picks on both families; minimum 4.
- Copies from the same run: 0 of 38. Every copy pair comes from two different runs.

### Per-episode consistency, NW E3 (k-free share test, 47 MAJOR episodes)
- Files: `episode_consistency_groups.json`, `episode_consistency_groups.csv`, `episode_residual_pairs_groups.csv`.
- Acceptance check: singleton groups reproduce the family-level V0 run exactly (89 tested, 58 over-dispersed, median
  χ²/df 4.034, distance −0.111, carrier +0.230).

| row | units | tested | over-dispersed (frac) | median χ²/df | distance coef (p) | carrier coef (p) | residual corr <3 / 3–6 / 6–10 / 10–20 / 20–40 km | Σ N |
|---|---|---|---|---|---|---|---|---|
| family level (V0) | 132 | 89 | 58 (0.65) | 4.03 | −0.111 (0.001) | +0.230 (0.001) | 0.31 / 0.11 / 0.16 / 0.06 / −0.03 | 15,900 |
| de-duplicated (S0b, keep ±2 s) | 132 | 73 | 50 (0.68) | 3.67 | −0.116 (0.001) | +0.269 (0.001) | 0.31 / 0.13 / 0.17 / 0.06 / −0.04 | – |
| **source groups (pre-registered)** | **93** | **54** | **36 (0.67)** | **4.30** | **−0.162 (0.001)** | **+0.166 (0.018)** | 0.37 / 0.12 / 0.20 / 0.08 / −0.05 | 13,011 |
| source groups, link s ≥ 0.2 | 88 | 49 | 33 (0.67) | 4.57 | −0.167 (0.001) | +0.161 (0.038) | 0.42 / 0.10 / 0.14 / 0.09 / −0.03 | 12,988 |
| source groups, link s ≥ 0.5 | 97 | 58 | 39 (0.67) | 4.07 | −0.162 (0.001) | +0.105 (0.107) | 0.48 / 0.11 / 0.18 / 0.07 / −0.04 | 13,015 |
| array-confirmed groups | 100 | 61 | 41 (0.67) | 3.97 | −0.135 (0.001) | +0.127 (0.047) | 0.44 / 0.18 / 0.13 / 0.06 / −0.04 | 12,980 |

- Group count = member headline repeats at P ≥ 0.75 minus the time-reversed control. A member's repeat is dropped if
  another member has a higher-cc repeat within ±2 s with the same arrival.
- Window = the union of the members' windows; position = the members' mean.
- Few pairs per distance bin at group level (10 at < 3 km), so the bin means are noisy.

### Counted once (`counted_once_by_group.json`)
NW headline data repeats with each event counted once per group:

| era | raw | counted once | dropped |
|---|---|---|---|
| E3 | 18,525 | 14,564 | 21.4 % |
| E2 | 6,311 | 5,006 | 20.7 % |
| E1 | 12,225 | 9,403 | 23.1 % |

Variant-D winners of NW E3 headline rows: 4,894 NW, 23 Olympic A, 3 Arm B.

## Checks
- Headline rows used here = `episode_rates/counts_wa_nw_fid_date_half.parquet` n075 for data (37,061) and time-reversed
  control (34), 0 per-family differences. The scrambled control was not needed here and was not extracted.
  (`check_headline_vs_episode_rates.json`.)
- Pair counts are symmetric in meaning. The directional variant-D winner counts from `post_v2/dedup/flags/` are columns
  (`n_Dflag_i_by_j`, `n_Dflag_j_by_i`).

## Caveats
- The link rule uses one shared station's arrival. It proves the same event, not the same location (point 4).
- The array test uses each run's own stacks and detections. Different runs searched partly different days (the
  common-day count is a column).
- The Revision 1 thresholds were chosen after seeing the linked-pair numbers. The pre-registered answer is the 93
  groups; the 100 is a diagnostic.
- Connected components chain: G001 links families that are not each other's copies.
- At group level fewer units are testable (54 vs 89), so the over-dispersed fraction has a wider error.

## Figures
- `fig1_sharing_map.png` — NW families with E3 sharing links, coloured and weighted by s (s ≥ 0.3 bold), multi-family
  groups coloured.
- `fig2_sharing_vs_distance.png` — s vs distance with the chance level; binned mean s, chance level and linked fraction.
- `fig3_separation_vs_erh.png` — pair separation vs combined erh, linked vs unlinked.

## Files
Everything is in `lfe_extension_wa/results/nw_sharing/`:
- `pair_table.csv` — 8,646 pairs × 3 eras, all columns of PLAN §1;
- `nw_families_meta.csv` — run, erh / erz / picks, Royer match;
- `pair_counts_raw.parquet` — per pair, shift and era;
- `nw_headline_rows.parquet` — NW headline data and time-reversed rows;
- `dflag_winner_counts.csv`.
