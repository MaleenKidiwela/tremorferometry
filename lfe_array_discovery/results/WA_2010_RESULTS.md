# Washington N-Olympic 2010 — Bostock's published recipe, run blind (no catalogue truth)

Agent L, 2026-09-20.  The **same recipe and the same code** as the north reproduction (`results/REPRO_BOSTOCK.md` §0 parameter table, `results/NORTH_DISCOVERY_SUMMARY.md`): 6-s windows, 1–8 Hz, per-component cc averaged over **all** admitted stations, threshold = median + 8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s declustering, connected-component grouping at the 4-s node merge (≥3 members), Bungum–Husebye realignment before every stack, 12-s station-tapered templates cut 2 s before the stack envelope peak, network matched filter at 8 × MAD, **no waveform-similarity merge**.  Engine files `scripts/k_auto.py`, `k_group.py`, `k_grow.py`, `k_locate.py` are agent K's, unedited; the WA bed is supplied by `scripts/l_common.py` (station table, borehole azimuths, net-qualified cache index, station-hour health mask, day-shift null) and driven by `l_auto.py` / `l_run.py` / `l_locate.py` / `l_verify.py`.

**There is no Washington LFE catalogue.**  Bostock's ≥30 % one-to-one criterion cannot be applied, so nothing below is a rediscovery count and no family here is claimed to be real on the strength of a match.  What replaces truth is stated in §5: location with a shuffled-station control, repeat/co-timing statistics against their own base rates, and the identical chain run once on a site-shifted network.

## 1. The bed

* window **2010-07-20 … 2010-09-10**, sub-region N-Olympic (47.8–48.3 N, -123.7…-122.7).
* **96 three-component stations**, every probe-admitted station (`results/site_probe_2010.csv`) except the one that fails the station-level health bar (> 30 % vetoed station-hours): XG.DR09.  No station is withheld — the north's one real lever is station count (11 → 15 stations raised true pairs 8.4 % → 12.9 %, ledger §8).

| network | stations |
|---|---|
| CN | 1 |
| PB | 10 |
| UW | 2 |
| XG | 64 |
| XU | 19 |
| **total** | **96** |

* PNSN candidate (tremor) hours: **454** on **29** days (≥5 epicentres within 50 km of an active cell of the sub-region).
* hour-matched, earthquake-screened quiet hours: **40** (`results/wa_quiet_hours_2010.csv`, WA_V2_INPUTS §2.1).
* intraslab-earthquake S station terms available for **19** of the 96 stations; 77 carry 0.0 and are FLAGGED, never imputed (`results/x_station_corrections.csv`).
* per-station-hour health veto applied as a mask in every autocorrelation block (`results/wa_station_hour_veto_2010.csv`).

## 2. Pair stage (autocorrelation at 8 × MAD)

| | tremor hours | quiet hours |
|---|---|---|
| blocks run | 454 | 40 |
| blocks rejected (8×MAD bar not > 3 sd) | 0 | — |
| stations per block (median) | 96 | 78 |
| 8 × MAD in sd of the network function (median) | 4.90 | 5.20 |
| detections | 2,582,086 | 54,477 |
| detections per clean hour | 5,687.4 | 1,361.9 |

**Tremor-vs-quiet pair-count ratio at 8 × MAD = 4.18×** (top-2,000 tier: 4,405,286,343.61×).  The north's 11× is a *truth-based* ratio (12.85 % true pairs against a 1.17 % shifted-catalogue chance rate, ledger §8) and is **not the same quantity** — with no catalogue the only available contrast is tremor hours against catalogue-quiet hours, which is a weaker test because a quiet hour is quiet in tremor, not in noise.

8 × MAD sits at 4.90 sd of the network function, against 5.13–5.24 sd on every north hour — the rule is a measured bar, not a high one, on this bed too.

nmad distribution of the tremor-hour detections: p50 8.9, p90 11.0, p99 16.7, max 74.3.

| tier by nmad | nmad floor | hours represented | days |
|---|---|---|---|
| top 100 | 58.6 | 2 | 2 |
| top 500 | 37.0 | 5 | 5 |
| top 2,000 | 30.2 | 8 | 8 |
| top 10,000 | 20.5 | 63 | 22 |
| top 50,000 | 14.6 | 250 | 27 |

Every autocorrelation block is **one hour**, so both windows of every pair lie inside the same hour by construction and the pair-stage "fraction in a tremor-active hour" is 1.0 with no information in it.  The informative version — the fraction of **grown-template** detections that land in a PNSN tremor-active hour, against that quantity's own base rate — is in §5.

## 3. Grouping

* nodes (detections merged at the 4-s declustering interval): 45,142
* connected components: 4,990
* components with ≥3 members: **1,705**
* grown: **1,500** of 1,705 (growth budget 1,500; the GPU is shared, and the matched filter allocates one day-length accumulator row per template)

## 4. Growth (stack → template → matched filter, pass 0)

Pass 0 only: the north measured that iterating **loses** families (50 → 47 → 47 → 40, ledger §7) because an enriching template absorbs its neighbours' events, and the null must run the identical chain.

| | real chain | site-shift null |
|---|---|---|
| groups grown | 1,500 | 1,500 |
| templates | 1,499 | 1,500 |
| matched-filter detections | 204,403 | 26,377 |
| located (≥5 stations) | 1,498 | 1,500 |
| shuffled-station control passes | 356 | 82 |
| **CANDIDATE FAMILIES** | **50** | **17** |

## 5. Verification without truth

A grown template is a **CANDIDATE FAMILY** only if all of:

1. it locates on **≥5 stations** (stack-and-pick S, ≥4 MAD stack envelope peak) on the Crosson 1-D S travel times with the **intraslab-earthquake S station terms** of `results/x_station_corrections.csv`; a station without a usable term gets 0.0 and is counted, never imputed;
2. the **shuffled-station control passes** — the observed demeaned-residual rms beats the 5th percentile of 50 station-label permutations of its own picks;
3. **|Δz| ≤ 8 km** of the Slab2 interface under its epicentre;
4. **≥3 distinct days** of detections;
5. **tremor co-timing above the same-day base rate**: the fraction of its detections within ±30 min of a PNSN epicentre ≤60 km exceeds the fraction obtained by redrawing the same number of times uniformly on the same days (200 realisations).

**(c) The whole-chain site-shift null**, run once: each station's day series is replaced by the same station's data from a station-specific whole-day offset, restamped with the nominal day's `t0`, and the **identical** autocorrelation → grouping → stacking → matched filter → locator → screen is applied.  Inter-station coherence is destroyed; everything else — noise, glitches, diurnal structure, station geometry, thresholds — is preserved.

> **50 candidate families on the real network, 17 on the null.**

Fraction of grown detections inside a PNSN tremor-active hour: **0.567** against a base rate of 0.652 (null: 0.748 against 0.652).  The matched filter is run over the tremor days only, so this base rate is high by construction and the test is weak.

### Candidate families

| grp | lat | lon | depth km | Slab2 km | Δz km | stations | n_det | days | co-timing | base | rms s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 505 | 47.50 | -124.24 | 26 | 25.3 | +0.7 | 50 | 4281 | 17 | 0.025 | 0.017 | 4.57 |
| 629 | 47.50 | -124.08 | 28 | 26.2 | +1.8 | 25 | 506 | 13 | 0.468 | 0.134 | 4.52 |
| 1152 | 47.92 | -123.84 | 32 | 36.1 | -4.1 | 53 | 208 | 27 | 0.639 | 0.610 | 6.60 |
| 521 | 47.50 | -123.20 | 30 | 31.7 | -1.7 | 50 | 173 | 19 | 0.266 | 0.261 | 5.33 |
| 497 | 47.50 | -122.68 | 20 | 22.3 | -2.3 | 17 | 149 | 6 | 0.987 | 0.978 | 6.83 |
| 451 | 48.04 | -123.12 | 20 | 26.0 | -6.0 | 17 | 147 | 3 | 0.993 | 0.989 | 6.35 |
| 1040 | 48.79 | -123.20 | 20 | 25.1 | -5.1 | 23 | 144 | 3 | 0.993 | 0.986 | 3.25 |
| 519 | 47.53 | -122.84 | 20 | 15.4 | +4.6 | 38 | 143 | 15 | 0.741 | 0.721 | 4.72 |
| 438 | 47.53 | -122.76 | 20 | 21.4 | -1.4 | 28 | 137 | 10 | 0.985 | 0.957 | 3.48 |
| 1656 | 48.04 | -123.92 | 30 | 36.8 | -6.8 | 37 | 115 | 13 | 0.643 | 0.350 | 5.70 |
| 328 | 48.28 | -122.68 | 20 | 25.7 | -5.7 | 25 | 108 | 8 | 0.676 | 0.219 | 6.43 |
| 875 | 48.52 | -122.56 | 20 | 18.7 | +1.3 | 20 | 103 | 8 | 0.922 | 0.441 | 6.05 |
| 1029 | 48.13 | -123.12 | 20 | 13.6 | +6.4 | 16 | 98 | 4 | 1.000 | 0.980 | 4.66 |
| 1678 | 48.58 | -124.72 | 20 | 22.7 | -2.7 | 28 | 96 | 7 | 0.969 | 0.822 | 4.27 |
| 319 | 47.83 | -122.24 | 20 | 21.2 | -1.2 | 10 | 92 | 3 | 0.880 | 0.623 | 4.57 |
| 1685 | 48.79 | -123.48 | 20 | 18.2 | +1.8 | 20 | 83 | 4 | 1.000 | 1.000 | 5.84 |
| 1193 | 48.16 | -123.52 | 30 | 35.9 | -5.9 | 40 | 72 | 19 | 0.889 | 0.793 | 7.37 |
| 189 | 47.50 | -122.28 | 28 | 25.0 | +3.0 | 28 | 72 | 7 | 0.958 | 0.908 | 3.13 |
| 181 | 48.19 | -122.28 | 20 | 26.8 | -6.8 | 21 | 64 | 6 | 0.875 | 0.236 | 5.18 |
| 537 | 48.52 | -122.24 | 20 | 21.2 | -1.2 | 39 | 53 | 14 | 0.170 | 0.056 | 6.01 |
| 410 | 47.50 | -122.76 | 20 | 19.4 | +0.6 | 48 | 52 | 12 | 0.923 | 0.879 | 6.12 |
| 75 | 47.50 | -122.68 | 20 | 22.3 | -2.3 | 25 | 52 | 8 | 0.923 | 0.796 | 2.86 |
| 901 | 48.16 | -124.60 | 20 | 22.3 | -2.3 | 17 | 48 | 7 | 0.896 | 0.356 | 6.64 |
| 1399 | 47.50 | -122.56 | 20 | 23.3 | -3.3 | 18 | 46 | 10 | 0.370 | 0.268 | 7.50 |
| 1063 | 47.92 | -122.56 | 20 | 21.0 | -1.0 | 8 | 38 | 5 | 0.079 | 0.078 | 5.24 |
| 32 | 47.56 | -122.24 | 24 | 20.5 | +3.5 | 23 | 35 | 7 | 0.943 | 0.831 | 4.96 |
| 223 | 47.50 | -122.80 | 20 | 19.8 | +0.2 | 32 | 29 | 7 | 0.931 | 0.902 | 6.81 |
| 179 | 47.50 | -122.64 | 20 | 22.0 | -2.0 | 35 | 27 | 12 | 0.593 | 0.506 | 5.66 |
| 468 | 47.77 | -122.24 | 26 | 25.4 | +0.6 | 29 | 27 | 7 | 0.852 | 0.518 | 6.65 |
| 1035 | 47.50 | -123.00 | 20 | 16.5 | +3.5 | 27 | 27 | 13 | 0.333 | 0.278 | 6.72 |
| 355 | 47.50 | -122.76 | 20 | 19.4 | +0.6 | 51 | 23 | 7 | 0.913 | 0.877 | 3.55 |
| 1344 | 47.50 | -122.92 | 20 | 14.2 | +5.8 | 19 | 22 | 8 | 0.091 | 0.078 | 7.49 |
| 271 | 47.50 | -123.04 | 20 | 20.9 | -0.9 | 21 | 21 | 5 | 1.000 | 0.983 | 5.74 |
| 1476 | 48.79 | -123.52 | 20 | 24.9 | -4.9 | 21 | 18 | 5 | 0.722 | 0.712 | 6.36 |
| 1598 | 47.50 | -122.88 | 20 | 16.5 | +3.5 | 14 | 17 | 9 | 0.353 | 0.295 | 7.04 |
| 99 | 47.50 | -123.40 | 44 | 41.0 | +3.0 | 24 | 17 | 6 | 0.882 | 0.859 | 3.90 |
| 8 | 48.55 | -123.60 | 20 | 26.0 | -6.0 | 14 | 16 | 7 | 0.500 | 0.374 | 7.09 |
| 305 | 47.80 | -122.32 | 20 | 21.7 | -1.7 | 22 | 16 | 7 | 0.562 | 0.503 | 7.92 |
| 132 | 47.50 | -122.84 | 20 | 21.9 | -1.9 | 32 | 15 | 5 | 0.933 | 0.868 | 6.63 |
| 429 | 47.56 | -123.04 | 34 | 36.2 | -2.2 | 37 | 14 | 6 | 0.929 | 0.914 | 3.97 |
| 1613 | 47.50 | -122.92 | 20 | 14.2 | +5.8 | 21 | 12 | 7 | 0.250 | 0.229 | 4.65 |
| 369 | 47.59 | -122.80 | 20 | 19.3 | +0.7 | 25 | 10 | 5 | 0.800 | 0.780 | 6.32 |
| 569 | 48.22 | -122.96 | 20 | 27.3 | -7.3 | 30 | 10 | 4 | 0.900 | 0.880 | 7.13 |
| 1471 | 48.79 | -122.92 | 58 | 58.9 | -0.9 | 34 | 9 | 3 | 0.778 | 0.700 | 8.15 |
| 1094 | 47.50 | -124.00 | 20 | 27.0 | -7.0 | 23 | 9 | 4 | 0.111 | 0.062 | 4.09 |
| 1227 | 47.80 | -122.24 | 20 | 20.3 | -0.3 | 22 | 8 | 6 | 0.250 | 0.241 | 8.00 |
| 1506 | 47.77 | -122.32 | 20 | 23.1 | -3.1 | 25 | 8 | 4 | 0.250 | 0.168 | 6.50 |
| 844 | 48.28 | -124.80 | 32 | 24.2 | +7.8 | 37 | 7 | 4 | 0.143 | 0.113 | 6.98 |
| 984 | 48.01 | -124.16 | 20 | 26.1 | -6.1 | 18 | 7 | 6 | 0.714 | 0.565 | 5.73 |
| 486 | 47.86 | -124.16 | 20 | 26.9 | -6.9 | 12 | 6 | 4 | 0.500 | 0.259 | 6.52 |

The null produced 17 "family" (families) passing the same screen; they are listed in `results/l_wa2010_null/verify_it0.csv`.

### (d) Chestler comparison — descriptive only

`../codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv` is, as its filename says, **not a catalogue**; it is used here only to say where the candidates sit relative to previously published Olympic LFE coordinates, never as truth.

* candidates: median **32.8 km** to the nearest of 69 derived coordinates (closest 1.0 km);
* 1,000 equal-sized random draws from the sub-region's active cells: median 8.1 km (5th percentile 6.7 km);
* fraction of random draws at least as close: **1.000**.

## 6. Files

| | |
|---|---|
| pair stage | `results/l_wa2010/auto/`, `pair_all.parquet`, `pair_pairstats.json` |
| quiet control | `results/l_wa2010/auto_quiet/`, `pair_quiet.parquet` |
| groups | `results/l_wa2010/gb_nodes.csv`, `gb_summary.csv` |
| templates + detections | `results/l_wa2010/grow/tpl_it0.npz`, `det_it0.parquet` |
| locations | `results/l_wa2010/loc_it0.csv` |
| screen | `results/l_wa2010/verify_it0.csv` |
| site-shift null (whole chain) | `results/l_wa2010_null/` (offsets in `null_offsets.json`) |
| engine parameters as run | `results/l_wa2010/auto/_params.json` |

```json
{
 "win_s": 6.0,
 "hop_s": 0.5,
 "nmad": 8.0,
 "min_lag_s": 6.0,
 "decluster_s": 4.0,
 "n_stations": 96,
 "hours": "tremor",
 "n_blocks": 454,
 "nmad_in_sd_median": 4.8969144330472965
}
```

Null day offsets (station → whole days): 96 stations, range 1–52 days, seed 11.


## Addendum (orchestrator, 2026-09-21 07:45Z) — grid-edge reading of the candidate table
46 of the 50 real "candidates" and 14 of the 17 null "candidates" have locations pinned to a search-grid boundary (lat 47.50 or ≥48.75,
lon ≥ −122.30 or ≤ −124.75, depth 20 or 58 km), i.e. their locations are unconstrained and the Slab2 clause is meaningless there (same
artefact agent K flagged in the north at the 18-km floor). Interior candidates: **REAL 4 (grp 1152, 1656, 1193, 429) vs NULL 3.**
With ≥5 days and a tremor co-timing excess ≥0.05: **REAL 2 (grp 1656: 48.04 N, −123.92, 30 km, 13 days, co-timing 0.64 vs 0.35;
grp 1193: 48.16 N, −123.52, 30 km, 19 days, 0.89 vs 0.79) plus grp 1152 marginal (47.92 N, −123.84, 32 km, 27 days, 0.64 vs 0.61,
p 0.045) vs NULL 0.** The headline "50 vs 17" is therefore NOT the evidence; the evidence is 2–3 interior repeating sources on the
northern-Olympic interface against 0 in the null — suggestive, unconfirmed. Confirmation = re-fire of the unchanged templates on the
sealed 2011 window (with the 17 null templates as the false-confirmation control) and re-location on a wider grid.
