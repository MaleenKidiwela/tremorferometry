# Static Vp/Vs pilot from LFE P and S picks (per-family Wadati κ) — PLAN

Written 2026-10-12 ~7:33 pm PT, BEFORE any κ was computed. Spec: `../MERLIN_TOMO.md` (§4 pilot). Scripts
`svi_science/scripts/s11_tomo_*.py`; tables `tables/`, figures `figures/`, logs `logs/`.

## 0. Data sets (all families kept; quality and class are columns, never filters)
- **A — spec set:** `lfe_array_discovery/results/ab_royer/NW/hypo_28/picks.csv` (run 1, 169 groups × 25 stations).
  Column `eye_class` from `NW/hypo_28/review/eye_class.csv` (33 `lfe`, 136 `noise`). The noise groups are a free
  second control: their κ should not look like the real families'.
- **B — catalogue NW:** the 141 `Royer_NW_WA` families of `data/catalogue/families_v1_3.csv`. Picks recomputed with the
  frozen picker `a4f_hypo.make_picks` (Hilbert rotation, zero-phase extremum, ratio rP/rS; S = largest horizontal
  extremum, P = largest vertical extremum 1.5–14 s before S) on each family's `stack_file` (hypo_28 stacks).
  Sanity: the same code on `NW/hypo_28/stacks_loc.npz` must reproduce set A's picks.csv exactly.
- **C — catalogue SW (step 2):** the 81 `Royer_SW_WA` families, same recipe (their SW_mute/SW_lag15 hypo_28 stacks;
  the SW_mute/hypo and SW_lag15/hypo picks.csv on disk are the older all-station runs, kept as a column).

## 1. Per-family κ (MERLIN_TOMO §4.1)
- Station rows with finite tP and tS and rP ≥ 1.2 AND rS ≥ 1.2; family measured when ≥ 6 such stations.
- OLS (tS − tP) = a + b·tP; **κ = b + 1**; rms = residual RMS (s).
- Bootstrap over stations (1,000 resamples with replacement, refit; draws with < 3 distinct stations skipped);
  σκ = SD of the bootstrap κ (16–84 % half-width as a column).
- Columns (not used for the decision): κ from OLS tS on tP; κ from orthogonal (total least squares) fit of tS on tP
  (the OLS form has a correlated-error bias because tP noise enters both axes); κ at rP, rS ≥ 1.5.
- **Acceptance (pilot gate G1):** ≥ 50 families with σκ ≤ 0.03 AND rms ≤ 0.25 s ("good" families).

## 2. Controls
- **G2 scrambled control (pairing scramble).** Within each family, the station labels of the S picks are permuted
  relative to the P picks (each S time paired with another station's P time), then κ is refit; 1,000 scrambles.
  (Permuting whole (tP, tS) rows between stations leaves a Wadati fit unchanged, so the S-to-P pairing is the
  only station-label scramble that tests anything here.)
  - G2a (literal §4.2): between-family spread of κ = 1.4826·MAD over the good families; real vs the scrambled
    control's distribution of the same statistic over the same families (one scramble per family per draw,
    1,000 draws); literal pass = real spread above the control's 95th percentile. Written down now: the scramble
    breaks the P–S pairing, so its κ collapse towards ~0 with a lever-arm-dependent scatter; the literal
    comparison may be ill-posed in direction, and is reported as stated.
  - G2b (information): per family, p = fraction of scrambles whose rms ≤ the real rms. Pass if ≥ 80 % of good
    families have p < 0.05, and the good-family count in the scrambled control (σκ ≤ 0.03, rms ≤ 0.25 s) is
    ≤ 5 % of the real count.
  - G2c (between-family spread is signal, not noise): χ²/ν of good-family κ about their weighted mean with the
    bootstrap σκ; pass if p(χ²) < 0.05. Column: split-half reproducibility for families with ≥ 10 stations
    (random station halves, corr(κ_half1, κ_half2) over families, 200 splits) vs the same in the scrambled control.
- **G3 leave-one-station-out:** per good family max |Δκ| over dropped stations; pass if ≥ 80 % of good families
  have max |Δκ| ≤ 0.05.
- **G4 site-term confounding:** (i) joint fit (tS − tP) = a_f + b_f·tP + c_s (station terms, sum-to-zero) over all
  measured families; pass if median |κ_joint − κ| ≤ 0.02 over good families and Spearman(κ_joint, κ) ≥ 0.8;
  (ii) regress κ of good families on station-presence indicators (ridge-free OLS on stations used by ≥ 5 good
  families); R² compared with the R² for κ randomly permuted across families (1,000); pass if real R² is below the
  permuted 95th percentile (the station set does not explain κ).
- **Fail signature (infeasible as is):** rms > 0.3 s typical (median over measured families) or G3 fails.
- **Adaptation (run if G1 or the fail signature fails, before any conclusion):** stack the vertical stacks of
  families within 3 km (catalogue epicentres; single-linkage clusters at 3 km), each station's family stacks
  aligned on that family's S pick (rS ≥ 1.2) before averaging; S is re-picked on the aligned horizontal stack and
  P on the stacked vertical with the same picker; κ per cluster with the same rules and gates. The S-alignment
  is an amendment-free part of the adaptation (the alignment uses only S, so it cannot manufacture a P–S slope;
  the scrambled control is rerun on the clusters).

## 3. Mapping (columns; no pass rule)
- κ vs catalogue depth (Hyp2000/catalogue `depth_km`; Slab2 depth as a column), vs along-strike distance
  (projection on azimuth 330° through 47.0 N, 123.0 W), vs firing rate (n_detections / era days).
- Spearman ρ with a permutation p (10,000) for each. ETS segmentation: all NW and SW WA families lie inside one
  Brudzinski & Allen (2007) segment (northern, ~14 months); no finer segmentation table is on disk → only NW vs SW
  is compared. Merrill et al. 2020 (if a table is found) as a column.

## 4. SVI (step 3)
- Provenance of `sp` in the template npz (`bw_gate_common.py:118–130` reads `z['sp']`): trace to its writer.
- P measured on the vertical, S on the horizontals of stacked family templates (boreholes B009/B010/B011/B926
  first, then core and POLARIS stations), same picker as WA, then κ with the same §1–§2 rules. Stack sources:
  S10 pilot/multi windows and grand stacks; exact recipe written as an amendment before the first SVI κ.

## 5. Step 4 (only if G1–G3 pass): plan, do not run, the earthquake-pick fetch (USGS phase-data / NRCan
`includearrivals`) for a differential Wadati of the 25–45 km layer per ~30-km zone; count requests and time.

## Amendments
(time-stamped below; none before the first number)

**Amendment 1 (2026-10-12 ~7:37 pm PT, after the first numbers: G1 fails in A/B/C; median rms 1.4–1.7 s, σκ 0.26–0.61;
only 17/61/40 families have ≥ 6 stations with rP, rS ≥ 1.2).** Diagnosis from record sections
(`figures/diag_recsec_A.png`, 4 eye-real families): P and S are visible and consistent at stations ≲ 45 km;
the rms is driven by a few far-station picks on noise. Added BEFORE computing them, as variants (the §1 OLS result
stays the headline of the spec test):
- V1 robust: iterative clipping of stations with |residual| > max(3 × 1.4826·MAD, 0.15 s), refit, until stable;
  family measured if ≥ 6 stations remain; bootstrap, LOSO and scrambled control on the kept stations.
- V2 distance: stations ≤ 60 km epicentral distance from the family location (catalogue / eye location), then V1.
- The adaptation (§2, cluster stacks aligned on S) runs regardless, with OLS and V1/V2 columns.
Same gates G1–G4 for every variant; a variant is only "the" result if it passes G1–G3 AND its scrambled control
gives ≤ 5 % as many good families.

**Amendment 2 (2026-10-12 ~7:40 pm PT, after V1/V2: still G1 fail, rms 0.67–1.07 s after clipping).** The
records show why: well-picked stations are the near ones (≲ 35 km), where tP spans only ~0.5–1 s for a 35–45 km
deep source, so the Wadati lever arm is tiny; far stations (40–90 km) carry S but rarely a clean P. The
adaptation therefore targets far-station P. Implementation (`s11_tomo_cluster.py`), fixed now:
- clusters = single linkage of catalogue epicentres at 3 km (5 km column), sets B and C separately;
- each member stack (3 × 44 s) unit-normalised per station; origin alignment by a robust two-way fit of the
  confident S picks (rS ≥ 1.2) within the cluster, tS(f, s) = τ_s + δ_f → member f shifted by −δ_f at every
  station (S1: preserves the moveout, reaches stations where that member has no S pick); S2 column = per-station
  S alignment where the member has rS ≥ 1.2 at that station (else the S1 shift);
- cluster stack = mean of the aligned members; picks with the frozen picker; κ by OLS and V1/V2, same gates;
- clusters with one member are carried (identical to the family) and flagged `n_members`.

**Amendment 3 (2026-10-12 ~7:44 pm PT, after the cluster adaptation: no gain; 3-km clusters are mostly
singletons (B 79/106, C 47/54), and the 22-member 5-km cluster stack shows P only ≲ 10 km and S ≲ 30 km;
cluster rms 0.8–2.1 s).** Diagnosis written down: pick precision at the best stations is ~0.1–0.15 s (co-located
B013/B943 differ by 0.15–0.18 s on SNR 65 stacks: a zero-phase extremum of a ~4 Hz pulse jumps by half a cycle)
and the usable tP span per family is 1–2 s, so per-family σκ ≥ ~0.08 is a floor; G1 (σκ ≤ 0.03 per family) is out
of reach for these picks whatever the variant. Next estimator, fixed before computing it:
- **Pooled zone Wadati:** (tS − tP)_fs = a_f + (κ_z − 1)·tP_fs with one κ per zone, family intercepts a_f;
  rows rP, rS ≥ 1.2; a family enters with ≥ 4 such stations; V1 clipping on the pooled residuals (per family
  median-centred); zone measured with ≥ 5 families. Zones: whole NW (B), whole SW (C), 30-km along-strike bins,
  depth bins (< 30, 30–36, 36–42, > 42 km). σκ by bootstrap over families (1,000) and, as a column, over stations.
- Controls: pairing-scrambled control (S picks permuted across stations within each family, 1,000; its pooled κ
  distribution); eye-noise groups of set A as a second control zone; leave-one-station-out over the zone
  (|Δκ| ≤ 0.05); joint station terms (|Δκ| ≤ 0.02).
- Pass for a zone κ: σκ ≤ 0.05, |real − scrambled median| > 4 σ of the scrambled distribution, LOSO and
  station-term rules met. κ outside 1.4–2.6 is flagged unphysical (column).

(Time correction, ~7:49 pm PT: the first versions of this file carried clock times 7:40/7:58/8:05/8:20 pm that were
guesses, not read from the clock; replaced by approximate times (Amendments 1–3 lie between 7:33 and ~7:40 pm; ±10 min). Order of events is exact.)

**Amendment 4 (2026-10-12 ~7:49 pm PT, after the pooled zone run).** Pooled κ is 0.66–1.10 in almost every zone
(physically impossible: S−P must grow with tP), with a scrambled control at ~0.0. This is the errors-in-variables
signature: P-pick noise enters both axes of (tS − tP) vs tP, so κ_obs ≈ κ·R with R = var(tP_true)/var(tP_obs);
κ ≈ 0.8 implies P noise comparable to the 1–2 s tP span. Two SW zones (along 0–30 km, depth 36–42) give
κ 1.85 ± 0.01–0.02 with rms 0.26–0.35 s after clipping (station-term rule fails, Δ −0.08/−0.11). Next, fixed now:
- Pick-error variances measured from co-located stations (same site < 1 km: B013/B943/W030, B007/B04A/LRIV,
  DOSE/C04A, W060/SMW and any SW pairs < 1 km): σ_P² = var(ΔtP)/2, σ_S² = var(ΔtS)/2 (robust: 1.4826·MAD),
  on rows with rP, rS ≥ 1.2.
- Estimators (pooled with family intercepts, same zones, controls, gates): E1 inverse regression tP on tS,
  κ = 1/slope (P noise in the dependent variable does not bias it; S noise does, by R_S); E2 Deming regression of
  tS on tP with δ = σ_S²/σ_P² from the co-located pairs; E3 = E1 corrected by the measured R_S. The original OLS
  stays in the table. Pass rules as Amendment 3; in addition E1 and E2 must agree within 2σ.

**Amendment 5 (2026-10-12 ~7:51 pm PT, after `tables/colocated_pick_errors.csv` and `figures/diag_sp_vs_distance.png`).**
Co-located pick errors are small (robust σ_P 0.01–0.09 s, σ_S 0.03–0.16 s at the B013/B943/W030 site) but with
gross outliers (6–50 % of pairs differ by > 0.5 s), so Gaussian errors-in-variables cannot explain κ < 1. The
distance plot shows the cause: a band of rows at tS − tP ≈ 1.5–2.5 s at all distances — the picker's P is the
largest vertical extremum 1.5–14 s before S, and it often lands on S energy arriving on the vertical just before
the S pick (plus early noise picks at tS − tP 10–13 s at short range). Both bias the Wadati slope down.
Adaptation (V3 picker, κ-free, fixed now): P characteristic function = |Hilbert-rotated Z| × w(t), with
w = E_Z/(E_Z + E_H) from 0.5-s running energies of the unrotated stack (P is vertically polarised at these
incidence angles, S-on-vertical leakage comes with larger horizontal energy); pick = maximum of the CF in
[tS − 14, tS − 1.5] s; rP = CF maximum / largest CF outside ±0.5 s. S unchanged. Same §1 per-family, Amendment 3
pooled zones (OLS + E1 inverse regression), all controls; run on sets A, B, C and the 3-km clusters.
The E1/E2 estimators of Amendment 4 are computed as columns in the pooled script on all pick sets.

**Amendment 6 — SVI recipe (2026-10-12 7:52 pm PT, before any SVI pick).**
- `sp` provenance (verified in code): `lfe_catalog_repair/scripts/e5r_templates.py:172` writes
  `SP = (ttS − ttP) + (S pick − nominal) − 0.5 s`, with ttS, ttP from the Crosson 1-D model at the catalogue
  location; `bw_gate_g3_templates.py` reuses the same builder (`T.do_station`) for the era/POLARIS templates.
  So `sp` is a MODELLED S−P plus a measured S-pick refinement and a placement constant; it contains no P
  measurement and is not used for κ (it would return the model's Vp/Vs).
- SVI-1: `lfe_catalog_repair/results/e5r_templates_2012.npz` stacks (116 families × 15 stations, 2012, 1–8 Hz,
  20 Hz, 40-s windows centred on the catalogue-origin + model ttS; per-event cc alignment shifts are ≈ 0,
  `e5r_template_qc_2012.csv`). Time frame per station: t = ttS − 20 s + i/20 relative to the catalogue origin
  (the model only places the window; picks are measured inside it). Upsampled to 40 Hz (resample_poly 2/1) so
  the frozen picker runs unchanged: S = frozen picker, P = V3 (headline) and v0 (column). All §1 per-family and
  Amendment 3–5 pooled analyses (OLS, E1) and controls; locations from `e5r_families_2012.csv`.
- SVI-2: S10 pilot grand stacks (`s10_dvv/pilot/tables/grandstacks.npz`, 3 families × up to 21 stations incl.
  POLARIS 2003–06; tau −14…+10 s at 40 Hz relative to the cut offset `offsets_sp.csv:offset_cut`); time frame
  offset_cut + tau; per (family, station) the epoch with most rows; same pickers; per-family κ only.

**Amendment 7 (2026-10-12 7:57 pm PT, after the pooled E1 runs).** E1 (tP on tS) with clipping gives
κ ≈ 1.85–1.96 at rms 0.16–0.22 s wherever the clipping converges on a coherent subset (NW all 1.95 ± 0.03,
SW all 1.90 ± 0.015, eye-real NW 1.89–1.91), but in other zones it converges on a mixture (rms ~2 s, κ 1.1–1.3).
The clipped answer depends on where the iteration starts, so a start-free estimator is added, fixed now:
- **Inlier profile:** for κ on a grid 1.30–2.60 (step 0.005), c = 1/κ, family intercept = median over its rows of
  (tP − c·tS), residual r = tP − a_f − c·tS; score(κ) = number of rows with |r| ≤ 0.25 s (families with ≥ 4
  rows). κ̂ = arg max; σ by bootstrap over families (1,000); scrambled control = same profile after the
  pairing scramble (200 draws): peak score of the data must exceed the control's 99th percentile; LOSO over
  stations |Δκ| ≤ 0.05; joint station terms on the inlier rows |Δκ| ≤ 0.02 (as Amendment 3).
- On the inlier rows at κ̂: E1 OLS, E3 = E1 corrected for S-pick noise (κ·R_S with R_S = 1 − σ_S²/var_f(tS),
  σ_S from co-located pairs), E2 Deming (δ = σ_S²/σ_P²). Pass also needs E1/E2/E3 within 2σ of κ̂.
- Zones as Amendment 3; pick sets AV3 (eye lfe/noise), BV3, CV3, S1V3, and v0 picks A, B, C as columns.

**Amendment 8 (2026-10-12 8:04 pm PT, coordinator request after RESULTS v1; before any new stack or pick).**
Clean P exists only within ~35 km, so use the SVI 2003–2006 era, when POLARIS stations sat 1–9 km from many families.
- Detections: `data/catalogue/svi_catalogue_v1_4_backward.parquet`, series data, counted_once, P ≥ 0.75,
  date ≤ 2006-06-30 (39,522 rows, 125 families, 391 days; all era detector). Family = e5r bkey; locations from
  `e5r_families_2012.csv`.
- Stations (20-Hz 1–8 Hz caches `lfe_extension/cache20hz_backward`): C8 KHVB/MGCB/SHDB/SHVB/TWBB,
  PO KELB/KLNB/SILB/SSIB/TSJB/TWGB/TWKB, CN LZB/PGC/VGZ/SNB/NLLB, PB B001/B003–B007/B009–B011 (2005–06 only).
- Stacks: per (family, station) linear mean of 3C windows [t + ttS − 16, t + ttS + 8] s (ttS Crosson 1-D from the
  catalogue location: it only places the window; the time frame is t + offset + i/fs, common to all stations of a
  family); each window RMS-normalised; windows with RMS > 5 × the (family, station) median dropped. Upsampled
  to 40 Hz; S = frozen picker, P = V3 picker.
- **Headline rows: stations with hypocentral distance ≤ 40 km (catalogue location; geometry selects rows only).**
  Columns: ≤ 60 km and all stations. Per-family κ (§1, OLS and tP-on-tS, clipped), pooled zone (Amendment 3,
  tP-on-tS) and inlier profile (Amendment 7); families enter the pooled fits with ≥ 3 rows in the ≤ 40 km set
  (fewer stations exist there). Same controls (scrambled control, LOSO, station terms) and pass rules.
- WA borehole-only column: PB boreholes (B0xx, B9xx) within 35 km epicentral, from picks_BV3/CV3 and B/C;
  pooled tP-on-tS and inlier profile with ≥ 3 rows per family.
- The earthquake-pick fetch stays NOT run (user decides).

## 6. Earthquake picks and the differential Wadati (written 2026-10-13 ~9:35 am PT, BEFORE any pick was fetched)
User OK 2026-10-13 ("do the earthquake download"). Workspace `tomo/eq/` (RESULTS.md, figures, tables, logs).
Scripts `svi_science/scripts/s11_eq_{select,fetch,parse,wadati,diff}.py`, runner `s11_eq_fetch.sh`.
Raw cache `data/eq_picks/raw/<source>/`, parsed picks `data/eq_picks/picks.parquet`.

**6.1 Selection** (`s11_eq_select.py` → `eq/selection.csv`, `eq/selection_counts.csv`). Source table
`g2_local/tables/eq_local.parquet`, rows `type == earthquake` and `removed == ""` (blasts, volcano circles out;
ComCat/NRCan duplicates are already merged: `dup_of` holds the other catalogue's id). LFE area = locations of the
families used for LFE κ (picks_BV3 = WA NW, picks_CV3 = WA SW, picks_EallV3 = SVI era; families deeper than 20 km).
Each event gets the region of its nearest family and a 30-km along-strike zone (azimuth 330° through 47.0 N,
123.0 W, the LFE zone coordinate).
- Tier 1 crustal (cls crustal, depth ≤ 25 km), epicentre ≤ 50 km from a family: SVI M ≥ 1.0 (all);
  WA M ≥ 1.5, at most 120 per (region, zone), largest magnitude first (keeps zones balanced).
- Tier 1 intraslab (cls intraslab, i.e. below Slab2 + 5 km), epicentre ≤ 50 km from a family, M ≥ 1.0. Their
  rays to stations above the LFE area cross the LFE source layer.
- Tier 2 layer events (25–45 km, overriding_deep / near_interface), epicentre ≤ 30 km, M ≥ 1.5, at most 50 per
  (region, zone): a column only (same depth range as the LFEs, with a known origin time). Fetched last.
- Download order: tier, then rank within zone (round robin over zones), so a partial download stays balanced.
- **Counts:** 2,692 events. ComCat-reachable: tier 1 **1,884** (crustal 1,328 incl. SVI 242; intraslab 556 incl.
  SVI 87), tier 2 389. 419 NRCan-only events (mostly SVI) have no public pick source (below). Per-zone counts in
  `eq/selection_counts.csv`: WA NW zones −30…180 km 13–120 crustal and 3–94 intraslab each; WA SW −120…60 km
  19–120 crustal and 4–42 intraslab; SVI 120–270 km 12–103 crustal (ComCat) and 2–23 intraslab.

**6.2 Sources (tested 2026-10-13 ~9:20 am PT).**
- USGS ComCat: event detail (`fdsnws/event/1/query?eventid=…&format=geojson`) → `phase-data` product QuakeML
  (picks + arrivals with phase, distance, timeResidual, timeWeight; pick time uncertainty). Prefer the `uw`
  (PNSN) product, else the first listed. ComCat `uw` phase-data IS the PNSN AQMS pick set, so a separate PNSN
  fetch adds nothing. For NRCan-primary rows the merged ComCat twin (`dup_of` = uw/us id) is fetched.
- **NRCan FDSN does not return arrivals:** `includearrivals=true` (also with `includeallorigins`) returns origins
  and magnitudes only (3 queries incl. 2024 M ≥ 3.5 events; 0 picks); `eventid` must be YYYYMMDD.HHmm. The
  NRCan web search offers HTML/text/QuakeML of the same origin-only service. So NRCan-only events are listed as
  `none_public` and not fetched. Next adaptation if SVI is short: pick P and S on CN/PO waveforms ourselves.
- Requests: 2 per event, ~4,550 in all. Response time observed 0.3–16 s, with intermittent 502 from CloudFront.
  ≤ 2 parallel workers, ≥ 0.5 s between a worker's requests, retry 6× with backoff 4, 8, 16, 32, 64, 128 s;
  a failure writes `<id>.fail` and is retried on the next run. Cached responses are skipped. Estimate ~1.5–3 h.

**6.3 Parsing** (`s11_eq_parse.py`, run with `python -I`). One row per (event, station, phase): event id,
origin (time, lat, lon, depth from the QuakeML preferred origin; catalogue values as columns), net, sta, channel,
phase, pick time, tP/tS relative to origin, uncertainty, timeWeight, timeResidual, distance (km), azimuth,
evaluation mode. Station key = net.sta (P on Z, S on horizontals may use different channels).

**6.4 Earthquake κ per zone** (`s11_eq_wadati.py`). Rows = stations with both P and S of one event.
- Headline rows: epicentral distance ≤ 100 km (direct crustal phases; Pn crosses over beyond ~100–150 km),
  all pick weights. Columns: ≤ 60 km, ≤ 150 km, timeWeight > 0 only. No filter on catalogue residuals (they
  carry the locator's own Vp/Vs).
- An event enters a zone fit with ≥ 4 such stations; a zone is measured with ≥ 5 events.
- Estimators, exactly as the LFE κ: pooled fit with one κ per zone and per-event intercepts (the intercept absorbs
  origin-time error), E1 tP on tS with V1 clipping (headline, as the SVI-era LFE value), OLS (tS − tP) on tP
  (column), and the start-free inlier profile (Amendment 7; threshold 0.25 s; column 0.5 s since catalogue
  picks of crustal events at 100 km carry larger scatter).
- Gates, same as Amendments 3 and 7: σκ = max(event bootstrap, station bootstrap) ≤ 0.05; |z| vs the pairing-
  scrambled control (S picks permuted across stations within each event, 1,000 draws) > 4; leave-one-station-out
  |Δκ| ≤ 0.05; joint station terms |Δκ| ≤ 0.02; profile: peak above the scrambled control's 99th percentile and
  E1/E2/E3 within 2σ. Every zone is reported with its gate columns, passing or not.
- Classes per zone: crustal (headline for the differential), intraslab (cross-check), layer (tier 2 column).

**6.5 Differential κ of the layer between 25 km and the LFE source** (`s11_eq_diff.py`).
- LFE κ inputs (from `tables/`, no new LFE fits): SVI era clipped tP-on-tS 1.773 ± 0.010; WA SW 1.904 ± 0.015
  (clipped) and 1.91 ± 0.024 (profile); WA NW 1.95 ± 0.029 (V3 clipped) with 1.815–1.85 (frozen P) as the range.
  Per-zone LFE values from `pooled_*_inv.csv` / `profile_*_zones.csv` as columns, with their gate status.
- Headline (ray-weighted): 1-D Vp from the CVM v1.7 profile of the region (`south_VI.csv`, `royer_NW_WA.csv`,
  `royer_SW_WA.csv`). P rays traced (flat layers) from the zone's median LFE depth to the epicentral distances of
  the zone's LFE inlier rows. Model tS = Σ_layers κ_l · tP_l with κ = κ_crust above 25 km (the zone's crustal
  earthquake κ) and κ_layer from 25 km to the source. The same pooled tP-on-tS slope on the model times gives
  κ_LFE,model(κ_layer); κ_layer solves κ_LFE,model = κ_LFE,obs. Uncertainty: 2,000 draws of κ_LFE and κ_crust
  from their σ (Gaussian), 16–84 % range.
- Column (simple layered): f = P travel-time fraction below 25 km on the vertical ray;
  κ_layer = (κ_LFE − (1 − f)·κ_crust)/f.
- Zones: whole region first (where LFE κ passed), then 30-km zones where both an LFE and a crustal κ exist.
- Cross-check: predicted intraslab Wadati κ from (κ_crust, κ_layer, κ below the source 1.73–1.80 as a range) on the
  intraslab events' own geometry vs their observed κ; tier-2 layer events' κ vs LFE κ.
- Receiver functions: Audet/Bostock κ 2.35–2.8 for a thin low-velocity layer at the slab top. Reported as the
  thickness h of a κ_RF layer inside the 25 km–source layer that the measured κ_layer allows (with the rest at the
  CVM value), noting that LFE rays start at the interface and may barely cross a layer below the sources.
- Reading rule: κ_layer counts as measured where the crustal κ passes the gates and the LFE κ is one of the
  passing values; otherwise it is reported as a range with the failed gate named.

**§6 Amendment 1 (2026-10-13 ~10:55 am PT, coordinator request, after tomo/eq/RESULTS.md v1).**
(a) Zone-label bug in `s11_eq_diff.py` fixed (crustal zone labels are "+30", not "30"); per-zone ray-weighted
layer κ and the layer figure rerun, no other change.
(b) **Picker test (does the LFE picker bias κ upward?)** `s11_eq_pickertest.py`. Events: the selected WA NW / WA SW
earthquakes with picks, deep sources (25–45 km "layer" + intraslab; crustal as a column), whose origin day is in
the existing 40-Hz 1–8 Hz cache (`lfe_array_discovery/cache40hz`, the data the LFE stacks came from) at ≥ 4 of the
region's LFE stations (NW: picks_BV3 stations; SW: picks_CV3). Per event and station: 3C window [origin − 4 s,
origin + 40 s] (44 s, as the LFE stacks), cut with the a4f DayCache (invalid samples → row dropped), picked with the
frozen picker (`a4f_hypo.make_picks`) and the V3 P picker (`s11_tomo_repick.pick_p_v3`) exactly as on LFE stacks.
Variants: raw; noise-matched (real noise from the same station-day, 10 min before the origin, added until the
envelope SNR equals that station's median SNR in the LFE picks — the closest analogue of an LFE stack).
Rows: rP, rS ≥ 1.2 (as LFE), epicentral ≤ 100 km. κ: the same pooled clipped tP-on-tS fit for (1) picker picks,
(2) the catalogue picks on the same event-station rows, (3) catalogue picks on all of these events' rows.
Columns: per-row ΔtP, ΔtS (picker − catalogue), fraction |ΔtP| > 0.5 s.
**Decision rule (written now):** bias = κ(picker) − κ(catalogue, same rows), noise-matched V3 as headline (the WA
LFE headline used V3 P), frozen-P as column. bias ≥ +0.10 → the WA LFE excess (~0.15–0.2 over the earthquake
deep-path κ) is picking, not rock (fully if ≥ +0.15); |bias| ≤ 0.05 → the picker is not biased on earthquake
records, and the excess stays unexplained (emergent LFE waveforms could still bias it; next test: synthetic LFE-like
stacks); in between → partly. σ of the bias by bootstrap over events (paired).

**§6 Amendment 2 (2026-10-13 ~11:16 am PT, coordinator request; rules written before the runs).**
(a) Volume check (done first, `eq/tables/volume_split_kappa.csv`): the "25–45 km" events split into overriding-plate
events (median 17–22 km ABOVE Slab2, depth ~26 km) and near-interface events (|dz| ≤ 5 km, ~37–41 km, the LFE source
level); intraslab split at 15 km below Slab2. Same fit and gates.
(b) Synthetic LFE-like test `s11_eq_synth.py`: rows = every (family, station) of the WA NW (picks_BV3) and WA SW
(picks_CV3) families at 25–50 km with a known SNR; tP from the CVM v1.7 ray model (family depth, epicentral
distance), tS = κ·tP, κ ∈ {1.75, 1.90, 2.10}. Emergent wavelet: carrier 3 Hz, random phase, envelope t²·e^(−t/0.25 s);
P on Z (radial 0.3), S on the horizontals at 4× P (random azimuth) with S leakage on Z 0.3× S (column: 0); real noise
windows from the station's own cached days (1–8 Hz, 40 Hz) added until the envelope SNR equals that row's LFE-stack
SNR; noise-free column. 3 noise draws. Frozen and V3 pickers, rows rP, rS ≥ 1.2, pooled clipped tP-on-tS fit.
Report recovered − true. **Rule:** at LFE noise, if κ_true 1.75 is recovered ≥ 1.85 (bias ≥ +0.10) the WA LFE κ of
1.9–1.95 is consistent with a normal ~1.75–1.8 source region (picker bias); if |bias| ≤ 0.05 at all three κ, the
picker cannot manufacture 1.9 and the WA LFE excess is in the waveforms (real or source-side); in between → partly.
(c) SVI picker test `s11_eq_pickertest_svi.py`: the Amendment 1b test on SVI deep (layer + intraslab) earthquakes
2003–2013 from the 20-Hz backward cache (upsampled 2/1, as for the era LFE picks). Expectation if the method is
sound: picker − catalogue ≈ the WA raw-record value (+0.06) or ≈ 0; a much larger SVI value would say the eq test
itself is unreliable.

**§6 Amendment 3 (2026-10-13 ~11:31 am PT, coordinator request; written before the run).** `s11_eq_swidth.py`, existing
LFE stacks only (WA NW picks_BV3 / WA SW picks_CV3 stacks; SVI era stacks upsampled 2/1).
(a) Pulse width: per row with rS ≥ 1.2, horizontal Hilbert envelope; S peak = envelope maximum within ±0.5 s of the
S pick; FWHM by walking out to half maximum (cap ±5 s); "onset" = left half-max crossing. P the same on the vertical
envelope around tP (rows rP ≥ 1.2). Width vs epicentral distance: OLS slope (ms/km) with bootstrap over families,
Spearman ρ; per region. A broadening is "a finding" only if the slope is > 0 at 3σ in WA and larger than in SVI
at 2σ; it is then reported as a path property (attenuation/scattering) with the caveat that stacking jitter also
widens far-station pulses.
(b) κ refits, same clipped tP-on-tS fit and gates: (i) V3 picks, stations ≤ 50 km epicentral (column ≤ 40, ≤ 60);
(ii) half-max ONSET picks (P and S) instead of peak picks, all distances and ≤ 50 km. Rule: if the onset or ≤ 50 km
κ in WA falls to within 0.05 of the near-interface earthquake κ (NW 1.775, SW 1.735), the +0.17 excess is a
peak-picking/broadening effect; if it stays ≥ 1.88 it is not.

**§6 Amendment 4 (2026-10-14 ~12:33 pm PT, coordinator request; written before any of the new measurements below).**
Status first: the Amendment 3 run on the LFE stacks (`eq/tables/swidth_rows.csv`, `swidth_slopes.csv`, `swidth_kappa.csv`)
finished 10-13 ~11:42 am PT, before the pod restart (~12:20 pm PT), and was not written up. Those numbers have been
seen; they are reported under the Amendment 3 rules unchanged. Everything below is new and fixed now.
(a) Earthquake-record widths (`s11_eq_swidth2.py`): rows = the `raw` rows of `eq/tables/pickertest_rows.csv` (WA NW, WA SW
deep earthquakes at the LFE stations, 40-Hz cache) and `pickertest_svi_rows.csv` (SVI, 20-Hz cache upsampled 2/1),
re-cut with the same windows. Same width function as Amendment 3 (horizontal Hilbert envelope, peak within ±0.5 s of the
picker's S pick, FWHM to half maximum, cap ±5 s; P the same on the vertical; rows rS ≥ 1.2 / rP ≥ 1.2). Band 1–8 Hz
(the κ band) as headline; 2–8 Hz as a column, computed for the LFE stacks too. Quantities: wS, wP, wS − wP,
S rise time (envelope peak − half-max onset), pick − envelope-peak delay, and where catalogue picks exist the picker −
catalogue ΔtS, ΔtP (a direct onset-to-pick delay). Each vs epicentral distance: OLS slope (ms/km), bootstrap over
events/families (1,000), Spearman ρ. Scrambled control: distances shuffled among the rows of each event/family
(1,000) → slope distribution; p = fraction of shuffled |slope| ≥ real |slope|.
(b) Synthetics with the measured broadening (`s11_eq_synth_width.py`): rows, ray geometry, real noise, SNR matching,
pickers and fit exactly as Amendment 2b, κ_true 1.75 (headline) and 1.90, S leak 0.3, 3 noise draws + noise-free.
Envelope t²·e^(−t/τ) with τ = FWHM/c (c = FWHM/τ of that envelope, computed numerically), separately for P and S,
wavelet length 6 s. Variants: W0 constant τ = 0.25 s (repeats Amendment 2b); **W1 headline** = region trend wS(d), wP(d)
from the OLS fits on the LFE-stack rows (1–8 Hz); W2 = each row's own measured LFE wS, wP (trend where missing);
W3 = W1 with 3× the distance slopes (sensitivity); W2s = scrambled control, W2 widths permuted across rows of the
region (breaks the width–distance link). bias = recovered − true. Closure: wS of the noisy synthetics vs distance,
measured with the same width function, against the observed slope. Rule: W1 bias ≥ +0.10 → broadening explains most
of the WA excess; |W1 bias| ≤ 0.05 → it does not; between → partly.
(c) Distance-limited κ: the Amendment 3(b) values (epicentral ≤ 40 and ≤ 50 km; family minimum 3 and 4 rows) with the
Amendment 3 rule (within 0.05 of the near-interface earthquake κ, NW 1.775 / SW 1.735 → measurement effect; ≥ 1.88 →
not). Columns: hypocentral ≤ 50 / ≤ 60 km. Scrambled control: for each family, the same number of rows as it has
≤ 50 km (resp. ≤ 40 km) drawn at random from all its rows (distances shuffled within family), clipped tP-on-tS κ,
300 draws; the real distance-limited κ is read against this distribution (percentile), i.e. whether near stations give
a different κ from random station subsets of the same size.
Verdict rule overall: the WA excess is called a broadening/peak-picking measurement effect only if (b) gives W1 bias
≥ +0.10 OR (c) drops to within 0.05 of the earthquake value; LFE-vs-earthquake and WA-vs-SVI slope contrasts are
reported with 2σ as the threshold for "different".

**§6 Amendment 5 (2026-10-14 ~2:50 pm PT, coordinator request; written before any of the numbers below were computed).**
Two tests of where the WA LFE κ excess (NW 1.95, SW 1.90 vs near-interface earthquakes 1.775 / 1.735) comes from, on existing
data only. Script `s11_eq_netsplit.py` (≤ 3 workers, OMP_NUM_THREADS=1). Rows = `eq/tables/swidth2_lfe_rows.csv` (the V3 picks
of the Amendment 3/4 κ fits, rP, rS ≥ 1.2, tS > tP); fit = the same clipped pooled tP-on-tS (E1) with per-family intercepts and
the Amendment 3 gates (σ ≤ 0.05, |z| vs the pairing-scrambled control > 4, LOSO ≤ 0.05, |station-term shift| ≤ 0.02).
Sanity: the all-station fit must reproduce NW 1.950 / SW 1.904.
(a) **Network split.** Station classes from `a4f_common.NET`: PERMANENT = UW, US (PNSN/ANSS) and PB (PBO boreholes);
TEMPORARY = XU (CAFE), TA (USArray; C04A TA→XU), YW (FACES, FACB). (No XG Array-of-Arrays station is in the 25-station NW/SW
pools; noted, nothing to fit.) Columns: PB only, UW/US only, XU only. κ per class and region, families with ≥ 4 class rows
(≥ 3 column); also on the common families (≥ 3 rows in both classes) and at epicentral ≤ 50 km.
Scrambled control: class labels permuted across the region's stations (class sizes kept), both groups refit, 300 draws →
Δκ = κ_perm − κ_temp distribution; column: labels permuted within station-distance terciles (median epicentral distance).
Earthquake column: the picker κ on the WA raw earthquake records (`pickertest_rows.csv`, origin-relative windows, event
intercepts) per class. **Rule:** one class "carries the excess" if its κ is within 0.05 of the near-interface earthquake
value while the other class is ≥ 1.88, AND the real Δκ lies outside the scrambled control's 2.5–97.5 % range; both
classes ≥ 1.88 → the excess is not a network (timing/orientation/response) effect; otherwise "partly", with the numbers.
(b) **Per-station residuals vs the SVI-calibrated model.** Model tS − tP = (κ_cal − 1)·tP_model, κ_cal = 1.773 (the SVI LFE
κ that matches the SVI earthquakes), tP_model = CVM v1.7 direct up-going P ray time (`s11_eq_diff.ray_partials`, profiles
royer_NW_WA / royer_SW_WA) from the catalogue location (depth clipped to 25–50 km) to the station. Residual
r_abs = (tS − tP)_obs − model; r_rel = r_abs − family median (removes family location/origin error). Per station: median
r_abs, r_rel, n, median distance, class; map (stations coloured by median r_rel, shape by class) and r_rel vs distance by
class. Column: κ_cal = 1.775 (WA NW near-interface earthquakes). Information only (no pass rule); read with (a).
(c) **B013 check.** The ballistic round-3 step (B013 horizontal S −3.9 ms after 2015, S−P −4.5…−7.7 ms) and the S14 B013
2012-08…2013-05 CC period: (i) the dates of the data in the WA stacks (family era_start/era_end); (ii) κ with B013 tS
shifted by −5 ms and −50 ms (10×), with B013 dropped, and with the B013/B943/W030 site dropped (NW all and ≤ 50 km);
(iii) B013 − B943 and B013 − W030 ΔtS, ΔtP on the same families (median, robust SD). Rule: B013 "matters" only if any
(ii) change exceeds 0.02.
(d) **S-pick scatter correction for the near-station fits.** Sets: all stations (≥ 4 rows), epicentral ≤ 50 km and ≤ 40 km
(≥ 3 rows; ≥ 4 column), regions WA NW, WA SW, SVI (reference). On the rows kept by the E1 clipping:
- **Headline: instrumental-variable estimator (IV).** z = tP_model (as (b)); κ_IV = Σ z_d·tS_d / Σ z_d·tP_d with family-
  demeaned quantities. Pick errors are independent of the station geometry, so IV is not biased by P or S pick scatter and
  needs no σ_S. σ = max(family, station bootstrap, 1,000); LOSO.
- Columns: E3 = κ_E1·R_S, R_S = 1 − σ_S²/var_within(tS), σ_S from co-located pairs (Amendment 7 code; fallback 0.092 s)
  and 2·σ_S; E2 Deming (δ = σ_P²/σ_S²); SIMEX (Gaussian noise added to tS at λ = 0.5, 1, 1.5, 2 × σ_S², 200 draws per λ,
  pooled slope, quadratic extrapolation to λ = −1).
- Synthetic closure on `synth_rows.csv` (real WA geometry and LFE noise, κ_true 1.75 and 1.90, S leak 0.3, draws 1–3, V3 P):
  E1 and IV at all / ≤ 50 / ≤ 40 km; the synthetic IV instrument is computed from the family location perturbed by
  N(0, 5 km) per horizontal axis and N(0, 5 km) in depth (one draw per family), so it carries location error like the data.
- Scrambled control: per family, the same number of rows as it has ≤ 40 (≤ 50) km drawn at random from all its rows, IV κ,
  200 draws; real near-station IV κ read as a percentile.
- **Rule:** the near-station rise "persists" if IV(≤ 40 km) − IV(all) > 2σ_diff AND > 0.05 (σ_diff from a joint family
  bootstrap of both fits); "removed" if ≤ 2σ_diff or ≤ 0.05. The synthetic closure "explains" the E1 rise if the synthetic
  E1(≤ 40) − E1(all) at κ_true 1.75 is ≥ half the real E1 rise. The corrected all-station IV κ is compared with the
  near-interface earthquake κ (within 0.05 → the excess is removed by the correction).
Verdict wording is fixed by these rules; outputs `eq/tables/netsplit_*.csv`, `eq/figures/fig_netsplit_*.png`.
