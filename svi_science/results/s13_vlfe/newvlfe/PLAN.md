# newvlfe PLAN — can a template-free detector find NEW VLFEs? (pre-registered)

Written 2026-10-16 ~10:05 am PT, before any battery, control or separation number of this study.
Question (user 2026-10-16): "we are only detecting repeats and no new VLFEs." The catalogue's matched filter (bank 2011 ∪
2018 relay) only finds events that resemble templates near the C0/C1 patches. Test whether the template-free scan, now
vetted by the denser v2 network (12–14 broadbands ≤ 160 km + 5 far LHZ references), can find VLFEs elsewhere.

Disclosure (looked at before writing this): the raw counts of positive scan-v2 candidates per z threshold for the REAL
runs (2013: 31,006 at z ≥ 5 … 7,162 at z ≥ 20; 2018: 22,969 … 6,093) and how many real candidates pass the cheap gate
(location rule c + earthquake-time masks; 2013: 13,172; 2018: 7,310). No control run, no battery value, no ETS/between
split, no known/new split was looked at.

Words: "scrambled control" = the same scan run on data scrambled so no real event can line up: (A) station-day shift,
(B) polarity scramble of the grid templates. PT = Pacific time. "Known" / "new" = location class below.

## 1. Candidate pool (no template involvement)
- Years 2013 and 2018, scan v2 (the catalogue's v2 recipe: `templates_grid_v2.npz`, `bp_<Y>v2`, gain fallback g = 1 on Z,
  same grid, noise rule, MAD normalisation; `../scripts/scan_detect.py`).
- Real: the existing `data/scan_<Y>/runs/real_v2` (read only).
- Controls, run now with a copy of scan_detect.py in `newvlfe/scripts/nv_detect.py` writing to `newvlfe/data/runs/<Y>/`:
  - (A) station-day shift, shifts s = 3, 7, 11: station k of the v2 list uses day (d + s·k) mod NDAY for **all** scan
    stations (the original shifts only the 7 pilot stations, which would leave the 5–7 regional v2 stations aligned);
    the references are shifted by s·(k_ctrl + 30) so that no reference shares a station's offset.
  - (B) polarity scramble: 4 draws in one pass (the scan's own rule: |ρ| ≤ 0.3 against the true template).
  - So 3 + 4 = 7 control catalogues per year, 14 in total.
- Rows: positive candidates (sign > 0, the thrust sense), z ≥ 5, all z thresholds 5, 7, 10, 15, 20 reported.
- Location = the candidate's best grid node. Time = its tsec on the candidate's timeline (PGC frame for controls).

## 2. Battery (every column kept; tiers pre-set here)
Cheap gate (the catalogue's compute shortcut, exact for every tier below): c (z ≥ 5, node not on the grid edge, inside the
footprint polygon or ≤ 20 km from an LFE family) AND the e′ time part (no bw5/reg/reg3 second in [t0 − 60, t0 + 120] and no
head2 window of M ≥ 5.5). Rows failing the gate fail every tier; they stay in the table with cheap columns.
Light stage (exact for every tier: f is in every tier): z_loso and z at the best node; rows failing f
(z_loso ≥ 3 AND z_loso ≥ 0.5·z) fail every tier.
Full columns (`scan_columns2.cols` + `hg_r3.plane_wave2`, the catalogue's machinery, read-only imports) for the rest:
ratio_far / n_far over all live far LHZ references (COR, NEW, HAWA, OMAK, I02D, MSO, WVOR), per-reference ex_/ratio_,
GNW ratio, vr_thrust, vr_free, share_excess, max_share, isolation, HF ratios (hf_loud), pw2_n / pw2_vr / pw2_nfar250, m0.
- **B0 = battery′** (catalogue PLAN_CAT §3): a (n_far ≥ 2, ratio_far ≤ 0.551) ∧ b (not travelling: pw2_n ≥ 8 ∧ pw2_vr ≥ 0.5
  ∧ ≥ 2 kept stations ≥ 250 km) ∧ c ∧ e′ ∧ f.
- **B1 = B0 ∧ n_far ≥ 3 ∧ share_excess ≤ 0.25 (no single-station dominance) ∧ vr_thrust ≥ 0.12 ∧ z_loso ≥ 5.**
- **B2 = B1 ∧ not travelling over the dense ring** (pw2_n ≥ 8 ∧ pw2_vr ≥ 0.5, without the ≥ 250 km clause).
- Tiers = {B0, B1, B2} × z ≥ {5, 7, 10, 15, 20} (15 tiers).
- HF-loud (2–8 Hz envelope ≥ 10× the day median at the two nearest pilot stations) is a class column, never a cut.
- **Validation of the added tests (rank, don't discard):** each test added in B1/B2 must keep ≥ 80 % of the catalogue's
  production tier events of 2013 + 2018 (real, all periods; their stored battery columns in
  `catalogue/data/prod2013/rows_scored_v2.parquet`, `prodr2018/rows_scored_v2.parquet`) that pass B0. A test that fails
  this is removed from B1/B2 and kept as a column (reported).
- Events: candidates within ±150 s are already declustered by the scan; no further merge.

## 3. Separation test
- Periods (as production): ETS = any SVI episode of the year (`episodes_b.csv`, start − 2 d to end + 3 d); "other" = NW-only
  episodes (reported, not tested); between = the rest.
- Location class: distance from the best node to the nearest bank template location = H&G T1–T5 (Table 1), the 14
  vetted + near-miss events of clusters C0–C5 (`catalogue/tables/e3_events.csv`), and the 17 2018 relay events RC0–RC7
  (`prod2018_relay_events.csv`). (a) **known** ≤ 15 km; (b) **new** > 15 km. Column: the same at 30 km.
- Live time (per run, per period): seconds with ≥ 3 scan stations available (navail) and outside the e′ time masks.
- Rates per live year: real = count / live yr; control = total count over the 7 control catalogues / their summed live yr
  (pooled over 2013 + 2018 for the headline; per year as columns). Control floor when 0 events: 0.5 / summed control live yr.
- **Pass for (b) new, separately ETS and between:** some pre-set tier with (pooled 2013 + 2018) real count ≥ 3 AND real
  rate ≥ 3 × control rate AND control rate ≤ 1/yr. Columns: Poisson p of the real count against the control rate × real
  live time; whether the pass also holds at the next z threshold.
- (a) known is reported the same way (positive check: the scan-only pool should separate there in the ETS if the dense
  vetting works).

## 4. Every surviving "new" candidate (any tier that passes, else the B2 real "new" rows at the lowest passing-free tier
  listed as candidates) gets:
- raw waveform figure: 0.02–0.05 Hz vertical displacement, near (≤ 40 km), mid (40–150 km), far LHZ, −120 to +300 s;
- location (node lat/lon/depth), nearest LFE family and distance;
- LFE/tremor coincidence: any SVI P ≥ 0.75 or WA-NW selected LFE within ±10 min; patch-specific: SVI P ≥ 0.75 LFEs from
  families ≤ 20 km of the node within ±10 min vs 1,000 random same-UTC-day times (p); PNSN tremor within 100 km ±10 min;
- storm percentile (PGC Z 5–10 s, day median of 10-min RMS, percentile in the year);
- whether any production template detected it (bank T1–T5, C0–C5, + RC0–RC7 in 2018; MF ≥ 5.5 MAD within ±30 s; tier flag).

## 5. If (b) passes
Stack the surviving new events per location cluster (complete linkage 15 km) into candidate templates (catalogue recipe V4,
same stations), leave-one-out re-find and 8 + 8 scrambled-control rates under the catalogue acceptance rule (relay rule:
LOO ≥ 85 % and scrambled ≤ 1/yr). Report only; nothing is added to the production bank.
If (b) fails: report the separation numbers and the closest tier; no template stage.

Rules: ≤ 2 CPU workers (GPU for the scan runs), memory < 40 GB, data on disk only, nothing written under `catalogue/`
or into shared scan directories (control runs live in `newvlfe/data/runs/`), PT times, raw figures first, points coloured,
no error bars. Log: `NOTES.md`; daily note `notes/2026-10-16_Notes.md`.

## Amendment NV3 (2026-10-16 ~4:55 pm PT, before any 2011/2016/2019 number) — extend to 2011, 2016, 2019 (user approval relayed by the coordinator)
Unchanged: candidate pool (positive scan-v2 candidates, z ≥ 5), scrambled controls (4 polarity draws + station-day shifts
s = 3, 7, 11 of all stations, references at k + 30), cheap gate, light f stage, full columns, known/new split (15 km
from the 2018-bank template locations, 30 km column), periods (SVI episodes −2/+3 d = ETS; NW-only = other; rest =
between), live-time rule, LFE/tremor columns.
Changes (lessons since, `catalogue/MERLIN_CAT_REVIEW5.md`):
1. **W2b in every year.** The scan runs zero the station-component-months in `catalogue/tables/qc_exclusions.csv`
   (`nv_detect.py`, NV_W2B=1). Real = variant 0 of the W2b polarity pass (polBW), plus shift controls shA3W/shA7W/shA11W.
   2013 is re-run the same way (its W2b list has PGC 3 months and LZB 1), so the pooled result is all W2b. 2018 is
   already W2b (identical to the non-W2b run, since the scan noise rule had already excluded LRIV); the 2018 shift
   controls are reused. 2011 uses `data/scan/bp_2011v2` and the reference arrays `bp2_ref_v2` (the production rule for
   2011–12).
2. **Ladder instead of tiers.** Level = number passed of 8 checks: the catalogue battery′ a, b, c, e′, f plus the 3 added
   checks that passed their validation (n_far ≥ 3; share_excess ≤ 0.25; no dense-ring plane wave). n_checks_missing is
   shown: rows stopped at the gate or at f have unevaluated checks, so level 8 is exact and lower levels are lower bounds
   for them. vr_thrust (an SNR proxy, review 5 flaw 3), Q2 and Tier A′ are columns only, never gates.
3. **New columns** per fully evaluated row, with each station's own source day so shift controls are scrambled
   consistently. Emergence follows the catalogue Y9 definition: window peak in [t0, t0 + 180 s] over the median Hilbert
   envelope of ±30 min, 20–50 s Z. Columns:
   - emerg_med_near over non-far stations ≤ 150 km;
   - n_near_emerg_ge2;
   - emerg_max_over_second and emerg_lead_sta;
   - emerg_far_med;
   - t3_near_mid, t3_far_over_nearmax, t3_spearman.
   **Earthquake-wave flag** = emerg_far_med ≥ emerg_med_near OR an M ≥ 4.5 event (scan eq table) with origin in
   [t0 − 90 min, t0 + 1 min].
4. **Separation, same rule as before, on the ladder.** For each level ≥ 6 × z ≥ {5, 7, 10, 15, 20} × known/new ×
   ETS/between, pooled over 2011, 2013, 2016, 2018 and 2019 (and per year): real count, live years, control rate, ratio,
   Poisson p. A "separating" level needs real ≥ 3, ≥ 3× the control rate and controls ≤ 1/yr. Reported with and without
   earthquake-wave-flagged rows (the flag is a column, both counts shown).
5. **Blind eye labels.** Every real "new" row at level 8 (all periods, all five years) and the same number of control rows
   drawn at random from the control rows at level 8 (all types; matched by year where possible) get the same raw record
   section. Figures are named by a random ID, real and control mixed, and the key goes in `data/nv3_label_key.csv`
   (not opened until labelling is done). Each figure gets one label: wavefield, eq-wave, spike or candidate.
   Comparison = label fractions, real vs control.
6. **Any new source that holds up** (candidate labels concentrated in place and time, beyond the controls) gets the
   stacking and relay acceptance (LOO ≥ 85 %, scrambled ≤ 1/yr), as for the 2018 western group.
Resources: ≤ 2 CPU workers, one GPU process at a time (sequential), < 40 GB, no downloads. Queue `scripts/nv3_queue.sh`,
marker `NV3_DONE`; failures as FAIL lines in `logs/queue.log`.
- NV3 §5 amendment (2026-10-16 ~5:15 pm PT, after a code test on the first finished year and before any control,
  separation or label number). Disclosed: the 2016 real run alone has 185 "new" level-8 rows (all periods). Labelling every
  row by eye in all five years is not feasible, so the label set is a year-stratified random sample of 120 real "new"
  level-8 rows (all of them if there are fewer) plus the same number of control rows matched by year. Label fractions come
  with their sampling scatter. Also noted from the same test: an M ≥ 4.5 event within 90 min is common (the global rate
  is about 19 per day), so that half of the eq-wave flag has a high base rate. The flag stays a column, shown for real and
  controls alike.
