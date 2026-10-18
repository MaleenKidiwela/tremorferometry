# S13 PLAN — long-period (20–50 s) stacks on our LFE times

Written 2026-10-13 (Tue), before any stack was computed. Rules below are fixed before the pilot numbers exist.
Any later change is an Amendment with its PT time, written before the numbers it affects.

## 0. Question
Do our LFE detection times carry a coherent 20–50 s (0.02–0.05 Hz) signal, i.e. a very-low-frequency (VLF) moment that
is time-locked to the LFEs? If yes, the long-period moment per family per episode is a per-patch slip measure for ~20 ETS
per region (2003–2026), independent of the LFE count.

Background: Hutchison & Ghosh 2018 (VLFEs Mw 3.1–4.3, repeat, track GPS, continue between ETS; their 2011 templates 2–4
sit at 48.3–48.4 N, 123.45–123.6 W, 30–54 km, right under our SVI families); Ide 2016 stacked tremor at 0.02–0.05 Hz.

## 1. Expected size (from `scripts/s13_green.py`, whole-space + free-surface ×2, thrust 315/12/90 at 35 km)
- Peak 20–50 s displacement per 1e12 N m (Mw 1.9): 0.4–1.2 nm at 0–40 km epicentral, 0.2–0.5 nm at 60 km,
  0.1–0.4 nm at 100 km. SVI LFE magnitudes (Bostock mag_med) are 1.3–2.1, median 1.7 (M0 ≈ 4e11 N m).
- If each LFE carries ONLY its own high-frequency moment: ≈ 0.2–0.4 nm per LFE at the near stations.
- Band noise on a good vault, vertical, 20–50 s: ≈ 3–10 nm RMS (NLNM ≈ 2.7 nm); horizontals higher (tilt).
- Pooled SVI episode (1,400–2,900 LFEs at P ≥ 0.75): SNR ≈ 0.3/5 × √2000 ≈ 2.7 per station-component, ~6–8 σ over
  6 stations × 3 components with a matched template. Per family per episode (~50–100 LFEs): SNR ≈ 0.5 → per-family
  numbers need a VLF moment several × the LFE moment, or pooling (zones / several episodes).
- If LFEs ride on VLF slow events (slow-earthquake scaling; Ide 2016, Ghosh 2015) the VLF moment per LFE can be ≫ the
  HF moment and the stack is easy. The pilot measures which case holds.

## 2. Data
- Catalogue SVI: backward `data/catalogue/svi_catalogue_v1_4_backward.parquet` (2003–2013; rows `primary &
  keep_xfam8==1`, `series=='data'`) + forward `svi_science/results/s0/analysis_table.parquet` (2014+, run prod,
  keep_xfam8). Same loader as `s10_multi_inventory.load_rows`. Headline cut **P ≥ 0.75**; P ≥ 0.5 a column.
  Counted-once rule (`counted_once`) for the backward rows.
- Family locations: `data/catalogue/families_v1_3.csv`, region S_Vancouver_Is, `family_id == fam`.
- Origin time per LFE: t (frame time) + offset_S(sta) − r(sta)/3.95 km/s, median over stations with an offset in
  `lfe_extension/results/cut/svi_ext/offsets_used.csv`. (1–2 s error is irrelevant at 20–50 s.)
- Waveforms: raw-count continuous day files already on disk (inventory, step 1, `tables/inventory_*.csv`): CN
  broadband BH (40 Hz) / HH (100 Hz) and other broadband. Boreholes (PB EH*, 2 Hz geophones) and short-period EH* are
  excluded (no usable 20–50 s response).
- Response: StationXML (response level) from EarthScope FDSN, cached in `data/resp/`.
- Processing per station-day (`s13_prep.py`): merge (gaps ≤ 5 s interpolated, longer gaps → masked windows), demean,
  detrend, 5 % taper, decimate to 1 sps (FIR anti-alias), remove response to displacement with pre-filter
  (0.004, 0.008, 0.25, 0.45) Hz, rotate to N/E/Z with the StationXML azimuth/dip, store float32 1-sps day arrays
  in `data/disp1sps/NET.STA/YYYY/DOY.npz` (≈ 1 MB/station-day). Bandpass 0.02–0.05 Hz Butterworth 4-pole zero-phase
  on the whole day, then cut windows (no window-edge effects).

## 3. Pilot (pre-registered)
- Episodes: **2011-08** (2011-08-12 → 09-03, Ghosh 2015 / Hutchison & Ghosh 2018 period) and **2012-08**
  (2012-08-31 → 09-28). Window = episode ±2 d. Column: 2010-08b1.
- Families: the 5 families with the most P ≥ 0.75 counted-once rows summed over the two episodes (expected among
  10, 56, 43, 0, 3). Also **pooled** = all SVI families.
- Stations: every broadband station with episode days on disk within 150 km of the SVI family centroid.
- Stack: window −150…+250 s about origin, 1 sps, linear mean (headline) and P-weighted mean (column), per
  (unit, episode, station, component); unit = family or pooled.
- Window QC (identical for data and control): drop a window if its band RMS > 3 × the median band RMS of all
  windows (data + control) of that station-component-episode, or if > 5 % of its samples are masked.
- Scrambled controls (same N, same QC):
  (a) **random times within the same UTC days** (each LFE moved to a uniform random time on its own day), 500 draws;
  (b) **shifted ±1 day** (t ± 86,400 s), two stacks.
- Template fit: synthetic 20–50 s displacement for a thrust (strike 315, dip 12, rake 90) at the unit's
  count-weighted location (pooled: count-weighted centroid), all station-components jointly, one common time shift
  scanned ±8 s; M̂0 = Σ s·d / Σ s·s (N m per LFE, = mean time-locked VLF moment), VR = 1 − |d − M̂0 s|² / |d|².
  Each station-component is weighted by 1/σ² where σ = SD of its control stacks (a).
  The same fit is run on every control draw → z = (M̂0 − mean_ctrl) / SD_ctrl, p = rank.
- Alternatives (columns): strike-slip 315/90/0, steep dip-slip 315/80/90, free 5-element moment tensor
  (deviatoric LSQ, same weights).
- Neighbour count: mean number of other catalogue LFEs (any family, P ≥ 0.5) within ±10 s of each stacked LFE →
  M̂0 per LFE corrected = M̂0 / (1 + n_nb) (column).

**Pass rules (pilot):**
- **P1 detection (pooled):** M̂0 above the 99th percentile of control (a) in ≥ 1 of the 2 episodes, and above both
  ±1-day stacks.
- **P2 thrust sign and pattern:** in that episode M̂0 > 0, VR ≥ 0.2, and the thrust VR is higher than both
  alternative fixed mechanisms; among station-components with template amplitude above the median, ≥ 70 % have the
  stack sign of the template at the best lag.
- **P3 replication:** odd-day and even-day halves each give M̂0 > 0 with z ≥ 2.
- **P4 per family:** ≥ 2 of the 5 families give z ≥ 3 in ≥ 1 episode (else per-patch numbers need pooling into
  zones or episodes; this does not stop the scale-up).
- Pilot "works" = P1 + P2 + P3. A non-detection is reported with the minimum detectable M̂0 per LFE
  (3 × SD_ctrl of M̂0), its Mw, and the obstacle (noise level per station-component).

Cross-check (column, no rule): Hutchison & Ghosh 2018 Table 1 template times 2011-08-21 04:09:06, 08:43:29,
2011-08-24 10:47:38 UTC (templates 2–4): number of our LFEs within ±60 s vs the day's mean rate.

## 4. Scale-up (only if the pilot works; rules finalised as an Amendment before the scale-up runs)
1. All SVI families × all episodes 2003–2026 (episodes_b SVI majors + Bostock-era p3 episodes), then WA-NW
   (`lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet`, selected_v2, P_own_pnet ≥ 0.75).
2. VLF moment per family per episode (or per zone if P4 fails), with the control z.
3. Moment partition: VLF moment / LFE count and / summed LFE amplitude per episode.
4. Versus GNSS slip per episode (S1 `s1/episode_geodetic.csv` D, tierB_M0) and versus the WA "larger episodes read
   higher dv/v" hint (s11_depth/round7/RESULTS.md).
5. Between-ETS bursts (q5_bursts/bursts.csv): do they carry VLF moment (small slow slips)?
Data beyond disk: 1-sps LH or BH decimated, ETS ±30 d, only if the inventory shows gaps (rule: total ≤ 20 GB → fetch).

## 5. Figures
Line graphs: pooled and per-family stacks per station-component with the template overlaid and the control 95 % band
as thin grey lines (no error bars); per-episode M̂0 line graph coloured by z/weight.

## Amendment 1 — scale-up rules (written 2026-10-13 ~10:55 am PT, after the 2011-08 pilot numbers and before the
## 2012-08 pilot verdict; applies only if the pilot verdict says "works"; no scale-up number exists yet)
Scripts: `s13_units.py` (rows), `s13_prep.py --days-file`, `s13_stack2.py` (pilot physics, QC and fit imported unchanged
from `s13_stack.py`), `s13_scale_summary.py`; queue `s13_scale_queue.sh` (markers `STACK_scale_*_DONE`, `SCALE_QUEUE_DONE`).
- **Episodes.** SVI: Bostock-era majors (`p3_bostock/episodes_bostock_final.csv`, major, start < 2010) + `episodes_b.csv`
  SVI majors; WA-NW: `episodes_nw_final.csv` majors. Window = episode ±2 d. Built when ≥ 30 LFEs (P ≥ 0.75).
- **Rows.** SVI as the pilot. WA-NW: `catalogue_v2.parquet` selected_v2 & region NW, P = P_own_pnet, family location
  `wa_prod/setup/families.csv` (lat_use/lon_use/depth_km; missing depth → 35 km, flagged), origin = time + median
  (s_off − r/3.95 km/s) over the family's template stations.
- **Units per episode:** pooled (thrust + alternatives, 500 draws); odd/even-day halves, 3 along-strike zones
  (terciles of family along-strike position, N45W axis), every family with ≥ 20 LFEs (thrust only, 200 draws).
- **Between ETS:** SVI `s7/between_events.csv` (BE, 76, tremor-free) and `s7/between_broad_events.parquet`
  (BE_broad, 688); WA `a_wa_between/between_events.csv` (220) and `between_broad_events.parquet`; one pooled unit each
  (full), same controls (random time on the same day; ±1 day).
- **Stations.** Every broadband station on disk in `tables/inventory_bb_near.csv` within ~100 km (SVI) / ~110 km (WA)
  (lists in the queue); a station enters a unit when its Z covers ≥ 20 % of the needed samples and ≥ 5 windows pass QC.
  No downloads in the first pass (the near SVI/WA long-record stations PGC, VGZ, SNB, NLLB, LRIV, DOSE, GNW, WISH, SP2
  are on disk for nearly all years).
- **Tests (pre-registered):**
  - A1 detection per episode: pooled z ≥ 3. Report the count per region.
  - A2 moment partition: weighted mean VLF moment per LFE; χ²/dof of the episode values about it (σ = control SD);
    p < 0.01 → the partition changes between episodes.
  - A3 size: Spearman ρ, 5,000-permutation p, of the episode total (M̂0 × n) and of M̂0 per LFE against SVI GNSS D
    (S1, via `s11_depth/round6/svi/tables/r6_sizes.csv` G), tremor count T, duration Du, LFE count L; WA against T, Du,
    L (no per-WA-episode GNSS series exists). Pass p < 0.05 (total vs GNSS = SVI headline; total vs T = WA headline).
  - A4 per family: weighted mean M̂0 per LFE per family over episodes; χ²/dof between families (σ from controls).
  - A5 between ETS: pooled z ≥ 3 → the bursts carry time-locked VLF moment; ratio to the ETS weighted mean.
  - A6 (column, no rule): Spearman of the WA and SVI episode totals with the round-7 ΔD per major (`r7_hint_dD.csv`).
- Scrambled control and ±1-day stacks are reported for every unit; the absolute moment scale carries the whole-space
  ×2 uncertainty (relative values are the product).
- (cost) Units with > 3,000 LFEs (SVI 2003-02, 2004-07, 2005-09) use 200 control draws instead of 500 (same rule
  otherwise; p floor 1/201).

## Amendment 2 (2026-10-13 ~11:00 am PT, after the pilot verdict, before any scale-up number)
- §2 said "gaps ≤ 5 s interpolated" but the pilot prep code did not do it: every gap split the day and cost 600 s per
  edge. PGC 2012 has 6–12 s gaps every few minutes, so PGC dropped out of the 2012-08 pilot (1 % of samples left).
  Fix: per channel, consecutive traces with gaps ≤ 15 s are joined by linear interpolation (`S13_MAXGAP`, default 15 s;
  at 20–50 s periods a ≤ 15 s straight-line fill is a small error). Days without gaps are bit-identical to the pilot
  (checked on PGC 2011-08-20, max difference 0).
- The scale-up writes to `data/disp1sps_v2/` (env `S13_DISP`); pilot files in `data/disp1sps/` are kept unchanged.
  The pilot verdict stands as run; `scale_svi_2012-08` (with PGC) is the re-run of that episode.

## Amendment 3 — gap fill (user OK 2026-10-13 for the ~35–45 GB raw download; written ~11:15 am PT, before any fill number)
- Disk: 56 TB free (≥ 200 GB rule met).
- List `tables/fetch_fill.csv`: 1,884 station-days on the SVI/WA episode and between-ETS prep days (±1 d) that are inside
  a BH/HH channel epoch but not on disk. Breakdown: GOBB 706, SP2 211, PFB 209, YOUB 205, WISH 138, LZB 130,
  LRIV 115, VGZ 53, DOSE 49, NLLB 43, PGC 20, GNW 5. Size ≈ 28 GB raw.
- Fetch raw (≤ 3 parallel, resumable) → `data/raw/`; prep → 1-sps displacement in `data/disp1sps_v2/`.
  The new raw is kept for now; it can be removed later with an OK. Nothing existing is touched.
- After the first-pass scale-up finishes, every row file whose days touch a fetched day is re-stacked with tag
  `scalef_*` (CN.GOBB added to the SVI list). The summary prefers `scalef_*`; `scale_*` tables stay as the
  pre-fill column.

## Amendment 4 — channel guard (v1 written 2026-10-14 ~5:52 pm PT; REPLACED by v2 at ~6:25 pm PT, both before any guarded number was used)
Reason: the write-up found flat channels (control SD exactly 0: UW.LRIV N/Z 2019-04/05, UW.OPC WA 2009-04) that take
infinite weight and collapse the moment to ~0, and a near-dead channel (UW.SP2 E 2012–2019) holding 33–70 % of the WA
fit weight.
- **v1 (withdrawn):** drop σ ≤ 1e-3 × the component median, then cap any station-component at 0.3 of the weight. Withdrawn
  after Merlin's review the same evening: the cap is not a defect test — it fired in 847 of 1,623 first-pass units and
  mainly down-weighted the best quiet verticals (PGC Z, NLLB Z, BS11 Z), and 1e-3 × median did not catch SP2 E.
  Stacks run with v1 (`refit_*` partial, SVI `scalef_*`) are moved to `../_quarantine_2026-10-14_envfail/s13_pre_guard/`
  or superseded, never used.
- **v2 (in force, Merlin):** in each unit, drop a station-component if any of
  - n_used < 5 (windows kept by QC);
  - σ ≤ 0 or not finite (σ = SD of the scrambled-control stacks over the fit window, `sigma_ctrl`);
  - σ < 0.05 × the median σ of the same station's other two components;
  - σ < 0.02 × the median σ over all channels in the unit.
  Kept channels keep 1/σ²; **no weight cap**. The same drop set and weights apply to the data stack, every
  scrambled-control draw, the ±1-day stacks, the P-weighted stack and the free moment tensor. Every drop is logged with
  its reason (`stack2.log`, units `guard_dropped`, stacomp `drop_reason`).
- New columns: `max_share` / `max_share_ch` (largest single-channel share of the weighted template energy,
  Σ(tpl/σ)² over the fit window) and the leave-one-station-out moment range (`loso_m0_min`, `loso_m0_max`, with the
  station left out for each).
- Pass checks: SVI 2019-04 and WA 2009-04 / 2019-04 pooled moments no longer ~0; SP2 E dropped 2012–2019; where nothing
  is dropped, the guarded run reproduces the stored first-pass moment exactly (otherwise it is a refit bug: stop).
- Re-runs (SVI only; WA is withdrawn pending the origin-time fix, Merlin F1): fill-affected SVI row files re-stacked
  with v2 under `scalef_*` (v1 outputs moved to quarantine first); non-fill SVI row files with any v2 drop under
  `refit_*`. Summary preference `scalef_*` > `refit_*` > `scale_*`; a unit from a table without `guard_rule = A4v2`
  is used only if v2 would drop nothing in it (then its numbers are identical).

## Amendment 5 — WA origin fix + dedupe + depths, SVI P-bin / thinning, between-ETS jackknife, n_nb (written 2026-10-14 ~6:30 pm PT, before any number below exists; executes MERLIN_SCALEUP_REVIEW.md steps 2–6 and Q1c)
New code and files only; nothing the running fill/refit queues read is touched (`data/rows/`, `s13_units.py`,
`s13_stack2.py` unchanged by this amendment).
- **Code.** `scripts/s13_units_v2.py` → `data/rows_v2/`; `scripts/s13_stack2_v2.py` = `s13_stack2.py` physics, QC, fit
  and the Amendment 4 v2 guard (`s13_common.guard2`), restructured station-outer, plus ±k-day stacks, fixed neighbour
  counts and two extra lag grids as columns. Check before any result: with `--kmax 1` on `svi_2025-09` it reproduces
  the current `s13_stack2.py` on every unit (M0 to 3e-16 relative, z to 2e-15, same lags, same drops; test outputs in
  `../_quarantine_2026-10-14_envfail/s13_v2_test/`). Tags `wa2_*` (WA), `pbin_*` (SVI step 4), `svibe2` (SVI step 5).
- **WA rows (Merlin F1–F3).**
  - Origin `to = time + 0.5 s` (catalogue `time`, PRE of `wa_common.py`); `wa_shift` dropped.
  - Dedupe: all selected_v2 NW rows (all P_own_pnet ≥ 0.75), greedy by P: keep the highest P, drop any other row
    (any family) within ±2 s of a kept row. The post_v2 variant-C flag is joined as a column (count only).
  - Depth: family depth < 20 km or missing → Slab2 depth at the family lat/lon (`data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz`),
    flagged (`tables/wa2_families.csv`).
  - Same episodes, ±2 d windows, ≥ 30 LFE build rule, units (pooled, odd/even, 3 along-strike terciles recomputed on
    the deduped family set, families ≥ 20), same WA station list as the first pass (fill-queue `WA_STA`), same data dir.
  - Between ETS: WA_BE, WA_BE_broad with the same origin rule and the same ±2 s dedupe within each set.
- **WA controls.** Scrambled (random time on the same day; 500 draws pooled/between, 200 other units) and ±k-day
  shifted stacks, k = 1…5, both signs (10 per unit), same QC and weights. z_k = (M0_k − control mean) / control SD.
- **Prior expectation stated now (not a rule change).** The WA template offsets are not travel times (s_off vs
  hypocentral distance: slope 0.005 s/km over 274 families, all offsets −5…+20 s). The 417 "clear" template arrivals
  give S at `time` + 2.7 s at 20–40 km (r/Vs ≈ 8–10 s), i.e. a physical origin ≈ 5–8 s before `time`. Under
  `to = time + 0.5` the best lag is therefore expected near −5…−7 s, inside the ±8 s scan but possibly outside the
  ±4 s lag rule. Columns (thrust, each with its own scrambled-control z): `c6_` grid −14…+2 s, `wide_` grid −40…+40 s.
  The headline stays the ±8 s grid (Merlin).
- **WA pass/fail (Merlin steps 2–3, as written).**
  - W1 lag: mode of the pooled-unit best lag over the built episodes within ±4 s, and ≤ 10 % of pooled units at
    ±8 s. Fail → zone-1 grand stack scanned over ±40 s, residual reported, WA stops there.
  - W2 duplicates: mean `n_nb2` (neighbours within ±2 s, pooled rows, self excluded) < 0.05.
  - W3 control calibration: SD of z_k over all pooled units and all 10 shifts in 0.85–1.15. Outside → every WA σ is
    multiplied by that SD (z divided by it) before A1–A5. SD > 1.3 → day-block-shift control is the next step.
  - W4 zones: zones 2 and 3 each positive in more than half of the ≥ 300-LFE episodes. Fail at the right lag →
    GF-sign check (Okada half-space Z signs at DOSE/GNW per zone).
  - W5 halves: Spearman ρ(odd-day z, even-day z) over episodes > 0.4.
  - Then A1–A6 of Amendment 1 on the `wa2_` tables.
- **SVI step 4 (2003–05).** For 2004-07, 2005-09, 2012-08, 2013-09, the pooled rows split by P into [0.75, 0.90),
  [0.90, 0.97), [0.97, 1]; for 2004-07 also the top 3,000 by P and a random 3,000 (seed 20261017). Thrust only,
  200 draws, ±1-day, same stations and data dir as the first pass. Reading (Merlin's three outcomes, with numbers fixed now):
  - (a) dilution by weak/false detections: 2004-07 top-P bin and top-3,000 ≥ 3.7e12 (within 30 % of 5.3e12) while
    its low-P bin < 2e12, and the 2012-08 bins spread < 1.5× (max/min).
  - (b) genuine era difference (provisional): all 2004-07 bins and top/random-3,000 within 1–3e12 and
    |top − random| < 2σ.
  - (c) accompaniment scales with LFE strength: high/low bin ratio ≥ 1.5 in both a 2003–05 and a 2012–13 episode.
  - None of these → "mixed", reported as such.
  - Also the single-channel Z table 2003–2013 (pooled units, per station Z alone, at the unit's best lag, from the
    stored stacks; σ = channel control σ scaled by the unit's control-SD / white-noise-SD ratio).
- **SVI step 5 (between ETS).** BE_broad (688, P = P_same) re-stacked as: all; minus 2020-09-13/14; minus
  2015-10-08/10; minus the top-5 days by count; tau_next ≥ 30 d; 2014–2019; 2020–2026; plus the P-weighted column.
  Pass: all positive, "minus top-5" z ≥ 2, both year halves positive. Fail: z collapses without 2020-09 → report the
  2020-09-13 burst specifically. The ratio to the ETS value is quoted with its uncertainty (0.79 ± 0.22 in the first pass).
- **Step 6.** `n_nb` from v2 (pooled rows only, self excluded, ±2/10/20/60 s and the band-pass-kernel-weighted sum);
  for SVI episodes computed from the row files (no re-stack).
- **Q1c PGC-out.** From the stored guarded SVI stacks: refit every pooled unit with CN.PGC weights set to 0 (same lag
  grid), σ scaled as for the single-channel table; weighted mean 2008–2026 and all. Merlin's criterion: a change
  > 30 % from the all-station value changes the reading.
- **Addendum (2026-10-14 ~7:35 pm PT, coordinator decision ~7:15 pm PT; before any wa2_ number was looked at):
  UW.SP2 E instrument veto.** guard2 misses the near-dead SP2 E in 75 of 315 first-pass WA units (E/Z 0.1–0.27, above
  the 0.05 clause). It is now an instrument veto in `svi_science/instrument_vetoes.csv` (same mechanism as B013/B009),
  span set from the data, not 2012–2019: **2011-09-05 → 2019-05-22 inclusive.** Daily 20–50 s band RMS shows E/N
  0.7–1.7 until 2011-09-04, 0.38 on 2011-09-05 (dies that day), 0.004–0.04 from 2011-09-06 through 2019-05-21, and
  normal from 2019-05-23 (no data 05-22).
  - `s13_stack2_v2.py` blanks vetoed station-component-days before the band-pass. Data, scrambled controls and
    ±k-day stacks therefore lose the same windows. A component left with < 5 windows is dropped by guard2. Blanked
    component-days are logged in `logs/stack2_v2.log`. Stations without a veto are unchanged (checked bit-identical on
    UW.GNW; SP2 N/Z identical, E blank).
  - `s13_v2_summary.py` refits zero a vetoed component when the unit overlaps the window.
  - Installed with os.replace at 7:29:58 pm PT. wa2_between and wa2_2012-02 had started with the old code. They were
    killed at 7:30 pm PT and are re-run by `scripts/s13_v2_rerun.sh` after the main queue ends.
  - The 8 wa2_ stacks run before the patch (2007-01 … 2011-07) load no day inside the veto (last loaded day
    2011-09-03), so they are identical either way and are kept.
  - Other readers of `instrument_vetoes.csv` (s11_invert, s14_common, s12_era_check) apply a row to the whole
    station, all components. They will now also drop UW.SP2 cells 2011-09-05 → 2019-05-22 if SP2 is in their sets.

## Amendment 6 — WA per-family origin offsets, location-suspect flags, re-dedupe, re-stack; SVI net_mad terciles (written 2026-10-14 ~8:50 pm PT, before any offset, row or stack number below exists; executes MERLIN_WA_ORIGIN_REVIEW.md steps 1–5)
New files only: `scripts/s13_offsets_v3.py` (steps 1–2 → `tables/wa_family_offsets.csv`, `tables/wa3_*`),
`scripts/s13_units_v3.py` (step 3 → `data/rows_v3/`), `scripts/s13_v3_queue.sh` (steps 4–5), `scripts/s13_v3_summary.py`
(step 6). `data/rows/`, `data/rows_v2/`, the wa2_ tables and stacks are not touched. `s13_stack2_v2.py` gets one optional
argument, `--grids` (default = the current '', c6_, wide_ set, so every earlier call is unchanged); v3 runs use
`--grids ,wide_,fix2_` where `fix2_` is the single lag +2 s (no scan).
- **Step 1 — per-family offsets from the catalogue alone.**
  - Input: `catalogue_v2.parquet` rows with `selected_v2` & region NW (the rows the stacks use), catalogue `time`.
  - Pairs: for every ordered family pair (A, B), each A row takes the nearest B row within ±12 s; Δ_AB = median of
    t_B − t_A over those matches, n_AB = number of matches, IQR kept as a column. Pairs with n_AB ≥ 8 enter.
  - Model: t_phys = time_f + d_f. A co-detection gives t_B − t_A = d_A − d_B. Solve the d_f by weighted least
    squares, weight √n_AB per pair (Merlin), separately in each connected component of the pair graph.
  - Absolute anchor (clear picks, `tpl_library.csv` clear == True): per pick δ = s_off − 0.5 − r/Vs, r = hypocentral
    distance from the family location (lat_use, lon_use, catalogue depth_km) to the station (`setup/stations.csv`);
    δ_f = median over the family's clear picks. δ is the physical origin relative to `time + 0.5`.
    Anchor families: located, location_quality A or B, catalogue depth ≥ 20 km (so not Slab2-replaced), not flagged
    location-suspect in step 2, with ≥ 1 clear pick. Component constant c = median over its anchor families of
    (δ_f − d_f). A component with no anchor family is not anchored (its families are treated as off-network).
  - offset_f = d_f + c, so to_f = time + 0.5 + offset_f is the physical origin estimate.
  - Families off the network (or in an unanchored component): offset_f = δ_f if they have a clear pick and are
    not location-suspect (flag `clear_only`); otherwise offset_f = median offset of the anchored network families
    (flag `default`). The between-ETS sets use the same per-fid offsets (fids outside NW: same fallback).
  - **Vs choice.** Headline Vs = 3.95 km/s straight ray, the same constant as the SVI origins (`s13_common.VS_KMS`) and
    as Merlin's "VLF lag − δ = +1…+4 s" anchor, so the W1 lag window means the same thing in both regions. Column:
    USGS CVM v1.7 `royer_NW_WA_S` profile, travel time = r × (vertical-average slowness from the surface to the
    source depth) (exact for a straight ray in a layered model). Relative offsets do not depend on Vs; only c does.
    The CVM constant shift is reported, and the W1 lags under the CVM anchor are the Vs-3.95 lags + that shift.
  - Acceptance (Merlin): pair-network residual rms < 2 s **and** SD of (offset_f − δ_f) over the anchor families
    < 3 s. Failure: that SD > 4 s → stop WA at step 1 (rows not rebuilt). Between 3 and 4 s, or rms ≥ 2 s: proceed,
    WA headline capped at "present, size provisional". Robust column (not a substitute unless stated now): the same
    solve with only pairs of IQR ≤ 2 s; if the primary rms ≥ 2 s and this variant passes both acceptance tests, the
    variant is used.
  - Epoch test (Merlin change-my-mind): fam38–fam75 Δ per epoch 2007–2011 / 2012–2017 / 2018–2024 (≥ 8 matches each).
    Drift = max − min of the epoch medians. Drift > 1.5 s → offsets are solved per family-epoch (each epoch's pairs,
    anchored by the same rule within the epoch; a family-epoch without network/anchor falls back as above, using
    the family's all-epoch offset before δ_f). Reported for the 10 largest pairs too (column, no rule).
  - Output `tables/wa_family_offsets.csv`: fam, n_rows, component, d_f, offset_f (Vs 3.95), offset_cvm, delta_f,
    n_clear, anchor, flag (network / clear_only / default), epoch columns if used.
- **Step 2 — a-priori location-suspect flags** (catalogue only, before any stack):
  (i) ≥ 30 % of the family's selected rows have a co-detection (±12 s, raw times) with one other family whose location
  is > 15 km away (epicentral + depth, catalogue locations); (ii) catalogue depth < 20 km or missing.
  Duplicate clusters: connected components of the graph whose edges are family pairs with co-detection ≥ 30 % of
  either family's rows. Each cluster's Green's-function location = its best-located member: not suspect, then
  location_quality A > B > C, then most selected rows. If every member is suspect, members keep their own location.
  Families outside clusters: suspect-by-depth → Slab2 depth at their lat/lon (as Amendment 5); otherwise catalogue.
  Table `tables/wa3_location_flags.csv` (flags, cluster, cluster span km, GF lat/lon/dep, source).
- **Step 3 — rows v3** (`data/rows_v3/`): to_f = time + 0.5 + offset_f (family-epoch offset if step 1 says so);
  greedy ±2 s dedupe on to_f over all families, keep the highest P_own_pnet (same code as v2); GF location from step 2;
  along-strike zones = terciles recomputed on the GF locations of the kept families; same 69 NW majors, ±2 d
  windows, ≥ 30 LFE build rule, units (pooled, odd/even, zone1–3, families ≥ 20). Between: WA_BE, WA_BE_broad with
  the same offsets and the same ±2 s dedupe within each set. Prep any day not yet in `data/disp1sps_v2` (raw on disk,
  no download).
  Acceptance (Merlin): pooled n_nb2 = 0; mean pooled n_nb10 < 0.15; 25–35 % of rows removed by the dedupe. Failure of
  n_nb10 or the removed fraction does not stop the stacks; it caps the WA headline at "present, size provisional"
  and is reported.
- **Step 4 — re-stack** (tag `wa3_`): the built WA majors + `wa3_between` (WA_BE, WA_BE_broad) with
  `s13_stack2_v2.py` (guard2, UW.SP2 E veto), `--kmax 5`, headline ±8 s grid, columns `wide_` (±40 s) and `fix2_`
  (+2 s fixed), same WA station list (`s13_fill_queue.sh` WA_STA), `data/disp1sps_v2`, ≤ 6 parallel, OMP_NUM_THREADS=1.
  Pass rules (Merlin step 4), scored by `s13_v3_summary.py`:
  - **W1** mode of the pooled best lag (±8 s grid) over the built episodes in +1…+4 s, and ≤ 10 % of pooled units at
    ±8 s. **W1 is the headline gate**: pass → WA headline (A2 weighted mean per LFE and detections) stated with the
    rules below listed; fail → "WA present, size provisional", the fixed +2 s column given as a column.
  - W4 zones 2 and 3 each positive in > half of the ≥ 300-LFE episodes.
  - Size: pooled A2 weighted mean ≥ 3.5e12 N m per LFE (Merlin's expected 0.7–1.0× SVI). Fail → stated as "WA per-LFE
    lower than expected" with the WA/SVI ratio; it does not remove the headline.
  - W3 ±k-day z SD over pooled units × 10 shifts in 0.85–1.15; outside → σ scaled by that SD (as Amendment 5).
  - W2 (n_nb2 < 0.05), W5 (odd/even ρ > 0.4) and A1–A6 as in Amendments 1 and 5, reported.
  - Independent check (no rule): per-family lag table — catalogue-only offset_f gives a predicted lag in the
    **v2** stacks (origin time + 0.5) of offset_f + 2 s; observed = best lag of that family's grand stack over all v2
    episodes on the wide grid (joint fit, episodes weighted by n). Also the v3 family grand-stack lag (expected near
    +2). Families with ≥ 3 episodes and grand-stack |z| reported.
- **Step 5 — SVI strength test.** For 2004-07, 2005-09, 2012-08, 2013-09: pooled rows (as in `data/rows/`) with
  P ≥ 0.90, split into terciles of the catalogue `net_mad` within each episode (units nm_lo / nm_mid / nm_hi); thrust,
  200 draws, kmax 1, SVI station list, `data/disp1sps_v2`, neighbour pool = the episode's pooled rows; tag `nmt_svi_*`.
  Rule (Merlin): nm_hi / nm_lo ≥ 1.5 in ≥ 3 of 4 episodes → "per-LFE VLF moment scales with detection amplitude".
  Ratio within 0.67–1.5 in ≥ 3 of 4 → "does not scale with detection amplitude at fixed P" (Merlin's change-my-mind).
  Otherwise "mixed". nm_lo ≤ 0 with nm_hi > 0 counts as ≥ 1.5 (stated if it happens). P medians per tercile and the
  Spearman(P, net_mad) within P ≥ 0.9 are reported as composition checks.
- **Addendum 6a (2026-10-14 ~9:05 pm PT, after the first step-1 pass on the catalogue, before any row or stack of
  Amendment 6).** What the first pass showed (catalogue only): every co-detected pair is an exact template-alignment
  constant (pair residual rms 0.011 s, 302 of 325 ordered pairs with IQR ≤ 0.1 s; fam38–fam75 identical to 0.05 s in all
  three epochs, so the epoch rule does not fire). But selected rows alone split the network into 26 components; only
  11 could be anchored, and 59 families (33 % of rows) would get the default offset. Two changes, both fixed now:
  - (a) Pairs are counted on all NW catalogue rows with P_own_pnet ≥ 0.5 (the selected rows are a subset), keeping
    pairs with n ≥ 8 **and IQR ≤ 0.5 s** (a duplicate pair has a fixed offset; chance coincidences within ±12 s have
    an IQR of several seconds). The selected-rows-only network is kept as a column, and its offsets must agree with
    the primary on shared families (reported).
  - (b) Duplicate clusters: the co-detection flag (i) is symmetric, so it flags every member of a far-apart cluster
    and "not suspect" would never pick a best member. The best-located member is chosen among members **not suspect by
    depth** (ii), by location_quality A > B > C, then most selected rows. Anchor families still exclude both flags.
  - The acceptance SD is quoted two ways: over all anchor families, and over anchors in components with ≥ 2 anchors
    using a leave-one-out constant (a single-anchor component matches its anchor by construction). The rule applies
    to the leave-one-out SD.

## Amendment 7 — between-ETS stacks at candidate times (Merlin's scan fallback; written 2026-10-15 (Thu) ~2:15–2:35 am PT, before any row, prep or stack number below exists)
Source: `MERLIN_SCAN_NOGO_REVIEW.md` last paragraph ("If still failing: stack at candidate times …"), coordinator track B
(2026-10-15 ~1:50 am PT). Same machinery as the scale-up: `s13_stack2_v2.py` unchanged (guard2, instrument vetoes,
scrambled control = random time on the same day, ±k-day stacks), `s13_prep.py` (v2 recipe) into `data/disp1sps_v2`.
New files only: `scripts/s13_units_be3.py` → `data/rows_be3/`, `scripts/s13_be3_env.py` (PNSN origins), queue
`scripts/s13_be3_queue.sh`, summary `scripts/s13_be3_summary.py`; tags `be3_*`.
Facts used to write this (catalogue counts only, no stack): between-ETS (outside SVI majors ±2 d, 2003-01-01 → 2026-09-15)
data rows before the cross-family dedupe: P ≥ 0.75 11,173 (648 days); 0.5–0.75 14,336 (2,332 days); 0.3–0.5 23,120
(4,283 days); 0.1–0.3 61,378 (6,620 days). PNSN tremor epochs inside the SVI footprint polygon between ETS (2010+): 12,656
on 403 days, 10,926 with no P ≥ 0.75 row within ±5 min.

- **Question.** Do between-ETS LFE candidates of lower confidence, and between-ETS tremor with no catalogued LFE, carry the
  20–50 s thrust accompaniment that catalogued LFEs carry in ETS? Output: VLF moment per candidate by class and the
  between-ETS VLF:LFE ratio vs the ETS value (5.39e12 N m per LFE, 2008–26 weighted mean ≈ 12 × 4.5e11).
  **Scope:** this sees VLF moment only at times with 2–8 Hz activity (a candidate or a tremor window); it cannot see
  VLF slip with no 2–8 Hz activity.
- **Between ETS** = outside every SVI major [start − 2 d, end + 2 d] (Bostock-era majors < 2010 + `episodes_b.csv` SVI
  majors, as Amendment 1), 2003-01-01 → 2026-09-15.
- **Candidates.** `s13_common.svi_rows()` (backward primary & keep_xfam8; forward prod & keep_xfam8), series 'data',
  between ETS, P ≥ 0.1; origin to = t + family shift (`origin_shift_svi`, as all SVI stacks). Cross-family dedupe: greedy
  ±2 s on to over all of them, keep the highest P (Amendment 5 code `nms_dedupe`). Classes: C075 [0.75, 1], C050
  [0.50, 0.75), C030 [0.30, 0.50), C010 [0.10, 0.30).
- **Day sample (cost).** Calendar split into 7-day blocks from 2003-01-01; a block is eligible if it holds ≥ 1 between-ETS
  day; a random 20 % of eligible blocks (numpy default_rng(20261020)) is drawn once; the class units use the candidates on
  the between-ETS days of the drawn blocks only (an unbiased sample of between-ETS time; expected ~1,450 days). Prep =
  drawn days ±1 (for the ±1-day stacks).
- **PNSN tremor windows (unit PNSN_noLFE).** `catalogs/pnsn_tremor_cascadia_full.csv` epochs inside
  `tables/scan_footprint_polygon.csv`, between ETS (all such days, not only drawn blocks; 2010+), with no catalogued LFE:
  no row with P ≥ 0.75 within ±300 s of the epoch time (column unit PNSN_noCand: no row with P ≥ 0.1 within ±300 s).
  - Origin = peak of the network 2–8 Hz envelope inside [t − 150, t + 300] s.
  - The envelope comes from the `scan_hfenv.py` 1-s RMS of each station's vertical, for CN.PGC, CN.VGZ, CN.SNB, CN.NLLB,
    UW.LRIV. Each envelope is divided by its day median and shifted by −r/3.95 km/s, with r the hypocentral distance
    from the PNSN epicentre at Slab2 depth to the station. The ≥ 2 available stations are averaged.
  - Dedupe: greedy ±60 s on the origin, keeping the larger envelope peak.
  - Source location = the PNSN epicentre at Slab2 depth; P = 1.
  - Prep: those days ±1.
- **Anchor.** BE_broad (688 inter-ETS tremor-burst LFEs 2014–2026, P_same, origins as `svi_be_v2` BE_all) re-stacked in
  the same machinery (unit BE_broad); its days are already prepped.
- **Units** (all `full`: thrust + competitors, free MT; 500 draws if ≤ 3,000 rows else 200, as Amendment 1): C075, C050,
  C030, C010, PNSN_noLFE, PNSN_noCand, BE_broad. Column units (thrust only, 200 draws):
  - per class × era 2003–08 / 2009–13 / 2014–26;
  - per class × "on a PNSN tremor day in the footprint" vs not (2010+ only).
- **Stations:** CN.PGC, CN.VGZ, CN.SNB, CN.NLLB, UW.LRIV (the long-record near stations that carry the SVI stack weight).
  The stacker's own rule drops a station-unit with < 20 % Z coverage. `--kmax 1` (±1 day), headline ±8 s lag grid (c6_ and
  wide_ columns as usual).
- **Readout (fixed now).** Per unit: M0 per candidate (N m), z vs scrambled control, ±1-day moments, VR thrust/free.
  - "VLF accompaniment present" if z ≥ 3 and both ±1-day stacks have |z| < 3.
  - Ratio to the ETS value: ρ = M0 / 5.39e12 with σ_ρ = control SD / 5.39e12. VLF:LFE ratio R = 12 ρ, assuming each
    candidate's LFE moment is 4.5e11 N m (stated, not measured).
  - Wording: "ETS-like" if |ρ − 1| < 2σ_ρ and σ_ρ ≤ 0.35; "below ETS" if ρ + 2σ_ρ < 1; "above ETS" if ρ − 2σ_ρ > 1;
    otherwise "undetermined" (with the 2σ range).
  - Mixture reading (column, no rule): if a fraction f of a class are real LFEs with the ETS accompaniment and the rest
    carry none, then ρ = f. Reported next to the class's mean P.
- **Raw figures first:** stacks of the strongest station verticals per unit (black, thrust template red, control band
  grey, ±1 day blue); a class plot of ρ coloured by z, no error bars.
- **Compute:** ≤ 3 CPU workers (track A of the same coordinator turn holds ≤ 5), OMP_NUM_THREADS=1, < 100 GB; nohup
  queue with markers `PREP_be3_DONE`, `STACK_be3_*_DONE`, `BE3_QUEUE_DONE`.
- **WA (only if quick):** the same class stacks on WA-NW with the Amendment 6 per-family offsets. Not started in this turn
  unless track A and the SVI stacks finish early (needs a WA prep of the drawn days); stated as such otherwise.
