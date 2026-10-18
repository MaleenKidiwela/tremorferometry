# G — Do large distant or regional earthquakes trigger deep LFEs under SVI and Washington? PLAN (pre-registered)

Written Sun 2026-10-11 ~5:45 am PT, **before any G statistic was computed**. Looked at while writing: file schemas
and date ranges of the catalogues, the two ComCat query results (row counts only: 356 global M ≥ 7, 139 regional
M ≥ 5.5), the raw-trace inventory (PGC broadband 40 Hz on disk 2003–2026, PB boreholes 100 Hz) and the station
metadata file. No LFE count near any earthquake was looked at. Changes after this point go to "Revisions" with a PT
time and the reason. Science only: no catalogue file is changed. Scripts `svi_science/scripts/g_trigger_*.py`;
outputs here (`tables/`, `figures/`). CPU ≤ 8 cores, RAM < 60 GB. "Scrambled control" = the same detector on the
same day with each station site's data shifted to another time of that day. Times shown to the reader are PT;
computations use UTC.

## 0. Data and usage rules

**SVI LFEs 2003-01-01 … 2026-09-15** = `results/r23/tables/events_2003_2026.parquet` (built under the
RESULTS_BACKWARD usage rules: primary detector, era rows P ≥ 0.90, frozen rows P ≥ 0.75, counted once; forward
2014+ rows `keep_xfam8` and P_same ≥ 0.75; data and scrambled-control rows through the same steps). Network
segment `seg`: CN4f, ERA (2003-02 … 2006-06), B2f, B3f, S15f, FWD. Multi-site support `n_z2` (sites with z ≥ 2):
backward rows carry it; forward rows get it joined from `lfe_extension/results/refine/feat/prod/<date>.parquet`
(key series, fam, t). Observable family-days: `r23/tables/famdays_observable.parquet`. Episodes 2003–2026:
`r23/tables/episodes_2003_2026.csv` (145; pre-2010 = 6 Bostock windows + 39 envelope minors, the R23 fallback,
flagged `validated == False`).

**WA LFEs (WA-NW) 2007-01-01 … 2026-09**: `lfe_extension_wa/results/wa_prod/conf/P/prod/P_*.parquet`, net core,
reg NW, P_own_WAVE_PNET ≥ 0.75 (the post_v2 headline), counted once = dedup variant D False
(`post_v2/dedup/flags`); series data, scrambled control c1 and time-reversed control rev (column). `n_z2` joined
from `mf/prod/det_<date>.parquet` (key date, net, series, fid, k). Eras E1 2007–2011, E2 2012–2014, E3 2015–2026
(fixed core network per era). Observable family-days = the NW family rows of `post_v2/dedup/rates_daily_v2_dedup.parquet`.
Episodes: `episode_rates/episodes_nw_final.csv` (118). Before 2009-08-06 only the 3 WECC major windows are listed
→ between-episode status unknown 2007-01 … 2009-08-05: those earthquakes are a separate column, not headline.
SW / Arm B / Olympic A are not used (no usable between-episode headline).

**Tremor**: PNSN epicentres (`catalogs/pnsn_tremor_cascadia_full.csv`, 2009-08 … 2026-05-27) and the pre-2010 env
flag only through the episode lists above, plus a descriptive column: PNSN epicentres within 100 km of the region
reference in each earthquake window vs its random windows.

**Earthquakes**: USGS ComCat FDSN (queried 2026-10-11 ~5:40 am PT, `dl/`): global M ≥ 7.0 2003-01-01 … 2026-10-11 and
M ≥ 5.5 within 1,500 km of (48.6 N, 123.6 W). Merged, duplicates removed by id. **Clusters:** earthquakes whose
origin times are within 6 h of each other form one trial, represented by the one with the largest peak dynamic
stress at the region (both listed). Reference points: SVI = mean lat/lon of the 115 families (`s0/families.csv`);
WA-NW = mean of the NW families' primary positions (`a_wa_between/families_wa.csv`). Distance Δ = great circle.

## 1. Arrival times, windows

- t_P = origin + first-arrival time (TauP iasp91, source depth from ComCat) at the reference point.
- Rayleigh slow end t_R = origin + Δ_km / 2.5 km/s. Love/Rayleigh start t_L = origin + Δ_km / 4.5 km/s.
- **W_inst (headline instantaneous window)** = [t_P, t_R + 60 min]. Columns: **W_surf** = [t_L, t_R + 30 min]
  (surface waves only), **W_6h** = [t_P, t_P + 6 h] ("hours"), **W_pre** = the same length as W_inst ending at t_P.
- Superposed-epoch view: counts in 30-min bins from t_P − 24 h to t_P + 24 h (figure; not a test).

## 2. Peak dynamic stress at the stations

- **Measured** (preferred): CN.PGC broadband from raw files on disk (`lfe_extension/traces/CN.PGC`,
  `traces_backward/CN.PGC`), response removed with `codex-LFE/results/full_history_metadata_v1/CN.xml` to velocity,
  band 0.02–2 Hz, 3 components; PGV = max |v| over components within [t_P, t_R + 30 min]. Strain ε = PGV / c with
  c = 3.5 km/s; stress σ = μ ε, μ = 30 GPa (Hill & Prejean convention). QC: ≥ 95 % of the window present on all 3
  components, no clipping (|counts| < 0.95 × 2^23), else not measured. Also the 2–8 Hz PGV and band RMS (for the
  eye check; not used for σ).
- **Estimated** (fallback, and a column for every event): surface-wave magnitude relation
  log10(A/T) = M − 1.66 log10(Δ°) − 3.3 (A in µm, T = 20 s; IASPEI Ms; as used by van der Elst & Brodsky 2010),
  M = ComCat magnitude (Mw ≈ Ms; saturation above ~8 noted), PGV ≈ 2π A/T, same ε and σ. Δ° < 2° is evaluated at 2°
  (relation outside its range; such events flagged `near`). The estimate is calibrated against the measured PGC values
  (median log ratio and scatter reported) but not rescaled.
- **σ used** = measured where it passes QC, else estimated; column `sigma_src`. WA-NW: measured PGC value is used
  for Δ ≥ 1,000 km (the ~150 km offset changes nothing at teleseismic range; `sigma_src = PGC_proxy`), the estimate
  at the WA reference otherwise.
- Stress tiers: **≥ 1 kPa (headline "strong")**, 0.1–1 kPa, < 0.1 kPa; ≥ 10 kPa listed.

## 3. Test 1 — instantaneous triggering (hours)

**Eligible (headline):** trials whose W_inst lies entirely outside every listed episode ± 2 d of that region, on an
observable day, inside the region's catalogue span with an episode list (SVI 2003-01-01 … 2026-05-27, the end of the
PNSN-based episode list; WA-NW 2009-08-06 … 2026-05-27). Trials during episodes → the ETS column (below).

**Per trial counts in W:** counted-once LFEs (data), scrambled control (c1), data with n_z2 ≥ 2, and the
artefact gauge: ALL candidates with any P (data and c1; SVI from `s0/analysis_table.parquet` + the backward
catalogue `primary` rows; WA from the P files) — this shows whether the detector itself reacts to earthquake energy.

**Random-date control:** for each trial, 2,000 windows of the same length starting at the same UTC time of day on
random days drawn from: the same network segment (SVI `seg`; WA era), outside every listed episode ± 2 d, observable,
and not within ±3 d of any listed earthquake. Seed 0.

**Statistics (stacked over eligible trials):** S = Σ N_data(W); X = Σ (N_data − N_c1); S_mz = Σ N_data,n_z2≥2.
Random stacks: one random window per trial, 10,000 stacks → p_R (S), p_S (X, the data increase beyond the scrambled
control's increase), p_mz (S_mz). One-sided, p = (1 + #{random ≥ observed}) / 10,001.

**PASS (per region, per stack):** p_R < 0.01 AND p_S < 0.01 AND p_mz < 0.01, and the eye check (§6) of the top
trials does not show the excess to be earthquake energy. Stacks: **headline = strong (σ ≥ 1 kPa)**; columns = all
eligible, 0.1–1 kPa, < 0.1 kPa, regional (Δ < 1,500 km) vs teleseismic, W_surf, W_6h, per network segment / era.
The same artefact-gauge numbers are reported for every stack.

**ETS column:** trials whose W_inst falls inside an episode: N in W vs a local background = all windows of the same
length tiling [t_P − 24 h, t_P) and [W end, W end + 24 h); per trial one random tile → 10,000 stacks → p_local;
scrambled control the same way. Same pass rule with p_local in place of p_R.

**Individual earthquakes named as "clearly triggered"** only if all hold: (a) N_data greater than every one of its
2,000 random windows; (b) Poisson p(X ≥ N_data | μ) < 0.01 / n_trials, μ = mean of its random windows (floor
1/2,000); (c) ≥ 2 of its LFEs have n_z2 ≥ 2; (d) N_c1 ≤ the 99th percentile of its random c1 windows and the
all-candidate c1 gauge not above its 99th percentile; (e) eye check passes. For ETS-column trials the local
background replaces the random windows.

## 4. Test 2 — delayed triggering (days 1–10)

- Day 0 = UTC day of t_P. Per trial: rate_after = Σ (data − c1) counted-once LFEs on days +1 … +10 / Σ observable
  family-days of those days; rate_before the same on days −10 … −1. Only between-episode days (outside every listed
  episode ± 2 d) count; a trial needs ≥ 5 such days on each side and day 0 between episodes.
- Statistics (mean over eligible trials): D_after = mean rate_after; D_diff = mean (rate_after − rate_before).
- Random-date control: pseudo-trials at random days of the same segment/era with the same eligibility, one per
  trial, 10,000 stacks → p_R(after), p_R(diff). Scrambled control: the same statistics on c1 alone must not pass
  (c1 p ≥ 0.01) and the data-only (not control-subtracted) statistics are columns. Multi-site: the same with
  n_z2 ≥ 2 rows.
- **PASS:** p_R(after) < 0.01 AND p_R(diff) < 0.01 AND multi-site p < 0.01 AND c1 alone not passing. Stacks:
  strong (headline), all, per segment/era; windows +1…+3 and +1…+10 (+1…+3 column).
- Between-episode rates are a control-subtracted bracket; nothing is compared across network segments (each trial
  is compared only with its own segment's random dates).

## 5. Test 3 — dose-response, and ETS onsets after big earthquakes

- **Dose-response (instantaneous and delayed):** Poisson model N_i ~ Poisson(μ_i (1 + β σ_i)), μ_i = trial i's
  random-window mean; β by maximum likelihood (β ≥ 0), profile-likelihood 99 % upper bound → "≤ X % extra LFEs per
  kPa". Spearman ρ(σ_i, (N_i + 0.5)/(μ_i + 0.5)) with a permutation p (10,000). Also the absolute 99 % upper bound on
  extra LFEs per strong trial: (S − 1st percentile of the random stacks)/n. **Response scales with stress** only if
  ρ > 0 with p < 0.01 AND β > 0 with the 99 % lower bound > 0. Delayed: the same with rate_after.
- **ETS onsets after big earthquakes:** onsets = list start of majors (SVI r23 majors 2003–2026; WA-NW majors
  2009-08 … 2026-05; the 3 WECC windows a column). Statistic K10 = number of onsets with ≥ 1 trial in the 10 days
  before (K30 for 30 days, column). Earthquake sets: strong (σ ≥ 1 kPa, headline) and all. Controls: (i) circular
  shift of all trial times by one random offset within the span (keeps clustering), (ii) uniform random dates;
  10,000 each. **Supported** if K10 > both at p < 0.01. (The scrambled control does not apply to onset dates; this
  test has no multi-site leg — stated, not substituted.)

## 6. Guards against the earthquake's own energy

1. Scrambled control in the same windows (p_S above; artefact gauge with all candidates).
2. Multi-site: n_z2 ≥ 2 subset (p_mz).
3. **Eye check:** 2–8 Hz record sections (SVI raw traces: PGC, VGZ/SNB/LZB/NLLB where present, PB boreholes,
   LRIV) for the top 5 trials per region by (N_data − μ) in W_inst plus every trial meeting (a)–(d) in §3, plus
   Haida Gwaii 2012 M7.8, Tohoku 2011, Chile 2010 and the largest regional event, as reference panels: envelope record
   section ± the window, P / S / Love / Rayleigh marks, LFE detections (data and c1) marked, and the broadband
   0.02–2 Hz PGC trace on top. Verdict per panel: "LFE-like bursts with moveout / modulated by surface waves",
   "earthquake coda only", or "nothing". The eye check can only veto, never pass.
4. **Literature (from memory, clearly labelled, not verified here):** Rubinstein et al. 2007 (Nature): 2002 Denali
   Mw 7.9 surface waves triggered tremor on Vancouver Island/WA, tremor pulses in phase with Love waves; Rubinstein
   et al. 2009 (JGR, Vancouver Island 1997–2008 survey): of the largest teleseisms, only Denali clearly triggered
   tremor on VI, with a stress threshold of order a few kPa reported. Other reports of Cascadia triggered tremor
   (e.g. 2012 Haida Gwaii, 2010 Chile, 2011 Tohoku) are not reliably recalled and will be stated as such. Denali (2002)
   is before our catalogue.

## 7. Regions and eras

SVI 2003–2026 and WA-NW 2009–2026 (2007-01 … 2009-08 column). Every stack is also given per network segment (SVI)
/ era (WA); random dates are always drawn inside the trial's own segment/era. No pooling of rates across segments;
the pooled stack is a sum of within-segment comparisons.

## 8. Outputs

`tables/eq_list.csv` (every earthquake, Δ, times in UTC and PT, σ measured/estimated, cluster, eligibility),
`tables/inst_trials_<region>.csv`, `tables/inst_stacks.csv`, `tables/delayed_*.csv`, `tables/dose.csv`,
`tables/ets_onsets.csv`; figures `figures/g_*.png`; RESULTS.md (plain words first).

## Revisions
- **R1 (Sun 10-11 ~5:55 am PT, before any LFE count was looked at).** `CN.xml` holds only the overall sensitivity
  for CN.PGC (no response stages), so §2 "response removed" is replaced by dividing by the epoch's sensitivity
  (flat-velocity broadband in 0.02–2 Hz; STS-2 / Trillium class). Clip test as written (|counts| ≥ 0.95 × 2^23).
- **R2 (same time).** A cluster trial uses the union of its members' windows (earliest P to latest Rayleigh + 60 min);
  σ = the members' maximum, as written. Random windows for the delayed test are all eligible same-segment/era days
  (not a 2,000 draw); the ±3-d exclusion around every listed earthquake is applied there too. The PNSN tremor column
  starts 2010-01-06 (the separate 2009 pull is not used). The dose bound is profiled on a β grid (0, 10^-5 … 10^3).
- **R3 (~6:05 am PT, after the tests ran).** ETS-column trials have only 24 local tiles for a ~2-h window, so rule (a)
  "greater than every tile" is a weak bar (smallest possible empirical p 0.04) and rule (b)'s Poisson law ignores
  burstiness. One trial met (a)–(d) (WA, Kermadec 2019-06-15). A **post-hoc** chance check was added
  (`g_trigger_posthoc.py`, onset definition fixed in its docstring before computing): how often does a burst *start*
  after ≥ 72 h of quiet, and how many starts are expected inside all the trial windows by chance. Labelled post-hoc.
- **R4 (~6:05 am PT).** Eye-check list extended with SVI trial 293 (2018-10-21 off-Vancouver-Island M6.5/6.8/6.5
  sequence), because the post-hoc onset list placed an SVI LFE start inside it; and with SVI 144 (2011-01-01 Argentina,
  the only between-episode SVI trial with N > 2 and empirical p < 0.01).
- **R5 (~6:12 am PT).** 4 earthquakes (incl. Kermadec 2019-06-15) had no PGC measurement only because of a
  100.0000076-Hz sampling-rate jitter that blocked the trace merge (a code fault, not data). Fixed (rates rounded to
  1 mHz) and every step re-run (PGV → tests → post-hoc → figures); the earlier tables were overwritten by the re-run.
- **R6 (~6:30 am PT).** The §2 clip test (|counts| ≥ 0.95 × 2^23) was wrong: PGC counts reach 4.3 × 10^7 without
  flat-topping (Haida Gwaii 2012), so the digitizer range exceeds 2^23. Replaced by a flat-top test (≥ 3 samples at
  the maximum |count|). 41 events had been sent to the estimate by the wrong test; all steps re-run.
- **R7 (~6:40 am PT).** After R6, 18 windows that cross UTC midnight read ~62 kPa. Cause: a 1-s glitch (~4 × 10^7
  counts, BHE) at 00:00:00–00:00:01 UTC at the end of the PGC day files, and day files that overlap the next day by
  ~1 min. Fixed by removing each piece's median, keeping the earlier piece where files overlap, and blanking
  midnight ± 3 s; all steps re-run. Remaining measured/estimated ratios outside 0.1–10: Haiti 2010 (×12), Illapel
  M7.0 aftershock (its window holds the M8.3 waves) and 2018-07-24 Oregon offshore M5.6 (×0.06); kept as measured.
