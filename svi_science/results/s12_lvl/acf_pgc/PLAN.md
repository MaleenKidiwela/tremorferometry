# PLAN — S12c: post-ETS LVL change under PGC from ambient-noise autocorrelation

Written Tue 2026-10-13 ~10:50 am PT, **before any ACF was measured**. Changes after that are dated and labelled in NOTES.md.

## 0. Question
Receiver functions at PGC show a short post-ETS change in the low-velocity layer (LVL) (Gosselin et al. 2020, reproduced in
`../rf_gosselin`: SV CC drop −0.12, days 0–15 after the ALBH midpoint, detectable to ~60 d, LVL-top contrast up = Vs up).
Teleseisms are sparse. Can continuous single-station noise correlations (ACF) follow the same change day by day?

## 1. Expected arrivals (step 1, done: `scripts/acf00_predict.py`, `tables/predicted_lags.csv`)
Observed RF Ps delays (from the RF depth picks, p = 0.06 s/km): PGC top 4.29 s, bottom 5.54 s; VGZ 3.35 / 5.03 s;
SNB 5.28 / 6.39 s. Zero-offset two-way times in 4 models that all honour these Ps delays
(RF model Vp/Vs 1.75; CVM v1.7 south_VI rescaled; LVL Vp/Vs 2.0 and 2.3):

| station | interface | PP (s) | PS = P down, S up (s) | SS (s) |
|---|---|---|---|---|
| PGC | LVL top | 9.7–10.9 | 13.9–15.0 | 18.0–19.1 |
| PGC | LVL bottom | 12.8–14.1 | 18.1–19.4 | 23.4–24.6 |
| VGZ | top | 7.2–8.5 | 10.4–11.7 | 13.6–14.9 |
| VGZ | bottom | 11.0–12.8 | 15.9–17.6 | 20.7–22.4 |
| SNB | top | 12.4–13.4 | 17.4–18.5 | 22.5–23.5 |
| SNB | bottom | 15.0–16.2 | 21.1–22.3 | 27.2–28.3 |

(The brief's 9–11 s PP and 16–19 s SS match the LVL **top**; the bottom is ~3 s (PP) and ~5 s (SS) later.)

**Forward target** (`tables/predicted_targets.json`). +2.5 % Vs in a 3.5 km layer at the LVL top (Gosselin size,
+0.05–0.1 km/s ≈ 1.7–3.3 %), Vs in the layer 3.0–4.1 km/s:
- everything reflected **at or below the LVL bottom** arrives earlier by **SS −40 to −60 ms** (Vs 3.0 → −58 ms),
  **PS −25 to −35 ms**, **PP −12 ms** if δVp = δVs/2 (0 if Vp unchanged);
- the LVL-**top** reflection does not move, but its amplitude changes (impedance contrast up ~2.5 %, i.e. ~+10–25 %
  of a contrast of ~0.1–0.2);
- as an apparent dv/v in the SS-bottom window (~24 s): **+0.17 to +0.25 %**; PP-bottom (~13.5 s): +0.09 %.
- The first 15 days carry the full change; it decays by ~60 d (RF bins).

## 2. Data (step 2) — all on disk, no download needed
Raw day-mseed, 3-C, already on disk for all three stations (from the LFE-extension work):
`lfe_extension/traces/` (2014–2026, HH 100 Hz) and `lfe_extension/traces_backward/` (2003–2013, BH 40 Hz / HH).
PGC 2003–2026 ~8,000 days (gaps: 2010 181 d, 2012 297 d, 2017 144 d); VGZ ~7,900; SNB 2003–2021-03 ~6,000.
Download size needed: 0 GB. Inventory in NOTES.md.

## 3. Preprocessing (fixed; `scripts/acf01_daily.py`)
Day → decimate to 10 Hz → rotate N,E to R = slab dip direction (RF dipdir: PGC 53°, VGZ 62°, SNB 58°), T = R + 90° →
per band (0.5–2 Hz and 1–4 Hz): bandpass → running-absolute-mean normalisation (window 1/f1, weight shared by the 3
components) → 1-h windows, dropping windows with gaps > 1 % or RMS > 4× / < 0.25× the day's median →
spectral whitening by a boxcar-smoothed spectrum (0.4 / 0.8 Hz wide, common to the 3 components; wide enough to keep
the spectral ripple from echoes at lags > 2.5 s) → ZZ, RR, TT autocorrelations (0–30 s) and ZR, ZT, RT cross-correlations
(±30 s), each normalised by the zero-lag autos → daily mean. Stored per station-year in `data/<NET.STA>/<year>.npz`.
No instrument correction: the ACF is zero-phase and whitening removes the amplitude response; the three components
share one sensor, so cross-component phase cancels too.

## 4. Pilot (step 3): PGC (+ VGZ, SNB) 2018-01 … 2021-12 (ETS 2018-07, 2019-08, 2019-11, 2020-11 weak, 2021-02, 2021-10)

### 4a. Visibility and stability (decided before measuring any time change)
- Long stack = mean of all pilot days. Plot 0–30 s per component and band with the predicted bands from §1.
- **Stability**: four yearly stacks (2018, 2019, 2020, 2021); mean pairwise CC of the yearly stacks inside each
  window. Same number in a late control window (27–30 s) for comparison.
- **Peak**: lag of the largest envelope peak inside each predicted band, its size relative to the median envelope in
  ±4 s around it, and whether the same peak (±0.3 s) is the largest one in ≥ 3 of the 4 yearly stacks.
- **Pass rule P1 (= "the pilot shows the arrival", the scale-up condition):** at PGC, for at least one band and
  component group (ZZ for PP; RR/TT for SS; ZR for PS), the LVL window (§4b) has yearly-stack CC ≥ 0.5 AND a peak in the
  predicted top or bottom band that is the largest in ≥ 3 of 4 years. If the windows are stable (CC ≥ 0.5) but no
  distinct peak shows, the differential measurement is still possible; report it as "stable coda, no distinct LVL echo"
  and still run §4b–§5 on the pilot, but scale-up then needs the pilot noise level (P2) to be within 2× of the target.

### 4b. Measurement (differential, immune to shallow change)
Windows for PGC (VGZ/SNB shifted by their own §1 bands), lag in s:

| arrival | pre window (above the LVL) | LVL window (at/after the LVL bottom) |
|---|---|---|
| SS (RR, TT) | 13.0–17.5 | 22.5–26.0 |
| PS (ZR, positive and negative lags measured separately) | 10.0–13.5 | 17.5–20.5 |
| PP (ZZ) | 6.0–9.5 | 12.5–15.5 |

For a test stack s and a reference r: time shift δt_W (cross-correlation, cosine-interpolated sub-sample peak) and
CC_W in each window. **Differential delay** Δt = δt_LVL − (δt_pre / t̄_pre) · t̄_LVL (removes any change that stretches
all lags alike, i.e. shallow/whole-column change). Also reported: stretch dv/v per window (`src/tremorferometry/dvv.stretch_dvv`)
and ΔCC = CC_LVL − CC_pre (decorrelation, the ACF analogue of the RF CC drop). Sign: a Vs **increase** in the LVL gives
Δt < 0 (earlier) in the LVL window.
**Control windows:** the same differential between two windows that are both above the LVL top
(SS: 9–12.5 vs 13–17.5; PP: 3–5.5 vs 6–9.5) must show no ETS effect.

### 4c. Superposed epoch on the ALBH ETS midpoints (step 4)
Midpoints = `../rf_gosselin/tables/ets_mid_SVI_ALBH.csv` (same list as the RF work; all 30 in the full record).
- Relative days k = −60 … +120. For each k: stack the daily ACFs at mid + k over episodes (daily, and with a ±3 d and
  ±7 d boxcar for the curve). Reference = stack of all days k ∈ [−60, −10].
- **Epoch metric (primary):** Δt(post) = Δt measured on the stack of days 0–15 against the stack of days −30…−1
  (the same windows as the RF test). Secondary: days 15–30, 30–60 (the RF decay bins).
- **Scrambled control:** 2,000 draws; each midpoint replaced by a random date with the same day of year ±30 d in a
  random year with data, not within ±60 d of any real midpoint; same metric. pct = fraction of draws with Δt at least as
  negative (one-sided, Vs-increase direction). Also the two-sided pct.
- **Primary channel (fixed now, chosen on stability only, which is blind to ETS):** the band and component group
  with the highest yearly-stack CC in the PGC SS LVL window; RR and TT measured separately and averaged.
  All other channels are secondary and are reported with the number of channels tried.
- **Pass rule P3 (detection, full record):** primary-channel Δt(post 0–15) < 0 with scrambled pct ≤ 0.05, and the control
  windows within their scrambled 90 % band. P2 (sensitivity): scrambled-control σ of Δt(post) scaled to the full record
  (×√(N_pilot/N_full)) ≤ half the target (≤ 25 ms SS). A non-detection is reported with the target beside it and the
  scrambled σ, i.e. as a bound.
- Comparison with the RF curve (`../rf_gosselin` recovery tables when they land): same bins, side by side.

## 5. Scale-up (step 5, only if P1 passes)
2003–2026 for PGC, VGZ, SNB (preprocessing ~2.5 s per station-day per CPU → ~2.5 h at 8 CPUs). Same §4b–4c with all
30 episodes. Recovery time constant: fit Δt(k) = A·exp(−k/τ) for k ≥ 0 to the ±3 d curve, bootstrap over episodes;
also the scrambled control for A.

## 6. Plots and words
Line graphs; points coloured by CC (no error bars); scrambled control as thin grey 5/95 % lines. Lead with what the data
show. "Scrambled control", never "null". Times in PT.

## 7. Resources
≤ 8 CPUs, < 40 GB RAM, nohup jobs with DONE markers, resume block in `scripts/resume_after_reboot.sh`. Kill by PID only.

## 8. Amendment A (Tue 2026-10-13 ~10:55 am PT, after the Merlin review, BEFORE any time-change measurement)
Only the long-stack shape of PGC 2018 had been looked at (test of the plotting script). Changes:
1. **Tremor hours are masked in the primary result.** Tremor under PGC is an ETS-locked noise source at the
   interface; its multiples land at LVL lags and move with the front, which the scrambled control cannot catch.
   Primary = tremor-free hours: drop an hour if a PNSN tremor epicentre lies within 60 km of the station within
   ±30 min of the hour (catalog 2010-08+), OR the hour's Z RMS at 2–4.5 Hz is above the 90th percentile of hours within
   ±15 days of the same day of year (works before 2010 too). All-hours is reported as secondary. Acceptance: tremor-free
   and all-hours Δt agree within the scrambled σ → source-independent; if they disagree, report a source effect, not LVL.
   Source diagnostics added: LVL-window envelope amplitude days 0–15 vs reference, and positive- vs negative-lag
   asymmetry of ZR/ZT.
2. **Two forms of the differential.** A change in a layer above both reflectors delays every deeper echo by the same
   additive δt. So the primary is the plain difference **Δt = δt_LVL − δt_pre**; the stretch-scaled form of §4b is
   reported beside it, and the gap between them is the systematic.
3. **P1 honesty.** 2 bands × ~6 component groups × top/bottom ≈ 24 trials, so one "stable peak" can pass by chance.
   P1 now also needs the same echo (same component group, same band, the same polarity) at VGZ or SNB at their own
   predicted lags (the station offset is the fixed prediction). The number of trials is reported.
4. Delay measurement on stacks FFT-upsampled ×10 (100 Hz) before the CC peak (already in `acf_common.measure`).

## 9. Amendment B (Tue 2026-10-13 ~11:40 am PT, coordinator request; written before any amplitude was measured)
Gosselin measured a change in the Vs CONTRAST at the LVL boundary (converted-wave amplitude), not a delay. Add an echo
AMPLITUDE and waveform-similarity stage (`scripts/acf05_amp.py`, queue `scripts/acf_amp_queue.sh`, after the full run):
- Echo windows = predicted LVL top and bottom bands ±0.3 s (§1) for S–S (RR, TT), P–S (ZR±), P–P (ZZ); reference window
  = the §4b "pre" window (above the LVL). Long-term echo L = stack of all days of the record.
- Metrics for a stack s: **projection amplitude** a_W = ⟨s, L⟩_W / ⟨L, L⟩_W; **relative echo strength**
  A = a_echo / a_ref (primary; immune to anything that scales the whole ACF); RMS ratio rms_echo / rms_ref (secondary);
  **similarity** CC(s, L) in the echo window, and dCC = CC_echo − CC_ref.
- Epoch: days 0–15, 15–30, 30–60, 60–120 vs −30…−1 (ratio A_post / A_pre − 1, CC difference); ±3-d curve from −60 to +120
  against days −60…−10; the same scrambled control (2,000 draws for bins, 200 for the curve).
- **Primary: PGC, 0.5–2 Hz, S–S TOP echo, tremor-free, relative strength change days 0–15; two-sided** (sign of the
  reflection change for a faster LVL depends on whether its top or bottom contrast moves; Gosselin report the top
  contrast up and the bottom down at VGZ). Pass: two-sided scrambled pct ≤ 0.05 and the same sign with all hours.
  Expected size: +2–3 % Vs in the LVL changes a ~10–20 % Vs contrast by ~15–25 % of the echo amplitude.

## 10. Amendment C (Wed 2026-10-14 ~2:50 pm PT) — confirming test of the "weeks before ETS" lead. Written BEFORE running.
**The lead (found post hoc, so it cannot confirm itself).** PGC 0.5–2 Hz, days −21…−1 vs −60…−25 (`acf07`, NOTES 10-14):
S–S LVL-bottom echo strength −19.5 % tremor-free (σ 9.8 %, pct 0.014) / −14 % all hours; P–P Δt +3.7 ms (pct 0.996 in the
Vs-up direction, i.e. +) and P–P top-echo strength +24 % (pct 0.003). In the RF (`../rf_gosselin` §6–8), CC_SV is ~+0.04–0.07
above its long-term mean 1–3 weeks before the midpoint at PGC (outside the band), VGZ, SNB, LZB, PFB, NLLB, YOUB and LON.
**Independent data for the test** (not used to find the lead): ACF at VGZ (and SNB if usable); RF at the stations whose
pre-ETS window was never computed or read.
**Same-month scrambled control (both parts):** each real midpoint is moved to a random other year with data, same calendar
date ± U{−15…+15} d, rejected if < 60 d from any real midpoint or outside the station's data (ACF: also < 30 d from a date
already drawn). ACF 2,000 draws, RF 20,000.

**C1 — ACF (primary).** Script `scripts/acf08_lead_test.py` (copy of acf07's loading, mask, windows, `differential` and echo
features; only the bins and the ±15 d control differ). Station **VGZ**, **0.5–2 Hz**, **tremor-free** hours, metric = S–S
LVL-bottom echo strength change A_bottom (RR, TT averaged), **days −21…−1 vs −60…−25**.
Pass: A_bottom < 0 with one-sided same-month pct ≤ 0.05, AND A_bottom < 0 with all hours.
Power: VGZ's scrambled σ for the bottom-echo strength was ~9 % (0–15 vs −30…−1 bins), so a PGC-size −20 % is ~2σ:
a real effect of that size passes about half the time. A fail is not proof of absence.
SNB: usable only if its 0.5–2 Hz tremor-free S–S LVL-window yearly-stack CC (RR and TT, `tables/visibility_full_tremorfree.csv`,
computed blind to ETS) is ≥ 0.5. It is 0.17 / 0.11, so SNB is **not usable** and is reported descriptively, not scored.
Secondary, reported and not scored: VGZ P–P Δt (lead sign +) and P–P top-echo strength (lead sign +); VGZ 1–4 Hz same three
metrics; PGC rerun with the ±15 d control (not independent; shows how the control choice moves the lead's pct).
0.2–1 Hz: the stored ACFs hold only 0.5–2 and 1–4 Hz (`acf01_daily.py` bands), so this band would need a new pass over
the raw day files (~2 h per station at 3 CPUs) and new windows; **not part of this test**.

**C2 — RF pre-window on unused stations.** Script `../rf_gosselin/scripts/rf07_samemonth.py`. Metric
H = mean CC_SV (`cc_sv_clean`) in [−21, 0) d minus mean in [−60, −25) d, days from the nearest real midpoint.
Primary set (local tremor ≥ 40 events per episode within 30 km, RF §9): **GOBB, MGB** (ALBH midpoints) and **GNW, SP2**
(SEAT midpoints). Secondary, not scored: OFR, WISH, NLWA (no local tremor: a pre-ETS high there would point to a non-LVL
cause). Honesty: these stations' 15-day Δ and recovery figures were already read (RF §7–8 "nothing outside the band near
the ETS"; GNW's 15-d Δ is +0.075, which if anything argues against a pre-ETS high there), but H was never computed for them.
Pass: Stouffer-combined one-sided same-month p (H > 0) over the 4 primary stations ≤ 0.05, AND H > 0 at ≥ 3 of 4.
PGC's H with the same control is also given, as the reference the lead came from.

**Overall reading (fixed now):** lead **confirmed** if C1 and C2 both pass; **partly supported** if one passes; **not
confirmed** if neither (with the power caveat above).
