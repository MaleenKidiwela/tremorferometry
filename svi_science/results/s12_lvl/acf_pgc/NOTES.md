# NOTES — S12c acf_pgc (decision log, newest last; times PT)

## Tue 2026-10-13
- ~10:45 am. Read rf_gosselin RESULTS/target.json (PGC ΔCC_SV −0.12 days 0–15, to ~60 d; LVL top 35.5 km, bottom 47 km
  by RF depth picks), lfe/RESULTS (no LVL change in LFE coda; no LVL-base reflection at 2–8 Hz in LFE records), and the
  old ACF script (`archive/superseded_or_support/scripts/autocorr_dvv.py`: Z only, 2–8 Hz, 0–10 s lag, whitened, stretch
  2–4 s vs all-time mean; the B003 "expansion" from it was shallow). New here: lags to 30 s, 3-C, 0.5–2/1–4 Hz,
  differential windows, epoch stacking.
- Step 1 (`acf00_predict.py`). Recovered Ps delays from the RF depth picks (p 0.06): PGC 4.29/5.54 s, VGZ 3.35/5.03 s,
  SNB 5.28/6.39 s. (The brief quoted VGZ Ps 3.8–4.8 s, which are the raw RF pulse times in rf NOTES; the depth picks give
  3.35/5.03. Either way SS ≈ 4.5 × t_Ps.) PGC: PP top 9.7–10.9 / bottom 12.8–14.1 s; SS top 18.0–19.1 / bottom 23.4–24.6 s.
  Target for +2.5 % Vs in 3.5 km: SS −40…−58 ms, PS −25…−35 ms, PP −12 ms (if δVp = δVs/2).
- Step 2 inventory (no download needed; 0 GB). Raw 3-C day-mseed on disk:
  | station | 2003–2013 (`lfe_extension/traces_backward`) | 2014–2026 (`lfe_extension/traces`) |
  |---|---|---|
  | PGC | 365/yr except 2010: 181, 2011: 314, 2012: 297, 2013: 305 | 365/yr except 2017: 144, 2026: 258 (to Sep) |
  | VGZ | 300–366/yr (2003: 300, 2011: 272) | 365/yr (2017: 311) |
  | SNB | 256–366/yr | 2014–2020 full, 2021: 86 (closed 2021-03) |
  Also 20 Hz 1–8 Hz caches (`lfe_extension/cache20hz*`) — not used (1 Hz high-cut would hurt the 0.5–2 Hz band).
  If ever needed: a 3-C 20 sps download 2003–2026 is ~3 × 8,400 d × 3 ch × 1.7 M samples × ~1.5 B ≈ 2 GB per station.
- Preprocessing timing: ~2.5 s per station-day per CPU; 8 workers → ~70 s per station-year.
- ~10:52 am. First run wrote daily means; switched to HOURLY storage (needed for tremor-hour masking). v0 file moved
  to `data_v0_daily/` (not deleted). Pilot queue relaunched (PGC, VGZ, SNB 2018–2021).
- ~10:55 am. Merlin review → PLAN amendment A (tremor-free hours primary; plain additive Δt primary, stretch-scaled
  beside it; P1 needs the echo at a second station; ×10 upsampling). Applied before any time-change measurement.
- First look (PGC 2018, long stack only, while testing the plot script): 0.5–2 Hz RR/TT show their largest envelope
  peaks at 18.2 s (predicted SS top 18.0–19.1) and 23.8 s (predicted SS bottom 23.4–24.6), ×2.3–2.8 the local median.
  ZZ: 12.9 s ×2.5 (PP bottom 12.8–14.1). Provisional, one year, one station.
- ~11:00 am. Debug run of acf03_epoch.py on 2018–2019 (3 ETS, 20 scrambled draws) to check the code; no plan change
  followed from it. Noted: SS window scrambled σ of Δt(0–15 d) is only ~1–3 ms (stable coda), PS (ZR, ZT) σ 20–50 ms.
- ~11:03 am. PGC pilot visibility (2018–2021, `tables/visibility_pilotPGC_*.csv`). 0.5–2 Hz: yearly-stack CC 0.89–0.99
  in every window, including the late 27–30 s control (0.94) — the whole ACF is a very stable ringing coda, so
  stability alone does not single out an LVL echo. Largest envelope peaks inside the predicted bands, same lag in 4/4
  years: RR top 18.2 s (×2.6), RR bottom 23.8 s (×2.8), ZZ bottom 12.9 s (×2.7), ZR± bottom 17.9 s (×1.7–1.8).
  The ratios (×2–3 of the local median envelope) are not much above other peaks in the same trace (e.g. ZZ 21 s,
  24.5 s), so "an echo at the predicted lag" is plausible but not unique. 1–4 Hz is less stable (CC 0.4–0.9).
  Tremor-free vs all hours: no visible difference in the long stacks.
- **Primary channel fixed (PLAN §4c rule, stability only): 0.5–2 Hz, SS group (RR, TT averaged), tremor-free hours**
  (SS LVL-window yearly CC RR 0.90 / TT 0.97 at 0.5–2 Hz vs 0.78 / 0.45 at 1–4 Hz).
- ~11:05 am. Pilot analysis queue launched (`scripts/acf_pilot_ana.sh`, waits for PILOT_PRE_DONE).
- ~11:12 am. P1 check with VGZ/SNB (`tables/visibility_pilot_*.csv`): the 0.5–2 Hz RR LVL-bottom peak shows at all three
  stations at their own predicted lags, same lag in 4/4 years: PGC 23.8 s (pred 23.4–24.6, ×2.8), VGZ 21.0 s
  (20.7–22.4, ×2.7), SNB 28.6 s (27.2–28.3 + 0.3 margin, ×3.3). **P1 passes → scale-up launched** (`acf_full_queue.sh`,
  gated on PILOT_ANA_DONE to keep ≤ 8 CPUs).
  Honesty check (post hoc, `acf02b_peak_rank.py`, `tables/peak_rank_pilot.csv`): a band of the same width slid over
  6–30 s finds a peak at least as strong (and as repeatable) 15 % (PGC), 22 % (VGZ), 6 % (SNB) of the time. So each
  station alone is weak; jointly ~0.2 % if independent (~2 % after ~10 trials). Plausible LVL-bottom SS echo, not proven.
- ~11:16 am. First PGC pilot epoch run (tag `pilot`, 6 ETS): the CC-peak delay cycle-skipped (scrambled draws with
  ±100 ms jumps; RR scrambled sd 70 ms). Fixes, made on the scrambled-draw behaviour only: maxshift 0.5 → 0.2 s, then
  the primary delay estimator replaced by a linearised (Taylor) delay that cannot skip cycles (synthetic test recovers
  0, 20, 50, −60, 150 ms exactly); CC-peak delay kept as `ddt_ccpk`; robust scrambled σ (1.4826·MAD) added.
  The `pilot` tag result for the record: primary SS (0.5–2 Hz, tremor-free) Δt(0–15 d) = −0.3 ms.
  Gain-proof RMS tremor proxy (divide by a 61-day running median) added for the multi-era full record.
  Pilot analysis rerun as tag `pilot2`.
- ~11:20 am. Full-record queue restructured (4 epoch runs in parallel; acf04_compare at the end) and relaunched
  (PID in logs/full_queue.pid; preprocessing started 11:26 am; ETA ~2 pm PT).
- ~11:26 am. Pilot2 done for all 3 stations (`tables/summary_pilot2.csv`, `figures/compare_rf_pilot2_CN.PGC.png`).
  Primary: PGC 0.5–2 Hz S–S tremor-free Δt(0–15 d) = −0.4 ms, σ 1.0 (target −56). 1–4 Hz PGC: right sign in S–S/P–S/P–P
  at ~2σ (candidate; P–P control also moves). SNB 0.5–2 Hz S–S −11 ms but control −8 ms and amplitude +63 % → source.
  RESULTS.md pilot section written.
- Note on the RF curve in compare figures: plotted at its `centre` column as given by rf05 (w15 window); its minimum sits
  at ~+25 to +35 d on that axis, so check rf05's centre convention before reading timing differences.
- ~11:40 am. Coordinator: Gosselin measured a CONTRAST (amplitude) change, so a timing null leaves it open. PLAN
  amendment B written first; `acf05_amp.py` (relative echo strength by projection on the long-term echo ÷ the same above
  the LVL; RMS ratio; CC with long-term echo; LVL top and bottom; S–S, P–S, P–P; same epoch bins, curve and scrambled
  control). Debug run on the pilot (10 draws) only to test code: scrambled σ of the S–S top strength ~3–10 % for
  6 ETS. Queue `acf_amp_queue.sh` (waits FULL_QUEUE_DONE) launched; `acf06_report.py` then appends an auto table to
  RESULTS.md. Resume block 17 updated.

## Wed 2026-10-14
- ~1:30 pm PT. Full record + amplitude stage done (AMP_QUEUE_DONE); auto table appended by acf06_report.py.
- ~1:40 pm PT. POST HOC (after the pre-set results were known; coordinator request from rf_gosselin §6/§11): new
  `scripts/acf07_posthoc_rfwindow.py` — same data/mask/estimators/scrambled rule, other day bins (0–24, 6–15 vs −30…−1;
  0–24 and −21…−1 vs −60…−25; 0–24 vs −21…−1). PGC both bands, both masks, 2,000 draws, ~2.5 min each.
  Out `tables/posthoc_rfwindow_CN.PGC_*.csv`. Pre-set 0–15 values reproduce exactly.
- ~1:55 pm PT. Found three VGZ epoch runs (b05_2 tremorfree, b1_4 none/tremorfree) had started 12:2x pm, before VGZ
  2024–26 was recomputed (12:41–12:43 pm) → 20 ETS. Rerun under tag `fullvgz2` (22 ETS); originals kept. Same reading.
- Found: `P3_pass` in acf04_compare.py ignores the control-window condition (not fixed; RESULTS reading scores P3 by hand).
- ~2:00 pm PT. '### Reading' written in RESULTS.md; top block items 1 and 3 edited (dated).
- ~2:50 pm PT. PLAN §10 Amendment C written (lead confirming test: C1 ACF VGZ 0.5–2 Hz tremor-free S–S bottom-echo strength,
  days −21…−1 vs −60…−25, same-month ±15 d control; C2 RF pre-window H on GOBB/MGB/GNW/SP2). SNB declared unusable from its
  ETS-blind yearly-stack CC (0.17/0.11 < 0.5). 0.2–1 Hz not in the stored ACFs → not tested.
- ~2:52 pm PT. `acf08_lead_test.py` (acf07 copy: lead bins + same-month ±15 d control, 30-d spacing between drawn dates) and
  `acf_lead_queue.sh` (3 lanes, OMP 1); 7 runs, ~1 min each, LEAD_QUEUE_DONE 2:56 pm. 0–15 d values reproduce acf03/acf05.
- Result: VGZ bottom echo +5.3 % (σ 8.4, lead sign −) → C1 fail; C2 fail (combined p 0.90) → lead NOT confirmed (RESULTS "Lead test").
  PGC lead with the ±15 d control: pct 0.006 (bottom), 0.002 (P–P Δt), 0.001 (P–P top).
