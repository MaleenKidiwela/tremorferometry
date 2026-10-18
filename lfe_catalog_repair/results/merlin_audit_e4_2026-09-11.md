# Merlin re-audit of E4 — 2026-09-11 ~04:45 (Fable 5.1; read-only; verbatim summary)

**Verdict.** Step 1's core finding STANDS and is the most important result E4 produced: with Z,N,E fixed and an empirical moveout,
Bostock's 2010 events register on our six-station ring at a raw network median 5.9 MAD (8-s template), 51% >=6 and 21% >=8, with the
rolled control at >=8 firing on 0.02% — the template route is open. But two downstream conclusions that drive the fleet decision are wrong
at the root: (i) the "false alarms per template per day" in RESULTS.md are TOTAL DETECTIONS, not false alarms; read off the quiet-day rate
and the only FA-viable threshold is ~8 MAD — "4-MAD-class threshold suffices in three regions" and "implied Bostock threshold 3.6 MAD" are
VOID; (ii) the prediction table compares the EXCESS statistic to the threshold when detection is on the RAW statistic -> every Mc ~0.5 mag
too pessimistic at 8 MAD; corrected, the model reproduces Bostock's 21-station array at 8 MAD with no threshold gymnastics. Decision: fix FA
accounting and raw/excess first, then DOWNLOAD Sept-2012 southern VI on the FULL 2012 station set (~17 sites, not 11) as the from-scratch
replication of Bostock's 28k-event catalog.

## Stands / falls
Stands: step 0; step 1 registration, held-out moveout (2.83 vs 3.17), 3-C > Z (direction), longer templates win; step 2 q-control and
borehole transfer (B011 0.85x, B926 1.33x); step 4 (seed only); step 6 station facts; the 21->7-station Bostock density drop (verified:
peak-day 4,670-7,783 in 2003-05 vs 2,551-3,640 in 2010-13; modal M 1.6-1.7 -> 1.7-1.8).
Falls: "FA 0-4/template/day at 8 MAD"; "4 MAD affordable"; "implied 3.6 MAD = D3's 3.5"; "Mc at 8 MAD none reaches Thomas density";
"0.078 anchor wrong by 7-13x" (coherent vs incoherent = two quantities); "injection 4-6x above real -> fine-tuning reopened"; "main front
reproduced at 3.3 km/day"; "secondary fronts not reproduced" (detector disabled by its own front model); "tidal growth reproduced"
(clustering artifact); "Sept-2012 pull = 11 stations".

## Flaws
F1 FA = total detections (FATAL for step 7, MAJOR for step 1 headline). e4_template_test.py:391-392, 452-459 count every declustered peak
   above T·MAD; e4_step1_model.py:62,65 report the MEDIAN over family-days. L=8: mean detections >=8 = 8.0/template/day, rolled control 1.8
   (medians 2/0 as published); 751 family-days: 6,037 detections >=8, 3,470 Bostock-matched, 1,342 rolled. QUIET days (doy 228-229,
   248-251): 0.1-1.4 detections/template/day at 8 MAD = rolled control -> noise FA at 8 MAD ~0.6/template/day (~10x Gaussian, heavy tails);
   at 6 MAD 45-115/template/day. 130 families: 8 MAD -> ~80 FA/day (precision >90% on ETS days); 6 MAD -> ~10,000 FA/day vs <=3,000 real;
   4 MAD -> >95% noise. Operating point 8 MAD (7.5 at best). Precision at 8 MAD on the 2010 ring ~92%; extras 1,200-2,100, ETS-concentrated.
F2 Prediction table EXCESS vs RAW (MAJOR, arithmetic). e4_prediction_table.py:139-140 uses log10(8/e_net_ref) with e = median-at-events
   minus random max-over-scan median (2.73). Raw median = excess + 2.73. In-sample contradiction: ring M1.8-2.1 L=8 excess 3.17, raw 5.89,
   frac>=8 0.206 — model predicts ~0. Every Mc too high by log10(T/(T-2.73))/0.38 = 0.48 mag at 8 MAD (1.31 at 4). CORRECTED 8-MAD
   equal-weight Mc: Puget 2023 1.71, N-Oregon 1.73, S-VI 2023 1.94, C-Oregon 1.98, NorCal 2025 1.95, N-VI 2.28 -> ~50%, 50%, ~17%, ~15%,
   ~17%, ~2% of Bostock's 2012 density. Bostock 21-station sanity: 2.10 -> 1.73 = his modal magnitude — the model reproduces his array.
F3 "Implied 3.6 MAD = D3's 3.5" false confirmation: inverted with p10=1.43 as "complete" (mode is 1.75) on 21-station geometry (his 2010-13
   catalog used 7); D3's 3.5 is a predicted network SIGNAL, a different quantity. A literal 3.6 threshold would make his catalog >90% noise;
   our data bound his 2010 false fraction at <=~8% (4% of his events below the random median 2.73 on our ring).
F4 Secondary-front detector disabled by its own main-front model: harness:286 passes first-EVENT Theil-Sen (line 229; -1.18 km/day,
   rho -0.12) into detect_secondary, which gates candidates to 1-10 km behind a front that never leaves the origin while the real one runs
   +35 -> -45 km. 0 vs 0.8 is uninformative. "Analysis bottleneck" NOT established either way.
F5 Main-front "3.3 km/day" is window-dependent: days 0-27 -3.32/-3.43; propagation phase (days 6-21) -4.53/-4.78; Mc2.0 days 6-21 -5.24.
   Thomas: 80 km in 11 days = 7.3 km/day (paper quotes ~5 and ~10); leading edge +35 km day 7 -> -45 day 20 ~6 km/day. PLAN GOAL -> "~7".
F6 Tidal "growth" = clustering artifact on a placeholder tide: Schuster assumes independent events; 28k events in 2,047 bursts; matched-n
   per-family test shows R_early 0.68 > R_late 0.42 (opposite). Withdraw.
F7 "0.078 wrong by 7-13x" conflates coherent and incoherent amplitude; FINE-TUNING CLOSED. e4_inject.py:60-85 measures S-window POWER
   excess (incoherent 0.57-1.46); the matched filter sees the coherent part: PGC excess 2.16 MAD x MAD-scale 0.046 ~ 0.10 in cc = 1.3-2.5x
   the 0.078 = step 1's 1.4-2.2x. A real LFE is ~15-25% template-coherent by amplitude. Anchored on coherent amplitude (0.1-0.2), injected
   recall 0.13-0.45 matches real 0.11-0.26 within 2x -> E3-audit acceptance PASSES.
F8 2012 southern-VI network is ~17 sites (f1_station_years.csv, 2012, 3-C, 366 d): C8 BPCB, PA01-03, SOKB (HH 100); PB B001, B003-B007,
   B009-B011, B926 (EH 100); CN NLLB, UW LRIV (BH 40); + ring 6. Co-located clusters PGC/B009/B010/B011/PA01-03; LRIV/B005-7. Rough
   SNR-weighted corrected prediction: raw ~7.5±1 at M1.95 -> Mc(8) ~2.0 -> ~25-35% of his 2012 events. ~12-17 GB at 40 Hz for 52 days.
F9 Magnitude slope 0.38 and 3-C/Z gain 3.4-6x inflated by the excess-over-max-over-scan nonlinearity and survivorship; incoherent-power
   slope 10^(0.46-0.58 dM) more honest; report recall-by-magnitude CURVES, not one Mc. 1.9x measured/predicted is a lower bound (~2.5-3.5x).
F10 Families 1 and 15 have ~0 excess at every station and drag n_sta=5 recall to 8%; blocklist or rebuild.
F11 SNR weighting right in direction but small (+17% MAD ~0.2 mag); station count and template length are the bigger levers; weights from
   build-half stack S-SNR only; rolled control recomputed with weights.
F12 Longer templates raise total detections at 8 MAD (1.2 -> 11/template/day, 2 -> 12 s) and the rolled control (0.2 -> 3.0); choose L on
   quiet-day FA-matched recall.

## Answers (short)
A accept "open"; 2.73 = extreme-value bias of max over 41 correlated samples (Gaussian ~2.77), like-for-like as a detection statistic; report
recall by n_sta. B 3.6 MAD is an artifact; operating point ~8 MAD (quiet-day 0.6 FA/template/day); budgets: front tracking <=10-20% FA on
peak days; re-rupture/burst metrics on bursts of >=3-4 detections within ~1 h; tidal burst-level, <=5% FA; with 8-12 s templates, weights
and 6-10 stations expect 30-50% recall. C verified; limits: his templates carry 1,400-2,000 stacked events vs our 100-300; his 7-station
catalog (Mc ~1.75) beats our corrected 7-station prediction (~2.2) by ~1.5x in MAD. D F7; fine-tune deferred/closed. E F4-F6. F fix F2;
S-VI post-2023 at 8 MAD ~17% equal-weight / ~35% SNR-weighted — adequate for the front, marginal for burst metrics.

## Prescription
0 Correct the record (RESULTS E4, NOTES decisions, PLAN GOAL "~7 km/day, 80 km in 11 d"): FA columns -> "detections/template/day" + add
  quiet-day FA; strike 3.6-MAD inversion, 4-MAD conclusion, "front reproduced at 3.3", "tidal growth reproduced", "secondary fronts not
  reproduced".
1 Detection-budget table on the 2010 ring (GPU minutes): thresholds 5,6,7,8,10; per family-day total declustered detections, rolled-control
  detections, detections within ±2 s of a Bostock event of that family; quiet (<=50 Bostock events/day) vs ETS days; per template length.
  Acceptance: quiet-day FA at 8 MAD <=1/template/day, at 7 <=10; ETS-day precision at 8 MAD >=85%. Failure: quiet-day FA at 8 >=5 ->
  operating point 9-10 MAD; re-derive downstream.
2 Fix the prediction model: raw = excess + 2.73 vs T; recall via measured log-spread (p10/50/90 = 3.5/6.1/9.6 at fixed M); SNR-weighted
  option sqrt(sum e²)·g; drop/label 4-MAD column; per-family actual station set in-sample. Acceptance: reproduces ring raw median 5.9±0.5
  and frac>=8 0.21±0.05 at M1.8-2.1; Bostock 21-station Mc(8) in 1.6-1.9; 7-station Mc(8) reported vs his 2010 mode 1.75 (gap = our
  implementation-vs-Bostock deficit). Predict the actual 2012 set; redo region-year table at 8 MAD only.
3 SNR-weighted network, 2010 held-out: weights ∝ 10^(build-half stack S-SNR dB/20), capped 3x median; rolled control with weights;
  blocklist families 1, 15. Acceptance: held-out excess at 8 s >=3.3 (from 2.64), frac>=8 >=0.22 at M1.8-2.1, rolled unchanged, quiet-day
  FA not worse.
4 Thomas harness fixes + re-run decimation (zero download): (a) main front from each family's first BURST (>=4 events within 10-min gaps)
  after 08-30, Theil-Sen over families with first burst inside the propagation window; centroid + leading-edge speeds with windows;
  reference 7.3 km/day, accept 5-10. (b) Secondary detector: front position from the interpolated 6-h leading-edge track; 3-h sliding
  window, 30-min step; event-level times; >=8 events from >=4 families spanning >=5 km; Theil-Sen speed 1-10 km/h with bootstrap CI
  excluding 0; |Spearman| >=0.6; merge overlaps; "behind the front" diagnostic not gate; nulls = circular per-family shift + within-day
  shuffle. Validate vs Thomas's four (dates/positions from PMC text; family 141 hosted six on 12-16 Sept). Acceptance >=3 of 4 with null
  <=1/episode. (c) Tides: drop Schuster p; one phase per burst start; time-shift null (1-30 d, >=200 draws); placeholder labelled; no
  claim without a Cascadia interface tidal-stress model. (d) Re-rupture = re-activation by a burst after >=12 h quiescence. (e) Re-run
  Mc 1.7-2.1 -> completeness requirement per metric.
5 DOWNLOAD RULING — Sept-2012 southern VI: APPROVED, conditional on steps 1-4 (one session) and step 2 predicting >=20% recall of Bostock's
  2012 events at 8 MAD on the actual 2012 set. Station set: every 3-C station within 60 km of the family centroid operating
  2012-08-15..10-05 (ring 6 + SOKB, BPCB, PA01-03, NLLB, LRIV, B001, B003-B007, B009-B011, B926), 40 Hz, ~12-17 GB. Test = from-scratch
  SNR-weighted network matcher vs his 28,149 events on the same (superset) geometry: templates for ALL stations from the 2012 build half
  (alternate 2-h blocks, 60-s guard), SeisBench stacks only as q-control; 8 and 12 s; 8 MAD. Pre-registered pass: recall of his test-half
  events within ±2 s >=25% overall, >=60% at M>=2.0, >=35% at M1.8-2.0, >=20% at M1.6-1.8; rolled control <=5% of our ETS-day detections;
  pre-SSE quiet days (08-15..08-29) <=1 detection/template/day; extras >=90% on days with >=100 Bostock events with median relative
  amplitude below matched events; Thomas metrics on our catalog vs his with the same estimators: front speed ±20%, re-rupture fraction
  ±0.1, secondary-front count ±1 (if 4b validated). Failure signatures: recall <15% with chance <=5% -> sensitivity short (compare
  per-station registration to prediction); chance >15% -> FA model wrong, raise threshold, add EQ/glitch vetoes; recall fine but metrics
  diverge -> family-dependent completeness.
6 T3 discovery on 2010 (on disk) once step 5's run is launched (parallel GPU); NorCal has independent truth on disk (Shelly MTJ). T4
  chaining after T3 passes. DL FINE-TUNING: deferred indefinitely (closed by F7).

## What would change Merlin's mind
Quiet-day FA at 8 MAD >=5/template/day -> operating point 9-10 MAD, 2012 replication likely fails 25%. Corrected model failing the
in-sample frac>=8 check -> raw/excess diagnosis wrong. Fixed secondary detector recovering >=3 of Thomas's four AND losing them at Mc 1.8
-> secondary fronts need >=35% of Bostock density; S-VI post-2023 cannot deliver -> replan the science before the fleet.
