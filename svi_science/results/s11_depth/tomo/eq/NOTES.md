# tomo/eq decision log

- **2026-10-13 ~9:35 am PT** — §6 written before the earthquake pick fetch (user OK).
- **2026-10-13 ~10:55 am PT** — §6 Amendment 1: zone-label fix; picker test on WA earthquake records.
- **2026-10-13 ~11:16 am PT** — §6 Amendment 2: volume split, synthetic emergent pulses, SVI picker calibration.
- **2026-10-13 ~11:31 am PT** — §6 Amendment 3: S/P width vs distance on LFE stacks; onset and ≤ 50 km κ refits.
  Ran ~11:32–11:42 am PT; pod restart ~12:20 pm PT killed the agent before the write-up. Outputs on disk.
- **2026-10-14 ~12:33 pm PT** — §6 Amendment 4 written before measuring: earthquake-record widths, synthetics with the
  measured broadening (W0–W3, W2s scrambled control), distance-shuffle scrambled control for the ≤ 40/50 km κ.
  Decision: use the 1–8 Hz κ band as headline (the band the picks come from), 2–8 Hz as a column; epicentral distance
  as in Amendment 3 (hypocentral ≤ 50 km leaves only epicentral ≲ 30 km at 40 km depth: column).
- **2026-10-14 12:37 pm PT** — launched `s11_eq_swidth2.py` (12:34, 4 workers for the cut, then 1 process) and
  `s11_eq_synth_width.py` (12:37, 4 workers; ~109k synthetic picks, est. ~30–40 min). Logs `logs/swidth2.log`,
  `logs/synthw.log`; PIDs in `swidth2.pid`, `synthw.pid`. Resume if killed: rerun the same command
  (`OMP_NUM_THREADS=1 nohup python svi_science/scripts/<script> --workers 4`); each script skips stages whose row
  table already exists (`swidth2_lfe_rows.csv`, `swidth2_eq_rows.csv`, `synthw_rows.csv`).
  Decisions while coding: the synthetic width trend uses only rows with finite widths; the per-row P width (W2) is
  used only where that LFE row had rP ≥ 1.2 (else the trend), the same rows the trend was fitted on.
- **2026-10-14 12:50 pm PT** — both runs done (~12:50). Rules applied as pre-written: Amendment 3 width rule not met (WA S − P
  broadening = SVI's); Amendment 4b W1 bias +0.002–0.015 (|bias| ≤ 0.05 → broadening does not explain); Amendment 3/4c
  distance limit: ≤ 50 km κ stays ≥ 1.96 (rises, not drops). Verdict: WA excess is not a broadening/peak-picking effect.
  Proposed next (not run): station-set split (boreholes + permanent vs temporary arrays); E3 S-noise correction on ≤ 50 km.
- **2026-10-14 ~2:50 pm PT** — §6 Amendment 5 written before computing: network split (permanent UW/US/PB vs temporary
  XU/TA/YW) with a class-label scrambled control, per-station residual maps vs the SVI-calibrated model (κ 1.773, CVM v1.7
  rays), B013 step check, and the S-pick scatter correction for near stations (headline instrumental-variable estimator with
  model tP as instrument; E3, Deming, SIMEX columns; synthetic closure; random-subset scrambled control).
  Decision: IV chosen as headline because it needs no σ_S (the co-located σ_S 0.04–0.11 s is from one site only).
- **2026-10-14 2:48 pm PT** — ran `s11_eq_netsplit.py --parts b,d,s,a --workers 3` (part c earlier, ~2:46); done ~2:54 pm.
  Log `logs/netsplit.log`, PID file `netsplit.pid`. Resume: rerun the same command (no stage caching; ~6 min).
  Found while running: the clipped E1 fit collapses (rms ≥ 1.5 s, κ 1.1–1.4) on the NW permanent class, on random half-networks
  (scrambled control median 1.39/1.47) and without the B013/B943/W030 site (40 % of kept NW rows). Reported with the
  start-free profile κ as the column that stays defined; the pre-registered rules were applied to E1 as written and the
  collapse is named. Co-located σ_S for NW from the Amendment 7 code = 0.131 s (σ_P 0.157); SW/SVI fallback 0.092/0.039.
  POST HOC (after seeing LOSO stations): drop tests for B013 site, B001, S020 (`tables/netsplit_posthoc_drop.csv`).
  Rules: (a) not a network effect; (c) B013 step ≤ 0.0002 (rule triggered only by a collapsed fit); (d) NW ≤ 40 km rise
  persists (+0.20 ± 0.06, fails σ/LOSO gates), ≤ 50 km and SW removed; synthetic scatter explains < 10 %; corrected
  all-station κ unchanged → excess not removed.
