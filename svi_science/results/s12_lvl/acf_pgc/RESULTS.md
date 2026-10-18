# RESULTS — S12c acf_pgc (pilot written Tue 2026-10-13 ~11:30 am PT; full record + reading Wed 2026-10-14 ~2 pm PT)

Rules: `PLAN.md` (§1–7 and amendment A written before any time change was measured). Decision log: `NOTES.md`.
Δt = delay of the window at/after the LVL-bottom echo minus the delay of a window above the LVL (ms). A Vs increase
in the LVL makes Δt negative. "pct" = fraction of 2,000 scrambled-control draws (ETS dates moved to other years, same
season, > 60 d from any real midpoint) with Δt at least as negative. σ = robust scrambled spread (1.48 × MAD).

## What the data show (short)
1. **The ACF has stable late arrivals at the predicted LVL lags, at all three stations.** At 0.5–2 Hz, the radial
   autocorrelation (RR) shows a repeatable peak (same lag in 4 of 4 years) where an S wave reflected at the LVL bottom
   should arrive: PGC 23.8 s (predicted 23.4–24.6), VGZ 21.0 s (20.7–22.4), SNB 28.6 s (27.2–28.3). The whole ACF is very
   stable (year-to-year CC ~0.9), so one peak alone could be chance (15 % / 22 % / 6 % of same-width bands elsewhere in
   the trace hold a peak as strong). Three stations agreeing at their own predicted lags is the stronger argument
   (~0.2 % by chance if independent). Pass rule P1 met → the full record (2003–2026) was launched.
   *[Edit 2026-10-14, full record: the peak holds at PGC (23.8 s in 20/24 years) and VGZ (21.0 s in 24/24) but not at
   SNB (28.0 s in 8/19 years; its ACF is unstable over 2003–2021), so the three-station argument is now two stations.]*
2. **Pilot 2018–2021 (6 ETS), primary channel (PGC, 0.5–2 Hz, S–S, tremor-free hours): no post-ETS shift of that
   echo.** Days 0–15 vs days −30…−1: **Δt = −0.4 ms**, scrambled σ 1.0 ms (pct 0.50). A Gosselin-size change
   (+2–3 % Vs in 3.5 km) predicts **−40 to −58 ms** if the window is all echo; even if the echo carries only 20 % of the
   window, −8 to −12 ms. So the pilot already sets |Δt| ≲ 2–3 ms, about 20× below the full-echo target.
   The ±3-day curve stays within ±3 ms from −60 to +120 days (`figures/compare_rf_pilot2_CN.PGC.png`).
3. **1–4 Hz at PGC: a small shift of the right sign in all three phases, ~2σ each.** Days 0–15: S–S −5.7 ms (σ 2.7,
   pct 0.03), P–S (ZR) −22 ms (σ 12, pct 0.02; −28 ms at days 15–30, pct 0.00), P–P −3.3 ms (σ 2.6, pct 0.03). All-hours
   gives the same sign, smaller (−3.1, −13.6, −2.7). But the P–P control (two windows both above the LVL) also moves
   (+5.7 ms, pct 1.00, opposite sign), and the P–S curve wanders by ±40 ms at other times too. **Candidate only**;
   the full record (≈30 ETS, σ ~2× smaller) decides.
   *[Edit 2026-10-14, full record (23 ETS): **not confirmed**. S–S +11.8 ms (opposite sign; driven by the window above the
   LVL), P–S −9.0 ms (σ 10, pct 0.14), P–P −2.1 ms (σ 3.8, pct 0.25). The primary 0.5–2 Hz result of item 2 holds:
   +1.2 ms, σ 1.4. See '### Reading' at the end.]*
4. VGZ: nothing (0.5–2 Hz S–S −0.1 ms, σ 1.4). SNB: 0.5–2 Hz S–S −11 ms (pct 0.01), but its above-LVL control moves by
   −8 ms (pct 0.008) and the LVL-window amplitude rises 63 % → a change in the noise field, not a deep velocity change.

## Pilot table (days 0–15 vs −30…−1, tremor-free hours; `tables/summary_pilot2.csv`)

| station | band | phase | Δt (ms) | σ (ms) | pct | control (ms) | target (ms) |
|---|---|---|---|---|---|---|---|
| **PGC** | **0.5–2** | **S–S (primary)** | **−0.4** | **1.0** | **0.50** | +1.3 | −56 |
| PGC | 0.5–2 | P–S | +15.5 | 27 | 0.53 | +29 | −34 |
| PGC | 0.5–2 | P–P | +2.3 | 1.6 | 0.97 | +3.3 | −12 |
| PGC | 1–4 | S–S | −5.7 | 2.7 | 0.03 | +1.9 | −56 |
| PGC | 1–4 | P–S | −21.9 | 12 | 0.02 | +12 | −34 |
| PGC | 1–4 | P–P | −3.3 | 2.6 | 0.03 | +5.7 | −12 |
| VGZ | 0.5–2 | S–S | −0.1 | 1.4 | 0.19 | −0.7 | −58 |
| SNB | 0.5–2 | S–S | −11.2 | 1.4 | 0.01 | −8.3 | −52 |

## Not yet known
- Whether the 23.8 s PGC peak is mostly LVL-bottom S–S energy (it lies at the predicted lag at three stations, but the
  ACF rings everywhere). If it is only a small part of the window, the bound in item 2 loosens in proportion.
- Full record 2003–2026 (queue `scripts/acf_full_queue.sh`; ETA ~2 pm PT): same tests on ~30 ETS at all three
  stations, recovery-time fit, RF-vs-ACF figure (`figures/compare_rf_full_<STA>.png`).

## Next adaptation if the full record is also flat
Measure the LVL-**top** echo amplitude (the RF sees a contrast change, not a delay; the ACF analogue is the
amplitude of the 18.2 s S–S top echo, which this test does not use as primary), and go to 0.2–1 Hz, closer to the RF band.

## Echo amplitude / similarity stage (PLAN amendment B, queued ~11:40 am PT)
Gosselin measured a change in the boundary Vs contrast, not a delay. The pilot timing result rules out a bulk Vs change
of their size inside the layer (if the 23.8 s peak is the LVL echo), but a contrast change is still open. Queued after the
full run: relative echo strength and similarity with the long-term echo, LVL top and bottom, PGC/VGZ/SNB, same epoch and
scrambled control. Primary: PGC 0.5–2 Hz S–S top echo, tremor-free, days 0–15, two-sided. Expected size ~15–25 % of the
echo amplitude. The numbers are appended below automatically when `AMP_QUEUE_DONE` appears.

## Full record 2003–2026 — AUTO-GENERATED by scripts/acf06_report.py (Wed 2026-10-14 13:30 PT); reading to be added

### Timing, days 0–15 vs −30…−1 (Δt ms; σ robust scrambled; pct one-sided, Vs-increase direction)

| station | band | mask | phase | Δt | σ | pct | control | target |
|---|---|---|---|---|---|---|---|---|
| CN.PGC | b05_2 | none | SS | +1.4 | 1.3 | 0.852 | -2.0 | -56 |
| CN.PGC | b05_2 | none | PS | +25.8 | 22.5 | 0.865 | -7.5 | -34 |
| CN.PGC | b05_2 | none | PP | +0.3 | 1.7 | 0.682 | -0.5 | -12 |
| CN.PGC | b05_2 | tremorfree | SS | +1.2 | 1.4 | 0.820 | -1.7 | -56 |
| CN.PGC | b05_2 | tremorfree | PS | +17.3 | 22.8 | 0.751 | -4.8 | -34 |
| CN.PGC | b05_2 | tremorfree | PP | -0.1 | 1.7 | 0.606 | +0.7 | -12 |
| CN.PGC | b1_4 | none | SS | +9.2 | 3.7 | 0.994 | -8.2 | -56 |
| CN.PGC | b1_4 | none | PS | -4.7 | 9.8 | 0.247 | -3.0 | -34 |
| CN.PGC | b1_4 | none | PP | -2.7 | 3.6 | 0.176 | +7.6 | -12 |
| CN.PGC | b1_4 | tremorfree | SS | +11.8 | 3.8 | 0.999 | -8.8 | -56 |
| CN.PGC | b1_4 | tremorfree | PS | -9.0 | 10.0 | 0.141 | +0.9 | -34 |
| CN.PGC | b1_4 | tremorfree | PP | -2.1 | 3.8 | 0.249 | +9.3 | -12 |
| CN.SNB | b05_2 | none | SS | -21.3 | 15.1 | 0.097 | +29.8 | -52 |
| CN.SNB | b05_2 | none | PS | -4.4 | 19.1 | 0.459 | +17.5 | -31 |
| CN.SNB | b05_2 | none | PP | -0.9 | 3.1 | 0.344 | -5.6 | -11 |
| CN.SNB | b05_2 | tremorfree | SS | -24.2 | 15.5 | 0.085 | +40.8 | -52 |
| CN.SNB | b05_2 | tremorfree | PS | -16.4 | 18.5 | 0.229 | +27.7 | -31 |
| CN.SNB | b05_2 | tremorfree | PP | -2.4 | 3.2 | 0.190 | -4.8 | -11 |
| CN.SNB | b1_4 | none | SS | -24.4 | 8.1 | 0.014 | +36.0 | -52 |
| CN.SNB | b1_4 | none | PS | -3.7 | 12.0 | 0.360 | -6.1 | -31 |
| CN.SNB | b1_4 | none | PP | +25.2 | 8.9 | 0.987 | -3.4 | -11 |
| CN.SNB | b1_4 | tremorfree | SS | -26.3 | 7.9 | 0.010 | +43.2 | -52 |
| CN.SNB | b1_4 | tremorfree | PS | -10.8 | 11.9 | 0.186 | -10.4 | -31 |
| CN.SNB | b1_4 | tremorfree | PP | +30.3 | 9.3 | 0.993 | -1.0 | -11 |
| CN.VGZ | b05_2 | none | SS | +1.8 | 1.1 | 0.914 | -5.2 | -58 |
| CN.VGZ | b05_2 | none | PS | +20.8 | 12.1 | 0.941 | -23.3 | -36 |
| CN.VGZ | b05_2 | none | PP | -3.0 | 2.0 | 0.037 | +0.8 | -13 |
| CN.VGZ | b05_2 | tremorfree | SS | +0.8 | 1.1 | 0.809 | -2.6 | -58 |
| CN.VGZ | b05_2 | tremorfree | PS | +17.7 | 11.8 | 0.936 | -18.2 | -36 |
| CN.VGZ | b05_2 | tremorfree | PP | -1.7 | 1.8 | 0.110 | -0.2 | -13 |
| CN.VGZ | b1_4 | none | SS | +5.9 | 5.0 | 0.846 | -2.5 | -58 |
| CN.VGZ | b1_4 | none | PS | +11.3 | 4.9 | 0.971 | -3.6 | -36 |
| CN.VGZ | b1_4 | none | PP | +6.4 | 5.2 | 0.869 | +6.5 | -13 |
| CN.VGZ | b1_4 | tremorfree | SS | +0.4 | 4.4 | 0.462 | -5.1 | -58 |
| CN.VGZ | b1_4 | tremorfree | PS | -1.8 | 4.9 | 0.313 | -1.9 | -36 |
| CN.VGZ | b1_4 | tremorfree | PP | +4.4 | 4.9 | 0.778 | +2.9 | -13 |

Recovery fits recovery_fit_full_CN.PGC.csv (meaningful only where P3 passes): b05_2 SS: A +1.4 ms, τ 103 d; b1_4 SS: A +16.3 ms, τ 35 d; b1_4 PS: A -17.1 ms, τ 120 d

Recovery fits recovery_fit_full_CN.SNB.csv (meaningful only where P3 passes): b05_2 SS: A -82.5 ms, τ 40 d; b1_4 SS: A -50.4 ms, τ 48 d; b1_4 PS: A +1135.4 ms, τ 119 d

Recovery fits recovery_fit_full_CN.VGZ.csv (meaningful only where P3 passes): b05_2 SS: A -0.5 ms, τ 120 d; b1_4 SS: A -11.3 ms, τ 120 d; b1_4 PS: A -17.2 ms, τ 4 d

### Echo strength and similarity, days 0–15 vs −30…−1 (PLAN amendment B)
A = relative echo strength change (projection on the long-term echo, divided by the same for the window above
the LVL), in %; CC = change of CC with the long-term echo; σ robust scrambled; p2 = two-sided scrambled pct.

| station | band | mask | phase | A top % (σ, p2) | A bottom % (σ, p2) | CC top (σ, p2) | CC bottom (σ, p2) |
|---|---|---|---|---|---|---|---|
| CN.PGC | b05_2 | none | SS | -16.8 (16.2, 0.236) | +27.1 (10.3, 0.019) | -0.052 (0.042, 0.392) | -0.000 (0.031, 0.755) |
| CN.PGC | b05_2 | none | PS | +3.5 (13.2, 0.834) | -2.5 (11.4, 0.602) | -0.063 (0.021, 0.104) | -0.032 (0.029, 0.542) |
| CN.PGC | b05_2 | none | PP | -13.4 (9.1, 0.202) | -6.3 (4.5, 0.187) | -0.008 (0.018, 0.935) | -0.006 (0.004, 0.426) |
| CN.PGC | b05_2 | tremorfree | SS | -20.2 (18.5, 0.244) | +24.5 (11.1, 0.047) | -0.014 (0.044, 0.872) | +0.029 (0.033, 0.220) |
| CN.PGC | b05_2 | tremorfree | PS | +11.6 (13.0, 0.445) | -10.9 (12.2, 0.176) | -0.020 (0.022, 0.714) | -0.039 (0.029, 0.419) |
| CN.PGC | b05_2 | tremorfree | PP | -6.5 (9.2, 0.591) | -2.3 (5.2, 0.613) | -0.009 (0.018, 0.942) | -0.009 (0.005, 0.269) |
| CN.PGC | b1_4 | none | SS | +19.2 (9.7, 0.108) | +2.0 (13.5, 0.933) | +0.038 (0.035, 0.157) | -0.000 (0.042, 0.698) |
| CN.PGC | b1_4 | none | PS | +4.1 (19.3, 0.989) | +7.3 (14.6, 0.771) | +0.099 (0.061, 0.157) | -0.008 (0.072, 0.939) |
| CN.PGC | b1_4 | none | PP | -3.3 (10.7, 0.827) | -12.8 (13.5, 0.391) | -0.023 (0.043, 0.731) | -0.127 (0.051, 0.090) |
| CN.PGC | b1_4 | tremorfree | SS | +12.0 (10.2, 0.265) | +4.5 (13.9, 0.809) | +0.028 (0.037, 0.231) | -0.004 (0.042, 0.799) |
| CN.PGC | b1_4 | tremorfree | PS | +7.7 (20.2, 0.883) | +6.4 (15.7, 0.892) | +0.113 (0.065, 0.136) | -0.004 (0.078, 0.919) |
| CN.PGC | b1_4 | tremorfree | PP | +2.1 (11.0, 0.719) | -9.6 (13.3, 0.564) | -0.031 (0.044, 0.624) | -0.189 (0.052, 0.024) |
| CN.SNB | b05_2 | none | SS | -9.2 (27.0, 0.582) | +27.1 (46.9, 0.542) | -0.082 (0.061, 0.347) | -0.073 (0.111, 0.628) |
| CN.SNB | b05_2 | none | PS | -8.3 (22.2, 0.586) | +6.1 (17.9, 0.755) | -0.097 (0.034, 0.061) | -0.031 (0.030, 0.426) |
| CN.SNB | b05_2 | none | PP | +10.5 (18.3, 0.548) | +0.4 (6.7, 0.811) | -0.001 (0.019, 0.935) | -0.036 (0.045, 0.536) |
| CN.SNB | b05_2 | tremorfree | SS | -13.7 (27.3, 0.465) | +21.1 (47.5, 0.610) | -0.052 (0.063, 0.490) | -0.087 (0.112, 0.577) |
| CN.SNB | b05_2 | tremorfree | PS | -16.2 (21.6, 0.322) | +1.3 (17.5, 0.964) | -0.151 (0.033, 0.011) | -0.055 (0.029, 0.190) |
| CN.SNB | b05_2 | tremorfree | PP | +15.0 (18.8, 0.420) | +1.4 (7.2, 0.717) | -0.026 (0.020, 0.362) | -0.075 (0.046, 0.244) |
| CN.SNB | b1_4 | none | SS | -16.2 (26.7, 0.396) | +15.4 (36.1, 0.660) | -0.113 (0.066, 0.256) | -0.107 (0.104, 0.516) |
| CN.SNB | b1_4 | none | PS | -36.6 (15.6, 0.026) | +19.6 (26.0, 0.413) | -0.379 (0.085, 0.000) | -0.108 (0.107, 0.485) |
| CN.SNB | b1_4 | none | PP | +97.1 (30.2, 0.061) | +62.3 (59.2, 0.402) | +0.146 (0.125, 0.203) | +0.013 (0.203, 0.864) |
| CN.SNB | b1_4 | tremorfree | SS | -21.1 (26.9, 0.289) | +9.0 (35.6, 0.751) | -0.089 (0.067, 0.349) | -0.123 (0.105, 0.461) |
| CN.SNB | b1_4 | tremorfree | PS | -39.1 (15.2, 0.019) | +13.0 (25.8, 0.526) | -0.403 (0.082, 0.000) | -0.144 (0.106, 0.334) |
| CN.SNB | b1_4 | tremorfree | PP | +95.6 (30.1, 0.063) | +73.9 (59.8, 0.369) | +0.123 (0.126, 0.242) | -0.018 (0.208, 0.994) |
| CN.VGZ | b05_2 | none | SS | -1.7 (5.4, 0.652) | -3.8 (8.9, 0.822) | -0.016 (0.008, 0.189) | -0.062 (0.021, 0.081) |
| CN.VGZ | b05_2 | none | PS | +4.9 (4.9, 0.414) | +0.1 (5.2, 0.811) | -0.014 (0.007, 0.164) | -0.004 (0.005, 0.671) |
| CN.VGZ | b05_2 | none | PP | -3.9 (3.5, 0.112) | +0.4 (7.4, 0.866) | -0.001 (0.002, 0.682) | -0.002 (0.010, 0.917) |
| CN.VGZ | b05_2 | tremorfree | SS | -1.3 (5.3, 0.757) | +2.1 (9.0, 0.773) | -0.007 (0.007, 0.499) | +0.001 (0.022, 0.693) |
| CN.VGZ | b05_2 | tremorfree | PS | +2.3 (5.3, 0.825) | +1.1 (5.2, 0.938) | -0.009 (0.006, 0.248) | -0.002 (0.005, 0.856) |
| CN.VGZ | b05_2 | tremorfree | PP | +1.0 (3.6, 0.870) | -3.0 (6.7, 0.463) | +0.000 (0.002, 0.935) | +0.002 (0.010, 0.754) |
| CN.VGZ | b1_4 | none | SS | +10.3 (12.7, 0.498) | -21.4 (14.6, 0.145) | -0.080 (0.048, 0.222) | -0.311 (0.064, 0.010) |
| CN.VGZ | b1_4 | none | PS | +7.8 (6.1, 0.281) | -5.0 (7.6, 0.457) | -0.259 (0.028, 0.008) | -0.094 (0.037, 0.113) |
| CN.VGZ | b1_4 | none | PP | -8.0 (10.3, 0.720) | -9.7 (14.6, 0.726) | -0.140 (0.034, 0.038) | -0.111 (0.052, 0.171) |
| CN.VGZ | b1_4 | tremorfree | SS | +3.2 (12.2, 0.885) | -2.8 (12.5, 0.822) | -0.130 (0.052, 0.066) | -0.101 (0.054, 0.202) |
| CN.VGZ | b1_4 | tremorfree | PS | +1.6 (4.9, 0.784) | -10.2 (6.6, 0.151) | -0.158 (0.023, 0.021) | -0.099 (0.032, 0.090) |
| CN.VGZ | b1_4 | tremorfree | PP | -8.0 (10.3, 0.666) | -9.6 (13.9, 0.696) | -0.067 (0.029, 0.182) | -0.140 (0.046, 0.110) |

Figures: `figures/amp_full_<STA>_<band>_<mask>.png`, `figures/compare_rf_full_<STA>.png`, `figures/epoch_full_*.png`.

### Reading (written Wed 2026-10-14 ~2:00 pm PT; full record, 23 ETS at PGC, 22 at VGZ, 18 at SNB)

**What the data show**

1. **PGC timing: the LVL-bottom echoes do not move after the ETS.** Days 0–15 vs −30…−1, Δt in ms (σ = scrambled-control
   spread; target = Gosselin-size change, +2.5 % Vs in a 3.5 km layer):

   | PGC | S–S (RR, TT) | P–S (ZR±) | P–P (ZZ) |
   |---|---|---|---|
   | **0.5–2 Hz, tremor-free (primary = S–S)** | **+1.2 (σ 1.4)** | +17 (σ 23) | −0.1 (σ 1.7) |
   | 0.5–2 Hz, all hours | +1.4 (σ 1.3) | +26 (σ 23) | +0.3 (σ 1.7) |
   | 1–4 Hz, tremor-free | +11.8 (σ 3.8) | −9.0 (σ 10) | −2.1 (σ 3.8) |
   | 1–4 Hz, all hours | +9.2 (σ 3.7) | −4.7 (σ 10) | −2.7 (σ 3.6) |
   | target | −40 to −58 | −25 to −35 | −12 (0 if Vp unchanged) |

   The primary S–S shift is +1.2 ms, the wrong sign and inside the scrambled control (pct 0.82). A Gosselin-size change
   predicts −40 to −58 ms, 30–40× the scrambled σ. Even if the 23.8 s peak carried only 20 % of the window energy (target
   −8 to −12 ms), it would sit 6–9σ out. The ±3-day curve stays within ±3 ms from −60 to +120 d
   (`figures/epoch_full_CN.PGC_b05_2_tremorfree.png`). Later bins are flat too (15–30 d +0.9, 30–60 d +1.4, 60–120 d +0.3 ms).
   The stretch-scaled form gives +2.1 ms, i.e. the systematic between the two forms is ~1 ms. P–S at 0.5–2 Hz is too noisy
   to say anything (σ 23 ms ≈ the target). At 1–4 Hz none of the three phases moves in the Vs-increase direction beyond
   the control (see check 3).

2. **PGC echo strength and similarity: no change of the expected size can be seen, either way.** Primary (S–S **top** echo,
   0.5–2 Hz, tremor-free, days 0–15): **−20 % (σ 18.5 %, two-sided pct 0.24)**; all hours −17 % (σ 16 %). The expected
   size is ±15–25 %, i.e. about 1σ, so this test cannot detect it: only changes beyond about −50 % / +80 % (the
   scrambled 95 % band) would stand out.
   The similarity with the long-term echo does not change (top ΔCC −0.01, σ 0.04; bottom +0.03, σ 0.03).
   One cell stands out: the S–S **bottom** echo is **+25 % (σ 11 %, pct 0.047)** tremor-free and +27 % (pct 0.019) all
   hours. It is one of 12 cells with pct ≤ 0.05 among 144 in the 0–15 d table (about 7 expected by chance), and the
   post-hoc split below shows it is a LOW in the three weeks before the midpoint, not a rise after it.

3. **VGZ and SNB.** VGZ (22 ETS) 0.5–2 Hz S–S **+0.5 ms (σ 1.0, target −58)**; its echo strength is flat with a much
   smaller scrambled spread than at PGC (top −1 %, σ 5 %; bottom +2 %, σ 9 %), so at VGZ a 15–25 % change of the S–S top
   echo **is excluded at ~3–5σ**. At 1–4 Hz VGZ the LVL-window amplitude drops (×0.55–0.74, pct 0.001–0.04) and its
   similarity falls (ΔCC down to −0.31 all hours, −0.10 to −0.16 tremor-free): the noise field changes after ETS, smaller
   when tremor hours are removed. SNB: the ACF is not stable over the full record (yearly-stack CC 0.1–0.2 in the S–S
   windows; σ 8–15 ms). Its 1–4 Hz S–S −26 ms (pct 0.010) comes with the above-LVL control moving +43 ms (pct 0.99) and
   P–P moving +30 ms the other way: a change in the noise field, as in the pilot, not an LVL velocity change.

**Checks**

1. **Above-LVL control window (two windows both above the LVL top).** PGC primary: −1.7 ms, inside the scrambled 90 %
   band (pct 0.26). Passes. At 1–4 Hz the controls move: S–S −8.8 ms (pct 0.03), P–P +9.3 ms (pct 0.998).
2. **Tremor-free vs all hours.** PGC primary +1.2 vs +1.4 ms: same within 0.2 ms (σ 1.4) → the result does not depend on
   the tremor hours. Amplitude: −20 % vs −17 %, same sign.
3. **1–4 Hz pilot candidate (S–S −5.7, P–S −22, P–P −3.3 ms with 6 ETS): not confirmed with 23 ETS.** S–S is now
   **+11.8 ms** (opposite sign, pct 0.998), P–S **−9.0 ms (σ 10, pct 0.14)**, P–P **−2.1 ms (σ 3.8, pct 0.25)**.
   The S–S value is driven by the 13–17.5 s reference window (above the LVL), which moves 7–15 ms earlier around the
   midpoint while the LVL window hardly moves (`tables/curve_full_CN.PGC_b1_4_tremorfree.csv`, dt_pre vs dt_lvl); it
   starts 10–15 days before the midpoint. That is not the pattern of a velocity change at any one depth.
4. Multiple trials: in the timing table 2–5 of 36 cells per bin reach pct ≤ 0.05 in the Vs-increase direction (~2
   expected) and 4–8 in the other direction. At PGC the only one is 0.5–2 Hz P–P at 60–120 d (−2.4 ms, σ 1.4), far
   outside any RF window; the rest are at SNB and VGZ, mostly at 1–4 Hz where the controls and amplitudes move too.

**Post hoc (not pre-set): the RF recovery window** (`scripts/acf07_posthoc_rfwindow.py`, `tables/posthoc_rfwindow_CN.PGC_*.csv`,
same data, masks, estimators and scrambled control, 2,000 draws; only the day bins differ). RF at PGC: change from the
midpoint to ~+24 d, lowest at +6…+15 d, and ~half the contrast is a high 2–3 weeks before the midpoint.

| PGC 0.5–2 Hz tremor-free | S–S Δt (ms) | S–S top strength | S–S bottom strength |
|---|---|---|---|
| days 0–24 vs −30…−1 | +0.7 (σ 1.2) | −21 % (σ 18, pct 0.20) | +25 % (σ 10, pct 0.03) |
| days 6–15 vs −30…−1 | +2.0 (σ 1.6) | −18 % (σ 23) | +27 % (σ 14, pct 0.08) |
| days 0–24 vs −21…−1 (RF-style contrast) | +1.0 (σ 1.3) | −14 % (σ 18) | +27 % (σ 10, pct 0.03) |
| days 0–24 vs −60…−25 (reference without the pre-ETS weeks) | +0.3 (σ 1.3) | −18 % (σ 20) | **+2 % (σ 10, pct 0.94)** |
| days −21…−1 vs −60…−25 (pre-ETS weeks) | −0.8 (σ 1.2) | −3 % (σ 17) | **−20 % (σ 10, pct 0.03)** |

- The ACF delay is flat in the RF window too (|Δt| ≤ 2 ms against a −40 to −58 ms target).
- The S–S bottom-echo "+25 %" disappears when the reference avoids the three weeks before the midpoint: the bottom echo
  is ~20 % **weaker in the 3 weeks before** the midpoint and back to normal after it. That is the same time span as the RF
  pre-ETS high, but post hoc, one cell of many, and not a post-ETS change.
- 0.5–2 Hz P–P (ZZ) also moves in the pre-ETS weeks: Δt +3.7 ms (σ 1.6, pct 0.008) and top-echo strength +24 % (σ 8 %,
  pct 0.003), both masks. Its control moves the other way (−3.0 ms), so much of it is the above-LVL reference window.
  With ~430 post-hoc cells, one or two at pct ≤ 0.003 are expected by chance. A lead only; testing it needs a pre-set
  run (e.g. VGZ, or 0.2–1 Hz).
- 1–4 Hz: S–S +5 to +12 ms (wrong sign) in every RF window; P–S −1 to −13 ms (σ 9–14); P–P −0.1 to −4.5 ms (σ 3.5–5),
  none beyond the scrambled control in the Vs-increase direction.

**Scope of the no-change result.** What is bounded: the **two-way S (and P) travel time to the LVL-bottom reflector
(~47 km) relative to a reflector above the LVL top**, at **0.5–2 Hz**, in the **first 15 days after the ALBH midpoint
(also 0–24 d, post hoc)**, compared with **the 30 days before**, averaged over **23 ETS (2003–2026)** at **PGC**
(22 at VGZ). |Δt| ≲ 3 ms, i.e. a bulk Vs change in the LVL of more than ~0.1–0.2 % (or of ~2.5 % confined to ~0.2–0.3 km)
is excluded, **if** the 23.8 s peak is the LVL-bottom echo (plausible: it sits at the predicted lag at PGC and VGZ in
20/24 and 24/24 years, but the ACF rings everywhere). **Not bounded:** a change of the boundary **contrast** of
Gosselin's size at PGC (the amplitude test is ~1σ for 15–25 %; at VGZ it is excluded at 3–5σ), any change outside
0.5–2 Hz or faster than a few days, and the single-episode level (only the 23-episode average is tested). The RF and the
ACF need not disagree: a short-lived contrast change at the LVL top with little change in the layer's average Vs would
change RF amplitudes and leave the ACF delay flat.

**Pre-set pass rules**
| rule | result |
|---|---|
| P1 (echo visible; second station) | Met in the pilot. Full record: PGC RR 23.8 s in 20/24 years, VGZ 21.0 s in 24/24; **SNB weaker** (28.0 s in 8/19, yearly CC ~0.1–0.2) |
| P2 (scrambled σ ≤ 25 ms S–S) | **Pass**: σ 1.4 ms. (σ did not shrink as √N from the pilot's 1.0 ms: the 2003–2013 BH era and the gappy PGC 2012–2014 add spread) |
| P3 (primary Δt < 0, pct ≤ 0.05, control inside band) | **Fail**: +1.2 ms, pct 0.82, control inside → non-detection, bound as above |
| Amendment A (tremor-free vs all hours agree) | **Pass**: +1.2 vs +1.4 ms |
| Amendment B (S–S top strength, two-sided pct ≤ 0.05, same sign all hours) | **Fail**: −20 %, pct 0.24 (same sign all hours); test too insensitive at PGC for the expected 15–25 % |
| Recovery fit (§5) | Not meaningful (P3 fails); the fitted A = +1.4 ms is inside the control |

**Data and run checks (and what is broken)**
- Day counts (both bands, after the morning restart): no station-year has 0 days. PGC 143–366 d/yr (2010 181, 2012 285,
  2017 143, 2026 258 to Sep); VGZ 258–366 (2022 266, 2024 279); SNB 250–366 except 2021 (83 d, closed 2021-03).
  PGC 2013 and 2014 keep only 7–8 hours/day on average (other years 18–24): gappy data, not empty.
- **VGZ timing runs were stale**: b05_2 tremor-free and both b1_4 runs started 12:2x pm PT, before VGZ 2024–26 was
  recomputed (12:41–12:43 pm), so they used 20 ETS. Rerun at ~1:55 pm PT as tag `fullvgz2` (22 ETS;
  `tables/epoch_fullvgz2_CN.VGZ_*`): 0.5–2 Hz tremor-free S–S +0.5 ms (σ 1.0) [auto table: +0.8], 1–4 Hz S–S +1.6 / +6.3 ms
  (tremor-free / all), P–S −1.1 / +12.3 ms. Same reading. The auto table above still shows the 20-ETS values for those
  three VGZ rows. VGZ b05_2 all-hours and all amplitude runs used 22 ETS; PGC and SNB runs used complete data.
- The `P3_pass` column in `tables/summary_full.csv` (acf04_compare.py) checks only sign and pct, not the control window;
  e.g. SNB 1–4 Hz S–S is flagged True although its control is outside the band. Read P3 from this section.
- `figures/compare_rf_full_CN.PGC.png`: the RF curve's minimum sits at ~+25–35 d on that axis, but RF RESULTS §6 puts it
  at +9 d. The RF centre convention (noted 10-13) is still unresolved; do not read timing offsets from that figure.
- "All 30 ETS" in PLAN §4c: the ALBH list starts in 1994, so the 2003–2026 record holds 23 (PGC), 22 (VGZ), 18 (SNB).

## Lead test: "weeks before ETS" (PLAN §10, Amendment C; pre-registered Wed 2026-10-14 ~2:50 pm PT, run ~2:52–2:56 pm PT)
Script `scripts/acf08_lead_test.py` (acf07's code; only the bins and the control changed), queue `scripts/acf_lead_queue.sh`
(3 lanes), tables `tables/leadtest_<STA>_<band>_<mask>.csv`. Days −21…−1 vs −60…−25. Same-month scrambled control
(same calendar date ±15 d in another year, > 60 d from any real midpoint), 2,000 draws. Check: the pre-set 0–15 d values come
out as before (PGC S–S Δt +1.23 ms, top −20.2 %, bottom +24.5 %; VGZ S–S +0.47 ms).

**What the data show: VGZ does not repeat the PGC pattern.**

| 0.5–2 Hz, days −21…−1 vs −60…−25 | S–S bottom-echo strength (lead sign −) | P–P Δt, ms (lead sign +) | P–P top-echo strength (lead sign +) |
|---|---|---|---|
| PGC tremor-free (where the lead was found) | **−19.5 %** (σ 9.0, pct 0.006) | **+3.7** (σ 1.5, pct 0.002) | **+24 %** (σ 7.8, pct 0.001) |
| **VGZ tremor-free (C1, primary)** | **+5.3 %** (σ 8.4; pct for a drop 0.70) | +0.4 (σ 2.2; 0.49) | +7.5 % (σ 3.5; 0.09) |
| VGZ all hours | +1.3 % (σ 8.9; 0.60) | −1.5 (σ 2.0; 0.78) | +4.4 % (σ 3.3; 0.30) |
| VGZ 1–4 Hz tremor-free / all hours (secondary) | +7.2 % / +1.7 % | +9.2 / +10.5 (σ 4.7–5.4; 0.05 / 0.06) | −6.7 % / +5.4 % |
| SNB tremor-free / all hours (not usable, descriptive) | +41 % / +33 % (σ 34) | +4.3 / +2.4 (σ 2.8–2.9) | +3.5 % / +0.7 % |

(pct = one-sided in the lead's direction; σ = robust same-month spread.)

- **C1 fails.** VGZ's bottom echo is 5 % *stronger* in the three weeks before the midpoint, not 20 % weaker. A PGC-size drop would sit
  2.3σ out at VGZ (power ~75 % at the 5 % level; PLAN §10 said "about half", an underestimate), so this is a clear miss, though not proof of absence. The P–P parts of the lead are small and not significant
  at VGZ in the 0.5–2 Hz band (top-echo strength +7.5 %, pct 0.09). At 1–4 Hz the P–P delay goes the lead's way (+9 to +10 ms,
  pct 0.05–0.06). That band was not scored, and its controls moved in the full record as well.
- **With the tighter ±15-d control the PGC lead looks a little stronger** (bottom echo pct 0.006, was 0.014 with ±30 d), so the
  PGC lead is not a season effect. But PGC is where it was found, and it does not carry over to VGZ.
- **C2 (RF pre-window on GOBB, MGB, GNW, SP2) fails too:** H > 0 at 2 of 4, combined p 0.90; GNW goes the other way
  (`../rf_gosselin/RESULTS.md` §13).
- **Overall (rule fixed before running): the lead is not confirmed.** Neither the independent ACF station nor the independent RF
  stations show a change in the 3 weeks before the midpoint. The PGC ACF lead (bottom echo −20 %, P–P +3.7 ms) is still there with
  the same-month control, but it stays a single-station observation.
- Scope: 0.5–2 Hz (and 1–4 Hz, unscored) ACF at VGZ, 22 ETS 2003–2026; the 3 weeks before vs 4–8 weeks before the ALBH midpoint;
  episode average. 0.2–1 Hz was not tested: the stored ACFs hold only 0.5–2 and 1–4 Hz (new pass over the raw files ~2 h per
  station at 3 CPUs).
