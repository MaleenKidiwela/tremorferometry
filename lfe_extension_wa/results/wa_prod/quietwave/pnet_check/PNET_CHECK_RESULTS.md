# Are the PNET-selected detections real? Waveform check of P_own_WAVE_PNET ≥ 0.75 (core network), 2026-10-06

## Verdict in plain words

- **NW, truly quiet days (catalogue era 2009-08-06 to 2026-05-27, no PNSN tremor within 60 km): the PNET set is mostly real.**
  - The 678 detections (591 scorable) score **0.95 of the tremor-day level [90 % range 0.81–1.07]** on the
    outside-the-template-window test. Matching on network cc gives 0.98 [0.86–1.06].
  - That is a false share of about **5 % (0–19 %)**. The coda part alone gives 0.85 [0.67–1.01], or about 15 % false.
  - The control-based false shares (reversed 3.2 %, scrambled 1.5 %) are therefore **somewhat optimistic but not badly
    wrong**. A fair bracket for this set is ~5–15 % false, with an upper end of ~20–30 %.
- **The events PNET adds** (PNET ≥ 0.75, WAVE < 0.75; 237 on quiet days) are **about half real**: 0.70 [0.52–0.91],
  or 0.61 [0.47–0.79] after matching on network cc. The events both models pass (441) look fully real: 1.09 [0.94–1.21].
  So the false share of the PNET set comes almost entirely from the PNET-only part.
- **The events PNET rejects** (WAVE ≥ 0.75, PNET < 0.75; 694 on quiet days) are **mostly noise, but not all**:
  0.40 [0.31–0.51], or 0.34 [0.26–0.42] after matching. That is ~60–65 % noise. About a third, roughly 200 events,
  look like real LFEs that PNET drops.
- **The coda-only check matters for PNET.** PhaseNet sees ±20–40 s of context, so the P-wave window is not fully
  independent of PNET's selection. On quiet days, the PNET-selected sets score higher on the P-wave part than on the coda
  part (both: full 1.09 vs coda 0.95; PNET-only: 0.70 vs 0.56). The WAVE-only set does not show this gap (0.40 vs 0.39).
  This looks like a small selection leak through the P window, so the coda-only numbers are the safer reading for
  PNET sets.

## Key numbers

The score is the outside-window match to an independent tremor-day stack of the same family and station. It covers the
P wave [−8, 0) s and the coda [12, 24) s. "Same rule" divides by tremor-day detections passing the same rule
(1 = as real as those). "Matched" divides by the tremor-day detections both models pass, at the same network cc.
Random times are the zero reference. Ranges are 90 % day-block bootstrap.

| NW set | day class | n scored | same rule | matched | coda only (same rule) | random times |
|---|---|---|---|---|---|---|
| PNET ≥ 0.75 | truly quiet | 591 | **0.95 [0.81–1.07]** | 0.98 [0.86–1.06] | 0.85 [0.67–1.01] | 0.05 |
| both ≥ 0.75 | truly quiet | 384 | 1.09 [0.94–1.21] | 1.14 [1.02–1.23] | 0.95 [0.78–1.12] | 0.03 |
| PNET only (WAVE < 0.75) | truly quiet | 207 | 0.70 [0.52–0.91] | 0.61 [0.47–0.79] | 0.56 [0.26–0.82] | 0.07 |
| WAVE only (PNET < 0.75) | truly quiet | 599 | **0.40 [0.31–0.51]** | 0.34 [0.26–0.42] | 0.39 [0.17–0.59] | −0.01 |
| WAVE ≥ 0.75 (old headline) | truly quiet | 983 | 0.60 [0.47–0.70] | 0.65 [0.53–0.76] | 0.56 [0.39–0.69] | 0.01 |
| PNET ≥ 0.75 | unknown (no catalogue) | 1,479 | 1.27 [1.21–1.33] | 1.08 [1.02–1.14] | 1.21 | 0.00 |
| PNET only | unknown | 363 | 1.35 | 0.92 [0.83–0.99] | 1.35 | 0.02 |
| WAVE only | unknown | 557 | 0.98 [0.79–1.13] | 0.69 [0.58–0.80] | 0.99 | −0.03 |
| PNET only | tremor days | 6,706 | (= 1) | 0.83 [0.81–0.85] | — | 0.02 |
| WAVE only | tremor days | 4,498 | (= 1) | 0.79 [0.77–0.81] | — | 0.00 |

- **Unknown days:** all sets look real. The PNET-rejected events there are ~70 % real after matching.
- **NW tremor days:** both disagreement sets are mostly real (~0.8 of the agreed events at the same network cc). PNET's
  extra tremor-day events are as real as the ones it rejects.
- **Arm B, tremor days.** These are relative to the events both models pass (not matched on cc).
  - PNET ≥ 0.75: 0.92 [0.85–0.99].
  - PNET only: 0.78 [0.71–0.87] (217 events).
  - WAVE only: 0.64 [0.53–0.80] (98 events).
  - Random times: 0.00.
  - Both disagreement sets are mostly real. The PNET-only events look slightly more real than the PNET-rejected ones.
- **SW, tremor days: not resolvable as a share.** SW windows carry little outside-window signal. The raw mean match is
  0.009 [0.006–0.012] for PNET ≥ 0.75, against 0.000 at random times; NW and Arm B are at 0.04–0.05. The tremor-day
  yardstick is therefore too noisy for a ratio: random times give ratios of 0.1–1.0. What it does show: PNET-only SW
  windows (0.0097 [0.006–0.014]) are at least as LFE-like as the PNET-rejected ones (0.0063 [0.003–0.010]), and both
  are above random times.
- **Controls under PNET:**
  - The few scrambled controls PNET passes (9 on quiet days, 172 on tremor days) score 0.5–0.6 of real. This fits the
    earlier finding: a scrambled window that PNET passes usually carries one site's real LFE.
  - The reversed controls PNET passes (22 quiet) score ~0.2, and their window energy is huge. They are transients.

## Caveats

- **Tremor-day detections are treated as 100 % real.** On tremor days the controls give < 1 % false. This assumes
  quiet-day LFEs are as strong as tremor-day LFEs of the same family. Matching on network cc corrects part of that.
- **Reference choice.** The quiet-day WAVE ≥ 0.75 share is 0.60 here but was 0.43 in the first test. The reference
  stack differs: here both models must pass ≥ 0.85, two-fold on all tremor days; before it was WAVE ≥ 0.85 on 100 sampled
  even days. Treat about ±0.15 as the method's own spread. The ranking of the sets does not change.
- **Random-time offset.** On quiet days, random times of the PNET families score 0.05–0.07 rather than 0. This is small,
  and subtracting it would lower the PNET real share to ~0.90.

## Files and method

- Tables: `tables/real_share_summary.csv` (every set × rule × series × normaliser) and `tables/per_candidate.parquet`
  (one row per scored candidate). Raw tremor-day means: `tables/tremor_raw_outcc_means.csv`. Replication check:
  `tables/replication_check.json`.
- Replication: the recomputed network cc equals the stored `net_cc` for all 44,625 data and control rows
  (max |Δ| 5.8e-4).
- Scripts:
  - `scripts/pc_extract.py`: NW/SW/Arm B core candidates with both P columns, plus the v1 day class.
  - `scripts/pc_windows.py`: cuts the windows with `quietwave/scripts/qw_windows.do_day`, unchanged; 20 CPU processes,
    21 min.
  - `scripts/pc_analyze.py`: reference stacks and per-candidate scores.
  - `scripts/pc_summary.py`: real shares.
- Reference stacks: tremor-day data where both models pass (NW ≥ 0.85; SW and Arm B ≥ 0.5), two-fold by day parity.
  Tremor-day rows are scored against the stack from the other parity. Quiet and unknown rows are scored against the
  stack from all tremor days.
- The run used CPU only. No downloads, deletions or git, and the production queue markers were not touched.
  `work/` holds 8.4 GB of float16 windows. It can be regenerated and deleted once you OK it.
