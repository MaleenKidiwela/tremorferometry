# S14 NOTES (decision log, PT)

- 2026-10-14 1:40 pm PT — task received (coordinator): test waveform decorrelation as its own observable on LFE data on disk.
- 1:40–1:58 pm — inventory (PLAN §1). Per-episode stretch CC exists for SVI (8,946 cells) and WA-NW (18,428) in
  `s11_depth/.../episode_raw.parquet`, with n, ΣP, pre_rms_ratio, 2-parameter (shift) CC and c2-control CC. Day CC exists in
  round-7 front days (headline series) and S10 pilot daily. No direct-S CC, no successive-episode CC, no per-day noise.
- ~1:55 pm — user clarification via coordinator: the observable is the CC of the existing dv/v stretching measurement (the
  colour on our dv/v plots). Primary = existing CC columns; raw dv/v-and-CC side-by-side figure first, at the top of RESULTS.
- 1:59–2:20 pm — PLAN.md written before any CC-vs-time number. Expected |ΔCC| 0.005–0.04 (PLAN §2).
- 2:04 pm — step A done (`s14_raw.py`, RAW_DONE): figures A1–A3, raw "CC low & dv/v flat" list. First look: the episode CC
  series follow ln n closely (A3), and the co-located instruments share the CC drops (A2).
- 2:04 pm — step B pilot (`s14_store.py --pilot`): rebuilt episode stacks reproduce cc3c_25 exactly (median |Δcc| 0.00000) → PASS.
- 2:06 pm — step C pilot (`s14_tests.py --pilot`). Sanity: REG β(ln n) +0.20, β(ln noise) −0.19 (right signs), REG explains
  63 % of the z variance at 2–5 s. Scrambled-control centring rule (|median| < 0.2σ): C3 late−early FAILS (−1.3σ), C1 0.1–0.5σ.
- 2:07 pm — **Amendment 1 (after the pilot, before any full-run number):** (a) day CC also measured against a
  leave-day-out reference (`lo_cc_*`, used for C3; grand-stack CC kept as a column) to remove a possible self-inclusion
  bias; (b) permutation vectorised (same rule: shuffle within series inside ±45-day-of-year blocks); (c) C4 uses ≤ 300 draws.
- 2:08 pm — re-run pilot: (a) does not move the C3 scrambled offset (−1.45σ). Reading: the front-offset derangement keeps
  each episode's day order, so its "late − early" measures the generic "later days in an episode" pattern (≈ −0.07 in z
  after n and noise). The observed contrast (−0.05) equals that pattern → nothing front-locked. The scrambled control is
  therefore the right baseline (p values are taken about its median, as planned). Decision: machinery is sane
  (reproduction PASS, β signs right, offset understood); go to the full run. Logged as a deviation from the centring rule.
- 2:25 pm — full-run tests: first launch failed (RAW design row mismatch when covariates are missing; fixed), second was slow in C4 (python loop; vectorised, pilot M identical 0.2675). Relaunched `s14_queue_tests.sh` (log logs/queue_tests.log; markers TESTS_SVI_DONE, TESTS_WA_NW_DONE, QUEUE_TESTS_DONE). Stores: SVI 8,937 / WA 18,428 cells rebuilt, reproduction PASS (median |Δcc| 0).
- 2:30–2:35 pm — RESULT full run (TESTS_SVI_DONE, TESTS_WA_NW_DONE): no pre-registered comparison passes. C1 SVI +0.003 ± 0.005 z
  (p 0.54), WA −0.012 ± 0.004 (3–6 s p 0.001 but glitch-dropped p 0.66 and the c2 noise windows show the same → activity
  noise). C2 RF-like post − pre: SVI +0.053 ± 0.020 (p 0.013; flagged cells dropped −0.016, p 0.83), WA −0.029 ± 0.014 (p 0.10).
  C3 front-locked: none (p 0.81 / 0.99). C5 seasonal: SVI none; WA 0.0055 z p 0.004 (lead). C4 shared within station:
  M 0.28 / 0.24 vs scrambled p99 0.13 / 0.10 → passes "shared".
- 2:35–2:48 pm — POST HOC (labelled; after reading C4): station-episode CC term SD 0.068 / 0.066 z (27 % / 17 % of residual
  variance); dv/v station-episode term 0.00 % / 0.07 %; correlation 0.18 / 0.05; not regional (cross-site r median 0.12);
  cc2p (±50 ms shift allowed) keeps it (0.066) → not a whole-trace time shift. First reading "mostly timing" from the
  |shift| correlation was wrong (shifts are noisy: median 20 ms, p90 at the 50 ms rail) — corrected before writing.
  Largest: B013 2012-08…2013-05 (−0.34…−0.46 z, dv/v flat, B943 normal → B013 instrument, not in the glitch list);
  B009 2007-01…2009-03 (B009 instrument). RESULTS.md written ~3:05 pm.
- 2026-10-14 ~2:40–2:58 pm PT — INSTRUMENT CHECK (B013, B009; metadata + MUSTANG QC metrics only, no waveforms;
  evidence saved in `tables/instrument_meta/`).
  - FDSN StationXML: B013 EH* one epoch since 2007-01-06, B009 EH* one epoch since 2005-11-02 (HS-1-LT, same gain
    7.81e8, 100 Hz, same azimuths) → no documented sensor/response/orientation/gain/sample-rate change for either.
  - B013: EH1/EH2 daily RMS jumps on 2012-07-13 (to 4–60× B943; EH2 up to ~2,800 counts vs ~25–50), EHZ normal;
    data outage 2013-07-10…07-23; after it horizontals quiet but EH1 sits at ~0.2–0.3× B943 EH1 (was ~0.5× in 2008–12),
    lasting → likely the undocumented cause of the ballistic B013 horizontal-S step (round 3 dated it 2011-14 … 2015-18).
    DECISION: VETO 2012-07-13…2013-07-24 (horizontal fault; matches the S14 CC drop 2012-08 … B2013-05-24 and the return to
    normal at B2013-08-15). KEEP after 2013-07-24 (CC normal), flagged; S11 sensitivity = B013 series split at that date.
  - B009: no gaps/timing flags/availability change; but on 2007-01-08 EH1 RMS fell ~35 % relative to EHZ and EH2 overnight
    (winter EHZ/EH1 0.65–0.85 before and from 2009-11 on; 1.2–1.5 in 2007-01 … 2009-04; recovery May–Oct 2009).
    DECISION: VETO 2007-01-08…2009-10-01 (EH1 sensitivity fault; matches the S14 CC drop 2007-01 … B2009-03-09 and the
    normal B2009-10-10). End date uncertain within 2009-05…10.
  - Applied as flags (nothing deleted): `svi_science/instrument_vetoes.csv` (decision column veto / keep); S11
    `s11_invert.load_cells` drops cells with date_mid in a veto window (default ON; `S11_APPLY_VETO=0` reproduces the
    pre-veto tables, which are all S11/S14 tables on disk); `s14_common` adds the windows to the `glitch` flag.
  - IMPACT (`scripts/s14_veto_impact.py`, linear refits of existing cells; `tables/veto_impact_*.csv`): SVI 6 Q3 cells
    (0.1 % of weight): composite −0.000 → −0.002 ± 0.040 %, Q1 0.175 → 0.177 ± 0.052 %, Q2 +0.125 → +0.109 ± 0.111 %.
    WA-NW 16 cells (0.2 %): composite +0.012 → +0.008 ± 0.020 %, Q1 0.164 → 0.162 ± 0.072 %, Q2 −0.003 → −0.005 ± 0.089 %;
    majors 2012-08 / 2012-09 move ≤ 0.09σ. B013 split at 2013-07-24 (sensitivity): Q2 −0.005 → +0.052 ± 0.093 % (0.6σ),
    Q1 0.168, composite +0.009. No headline reading changes; gates/injections not rerun (would need the full queue).
