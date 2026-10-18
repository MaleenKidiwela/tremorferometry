# newvlfe NOTES (times PT)

- 2026-10-16 ~10:00 am — read catalogue RESULTS §1–17, PLAN_CAT, Merlin reviews (design F1, review 2/3, scan NO-GO), S13 scan
  pilots 1–3. Found: the 2013/2018 scan-v2 runs exist for REAL data only (`data/scan_<Y>/runs/real_v2`, ndraw 0, shift 0);
  no v2 scrambled controls exist → must be run. scan_detect's shift control moves only the 7 pilot stations; with v2 the 5–7
  regional stations would stay aligned → copy `scripts/nv_detect.py` shifts all stations; references shifted by k_ctrl + 30.
- Looked at before PLAN (disclosed there): real positive candidate counts per z, real gate counts.
- ~10:05 am — PLAN.md written. ~10:04 am scan control queue launched (GPU): polB (4 polarity draws), shA3/7/11, per year.
- ~10:07 am — battery queue launched (1 worker while the GPU queue runs, then 2). Test on 300 rows: light z_chk equals the
  scan z (check of the light-stage arithmetic); ~0.07 s per cols call.
- ~10:15 am — nv_summary.py / nv_survivors.py tested on partial data (2013 real only; no control numbers yet). Validation of
  the added tests on 93 production tier events (prodq2013 + prodr2018): n_far ≥ 3 1.00, share_excess ≤ 0.25 0.98,
  vr_thrust ≥ 0.12 0.85, dense-ring plane wave 1.00 → kept; **z_loso ≥ 5 keeps 0.58 → fails the 80 % rule, removed from
  B1/B2 (column only)**, as PLAN §2 prescribes.
- HF envelopes on disk cover only 31 days of 2013 and 22 of 2018, so the battery's hf_loud column is mostly empty; the
  survivors stage computes the 2–8 Hz envelope ratio from the raw day files for every listed new candidate (real and
  control), cached in `data/hfenv/`.
- Waveform figure: CN.LZB (never used by the scan) is drawn grey and excluded from the common scale.
- Post queue launched (`scripts/nv_post_queue.sh`, waits for BATTERY_DONE → summary → survivors → POST_DONE).
- 10:43 am — POST_DONE. No pre-set tier passes for "new" (best: ETS B1 z ≥ 5, 45 vs 32.0 expected, 1.41×, p 0.017).
- ~10:50–11:05 am — coordinator follow-up: waveform figure scale now clips stations > 10× the median peak (broken LRIV 2018);
  `scripts/nv_wave_extra.py` (top new-ETS real rows + patch-LFE rows + polB controls of similar z), `scripts/nv_rank.py`
  (rank curves per year and control type, LRIV-broken flag). Findings: the new-ETS excess is all z 5–7 and all 2018
  (2013: 9 vs 8.7); LRIV was broken (3.5–9.8 µm square waves) at 92 % of the 2018 new rows; LRIV-ok rows 12 vs 12.8. Shift
  controls pass far/near 87 % vs real 13 % (incoherent noise population); polB is the closer control between ETS but
  is raised in ETS. RESULTS.md written.
- ~10:48 am — W2b 2018 re-run launched (scripts/nv_w2b_queue.sh: nv_detect NV_W2B=1 zeroes catalogue qc_exclusions 2018 station-component-months = UW.LRIV N/E/Z; polBW pass with 4 draws, variant 0 = realW; batteries realW/polBW). Marker W2B2018_DONE.
- ~10:58 am — W2b re-count (`scripts/nv_w2b_recount.py`): identical to before, because the scan noise rule had already zeroed
  LRIV in all 2018 months (daycfg w_on = 0 for N/E/Z in real_v2, polB, polBW). LRIV attribution in RESULTS retracted.
  2018 new-ETS B1: 36 vs 18.0 expected from polB (+18; z 5–7 +13.8). Western group: 17 real vs 1.5 expected; 6 with
  patch-LFE p ≤ 0.05 vs 0.5 expected; none a production tier event.
- ~11:00–11:10 am — point 5 (coordinator): nv_stackw.py (stacks W6/W17, CC vs bank, amplitude pattern, stacked scan location),
  nv_stackw_coh.py (member vs LOO stack vs 300 random ETS times), nv_stackw_ref.py (selection-matched reference = 208
  between-ETS scan candidates at the same nodes). Not RC4-shifted; weakly self-similar. Relay acceptance queued
  (nv_wbank.py; W2b QC on via P.QC = True; outputs newvlfe/data/wbank; K.log redirected to newvlfe/logs). Marker STACKW_DONE.
- ~11:20 am — STACKW_DONE: western relay test LOO 2/17 (rule ≥ 85 %), scrambled 0.11/yr, 8 extra tier events (5 known-patch, also production; 1 WA; 1 NW-only episode; 1 between, W-only, vr 0.117) → not accepted; no 2013/forward run.
- 2026-10-16 ~4:55–5:15 pm — NV3 (user approval ~5:05 pm relayed): PLAN Amendment NV3 written first. Inputs checked: scan-v2
  real runs and bp_<Y>v2 exist for 2011 (data/scan), 2016 and 2019; nothing to build. nv_battery.py now handles 2011 paths
  (bp2_ref_v2) and adds the Y9 emergence / T3 columns with each station on its own source day. Queues: nv3_scan.sh (GPU,
  W2b, waits for ≥ 22 GB free GPU memory), nv3_bat.sh (1 worker while the GPU queue runs, else 2), nv3_post.sh (2018
  batteries re-run with NV3 columns; old 2018 files moved to data/pre_nv3/; then nv3_summary.py and nv3_labels.py). Test
  outputs from 2016-only runs moved to data/nv3_test/. Label set capped at 120 real + 120 control (disclosed in PLAN). The
  M ≥ 4.5-within-90-min part of the eq-wave flag has a high base rate (~19 such events per day worldwide). Marker NV3_DONE.
- NV3 label criteria (fixed before looking at the first label figure): **wavefield** = similar-amplitude oscillation on most
  near and mid stations through the whole window, no emergence in the fit window. **eq-wave** = a transient wavetrain with
  moveout across the network and the far stations comparably excited, or a dispersive surface-wave train. **spike** =
  one or two stations dominate with a glitch/transient while the others are quiet. **candidate** = a pulse in the fit
  window, largest at the near (red) stations, decaying with distance, with little at the far stations and a quieter
  background before and after. Labelled from 4-figure sheets (data/nv3_sheets/), key not opened.
- labels done (~18:02 PT), 238 figures, key not yet opened; data/nv3_mylabels.txt
- ~7:30 pm PT — NV3 written up (RESULTS NV3 section). Labels: real cand 5/119 vs control 1/119 (p 0.21); no new source; 2016-06-13 07:41 is level 8 but 'known' (9 km).
