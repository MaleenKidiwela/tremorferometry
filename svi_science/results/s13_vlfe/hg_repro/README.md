# hg_repro — Hutchison & Ghosh (2019) VLFE matched filter, reproduced literally, then our controls

Rules and paper-step table: `PLAN.md` (written before any number). Log: `NOTES.md`. Answers: `RESULTS.md`.
Paper: `notes/papers/Hutchison_Ghosh_2018_JGR_VLFE.pdf`; Table 1 in `../tables/hg2018_table1.csv`.

Scripts (`scripts/`):
- `hg_digitize_fig6.py`, `hg_digitize_fig7.py`: known-answer numbers digitised from the paper figures (approximate).
- `hg_fetch.py`: missing 2011–2013 data for the paper's 9 stations (TA LH; SP2 BH decimated to 1 sps), log `logs/fetch.log`.
- `hg_prep.py`: 1-sps displacement day files (v3 recipe of `../scripts/scan_prep_v3.py`, LH allowed) in `data/disp1sps/`.
- `hg_cont.py`: continuous 0.02–0.05 Hz arrays per station (`data/mf/cont_*.npy`) + 90-s window norms (`den_*.npy`).
- `hg_refs.py`: far-reference LHZ arrays (`data/mf/ref_*.npy`) from ../data/raw_ref2011, raw_ref2012 (read-only).
- `hg_mf.py`: templates, summed-CC matched filter (real, noise templates, station-day-shift and polarity-scrambled controls).
- `hg_analyze.py`: thresholds (literal 3rd-highest, 6×MAD, rate-calibrated), catalogues, known-answer check K1–K4.
- `hg_cmt.py`: whole-space grid-search CMT (deviatoric, VR, CLVD, Kagan). `hg_stackinv.py`: stack-and-invert + selection-bias controls.
- `hg_controls.py`: A1/A2 controls, far/near flag, clustering, survival.
- `hg_queue.sh`, `hg_queue2.sh`: the two chained queues (markers below).

Data: `data/raw/` (fetched day files), `data/resp/` (StationXML), `data/disp1sps/`, `data/mf/` (CC series/peaks).
Status: `ls svi_science/results/s13_vlfe/hg_repro/*DONE` — PREP_DONE CONT_DONE MF_REAL_DONE MF_NOISE_DONE ANALYZE_DONE
(queue 1, `logs/queue.log`), GF_DONE CONTROLS_DONE STACKINV_DONE (queue 2, `logs/queue2.log`).
Restart a stage: each script skips finished outputs; re-run `nohup bash scripts/hg_queue.sh "" "" > logs/queue.log 2>&1 &`
after PREP_DONE exists (the PID arguments are only used before PREP_DONE), then `nohup bash scripts/hg_queue2.sh "" ...`.

R2 (reverse-engineered reproduction, PLAN Amendments R2/R2a/R2b/R2c; results in RESULTS.md "Reverse-engineered
reproduction"): `hg_r2.py` (variants, chosen-setting MF + controls + noise templates), `hg_r2_calib.py` (thresholds fitted
to Table S1 2011-ETS counts, automatic variant choice), `hg_r2_heldout.py`, `hg_r2_match.py` (detection match),
`hg_r2_stage2.py`, `hg_r2_s1.py` (tests on Table S1's own detections), `hg_r2b.py` (self-detection, stacks, overlap,
grid), `hg_r2c.py` (moment budget); queues `r2_queue1.sh`, `r2_queue2.sh`, `r2b_queue.sh`. Table S1:
`tables/hg2019_tableS1.csv` (from notes/papers/HG2019_SI/). Markers R2_DONE R2B_DONE R2C_DONE.

R3 (individual between-ETS candidates, PLAN Amendment R3; RESULTS.md "Individual between-ETS candidates"): `hg_r3.py` (rows / events YEAR with SCAN_YEAR=YEAR / aux / agg; R3_TAG=ets = post-hoc ETS positive control), `hg_r3_summary.py` (pass rules, classes, tables r3_events/r3_pass_counts/r3_aggregate, figures r3_*); queues `r3_queue.sh`, `r3_events.sh`, `r3_ets_queue.sh`; data/r3/. Markers R3_EVENTS_DONE R3_AUX_DONE R3_AGG_DONE R3_ETS_DONE.
