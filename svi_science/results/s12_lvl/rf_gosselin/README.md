# S12 / rf_gosselin — receiver-function reproduction of Gosselin et al. (2020)

Goal: reproduce the post-ETS change in the low-velocity layer (LVL) seen with receiver functions by Gosselin et al.
2020 (Sci. Adv. 6:eaay5174) at CN PGC, SNB, VGZ; then extend it (2019-2026, full 1993-2026, more SVI + WA stations,
recovery curve) and write a calibrated target (`target.json`) for the LFE coda dv/v search.

- PLAN.md — recipe, data, scoring rule (written before results)
- NOTES.md — running log (bugs, decisions, numbers)
- RESULTS.md — results (after the queue)
- scripts/
  - rf00_events.py — ComCat M>=5.8 list, per-station 30-90° P arrivals (iasp91), channel epochs
  - rf01_download.py — P-150..P+150 s windows (CN from NRCAN/CNDC FDSN, UW/US from EarthScope), <=3 parallel, resumable
  - rf02_ets_gps.py — ETS midpoints from GPS east (paper recipe), ALBH for SVI, SEAT/SC04 for WA
  - rf03_rf.py — P-SV-SH receiver functions (rfpy conventions), 0.05-0.5 and 0.05-0.7 Hz
  - rf04_analysis.py — depth conversion, harmonic decomposition, CC, LVL dVs, ETS bins, scrambled control, Fig. 2 analog
  - rf05_recovery.py — recovery curve (sliding windows) and slip-front timing from PNSN tremor
  - rf06_align_check.py — quick check: ΔCC_SV aligned on GPS midpoint vs local tremor (2026-10-14)
  - rf07_samemonth.py — same-month scrambled control, seasonal regression, pre-ETS window H (2026-10-14, PLAN Amendment 2)
  - rf_queue.sh — the processing queue (DONE markers)
- data: `data/rf_teleseis/` (events, arrivals, raw/, rf/, gps/)
- restart after a reboot: block 14 / 14b of `scripts/resume_after_reboot.sh`

How to check progress:
`tail -n 2 logs/download_*.log; cat logs/queue.log; ls *_DONE ../../../../data/rf_teleseis/raw/*.DONE`
