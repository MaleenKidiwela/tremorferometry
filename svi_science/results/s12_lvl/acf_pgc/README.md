# S12c — post-ETS LVL change under PGC from ambient-noise autocorrelation (acf_pgc)

Started Tue 2026-10-13 ~10:45 am PT. Question, rules and pass criteria: `PLAN.md` (written before measuring; amendment A
dated). Decision log and inventory: `NOTES.md`. Results in plain words: `RESULTS.md`.

Scripts (`scripts/`):
- `acf00_predict.py` — predicted zero-offset PP / PS / SS lags of LVL top and bottom (PGC, VGZ, SNB) from the RF Ps
  delays + CVM v1.7; forward target delays → `tables/predicted_*`.
- `acf01_daily.py` — hourly whitened auto/cross-correlations (ZZ RR TT, ZR ZT RT; 0.5–2 and 1–4 Hz; 0–30 s) from the raw
  day-mseed already on disk → `data/<NET.STA>/<year>.npz` (~190 MB per station-year).
- `acf_common.py` — loading, tremor-hour masks (PNSN catalog 60 km ±30 min + 2–4.5 Hz RMS > DOY-matched 90th pct),
  stacks, window shift/CC measurement (×10 FFT upsampling), differential delay.
- `acf02_visibility.py` — long stacks, yearly stability, LVL-peak test (P1) → `tables/visibility_*`, `figures/visibility_*`.
- `acf03_epoch.py` — superposed epoch on ALBH ETS midpoints, scrambled control, day-by-day curve →
  `tables/epoch_*`, `tables/curve_*`, `figures/epoch_*`.
- `acf07_posthoc_rfwindow.py` — post-hoc RF-window bins (2026-10-14). `acf08_lead_test.py` + `acf_lead_queue.sh` — pre-registered lead test with a same-month control (PLAN §10).
- `acf_pilot_queue.sh`, `acf_full_queue.sh` — nohup queues; DONE markers in this folder.

`data_v0_daily/` = first PGC 2018 run in daily-mean format (superseded by the hourly format; kept).
