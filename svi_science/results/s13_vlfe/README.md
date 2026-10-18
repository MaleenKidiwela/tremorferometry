# S13 — VLFE-band (20–50 s) stacks on LFE detection times

Question: do our LFE times carry a time-locked very-low-frequency moment, and how does it vary per family, per
episode, between ETS, and with GNSS slip? Rules: `PLAN.md` (+ Amendments). Log: `NOTES.md`. Answers: `RESULTS.md`.

Scripts (`scripts/`):
- `s13_inventory.py`: continuous data on disk. Output: `tables/inventory_station_year.csv`, `tables/inventory_bb_near.csv`.
- `s13_green.py`: whole-space moment-tensor Green's function (near, intermediate and far field; free surface ×2).
- `s13_common.py`: paths, SVI rows and family loaders, origin shift, station coordinates.
- `s13_fetch.py`: resumable EarthScope day-file downloads, ≤ 3 parallel, into `data/raw/`.
- `s13_prep.py`: raw → 1-sps displacement N/E/Z day arrays, in `data/disp1sps/` (pilot) or `data/disp1sps_v2/` (scale-up).
- `s13_stack.py`: pilot stacker and fit. `s13_pilot_queue.sh` and `s13_pilot_figs.py` make the pilot figures and verdict.
- `s13_units.py`: per-episode and between-ETS row files. `s13_stack2.py`: scale-up stacker.
  `s13_scale_summary.py`: tests A1–A6. `s13_scale_queue.sh`: the scale-up queue.
- `s13_hg_crosscheck.py`: our LFEs around the Hutchison & Ghosh 2018 VLFE template times.

Data:
- `data/resp/`: StationXML responses.
- `data/raw/`: downloaded day files (2012 pilot gap: LZB, YOUB; PFB has no data at EarthScope).
- `data/stacks/*.pkl`: stacks.
- `data/rows/`: unit rows.

Check status: `ls svi_science/results/s13_vlfe/*DONE`; `tail svi_science/results/s13_vlfe/logs/scale_queue.log`.
Restart after a reboot: `bash scripts/resume_after_reboot.sh` (block 16).
