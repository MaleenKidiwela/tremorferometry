# lfe_oregon — Thomas & Bostock (2015) central-Oregon reproduction

Paper: notes/1-s2.0-S0040195115003807-main.pdf (Tectonophysics 658:111). Summary: notes/CATALOGUE_SUMMARY.md §7.
Target: 18 families, >15,000 detections 2009-08-03..09-17, 13 located (4 free depths 31–34 km), front 5 km/day south.
Method (not yet built): Savard & Bostock cross-station detection on ALP0/ALS0/BWO0, network cross-correlation 8×MAD,
phase-weighted stacking, 10–15 iterations, ≥100 detections in 3 days, 10 % / 5 s dedupe. 1–10 Hz, 40 sps, 24 h files.

Data: scripts/fetch.py → raw/<NET.STA>/<date>.mseed (30 stations × 46 days), manifest results/fetch_manifest.csv,
marker results/FETCH_DONE. Station list notes/thomas2015_candidate_stations.csv.

Status 2026-09-28: 1–10 Hz cache built (cache40/, scripts/or_cache40.py). Cross-station detector = scripts/xs_detect.py
(port of Savard's MATLAB, notes/oregon_sources/CS-quake-detection). Chain = scripts/xs_chain_queue.sh (not yet run).
Recipe + silent points: METHOD_FROM_PAPER.md. Log: NOTES.md.
