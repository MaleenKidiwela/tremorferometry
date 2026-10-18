# newvlfe — can a template-free detector find NEW VLFEs? (S13, 2026-10-16)

User (2026-10-16): "we are only detecting repeats and no new VLFEs." The catalogue (`../catalogue/`, read-only here) detects
with a template bank (2011 ∪ 2018 relay) that only finds events resembling templates near the C0/C1 patches. This study
re-tests the template-free scan (v2 network: 12–14 broadbands ≤ 160 km + 5 far LHZ references) as a detector, with the
catalogue's battery plus stricter far/near, plane-wave, single-station and mechanism tests, against scrambled controls,
separately for candidates near template patches ("known") and > 15 km away ("new"), in ETS and between ETS (2013, 2018).

Files: `PLAN.md` (pre-registered rules), `NOTES.md` (log), `RESULTS.md` (answer), `scripts/` (nv_detect.py = scan
detector copy with all-station shift; nv_battery.py; nv_summary.py; nv_survivors.py; queue scripts), `data/` (control scan
runs, battery parquet files, live time), `tables/`, `figures/`. Markers: `SCANCTRL_DONE`, `BATTERY_DONE`.
