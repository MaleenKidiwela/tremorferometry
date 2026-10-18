---
name: pnsn-catalog-cut
description: The pnsn_tremor_2014-2026 catalog is download-cut at lat 47.5N; use cascadia_full master instead
metadata: 
  node_type: memory
  type: reference
  originSessionId: 5cecf374-37e4-4bcc-aea9-79b2b6cdc569
---

The PNSN tremor catalogs `catalogs/pnsn_tremor_2010.csv` and `pnsn_tremor_2014-2026.csv`
were fetched with `00b_fetch_pnsn_tremor.py --bbox 47.5 ...` so they are **clipped at
lat 47.50°N** — they contain ZERO tremor south of 47.5. This silently starves any station
near/south of that line (e.g. HDW at 47.649 looked like it had no southern families when
really the catalog just lacked the data).

**Use the master instead:** `catalogs/pnsn_tremor_cascadia_full.csv` — 748k rows,
lat 39.5–50.9°N, lon −128 to −121, 2010–2026 (full Cascadia, OR→BC). It was pulled earlier
for tremor-density plots.

When defining a discovery box, **center it on the station** and pick ±km (per user
preference) rather than an arbitrary geographic bbox; plot it first with
`scripts/plot_station_box.py`. For HDW, 43% of tremor within 60 km is *south* of the station.
Related: [[gpu-family-discovery]], [[nisqually-precursor-project]].
