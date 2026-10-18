---
name: family-selection-rule
description: Family selection = coverage-balanced + top-10% highest-SNR repeaters added back; how to add families efficiently
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**Family selection for dv/v (user-set 2026-06-05): coverage-balanced selection, THEN add back the top ~10% highest-SNR families of the eligible pool.** Don't let coverage-balancing drop strong repeaters — high-SNR repeating families are the cleanest dv/v sources, so keep them all.

**Why:** `select_coverage_families.py` picks top-2 by SNR per az×dist cell (so it DOES prioritize SNR within a cell), but caps per cell for spatial spread — which drops a 3rd/4th strong family in a dense cell. The user wants those kept. Pure SNR-ranking alone is wrong too (families cluster in one azimuth → spatially biased cross-patch median). So: **coverage set ∪ top-10%-by-SNR of the pool** (pool = families with SNR≥floor AND dist≤rings-max).

**How (recipe):**
```
elig = pool[(snr>=floor) & (dist<=135)]
top10 = elig.nlargest(ceil(0.10*len(elig)), 'snr')   # top 10% by SNR
final = coverage_ids ∪ top10_ids
```
Mostly matters for SPARSE stations (B026: 35→45). Rich stations (B011/B013, 70+ fam) already swept up the high-SNR ones.

**EFFICIENCY (user feedback, learned the hard way):** to ADD families to an already-densified station, densify ONLY the new families (separate `--out-prefix mf_<sta>extra_`) and concat into `mf_<sta>_all` — **do NOT delete the base mf_<sta>_YYYY.csv and re-densify the whole set** (the densify-resume skips existing year files, so I deleted them to force the new set in, then wastefully re-densified all 45 at B026). Better: select coverage+top10% UPFRONT so the first densify already includes them (zero re-densify). See [[margin-station-status]], [[glitch-days-false-dvv-drop]].
