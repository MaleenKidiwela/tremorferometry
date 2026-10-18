---
name: pandas-datetime64us-trap
description: "Recurring bug — pandas to_datetime / parquet timestamps are datetime64[us] on this pod, so .astype(\"int64\")/1e9 mis-scales times 1000x (empty joins, epoch ~1e6)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 03499166-a677-42d8-9397-e2cf5e28621f
  modified: 2026-09-10T19:42:30.856Z
---

On this pod pandas returns `datetime64[us]` (not ns) from `to_datetime` and from parquet (e.g. the Lin LFE parquet
`data/detections_lin_ets_2010_vi.parquet`, `lfe_features/data/lin_precise_b011.csv` joins). `.astype("int64") / 1e9`
therefore gives epoch ~1.27e6 instead of 1.27e9 and time joins come back EMPTY. Hit twice: lfe_features DECISIONS E23
(2026-06) and lfe_catalog_repair D3 (2026-09-10).

**How to apply:** always `pd.to_datetime(x).astype("datetime64[ns]")` (or use `.view("int64")` after `astype("datetime64[ns]")`,
or `(x - pd.Timestamp(0)) / pd.Timedelta("1s")`) before any epoch arithmetic; put this in every agent spec that touches
time columns. Related: [[orchestrate-with-opus-subagents]].
