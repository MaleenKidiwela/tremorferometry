#!/usr/bin/env python
"""Cap the reused old-pipeline detections (32 NLLB families) to top-100/family-day by cc,
to match the new densify's per-(family,day) cap. Memory-safe chunked pass over ~107M rows.
Stacks are rebuilt from raw waveforms downstream, so resampler stays consistent."""
import pandas as pd

SRC = 'data/mf_nllb_old32_raw.csv'
OUT = 'data/mf_nllb_old32_capped.csv'
TOPN = 100

parts = []
nrows = 0
for chunk in pd.read_csv(SRC, usecols=['template', 'time', 'cc', 'station'],
                         dtype={'cc': 'float32'}, chunksize=15_000_000):
    nrows += len(chunk)
    chunk['day'] = chunk['time'].str.slice(0, 10)
    chunk = chunk.sort_values('cc', ascending=False).groupby(['template', 'day'], sort=False).head(TOPN)
    parts.append(chunk)
    print(f'  read {nrows:,} rows; chunk kept {len(chunk):,}', flush=True)

cap = pd.concat(parts, ignore_index=True)
# final re-cap (fix any (template,day) split across chunk boundaries)
cap = cap.sort_values('cc', ascending=False).groupby(['template', 'day'], sort=False).head(TOPN)
cap = cap.drop(columns='day')
cap.to_csv(OUT, index=False)
print(f'[done] {nrows:,} raw -> {len(cap):,} capped rows, {cap.template.nunique()} families -> {OUT}', flush=True)
