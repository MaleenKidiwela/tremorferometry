#!/usr/bin/env python
"""Incrementally splice one year's daily stacks into an existing long-window npz WITHOUT rebuilding
the whole record. Usage: merge_year_into_npz.py <full_npz> <year_only_npz> <YEAR>
Drops the full npz's rows for YEAR (the old truncated ones), appends all rows from the year-only npz."""
import sys, numpy as np

full_p, year_p, YEAR = sys.argv[1], sys.argv[2], sys.argv[3]
F = dict(np.load(full_p, allow_pickle=True))
Y = np.load(year_p, allow_pickle=True)
keep = np.array([not str(d).startswith(YEAR) for d in F["dates"]])
n_old = int((~keep).sum())
out = {}
for k in ("stacks", "patches", "dates", "n_det"):
    out[k] = np.concatenate([F[k][keep], Y[k]], axis=0)
out["t"] = F["t"]; out["fs"] = F["fs"]
np.savez(full_p, **out)
print(f"{full_p}: dropped {n_old} old {YEAR} rows, added {len(Y['dates'])} -> total {len(out['dates'])}")
