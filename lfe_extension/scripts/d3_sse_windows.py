#!/usr/bin/env python
"""Track D, D3 -- the SSE window list for the science QC (Merlin Step 8).

Definition (pre-registered, from the brief): a day is ACTIVE when the PNSN tremor catalog has >= 50
epicentres inside the family box (the Bostock 2012 family lat/lon range padded 0.2 deg -- the same box
Track C uses for `tremor_flag`); active days are clustered with gaps <= 20 days; a cluster is an SSE
when it spans >= 10 days.  The window is the cluster's first..last active day, padded by --pad days so
that the onset and the tail of the episode are inside it."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-per-day", type=int, default=50)
    ap.add_argument("--max-gap-d", type=int, default=20)
    ap.add_argument("--min-span-d", type=int, default=10)
    ap.add_argument("--pad", type=int, default=3)
    ap.add_argument("--t0", default="2014-01-01"); ap.add_argument("--t1", default="2026-09-30")
    ap.add_argument("--out", default=f"{EXT}/results/trackD/sse_windows.csv")
    a = ap.parse_args()
    F = pd.read_csv(f"{REP}/results/e5r_families_2012.csv")
    lat0, lat1 = F.lat.min() - 0.2, F.lat.max() + 0.2
    lon0, lon1 = F.lon.min() - 0.2, F.lon.max() + 0.2
    tr = pd.read_csv(PNSN); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    box = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1) & tr.time.notna()].copy()
    box["date"] = box.time.dt.normalize()
    c = box.groupby("date").size().rename("n").reset_index()
    c = c[(c.date >= a.t0) & (c.date <= a.t1)]
    act = c[c.n >= a.min_per_day].sort_values("date").reset_index(drop=True)
    rows = []
    if len(act):
        grp = (act.date.diff().dt.days.fillna(0) > a.max_gap_d).cumsum()
        for g, sub in act.groupby(grp):
            span = (sub.date.max() - sub.date.min()).days + 1
            if span < a.min_span_d:
                continue
            s = sub.date.min() - pd.Timedelta(days=a.pad)
            e = sub.date.max() + pd.Timedelta(days=a.pad)
            tot = int(c[(c.date >= s) & (c.date <= e)].n.sum())
            rows.append(dict(sse_id=f"SSE{sub.date.min():%Y%m}", t0=str(s.date()), t1=str(e.date()),
                             peak_date=str(sub.loc[sub.n.idxmax(), "date"].date()),
                             n_active_days=len(sub), span_days=span + 2 * a.pad,
                             n_tremor_in_box=tot, peak_per_day=int(sub.n.max())))
    W = pd.DataFrame(rows)
    # disambiguate two episodes that start in the same month
    if len(W) and W.sse_id.duplicated().any():
        W["sse_id"] = [f"{r.sse_id}{chr(97 + i)}" if n > 1 else r.sse_id
                       for (r, n, i) in zip(W.itertuples(), W.groupby("sse_id").sse_id.transform("size"),
                                            W.groupby("sse_id").cumcount())]
    W.to_csv(a.out, index=False)
    print(f"[box] lat {lat0:.2f}..{lat1:.2f}  lon {lon0:.2f}..{lon1:.2f}; "
          f"{len(box)} tremor detections in box, {len(act)} active days (>= {a.min_per_day}/day)")
    print(W.to_string(index=False))
    print(f"[write] {a.out}  ({len(W)} SSE windows)")


if __name__ == "__main__":
    main()
