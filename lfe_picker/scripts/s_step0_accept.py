#!/usr/bin/env python
"""P3D STEP 0 acceptance: SSE-window delivery and per-station-day coverage.

Reads the downloader's manifest READ-ONLY.  Accept: every SSE station-day has coverage_frac >= 0.9
or is MASKED (written to results/p3d/masked_station_days.csv and excluded from every later step).
FAIL (report before continuing): SYMB delivery < 80 % of the SSE days it could have delivered.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

COV_MIN = 0.9


def main():
    M = pd.read_csv(f"{S.P3D}/download_manifest.csv", low_memory=False)
    ST = S.station_table()
    W = S.sse_windows()
    days = {}
    for r in W.itertuples():
        for t in pd.date_range(r.start, r.end, freq="D"):
            days[t.strftime("%Y-%m-%d")] = r.id
    M = M[M.date.isin(days)].copy()
    M["episode"] = M.date.map(days)
    M["cov"] = pd.to_numeric(M["coverage_frac"], errors="coerce")
    rows, masked = [], []
    for sta, g in M.groupby("sta"):
        if sta not in ST.index:
            continue
        inep = g[g.status != "no_epoch"]
        data = inep[inep.status == "data"]
        bad = data[~(data["cov"] >= COV_MIN)]
        for r in bad.itertuples():
            masked.append(dict(sta=sta, date=r.date, episode=r.episode, coverage_frac=r.cov,
                               reason="coverage_frac < 0.9"))
        rows.append(dict(sta=sta, site=ST.site[sta], role=ST.role[sta],
                         sse_days_in_epoch=len(inep), delivered=len(data),
                         frac_delivered=len(data) / max(len(inep), 1),
                         cov_ge_0p9=int((data["cov"] >= COV_MIN).sum()),
                         masked=len(bad),
                         frac_usable=float((data["cov"] >= COV_MIN).sum() / max(len(inep), 1)),
                         cov_median=float(data["cov"].median()) if len(data) else np.nan))
    D = pd.DataFrame(rows).sort_values("frac_usable", ascending=False)
    D.to_csv(f"{S.P3D}/step0_delivery.csv", index=False)
    pd.DataFrame(masked).to_csv(f"{S.P3D}/masked_station_days.csv", index=False)
    print(D.to_string(index=False))
    sym = D[D.sta == "SYMB"]
    ok = bool(len(sym)) and float(sym.frac_delivered.iloc[0]) >= 0.80
    v = dict(symb_frac_delivered=float(sym.frac_delivered.iloc[0]) if len(sym) else np.nan,
             symb_gate="PASS" if ok else "FAIL",
             n_masked_station_days=len(masked),
             masked_rule=f"coverage_frac < {COV_MIN} -> excluded from every later step")
    json.dump(v, open(f"{S.P3D}/step0_accept.json", "w"), indent=1, default=str)
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
