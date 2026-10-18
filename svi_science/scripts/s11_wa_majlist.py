#!/usr/bin/env python
"""s11_wa_majlist -- PLAN §5: the Q3 major list for a region run, fixed by COUNTING headline cells only (no dv/v value is
read): every major == True unit with >= 20 headline cells (n >= 20, cc3c >= 0.35, causality > 1.5, series >= 3 cells).
Column 'cotimed_svi': unit window within +-30 d of one of the SVI 15 majors. Writes C.TB/maj_list.csv (read by s11_common).
    S11_REGION=WA_NW python s11_wa_majlist.py
"""
from __future__ import annotations
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I

SVI_MAJ15 = ["2003-02", "2004-07", "2005-09", "2008-05", "2010-08b1", "2011-08", "2012-08", "2013-09", "2015-12", "2018-06",
             "2019-08", "2019-11", "2020-11", "2023-04", "2024-08"]
assert C.REGION != "SVI"
cells = I.prepare("raw", C.HEAD, "", nmin=20, ccmin=0.35)
n = cells.groupby("episode").size()
U = C.episodes().reset_index()
M = U[U.cls == "major"].copy()
M["n_cells"] = M.episode.map(n).fillna(0).astype(int)
SV = pd.read_csv(f"{C.ROOT}/svi_science/results/s10_dvv/multi/tables/episode_units.csv").set_index("episode").loc[SVI_MAJ15]
def cot(r):
    for e, s in SV.iterrows():
        if pd.Timestamp(r.w0) <= pd.Timestamp(s.end) + pd.Timedelta(days=30) and pd.Timestamp(r.w1) >= pd.Timestamp(s.start) - pd.Timedelta(days=30):
            return e
    return ""
M["cotimed_svi"] = M.apply(cot, axis=1)
L = M[M.n_cells >= 20][["episode", "start", "end", "n_cells", "cotimed_svi"]]
L.to_csv(f"{C.TB}/maj_list.csv", index=False)
M[["episode", "start", "end", "n_cells", "cotimed_svi"]].to_csv(f"{C.TB}/maj_all_counts.csv", index=False)
print(f"{C.REGION}: {len(L)} of {len(M)} majors with >= 20 headline cells ({len(cells)} cells); co-timed with SVI majors: "
      f"{int((L.cotimed_svi != '').sum())}")
open(f"{C.OUT}/MAJLIST_DONE", "w").write(C.pt_now() + "\n")
