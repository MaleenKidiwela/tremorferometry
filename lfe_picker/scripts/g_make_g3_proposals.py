#!/usr/bin/env python
"""Gate-3 validation sets on a LONG span (the 2012 truth window is 52 d, on which a circular shift of
>= 30 d does not exist, so Gate 3 is structurally UNTESTABLE there).

  g3_pos            60 families of the lfe_extension catalog v1 (cat_core, 2014-01..2026-09), members
                    = the catalog's own detections, cell = the Bostock family position
  g3_neg_random     the same families / N, members drawn uniformly over the same span
  g3_neg_activity   the same families / N, members drawn only on days with PNSN tremor ANYWHERE within
                    300 km -- Merlin F13's blind spot: an activity-covarying false-alarm process
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402

N_FAM, N_MEM, SEED = 60, 300, 20260917
T0, T1 = np.datetime64("2014-01-01"), np.datetime64("2026-09-15")


def main():
    rng = np.random.default_rng(SEED)
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv").set_index("bkey")
    C = pd.read_parquet(f"{G.EXT}/results/trackD/lfe_extension_catalog_v1.parquet",
                        columns=["family_id", "t", "cat_core"])
    C = C[C.cat_core]
    C = C.rename(columns={"family_id": "fam"})
    cnt = C.groupby("fam").size()
    fams = [int(b) for b in cnt[cnt >= N_MEM].index if int(b) in F.index]
    rng.shuffle(fams)
    fams = sorted(fams[:N_FAM])
    days = np.arange(T0, T1 + np.timedelta64(1, "D"))
    Tr = pd.read_csv(f"{G.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    Tr["day"] = pd.to_datetime(Tr.time).values.astype("datetime64[D]")
    sets = {"g3_pos": [], "g3_neg_random": [], "g3_neg_activity": []}
    for b in fams:
        fr = F.loc[b]
        t = np.sort(rng.choice(C.t[C.fam == b].to_numpy(float), size=N_MEM, replace=False))
        base = dict(lat=float(fr.lat), lon=float(fr.lon), dep=float(fr.dep),
                    offset_search_s=0.0, truth_bkey=int(b), drop_bkey=int(b))
        sets["g3_pos"] += [dict(pid=f"g3pos_f{b}", t=float(x), kind="g3_pos", **base) for x in t]
        d = rng.choice(days, size=N_MEM)
        tt = (d.astype("datetime64[s]").astype(np.int64) + rng.integers(0, 86400, N_MEM))
        sets["g3_neg_random"] += [dict(pid=f"g3rnd_f{b}", t=float(x), kind="g3_neg_random", **base)
                                  for x in tt]
        dd = G.epi_km(Tr.lat.values, Tr.lon.values, fr.lat, fr.lon)
        act = np.unique(Tr.day.values[dd <= 300.0])
        act = act[(act >= T0) & (act <= T1)]
        d = rng.choice(act, size=N_MEM)
        tt = (d.astype("datetime64[s]").astype(np.int64) + rng.integers(0, 86400, N_MEM))
        sets["g3_neg_activity"] += [dict(pid=f"g3act_f{b}", t=float(x), kind="g3_neg_activity", **base)
                                    for x in tt]
    for k, v in sets.items():
        D = pd.DataFrame(v)
        D.to_parquet(f"{G.GATES}/proposals_{k}.parquet", index=False)
        print(k, D.pid.nunique(), len(D))


if __name__ == "__main__":
    main()
