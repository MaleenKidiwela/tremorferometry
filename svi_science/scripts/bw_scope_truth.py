#!/usr/bin/env python
"""bw_scope_truth -- SVI backward scope, step 3: Bostock (2015) 2003-2013 truth per year.

Maps every Bostock event to its family by exact (lat, lon, dep) match to lfe_catalog_repair e5r_families_2012
(fid = bkey, same rule as the lfe_picker referee freeze), keeps the 115 SVI families (svi_science s0 families.csv),
and tabulates per year: events, catalogue days (>= 1 event), days with >= 50 events ("ETS-active"), episodes
(runs of catalogue days separated by > 10 days), families active, first/last day.  Also per (year, family) counts.

Out: svi_science/results/backward_scope/truth_by_year.csv, truth_episodes.csv, truth_family_year.csv
"""
import pandas as pd, numpy as np

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
B = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
E = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv")
F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
key = lambda d: d.lat.round(2).astype(str) + "|" + d.lon.round(2).astype(str) + "|" + d.dep.round(2).astype(str)
m = dict(zip(key(E), E.bkey))
B["fid"] = key(B).map(m)
print("unmapped", B.fid.isna().sum(), "of", len(B))
B = B[B.fid.isin(F.fid)].copy()
B["t"] = pd.to_datetime(B[["year", "month", "day", "hour", "min"]].rename(columns={"min": "minute"})) + \
    pd.to_timedelta(B.sec, unit="s")
B["date"] = B.t.dt.floor("D")
days = B.groupby("date").size().rename("n").reset_index().sort_values("date")
gap = days.date.diff().dt.days.fillna(999)
days["ep"] = (gap > 10).cumsum()
EP = days.groupby("ep").agg(start=("date", "min"), end=("date", "max"), n_days=("date", "size"), n_events=("n", "sum")).reset_index(drop=True)
EP["year"] = EP.start.dt.year
EP.to_csv(f"{OUT}/truth_episodes.csv", index=False)
rows = []
for y, g in B.groupby("year"):
    d = days[days.date.dt.year == y]
    rows.append(dict(year=y, n_events=len(g), n_cat_days=len(d), n_days_ge50=int((d.n >= 50).sum()),
                     n_episodes=int(EP[EP.year == y].shape[0]), n_families_active=g.fid.nunique(),
                     n_families_ge20=int((g.groupby("fid").size() >= 20).sum()),
                     first=str(d.date.min().date()), last=str(d.date.max().date())))
T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/truth_by_year.csv", index=False)
FY = B.groupby(["fid", "year"]).size().unstack(fill_value=0)
FY.to_csv(f"{OUT}/truth_family_year.csv")
print(T.to_string()); print(EP.to_string())
