#!/usr/bin/env python
"""SH_CHECKS CHECK 2b: day-level coincidence of LFE counts with PNSN tremor days (30 km circle; also 15/60 km)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_checks_pnsn import hav, RES, OUT, PN, LASTPN
from sh_v3_common import his_catalogue
F = pd.read_csv(f"{RES}/families_south.csv"); F = F[F.source == "Shelly2025"]; la0, lo0 = F.lat.mean(), F.lon.mean()
P = pd.read_csv(PN, parse_dates=["time"]); P = P[(P.time >= "2018-01-01") & (P.time <= LASTPN)]
P["d"] = hav(la0, lo0, P.lat.values, P.lon.values); P["day"] = P.time.dt.floor("D")
days = pd.date_range("2018-01-01", "2026-05-27", freq="D")
D = pd.DataFrame(index=days)
for r in (15, 30, 60): D[f"pn{r}"] = P[P.d <= r].groupby("day").size().reindex(days).fillna(0)
S = his_catalogue(); D["his"] = S.groupby(S.t.dt.floor("D")).size().reindex(days)
D.loc[D.index >= "2024-12-01", "his"] = np.nan; D.loc[D.index < "2018-01-01", "his"] = np.nan
D["his"] = D.his.where(D.index >= "2024-12-01", D.his.fillna(0)).where(D.index < "2024-12-01")
C = pd.read_parquet(f"{RES}/cut/sh_ext/conf_P_by_det_id.parquet"); C["day"] = pd.to_datetime(C.date)
dp = pd.to_datetime(pd.read_csv(f"{RES}/sh_ext2026/days_processed.csv").day)
D["oursP"] = C.groupby("day").P_same.sum().reindex(days).fillna(0).where(days.isin(dp))
DC = pd.read_csv(f"{RES}/sh_ext2026/daily_counts.csv", parse_dates=["day"]).set_index("day")
D["oursS1"] = DC[[str(i) for i in range(1, 28)]].sum(1).reindex(days).where(days.isin(dp))
rows = []
per = [("his 2018-01..2023-11", "his", "2018-01-01", "2023-11-30"), ("his join 2023-12..2024-11", "his", "2023-12-01", "2024-11-30"),
       ("ours ΣP join 2023-12..2024-11", "oursP", "2023-12-01", "2024-11-30"), ("ours S1 join", "oursS1", "2023-12-01", "2024-11-30"),
       ("ours ΣP 2024-12..2026-05-27", "oursP", "2024-12-01", "2026-05-27"), ("ours S1 2024-12..2026-05-27", "oursS1", "2024-12-01", "2026-05-27"),
       ("ours ΣP 2025", "oursP", "2025-01-01", "2025-12-31"), ("ours ΣP 2026-01..05", "oursP", "2026-01-01", "2026-05-27")]
for r in (30, 15, 60):
    for lab, c, a, b in per:
        X = D.loc[a:b, [c, f"pn{r}"]].dropna(); y = X[c]; on = X[f"pn{r}"] > 0
        rows.append(dict(radius_km=r, series=lab, n_days=len(X), n_pnsn_days=int(on.sum()), frac_pnsn_days=on.mean(),
                         mean_on=y[on].mean(), mean_off=y[~on].mean(), ratio_on_off=y[on].mean() / y[~on].mean(),
                         share_of_lfe_on_pnsn_days=y[on].sum() / y.sum(),
                         p90_lfe_day_frac_pnsn=(X[f"pn{r}"][y >= y.quantile(0.9)] > 0).mean()))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/pnsn_daily_coincidence.csv", index=False)
pd.set_option("display.width", 250); print(R.round(3).to_string())
D.to_csv(f"{OUT}/pnsn_daily.csv")
