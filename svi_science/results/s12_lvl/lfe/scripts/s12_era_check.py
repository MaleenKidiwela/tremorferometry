#!/usr/bin/env python
"""s12_era_check -- post hoc (2026-10-14): SVI recharge slope by era (2003-2013 / 2014-2026) as in s12_recharge2, plus
(i) a scrambled control per era (time since ETS shuffled within series, 300 draws) and (ii) the same fits with the B009
instrument-veto window (2007-01-08..2009-10-01, svi_science/instrument_vetoes.csv) removed. Writes tables/svi/r2_era_check.json.
    S11_REGION=SVI python s12_era_check.py
"""
import os, sys, json
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s12_recharge2 as R2
from s12_epoch import C, TB, SUB

F = pd.read_parquet(f"{TB}/b_days.parquet"); F["date"] = pd.to_datetime(F.date)
F["site"] = [C.SITE.get(n, n) for n in F.netsta]
ph = 2 * np.pi * F.date.dt.dayofyear.to_numpy() / 365.25
F["cs"], F["sn"], F["cs2"], F["sn2"] = np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph)
F["ty"] = F.prev / 365.25
B = F[F.betw & (F.prev < 2000)].copy(); q = "dvv_2-5"; d = R2.prep(B, q)
V = pd.read_csv(f"{C.ROOT}/svi_science/instrument_vetoes.csv"); V = V[V.decision.str.strip() == "veto"]
veto = np.zeros(len(d), bool)
for r in V.itertuples():
    veto |= ((d.netsta == r.netsta) & (d.date >= r.start) & (d.date <= r.end)).to_numpy()
out = dict(n_vetoed_days=int(veto.sum()), n_days=len(d))
for lab, dd in (("asis", d), ("veto", d[~veto])):
    for era, m in (("2003-2013", dd.date.dt.year <= 2013), ("2014-2026", dd.date.dt.year >= 2014), ("all", dd.date.dt.year > 0)):
        f = R2.fit(dd[m], q)
        if lab == "asis":
            f.update({k: v for k, v in R2.scramble(dd[m], q, ["key"], n=300).items() if k != "b"})
        out[f"{lab}_{era}"] = f
        print(lab, era, {k: round(v, 4) if isinstance(v, float) else v for k, v in f.items()}, flush=True)
json.dump(out, open(f"{TB}/r2_era_check.json", "w"), indent=1, default=float)
