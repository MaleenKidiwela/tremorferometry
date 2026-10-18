#!/usr/bin/env python
"""s12_recharge_check -- robustness of the task-B recharge slope (PLAN §6 B), written after the WA B number was seen
(labelled post-hoc checks, NOTES.md amendment A1). Day-level between-ETS dv/v (2-5, 4-6 s) from novel_days.parquet:
  * headline refit (linear in time since the previous major front, per-site annual cycle, series offsets, cc^2 weights);
  * + per-site linear calendar trend;  * leave-one-site-out;  * scrambled control: within each series the 'time since
    the previous front' values are shuffled across days of different inter-ETS intervals (500 draws);
  * binned curve: months since the previous front (2-4, 4-6, 6-9, 9-12, 12-18), per-site annual removed.
    S11_REGION=SVI|WA_NW python s12_recharge_check.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s12_epoch import C, OUT, SUB, TB, log
from s12_novel import within, wls

rng = np.random.default_rng(5)


def design(d, trend=False):
    cols = ["ty"]
    for s in sorted(d.site.unique()):
        d[f"c_{s}"] = d.cs * (d.site == s); d[f"s_{s}"] = d.sn * (d.site == s); cols += [f"c_{s}", f"s_{s}"]
        if trend:
            d[f"t_{s}"] = (d.year - 2015) * (d.site == s); cols.append(f"t_{s}")
    return cols


def fit(d, q, w, trend=False):
    cols = design(d, trend)
    G = within(d, [q] + cols, pd.Series(w, index=d.index))
    return wls(G[q].to_numpy(), G[cols].to_numpy(), w)[0]


def main():
    F = pd.read_parquet(f"{TB}/novel_days.parquet")
    F = F[F.betw & (F.prev < 2000)].copy()
    F["ty"] = F.prev / 365.25; F["year"] = F.date.dt.year + F.date.dt.dayofyear / 365.25
    F["iv"] = (F.date - pd.to_timedelta(F.prev, "D")).dt.strftime("%Y-%m")
    rows, curve = [], []
    for q in ("dvv_2-5", "dvv_4-6"):
        wn = q.split("_")[1]
        d = F[np.isfinite(F[q]) & (F[f"dvv_cc{wn}"] > 0.3)].copy()
        w = (d[f"dvv_cc{wn}"] ** 2).to_numpy()
        b0 = fit(d.copy(), q, w); bt = fit(d.copy(), q, w, trend=True)
        loso = {s: fit(d[d.site != s].copy(), q, w[(d.site != s).to_numpy()]) for s in sorted(d.site.unique())}
        sc = []
        for k in range(500):
            e = d.copy()
            e["ty"] = e.groupby("key").ty.transform(lambda x: rng.permutation(x.to_numpy()))
            sc.append(fit(e, q, w))
        sc = np.array(sc)
        r = dict(obs=q, slope=b0, slope_trend=bt, scr_sd=float(sc.std()), scr_p=float((np.sum(np.abs(sc) >= abs(b0)) + 1) / 501),
                 loso_min=float(min(loso.values())), loso_max=float(max(loso.values())), n_days=len(d),
                 n_series=d.key.nunique(), n_intervals=d.iv.nunique(), **{f"loso_{k}": v for k, v in loso.items()})
        rows.append(r); log(f"[B-check {SUB}] {r}")
        # binned curve with the per-site annual cycle removed
        cols = design(d)[1:]
        G = within(d, [q] + cols, pd.Series(w, index=d.index))
        ann = wls(G[q].to_numpy(), G[cols].to_numpy(), w)
        res = G[q].to_numpy() - G[cols].to_numpy() @ ann
        mo = d.prev.to_numpy() / 30.44
        for a, b in ((2, 4), (4, 6), (6, 9), (9, 12), (12, 18)):
            m = (mo >= a) & (mo < b)
            if m.sum() > 10:
                curve.append(dict(obs=q, months=f"{a}-{b}", mid=(a + b) / 2, value=float(np.sum(w[m] * res[m]) / w[m].sum()),
                                  weight=float(w[m].sum()), n_days=int(m.sum())))
    pd.DataFrame(rows).to_csv(f"{TB}/novel_B_check.csv", index=False)
    pd.DataFrame(curve).to_csv(f"{TB}/novel_B_curve.csv", index=False)


if __name__ == "__main__":
    main()
