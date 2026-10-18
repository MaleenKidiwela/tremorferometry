#!/usr/bin/env python
"""P4-0b growth driver: quiet-day thresholds -> growth pass -> survival.

  quiet     run the bank over the quiet days at the 8.0 floor and derive T_seed per template by
            the FROZEN FA-only rule (c_tiers_by_epoch.py): smallest T on a 0.25 grid, >= 8.0, with
            <= 0.1 false hits per template-day.
  grow      run the bank over the pooled episode days at each template's own T_seed.
  survive   apply Merlin's survival clauses (>= 20 independent hits; >= 3 distinct days and >= 2
            episodes; <= 20 % of hits within +-10 min of the seed).  The N-moveout ratio is added
            by `--nmoveout-hits`.

  python scripts/v_run.py quiet --bank ... --out ...
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import v_engine as E                                                        # noqa: E402

OUT = V.OUT


def _cap(days):
    """V_MAX_DAYS=<n> truncates every day list -- SMOKE TESTS ONLY; never set in the queue."""
    k = int(os.environ.get("V_MAX_DAYS", "0") or 0)
    return days[:k] if k else days


def _quiet_days():
    p = f"{OUT}/quiet_days.csv"
    if os.path.exists(p):
        Q = pd.read_csv(p)
        return _cap([pd.Timestamp(x) for x in Q.date])
    Q = V.quiet_days()
    Q.to_csv(p, index=False)
    return _cap(list(Q.date))


def _episode_days():
    ED = V.episode_days()
    A = V.availability()
    usable = {t for t, u in zip(A.date, A.usable) if u}
    return _cap([t for t in sorted(ED) if t in usable])


def _n_obs_days(bank_stas, days):
    idx = V.cache_idx()
    n = 0
    for d in days:
        av = {V.SITE_OF[s] for s in bank_stas if (s, int(d.year), int(d.dayofyear)) in idx}
        if len(av) >= V.MIN_SITES:
            n += 1
    return n


def stage_quiet(bank, out):
    tids, stas, W, off, meta = V.load_bank(bank)
    days = _quiet_days()
    H = E.grow(bank, days, out + "_quiethits.parquet", thr=None, floor=V.T_FLOOR, tag="/quiet")
    nday = _n_obs_days(stas, days)
    grid = np.arange(V.T_FLOOR, V.T_MAX + 1e-9, V.T_GRID)
    rows = []
    for bi, tid in enumerate(tids):
        v = H[H.tpl == bi].net_mad.to_numpy() if len(H) else np.zeros(0)
        fa = np.array([(v >= t).sum() / max(nday, 1) for t in grid])
        ok = np.nonzero(fa <= V.FA_PER_TPL_DAY)[0]
        T = float(grid[ok[0]]) if len(ok) else float(V.T_MAX)
        rows.append(dict(tpl=bi, tid=tid, T_seed=T, n_quiet_hits=int(len(v)),
                         quiet_days=nday, fa_at_T=float(fa[ok[0]]) if len(ok) else float(fa[-1]),
                         saturated=bool(len(ok) == 0)))
    Tt = pd.DataFrame(rows)
    Tt.to_csv(out + "_thresholds.csv", index=False)
    V.log(f"[quiet] {len(Tt)} templates, {nday} quiet days; T_seed median "
          f"{Tt.T_seed.median():.2f} range {Tt.T_seed.min():.2f}-{Tt.T_seed.max():.2f}; "
          f"saturated {int(Tt.saturated.sum())}")
    return Tt


def stage_grow(bank, out, thr_csv=""):
    tids, stas, W, off, meta = V.load_bank(bank)
    thr_csv = thr_csv or (out + "_thresholds.csv")
    T = pd.read_csv(thr_csv).set_index("tpl").T_seed.reindex(range(len(tids))).to_numpy()
    days = _episode_days()
    E.grow(bank, days, out + "_hits.parquet", thr=T, floor=V.T_FLOOR, tag="/grow")


def stage_survive(bank, out, seeds_path, nmoveout_out="", nmoveout_bank=""):
    tids, stas, W, off, meta = V.load_bank(bank)
    H = pd.read_parquet(out + "_hits.parquet")
    Sd = pd.read_parquet(seeds_path).set_index("seed")
    ED = V.episode_days()
    NM = None
    if nmoveout_out and os.path.exists(nmoveout_out + "_hits.parquet") and nmoveout_bank:
        NH = pd.read_parquet(nmoveout_out + "_hits.parquet")
        ntids, *_ = V.load_bank(nmoveout_bank)
        cnt = {}
        for bi, g in NH.groupby("tpl"):
            m = V.independent(g.sort_values("t").t.to_numpy())
            cnt[ntids[int(bi)]] = int(m.sum())
        NM = cnt
    rows = []
    indep_rows = []
    for bi, tid in enumerate(tids):
        g = H[H.tpl == bi].sort_values("t") if len(H) else H
        if len(g) == 0:
            rows.append(dict(tpl=bi, tid=tid, n_hits=0, n_indep=0, n_days=0, n_episodes=0,
                             frac_near_seed=np.nan, n_moveout=NM.get(tid, np.nan) if NM else np.nan,
                             moveout_ratio=np.nan, survives=False))
            continue
        m = V.independent(g.t.to_numpy())
        gi = g[m].copy()
        gi["tid"] = tid
        indep_rows.append(gi)
        eps = {ED.get(pd.Timestamp(d), "") for d in gi.date}
        t_seed = float(Sd.t[tid]) if tid in Sd.index else np.nan
        near = float(np.mean(np.abs(gi.t.to_numpy() - t_seed) <= V.NEAR_SEED_S)) \
            if np.isfinite(t_seed) else np.nan
        nm = NM.get(tid, np.nan) if NM else np.nan
        ratio = (len(gi) / max(nm, 1)) if (NM and np.isfinite(nm)) else np.nan
        surv = bool(len(gi) >= V.SURV_MIN_INDEP and gi.date.nunique() >= V.SURV_MIN_DAYS
                    and len(eps - {""}) >= V.SURV_MIN_EPISODES
                    and (not np.isfinite(near) or near <= V.SURV_MAX_NEAR_SEED)
                    and (NM is None or (np.isfinite(ratio) and ratio >= V.NMOVEOUT_RATIO)))
        rows.append(dict(tpl=bi, tid=tid, n_hits=int(len(g)), n_indep=int(len(gi)),
                         n_days=int(gi.date.nunique()), n_episodes=int(len(eps - {""})),
                         frac_near_seed=near, n_moveout=nm, moveout_ratio=ratio, survives=surv))
    R = pd.DataFrame(rows)
    R.to_csv(out + "_survival.csv", index=False)
    if indep_rows:
        pd.concat(indep_rows, ignore_index=True).to_parquet(out + "_indep.parquet", index=False)
    V.log(f"[survive] {len(R)} templates -> {int(R.survives.sum())} survivors "
          f"(>=20 indep {int((R.n_indep>=V.SURV_MIN_INDEP).sum())}, "
          f">=3 days {int((R.n_days>=V.SURV_MIN_DAYS).sum())}, "
          f">=2 episodes {int((R.n_episodes>=V.SURV_MIN_EPISODES).sum())})")
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["quiet", "grow", "survive"])
    ap.add_argument("--bank", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--thr", default="")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--nmoveout-out", default="")
    ap.add_argument("--nmoveout-bank", default="")
    a = ap.parse_args()
    if a.stage == "quiet":
        stage_quiet(a.bank, a.out)
    elif a.stage == "grow":
        stage_grow(a.bank, a.out, a.thr)
    else:
        stage_survive(a.bank, a.out, a.seeds, a.nmoveout_out, a.nmoveout_bank)


if __name__ == "__main__":
    main()
