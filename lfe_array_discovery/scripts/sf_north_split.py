#!/usr/bin/env python
"""STAGE 1 -- calibrate the CONFIRMATION rule on the north bed (agent R's never-executed brief).

The Washington confirmation rules used so far (>= 3 tremor days in the sealed year, quiet-day FA
<= 0.1/template-day, an INTERIOR location within 10 km of the discovery location) were never
measured against a known catalogue.  Here they are, on the only bed with truth: the 24-day 2005
POLARIS ETS window is split into calendar halves H1 (first 12 days) and H2 (last 12), every
template is RE-STACKED FROM ITS H1 DETECTIONS ONLY and re-fired, unchanged, on H2 and on the 20
literature-quiet 2005 days.  That is exactly what a Washington confirmation does -- a template
built on one window, fired on another it never saw -- with the difference that here we know which
templates belong to Bostock's own families.

Populations: the 177 templates of the 58 blind-rediscovered KEPT families (`kept_family` in
results/repro/classify15.csv) and a matched-size random sample of the 1,201 unlocatable/unmatched
templates (the noise population of results/SCREEN_VALIDATION_NORTH.md).

Nothing in scripts/k_*.py is edited: `k_grow.build_templates`, `k_grow.matched_filter` (8 x MAD)
and `k_locate` are called verbatim.

  python scripts/sf_north_split.py --stage select [--limit 3]
  python scripts/sf_north_split.py --stage h1
  python scripts/sf_north_split.py --stage refire --which h2|quiet
  python scripts/sf_north_split.py --stage h2stack
  python scripts/sf_north_split.py --stage locate --half 1|2
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402
import k_grow as KG                                                         # noqa: E402
import k_locate as KL                                                       # noqa: E402
import sf_common as SF                                                      # noqa: E402

S15 = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC",
       "B003", "B005", "JRBC", "LZB", "VGZ", "SNB"]
TPLF = "results/repro/grow15/tpl_it0.npz"
DETF = "results/repro/grow15/det_it0.parquet"
CLS = "results/repro/classify15.csv"
OUT = f"{SF.RES}/sf_final/north"
SEED = 20260923


def days():
    d = [x.strftime("%Y-%m-%d") for x in pd.date_range("2005-08-29", "2005-09-21")]
    return d[:12], d[12:]


def quiet_days_usable(dc):
    use, rows = [], []
    for d in K.QUIET_DAYS:
        n = len(dc.day(d))
        rows.append(dict(date=d, n_sta=n, used=bool(n >= 4)))
        if n >= 4:
            use.append(d)
    pd.DataFrame(rows).to_csv(f"{OUT}/quiet_days.csv", index=False)
    return use


# ---------------------------------------------------------------- select ----------------------
def stage_select(a):
    C = pd.read_csv(CLS)
    true = C[C.cls == "kept_family"].copy()
    noise_pool = C[C.cls == "unlocatable"].copy()
    rng = np.random.default_rng(SEED)
    k = min(len(true), len(noise_pool))
    idx = rng.choice(len(noise_pool), size=k, replace=False)
    noise = noise_pool.iloc[np.sort(idx)].copy()
    true["pop"] = "true"
    true["bkey"] = true.kept_bkey
    noise["pop"] = "noise"
    noise["bkey"] = -1
    S = pd.concat([true, noise], ignore_index=True)[
        ["grp", "pop", "bkey", "n_det", "n_days", "located", "lat", "lon", "dep"]]
    if a.limit:
        # smoke: the a.limit largest TRUE families plus the same number of noise templates
        t = S[S["pop"] == "true"].nlargest(a.limit, "n_det")
        nz = S[S["pop"] == "noise"].nlargest(a.limit, "n_det")
        S = pd.concat([t, nz], ignore_index=True)
    S = S.sort_values("grp").reset_index(drop=True)
    S.to_csv(f"{OUT}/set{'_smoke' if a.limit else ''}.csv", index=False)
    K.log(f"[select] {len(S)} templates: {(S['pop'] == 'true').sum()} true "
          f"({S[S['pop'] == 'true'].bkey.nunique()} families), "
          f"{(S['pop'] == 'noise').sum()} noise; seed {SEED}")


def setfile(a):
    return f"{OUT}/set{'_smoke' if a.limit else ''}.csv"


def sfx(a):
    return "_smoke" if a.limit else ""


# ---------------------------------------------------------------- H1 stacks -------------------
def stage_h1(a):
    S = pd.read_csv(setfile(a))
    grps = set(int(g) for g in S.grp)
    D = pd.read_parquet(DETF)
    H1, _ = days()
    d1 = D[D.grp.isin(grps) & D.date.isin(H1)]
    groups, info = {}, []
    for g in sorted(grps):
        gg = d1[d1.grp == g]
        t = gg.nlargest(KG.MAX_STACK, "cc").t.to_numpy() if len(gg) else np.array([])
        info.append(dict(grp=g, n_det_h1=len(gg),
                         n_days_h1=int(gg.date.nunique()) if len(gg) else 0,
                         n_stack_h1=len(t) if len(t) >= 3 else 0))
        if len(t) >= 3:
            groups[g] = np.sort(t)
    pd.DataFrame(info).to_csv(f"{OUT}/h1_stackinfo{sfx(a)}.csv", index=False)
    K.log(f"[h1] {len(groups)}/{len(grps)} templates have >= 3 H1 detections to re-stack")
    dc = K.DayCache(S15)
    TPL, _ = KG.build_templates(dc, S15, groups)
    SF.save_tpl(f"{OUT}/tpl_h1{sfx(a)}.npz", TPL)
    K.log(f"[h1] {len(TPL)} H1-only templates written")


# ---------------------------------------------------------------- re-fire ---------------------
def stage_refire(a):
    TPL = SF.load_tpl(f"{OUT}/tpl_h1{sfx(a)}.npz")
    dc = K.DayCache(S15)
    H1, H2 = days()
    if a.which == "h2":
        dates = H2
    else:
        dates = quiet_days_usable(dc)
    K.log(f"[refire {a.which}] {len(TPL)} templates x {len(dates)} days")
    D = KG.matched_filter(dc, S15, TPL, dates)
    D.to_parquet(f"{OUT}/det_{a.which}{sfx(a)}.parquet", index=False)
    json.dump(dict(which=a.which, n_templates=len(TPL), n_dates=len(dates), dates=dates,
                   n_det=int(len(D)),
                   per_tpl_day=float(len(D) / max(len(TPL) * len(dates), 1))),
              open(f"{OUT}/det_{a.which}{sfx(a)}_summary.json", "w"), indent=1)
    K.log(f"[refire {a.which}] {len(D)} detections "
          f"({len(D) / max(len(TPL) * len(dates), 1):.3f} per template-day)")


# ---------------------------------------------------------------- H2 stacks -------------------
def stage_h2stack(a):
    D = pd.read_parquet(f"{OUT}/det_h2{sfx(a)}.parquet")
    groups, info = {}, []
    for g, gg in D.groupby("grp"):
        t = gg.nlargest(KG.MAX_STACK, "cc").t.to_numpy()
        info.append(dict(grp=int(g), n_det_h2=len(gg), n_days_h2=int(gg.date.nunique()),
                         n_stack_h2=len(t) if len(t) >= 3 else 0))
        if len(t) >= 3:
            groups[int(g)] = np.sort(t)
    pd.DataFrame(info).to_csv(f"{OUT}/h2_stackinfo{sfx(a)}.csv", index=False)
    K.log(f"[h2stack] {len(groups)} templates re-stack from their own H2 detections")
    dc = K.DayCache(S15)
    TPL, _ = KG.build_templates(dc, S15, groups)
    SF.save_tpl(f"{OUT}/tpl_h2{sfx(a)}.npz", TPL)
    K.log(f"[h2stack] {len(TPL)} H2-only templates written")


# ---------------------------------------------------------------- locate ----------------------
def grid():
    return (np.arange(47.9, 49.31, 0.03), np.arange(-124.9, -122.79, 0.04),
            np.arange(18.0, 56.1, 2.0))


def is_edge(la, lo, dp, lats, lons, deps, tol=1e-6):
    return bool(abs(la - lats[0]) < tol or abs(la - lats[-1]) < tol
                or abs(lo - lons[0]) < tol or abs(lo - lons[-1]) < tol
                or abs(dp - deps[0]) < tol or abs(dp - deps[-1]) < tol)


def stage_locate(a):
    TPL = SF.load_tpl(f"{OUT}/tpl_h{a.half}{sfx(a)}.npz")
    P = {g: {s: (TPL[g][s][1] + KL.TPL_PRE, TPL[g][s][2]) for s in TPL[g]} for g in TPL}
    tstas = sorted({s for g in P for s in P[g]})
    lats, lons, deps = grid()
    K.log(f"[loc h{a.half}] {len(P)} templates, grid {lats.size}x{lons.size}x{deps.size}, "
          f"{len(tstas)} stations")
    NODE, TTg = KL.tt_table(tstas, lats, lons, deps)
    si = {s: i for i, s in enumerate(tstas)}
    rng = np.random.default_rng(3)
    rows = []
    for g in sorted(P):
        ss = [s for s, (o, h) in P[g].items() if h >= KL.MIN_PEAK and s in si]
        base = dict(grp=int(g), n_sta=len(ss))
        if len(ss) < KL.MIN_STA:
            rows.append({**base, "located": False})
            continue
        obs = np.array([P[g][s][0] for s in ss])
        idx = np.array([si[s] for s in ss])
        i, r, _ = KL.fit(obs, TTg, idx)
        rs = np.array([KL.fit(obs, TTg, idx[rng.permutation(len(idx))])[1]
                       for _ in range(KL.N_SHUF)])
        la, lo, dp = NODE[i]
        rows.append({**base, "located": True, "lat": float(la), "lon": float(lo),
                     "dep": float(dp), "rms": r,
                     "ctrl_pass": bool(r < np.quantile(rs, .05)),
                     "edge": is_edge(la, lo, dp, lats, lons, deps)})
    R = pd.DataFrame(rows).sort_values("grp")
    R.to_csv(f"{OUT}/loc_h{a.half}{sfx(a)}.csv", index=False)
    K.log(f"[loc h{a.half}] {len(R)} templates; located {int(R.located.sum())}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["select", "h1", "refire", "h2stack", "locate"])
    ap.add_argument("--which", default="h2", choices=["h2", "quiet"])
    ap.add_argument("--half", type=int, default=1, choices=[1, 2])
    ap.add_argument("--limit", type=int, default=0, help="smoke test: this many true families")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    {"select": stage_select, "h1": stage_h1, "refire": stage_refire,
     "h2stack": stage_h2stack, "locate": stage_locate}[a.stage](a)


if __name__ == "__main__":
    main()
