#!/usr/bin/env python
"""SCREEN-VALIDATION data stage (prefix `q_sv_`).  Produces the two things the 2005 POLARIS bed
is missing before the five-clause Washington screen can be applied to it:

  (a) QUIET-DAY detections -- the unchanged K_REPRO3 15-station templates (grow15/tpl_it0.npz)
      re-fired, at the identical 8 x MAD rule, over the 20 literature-quiet 2005 days
      (2005-06-25..2005-07-14, k_common.QUIET_DAYS).  Nothing about the templates changes.
  (e) FIRST-HALF and SECOND-HALF re-stacks of every grown template on its own ETS detections,
      re-located on k_locate's own grid with the same 50-realisation shuffled-station control.
      k_grow.build_templates is used verbatim, exactly as n_relocate.py does for the 2008 bed.

  python scripts/q_sv_data.py --stage quiet
  python scripts/q_sv_data.py --stage restack --half 1|2
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import k_common as K                                                        # noqa: E402
import k_grow as KG                                                         # noqa: E402
import k_locate as KL                                                       # noqa: E402

S15 = ["GLBC", "TWKB", "MGCB", "TSJB", "TWBB", "PGC", "SILB", "SSIB", "LCBC",
       "B003", "B005", "JRBC", "LZB", "VGZ", "SNB"]
TPLF = "results/repro/grow15/tpl_it0.npz"
DETF = "results/repro/grow15/det_it0.parquet"
OUT = "results/q_screenval"


def load_tpl(path):
    z = np.load(path, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        TPL.setdefault(int(g), {})[s] = (np.asarray(z[k]), float(v[0]), float(v[1]))
    return TPL


def stage_quiet(a):
    TPL = load_tpl(TPLF)
    dates = list(K.QUIET_DAYS)
    dc = K.DayCache(S15)
    avail = []
    for d in dates:
        n = len(dc.day(d))
        avail.append(dict(date=d, n_sta=n, used=bool(n >= 4)))
        K.log(f"[quiet-cache] {d} {n} stations")
    A = pd.DataFrame(avail)
    A.to_csv(f"{OUT}/quiet_day_stations.csv", index=False)
    use = [r.date for r in A.itertuples() if r.used]
    K.log(f"{len(TPL)} templates, {len(use)} usable quiet days of {len(dates)}")
    D = KG.matched_filter(dc, S15, TPL, use)
    D.to_parquet(f"{OUT}/det_quiet15.parquet", index=False)
    K.log(f"[quiet] {len(D)} detections, {D.grp.nunique() if len(D) else 0} templates firing; "
          f"network rate {len(D) / max(len(TPL) * len(use), 1):.4f} per template-day")


def grid():
    return (np.arange(47.9, 49.31, 0.03), np.arange(-124.9, -122.79, 0.04),
            np.arange(18.0, 56.1, 2.0))


def is_edge(la, lo, dp, lats, lons, deps, tol=1e-6):
    return bool(abs(la - lats[0]) < tol or abs(la - lats[-1]) < tol
                or abs(lo - lons[0]) < tol or abs(lo - lons[-1]) < tol
                or abs(dp - deps[0]) < tol or abs(dp - deps[-1]) < tol)


def stage_restack(a):
    import p_common as P_
    sl = P_.slab2_interp()      # the SAME slab2 k_locate used for loc15_it0
    D = pd.read_parquet(DETF)
    groups, skipped = {}, []
    for g, gg in D.groupby("grp"):
        if a.half:
            gg = gg.sort_values("t")
            h = len(gg) // 2
            gg = gg.iloc[:h] if a.half == 1 else gg.iloc[h:]
        t = gg.nlargest(KG.MAX_STACK, "cc").t.to_numpy()
        if len(t) >= 3:
            groups[int(g)] = np.sort(t)
        else:
            skipped.append(int(g))
    K.log(f"[half{a.half}] {len(groups)} groups to re-stack, {len(skipped)} with <3 detections")
    nd = D.groupby("grp").size(); ndays = D.groupby("grp").date.nunique()
    dc = K.DayCache(S15)
    TPL, info = KG.build_templates(dc, S15, groups)
    pref = f"{OUT}/reloc_h{a.half}"
    np.savez_compressed(f"{pref}_tpl.npz",
                        **{f"{g}_{s}": TPL[g][s][0].astype(np.float32) for g in TPL for s in TPL[g]},
                        meta=np.array(json.dumps({f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                                  for g in TPL for s in TPL[g]})))
    P = {g: {s: (TPL[g][s][1] + KL.TPL_PRE, TPL[g][s][2]) for s in TPL[g]} for g in TPL}
    tstas = sorted({s for g in P for s in P[g]})
    lats, lons, deps = grid()
    K.log(f"[half{a.half}] grid {lats.size}x{lons.size}x{deps.size}, {len(tstas)} stations")
    NODE, TTg = KL.tt_table(tstas, lats, lons, deps)
    si = {s: i for i, s in enumerate(tstas)}
    rng = np.random.default_rng(3)
    rows = []
    for g in sorted(P):
        ss = [s for s, (o, h) in P[g].items() if h >= KL.MIN_PEAK and s in si]
        base = dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)), n_days=int(ndays.get(g, 0)),
                    n_stack=len(groups.get(g, [])))
        if len(ss) < KL.MIN_STA:
            rows.append({**base, "located": False}); continue
        obs = np.array([P[g][s][0] for s in ss]); idx = np.array([si[s] for s in ss])
        i, r, _ = KL.fit(obs, TTg, idx)
        rs = np.array([KL.fit(obs, TTg, idx[rng.permutation(len(idx))])[1]
                       for _ in range(KL.N_SHUF)])
        la, lo, dp = NODE[i]
        try:
            sz = float(np.asarray(sl(np.array([la]), np.array([lo]))).ravel()[0])
        except Exception:
            sz = np.nan
        rows.append({**base, "located": True, "lat": float(la), "lon": float(lo),
                     "dep": float(dp), "rms": r, "rms_shuf_med": float(np.median(rs)),
                     "rms_shuf_p05": float(np.quantile(rs, .05)),
                     "ctrl_pass": bool(r < np.quantile(rs, .05)), "slab2": sz,
                     "dz": (float(dp) - sz) if np.isfinite(sz) else np.nan,
                     "edge": is_edge(la, lo, dp, lats, lons, deps)})
    for g in skipped:
        rows.append(dict(grp=g, n_sta=0, n_det=int(nd.get(g, 0)), n_days=int(ndays.get(g, 0)),
                         n_stack=0, located=False))
    R = pd.DataFrame(rows).sort_values("grp")
    R.to_csv(f"{pref}.csv", index=False)
    Lo = R[R.located == True]                                               # noqa: E712
    K.log(f"[half{a.half}] {len(R)} templates; located {len(Lo)}; ctrl_pass "
          f"{int(Lo.ctrl_pass.sum()) if len(Lo) else 0}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["quiet", "restack"])
    ap.add_argument("--half", type=int, default=2, choices=[0, 1, 2])
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    (stage_quiet if a.stage == "quiet" else stage_restack)(a)
