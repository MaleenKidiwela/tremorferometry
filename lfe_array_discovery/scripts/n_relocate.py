#!/usr/bin/env python
"""RE-STACK a template on its own detections and RE-LOCATE it.  Prefix `n_`.

Agent M's `m_relocate.py` with the bed adapter swapped to `n_common` (Puget/S-Olympic 2008 and
the 2007-01 confirmation) and one option added: `--half 1|2` restacks from the FIRST or SECOND
half of the template's detections in time, which is the stability clause of the 2008 screen
(re-stack from the second half, re-locate, the location must reproduce within 10 km).

Why re-stack.  `l_locate`/`k_locate` pick S from the TEMPLATE's own per-station stack, so locating
the same template file twice gives the same answer by construction.  The 2011 question -- "does
this family locate in the same place from 2011 data?" -- is only answerable from a stack built at
the 2011 detection times.  The stacking code is `k_grow.build_templates` VERBATIM (Bungum-Husebye
realignment, 12 s cut 2 s before the station's stack-envelope peak, >= 4 stations, >= 3 MAD peak);
the member set is the same `nlargest(MAX_STACK, cc)` rule `k_grow` itself uses between iterations.
Both years go through this identical path, so 2010 and 2011 locations are comparable.

The grid is `n_locate`'s, imported from it so the two can never drift: 46.5-48.8 N /
-125.2..-121.8 / 15-60 km at 0.03 deg, 0.04 deg, 2 km, with the same `edge` flag.

Everything else is k_locate: >= 5 picked stations at >= 4 MAD, demeaned-residual rms on the
Crosson 1-D S times with the intraslab S station terms, and the 50-realisation SHUFFLED-STATION
control.

  python scripts/n_relocate.py --det <det.parquet> --grps <csv|list> --out <prefix> [--half 2]
  N_BED=2007 python scripts/n_relocate.py ...                      # the confirmation bed
  N_NULL_SHIFT=1 python scripts/n_relocate.py ...                  # the site-shift null chain
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import k_grow as KG                                                         # noqa: E402
import k_locate as KL                                                       # noqa: E402
import n_common as LC                                                       # noqa: E402
import n_locate as NL                                                       # noqa: E402

grid = NL.grid
is_edge = NL.is_edge


def tt_table(stas, lats, lons, deps, corr, ST):
    """(n_node, n_sta) predicted S -- k_locate.tt_table's computation, distance table per depth."""
    import s_common as S
    import p_vel1d as V
    LA, LO, DP = np.meshgrid(lats, lons, deps, indexing="ij")
    LA, LO, DP = LA.ravel(), LO.ravel(), DP.ravel()
    dgrid = np.arange(0.0, 400.1, 0.5)
    tabs = {float(dp): np.array([V.tt(float(x), float(dp), "S") for x in dgrid])
            for dp in np.unique(DP)}
    out = np.zeros((LA.size, len(stas)), np.float32)
    for k, s in enumerate(stas):
        d = S.epi_km(LA, LO, float(ST.lat[s]), float(ST.lon[s]))
        for dp, tab in tabs.items():
            m = DP == dp
            out[m, k] = np.interp(d[m], dgrid, tab) + float(corr.get(s, 0.0))
    return np.column_stack([LA, LO, DP]), out


def read_grps(spec):
    if not spec:
        return None
    if os.path.exists(spec):
        V = pd.read_csv(spec)
        if "candidate" in V.columns:
            V = V[V.candidate.astype(bool)]
        return sorted(int(x) for x in V.grp)
    return sorted(int(x) for x in spec.split(","))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True)
    ap.add_argument("--half", type=int, default=0, choices=[0, 1, 2],
                    help="0 = all detections; 1/2 = the first/second half in time")
    ap.add_argument("--grps", default="")
    ap.add_argument("--out", required=True, help="output prefix; writes <p>_tpl.npz and <p>.csv")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    tag = a.tag or os.path.basename(a.out)

    import n_common as B
    K = B.activate()                         # honours N_BED and N_NULL_SHIFT
    import g_common as G
    stas = list(K.CLEAN11)
    log = B.log
    log(f"[{tag}] bed {B.BED}, {len(stas)} stations, half={a.half}, null_shift="
        f"{os.environ.get('N_NULL_SHIFT', '0')}")

    D = pd.read_parquet(a.det)
    grps = read_grps(a.grps)
    if grps is not None:
        D = D[D.grp.isin(grps)]
    days = set(K.bed_days())
    D = D[D.date.isin(days)]
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
    log(f"[{tag}] {len(groups)} groups to re-stack ({len(D)} detections); "
        f"{len(skipped)} with < 3 detections: {skipped[:20]}")

    nd = D.groupby("grp").size()
    ndays = D.groupby("grp").date.nunique()

    rows = []
    if groups:
        dc = K.DayCache(stas)
        TPL, info = KG.build_templates(dc, stas, groups)
        np.savez_compressed(f"{a.out}_tpl.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(json.dumps(
                                {f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                 for g in TPL for s in TPL[g]})))
        info.to_csv(f"{a.out}_stackinfo.csv", index=False)

        P = {g: {s: (TPL[g][s][1] + KL.TPL_PRE, TPL[g][s][2]) for s in TPL[g]} for g in TPL}
        tstas = sorted({s for g in P for s in P[g]})
        corr, flagged = LC.station_corrections()
        ST = LC.station_table()
        lats, lons, deps = grid()
        log(f"[{tag}] grid {lats.size}x{lons.size}x{deps.size} = "
            f"{lats.size*lons.size*deps.size} nodes, {len(tstas)} template stations; "
            f"intraslab S terms for {len([s for s in tstas if s not in flagged])}, "
            f"0.0-FLAGGED for {len([s for s in tstas if s in flagged])}")
        NODE, TTg = tt_table(tstas, lats, lons, deps, corr, ST)
        si = {s: i for i, s in enumerate(tstas)}
        rng = np.random.default_rng(3)
        for g in sorted(P):
            ss = [s for s, (o, h) in P[g].items() if h >= KL.MIN_PEAK and s in si]
            base = dict(grp=g, n_sta=len(ss), n_det=int(nd.get(g, 0)),
                        n_days=int(ndays.get(g, 0)), n_stack=len(groups.get(g, [])))
            if len(ss) < KL.MIN_STA:
                rows.append({**base, "located": False})
                continue
            obs = np.array([P[g][s][0] for s in ss])
            idx = np.array([si[s] for s in ss])
            i, r, _ = KL.fit(obs, TTg, idx)
            rs = np.array([KL.fit(obs, TTg, idx[rng.permutation(len(idx))])[1]
                           for _ in range(KL.N_SHUF)])
            la, lo, dp = NODE[i]
            sz = float(np.asarray(G.slab2_depth(np.array([la]), np.array([lo]))).ravel()[0])
            rows.append({**base, "located": True, "lat": float(la), "lon": float(lo),
                         "dep": float(dp), "rms": r,
                         "rms_shuf_med": float(np.median(rs)),
                         "rms_shuf_p05": float(np.quantile(rs, .05)),
                         "ctrl_pass": bool(r < np.quantile(rs, .05)),
                         "slab2": sz,
                         "dz": (float(dp) - sz) if np.isfinite(sz) else np.nan,
                         "edge": is_edge(la, lo, dp, lats, lons, deps),
                         "n_corr": int(sum(1 for s in ss if s not in flagged))})
    for g in skipped:
        rows.append(dict(grp=g, n_sta=0, n_det=int(nd.get(g, 0)),
                         n_days=int(ndays.get(g, 0)), n_stack=0, located=False))
    if grps is not None:
        for g in sorted(set(grps) - {r["grp"] for r in rows}):
            rows.append(dict(grp=g, n_sta=0, n_det=0, n_days=0, n_stack=0, located=False))
    R = pd.DataFrame(rows).sort_values("grp")
    if "edge" in R.columns:
        R["edge"] = R.edge.map(lambda v: bool(v) if v is not None and v == v else None)
    R.to_csv(f"{a.out}.csv", index=False)
    Lo = R[R.located == True]                                               # noqa: E712
    ed = Lo.edge.astype(bool) if len(Lo) else pd.Series(dtype=bool)
    log(f"[{tag}] {len(R)} templates; located {len(Lo)}; control pass "
        f"{int(Lo.ctrl_pass.astype(bool).sum()) if len(Lo) else 0}; "
        f"on a grid edge {int(ed.sum())}; interior {int((~ed).sum())}")


if __name__ == "__main__":
    main()
