#!/usr/bin/env python
"""Rebuild the LONG per-station stacks (44 s, starting 12 s before the node time) of every grown
family from its OWN detections.  Prefix `u_`.

Why: `k_grow.build_templates` throws the long stack away and keeps only a 12-s cut starting 2 s
before the stack's envelope maximum.  Every downstream location therefore sees ONE number per
station -- that envelope-maximum offset -- and nothing else.  A P arrival is 3-7 s BEFORE the S
here, i.e. outside the cut entirely, and an onset cannot be picked from a window that starts 2 s
before the peak by construction.  This script re-runs the identical stacking code
(k_grow.build_templates's alignment loop, verbatim in structure: linear stack, Bungum-Husebye
re-shift on the network cc, MAX_STACK=200 strongest detections) and SAVES THE LONG STACK.

Detection times in `det_it*.parquet` are NODE-equivalent times (k_grow.matched_filter shifts each
station's cc row back by that station's template offset before summing), so the same
`cut(sta, t - LONG_PRE, LONG_LEN)` convention applies unchanged.

Nothing in the detection or growth chain is modified; this is a read-only re-stack.
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
import uu_common as U                                                        # noqa: E402

ALIGN_ROUNDS = 2
PAD_S = 1.0


def norm3(c):
    n = np.sqrt((c ** 2).sum())
    return c / n if n > 0 else c


def build_long_stacks(dc, stas, groups, fs=U.FS, align_rounds=ALIGN_ROUNDS, pad_s=PAD_S,
                      log=U.log):
    nL = U.NLONG
    nP = int(round(pad_s * fs))
    byd = {}
    for g, t in groups.items():
        for x in t:
            byd.setdefault(pd.Timestamp(float(x), unit="s").strftime("%Y-%m-%d"),
                           []).append((g, float(x)))
    ST = {}
    for r in range(align_rounds + 1):
        acc, cnt = {}, {}
        for date in sorted(byd):
            dc.day(date)
            for g, x in byd[date]:
                seg = {}
                for s in stas:
                    c = dc.cut(s, x - U.LONG_PRE - pad_s, U.LONG_LEN + 2 * pad_s, date)
                    if c is None or (c ** 2).sum() <= 0:
                        continue
                    seg[s] = norm3(c)
                if not seg:
                    continue
                k = 0
                if r > 0 and g in ST:
                    ccs = np.zeros(2 * nP + 1)
                    nn = 0
                    for s, c in seg.items():
                        rv = ST[g].get(s)
                        if rv is None:
                            continue
                        ref = rv[0]
                        rn = np.sqrt((ref ** 2).sum())
                        if rn <= 0:
                            continue
                        V = np.lib.stride_tricks.sliding_window_view(c, nL, axis=1)
                        V = V[:, :2 * nP + 1, :]
                        num = np.einsum("clt,ct->l", V, ref)
                        den = np.sqrt(np.einsum("clt,clt->l", V, V))
                        ccs += np.where(den > 0, num / (rn * np.maximum(den, 1e-30)), 0.0)
                        nn += 1
                    if nn:
                        k = int(np.argmax(ccs)) - nP
                for s, c in seg.items():
                    i0 = nP + k
                    if i0 < 0 or i0 + nL > c.shape[1]:
                        continue
                    A = acc.setdefault((g, s), np.zeros((3, nL)))
                    A += c[:, i0:i0 + nL]
                    cnt[(g, s)] = cnt.get((g, s), 0) + 1
        ST = {}
        for (g, s), A in acc.items():
            if cnt[(g, s)] >= 2:
                ST.setdefault(g, {})[s] = (A / cnt[(g, s)], cnt[(g, s)])
        if log:
            log(f"[stack] round {r}: {len(ST)} groups")
    return ST


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", default="north")
    ap.add_argument("--det", required=True)
    ap.add_argument("--stas", default="")
    ap.add_argument("--tpl", default="", help="take the station set from this template npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-groups", type=int, default=0)
    a = ap.parse_args()
    B = U.bed(a.bed)
    if a.bed == "north":
        stas = a.stas.split(",") if a.stas else list(B.m.ALL2005)
        have = set(B.station_table().index)
        stas = [s for s in stas if s in have]
    else:
        B.m.activate()
        if a.tpl:
            zt = np.load(a.tpl, allow_pickle=True)
            stas = sorted({k.split("_", 1)[1] for k in zt.files if k != "meta"})
        else:
            stas = a.stas.split(",") if a.stas else list(B.m.stations())
        have = set(B.station_table().index)
        drop = [s for s in stas if s not in have]
        if drop:
            U.log(f"{len(drop)} template stations have no coordinates and are skipped: {drop}")
        stas = [s for s in stas if s in have]
    D = pd.read_parquet(a.det)
    groups = {}
    for g, gg in D.groupby("grp"):
        t = gg.nlargest(U.MAX_STACK, "cc").t.to_numpy()
        if len(t) >= 3:
            groups[int(g)] = np.sort(t)
    if a.max_groups:
        order = sorted(groups, key=lambda g: -len(groups[g]))[:a.max_groups]
        groups = {g: groups[g] for g in order}
    U.log(f"bed={a.bed} {len(stas)} stations, {len(groups)} groups, "
          f"{sum(len(v) for v in groups.values())} stacked detections")
    dc = B.DayCache(stas)
    ST = build_long_stacks(dc, stas, groups)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    payload, meta = {}, {}
    for g in ST:
        for s, (A, n) in ST[g].items():
            payload[f"{g}_{s}"] = A.astype(np.float32)
            meta[f"{g}_{s}"] = int(n)
    np.savez_compressed(a.out, meta=np.array(json.dumps(meta)), **payload)
    U.log(f"wrote {a.out}: {len(ST)} groups, {len(payload)} station stacks")


if __name__ == "__main__":
    main()
