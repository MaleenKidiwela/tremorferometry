#!/usr/bin/env python
"""bw_gate_retime -- re-time Bostock catalogue events onto the FROZEN detector time axis with the frozen 2012 templates.

Why: Bostock's 2004 and 2005 event times are TRUNCATED to whole seconds (every 'sec' value is an integer in the
catalogue; 2003 and 2007-2013 carry 0.05-s precision).  At the never-changed PGC the frozen template then sits +0.45 s
late with ~0.3 s scatter (bw_gate_g3_templates.py check), which smears any stack built at his times.

retime(): per event, t0 = catalogue t_origin (+0.5 s when the year is truncated), then the network cc of the frozen
12-s templates of the TIMING stations (never the station under test) at their frozen placements, weighted 3 w as in
the detector, is scanned over +-0.6 s (20 Hz); the event moves to the argmax.  Columns added: t_cat, lag, net_cc.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G

TRUNC_YEARS = (2004, 2005)
_FZ = {}


def frozen():
    if not _FZ:
        R, C = G.patch_forward()
        fz = C.frozen(); L = float(fz["template_length_s"]); tpath, _ = C.resolve_templates(fz, "")
        tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpath, L)
        _FZ.update(tpl=tpl, off=C.placement(tpl, sp, F, mv, C.station_table(), L), wgt=wgt, L=L)
    return _FZ


def retime(ev, year, caches, stations, half=0.6):
    import e5r_templates as T
    from numpy.lib.stride_tricks import sliding_window_view as swv
    Z = frozen(); tpl, off, wgt = Z["tpl"], Z["off"], Z["wgt"]
    T.CACHES = caches; T.YEAR = year
    D = {s: T.load_days(G.NETSTA[s], year) for s in stations}
    h = int(round(half * 20)); n = int(round(Z["L"] * 20))
    ev = ev.copy().reset_index(drop=True)
    ev["t_cat"] = ev.t_origin
    t0 = ev.t_origin.to_numpy() + (0.5 if year in TRUNC_YEARS else 0.0)
    acc = np.zeros((len(ev), 2 * h + 1)); wsum = np.zeros(len(ev))
    for s in stations:
        if not D[s]:
            continue
        for f, idx in ev.groupby("fam").indices.items():
            k = (float(f), s)
            if k not in tpl or k not in off:
                continue
            w = np.asarray(tpl[k][Z["L"]], float); wn = np.sqrt((w * w).sum())
            Sw, keep = T.cut_windows(D[s], t0[idx] + off[k] - half, n + 2 * h)
            if Sw is None:
                continue
            V3 = swv(Sw, n, axis=2)                                    # (N, 3, 2h+1, n)
            num = np.einsum("ackn,cn->ak", V3, w)
            den = np.sqrt(np.einsum("ackn,ackn->ak", V3, V3)) * wn
            ww = 3.0 * wgt.get(k, 1.0)
            acc[idx[keep]] += ww * num / np.maximum(den, 1e-20); wsum[idx[keep]] += ww
    ok = wsum > 0
    cc = np.where(ok[:, None], acc / np.maximum(wsum, 1e-20)[:, None], np.nan)
    j = np.where(ok, np.nanargmax(np.nan_to_num(cc, nan=-9), axis=1), h)
    ev["lag"] = (j - h) / 20.0 + (0.5 if year in TRUNC_YEARS else 0.0)
    ev["net_cc"] = np.where(ok, cc[np.arange(len(ev)), j], np.nan)
    ev["t_origin"] = t0 + (j - h) / 20.0
    ev["retimed"] = ok
    return ev
