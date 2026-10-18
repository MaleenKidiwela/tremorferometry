#!/usr/bin/env python
"""P4-0b iteration and screen (Merlin P4-2 steps 5 and 6).

  iterate  align_stack (2 iterations, +-1.5 s -- the frozen recipe in s_common.align_stack) on the
           INDEPENDENT hits at each detection site -> re-pick S on the aligned stack at the
           HIGHEST-S-SNR site -> that single global shift is applied to every site -> re-cut 12-s
           templates from S - 2 s.  Exactly two iterations are run by the queue; no third.
  screen   CHEAP SCREEN, NOT EVIDENCE (Merlin F3): alignment-free mean stack of the independent
           hits at each site, S excess in dB (s_common.s_excess_db).  Reported, never used as
           acceptance.
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
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402
import v_bank as B                                                          # noqa: E402

CUT_PRE_S = 8.0
CUT_LEN_S = 24.0
MAX_STACK = 400


def _cut_hits(sta, times, off_s, idx, n, pre):
    segs = []
    for t in times:
        w = B._cut(sta, float(t) + off_s, n, idx)
        if w is not None:
            segs.append(w)
    return np.stack(segs) if segs else None


def iterate(bank, hits_indep, survivors, out_bank, out_diag, use_all=False, min_indep=10):
    """use_all: take every template with >= min_indep independent hits instead of the SURVIVORS.
    The pre-registered survival clauses are applied ONCE, at pass 0; iterations 1 and 2 simply
    re-stack and re-grow (Merlin P4-2 step 5: "iterate exactly twice")."""
    tids, stas, W, off, meta = V.load_bank(bank)
    Hs = pd.read_parquet(hits_indep)
    SV = pd.read_csv(survivors)
    keep = (SV[SV.n_indep >= min_indep].tid.tolist() if use_all else SV[SV.survives].tid.tolist())
    tid2i = {t: i for i, t in enumerate(tids)}
    n_cut = int(round(CUT_LEN_S * V.FS))
    nom = int(round(CUT_PRE_S * V.FS))
    n = int(round(V.TPL1_LEN_S * V.FS))
    idx = V.cache_idx()
    pre_in = float(meta.get("pre_s", V.TPL0_PRE_S))
    Wn = np.full((len(keep), len(stas), 3, n), np.nan, np.float32)
    offn = np.full((len(keep), len(stas)), np.nan)
    diag = []
    tids_out = []
    for bi, tid in enumerate(keep):
        src = tid2i[tid]
        g = Hs[Hs.tid == tid].sort_values("net_mad", ascending=False).head(MAX_STACK)
        if len(g) < 10:
            continue
        times = g.t.to_numpy(float)
        stacks, picks, snr = {}, {}, {}
        for si, s in enumerate(stas):
            o = off[src, si]
            if not np.isfinite(o):
                continue
            segs = _cut_hits(s, times, o + pre_in - CUT_PRE_S, idx, n_cut, CUT_PRE_S)
            if segs is None or len(segs) < 10:
                continue
            st, _, _ = S.align_stack(segs, iters=V.ALIGN_ITERS, lag=V.ALIGN_LAG)
            i_s = S.pick_s(st, nom, search_s=4.0)
            stacks[si] = st
            picks[si] = i_s
            snr[si] = float(S.s_snr_db(st, i_s))
            diag.append(dict(tid=tid, sta=stas[si], n_stacked=len(segs), s_snr_db=snr[si],
                             s_excess_db=float(S.s_excess_db(st, i_s)),
                             pick_shift_s=float((i_s - nom) / V.FS)))
        if not stacks:
            continue
        best = max(snr, key=lambda k: (snr[k] if np.isfinite(snr[k]) else -1e9))
        dt = (picks[best] - nom) / V.FS
        bi2 = len(tids_out)
        got = 0
        for si, st in stacks.items():
            a = nom + int(round(dt * V.FS)) - int(round(V.TPL1_PRE_S * V.FS))
            if a < 0 or a + n > st.shape[1]:
                continue
            Wn[bi2, si] = V.norm_channels(st[:, a:a + n])
            base = off[src, si] + pre_in          # = predicted S delay at this site
            offn[bi2, si] = base + dt - V.TPL1_PRE_S
            got += 1
        if got == 0:
            continue
        tids_out.append(tid)
        diag.append(dict(tid=tid, sta="__global__", n_stacked=len(g),
                         s_snr_db=snr[best], s_excess_db=np.nan, pick_shift_s=float(dt)))
    Wn = Wn[:len(tids_out)]
    offn = offn[:len(tids_out)]
    ok = np.isfinite(Wn).all(axis=(2, 3))
    nsite = np.array([len({V.SITE_OF[stas[i]] for i in np.nonzero(ok[b])[0]})
                      for b in range(len(tids_out))])
    sel = np.nonzero(nsite >= V.MIN_SITES)[0]
    pd.DataFrame(diag).to_csv(out_diag, index=False)
    V.save_bank(out_bank, [tids_out[i] for i in sel], stas, Wn[sel], offn[sel],
                dict(kind="iterated", len_s=V.TPL1_LEN_S, pre_s=V.TPL1_PRE_S,
                     parent=os.path.basename(bank), n_in=len(keep)))
    V.log(f"[iterate] {len(keep)} survivors -> {len(sel)} iterated templates -> {out_bank}")


def screen(bank, hits_indep, out_csv):
    tids, stas, W, off, meta = V.load_bank(bank)
    Hs = pd.read_parquet(hits_indep)
    n_cut = int(round(CUT_LEN_S * V.FS))
    nom = int(round(CUT_PRE_S * V.FS))
    idx = V.cache_idx()
    pre = meta.get("pre_s", V.TPL1_PRE_S)
    rows = []
    for bi, tid in enumerate(tids):
        g = Hs[Hs.tid == tid].sort_values("net_mad", ascending=False).head(MAX_STACK)
        if len(g) == 0:
            continue
        times = g.t.to_numpy(float)
        n_ok = n_reg = 0
        per = {}
        for si, s in enumerate(stas):
            o = off[bi, si]
            if not np.isfinite(o):
                continue
            segs = _cut_hits(s, times, o + pre - CUT_PRE_S, idx, n_cut, CUT_PRE_S)
            if segs is None or len(segs) < 10:
                continue
            st = segs.mean(0)                       # ALIGNMENT-FREE
            ex = float(S.s_excess_db(st, nom))
            per[s] = ex
            n_reg += 1
            n_ok += int(ex >= V.SCREEN_DB)
        rows.append(dict(tid=tid, n_registering=n_reg, n_ge3db=n_ok,
                         frac=n_ok / max(n_reg, 1),
                         passes=bool(n_ok >= 3 or (n_reg and n_ok / n_reg >= 0.5)),
                         **{f"ex_{k}": v for k, v in per.items()}))
    R = pd.DataFrame(rows)
    R.to_csv(out_csv, index=False)
    V.log(f"[screen] {len(R)} templates, {int(R.passes.sum()) if len(R) else 0} pass the SCREEN "
          f"(cheap, not evidence)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["iterate", "screen"])
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--survivors", default="")
    ap.add_argument("--out-bank", default="")
    ap.add_argument("--out-diag", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--use-all", action="store_true")
    a = ap.parse_args()
    if a.cmd == "iterate":
        if os.path.exists(a.out_bank):
            V.log(f"[iterate] {a.out_bank} exists -- skip")
            return
        iterate(a.bank, a.indep, a.survivors, a.out_bank, a.out_diag, use_all=a.use_all)
    else:
        screen(a.bank, a.indep, a.out)


if __name__ == "__main__":
    main()
