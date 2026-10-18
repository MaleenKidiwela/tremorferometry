#!/usr/bin/env python
"""P4-0b template banks.

  pass0      20-s templates from predicted S - 6 s at the seed's moveout (Merlin P4-2 step 3)
  nmoveout   the SAME seeds re-cut with a WRONG cell's per-station timing (predicted S differing
             >= 5 s at >= half the stations) -- the per-template moveout-scramble null
  ref        the arm-(i) REFERENCE: each eligible kept family's own 12-s truth-stack template,
             built by align_stack over ALL its hp events at each site and cut from S - 2 s, run
             through the identical threshold/growth machinery.

  DOCUMENTED DEVIATION.  Merlin's arm-(i) reference is "the family's own transferred template".
  The P3D transfer bank on disk (results/p3d/bank_transfer_v1*/bank_production.npz) covers only
  26-184 (family, station) pairs and NONE at PGC/VGZ/ALBHB, so it cannot supply a 7-site reference
  for the 18 eligible families.  The truth-stack template uses all of the family's hp events and is
  therefore AT LEAST as strong as a transferred 2012 template: the relative clauses
  (recall >= 0.5 x ref, purity >= 0.8 x ref) are conservative against it.  Absolute recall/purity
  are reported alongside so the bar can be re-read.

  python scripts/v_bank.py pass0 --seeds results/p4/p4_0b/seeds_given.parquet --out .../bank_given0.npz
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

OUT = V.OUT


def _cut(sta, t_abs, n, idx):
    """(3, n) float32 window at an absolute start time, each channel normalised to unit L2."""
    ts = pd.to_datetime(float(t_abs), unit="s")
    v = S.load_day_zne(sta, int(ts.year), int(ts.dayofyear), idx)
    if v is None:
        return None
    t0, X = v
    i0 = int(round((t_abs - t0) * V.FS))
    if i0 < 0 or i0 + n > X.shape[1]:
        return None
    seg = X[:, i0:i0 + n].astype(np.float64)
    if not np.isfinite(seg).all() or (seg ** 2).sum() <= 0:
        return None
    if (seg ** 2).sum(1).min() <= 0:
        return None
    return V.norm_channels(seg)


def build_pass0(seeds, moveout, out_path, pre_s=V.TPL0_PRE_S, len_s=V.TPL0_LEN_S, tag=""):
    """moveout: {seed_id: {sta: tt}}.  off = tt - pre_s (template start minus origin)."""
    idx = V.cache_idx()
    stas = V.DET_STAS
    n = int(round(len_s * V.FS))
    B = len(seeds)
    W = np.full((B, len(stas), 3, n), np.nan, np.float32)
    off = np.full((B, len(stas)), np.nan)
    tids = []
    for bi, r in enumerate(seeds.itertuples()):
        mv = moveout[r.seed]
        tids.append(r.seed)
        for si, s in enumerate(stas):
            o = float(mv[s]) - pre_s
            off[bi, si] = o
            w = _cut(s, float(r.t) + o, n, idx)
            if w is not None:
                W[bi, si] = w
    ok = np.isfinite(W).all(axis=(2, 3))
    nsite = np.array([len({V.SITE_OF[stas[i]] for i in np.nonzero(ok[b])[0]}) for b in range(B)])
    meta = dict(kind=f"pass0{tag}", len_s=len_s, pre_s=pre_s, n_seeds=B,
                n_sta_median=float(np.median(ok.sum(1))), n_site_median=float(np.median(nsite)),
                n_with_ge3_sites=int((nsite >= V.MIN_SITES).sum()))
    keep = np.nonzero(nsite >= V.MIN_SITES)[0]
    V.save_bank(out_path, [tids[i] for i in keep], stas, W[keep], off[keep], meta)
    V.log(f"[bank{tag}] {len(keep)}/{B} seeds with >= {V.MIN_SITES} sites -> {out_path}")
    return out_path


def wrong_cells():
    """Coarse cell grid over the north box at the family-field depth, for the moveout scramble."""
    lat0, lat1 = 48.0, 49.0
    lon0, lon1 = -124.5, -123.0
    _, _, dep, _ = V.region_centre()
    cells = []
    for la in np.arange(lat0, lat1 + 1e-9, 0.1):
        for lo in np.arange(lon0, lon1 + 1e-9, 0.15):
            for dz in (-8.0, 0.0, 8.0):
                cells.append((float(la), float(lo), float(np.clip(dep + dz, 15, 55))))
    return cells


def build_nmoveout(seeds, moveout, out_path, seed_rng=20260920):
    """Per template: a WRONG cell whose predicted S differs >= 5 s at >= half the stations."""
    ST = V.station_table()
    cells = wrong_cells()
    TT = [V.cell_tt(la, lo, dz, ST=ST) for la, lo, dz in cells]
    stas = V.DET_STAS
    mv2, chosen = {}, []
    for r in seeds.itertuples():
        true = np.array([moveout[r.seed][s] for s in stas])
        h = int(hashlib.sha256(f"{seed_rng}:{r.seed}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(h)
        cand = [k for k, tt in enumerate(TT)
                if (np.abs(np.array([tt[s] for s in stas]) - true) >= V.NMOVEOUT_MIN_DT).sum()
                >= int(np.ceil(len(stas) / 2))]
        if not cand:
            continue
        k = int(rng.choice(cand))
        mv2[r.seed] = TT[k]
        d = np.abs(np.array([TT[k][s] for s in stas]) - true)
        chosen.append(dict(seed=r.seed, cell_lat=cells[k][0], cell_lon=cells[k][1],
                           cell_dep=cells[k][2], n_sta_ge5s=int((d >= V.NMOVEOUT_MIN_DT).sum()),
                           median_dt=float(np.median(d))))
    sel = seeds[seeds.seed.isin(mv2)].copy()
    pd.DataFrame(chosen).to_csv(out_path.replace(".npz", "_cells.csv"), index=False)
    return build_pass0(sel, mv2, out_path, tag="_nmoveout")


def build_ref(out_path):
    """Truth-stack 12-s templates for the eligible kept families (the arm-(i) reference)."""
    keep, D, _ = V.eligible_families()
    idx = V.cache_idx()
    stas = V.DET_STAS
    FT = V.family_tt(stas)
    n_cut = int(round(24.0 * V.FS))
    pre_cut = 8.0                       # cut S-8 .. S+16 for the alignment, then window 12 s
    n = int(round(V.TPL1_LEN_S * V.FS))
    W = np.full((len(keep), len(stas), 3, n), np.nan, np.float32)
    off = np.full((len(keep), len(stas)), np.nan)
    diag = []
    for bi, b in enumerate(keep):
        g = D[D.bkey == b]
        for si, s in enumerate(stas):
            tt = FT[(b, s)]
            segs = []
            for t in g.t.to_numpy(float):
                w = _cut(s, t + tt - pre_cut, n_cut, idx)
                if w is not None:
                    segs.append(w)
            if len(segs) < 10:
                continue
            st, _, _ = S.align_stack(np.stack(segs), iters=V.ALIGN_ITERS, lag=V.ALIGN_LAG)
            i_s = S.pick_s(st, int(round(pre_cut * V.FS)), search_s=4.0)
            a = i_s - int(round(V.TPL1_PRE_S * V.FS))
            if a < 0 or a + n > st.shape[1]:
                continue
            W[bi, si] = V.norm_channels(st[:, a:a + n])
            off[bi, si] = tt + (i_s - int(round(pre_cut * V.FS))) / V.FS - V.TPL1_PRE_S
            diag.append(dict(bkey=int(b), sta=s, n_stacked=len(segs),
                             s_snr_db=float(S.s_snr_db(st, i_s)),
                             s_excess_db=float(S.s_excess_db(st, i_s)),
                             pick_shift_s=float((i_s - int(round(pre_cut * V.FS))) / V.FS)))
    ok = np.isfinite(W).all(axis=(2, 3))
    nsite = np.array([len({V.SITE_OF[stas[i]] for i in np.nonzero(ok[b])[0]})
                      for b in range(len(keep))])
    sel = np.nonzero(nsite >= V.MIN_SITES)[0]
    pd.DataFrame(diag).to_csv(out_path.replace(".npz", "_diag.csv"), index=False)
    V.save_bank(out_path, [f"REF{keep[i]}" for i in sel], stas, W[sel], off[sel],
                dict(kind="truth_stack_reference", len_s=V.TPL1_LEN_S, pre_s=V.TPL1_PRE_S,
                     families=[int(keep[i]) for i in sel],
                     note="align_stack over ALL hp events per family/site; S re-picked per site"))
    V.log(f"[ref] {len(sel)}/{len(keep)} families -> {out_path}")
    return out_path


def moveout_for(seeds):
    """{seed_id: {sta: predicted S delay}} -- family cell for arm (i), region centre for arm (ii)."""
    stas = V.DET_STAS
    if "bkey" in seeds.columns:
        FT = V.family_tt(stas)
        return {r.seed: {s: FT[(int(r.bkey), s)] for s in stas} for r in seeds.itertuples()}
    lat, lon, dep, _ = V.region_centre()
    cst = V.cell_tt(lat, lon, dep)
    return {r.seed: cst for r in seeds.itertuples()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["pass0", "nmoveout", "ref"])
    ap.add_argument("--seeds", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    if os.path.exists(a.out):
        V.log(f"[bank] {a.out} exists -- skip")
        return
    if a.cmd == "ref":
        build_ref(a.out)
        return
    sd = pd.read_parquet(a.seeds)
    if a.limit:
        sd = sd.head(a.limit)
    mv = moveout_for(sd)
    if a.cmd == "pass0":
        build_pass0(sd, mv, a.out)
    else:
        build_nmoveout(sd, mv, a.out)


if __name__ == "__main__":
    main()
