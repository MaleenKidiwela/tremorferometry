#!/usr/bin/env python
"""Shared helpers for agent S's FINAL assessment chain.  Prefix `sf_`.

Nothing here is a screen clause.  Two things live in this file:

  * `dedupe()` -- the ONE family-acceptance rule that is positively attested in the published
    Cascadia LFE literature (results/SCREEN_VALIDATION_NORTH.md sec.4):
      - Bostock et al. 2012 (via Royer, Thomas & Bostock 2015): 130 templates "culled to 93
        geographically independent LFE families";
      - Plourde et al. 2015: "the templates were considered to originate from the same family if
        > 10 % of detections were shared ... A correlation coefficient > 0.3 was considered to
        suggest templates originated from the same family";
      - Thomas & Bostock 2015: "no family shares more than 10 % of detection times";
      - Chestler & Creager 2017: R_com > 0.5, keep the family with the higher median cc.
    Implemented literally: two templates are the SAME family if their declustered hit sets share
    >= SHARE_MIN of the smaller set's detections OR their stacks cross-correlate >= CC_MIN over
    the stations they have in common; the survivor is the one with more detections (tie-break:
    higher median cc), which is Chestler's "nearly always the family with fewer detections" rule
    written the way she describes its effect.
  * small geometry / template I/O helpers.

No k_* file is edited or imported for anything but its published behaviour.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
OUT = f"{RES}/sf_final"
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)

SHARE_MIN = 0.10        # Plourde 2015 / Thomas & Bostock 2015: > 10 % of detections shared
CC_MIN = 0.30           # Plourde 2015: cc > 0.3 suggests the same family
DECL_S = 4.0            # Bostock 2015's declustering interval -- the shared-time bin
MAX_LAG_S = 2.0         # lag search for the stack cross-correlation


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


def ensure_out():
    os.makedirs(OUT, exist_ok=True)
    return OUT


# ---------------------------------------------------------------- template I/O ----------------
def load_tpl(path, grps=None):
    """{g: {sta: (tpl (3,n) float32, off_s, peak_mad)}} exactly as k_grow saved it."""
    z = np.load(path, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k, v in meta.items():
        g, s = k.split("_", 1)
        g = int(g)
        if grps is not None and g not in grps:
            continue
        TPL.setdefault(g, {})[s] = (np.asarray(z[k], np.float32), float(v[0]), float(v[1]))
    return TPL


def save_tpl(path, TPL):
    np.savez_compressed(
        path,
        **{f"{g}_{s}": TPL[g][s][0].astype(np.float32) for g in TPL for s in TPL[g]},
        meta=np.array(json.dumps({f"{g}_{s}": [float(TPL[g][s][1]), float(TPL[g][s][2])]
                                  for g in TPL for s in TPL[g]})))


# ---------------------------------------------------------------- geometry --------------------
def epi_km(la1, lo1, la2, lo2):
    import s_common as S
    return np.asarray(S.epi_km(la1, lo1, la2, lo2), float)


def hyp_km(la1, lo1, dp1, la2, lo2, dp2):
    d = epi_km(la1, lo1, la2, lo2)
    return np.sqrt(d ** 2 + (np.asarray(dp1, float) - np.asarray(dp2, float)) ** 2)


# ---------------------------------------------------------------- shared detections -----------
def shared_counts(det, grps):
    """(n, n) count of shared declustered detection times, and the per-group detection count.

    Times are binned at the 4-s declustering interval in TWO bin systems offset by half a bin, and
    the larger of the two counts is used, so any two detections within 2 s are counted as shared
    in at least one system (and nothing further apart than 4 s ever is)."""
    from scipy import sparse
    gi = {g: i for i, g in enumerate(grps)}
    d = det[det.grp.isin(gi)]
    rows = np.array([gi[int(g)] for g in d.grp], np.int64)
    t = d.t.to_numpy(np.float64)
    n = len(grps)
    S = None
    for off in (0.0, DECL_S / 2):
        b = np.floor((t + off) / DECL_S).astype(np.int64)
        ub, bi = np.unique(b, return_inverse=True)
        M = sparse.csr_matrix((np.ones(len(rows), np.float32), (rows, bi)),
                              shape=(n, len(ub)))
        M.data[:] = 1.0
        M.sum_duplicates()
        M.data[:] = 1.0                       # binary: a group is in a bin or it is not
        C = (M @ M.T).toarray()
        S = C if S is None else np.maximum(S, C)
    ndet = np.array([int((d.grp == g).sum()) for g in grps], np.int64)
    return S.astype(np.int64), ndet


# ---------------------------------------------------------------- template cc -----------------
def stack_cc_matrix(TPL, grps, fs=20.0, dev=None, chunk=None):
    """(n, n) network-mean normalised cross-correlation of the stacks, maximised over +-2 s.

    Per station the two 3-component templates are normalised jointly over the compared window and
    correlated; the pair value is the mean over the stations the two templates have in common
    (Plourde 2015 correlates template waveforms; the lag search absorbs the different stack
    onsets).  Runs on the GPU when one is available, otherwise on the CPU."""
    import torch
    if dev is None:
        dev = "cuda" if torch.cuda.is_available() else "cpu"
    n = len(grps)
    stas = sorted({s for g in grps for s in TPL.get(g, {})})
    nt = min(TPL[g][s][0].shape[1] for g in grps for s in TPL[g])
    ml = int(round(MAX_LAG_S * fs))
    c0, c1 = ml, nt - ml
    if c1 - c0 < int(fs):
        ml = max(int(nt // 6), 1)
        c0, c1 = ml, nt - ml
    acc_best = torch.full((n, n), -2.0, device=dev)
    X, HAVE = {}, {}
    for s in stas:
        A = np.zeros((n, 3, nt), np.float32)
        h = np.zeros(n, bool)
        for i, g in enumerate(grps):
            v = TPL.get(g, {}).get(s)
            if v is not None:
                A[i] = v[0][:, :nt]
                h[i] = True
        X[s] = torch.from_numpy(A).to(dev)
        HAVE[s] = torch.from_numpy(h).to(dev)

    def unit(t):
        v = t.reshape(t.shape[0], -1)
        return v / v.norm(dim=1, keepdim=True).clamp_min(1e-20)

    for k in range(-ml, ml + 1):
        num = torch.zeros((n, n), device=dev)
        cnt = torch.zeros((n, n), device=dev)
        for s in stas:
            A = unit(X[s][:, :, c0:c1])
            B = unit(X[s][:, :, c0 + k:c1 + k])
            m = (HAVE[s][:, None] & HAVE[s][None, :]).float()
            num += (A @ B.T) * m
            cnt += m
        cc = num / cnt.clamp_min(1.0)
        cc = torch.where(cnt > 0, cc, torch.full_like(cc, -2.0))
        acc_best = torch.maximum(acc_best, cc)
    C = acc_best.cpu().numpy()
    C = np.maximum(C, C.T)                    # the lag search is one-sided per entry
    np.fill_diagonal(C, 1.0)
    return C


# ---------------------------------------------------------------- the dedupe ------------------
def dedupe(det, TPL, grps=None, share_min=SHARE_MIN, cc_min=CC_MIN, do_cc=True):
    """The literature dedupe.  Returns (keep_list, DataFrame of every template's fate, pairs)."""
    grps = sorted(TPL) if grps is None else sorted(grps)
    S, ndet = shared_counts(det, grps)
    med_cc = det.groupby("grp").cc.median()
    n = len(grps)
    mn = np.minimum(ndet[:, None], ndet[None, :])
    share = S / np.maximum(mn, 1)
    np.fill_diagonal(share, 0.0)
    if do_cc and n > 1:
        C = stack_cc_matrix(TPL, grps)
        np.fill_diagonal(C, 0.0)
    else:
        C = np.zeros((n, n), np.float32)
    same = (share >= share_min) | (C >= cc_min)
    np.fill_diagonal(same, False)
    order = sorted(range(n), key=lambda i: (-int(ndet[i]),
                                            -float(med_cc.get(grps[i], 0.0))))
    keep, absorbed = [], {}
    for i in order:
        hit = None
        for j in keep:
            if same[i, j]:
                hit = j
                break
        if hit is None:
            keep.append(i)
        else:
            absorbed[i] = hit
    rows = []
    for i in range(n):
        rows.append(dict(grp=int(grps[i]), n_det=int(ndet[i]),
                         med_cc=float(med_cc.get(grps[i], np.nan)),
                         kept=bool(i in keep),
                         merged_into=(int(grps[absorbed[i]]) if i in absorbed else -1),
                         n_partners=int(same[i].sum()),
                         max_share=float(share[i].max()) if n > 1 else 0.0,
                         max_cc=float(C[i].max()) if n > 1 else 0.0))
    D = pd.DataFrame(rows).sort_values("grp").reset_index(drop=True)
    pairs = dict(n_templates=n, n_kept=len(keep),
                 n_pairs_same=int(same.sum() // 2),
                 n_pairs_share=int(((share >= share_min).sum() - 0) // 2),
                 n_pairs_cc=int(((C >= cc_min).sum()) // 2),
                 share_min=share_min, cc_min=cc_min,
                 max_cc_overall=float(C.max()) if n > 1 else 0.0,
                 max_share_overall=float(share.max()) if n > 1 else 0.0)
    return [int(grps[i]) for i in keep], D, pairs
