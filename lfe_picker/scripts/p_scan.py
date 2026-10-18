#!/usr/bin/env python
"""Continuous per-station characteristic function over whole station-days, + peak declustering.

The scan reproduces the TRAINING input distribution exactly: overlapping windows of the model's own
length, each peak-normalised over its three components, with only the central third of each window's
prediction kept (so every output sample comes from a window in which it sat away from the edge).
"""
from __future__ import annotations
import os, sys
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import p_models as M

FS = P.FS


def load_model(path, dev="cuda"):
    z = torch.load(path, map_location=dev, weights_only=False)
    m = M.LinUNet(z["fac"]) if z["arch"] == "lin" else M.PhaseNetS()
    m.load_state_dict(z["state"])
    m.to(dev).eval()
    return m, z


@torch.no_grad()
def cf_day(model, meta, x, dev="cuda", batch=1024):
    """x (3, N) float16/float32 -> CF (N,) float32 in [0,1].  Zero where the data are invalid."""
    nw = int(round(meta["win_s"] * FS))
    step = nw // 3
    N = x.shape[1]
    starts = np.arange(0, max(1, N - nw + 1), step, dtype=np.int64)
    if starts[-1] != N - nw:
        starts = np.append(starts, N - nw)
    out = np.zeros(N, np.float32)
    cov = np.zeros(N, np.bool_)
    X = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
    for i0 in range(0, len(starts), batch):
        sl = starts[i0:i0 + batch]
        idx = torch.from_numpy(sl).to(dev)[:, None] + torch.arange(nw, device=dev)[None, :]
        w = X[:, idx.reshape(-1)].reshape(3, len(sl), nw).permute(1, 0, 2).float()
        m = w.abs().amax(dim=(1, 2), keepdim=True).clamp_min(1e-12)
        w = w / m
        with torch.amp.autocast("cuda", dtype=torch.bfloat16):
            if meta["arch"] == "lin":
                p = model(M.zen2inp(M.upsample2(w))).float()
                p = p.reshape(len(sl), nw, 2).amax(-1)          # 40 Hz -> back to 20 Hz
            else:
                p = model(w).float()
        p = p.cpu().numpy()
        a, b = nw // 3, 2 * (nw // 3)
        for k, s in enumerate(sl):
            lo, hi = (0 if s == 0 else a), (nw if s == starts[-1] else b)
            out[s + lo:s + hi] = np.maximum(out[s + lo:s + hi], p[k, lo:hi])
            cov[s + lo:s + hi] = True
    out[~cov] = 0.0
    dead = (np.asarray(x, np.float32) == 0).all(0)
    out[dead] = 0.0
    return out


def decluster(v, sep_s=2.0, keep=20000):
    """(times_idx, values) of local maxima at least sep_s apart, the strongest `keep` kept."""
    n = int(round(sep_s * FS))
    t = torch.from_numpy(v)[None, None]
    mx = torch.nn.functional.max_pool1d(t, kernel_size=2 * n + 1, stride=1, padding=n)[0, 0].numpy()
    m = (v >= mx - 1e-9) & (v > 0)
    idx = np.flatnonzero(m)
    if len(idx) > keep:
        o = np.argsort(v[idx])[::-1][:keep]
        idx = np.sort(idx[o])
    return idx, v[idx]


def one2one(pred_t, obs_t, tol=1.0):
    """greedy one-to-one match of observed arrivals to predicted peaks; returns matched mask."""
    if len(pred_t) == 0 or len(obs_t) == 0:
        return np.zeros(len(obs_t), bool)
    pred_t = np.sort(pred_t)
    used = np.zeros(len(pred_t), bool)
    hit = np.zeros(len(obs_t), bool)
    j = np.searchsorted(pred_t, obs_t)
    for i, o in enumerate(obs_t):
        for k in range(max(0, j[i] - 2), min(len(pred_t), j[i] + 2)):
            if not used[k] and abs(pred_t[k] - o) <= tol:
                used[k] = True
                hit[i] = True
                break
    return hit
