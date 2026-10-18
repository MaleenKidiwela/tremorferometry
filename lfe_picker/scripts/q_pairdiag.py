#!/usr/bin/env python
"""P3 STEP 1a — diagnostic: is the pair statistic working at all, and how much of it does the
template-free window placement cost?

Compares, on one dev day, the SAME-FAMILY pair cc against the cross-family pair cc under four
window configurations, at the PGC site and summed over the 7 discovery sites:

  A  24 s at the grid-centre constant   (exactly what q_pairlink does — template-free)
  B  24 s at the family's OWN predicted S  (oracle placement, still 24 s)
  C   8 s at the family's OWN predicted S  (oracle placement + oracle duration)
  D   4 s at the family's OWN predicted S

If A is flat but C/D separate, the kill test is about window DILUTION, not about the machinery.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import q_pairlink as Q

dev = "cuda"


def pairs_cc(store, B, tt, win_s, lag_s=3.0):
    """network pair cc (max over lag of the site-summed MAD-z), (n, n) numpy."""
    nw, nl = int(win_s * P.FS), int(lag_s * P.FS)
    nb = nw + 2 * nl
    N = len(B)
    accs = None
    for si, s in enumerate(Q.STAS):
        W = np.zeros((N, 3, nb), np.float32)
        ok = np.zeros(N, bool)
        for doy, g in B.groupby("doy"):
            if not store.has(s, 2012, int(doy)):
                continue
            k = store.idx(s, 2012, int(doy))
            n = int(store.nvalid[k])
            X = np.asarray(store.get_i(k)[:, :n], np.float32)
            t0 = float(store.t0[k])
            for row, t, b in zip(g.index.to_numpy(), g.t_origin.to_numpy(), g.bkey.to_numpy()):
                a = int(round((t + tt(int(b), s) - lag_s - t0) * P.FS))
                if a < 0 or a + nb > n:
                    continue
                W[row, :, :] = X[:, a:a + nb]
                ok[row] = True
        Wd = torch.from_numpy(W).to(dev)
        A = Wd[:, :, nl:nl + nw].reshape(N, 3 * nw)
        A = A / A.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
        U = Wd.unfold(2, nw, 1).permute(0, 2, 1, 3).reshape(N * (2 * nl + 1), 3 * nw)
        U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
        r = (A @ U.T).float()
        v = r.reshape(-1)
        v = v[torch.randint(0, v.numel(), (min(400_000, v.numel()),), device=dev)]
        m = v.median(); md = (v - m).abs().median().clamp_min(1e-6) * 1.4826
        z = ((r - m) / (md * Q.NMEM[Q.SITE_OF[s]])).reshape(N, N, 2 * nl + 1)
        accs = z if accs is None else accs + z
        del Wd, A, U, r, z
        torch.cuda.empty_cache()
    return accs.amax(-1).cpu().numpy(), ok


def main():
    B = Q.candidates()
    doy = 251
    B = B[B.doy == doy].reset_index(drop=True)
    # keep the 12 most prolific families of the day, to bound n^2
    top = B.bkey.value_counts().head(12).index.tolist()
    B = B[B.bkey.isin(top)].reset_index(drop=True)
    print(f"[diag] doy {doy}: {len(B)} events, {B.bkey.nunique()} families")
    store = P.DayStore()
    cst, _ = Q.station_constants()
    ftt = P.fam_tt(Q.STAS)
    conf = [("A grid-centre 24 s", 24.0, lambda b, s: cst[s]),
            ("B own-cell 24 s", 24.0, lambda b, s: ftt[(b, s)]),
            ("C own-cell 8 s", 8.0, lambda b, s: ftt[(b, s)]),
            ("D own-cell 4 s", 4.0, lambda b, s: ftt[(b, s)])]
    bk = B.bkey.to_numpy()
    tt0 = B.t_origin.to_numpy()
    same = bk[:, None] == bk[None, :]
    iu = np.triu_indices(len(B), 1)
    # pairs closer in time than window+lag SHARE SAMPLES -> cc ~ 1 trivially (same rule as q_pairlink)
    keep = np.abs(tt0[iu[0]] - tt0[iu[1]]) >= 27.0
    rows = []
    for name, ws, tt in conf:
        M, ok = pairs_cc(store, B, tt, ws)
        v = M[iu][keep]
        sm = same[iu][keep]
        rows.append(dict(config=name, n_same=int(sm.sum()), n_cross=int((~sm).sum()),
                         same_med=float(np.median(v[sm])), cross_med=float(np.median(v[~sm])),
                         same_p99=float(np.quantile(v[sm], .99)),
                         cross_p99=float(np.quantile(v[~sm], .99)),
                         same_max=float(v[sm].max()), cross_max=float(v[~sm].max()),
                         cross_p9999=float(np.quantile(v[~sm], .9999)),
                         n_same_above_cross_p9999=int(
                             (v[sm] > np.quantile(v[~sm], .9999)).sum()),
                         enrich_at_cross_p9999=float(
                             (v[sm] > np.quantile(v[~sm], .9999)).mean() / 1e-4)))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(f"{P.RES}/p3/pairlink_window_diag.csv", index=False)


if __name__ == "__main__":
    main()
