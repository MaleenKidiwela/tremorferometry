#!/usr/bin/env python
"""P3 STEP 1a audit — is the WINDOW machinery sound, or is the pair statistic the thing that fails?

Same events, same windows, same MAD normalisation, same site sum, same max-over-lag.  Two
statistics:
  (i)  EVENT vs EVENT   — what q_pairlink measures (template-free linking);
  (ii) EVENT vs the family's own cc-aligned STACK of its other events (leave-one-out)  — what a
       matched filter measures.

If (ii) separates and (i) does not, the data, the time base and the code are fine and the pair
statistic is simply too weak: single-event-to-single-event cc costs sqrt(signal fraction) relative
to event-to-stack cc, which is the whole reason template matching exists.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P
import q_pairlink as Q
import g_common as G

dev = "cuda"


def main(doy=251, win_s=8.0, lag_s=3.0, nfam=12):
    B = Q.candidates()
    B = B[B.doy == doy].reset_index(drop=True)
    top = B.bkey.value_counts().head(nfam).index.tolist()
    B = B[B.bkey.isin(top)].reset_index(drop=True)
    N = len(B)
    bk = B.bkey.to_numpy()
    print(f"[audit] doy {doy}  {N} events  {len(top)} families  win {win_s}s")
    store = P.DayStore()
    ftt = P.fam_tt(Q.STAS)
    nw, nl = int(win_s * P.FS), int(lag_s * P.FS)
    nb = nw + 2 * nl
    accP = torch.zeros(N, N, 2 * nl + 1, device=dev)          # event vs event
    accS = torch.zeros(N, len(top), 2 * nl + 1, device=dev)   # event vs family stack (LOO)
    for s in Q.STAS:
        W = np.zeros((N, 3, nb), np.float32)
        k = store.idx(s, 2012, int(doy)); n = int(store.nvalid[k])
        X = np.asarray(store.get_i(k)[:, :n], np.float32); t0 = float(store.t0[k])
        for row, (t, b) in enumerate(zip(B.t_origin.to_numpy(), bk)):
            a = int(round((t + ftt[(int(b), s)] - lag_s - t0) * P.FS))
            if 0 <= a and a + nb <= n:
                W[row] = X[:, a:a + nb]
        Wd = torch.from_numpy(W).to(dev)
        A = Wd[:, :, nl:nl + nw].reshape(N, 3 * nw)
        A = A / A.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
        U = Wd.unfold(2, nw, 1).permute(0, 2, 1, 3).reshape(N * (2 * nl + 1), 3 * nw)
        U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
        r = (A @ U.T).float()
        v = r.reshape(-1)[torch.randint(0, r.numel(), (400_000,), device=dev)]
        m = v.median(); md = (v - m).abs().median().clamp_min(1e-6) * 1.4826
        accP += ((r - m) / (md * Q.NMEM[Q.SITE_OF[s]])).reshape(N, N, 2 * nl + 1)
        # --- family stacks (cc-aligned, leave-one-out) on the same windows
        Sk = np.zeros((len(top), 3 * nw), np.float32)
        Wn = Wd[:, :, nl:nl + nw].cpu().numpy()
        for fi, f in enumerate(top):
            sel = np.flatnonzero(bk == f)
            Ssub = Wn[sel]
            rms = np.sqrt((Ssub ** 2).mean((1, 2), keepdims=True))
            Ssub = Ssub / np.maximum(rms, 1e-12)
            st, _, _ = G.align_stack(Ssub.astype(np.float64), iters=2, lag=20)
            Sk[fi] = st.reshape(-1)
        Sd = torch.from_numpy(Sk).to(dev)
        Sd = Sd / Sd.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
        rs = (Sd @ U.T).float()
        accS += ((rs - m) / (md * Q.NMEM[Q.SITE_OF[s]])).reshape(len(top), N,
                                                                 2 * nl + 1).permute(1, 0, 2)
        del Wd, A, U, r, rs
        torch.cuda.empty_cache()
    Mp = accP.amax(-1).cpu().numpy()
    Ms = accS.amax(-1).cpu().numpy()
    iu = np.triu_indices(N, 1)
    tt0 = B.t_origin.to_numpy()
    keep = np.abs(tt0[iu[0]] - tt0[iu[1]]) >= 27.0
    same = (bk[:, None] == bk[None, :])[iu][keep]
    v = Mp[iu][keep]
    own = np.array([list(top).index(b) for b in bk])
    s_own = Ms[np.arange(N), own]
    s_oth = np.array([Ms[i, [j for j in range(len(top)) if j != own[i]]].max() for i in range(N)])
    rows = dict(
        pair_same_med=float(np.median(v[same])), pair_cross_med=float(np.median(v[~same])),
        pair_same_p99=float(np.quantile(v[same], .99)),
        pair_cross_p999=float(np.quantile(v[~same], .999)),
        pair_n_same=int(same.sum()), pair_n_cross=int((~same).sum()),
        pair_cross_p9999=float(np.quantile(v[~same], .9999)),
        pair_enrich_at_cross_p9999=float((v[same] > np.quantile(v[~same], .9999)).mean() / 1e-4),
        stack_own_med=float(np.median(s_own)), stack_other_med=float(np.median(s_oth)),
        stack_own_p10=float(np.quantile(s_own, .1)),
        stack_frac_own_above_other_p999=float((s_own > np.quantile(s_oth, .999)).mean()))
    for k, val in rows.items():
        print(f"  {k:42s} {val:.3f}")
    pd.DataFrame([rows]).to_csv(f"{P.RES}/p3/pairlink_stack_audit.csv", index=False)


if __name__ == "__main__":
    main()
