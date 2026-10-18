#!/usr/bin/env python
"""Per-cell +1 % recovery table for a daily-inversion run (RESDIR env; default production res_catalog).
T_cell(t) = self-recovery of a unit injection in cell j through the production per-day estimator day_solve
(joint [m, s], run's lam_f/lam_s, c = 0). day_solve is linear in d, so all cells are injected at once:
R_t = Minv[:Mf] . Aw^T . (w * Gc[rows])  ->  column j = model recovered from injecting 1 in cell j.
Days: mid-June of every year 2011-2026 (gate-passing nearest day). T_cell = median over days of diag(R_t).
Outputs (RESDIR): cell_recovery.csv (cell, lat, lon, depth, T_cell, T_min, T_max, ratio_spread, crosstalk_sum),
                  cell_recovery.npz (R [days x M x M] float32, days, cell).
NOT the quantity the interface explorer applies: the explorer's per-cell T is the pattern-weighted mean of the 20 PATCH
transmissions (125-km Gaussian blobs, transmission_T.npz) -> cell_T_patchmapped.csv. Single-cell self-recovery is much
smaller (production median 0.04, G_diff 0.77) because a 0.1-deg cell is far below the resolution scale; diagnostic only.
"""
import os, sys, numpy as np, pandas as pd
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_v, "2")
sys.path.insert(0, "fault_tomography/inversion"); import invert_dvv_4d_daily as inv
from numpy.linalg import inv as npinv
from concurrent.futures import ProcessPoolExecutor
RD = inv.OUT; z = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); LAMF, LAMS = float(z["lam_f"]), float(z["lam_s"])
days = inv.days; Mf = inv.Mf
assert list(z["cell"]) == list(inv.cellids), "cell order mismatch between npz and G.npz"
rep = sorted({int(days.get_indexer([pd.Timestamp(f"{y}-06-15")], method="nearest")[0]) for y in range(2011, 2027)})

def recover(t):
    idx = inv.BY_DAY[t]; rows = inv.ROW[idx]; w = np.sqrt(inv.Wr[idx]); stj = inv.STJ[idx]; ust = np.unique(stj)
    colmap = {j: k for k, j in enumerate(ust)}; ns = len(ust)
    S = np.zeros((len(idx), ns)); S[np.arange(len(idx)), [colmap[j] for j in stj]] = 1.0
    A = np.hstack([inv.Gc[rows], S]); Aw = A * w[:, None]; AtA = Aw.T @ Aw
    reg = np.zeros_like(AtA); reg[:Mf, :Mf] = LAMF ** 2 * inv.LtL; reg[np.arange(Mf), np.arange(Mf)] += 1e-4
    si = np.arange(Mf, Mf + ns); reg[si, si] += LAMS ** 2
    Minv = npinv(AtA + reg)
    R = Minv[:Mf] @ (Aw.T @ (w[:, None] * inv.Gc[rows]))          # (Mf x Mf): column j = recovered model for unit in cell j
    return t, R.astype(np.float32)

if __name__ == "__main__":
    with ProcessPoolExecutor(min(16, len(rep))) as ex: res = dict(ex.map(recover, rep))
    R = np.stack([res[t] for t in rep]); diag = np.array([np.diag(r) for r in R])
    T = np.median(diag, 0); lo, hi = diag.min(0), diag.max(0); spread = (hi - lo) / np.maximum(np.abs(T), 1e-6)
    leak = np.array([np.abs(r - np.diag(np.diag(r))).sum(0) for r in R]).mean(0)
    tab = pd.DataFrame(dict(cell=list(z["cell"]), lat=np.asarray(z["clat"], float), lon=np.asarray(z["clon"], float), depth_km=np.asarray(z["depth"], float),
                            T_cell=T, T_min=lo, T_max=hi, ratio_spread=spread, crosstalk_sum=leak))
    tab.to_csv(f"{RD}/cell_recovery.csv", index=False)
    np.savez_compressed(f"{RD}/cell_recovery.npz", R=R, days=days[rep].astype(str).values, cell=np.array(list(z["cell"]), object), lam_f=LAMF, lam_s=LAMS)
    dep = tab.depth_km.values; deep = dep >= 32
    print(f"{RD}: {Mf} cells, {len(rep)} days | T_cell median {np.median(T):.3f} (p10-p90 {np.percentile(T,10):.2f}-{np.percentile(T,90):.2f}); deep {np.median(T[deep]):.3f}, up-dip {np.median(T[~deep]):.3f} | time spread median {np.median(spread):.3f} | cells T<0.2: {(T<0.2).sum()}")
    print(f"wrote {RD}/cell_recovery.csv, cell_recovery.npz")
