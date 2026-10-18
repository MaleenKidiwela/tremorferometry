#!/usr/bin/env python
"""Merlin step-3 acceptance test (I2): on a REDUCED problem (3 months x ~30 best-covered stations),
solve the daily fixed-effects joint inversion two ways and compare:
  (a) BCD exactly as invert_dvv_4d_daily.py does it (per-day solves + closed-form c + per-sweep gauge);
  (b) DIRECT sparse joint solve of the full normal equations (m all-days + s all-station-days + c).
Both gauged identically post-hoc (mean_t m=0 per cell, mean_t s=0 per station-era, absorbed into c).
ACCEPT: max field discrepancy << per-frame sigma_m (report both); monotone BCD objective.
FAILURE signature: objective plateaus while c still walks -> era x site gauge (B2) missing.
"""
import os
import numpy as np, pandas as pd
from numpy.linalg import inv
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

OUT = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
T0, T1 = pd.Timestamp("2021-01-01"), pd.Timestamp("2021-03-31")
NSTA, LAM_F, LAM_RATIO = 30, 0.5, 0.125
SWEEPS = 40

gz = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = gz["G"].astype(np.float64)
fcap = gz["captured"].astype(np.float64)
cellids = gz["cell"]; cxy = gz["cxy"]
keep = fcap > 0
Gc_full = (G * fcap[:, None])[keep]
ptag, pcell = gz["pair_tag"][keep], gz["pair_cell"][keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
cellpos = {c: i for i, c in enumerate(cellids)}

fd = pd.read_parquet(f"{OUT}/family_days.parquet")
fd["date"] = pd.to_datetime(fd.date)
fd = fd[(fd.date >= T0) & (fd.date <= T1)]
top = fd.groupby("tag").size().sort_values(ascending=False).head(NSTA).index
fd = fd[fd.tag.isin(top)].copy()
fd["grow"] = [row_of.get((t, c), -1) for t, c in zip(fd.tag, fd.cell)]
fd = fd[fd.grow >= 0].copy()
fams = pd.read_csv(f"{OUT}/families_daily.csv")
sig_map = dict(zip(fams.tag + "|" + fams.fam + "|e" + fams.era.astype(str), fams.sigma))
fd["fe"] = fd.tag + "|" + fd.fam + "|e" + fd.era.astype(str)
fd["sig"] = fd.fe.map(sig_map).fillna(0.2).clip(lower=0.02)
fd["w"] = 1.0 / fd.sig ** 2

# restrict cells to those illuminated by the subset; slice kernel + rebuild Laplacian on the subset
cells_sub = sorted({c for c in fd.cell.unique()})
csub_idx = np.array([cellpos[c] for c in cells_sub]); Mc = len(csub_idx)
Gs = Gc_full[:, csub_idx]
_, nbr = cKDTree(cxy[csub_idx]).query(cxy[csub_idx], k=min(7, Mc))
L = np.zeros((Mc, Mc))
for i in range(Mc):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L

fe_ids = pd.Index(sorted(fd.fe.unique())); fd["fei"] = fe_ids.get_indexer(fd.fe)
sta_ids = pd.Index(sorted(fd.tag.unique())); fd["staj"] = sta_ids.get_indexer(fd.tag)
days = pd.Index(sorted(fd.date.unique())); fd["ti"] = days.get_indexer(fd.date)
P, NS, T = len(fe_ids), len(sta_ids), len(days)
ROW = fd.grow.values; TI = fd.ti.values; FEI = fd.fei.values; STJ = fd.staj.values
ERA = fd.era.values.astype(int)
DVV = fd.dvv.values.astype(np.float64); Wr = fd.w.values.astype(np.float64)
g1 = fd.groupby("fei")
fe_st = g1.staj.first().values; fe_era = g1.era.first().values.astype(int)
fe_row = g1.grow.first().values
fe_wsum = np.bincount(FEI, weights=Wr, minlength=P)
BY_DAY = [np.where(TI == t)[0] for t in range(T)]
lam_s = LAM_F * LAM_RATIO
N = len(fd)
print(f"reduced problem: rows {N:,} | frames {T} | cells {Mc} | stations {NS} | family-eras {P}")

# station-day site columns present in data
sd = pd.DataFrame({"ti": TI, "staj": STJ}).drop_duplicates().sort_values(["ti", "staj"])
sd["col"] = np.arange(len(sd))
sdmap = {(t, j): k for t, j, k in zip(sd.ti, sd.staj, sd.col)}
NSD = len(sd)
SDCOL = np.array([sdmap[(t, j)] for t, j in zip(TI, STJ)])


def gauge_all(M, S, c):
    """identical post-hoc gauge for both solutions: mean_t m=0 per cell; mean_t s=0 per (station,era)."""
    mu = M.mean(0)
    M = M - mu[None, :]
    c = c + Gs[fe_row] @ mu
    for j in range(NS):
        for e in np.unique(fe_era[fe_st == j]):
            tid = np.unique(TI[(STJ == j) & (ERA == e)])
            if not len(tid): continue
            vals = S[tid, j]; vals = vals[np.isfinite(vals)]
            if not len(vals): continue
            b = float(vals.mean())
            S[tid, j] -= b
            c[(fe_st == j) & (fe_era == e)] += b
    return M, S, c


# ---------------- (a) BCD ----------------
def day_solve(t, c):
    idx = BY_DAY[t]
    w = np.sqrt(Wr[idx]); d = DVV[idx] - c[FEI[idx]]
    stj = STJ[idx]; ust = np.unique(stj)
    cm = {j: k for k, j in enumerate(ust)}; ns = len(ust)
    Sm = np.zeros((len(idx), ns)); Sm[np.arange(len(idx)), [cm[j] for j in stj]] = 1.0
    A = np.hstack([Gs[ROW[idx]], Sm]); Aw = A * w[:, None]
    AtA = Aw.T @ Aw
    reg = np.zeros_like(AtA)
    reg[:Mc, :Mc] = LAM_F ** 2 * LtL; reg[np.arange(Mc), np.arange(Mc)] += 1e-4
    si = np.arange(Mc, Mc + ns); reg[si, si] += lam_s ** 2
    mfull = inv(AtA + reg) @ (Aw.T @ (w * d))
    return mfull[:Mc], ust, mfull[Mc:]


def c_update(M, S):
    pred = np.einsum("ij,ij->i", Gs[ROW], M[TI]) + S[TI, STJ]
    return np.bincount(FEI, weights=Wr * (DVV - pred), minlength=P) / np.maximum(fe_wsum, 1e-12)


def objective(M, S, c):
    pred = np.einsum("ij,ij->i", Gs[ROW], M[TI]) + S[TI, STJ] + c[FEI]
    reg = LAM_F ** 2 * float(np.sum((M @ L.T) ** 2)) + 1e-4 * float(np.sum(M ** 2)) + \
        lam_s ** 2 * float(np.nansum(S ** 2))
    return float(np.sum(Wr * (DVV - pred) ** 2)) + reg


c = c_update(np.zeros((T, Mc)), np.zeros((T, NS)))
objs = []
for it in range(SWEEPS):
    M = np.zeros((T, Mc)); S = np.zeros((T, NS))
    for t in range(T):
        m, ust, sv = day_solve(t, c)
        M[t] = m; S[t, ust] = sv
    cn = c_update(M, S)
    dC = float(np.max(np.abs(cn - c)))
    c = cn
    objs.append(objective(M, S, c))
    if it and objs[-2] - objs[-1] < 1e-10 * objs[-1] and dC < 1e-6:
        break
mono = all(objs[i + 1] <= objs[i] + 1e-9 * objs[i] for i in range(len(objs) - 1))
Mb, Sb, cb = gauge_all(M.copy(), S.copy(), c.copy())
print(f"BCD: {len(objs)} sweeps, monotone objective: {mono}, final obj {objs[-1]:.6e}")

# ---------------- (b) direct sparse joint solve ----------------
# unknowns: [m (T*Mc) | s (NSD station-days) | c (P)]
nM = T * Mc
rows_i = np.repeat(np.arange(N), Mc)
cols_m = (TI[:, None] * Mc + np.arange(Mc)[None, :]).ravel()
vals_m = Gs[ROW].ravel()
Am = sp.csr_matrix((vals_m, (rows_i, cols_m)), shape=(N, nM))
As = sp.csr_matrix((np.ones(N), (np.arange(N), SDCOL)), shape=(N, NSD))
Ac = sp.csr_matrix((np.ones(N), (np.arange(N), FEI)), shape=(N, P))
A = sp.hstack([Am, As, Ac]).tocsr()
Wd = sp.diags(Wr)
AtWA = (A.T @ Wd @ A).tocsr()
# regularization: Laplacian per day on m, ridge on s; c unregularized (+tiny numeric)
LtLs = sp.csr_matrix(LtL)
regm = sp.block_diag([LAM_F ** 2 * LtLs + 1e-4 * sp.eye(Mc)] * T)
reg = sp.block_diag([regm, lam_s ** 2 * sp.eye(NSD), 1e-9 * sp.eye(P)])
AtWd = A.T @ (Wr * DVV)
x = spsolve((AtWA + reg).tocsc(), AtWd)
Mdir = x[:nM].reshape(T, Mc)
Sdir = np.zeros((T, NS))
Sdir[sd.ti.values, sd.staj.values] = x[nM:nM + NSD]
cdir = x[nM + NSD:]
Md, Sd, cd = gauge_all(Mdir, Sdir, cdir)

# ---------------- compare ----------------
dm = np.abs(Mb - Md); ds = np.abs(Sb - Sd); dc = np.abs(cb - cd)
# reference scale: per-frame posterior sigma from one representative day
t0 = T // 2
idx = BY_DAY[t0]; w = np.sqrt(Wr[idx])
ust = np.unique(STJ[idx]); ns = len(ust)
cm = {j: k for k, j in enumerate(ust)}
Sm = np.zeros((len(idx), ns)); Sm[np.arange(len(idx)), [cm[j] for j in STJ[idx]]] = 1.0
A1 = np.hstack([Gs[ROW[idx]], Sm]); Aw = A1 * w[:, None]
AtA = Aw.T @ Aw
reg1 = np.zeros_like(AtA); reg1[:Mc, :Mc] = LAM_F ** 2 * LtL
reg1[np.arange(Mc), np.arange(Mc)] += 1e-4
reg1[np.arange(Mc, Mc + ns), np.arange(Mc, Mc + ns)] += lam_s ** 2
Minv = inv(AtA + reg1)
sig_m = float(np.median(np.sqrt(np.clip(np.diag(Minv @ AtA @ Minv)[:Mc], 0, None))))
print(f"\nmax|dM| {dm.max():.2e}%  median {np.median(dm):.2e}%   (median per-frame sigma_m {sig_m:.3f}%)")
print(f"max|dS| {ds.max():.2e}%  max|dc| {dc.max():.2e}%")
ok = dm.max() < 0.1 * sig_m
print(f"\nACCEPT (max|dM| < 0.1 sigma_m): {'PASS' if ok else 'FAIL'} | monotone: {'PASS' if mono else 'FAIL'}")
