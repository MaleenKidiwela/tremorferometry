#!/usr/bin/env python
"""DAILY fixed-effects 4-D joint inversion (Merlin-audited design + rulings, 2026-08-22 —
notes/2026-08-22_daily_joint_inversion_design.md).

Model (data are eps = -(dv/v); G.npz stores -K, so eps = G·m AS STORED — sign-corrected 2026-08-22):

    d_{f,t} = K_f^T m(·,t) + s_{j(f)}(t) + c_{f,e(t)} + noise

  m(x,t)   interface dbeta/beta movie (cells x days)          — the field
  s_j(t)   per-station-DAY site term (absorbs station-common: shallow, weather, seasonal)
  c_{f,e}  per-FAMILY x instrument-era offset (measurement-reference constant; Merlin B1 ruling)

Gauge (the reference window R = 2019-2024 IS the gauge): after every sweep, exactly reparameterize so
  mean_{t in R} m(x,t) = 0 per cell,  mean_{t in era ∩ R} s_j(t) = 0 per (station, era)   [B2]
absorbed into c — data fit invariant. Eras/cells with no R coverage gauge on their own span (flagged).

Solver: block coordinate descent (exact alternating minimization of a strictly convex quadratic — I2):
  (i) given c: per-day independent WLS solves for [m_t, s_t] (Laplacian on m, ridge on s), parallel;
  (ii) given m,s: closed-form c = per-(family,era) weighted mean residual (I3: gate-passing days only).
Stop: relative objective decrease < 1e-8 AND max|dc|/sigma_c < 0.1. Validate first with
validate_daily_bcd.py (reduced problem, BCD vs direct joint solve) — Merlin step 3.

Lambda: oracle on representative single days spanning roster growth (2013/2018/2023), noise draws =
iid per-row + STATION-COMMON component (I5b), 100-km checkerboard recovery; ONE lambda for the whole
record (B4 — per-epoch lambda makes the response jump; amplitude evolution needs the separate
transmission normalization T(y,patch) before ANY slow-field claim).

Outputs (RESDIR): inversion_4d_daily.npz — MODEL (days x cells, TRUE dbeta/beta), SIGM (final-sweep
posterior), SITES (days x stations, TRUE near-station dv/v = -s, matching the annual product, I7a),
day list, gauge flags; daily_offsets.csv — c per family-era in TRUE-dv/v sign with N_eff-deflated
sigma (the step-scan readout). Movie temporal resolution is 30 d HARD (I6): value at t = trailing
[t-30d, t] stack (~15-d centroid lag); sub-30-d structure = artifact.
"""
import argparse, os
# Pin BLAS to 1 thread/process BEFORE numpy loads: pod is capped at 32 CPUs but advertises 176;
# unpinned BLAS thrashed the pod on the 07:31 launch (parallelism = the 28 process workers).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
from concurrent.futures import ProcessPoolExecutor

OUT = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
MIN_ROWS, MIN_SITES = 50, int(os.environ.get("MIN_SITES", "10"))   # env override for subset runs (2026-09-02)
R_WIN = (pd.Timestamp("2019-01-01"), pd.Timestamp("2024-12-31"))
LAM_DAYS = ["2013-06-15", "2018-06-15", "2023-06-15"]
MAXSWEEP, TOL_OBJ, TOL_C = 30, 1e-8, 0.1

# ---------------- kernel ----------------
gz = np.load(f"{OUT}/G.npz", allow_pickle=True)
G = gz["G"].astype(np.float64)                 # AS STORED (=-K): eps = G.m — NO flip (2026-08-22)
fcap = gz["captured"].astype(np.float64)
cellids = gz["cell"]; Mf = G.shape[1]; cxy = gz["cxy"]
keep = fcap > 0
Gc = (G * fcap[:, None])[keep]
ptag, pcell = gz["pair_tag"][keep], gz["pair_cell"][keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
cellpos = {c: i for i, c in enumerate(cellids)}
pcell_idx = np.array([cellpos[c] for c in pcell])
_, nbr = cKDTree(cxy).query(cxy, k=7)
L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T @ L

# ---------------- data ----------------
fd = pd.read_parquet(f"{OUT}/family_days.parquet")
fams = pd.read_csv(f"{OUT}/families_daily.csv")
fd["grow"] = [row_of.get((t, c), -1) for t, c in zip(fd.tag, fd.cell)]
n_unmapped = int((fd.grow < 0).sum())
fd = fd[fd.grow >= 0].copy()
fd["fe"] = fd.tag + "|" + fd.fam + "|e" + fd.era.astype(str)
sig_map = dict(zip(fams.tag + "|" + fams.fam + "|e" + fams.era.astype(str), fams.sigma))
fd["sig"] = fd.fe.map(sig_map).fillna(0.2).clip(lower=0.02)
fd["w"] = 1.0 / fd.sig ** 2
fd["date"] = pd.to_datetime(fd.date)

gd = fd.groupby("date").agg(nrows=("dvv", "size"), nsta=("tag", "nunique"))
good = gd[(gd.nrows >= MIN_ROWS) & (gd.nsta >= MIN_SITES)].index
n_gated = gd.index.size - good.size
fd = fd[fd.date.isin(good)].copy()

fe_ids = pd.Index(sorted(fd.fe.unique())); fd["fei"] = fe_ids.get_indexer(fd.fe)
sta_ids = pd.Index(sorted(fd.tag.unique())); fd["staj"] = sta_ids.get_indexer(fd.tag)
days = pd.Index(sorted(fd.date.unique())); fd["ti"] = days.get_indexer(fd.date)
P, NS, T = len(fe_ids), len(sta_ids), len(days)
inR = np.asarray((days >= R_WIN[0]) & (days <= R_WIN[1]))

ROW = fd.grow.values; TI = fd.ti.values; FEI = fd.fei.values; STJ = fd.staj.values
ERA = fd.era.values.astype(np.int16)
DVV = fd.dvv.values.astype(np.float64); Wr = fd.w.values.astype(np.float64)
g1 = fd.groupby("fei")
fe_st = g1.staj.first().values; fe_era = g1.era.first().values.astype(np.int16)
fe_row = g1.grow.first().values
fe_wsum = np.bincount(FEI, weights=Wr, minlength=P)
fe_n = np.bincount(FEI, minlength=P)
fe_span = (g1.date.max() - g1.date.min()).dt.days.values
fe_neff = np.clip(fe_span / 75.0, 1.0, None)                   # decorrelation ~75 d (I5a)
sigma_c = (1.0 / np.sqrt(np.maximum(fe_wsum, 1e-12))) * np.sqrt(np.maximum(fe_n / fe_neff, 1.0))
scom = fd.groupby(["staj", "ti"]).dvv.mean().groupby("staj").std()
SIG_COMMON = scom.reindex(range(NS)).fillna(scom.median()).values
BY_DAY = [np.where(TI == t)[0] for t in range(T)]
N_ROWS = len(fd)
del fd, g1, gd, sig_map                        # free the 30M-row frame before workers fork
import gc; gc.collect()


def predict(M, S):
    """Row predictions Gc[ROW]·M[TI] + S[TI,STJ], computed per-day.
    NEVER materialize Gc[ROW] / M[TI] whole: that is a 30.4M x Mf float64 = 212 GB
    allocation — it OOM-killed the pod on 2026-08-22."""
    pred = np.empty(N_ROWS)
    for t in range(T):
        idx = BY_DAY[t]
        pred[idx] = Gc[ROW[idx]] @ M[t] + S[t, STJ[idx]]
    return pred


def day_solve(t, c, lam_f, lam_s, want_cov=False, d_override=None):
    idx = BY_DAY[t]
    rows = ROW[idx]; w = np.sqrt(Wr[idx])
    d = (DVV[idx] if d_override is None else d_override) - c[FEI[idx]]
    stj = STJ[idx]; ust = np.unique(stj)
    colmap = {j: k for k, j in enumerate(ust)}; ns = len(ust)
    S = np.zeros((len(idx), ns))
    S[np.arange(len(idx)), [colmap[j] for j in stj]] = 1.0
    A = np.hstack([Gc[rows], S]); Aw = A * w[:, None]
    AtA = Aw.T @ Aw
    reg = np.zeros_like(AtA)
    reg[:Mf, :Mf] = lam_f ** 2 * LtL; reg[np.arange(Mf), np.arange(Mf)] += 1e-4
    si = np.arange(Mf, Mf + ns); reg[si, si] += lam_s ** 2
    Minv = inv(AtA + reg)
    mfull = Minv @ (Aw.T @ (w * d))
    if want_cov == "full":                            # full m-block posterior (B4 noise propagation)
        C = Minv @ AtA @ Minv
        return mfull[:Mf], ust, mfull[Mf:], C[:Mf, :Mf]
    sigm = np.sqrt(np.clip(np.diag(Minv @ AtA @ Minv)[:Mf], 0, None)) if want_cov else None
    return mfull[:Mf], ust, mfull[Mf:], sigm


def _worker(argt):
    t, c, lam_f, lam_s, want_cov = argt
    return day_solve(t, c, lam_f, lam_s, want_cov)


def sweep(c, lam_f, lam_s, want_cov=False, workers=1):
    M = np.zeros((T, Mf)); S = np.full((T, NS), np.nan)
    SG = np.full((T, Mf), np.nan, np.float32) if want_cov else None
    tasks = [(t, c, lam_f, lam_s, want_cov) for t in range(T)]
    if workers > 1:
        with ProcessPoolExecutor(workers) as ex:
            res = list(ex.map(_worker, tasks, chunksize=8))
    else:
        res = [_worker(a) for a in tasks]
    for t, (m, ust, sv, sg) in enumerate(res):
        M[t] = m; S[t, ust] = sv
        if want_cov: SG[t] = sg
    return M, S, SG


def c_update(M, S):
    r = DVV - predict(M, S)
    return np.bincount(FEI, weights=Wr * r, minlength=P) / np.maximum(fe_wsum, 1e-12)


def gauge(M, S, c):
    """Exact reparameterization (I1): applied EVERY sweep. Data fit invariant."""
    Rm = inR if inR.sum() >= 30 else np.ones(T, bool)
    mu = np.nanmean(np.where(np.isfinite(M[Rm]), M[Rm], 0.0), axis=0)   # cell level over R
    M = M - mu[None, :]
    c = c + Gc[fe_row] @ mu
    # site gauge per (station, era): era membership of a station's days from its rows
    for j in range(NS):
        for e in np.unique(fe_era[fe_st == j]):
            tid = np.unique(TI[(STJ == j) & (ERA == e)])
            if not len(tid): continue
            sel = np.zeros(T, bool); sel[tid] = True
            selR = sel & Rm
            vals = S[selR, j] if selR.any() else S[sel, j]
            vals = vals[np.isfinite(vals)]
            if not len(vals):
                vals = S[sel, j][np.isfinite(S[sel, j])]
            if not len(vals): continue
            base = float(vals.mean())
            S[sel, j] -= base
            c[(fe_st == j) & (fe_era == e)] += base
    return M, S, c


def objective(M, S, c):
    pred = predict(M, S) + c[FEI]
    return float(np.sum(Wr * (DVV - pred) ** 2))


def checker(scale):
    return np.sign(np.sin(np.pi * cxy[:, 0] / scale) * np.sin(np.pi * cxy[:, 1] / scale))


def oracle_lambda():
    ck = checker(100.0)
    tgt = [int(days.get_indexer([pd.Timestamp(x)], method="nearest")[0]) for x in LAM_DAYS]
    best, bestc = 0.5, -1.0
    for lam in np.logspace(-2, 1, 10):
        cs = []
        for t in tgt:
            idx = BY_DAY[t]
            dcl = Gc[ROW[idx]] @ ck
            ill = np.zeros(Mf, bool); ill[np.unique(pcell_idx[ROW[idx]])] = True
            for r in range(3):
                rng = np.random.RandomState(1000 * t + r)
                noise = rng.randn(len(idx)) / np.sqrt(Wr[idx]) + \
                    (SIG_COMMON * rng.randn(NS))[STJ[idx]]       # iid + station-common (I5b)
                m, _, _, _ = day_solve(t, np.zeros(P), lam, lam * 0.125, d_override=dcl + noise)
                cs.append(np.corrcoef(ck[ill], m[ill])[0, 1])
        if np.nanmean(cs) > bestc:
            bestc, best = float(np.nanmean(cs)), float(lam)
    print(f"lambda oracle (days {LAM_DAYS}, iid+station-common noise): lam_f={best:.3f} "
          f"(100 km recovery {bestc:.2f})")
    return best


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=28)
    ap.add_argument("--lam", type=float, default=0.0, help="skip oracle, use this lam_f")
    ap.add_argument("--lam-ratio", type=float, default=0.125)
    ap.add_argument("--max-sweeps", type=int, default=MAXSWEEP)
    args = ap.parse_args()
    print(f"family-days {N_ROWS:,} | family-eras {P} | stations {NS} | frames {T} "
          f"(day-gated out {n_gated}) | unmapped-to-G rows dropped {n_unmapped:,}")

    lam_f = args.lam if args.lam > 0 else oracle_lambda()
    lam_s = lam_f * args.lam_ratio

    c = c_update(np.zeros((T, Mf)), np.zeros((T, NS)))           # init = per-fe weighted mean
    prev = np.inf
    prev_cn = None
    for it in range(args.max_sweeps):
        M, S, _ = sweep(c, lam_f, lam_s, workers=args.workers)
        cn = c_update(M, np.nan_to_num(S))
        # convergence on successive UN-gauged optima: comparing cn to the post-gauge c measures the
        # constant gauge shuttle b (ridge-optimal s has nonzero era-mean; gauge moves it into c every
        # sweep), which pinned dC at 2.320 for 16 sweeps on the 2026-08-22 run — not true dc motion.
        dC = float(np.max(np.abs(cn - prev_cn) / sigma_c)) if prev_cn is not None else np.inf
        prev_cn = cn.copy()
        c = cn
        M, S, c = gauge(M, S, c)
        obj = objective(M, np.nan_to_num(S), c)
        rel = (prev - obj) / max(obj, 1e-12)
        print(f"sweep {it}: obj {obj:.6e} rel {rel:.2e} max|dc|/sig_c {dC:.3f}", flush=True)
        if rel < TOL_OBJ and dC < TOL_C:
            break
        prev = obj
    M, S, SG = sweep(c, lam_f, lam_s, want_cov=True, workers=args.workers)
    M, S, c = gauge(M, S, c)

    for t in range(T):                                            # illumination mask per frame
        ill = np.zeros(Mf, bool); ill[np.unique(pcell_idx[ROW[BY_DAY[t]]])] = True
        M[t][~ill] = np.nan; SG[t][~ill] = np.nan

    ctab = pd.DataFrame(dict(fe=fe_ids, tag=sta_ids[fe_st], era=fe_era,
                             c_true_dvv=np.round(-c, 4), sigma_c=np.round(sigma_c, 4),
                             n=fe_n, neff=np.round(fe_neff, 1)))
    ctab.to_csv(f"{OUT}/daily_offsets.csv", index=False)
    np.savez_compressed(f"{OUT}/inversion_4d_daily.npz",
                        days=days.astype(str).values, MODEL=M.astype(np.float32), SIGM=SG,
                        SITES=(-S).astype(np.float32),           # TRUE near-station dv/v (I7a)
                        cell=cellids, clat=gz["cell_lat"], clon=gz["cell_lon"],
                        depth=gz["depth_km"], lam_f=lam_f, lam_s=lam_s, inR=inR,
                        sta_tags=sta_ids.values.astype(object))
    print(f"wrote {OUT}/inversion_4d_daily.npz + daily_offsets.csv | frames {T} | lam_f {lam_f:.3f}")
