#!/usr/bin/env python
"""s10_tide_phase -- S10 follow-up part A (rules: results/s10_dvv/tide_subday/PLAN.md, A1-A4).
Per (family, station, epoch) series: tidal / fake-cycle bin stacks in designs D1 (within-day differential), D2
(day-matched pooled) and D3 (simple pooled vs epoch grand stack), 2-4 and 3-5 s, for the observed binning and
1,000 random-shift control draws (control i). Scrambled control c2 through D1/D3 (control iii, 200 draws). Injection
of eps = a cos(theta_UDSS) (family 10). Outputs per-series sums S1 = sum CC^2 eps, S0 = sum CC^2 per bin; pooling,
fits and rules are in s10_tide_phase_summary.py.
    python s10_tide_phase.py --workers 6 [--fams 10,87,17] [--ndraw 1000]
"""
from __future__ import annotations
import os, sys, argparse, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
import scipy.sparse as sp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M
import s10_tide_common as T

NB = 8
SEED = 20261011
T_M2_H = 12.4206012
PHASE_VARS = ["th_udss", "th_sn", "th_fake12", "th_fake13.5", "th_fake10.5"]
VALUE_VARS = ["udss", "udss_rate", "sn", "sn_rate"]
DAY_MODE_VARS = ["th_udss", "udss"]          # control (i) column: one common offset per (family, day)
C2_VARS = ["th_udss", "udss"]
INJ_A = (0.001, 0.0003)
OUTD = f"{T.TB}/phase_series"

G = {}   # globals shared with forked workers


def phase_code(th):
    return (np.floor(((np.asarray(th) + 22.5) % 360.0) / 45.0).astype(np.int8)) % NB


def build_codes(nd):
    """observed + control-(i) bin codes for every row and variable; row 0 = observed."""
    R = pd.read_parquet(f"{T.TB}/rows_tide.parquet")
    R = R.sort_values("row_id").reset_index(drop=True)
    H = np.load(f"{T.TB}/hourly_3fam.npz")
    fams = list(H["fams"]); t0 = int(H["t0_ns"])
    fi = R.fam.map({f: k for k, f in enumerate(fams)}).to_numpy()
    tns = pd.to_datetime(R.t, unit="s").astype("datetime64[ns]").to_numpy().astype(np.int64)
    rng = np.random.default_rng(SEED)
    n = len(R)
    # (family, day) groups for the day-common mode
    fd = pd.factorize(R.fam.astype(str) + "|" + R.date.astype(str))[0]
    nfd = fd.max() + 1
    off_rep = rng.uniform(0, 1, (nd, n)).astype(np.float32)          # fraction of a cycle
    off_day = rng.uniform(0, 1, (nd, nfd)).astype(np.float32)[:, fd]
    series = {}
    arrs = {"udss": H["udss"].astype(np.float64), "sn": H["sn"].astype(np.float64)}
    for nm in ("udss", "sn"):
        A = arrs[nm]; D = np.full_like(A, np.nan); D[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / 2.0; arrs[nm + "_rate"] = D

    def val(nm, tshift_ns):
        A = arrs[nm]
        h = (tns[None, :] + tshift_ns - t0) / 3.6e12
        h0 = np.floor(h).astype(np.int64); w = h - h0
        return A[fi[None, :], h0] * (1 - w) + A[fi[None, :], h0 + 1] * w

    edges = {}
    data = (R.series == "data").to_numpy()
    for v in VALUE_VARS:
        x = R[v].to_numpy()
        e = np.full((len(fams), NB + 1), np.nan)
        for k, f in enumerate(fams):
            m = data & (fi == k)
            e[k] = np.quantile(x[m], np.linspace(0, 1, NB + 1))
        edges[v] = e
    xmean = {}
    for v in PHASE_VARS:
        th = R[v].to_numpy()
        codes = {"rep": np.vstack([phase_code(th)[None, :], phase_code(th[None, :] + 360.0 * off_rep)])}
        if v in DAY_MODE_VARS:
            codes["day"] = np.vstack([phase_code(th)[None, :], phase_code(th[None, :] + 360.0 * off_day)])
        series[v] = codes
    for v in VALUE_VARS:
        e = edges[v]
        def cod(x):
            out = np.zeros(x.shape, np.int8)
            for k in range(len(fams)):
                m = fi == k
                out[..., m] = np.clip(np.searchsorted(e[k, 1:-1], x[..., m], side="right"), 0, NB - 1)
            return out
        x0 = R[v].to_numpy()
        codes = {"rep": np.vstack([cod(x0)[None, :], cod(val(v, (off_rep * T_M2_H * 3.6e12).astype(np.int64)))])}
        if v in DAY_MODE_VARS:
            codes["day"] = np.vstack([cod(x0)[None, :], cod(val(v, (off_day * T_M2_H * 3.6e12).astype(np.int64)))])
        series[v] = codes
        # bin mean value (kPa or kPa/h) per family, observed data rows
        xm = np.full((len(fams), NB), np.nan)
        c0 = codes["rep"][0]
        for k in range(len(fams)):
            for b in range(NB):
                m = data & (fi == k) & (c0 == b)
                xm[k, b] = x0[m].mean() / 1e3 if m.any() else np.nan
        xmean[v] = xm
    return R, series, edges, xmean, fams


# ------------------------------------------------------------------ per-series engine
class Eng:
    def __init__(self, Z, P, day, tS, G1):
        self.n = len(P); self.P = P; self.day = day; self.nd = day.max() + 1
        self.seg, self.Bd, self.Bdu, self.BG = {}, {}, {}, {}
        # day all-bin stacks (P-weighted, full trace) -> banks
        Wm = sp.csr_matrix((P, (day, np.arange(self.n))), shape=(self.nd, self.n))
        Rd = np.asarray(Wm @ Z.astype(np.float64)) / np.asarray(Wm.sum(1))
        for k, win in T.WINS.items():
            idx = T.seg_idx(tS, win)
            self.seg[k] = Z[:, idx].astype(np.float64)
            Bu = np.stack([T.bank_u(Rd[d], tS, win)[0] for d in range(self.nd)])       # (nd, ne, K)
            self.Bdu[k] = Bu
            self.Bd[k] = np.swapaxes(T.norm_bank(Bu), 1, 2).copy()                     # (nd, K, ne)
            self.BG[k] = T.norm_bank(T.bank_u(G1, tS, win)[0]) if G1 is not None else None

    def run(self, code, designs=("D1", "D2", "D3")):
        """code (n,) bin per window -> dict design|win -> (S1, S0, S3, N) arrays of shape (NB,)"""
        out = {}
        n, nd = self.n, self.nd
        cell = self.day.astype(np.int64) * NB + code
        Wc = sp.csr_matrix((self.P, (cell, np.arange(n))), shape=(nd * NB, n))
        cnt = np.bincount(cell, minlength=nd * NB).reshape(nd, NB)
        Wb = sp.csr_matrix((self.P, (code.astype(np.int64), np.arange(n))), shape=(NB, n))
        cntb = np.bincount(code, minlength=NB)
        Wdb = np.asarray(Wc.sum(1)).reshape(nd, NB)                                        # sum P per day-bin
        for k in T.WINS:
            if "D1" in designs:
                C = np.asarray(Wc @ self.seg[k]).reshape(nd, NB, -1)
                C[cnt < 3] = 0.0
                e, c = T.measure_bank(C, np.swapaxes(self.Bd[k], 1, 2))
                e, c = T.unrail(e, c)
                out[f"D1|{k}"] = _sums(e, c, axis=0, ncell=np.isfinite(e).sum(0))
            Sb = np.asarray(Wb @ self.seg[k]); Sb[cntb < 3] = 0.0
            if "D2" in designs:
                ws = Wdb.sum(0); ws[ws == 0] = 1
                Bref = (Wdb.T @ self.Bdu[k].reshape(nd, -1)).reshape(NB, *self.Bdu[k].shape[1:]) / ws[:, None, None]
                e, c = T.measure_bank(Sb[:, None, :], T.norm_bank(Bref))
                e, c = T.unrail(e[:, 0], c[:, 0])
                out[f"D2|{k}"] = _sums(e[None], c[None], axis=0, ncell=np.isfinite(e).astype(int))
            if "D3" in designs and self.BG[k] is not None:
                e, c = T.measure_bank(Sb, self.BG[k])
                e, c = T.unrail(e, c)
                out[f"D3|{k}"] = _sums(e[None], c[None], axis=0, ncell=np.isfinite(e).astype(int))
        return out


def _sums(e, c, axis, ncell):
    w = np.where(np.isfinite(e), c ** 2, 0.0)
    return np.stack([np.nansum(w * np.nan_to_num(e), axis), w.sum(axis), np.nansum(w * np.nan_to_num(c), axis), ncell])


def run_series(fn):
    z = np.load(fn, allow_pickle=False)
    fam, ns, ep = int(z["fam"]), str(z["netsta"]), str(z["epoch"])
    tag = f"{fam}_{ns}_{ep}"
    of = f"{OUTD}/{tag}.npz"
    if os.path.exists(of):
        return tag, "cached"
    pos = G["pos"]
    ri = np.array([pos[r] for r in z["row_id"]])
    days, day = np.unique(z["date"], return_inverse=True)
    Z = z["Z"]; P = z["P"].astype(np.float64); tS = float(z["tS"])
    eng = Eng(Z, P, day, tS, z["G1"])
    res = {}
    nd1 = G["nd"] + 1
    for v, modes in G["codes"].items():
        for mode, CODES in modes.items():
            cc = CODES[:, ri]
            acc = None
            for k in range(nd1):
                o = eng.run(cc[k])
                if acc is None:
                    acc = {key: np.zeros((nd1,) + val.shape, np.float32) for key, val in o.items()}
                for key, val in o.items():
                    acc[key][k] = val
            for key, val in acc.items():
                res[f"{v}|{mode}|{key}"] = val
    # control iii: scrambled-control (c2) windows, own day references for D1, data grand stack for D3
    if len(z["row_id2"]) >= 6:
        ri2 = np.array([pos[r] for r in z["row_id2"]])
        d2u, day2 = np.unique(z["date2"], return_inverse=True)
        eng2 = Eng(z["Z2"], np.ones(len(ri2)), day2, tS, z["G1"])
        for v in C2_VARS:
            CODES = G["codes"][v]["rep"][: G["nd_c2"] + 1, ri2]
            acc = None
            for k in range(len(CODES)):
                o = eng2.run(CODES[k], designs=("D1", "D3"))
                if acc is None:
                    acc = {key: np.zeros((len(CODES),) + val.shape, np.float32) for key, val in o.items()}
                for key, val in o.items():
                    acc[key][k] = val
            for key, val in acc.items():
                res[f"c2:{v}|rep|{key}"] = val
    # injection (family 10): eps_i = a cos(theta_UDSS,i) on every data window; observed binning only
    if fam == 10:
        th = G["th_udss"][ri]
        for a in INJ_A:
            Zi = np.stack([M.stretch_trace(Z[i].astype(np.float64), tS, a * np.cos(np.radians(th[i]))) for i in range(len(Z))])
            G1i = (P @ Zi) / P.sum()
            engi = Eng(Zi, P, day, tS, G1i)
            o = engi.run(G["codes"]["th_udss"]["rep"][0, ri])
            for key, val in o.items():
                res[f"inj{a:g}|{key}"] = val[None].astype(np.float32)
    # leave-bin-out day reference (D1 column, observed, th_udss and udss)
    for v in ("th_udss", "udss"):
        code = G["codes"][v]["rep"][0, ri]
        for k, win in T.WINS.items():
            S1 = np.zeros(NB); S0 = np.zeros(NB)
            for d in range(len(days)):
                md = day == d
                for b in range(NB):
                    mc = md & (code == b)
                    if mc.sum() < 3 or (md & ~mc).sum() < 3:
                        continue
                    ref = (P[md & ~mc] @ Z[md & ~mc]) / P[md & ~mc].sum()
                    cur = (P[mc] @ Z[mc]) / P[mc].sum()
                    e, c = M.Ref(ref, tS, win).measure(cur)
                    if np.abs(e[0]) < T.RAIL:
                        S1[b] += c[0] ** 2 * e[0]; S0[b] += c[0] ** 2
            res[f"{v}|lbo|D1|{k}"] = np.stack([S1, S0])[None].astype(np.float32)
    np.savez_compressed(of, fam=fam, netsta=ns, epoch=ep, n=len(ri), n_days=len(days), **{k.replace("|", "~"): v for k, v in res.items()})
    return tag, "done"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--ndraw", type=int, default=1000); ap.add_argument("--nd_c2", type=int, default=200)
    ap.add_argument("--fams", default="10,87,17"); ap.add_argument("--only", default="")
    ap.add_argument("--outd", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    global OUTD
    if a.outd:
        OUTD = a.outd
    T.ensure(); os.makedirs(OUTD, exist_ok=True)
    log = S.Log(f"{T.LOGD}/phase{a.tag}.log")
    t0 = time.time()
    R, codes, edges, xmean, fams = build_codes(a.ndraw)
    G.update(codes=codes, nd=a.ndraw, nd_c2=a.nd_c2, pos={r: i for i, r in enumerate(R.row_id)},
             th_udss=R.th_udss.to_numpy())
    np.savez(f"{T.TB}/phase_bins_meta{a.tag}.npz", fams=np.array(fams), **{f"edges_{v}": e for v, e in edges.items()},
             **{f"xmean_{v}": x for v, x in xmean.items()})
    log(f"codes built ({a.ndraw} draws) {time.time() - t0:.0f}s")
    files = sorted(f"{T.CACHE}/{f}" for f in os.listdir(T.CACHE))
    for fam in [int(x) for x in a.fams.split(",")]:
        fl = [f for f in files if os.path.basename(f).startswith(f"{fam}_") and (a.only in f)]
        import multiprocessing as mp
        with mp.get_context("fork").Pool(a.workers) as p:
            for tag, st in p.imap_unordered(run_series, fl):
                log(f"  {tag}: {st} {time.time() - t0:.0f}s")
        open(f"{T.OUT}/PHASE{a.tag}_fam{fam}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
        log(f"family {fam} done {time.time() - t0:.0f}s")
    open(f"{T.OUT}/PHASE{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
