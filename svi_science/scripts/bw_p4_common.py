#!/usr/bin/env python
"""bw_p4_common -- shared pieces of P4 (thresholds + quiet flag) of the SVI backward extension.

Read-only re-implementation of the forward false-alarm rule (lfe_extension/scripts/c_tiers_by_epoch.py) on the
svi_ref_mf feature files, plus hf_pass computed exactly as c_run_month.py does (the feature files carry hfr / hf_fs,
not hf_pass).  Rules: svi_science/results/backward/p4/P4.md (written before computing).
"""
from __future__ import annotations
import os, sys, glob
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
LFE = f"{ROOT}/lfe_extension"
sys.path.insert(0, f"{LFE}/scripts")
import c_common as C                                      # noqa: E402

BW = f"{ROOT}/svi_science/results/backward"
P4 = f"{BW}/p4"
TAB = f"{P4}/tables"
REFINE = f"{LFE}/results/refine"
FAMS12 = None
GRID = np.round(np.arange(8.0, 15.0 + 1e-9, 0.25), 4)
PNSN_START = "2010-01-10"                                 # catalogue starts 2010-01-06; +-3 d rule needs 3 d before
EPISODES = f"{ROOT}/svi_science/results/backward_scope/truth_episodes.csv"


def log(m):
    import time
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def pnsn_quiet(t0="2003-01-01", t1="2026-09-30"):
    """per UTC date: PNSN tremor count in the family box, quiet (Track C rule) and post-SSE (c_tiers rule).
    Days before PNSN_START: quiet = post = NaN (no catalogue)."""
    import svi_ref_common as R
    F = pd.read_csv(R.FAMS12); lat0, lat1, lon0, lon1 = C.family_box(F)
    tr = pd.read_csv(C.PNSN, usecols=["time", "lat", "lon"]); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"]); tr = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1)]
    cnt = tr.time.dt.strftime("%Y-%m-%d").value_counts()
    days = pd.date_range(t0, t1, freq="D")
    td = pd.to_datetime(pd.Index(cnt[cnt > 0].index)).sort_values()
    Q = pd.DataFrame(dict(date=[str(d.date()) for d in days]))
    Q["n_tremor_day"] = Q.date.map(cnt).fillna(0).astype(int)
    di = days.values.astype("datetime64[D]").astype(np.int64); ti = td.values.astype("datetime64[D]").astype(np.int64)
    j = np.searchsorted(ti, di)
    dl = np.abs(di - ti[np.clip(j - 1, 0, len(ti) - 1)]); dh = np.abs(ti[np.clip(j, 0, len(ti) - 1)] - di)
    Q["quiet"] = (Q.n_tremor_day == 0) & (np.minimum(dl, dh) > 3)
    # post-SSE: 1..10 d after the last tremor day of an episode (episode ends when >= 3 tremor-free days follow)
    tdl = list(td)
    ends = [tdl[i] for i in range(len(tdl)) if i == len(tdl) - 1 or (tdl[i + 1] - tdl[i]).days >= 3]
    post = set()
    for e in ends:
        for k in range(1, 11):
            post.add((e + pd.Timedelta(days=k)).strftime("%Y-%m-%d"))
    post -= {t.strftime("%Y-%m-%d") for t in tdl}
    Q["post"] = Q.date.isin(post)
    pre = Q.date < PNSN_START
    Q["quiet"] = Q.quiet.astype(object); Q["post"] = Q.post.astype(object)
    Q.loc[pre, "quiet"] = np.nan; Q.loc[pre, "post"] = np.nan
    return Q


def bostock_near(dates, pad=3):
    """True where the date is within +-pad days of a Bostock ETS window (inclusive window ends)."""
    E = pd.read_csv(EPISODES)
    d = pd.to_datetime(pd.Series(list(dates)))
    out = np.zeros(len(d), bool)
    for r in E.itertuples():
        a = pd.Timestamp(r.start) - pd.Timedelta(days=pad); b = pd.Timestamp(r.end) + pd.Timedelta(days=pad)
        out |= ((d >= a) & (d <= b)).to_numpy()
    return out


_HFC = None


def hf_pass(D):
    """c_run_month.py's hf_pass: ~(isfinite(hfr) & hfr > hf_cut), hf_cut from results/a4_hf_veto_cuts.csv (PGC)."""
    global _HFC
    if _HFC is None:
        _HFC = C.hf_cuts()
    fz = C.frozen(); cut0 = float(fz.get("hf_veto", {}).get("cut", np.nan))
    if not len(D):
        return np.zeros(0, bool)
    keys = pd.DataFrame(dict(date=D.date.astype(str).to_numpy(), fs=D.hf_fs.to_numpy(float)))
    u = keys.drop_duplicates()
    cmap = {}
    for r in u.itertuples():
        cmap[(r.date, r.fs if np.isfinite(r.fs) else -1.0)] = C.hf_cut_for(_HFC, "PGC", r.date, r.fs, cut0)[0] if _HFC else cut0
    cut = np.array([cmap[(d, f if np.isfinite(f) else -1.0)] for d, f in zip(keys.date, keys.fs)], float)
    h = D.hfr.to_numpy(float)
    return ~(np.isfinite(h) & (h > cut))


def load_run(featdir, dates=None, cols=("series", "fam", "date", "t", "net_mad", "hfr", "hf_fs")):
    """candidates and family-days of one svi_ref_mf run directory (optionally only `dates`)."""
    ds = sorted(os.path.basename(p)[:10] for p in glob.glob(f"{featdir}/*.json"))
    if dates is not None:
        dset = set(dates); ds = [d for d in ds if d in dset]
    A, FD = [], []
    for d in ds:
        p = f"{featdir}/{d}.parquet"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                x = pd.read_parquet(p)
                if len(x):
                    A.append(x[[c for c in cols if c in x.columns]])
            except Exception as e:
                log(f"[load] {p}: {e}")
        q = f"{featdir}/{d}_fd.parquet"
        if os.path.exists(q) and os.path.getsize(q) > 0:
            FD.append(pd.read_parquet(q, columns=["series", "fam", "date", "observable"]))
    A = pd.concat(A, ignore_index=True) if A else pd.DataFrame(columns=list(cols))
    FD = pd.concat(FD, ignore_index=True) if FD else pd.DataFrame(columns=["series", "fam", "date", "observable"])
    FD["observable"] = FD.observable.astype(bool)
    return A, FD, ds


def fa_curve(A, FD, dates, grid=GRID, series="data", use_hf=True):
    """candidates per observable family-day on `dates`, at each T of grid -> (curve, n_family_days)."""
    dset = set(dates)
    n = int(((FD.series == series) & FD.observable & FD.date.isin(dset)).sum())
    if n == 0:
        return np.full(len(grid), np.nan), 0
    d = A[(A.series == series) & A.date.isin(dset)]
    if use_hf:
        d = d[hf_pass(d)]
    m = np.sort(d.net_mad.to_numpy(float))
    return np.array([float(len(m) - np.searchsorted(m, T, side="left")) / n for T in grid]), n


def forward_rule(A, FD, quiet_dates, post_dates=(), min_quiet=200):
    """c_tiers_by_epoch.py's rule -> dict(T_core, T_ext, n_quiet, n_post, FA at T and at 8.75 / 8.25, curve)."""
    q, nq = fa_curve(A, FD, quiet_dates)
    p_, npo = fa_curve(A, FD, post_dates) if len(post_dates) else (np.full(len(GRID), np.nan), 0)
    tc = next((float(T) for T, v in zip(GRID, q) if np.isfinite(v) and v <= 0.1), np.nan)
    te = next((float(T) for T, v, w in zip(GRID, q, p_)
               if np.isfinite(v) and v <= 1.0 and (not np.isfinite(w) or w <= 1.5)), np.nan)
    gi = lambda T: int(np.argmin(np.abs(GRID - T))) if np.isfinite(T) else 0
    status = "derived" if nq >= min_quiet else f"insufficient quiet family-days ({nq} < {min_quiet})"
    return dict(T_core=tc, T_ext=te, n_quiet_family_days=nq, n_post_family_days=npo, status=status,
                fa_quiet_at_T_core=float(q[gi(tc)]) if nq else np.nan,
                fa_quiet_at_T_ext=float(q[gi(te)]) if nq else np.nan,
                fa_post_at_T_ext=float(p_[gi(te)]) if npo else np.nan,
                fa_quiet_at_8p00=float(q[gi(8.0)]) if nq else np.nan,
                fa_quiet_at_8p25=float(q[gi(8.25)]) if nq else np.nan,
                fa_quiet_at_8p75=float(q[gi(8.75)]) if nq else np.nan,
                fa_quiet_at_9p00=float(q[gi(9.0)]) if nq else np.nan,
                fa_post_at_8p25=float(p_[gi(8.25)]) if npo else np.nan), q, p_


def ctl_rate(A, FD, T, dates):
    """G1: scrambled-control candidates >= T per observable control family-day (no HF veto)."""
    dset = set(dates)
    n = int(((A.series == "c1") & (A.net_mad >= T) & A.date.isin(dset)).sum())
    m = int(((FD.series == "c1") & FD.observable & FD.date.isin(dset)).sum())
    return n / max(m, 1)


def matched_T(A, FD, target, dates):
    """G1 (bw_gate_g1_score.matched_T): smallest T on a 0.01 grid >= 8.00 with control rate <= target."""
    grid = np.round(np.arange(8.0, 14.0, 0.01), 2)
    dset = set(dates)
    c = np.sort(A[(A.series == "c1") & A.date.isin(dset)].net_mad.to_numpy(float))
    m = int(((FD.series == "c1") & FD.observable & FD.date.isin(dset)).sum())
    for T in grid:
        if (len(c) - np.searchsorted(c, T, side="left")) / max(m, 1) <= target + 1e-12:
            return float(T)
    return float(grid[-1])


def up_to_grid(T, step=0.25):
    return float(np.ceil(np.round(T / step, 6)) * step) if np.isfinite(T) else np.nan
