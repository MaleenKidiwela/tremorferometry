#!/usr/bin/env python
"""SVI_REF R3 prep (CPU) -- the southern cleaning (lfe_south/scripts/sh_v2_common.py) ported to the SVI 20 Hz cache.

The SVI cache is 1-8 Hz Butterworth-4 zero-phase at the native rate, decimated to 20 Hz, gaps merged with
fill_value=0 BEFORE the filter (c_cache_month.py).  A gap / dead channel is therefore float residue
(1e-11 .. 1e-39) rather than exact zeros, so the sh_v2 'exactly zero or repeated' rule is kept AND extended by a
station-relative residue rule.  The band is NOT changed (no 2 Hz high-pass: that is Shelly's band, not a cleaning
step).  Per station-day, written to OUT/prep/<date>.npz (nothing else is touched):

  1. flat samples   exactly 0, exactly repeating the previous sample, or |x| <= 1e-6 x that channel's median |x|
                    (residue); runs >= 0.5 s are invalid, padded by 2 s (the sh_v2 FLAT_PAD_S)       -> <NS>|flat
  2. dead channels  plain daily median 2-s energy <= 1e-4 x the station's loudest channel (sh_v2 GUARD_STA_REL), or
                    the station's loudest < 1e-6 x the day's median station level (GUARD_NET_REL)     -> <NS>|dead
  3. burst muting   sh_v2 ratios/mute_from verbatim (K = 15, qmax 2, qlong 1.3, pad 3 s, 2-s windows, +-3 windows
                    moveout, 5 nearest other sites, sites < 3 km) on the 1-8 Hz 2-s energies          -> <NS>|mute
  4. line notching  sh_v2 line_scan/find_lines (hourly Welch, daily median prominence >= 4, narrow +-0.05 Hz,
                    2-9 Hz), applied by the matched filter as the sh_v2 0.06-Hz zero-phase notch       -> <NS>|lines
Status per station-day -> OUT/status/prep_<tag>.csv (guarded: > 3 % error/missing fails the stage).

  python scripts/svi_ref_prep.py --days-file OUT/days/prep_all.json --workers 10 --tag all
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings, traceback
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C
sys.path.insert(0, f"{R.SOUTH}/scripts")
import sh_v2_common as V                                 # functions only (ratios, mute_from, line_scan, find_lines, notch)

FS = R.FS
RESID_REL = 1e-6
PREP = f"{R.OUT}/prep"


def coords():
    T = C.station_table()
    ns = R.netsta_of()
    return {ns[s]: (float(T.loc[s, "lat"]), float(T.loc[s, "lon"])) for s in T.index if s in ns}


def cache_dirs_for(date):
    y = int(date[:4])
    if y in R.TRUTH:
        return R.TRUTH[y]["caches"]
    return [R.CACHE_R1, C.CACHE]


def find_cache(ns, date, dirs):
    d = pd.Timestamp(date)
    for c in dirs:
        p = f"{c}/{ns}/{d.year}/{d.dayofyear:03d}.npz"
        if os.path.exists(p):
            return p
    return None


def flat_runs(x):
    """(C, N) -> list per channel of padded invalid [a, b) sample intervals."""
    Cn, N = x.shape
    out = []
    runmin = int(0.5 * FS); w = int(V.FLAT_PAD_S * FS)
    for c in range(Cn):
        a_ = np.abs(x[c]); med = float(np.median(a_))
        fl = (x[c] == 0) | (a_ <= RESID_REL * med)
        fl[1:] |= x[c][1:] == x[c][:-1]
        iv = []
        if fl.any():
            d = np.diff(np.concatenate([[0], fl.astype(np.int8), [0]]))
            aa = np.flatnonzero(d == 1); bb = np.flatnonzero(d == -1)
            for i, j in zip(aa, bb):
                if j - i >= runmin:
                    iv.append((max(i - w, 0), min(j + w, N)))
        # merge overlaps
        m = []
        for i, j in sorted(iv):
            if m and i <= m[-1][1]: m[-1][1] = max(m[-1][1], j)
            else: m.append([i, j])
        out.append(np.array(m, np.int64).reshape(-1, 2))
    return out


def mask_from(N, iv):
    mc = np.ones(N, bool)
    for i, j in iv:
        mc[i:j] = False
    return mc


def prep_day(date, CO, write=True):
    d = pd.Timestamp(date)
    out = f"{PREP}/{date}.npz"
    status = []
    dirs = cache_dirs_for(date)
    RAW = {}
    for ns in sorted(CO):
        p = find_cache(ns, date, dirs)
        if p is None:
            continue
        try:
            z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"])
            if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * FS:
                status.append((date, ns, "short", "cache < 12 h")); continue
            RAW[ns] = (t0, x)
            status.append((date, ns, "ok", ""))
        except Exception as e:
            status.append((date, ns, "error", f"cache unreadable: {str(e)[:60]}"))
    S = sorted(RAW)
    save, rows = {}, []
    if not S:
        if write: np.savez_compressed(out + ".tmp.npz", stas=np.array([], dtype=str)); os.replace(out + ".tmp.npz", out)
        return status, rows
    CL = {}
    for ns in S:
        t0, x = RAW[ns]
        fr = flat_runs(x)
        mc = np.stack([mask_from(x.shape[1], iv) for iv in fr])
        X = np.where(mc, x, 0.0)
        CL[ns] = (t0, X, mc, fr)
    E = [V.energy_2s(CL[s][1], CL[s][2], CL[s][0], d) for s in S]
    with np.errstate(invalid="ignore"):
        PM = [np.nanmedian(e, axis=1) for e in E]
    DEADC = {}
    for s, pm in zip(S, PM):
        pm = np.where(np.isfinite(pm), pm, 0.0); top = pm.max()
        DEADC[s] = (pm <= V.GUARD_STA_REL * top) if top > 0 else np.ones(pm.size, bool)
    stop = np.array([np.nanmax(np.where(np.isfinite(pm), pm, 0.0)) for pm in PM])
    nfloor = V.GUARD_NET_REL * float(np.median(stop[stop > 0])) if (stop > 0).any() else np.inf
    for s, t_ in zip(S, stop):
        if t_ < nfloor:
            DEADC[s][:] = True
    E = [np.where(DEADC[s][:, None], np.nan, e) for s, e in zip(S, E)]
    if len(S) >= 2:
        Rr, Q, QL = V.ratios(E, S, CO)
        MU = V.mute_from(Rr, Q, QL)
    else:
        MU = np.zeros((len(S), V.NK), bool); Rr = np.full((len(S), V.NK), np.nan)
    for i, ns in enumerate(S):
        t0, X, mc, fr = CL[ns]
        ms = V.mute_samples(MU[i], t0, d, X.shape[1])
        mc2 = mc & ~ms[None, :]
        X2 = np.where(mc2, X, 0.0)
        r = V.line_scan(X2, mc2)
        L = np.array(V.find_lines(*r), float).reshape(-1, 3) if r is not None else np.zeros((0, 3))
        save[f"{ns}|mute"] = np.packbits(MU[i]); save[f"{ns}|lines"] = L; save[f"{ns}|dead"] = DEADC[ns]
        for c, iv in enumerate(fr):
            save[f"{ns}|flat{c}"] = iv
        live = np.isfinite(Rr[i])
        rows.append(dict(date=date, netsta=ns, nch=X.shape[0], N=X.shape[1], t0=t0,
                         flat_frac=float((~mc).mean()), dead_ch=int(DEADC[ns].sum()),
                         muted_frac=float((MU[i] & live).sum() / max(live.sum(), 1)), n_mute_win=int(MU[i].sum()),
                         n_lines=int(len(L)), lines=";".join(f"{int(c)}:{f0:.2f}" for c, f0, _ in L)))
    if write:
        np.savez_compressed(out + ".tmp.npz", stas=np.array(S), nk=V.NK, **save); os.replace(out + ".tmp.npz", out)
    return status, rows


class DayPrep:
    """matched-filter side: per station the cleaning sample mask and the line list."""
    def __init__(self, date, prep_dir=None):
        p = f"{prep_dir or PREP}/{date}.npz"
        self.date = date
        self.z = dict(np.load(p, allow_pickle=False)) if os.path.exists(p) else None
        self.stas = set(str(s) for s in self.z["stas"]) if self.z is not None else set()

    def ok(self):
        return self.z is not None

    def mask(self, ns, t0, N):
        """(N,) bool INVALID samples (flat runs on any channel | muted | a dead channel -> all), or None."""
        if ns not in self.stas:
            return None
        inv = np.zeros(N, bool)
        dead = self.z[f"{ns}|dead"]
        if dead.any():
            inv[:] = True; return inv
        c = 0
        while f"{ns}|flat{c}" in self.z:
            for i, j in self.z[f"{ns}|flat{c}"]:
                inv[max(int(i), 0):min(int(j), N)] = True
            c += 1
        mu = np.unpackbits(self.z[f"{ns}|mute"])[:V.NK].astype(bool)
        inv |= V.mute_samples(mu, t0, pd.Timestamp(self.date), N)
        return inv

    def lines(self, ns):
        return self.z[f"{ns}|lines"] if ns in self.stas else np.zeros((0, 3))

    def apply(self, ns, x, t0):
        """cleaned copy of x (C, N): invalid samples zeroed, lines notched over valid runs >= 60 s; + invalid mask."""
        inv = self.mask(ns, t0, x.shape[1])
        if inv is None:
            return x, None
        mc = np.broadcast_to(~inv, x.shape).copy()
        X = np.where(mc, x.astype(np.float64), 0.0)
        L = self.lines(ns)
        if len(L):
            X = V.notch(X, mc, L)
        return X.astype(np.float32), inv


_CO = {}
def _job(date):
    try:
        st, rows = prep_day(date, _CO)
        return date, st, rows, ""
    except Exception as e:
        return date, [(date, ns, "error", f"prep: {str(e)[:80]}") for ns in _CO], [], traceback.format_exc()[-800:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--days-file", required=True); ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--tag", default="all")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    R.ensure(PREP, f"{R.OUT}/status")
    days = json.load(open(a.days_file))
    todo = [d for d in days if not os.path.exists(f"{PREP}/{d}.npz")]
    _CO.update(coords())
    R.log(f"[prep {a.tag}] {len(todo)} of {len(days)} days to do, {a.workers} workers, {len(_CO)} stations")
    sp = f"{R.OUT}/status/prep_{a.tag}.csv"
    status = [tuple(x) for x in pd.read_csv(sp, dtype=str).fillna("").itertuples(index=False)] if os.path.exists(sp) else []
    rp = f"{R.OUT}/prep_summary_{a.tag}.parquet"
    rows = pd.read_parquet(rp).to_dict("records") if os.path.exists(rp) else []
    import multiprocessing as mp
    t0 = time.time()
    with mp.Pool(a.workers, maxtasksperchild=50) as pool:
        for n, (d, st, rw, tb) in enumerate(pool.imap_unordered(_job, todo, chunksize=1), 1):
            status += st; rows += rw
            if tb: R.log(f"  ERROR {d}: {tb}")
            if n % 50 == 0 or n == len(todo):
                R.log(f"  {n}/{len(todo)} {time.time() - t0:.0f}s")
                R.write_status("prep", a.tag, status); pd.DataFrame(rows).to_parquet(rp, index=False)
    R.write_status("prep", a.tag, status)
    if rows: pd.DataFrame(rows).to_parquet(rp, index=False)
    R.log(f"[prep {a.tag}] done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
