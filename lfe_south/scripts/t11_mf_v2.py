#!/usr/bin/env python
"""SH_V2 matched filter (2026-09-28) -- a new file built from t11_mf_aa.py (never edited).

Shelly et al. 2025 recipe, literally where the main text states it:
  band 2-9 Hz (Fig. 2) ; 6 s templates (Fig. 2a) ; network CORRELATION SUM (Fig. 4a) with each
  station's moveout (t11_mf_aa F2 fix: shift by tt_s + best_tau - pre) ; threshold k x DAILY MAD
  (Fig. 4a, k = 9; a sweep 6..12 is written as columns) ; keep the highest correlation sum in each
  FIXED 4 s window (sec. 3).
The statistic is the SUM over live channels (his words, "correlation sum"); the per-sample live
channel count and the mean cc (sum / live) are written beside it.

Cleaning (sh_v2_common.py, results/SH_V2_LOG.md): burst muting + line notching, exact block-local
window energies, relative gap guard, dead-channel and dead-station floors.  A channel-window enters
the sum only if every sample is valid and e > 1e-3 x the channel's reference energy.  Any |cc| > 1 +
1e-6 left after that is counted (col n_ccgt1 in the diag) and zeroed -- it should never happen.

Per day: each channel's spectrum and exact S1/S2 are computed ONCE and reused by all families.

  python scripts/t11_mf_v2.py --templates results/sh_v2/tpl_it0 --days-file results/aa_shelly/stage2_days_test.json \
      --truth shelly --t0 2018-01-01 --t1 2024-12-01 --out results/sh_v2/mf_it0_test.parquet
"""
import os, sys, time, argparse, warnings, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v2_common import *
try:
    import cupy as cp; XP = cp; GPU = True
except Exception:
    XP = np; GPU = False
def to_np(a): return cp.asnumpy(a) if GPU else a
SWEEP = [6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]


def decl_bin(nm, thr, gap_n):
    n = (nm.size // gap_n) * gap_n
    if n == 0: return np.array([], int)
    B = nm[:n].reshape(-1, gap_n); am = B.argmax(1); mx = B.max(1)
    ok = np.where(mx >= thr)[0]
    return ok * gap_n + am[ok]


class Chan:
    """one channel-day on the GPU: spectrum, exact window sums, window validity."""
    def __init__(self, x, mc, ref, L, nfft):
        N = x.size; self.N = N; self.L = L; self.nw = N - L + 1
        xg = XP.asarray(x)
        self.F = XP.fft.rfft(xg, nfft)
        starts = XP.arange(self.nw)
        self.S1 = box_sum(xg[None, :], L, starts, XP)[0]
        self.S2 = box_sum((xg * xg)[None, :], L, starts, XP)[0]
        inv = XP.asarray((~mc).astype(np.int32))
        cf = XP.concatenate([XP.zeros(1, dtype=XP.int64), XP.cumsum(inv)])
        ninv = cf[L:N + 1] - cf[:N - L + 1]
        var = L * self.S2 - self.S1 * self.S1
        self.ok = (ninv == 0) & (self.S2 > GUARD_REL * ref) & (var > 0) & (ref > 0)
        self.sd = XP.sqrt(XP.where(self.ok, var, 1.0))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--source", default="Shelly2025")
    ap.add_argument("--t0", required=True); ap.add_argument("--t1", required=True)
    ap.add_argument("--tmin", type=float, default=6.0); ap.add_argument("--gap", type=float, default=4.0)
    ap.add_argument("--out", required=True); ap.add_argument("--days-file", default="")
    ap.add_argument("--templates", required=True)
    ap.add_argument("--truth", default="", choices=["", "shelly"])
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    M = pd.read_csv(a.templates + ".csv"); M = M[M.source == a.source].copy(); M["fam"] = M.fam.astype(str)
    Z = np.load(a.templates + ".npz"); TPL = {k: Z[k].astype(np.float64) for k in M.key if k in Z}
    fams = sorted(M.fam.unique(), key=lambda s: int(s) if s.isdigit() else s)
    stas = sorted(M.netsta.unique())
    PRE = float(M.pre_s.iloc[0]); LEN = float(M.len_s.iloc[0]); L = int(round(LEN * FS))
    log(f"[mf-v2] {len(fams)} fam, {len(stas)} sta, {len(TPL)} tpl, GPU={GPU}, L={L}, pre={PRE}")
    TRU = None
    if a.truth == "shelly":
        S = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
        S.columns = [c.strip('%') for c in S.columns]
        S["t"] = (pd.to_datetime(dict(year=S.year, month=S.month, day=S.day, hour=S.hour, minute=S.minute,
                  second=S.second.astype(int))).astype("datetime64[ns]") + pd.to_timedelta(S.second % 1, unit="s"))
        TRU = S.rename(columns={"idnum": "fam", "meancc": "cc"})[["fam", "t", "cc"]]
        TRU["fam"] = TRU.fam.astype(str); TRU["day"] = TRU.t.dt.floor("D")
    DAYS = None
    if a.days_file:
        dd = json.load(open(a.days_file)); DAYS = set()
        for k, v in dd.items():
            if k.split("|")[0] == a.source: DAYS |= {pd.Timestamp(x) for x in v}
    gap_n = int(round(a.gap * FS))
    rows = []; diag = []; tru = []; t00 = time.time(); nday = 0
    for d in pd.date_range(a.t0, a.t1, freq="D", inclusive="left"):
        if DAYS is not None and d not in DAYS: continue
        P = DayPrep(d)
        if P.z is None: continue
        DAT = {}
        for ns in stas:
            v = P.load(ns)
            if v is not None: DAT[ns] = v
        if not DAT: continue
        nday += 1
        tref = pd.Timestamp(d).timestamp()
        NS = int(round(86400 * FS))
        N = max(v[1].shape[1] for v in DAT.values())
        nfft = int(2 ** np.ceil(np.log2(N + L)))
        CH = {}
        for ns, (t0, X, mc, ref) in DAT.items():
            CH[ns] = [Chan(X[c], mc[c], float(ref[c]), L, nfft) for c in range(X.shape[0])]
        for f in fams:
            sub = M[(M.fam == f) & M.netsta.isin(DAT)]
            SUM = XP.zeros(NS); LIVE = XP.zeros(NS); nch = 0; ngt1 = 0; nst = 0
            for r in sub.itertuples():
                T = TPL.get(r.key)
                if T is None: continue
                t0 = DAT[r.netsta][0]
                chs = CH[r.netsta]
                if T.shape[0] != len(chs): continue
                sh = int(round((float(r.tt_s) + float(r.best_tau) - PRE - (t0 - tref)) * FS))
                used = 0
                for c, ch in enumerate(chs):
                    tp = T[c]
                    if not np.all(np.isfinite(tp)) or not np.any(tp): continue
                    tg = XP.asarray(tp)
                    st = float(np.sqrt(max(L * np.sum(tp * tp) - np.sum(tp) ** 2, 0.0)))
                    if st <= 0: continue
                    dot = XP.fft.irfft(ch.F * XP.conj(XP.fft.rfft(tg, nfft)), nfft)[:ch.nw]
                    cc = (L * dot - float(np.sum(tp)) * ch.S1) / (st * ch.sd)
                    ok = ch.ok & XP.isfinite(cc)
                    bad = ok & (XP.abs(cc) > 1 + 1e-6)
                    ngt1 += int(bad.sum()); ok = ok & ~bad
                    v = XP.where(ok, cc, 0.0); lv = ok.astype(XP.float64)
                    # SUM[j] += v[j + sh]   (j on the published-time axis from midnight)
                    if sh >= 0:
                        m_ = min(NS, ch.nw - sh)
                        if m_ > 0: SUM[:m_] += v[sh:sh + m_]; LIVE[:m_] += lv[sh:sh + m_]
                    else:
                        s_ = -sh; m_ = min(NS - s_, ch.nw)
                        if m_ > 0: SUM[s_:s_ + m_] += v[:m_]; LIVE[s_:s_ + m_] += lv[:m_]
                    nch += 1; used += 1
                nst += used > 0
            if nch == 0: continue
            SUM = to_np(SUM); LIVE = to_np(LIVE)
            ok = LIVE > 0
            if not ok.any(): continue
            med = np.median(SUM[ok]); mad = np.median(np.abs(SUM[ok] - med)) + 1e-30
            nm = np.where(ok, (SUM - med) / mad, -np.inf)
            keep = decl_bin(nm, a.tmin, gap_n)
            for i in keep:
                rows.append((f, tref + i / FS, float(nm[i]), float(SUM[i]), float(SUM[i] / LIVE[i]), int(LIVE[i]), nch, nst))
            mean = np.where(ok, SUM / np.maximum(LIVE, 1), 0)
            row = dict(fam=f, day=str(d.date()), nchan=nch, nsta=nst, live_med=float(np.median(LIVE[ok])),
                       live_frac=float(ok.mean()), mad=float(mad), n_ccgt1=ngt1,
                       nm_p99=float(np.percentile(nm[ok], 99)), nm_p999=float(np.percentile(nm[ok], 99.9)),
                       nm_max=float(nm[ok].max()), mean_max=float(mean[ok].max()), sum_max=float(SUM[ok].max()))
            for s_ in SWEEP: row[f"n{int(s_)}"] = int(len(decl_bin(nm, s_, gap_n)))
            if TRU is not None:
                th = TRU[(TRU.fam == f) & (TRU.day == d)]
                row["n_theirs"] = int(len(th))
                for e in th.itertuples():
                    j = int(round((e.t.value / 1e9 - tref) * FS))
                    if 0 <= j < NS:
                        lo = max(0, j - int(2 * FS)); hi = min(NS, j + int(2 * FS) + 1)
                        w = np.where(ok[lo:hi], nm[lo:hi], np.nan)
                        tru.append(dict(fam=f, day=str(d.date()), t=e.t, cc_theirs=e.cc,
                                        nm_at=float(nm[j]) if ok[j] else np.nan,
                                        nm_pm2s=float(np.nanmax(w)) if np.isfinite(w).any() else np.nan,
                                        dt_best=float((np.nanargmax(w) + lo - j) / FS) if np.isfinite(w).any() else np.nan,
                                        live_at=float(LIVE[j])))
            diag.append(row)
        del CH
        if GPU: cp.get_default_memory_pool().free_all_blocks()
        if nday % 20 == 0: log(f"  {d.date()} days={nday} det={len(rows)} {time.time() - t00:.0f}s")
    D = pd.DataFrame(rows, columns=["fam", "t_epoch", "net_mad", "ccsum", "meancc", "nlive", "nchan", "nsta"])
    if len(D): D["t"] = pd.to_datetime(D.t_epoch, unit="s").astype("datetime64[ns]")
    D.to_parquet(a.out)
    pd.DataFrame(diag).to_csv(a.out.replace(".parquet", "_diag.csv"), index=False)
    if tru: pd.DataFrame(tru).to_csv(a.out.replace(".parquet", "_truthstat.csv"), index=False)
    dg = pd.DataFrame(diag)
    log(f"[mf-v2] {len(D)} detections on {nday} days in {time.time() - t00:.0f}s -> {a.out}; "
        f"family-days {len(dg)}, |cc|>1 samples {int(dg.n_ccgt1.sum()) if len(dg) else 0}, "
        f"mean_max>1 family-days {int((dg.mean_max > 1).sum()) if len(dg) else 0}")
if __name__ == "__main__": main()
