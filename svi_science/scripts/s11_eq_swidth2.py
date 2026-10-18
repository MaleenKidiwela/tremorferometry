#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 4a, 4c): S/P envelope width vs distance on the LFE stacks AND on earthquake
records (same band, same picker), scrambled controls (distances shuffled within family/event), and the
distance-shuffle scrambled control for the distance-limited LFE κ.

  python s11_eq_swidth2.py [--workers 4]
writes tomo/eq/tables/swidth2_lfe_rows.csv, swidth2_eq_rows.csv, swidth2_slopes.csv, swidth2_dist_kappa.csv
"""
import argparse
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import hilbert, resample_poly, butter, sosfiltfilt  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_eq_pickertest as PT  # noqa: E402
import s11_eq_pickertest_svi as PS  # noqa: E402
import s11_eq_swidth as SW  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

A = PT.A
ROOT, EQ, TT = PT.ROOT, PT.EQ, PT.TT
FS = 40.0
PRE, LEN = PT.PRE, PT.LEN
SOS28 = butter(4, [2.0, 8.0], btype="band", fs=FS, output="sos")


def measure(st, tP, tS, off=0.0):
    """widths for one 3C record; tP/tS in the record's time frame (s from sample 0 = off)."""
    out = {}
    for band, x in (("", st), ("_28", sosfiltfilt(SOS28, st, axis=-1))):
        eh = np.sqrt(np.abs(hilbert(x[1])) ** 2 + np.abs(hilbert(x[2])) ** 2)
        ez = np.abs(hilbert(x[0]))
        wS, onS, pkS = SW.width(eh, (tS - off) * FS) if np.isfinite(tS) else (np.nan,) * 3
        wP, onP, pkP = SW.width(ez, (tP - off) * FS) if np.isfinite(tP) else (np.nan,) * 3
        out.update({f"wS{band}": wS, f"wP{band}": wP, f"riseS{band}": pkS - onS, f"riseP{band}": pkP - onP,
                    f"pkmS{band}": tS - (pkS + off), f"pkmP{band}": tP - (pkP + off),
                    f"tS_on{band}": onS + off, f"tP_on{band}": onP + off})
    return out


# ------------------------------------------------------------------ LFE stacks
def lfe_rows():
    xy = PR.coords()
    out, cache = [], {}
    for reg, ps in SW.SETS.items():
        P = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False)
        P = P[P.tS.notna() & (P.rS >= 1.2)].copy()
        if reg != "SVI":
            P["epi_km"] = [np.hypot((r.lon - xy[r.sta][1]) * 111.19 * np.cos(np.radians(r.lat)),
                                    (r.lat - xy[r.sta][0]) * 111.19) if r.sta in xy else np.nan for r in P.itertuples()]
        off = P.frame_offset.values if reg == "SVI" else np.zeros(len(P))
        for r, o in zip(P.itertuples(), off):
            if reg == "SVI":
                p = f"{TT}/era_stacks/{r.sta}.npz"
                if p not in cache:
                    cache[p] = np.load(p)
                st = resample_poly(cache[p][str(int(r.grp))].astype(np.float64), 2, 1, axis=-1)
            else:
                p = f"{ROOT}/{r.src}"
                if p not in cache:
                    cache[p] = np.load(p)
                st = cache[p][f"{int(r.grp)}_{r.sta}"].astype(np.float64)
            out.append(dict(reg=reg, fam=str(r.fam), sta=r.sta, epi_km=r.epi_km, depth_km=r.depth_km, snr=r.snr,
                            rP=r.rP, rS=r.rS, tP=r.tP, tS=r.tS, **measure(st, r.tP, r.tS, o)))
        print("LFE", reg, len(P), flush=True)
    return pd.DataFrame(out)


# ------------------------------------------------------------------ earthquake records
def _wa_date(args):
    date, g = args
    dc = A.DayCache(sorted(g.sta.unique()))
    dc.day(date)
    rows = []
    for r in g.itertuples():
        st = dc.cut(r.sta, r.t_or - PRE, LEN, date)
        if st is None:
            continue
        rows.append(dict(r._asdict(), **measure(st, r.tPv3, r.tS, -PRE)))
    return rows


def _svi_date(args):
    date, g = args
    ts = pd.Timestamp(date)
    rows = []
    for sta, gg in g.groupby("sta"):
        p = f"{PS.CACHE}/{sta}/{ts.year}/{ts.dayofyear:03d}.npz"
        if not os.path.exists(p):
            continue
        z = np.load(p)
        for r in gg.itertuples():
            st = PS.seg(z, r.t_or - PRE)
            if st is None:
                continue
            rows.append(dict(r._asdict(), **measure(st, r.tPv3, r.tS, -PRE)))
    return rows


def eq_rows(workers):
    E = pd.read_parquet(f"{ROOT}/data/eq_picks/picks.parquet", columns=["cc_id", "origin_time"]).drop_duplicates("cc_id")
    E["t_or"] = (E.origin_time - pd.Timestamp("1970-01-01")).dt.total_seconds()
    out = []
    for f, fn in (("pickertest_rows", _wa_date), ("pickertest_svi_rows", _svi_date)):
        R = pd.read_csv(f"{EQ}/tables/{f}.csv")
        R = R[(R.variant == "raw") & (R.rS >= 1.2)].drop(columns=[c for c in ("noise_a", "snr_target") if c in R])
        R = R.merge(E[["cc_id", "t_or", "origin_time"]], on="cc_id")
        R["date"] = R.origin_time.dt.strftime("%Y-%m-%d")
        R = R.drop(columns=["origin_time"])
        jobs = [(d, g) for d, g in R.groupby("date")]
        with Pool(workers) as pool:
            res = pool.map(fn, jobs, chunksize=1)
        out += [r for rr in res for r in rr]
        print("EQ", f, len(R), "measured", sum(len(r) for r in res), flush=True)
    Q = pd.DataFrame(out).drop(columns=["Index"], errors="ignore")
    Q["dtS_cat"] = Q.tS - Q.cat_tS
    Q["dtP_cat"] = Q.tPv3 - Q.cat_tP
    return Q


# ------------------------------------------------------------------ slopes + scrambled control
def slope_stats(d, col, gcol, rng, nb=1000, ns=1000):
    d = d.dropna(subset=[col, "epi_km"])
    if len(d) < 15 or d[gcol].nunique() < 3:
        return dict(n=len(d))
    x, y, g = d.epi_km.values, d[col].values, d[gcol].values
    b = np.polyfit(x, y, 1)
    groups = np.unique(g); idx = {k: np.where(g == k)[0] for k in groups}
    bs = []
    for _ in range(nb):
        ii = np.concatenate([idx[k] for k in rng.choice(groups, len(groups))])
        bs.append(np.polyfit(x[ii], y[ii], 1)[0])
    sc = []
    for _ in range(ns):                           # scrambled control: distances shuffled within each group
        xs = x.copy()
        for k, j in idx.items():
            xs[j] = x[rng.permutation(j)]
        sc.append(np.polyfit(xs, y, 1)[0])
    sc = np.array(sc)
    rho, p = spearmanr(x, y)
    return dict(n=len(d), n_groups=len(groups), slope_ms_per_km=1000 * b[0], intercept_s=b[1],
                slope_sd=1000 * float(np.std(bs)), z=b[0] / np.std(bs) if np.std(bs) > 0 else np.nan, rho=rho, p_rho=p,
                scr_slope_med=1000 * float(np.median(sc)), scr_slope_sd=1000 * float(np.std(sc)),
                p_scr=float((np.abs(sc - np.median(sc)) >= abs(b[0] - np.median(sc))).mean()),
                med_le30=float(np.median(y[x <= 30])) if (x <= 30).any() else np.nan,
                med_30_50=float(np.median(y[(x > 30) & (x <= 50)])) if ((x > 30) & (x <= 50)).any() else np.nan,
                med_ge60=float(np.median(y[x >= 60])) if (x >= 60).any() else np.nan, med_all=float(np.median(y)))


def slopes(L, Q, rng):
    S = []
    for src, D, gcol, rpcol in (("LFE", L, "fam", "rP"), ("EQ", Q, "cc_id", "rPv3")):
        for reg in ("WA_NW", "WA_SW", "SVI"):
            d = D[D.reg == reg].copy()
            for band in ("", "_28"):
                d[f"wSmP{band}"] = d[f"wS{band}"] - d[f"wP{band}"]
                cols = [(f"wS{band}", d.rS >= 1.2), (f"wP{band}", d[rpcol] >= 1.2),
                        (f"wSmP{band}", (d.rS >= 1.2) & (d[rpcol] >= 1.2)), (f"riseS{band}", d.rS >= 1.2),
                        (f"riseP{band}", d[rpcol] >= 1.2), (f"pkmS{band}", d.rS >= 1.2)]
                if src == "EQ" and band == "":
                    cols += [("dtS_cat", d.rS >= 1.2), ("dtP_cat", d[rpcol] >= 1.2)]
                for col, m in cols:
                    o = dict(src=src, reg=reg, col=col, **slope_stats(d[m], col, gcol, rng))
                    S.append(o)
                    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items()}, flush=True)
    return pd.DataFrame(S)


# ------------------------------------------------------------------ distance-limited κ + scrambled control
def dist_kappa(L, rng, ndraw=300):
    PO.NB, PO.NS, PO.EST = 1000, 300, "inv"
    K = []
    for reg in ("WA_NW", "WA_SW", "SVI"):
        d = L[(L.reg == reg) & (L.rP >= 1.2) & (L.rS >= 1.2)].copy()
        d["hyp_km"] = np.hypot(d.epi_km, d.depth_km)
        base = pd.DataFrame(dict(fam=d.fam.values, sta=d.sta.values, tP=d.tP.values, tS=d.tS.values,
                                 epi=d.epi_km.values, hyp=d.hyp_km.values)).dropna(subset=["tP", "tS", "epi"])
        base = base[base.tS > base.tP]
        for dcol, dmax, nmin in (("epi", 40, 3), ("epi", 40, 4), ("epi", 50, 3), ("epi", 50, 4), ("hyp", 50, 3),
                                 ("hyp", 60, 3), ("hyp", 60, 4), ("epi", 1000, 4)):
            PO.NSTA_MIN = nmin
            q = base[base[dcol] <= dmax].groupby("fam").filter(lambda g: len(g) >= nmin)
            row = dict(reg=reg, dist=dcol, dmax=dmax, nmin=nmin, n_fam_in=q.fam.nunique(), n_rows_in=len(q))
            if q.fam.nunique() >= 5:
                o, _ = PO.zone(q, rng, f"{reg}|{dcol}{dmax}", True)
                row.update({k: o[k] for k in ("n_fam", "n_rows", "kappa", "sig", "rms", "z_vs_scr", "loso_max",
                                              "d_site", "pass")})
                row["med_dist_used"] = float(q[dcol].median())
                if dmax < 1000:
                    # scrambled control: same per-family row counts drawn at random from all that family's rows
                    cnt = q.fam.value_counts()
                    pool = base[base.fam.isin(cnt.index)]
                    gi = {f: np.where(pool.fam.values == f)[0] for f in cnt.index}
                    ks, md = [], []
                    for _ in range(ndraw):
                        ii = np.concatenate([rng.choice(gi[f], cnt[f], replace=False) for f in cnt.index])
                        s = pool.iloc[ii]
                        x, y = PO.xy(s.tP.values, s.tS.values)
                        ks.append(PO.kap(PO.clip_fit(s.fam.values, x, y)[0]))
                        md.append(float(s[dcol].median()))
                    ks = np.array(ks, float)
                    row.update(scr_kappa_med=float(np.nanmedian(ks)), scr_kappa_p16=float(np.nanpercentile(ks, 16)),
                               scr_kappa_p84=float(np.nanpercentile(ks, 84)),
                               real_pct_in_scr=float(100 * np.nanmean(ks < row["kappa"])),
                               scr_med_dist=float(np.median(md)), scr_frac_measured=float(np.isfinite(ks).mean()))
            K.append(row)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    return pd.DataFrame(K)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    rng = np.random.default_rng(41)
    fl, fe = f"{EQ}/tables/swidth2_lfe_rows.csv", f"{EQ}/tables/swidth2_eq_rows.csv"
    if not os.path.exists(fl):
        lfe_rows().to_csv(fl, index=False)
    if not os.path.exists(fe):
        eq_rows(a.workers).to_csv(fe, index=False)
    L, Q = pd.read_csv(fl), pd.read_csv(fe)
    L["fam"] = L.fam.astype(str)
    slopes(L, Q, rng).to_csv(f"{EQ}/tables/swidth2_slopes.csv", index=False)
    dist_kappa(L, rng).to_csv(f"{EQ}/tables/swidth2_dist_kappa.csv", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
