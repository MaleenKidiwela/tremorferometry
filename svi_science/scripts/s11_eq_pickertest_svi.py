#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 2b): the WA picker test repeated on SVI deep earthquakes (calibration region,
where LFE κ 1.773 already agrees with earthquake κ). Records from the 20-Hz 1–8 Hz backward cache (the data the SVI
era stacks came from), upsampled 2/1 to 40 Hz exactly as for the era LFE picks; same pickers, rows, fits.

  python s11_eq_pickertest_svi.py [--workers 8]
writes tomo/eq/tables/pickertest_svi_rows.csv, pickertest_svi_same_events.csv
"""
import argparse
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import resample_poly  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_eq_pickertest as PT  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_era_stacks as ES  # noqa: E402

ROOT = PT.ROOT
EQ = PT.EQ
CACHE = f"{ROOT}/lfe_extension/cache20hz_backward"
STAS = sorted(os.listdir(CACHE))
PRE, LEN = PT.PRE, PT.LEN


def seg(z, t_abs):
    x, t0 = z["x"], float(z["t0"])
    n = int(LEN * 20)
    i0 = int(round((t_abs - t0) * 20))
    if i0 < 0 or i0 + n > x.shape[1]:
        return None
    s = x[:, i0:i0 + n].astype(np.float64)
    if not np.isfinite(s).all() or (s ** 2).sum() == 0:
        return None
    s = s - s.mean(1, keepdims=True)
    return resample_poly(s, 2, 1, axis=-1)


def do_date(args):
    date, g, xy, snr_t, snr_all, cat = args
    ts = pd.Timestamp(date)
    rows = []
    for s in STAS:
        p = f"{CACHE}/{s}/{ts.year}/{ts.dayofyear:03d}.npz"
        if not os.path.exists(p) or s not in xy.index:
            continue
        z = np.load(p)
        order = [str(o) for o in z["order"]]
        if len(order) < 3 or not order[0].endswith("Z"):
            continue
        la, lo = xy.loc[s, "lat"], xy.loc[s, "lon"]
        for _, e in g.iterrows():
            t_or = (e.origin_time - pd.Timestamp("1970-01-01")).total_seconds()
            epi = float(np.hypot((e.o_lon - lo) * 111.19 * np.cos(np.radians(la)), (e.o_lat - la) * 111.19))
            if epi > 100:
                continue
            st = seg(z, t_or - PRE)
            if st is None:
                continue
            nz = seg(z, t_or - 600 - PRE)
            base = dict(reg="SVI", cc_id=e.cc_id, kind=e.kind, mag=e.mag, depth=e.o_depth_km, sta=s, stakey=s,
                        epi_km=epi, cat_tP=cat.get((e.cc_id, s, "P"), np.nan), cat_tS=cat.get((e.cc_id, s, "S"), np.nan))
            pr = PT.pick(st)
            rows.append(dict(base, variant="raw", **{k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}))
            if nz is not None:
                sm, a = PT.noise_match(st, nz, snr_t.get(s, snr_all))
                pr = PT.pick(sm)
                rows.append(dict(base, variant="noise", noise_a=a,
                                 **{k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    rng = np.random.default_rng(13)
    f = f"{EQ}/tables/pickertest_svi_rows.csv"
    if not os.path.exists(f):
        P = pd.read_parquet(f"{ROOT}/data/eq_picks/picks.parquet")
        P = P[(P.reg == "SVI") & P.kind.isin(["layer", "slab"]) & (P.origin_time < "2014-01-01")]
        ev = P.drop_duplicates("cc_id")[["cc_id", "origin_time", "o_lat", "o_lon", "o_depth_km", "kind", "mag"]]
        c = P[P.phase1.isin(["P", "S"])].copy()
        c["man"] = (c["mode"] != "manual").astype(int)
        c = c.sort_values(["man", "t_rel"]).drop_duplicates(["cc_id", "stakey", "phase1"])
        L = pd.read_csv(f"{PT.TT}/picks_EallV3.csv", low_memory=False)
        ok = L[(L.rP >= 1.2) & (L.rS >= 1.2)]
        snr_t = ok.groupby("sta").snr.median().to_dict(); snr_all = float(ok.snr.median())
        xy = ES.station_coords()
        ev["date"] = ev.origin_time.dt.strftime("%Y-%m-%d")
        jobs = []
        for date, g in ev.groupby("date"):
            cc = c[c.cc_id.isin(g.cc_id)]
            jobs.append((date, g, xy, snr_t, snr_all, {(r.cc_id, r.stakey, r.phase1): r.t_rel for r in cc.itertuples()}))
        print("SVI dates", len(jobs), "events", len(ev), flush=True)
        with Pool(a.workers) as pool:
            res = pool.map(do_date, jobs, chunksize=1)
        R = pd.DataFrame([r for rr in res for r in rr])
        R.to_csv(f, index=False)
    R = pd.read_csv(f)
    C = pd.read_csv(f"{EQ}/tables/eq_kappa_rows.csv", low_memory=False)
    PO.NB, PO.NS, PO.EST, PO.NSTA_MIN = 1000, 300, "inv", 4
    out = []

    def fit(q, lab):
        q = q.groupby("fam").filter(lambda g: len(g) >= 4)
        if q.fam.nunique() < 5:
            return dict(label=lab, n_fam=q.fam.nunique(), note="< 5 events")
        o, _ = PO.zone(q, rng, lab, True)
        return dict(label=lab, **{k: o[k] for k in ("n_fam", "n_rows", "kappa", "sig", "rms", "z_vs_scr", "loso_max",
                                                    "d_site", "pass")})
    for var in ("raw", "noise"):
        for pk, tp, rp in (("V3", "tPv3", "rPv3"), ("frozen", "tP0", "rP0")):
            q = R[(R.variant == var) & (R[rp] >= 1.2) & (R.rS >= 1.2)]
            q = pd.DataFrame(dict(fam=q.cc_id, sta=q.sta, tP=q[tp], tS=q.tS)).dropna()
            q = q[q.tS > q.tP]
            o1 = fit(q, f"SVI|{var}|{pk}|picker")
            evs = set(q.groupby("fam").filter(lambda g: len(g) >= 4).fam)
            o2 = fit(C[C.fam.isin(evs) & (C.dist_km <= 100)][["fam", "sta", "tP", "tS"]], f"SVI|{var}|{pk}|catalogue_same_events")
            out += [o1, o2]
            print(o1, "\n", o2, flush=True)
    dq = R[R.variant == "raw"].dropna(subset=["cat_tP"])
    print("raw dtP median", float((dq.tPv3 - dq.cat_tP).median()), "dtS median", float((dq.tS - dq.cat_tS).median()))
    pd.DataFrame(out).to_csv(f"{EQ}/tables/pickertest_svi_same_events.csv", index=False)


if __name__ == "__main__":
    main()
