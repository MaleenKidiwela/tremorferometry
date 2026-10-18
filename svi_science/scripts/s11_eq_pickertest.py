#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 1b): run the LFE picker on earthquake records at the LFE stations and compare
the Wadati κ with the catalogue picks of the same event-station rows.

  python s11_eq_pickertest.py [--nb 1000] [--ns 300]
writes tomo/eq/tables/pickertest_rows.csv, pickertest_kappa.csv, pickertest_bias.csv
"""
import argparse
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import hilbert  # noqa: E402

ROOT = "/home/jovyan/tremorferometry"
EQ = f"{ROOT}/svi_science/results/s11_depth/tomo/eq"
TT = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
_cwd = os.getcwd(); os.chdir(f"{ROOT}/lfe_array_discovery")
import a4f_common as A  # noqa: E402
import a4f_hypo as H  # noqa: E402
os.chdir(_cwd)
import s11_tomo_repick as RP  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402

PRE, LEN = 4.0, 44.0
REG = {"WA_NW": "BV3", "WA_SW": "CV3"}


def env_snr(st):
    env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    md = np.median(env)
    return float(env.max() / md) if md > 0 else np.nan


def noise_match(sig, noise, target):
    """add a*noise (both unit-rms normalised) so that env_snr(sig + a*noise) ~ target (bisection on log a)."""
    s = sig / np.sqrt((sig ** 2).mean()); n = noise / np.sqrt((noise ** 2).mean())
    if not np.isfinite(target) or env_snr(s) <= target:
        return s, 0.0
    lo, hi = -3.0, 3.0
    for _ in range(25):
        m = 0.5 * (lo + hi)
        if env_snr(s + 10 ** m * n) > target:
            lo = m
        else:
            hi = m
    a = 10 ** (0.5 * (lo + hi))
    return s + a * n, a


def pick(st):
    p = H.make_picks({0: {"x": st}})
    if not len(p):
        return {}
    r = p.iloc[0].to_dict()
    out = dict(tS=r.get("tS", np.nan), rS=r.get("rS", np.nan), tP0=r.get("tP", np.nan), rP0=r.get("rP", np.nan),
               snr=r.get("snr", np.nan))
    if np.isfinite(out["tS"]):
        tp, rp, wp = RP.pick_p_v3(st, out["tS"])
        out.update(tPv3=tp, rPv3=rp)
    return out


def _do_date(args):
    reg, date, g, stas, snr_t, snr_all, cat = args
    dc = A.DayCache(stas)
    D = dc.day(date)
    rows = []
    if len(D) < 4:
        return rows
    for _, e in g.iterrows():
        t_or = (e.origin_time - pd.Timestamp("1970-01-01")).total_seconds()
        for s in D:
            st = dc.cut(s, t_or - PRE, LEN, date)
            if st is None:
                continue
            la, lo = A.COORD[s][0], A.COORD[s][1]
            epi = float(np.hypot((e.o_lon - lo) * 111.19 * np.cos(np.radians(la)), (e.o_lat - la) * 111.19))
            if epi > 100:
                continue
            nz = dc.cut(s, t_or - 600 - PRE, LEN, date)
            if nz is None:
                nz = dc.cut(s, t_or + 600, LEN, date)
            key = A.netsta(s)
            cP, cS = cat.get((e.cc_id, key, "P"), np.nan), cat.get((e.cc_id, key, "S"), np.nan)
            base = dict(reg=reg, cc_id=e.cc_id, kind=e.kind, mag=e.mag, depth=e.o_depth_km, sta=s,
                        stakey=key, epi_km=epi, cat_tP=cP, cat_tS=cS)
            pr = pick(st)
            rows.append(dict(base, variant="raw", **{k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}))
            if nz is not None:
                sm, a = noise_match(st, nz, snr_t.get(s, snr_all))
                pr = pick(sm)
                rows.append(dict(base, variant="noise", noise_a=a, snr_target=snr_t.get(s, snr_all),
                                 **{k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}))
    return rows


def build(rng, kinds, workers):
    from multiprocessing import Pool
    P = pd.read_parquet(f"{ROOT}/data/eq_picks/picks.parquet")
    P = P[P.reg.isin(list(REG)) & P.kind.isin(kinds)]
    ev = P.drop_duplicates("cc_id")[["cc_id", "origin_time", "o_lat", "o_lon", "o_depth_km", "reg", "kind", "mag"]]
    c = P[P.phase1.isin(["P", "S"])].copy()
    c["man"] = (c["mode"] != "manual").astype(int)
    c = c.sort_values(["man", "t_rel"]).drop_duplicates(["cc_id", "stakey", "phase1"])
    jobs = []
    for reg, ps in REG.items():
        L = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False)
        ok = L[(L.rP >= 1.2) & (L.rS >= 1.2)]
        snr_t = ok.groupby("sta").snr.median().to_dict()
        snr_all = float(ok.snr.median())
        stas = sorted(L.sta.unique())
        E = ev[ev.reg == reg].copy()
        E["date"] = E.origin_time.dt.strftime("%Y-%m-%d")
        for date, g in E.groupby("date"):
            yd = pd.Timestamp(date)
            have = [s for s in stas if os.path.exists(A.cache_path(s, yd))]
            if len(have) < 4:
                continue
            cc = c[c.cc_id.isin(g.cc_id)]
            cat = {(r.cc_id, r.stakey, r.phase1): r.t_rel for r in cc.itertuples()}
            jobs.append((reg, date, g, have, snr_t, snr_all, cat))
    print("dates to read", len(jobs), "events", sum(len(j[2]) for j in jobs), flush=True)
    with Pool(workers) as pool:
        res = pool.map(_do_date, jobs, chunksize=1)
    return pd.DataFrame([r for rr in res for r in rr])


def fit(Q, tp, ts, rng, lab):
    q = pd.DataFrame(dict(fam=Q.cc_id.values, sta=Q.sta.values, tP=Q[tp].values, tS=Q[ts].values)).dropna()
    q = q[(q.tS > q.tP)]
    q = q.groupby("fam").filter(lambda g: len(g) >= 4)
    if q.fam.nunique() < 5:
        return dict(label=lab, n_ev=q.fam.nunique(), note="< 5 events")
    PO.EST = "inv"; PO.NSTA_MIN = 4
    o, _ = PO.zone(q, rng, lab, True)
    return dict(label=lab, n_ev_in=q.fam.nunique(), **{k: o[k] for k in ("n_fam", "n_rows", "kappa", "sig", "rms",
                                                                       "z_vs_scr", "loso_max", "d_site", "pass")})


def paired_bias(Q, tp, rng, nb):
    """κ(picker) − κ(catalogue) on identical rows, bootstrap over events (no clipping: same rows both)."""
    q = Q.dropna(subset=[tp, "tS", "cat_tP", "cat_tS"])
    q = q[(q.tS > q[tp]) & (q.cat_tS > q.cat_tP)]
    q = q.groupby("cc_id").filter(lambda g: len(g) >= 4)
    if q.cc_id.nunique() < 5:
        return dict(n_ev=q.cc_id.nunique())

    def k(dd, a, b):
        xd, yd, sl = __import__("s11_tomo_profile").pooled_slope(dd.cc_id.values, dd[b].values, dd[a].values)
        return 1 / sl

    # robust: drop rows whose picker S or P is > 1 s from the catalogue in either fit? No — report both
    kp, kc = k(q, tp, "tS"), k(q, "cat_tP", "cat_tS")
    evs = q.cc_id.unique(); idx = {e: np.where(q.cc_id.values == e)[0] for e in evs}
    bs = []
    for _ in range(nb):
        pick_ = rng.choice(evs, len(evs))
        ii = np.concatenate([idx[e] for e in pick_])
        d = q.iloc[ii].copy(); d["cc_id"] = np.concatenate([[f"{j}"] * len(idx[e]) for j, e in enumerate(pick_)])
        bs.append(k(d, tp, "tS") - k(d, "cat_tP", "cat_tS"))
    return dict(n_ev=len(evs), n_rows=len(q), kappa_picker_ols_all=kp, kappa_cat_ols_all=kc, bias_all=kp - kc,
                bias_all_sd=float(np.std(bs)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nb", type=int, default=1000)
    ap.add_argument("--ns", type=int, default=300)
    ap.add_argument("--kinds", default="layer,slab")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    rng = np.random.default_rng(11)
    PO.NB, PO.NS = a.nb, a.ns
    f = f"{EQ}/tables/pickertest_rows.csv"
    if os.path.exists(f):
        R = pd.read_csv(f)
    else:
        R = build(rng, a.kinds.split(","), a.workers)
        R.to_csv(f, index=False)
    R["dtP0"] = R.tP0 - R.cat_tP; R["dtPv3"] = R.tPv3 - R.cat_tP; R["dtS"] = R.tS - R.cat_tS
    K, B = [], []
    for reg in REG:
        for cls, kinds in (("deep", ["layer", "slab"]), ("crust", ["crust"])):
            for var in ("noise", "raw"):
                Q = R[(R.reg == reg) & R.kind.isin(kinds) & (R.variant == var)]
                for pname, tp, rp in (("V3", "tPv3", "rPv3"), ("frozen", "tP0", "rP0")):
                    q = Q[(Q[rp] >= 1.2) & (Q.rS >= 1.2)]
                    lab = f"{reg}|{cls}|{var}|{pname}"
                    o1 = fit(q, tp, "tS", rng, lab + "|picker")
                    o2 = fit(q, "cat_tP", "cat_tS", rng, lab + "|catalogue_same_rows")
                    for o in (o1, o2):
                        o.update(reg=reg, cls=cls, variant=var, picker=pname)
                        K.append(o)
                    b = paired_bias(q, tp, rng, min(a.nb, 500))
                    dq = q.dropna(subset=[tp, "cat_tP"])
                    b.update(reg=reg, cls=cls, variant=var, picker=pname,
                             bias_clip=o1.get("kappa", np.nan) - o2.get("kappa", np.nan),
                             bias_clip_sd=float(np.hypot(o1.get("sig", np.nan), o2.get("sig", np.nan))),
                             med_dtP=float((dq[tp] - dq.cat_tP).median()) if len(dq) else np.nan,
                             frac_dtP_gt05=float(((dq[tp] - dq.cat_tP).abs() > 0.5).mean()) if len(dq) else np.nan,
                             med_dtS=float((q.tS - q.cat_tS).median()) if len(q) else np.nan,
                             frac_dtS_gt05=float(((q.tS - q.cat_tS).abs() > 0.5).mean()) if len(q) else np.nan,
                             n_rows_r12=len(q), n_rows_cat=int(q.cat_tP.notna().sum()))
                    B.append(b)
                    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in b.items()}, flush=True)
        # catalogue picks of the same events at all stations (≤ 100 km) for reference
    pd.DataFrame(K).to_csv(f"{EQ}/tables/pickertest_kappa.csv", index=False)
    pd.DataFrame(B).to_csv(f"{EQ}/tables/pickertest_bias.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(pd.DataFrame(K).round(3).to_string())


if __name__ == "__main__":
    main()
