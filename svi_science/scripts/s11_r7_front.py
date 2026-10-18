#!/usr/bin/env python
"""s11_r7_front -- S11 round 7, task 1 (PLAN.md §11.1): front timing (Merlin test 3).
Day stacks of every H series (round-3 near / round-4 far cc-weighted stores) are measured against the series' grand stack
(2-5 s joint 3C stretch); superposed-epoch dv/v against days since the slip front reached the family (S3 onsets, else the
S3 10 %-quantile definition on the S11 detection days); within-episode shape and absolute curves, bootstrap over majors,
500 front-scrambled draws (front times shifted to other episodes), the time-of-day scrambled control, the bias on ΔD.
    S11_REGION=... python s11_r7_front.py [--workers 12] [--test N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_dist_common as D
PART = sys.argv[sys.argv.index("--part") + 1] if "--part" in sys.argv else "near"
if PART == "far":                                  # far stores live on the far grid (-28..+10 s), as s11_r4_measure --set far
    D.set_grid(D.FAR_GRID)
import s11_r7_common as Q
import s11_r5_common as R5
import s11_r4_common as R4
import s11_r3_common as R3
import s11_r3_inject as J
from s11_measure import Ref3

RES = f"{C.ROOT}/svi_science/results"
WIN = C.WINS[C.HEAD]
BINS = [("le-1", -999, -1), ("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3", 3, 3), ("4-6", 4, 6), ("ge7", 7, 999)]
NBOOT, NSCR = 1000, 500
_G = {}


# ------------------------------------------------------------------------------ day measurement
def measure_series(key):
    st = _G["store"][key]
    St, Ws = st["St"].astype(np.float64), st["Ws"]
    G = np.full((3, St.shape[-1]), np.nan)
    for c in range(3):
        w = Ws[:, c]
        if w.sum() > 0 and np.isfinite(st["G1"][c]).all():
            G[c] = (w @ np.nan_to_num(St[:, c])) / w.sum()
    if not np.isfinite(G).all(1).any():
        return key, None
    Rf = Ref3(G, st["tS"], WIN)
    e, cc = Rf.measure(St, (0, 1, 2))
    out = dict(days=np.asarray(st["days"]).astype(str), eps=e, cc=cc, ws=Ws.mean(1), nd=st["Nd"] if st["Nd"] is not None else np.full(len(e), -1))
    if len(st["Sc"]):
        ec, ccc = Rf.measure(st["Sc"].astype(np.float64), (0, 1, 2))
        out.update(c_days=np.asarray(st["days"])[st["c2_day"]].astype(str), c_eps=ec, c_cc=ccc, c_ws=st["Cn"].mean(1))
    # bin-stack check is done later from the stacks themselves -> keep G and the store reference
    return key, out


def store_paths():
    near = [f"{R3.TB3}/store_ccw"]
    if C.REGION == "SVI" and os.path.isdir(f"{R3.TB3}/store_ccw_ext"):
        near.append(f"{R3.TB3}/store_ccw_ext")
    return near, [f"{R4.TB4}/store_ccw_far"]


def measure_all(H, workers, lg):
    import multiprocessing as mp
    keys = set(H.series)
    near, far = store_paths()
    stas = sorted(H.netsta.unique())
    res = {}
    for ns in stas:
        store = {}
        for sp in near + far:
            f = f"{sp}/{ns}.npz"
            if os.path.exists(f):
                for k, v in J.load_store_file(f).items():
                    if k in keys and k not in store:
                        store[k] = v
        ks = sorted(store)
        if not ks:
            continue
        _G.clear(); _G.update(store=store)
        with mp.get_context("fork").Pool(workers) as p:
            for k, o in p.imap_unordered(measure_series, ks, chunksize=2):
                if o is not None:
                    res[k] = o
        # bin-stack material: keep the day stacks only for the bin check (weighted per bin later) -> store per series
        for k in ks:
            if k in res:
                res[k]["_st"] = dict(St=store[k]["St"], Ws=store[k]["Ws"], G1=store[k]["G1"], tS=store[k]["tS"])
        lg(f"  {ns}: {len(ks)} series measured")
    lg(f"day measurement: {len(res)} of {len(keys)} H series found in the stores")
    return res


# ------------------------------------------------------------------------------ front times
def s11_major_for(name, w0, w1, majors):
    if name in majors:
        return name
    mid = pd.Timestamp(w0) + (pd.Timestamp(w1) - pd.Timestamp(w0)) / 2
    for e in majors:
        if pd.Timestamp(Q.U.start[e]) <= mid <= pd.Timestamp(Q.U.end[e]) + pd.Timedelta(days=1):
            return e
    return None


def fronts_s3(majors, fams, lg):
    """(A) S3 onsets -> {(fam, major): UTC timestamp}"""
    O = pd.read_csv(f"{RES}/s3/family_onsets.csv")
    K = pd.read_csv(f"{RES}/s3/episodes_kinematics.csv")
    EP = pd.read_csv(f"{RES}/s3/episodes_all.csv", parse_dates=["w0", "w1"])
    reg = "SVI" if C.REGION == "SVI" else "NW"
    O = O[O.region == reg]; K = K[K.region == reg].set_index("episode"); EP = EP[EP.region == reg].set_index("episode")
    if reg == "NW":
        src = pd.read_csv(f"{RES}/s3/sources_nw.csv").set_index("source")
        mem = {s: [int(x) for x in str(src.loc[s, "members"]).split(";") if x != "nan"] for s in src.index}
    out, mp_ = {}, {}
    for r in O.itertuples():
        if r.episode not in K.index or r.episode not in EP.index:
            continue
        e = s11_major_for(r.episode, EP.loc[r.episode, "w0"], EP.loc[r.episode, "w1"], majors)
        if e is None:
            continue
        mp_[r.episode] = e
        t = pd.Timestamp(EP.loc[r.episode, "w0"]) + pd.Timedelta(days=float(K.loc[r.episode, "lfe_t_first"]) + float(r.onset))
        fl = [int(r.unit)] if reg == "SVI" else mem.get(r.unit, [])
        for f in fl:
            if f in fams:
                out[(f, e)] = t
    lg(f"S3 fronts: {len(out)} (family, major) pairs over {len(set(e for _, e in out))} majors; S3 -> S11 episode map {mp_}")
    return out


def fronts_own(majors, fams, lg):
    """(B) date where the family's cumulative sumP in the major reaches 10 %"""
    A = pd.read_csv(f"{C.MTB}/active_days.csv")
    if C.REGION == "SVI" and os.path.exists(f"{C.TB}/ext/active_days.csv"):
        A = pd.concat([A, pd.read_csv(f"{C.TB}/ext/active_days.csv")], ignore_index=True).drop_duplicates(["fam", "date"])
    A = A[A.episode.isin(majors) & A.fam.isin(fams)]
    out = {}
    for (f, e), g in A.groupby(["fam", "episode"]):
        g = g.sort_values("date"); cs = np.cumsum(g.sumP.to_numpy()); tot = cs[-1]
        if tot <= 0:
            continue
        i = int(np.argmax(cs >= 0.1 * tot))
        out[(int(f), e)] = pd.Timestamp(g.date.iloc[i])
    lg(f"own fronts: {len(out)} (family, major) pairs over {len(set(e for _, e in out))} majors")
    return out


# ------------------------------------------------------------------------------ statistics
def bin_of(lag):
    b = np.full(len(lag), -1)
    for i, (_, lo, hi) in enumerate(BINS):
        b[(lag >= lo) & (lag <= hi)] = i
    return b


def pooled(D, val, k=None):
    """per-bin weighted means, C = late (lag >= 3) - early (0-1), slope b on lag clipped [-1, 10]; k = bootstrap counts per
    major index (None = 1)"""
    w = D.w.to_numpy() * (1.0 if k is None else k[D.mi.to_numpy()])
    v = val; lag = D.lag.to_numpy()
    out = {}
    b = D.bin.to_numpy()
    for i, (nm, _, _) in enumerate(BINS):
        m = b == i
        out[nm] = float(np.sum(w[m] * v[m]) / np.sum(w[m])) if np.sum(w[m]) > 0 else np.nan
    late = lag >= 3; early = (lag >= 0) & (lag <= 1)
    out["C"] = (float(np.sum(w[late] * v[late]) / np.sum(w[late])) - float(np.sum(w[early] * v[early]) / np.sum(w[early]))
                if np.sum(w[late]) > 0 and np.sum(w[early]) > 0 else np.nan)
    x = np.clip(lag, -1, 10); W = w.sum()
    if W > 0:
        xm = np.sum(w * x) / W; vm = np.sum(w * v) / W
        den = np.sum(w * (x - xm) ** 2)
        out["b"] = float(np.sum(w * (x - xm) * (v - vm)) / den) if den > 0 else np.nan
    else:
        out["b"] = np.nan
    return out


def scr_tables(front_tab, uniq):
    """per (fam, major): its own front offset from the episode start and the offsets of the same family in other majors"""
    by_f = {f: (g.major.to_numpy(), g.rel.to_numpy(float)) for f, g in front_tab.groupby("fam")}
    rel = {(r.fam, r.major): float(r.rel) for r in front_tab.itertuples()}
    r0 = np.array([rel.get(k, np.nan) for k in uniq])
    alts = []
    for f, e in uniq:
        m, rr = by_f.get(f, (np.zeros(0), np.zeros(0)))
        alts.append(rr[m != e])
    return r0, alts


def scr_shift(r0, alts, rng):
    sh = np.zeros(len(r0))
    for i, a in enumerate(alts):
        if len(a) and np.isfinite(r0[i]):
            sh[i] = r0[i] - a[rng.integers(0, len(a))]
    return sh


def shape_vals(D):
    g = D.grp.to_numpy(); w = D.w.to_numpy(); v = D.dvv.to_numpy()
    codes, _ = pd.factorize(g)
    m = np.bincount(codes, w * v) / np.maximum(np.bincount(codes, w), 1e-300)
    return v - m[codes]


def analyse(D, majors, lg, label, front_tab, scr=True, rng_seed=None):
    """D: major-day rows with columns series fam major mi lag bin w dvv base grp; -> summary dict + curve rows"""
    rng = np.random.default_rng(Q.SEED if rng_seed is None else rng_seed)
    D = D[(D.w > 0) & np.isfinite(D.dvv)].reset_index(drop=True)
    if len(D) < 20:
        return None, []
    shp = shape_vals(D)
    has_b = np.isfinite(D.base.to_numpy())
    absv = D.dvv.to_numpy() - np.nan_to_num(D.base.to_numpy())
    S0 = pooled(D, shp); A0 = pooled(D[has_b].reset_index(drop=True), absv[has_b])
    nm = len(majors)
    bs, ba = [], []
    for _ in range(NBOOT):
        k = np.bincount(rng.integers(0, nm, nm), minlength=nm).astype(float)
        bs.append(pooled(D, shp, k)); ba.append(pooled(D[has_b].reset_index(drop=True), absv[has_b], k))
    sig = lambda L, key: float(np.nanstd([x[key] for x in L], ddof=1))
    summ = dict(label=label, n_days=len(D), n_series=D.series.nunique(), n_fam_major=D.grp.nunique(), n_majors=D.major.nunique(),
                C=S0["C"], C_sig=sig(bs, "C"), b=S0["b"], b_sig=sig(bs, "b"), abs_C=A0["C"], abs_C_sig=sig(ba, "C"),
                abs_b=A0["b"], abs_b_sig=sig(ba, "b"))
    rows = []
    for nmb, _, _ in BINS:
        rows.append(dict(label=label, bin=nmb, shape=S0[nmb], shape_sig=sig(bs, nmb), absolute=A0[nmb], absolute_sig=sig(ba, nmb),
                         n_days=int((D.bin == [i for i, x in enumerate(BINS) if x[0] == nmb][0]).sum()),
                         weight=float(D.w[D.bin == [i for i, x in enumerate(BINS) if x[0] == nmb][0]].sum())))
    if scr and front_tab is not None:
        # front-scrambled control: (fam, major) front replaced by start(e) + (front(f, e') - start(e'))
        key_fm = list(zip(D.fam.to_numpy(), D.major.to_numpy()))
        uniq = sorted(set(key_fm)); ui = {k: i for i, k in enumerate(uniq)}; kidx = np.array([ui[k] for k in key_fm])
        r0, alts = scr_tables(front_tab, uniq)
        Cs, bsc, Cabs = [], [], []
        lag0 = D.lag.to_numpy()
        for _ in range(NSCR):
            shift = scr_shift(r0, alts, rng)
            lag_s = lag0 + shift[kidx]
            Ds = D.assign(lag=lag_s, bin=bin_of(lag_s))
            p = pooled(Ds, shp); Cs.append(p["C"]); bsc.append(p["b"])
            pa = pooled(Ds[has_b].reset_index(drop=True), absv[has_b]); Cabs.append(pa["C"])
        Cs, bsc, Cabs = np.array(Cs), np.array(bsc), np.array(Cabs)
        summ.update(scr_p_C=float((np.sum(np.abs(Cs) >= abs(S0["C"])) + 1) / (NSCR + 1)),
                    scr_p_b=float((np.sum(np.abs(bsc) >= abs(S0["b"])) + 1) / (NSCR + 1)),
                    scr_p_absC=float((np.sum(np.abs(Cabs) >= abs(A0["C"])) + 1) / (NSCR + 1)),
                    scr_C_sd=float(np.std(Cs)), scr_C_q95=float(np.percentile(np.abs(Cs), 95)), scr_b_sd=float(np.std(bsc)))
        summ["_scr_C"] = Cs
    lg(f"[{label}] days {len(D)}, series {D.series.nunique()}, majors {D.major.nunique()}: shape C {S0['C']:+.4f} ± {summ['C_sig']:.4f} %, "
       f"b {S0['b']:+.4f} ± {summ['b_sig']:.4f} %/d; absolute C {A0['C']:+.4f} ± {summ['abs_C_sig']:.4f}"
       + (f"; front-scrambled p {summ['scr_p_C']:.3f} (C), {summ['scr_p_b']:.3f} (b)" if "scr_p_C" in summ else ""))
    return summ, rows


def binstack_check(Dh, res, TB, part, lg):
    """headline fronts: per (series, bin) stack the days over majors and measure against the series grand stack"""
    chk = []
    for k, g in Dh.groupby("series"):
        o = res.get(k)
        if o is None:
            continue
        st = o["_st"]; dix = {d: i for i, d in enumerate(o["days"])}
        Gs = np.full((3, st["St"].shape[-1]), np.nan)
        for c in range(3):
            w = st["Ws"][:, c]
            if w.sum() > 0 and np.isfinite(st["G1"][c]).all():
                Gs[c] = (w @ np.nan_to_num(st["St"][:, c].astype(np.float64))) / w.sum()
        if not np.isfinite(Gs).all(1).any():
            continue
        Rf = Ref3(Gs, st["tS"], WIN)
        for b, gb in g.groupby("bin"):
            ii = [dix[d] for d in gb.date if d in dix]
            if len(ii) < 2:
                continue
            Sx = np.full((3, st["St"].shape[-1]), np.nan)
            for c in range(3):
                w = st["Ws"][ii, c]
                if w.sum() > 0:
                    Sx[c] = (w @ st["St"][ii, c].astype(np.float64)) / w.sum()
            e_, c_ = Rf.measure(Sx[None], (0, 1, 2))
            vv = gb[np.isfinite(gb.dvv) & (gb.w > 0)]
            chk.append(dict(series=k, bin=BINS[int(b)][0], n_days=len(ii), stack_dvv=-100 * float(e_[0]) if abs(e_[0]) < C.RAIL else np.nan,
                            stack_cc=float(c_[0]), daypool_dvv=float(np.sum(vv.w * vv.dvv) / vv.w.sum()) if len(vv) else np.nan))
    CK = pd.DataFrame(chk); CK.to_csv(f"{TB}/r7_front_binstack_check_{part}.csv", index=False)
    ok = CK.stack_dvv.notna() & CK.daypool_dvv.notna()
    if ok.sum() > 5:
        lg(f"bin-stack check: {int(ok.sum())} (series, bin) cells, corr {np.corrcoef(CK.stack_dvv[ok], CK.daypool_dvv[ok])[0, 1]:.2f}, "
           f"median (stack − day pool) {np.median(CK.stack_dvv[ok] - CK.daypool_dvv[ok]):+.4f} %")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--test", type=int, default=0)
    ap.add_argument("--part", default="near", choices=["near", "far", "analyse"])
    a = ap.parse_args()
    test = a.test > 0
    TB = Q.tbdir(test); lg = Q.log(f"front_{a.part}" + ("_test" if test else "")); t0 = time.time()
    H, inj = R5.load_H(lg, limit=a.test)
    majors = Q.majors()
    cs, good = R5.setup_cells(H, Q.LSEL)
    ad_ser = cs[good].groupby("series").a_d.median()
    far_ser = H.groupby("series").far.first()
    if a.part == "analyse":
        analyse_part(a, test, TB, lg, t0, H, inj, majors, ad_ser, far_ser); return
    Hp = H[H.far == (a.part == "far")].reset_index(drop=True)
    res = measure_all(Hp, a.workers, lg)
    unit = __import__("s11_paths").unit_map() if C.REGION == "SVI" else dict(zip(pd.read_csv(f"{C.MTB}/active_days.csv").date,
                                                                              pd.read_csv(f"{C.MTB}/active_days.csv").episode))
    rows, crow = [], []
    for k, o in res.items():
        f, ns, ep = k.split("|")
        for kind, dd, e_, cc_, ws_ in (("data", o["days"], o["eps"], o["cc"], o["ws"]),
                                       ("ctrl", o.get("c_days", np.zeros(0)), o.get("c_eps", np.zeros(0)), o.get("c_cc", np.zeros(0)), o.get("c_ws", np.zeros(0)))):
            if not len(dd):
                continue
            ok = np.isfinite(e_) & (np.abs(e_) < C.RAIL)
            df = pd.DataFrame(dict(series=k, fam=int(f), netsta=ns, kind=kind, date=dd, dvv=np.where(ok, -100 * e_, np.nan), cc=cc_,
                                   ws=ws_))
            df["w"] = df.ws * np.clip(np.nan_to_num(df.cc), 0, None) ** 2
            df["episode"] = df.date.map(unit)
            rows.append(df)
    DAY = pd.concat(rows, ignore_index=True)
    DAY["far"] = DAY.series.map(far_ser).astype(bool)
    DAY["a_d"] = DAY.series.map(ad_ser)
    DAY.to_parquet(f"{TB}/r7_front_days_{a.part}.parquet", index=False)
    lg(f"day table ({a.part}): {len(DAY)} rows ({(DAY.kind == 'data').sum()} data, {(DAY.kind == 'ctrl').sum()} control); railed/missing "
       f"{int(DAY.dvv.isna().sum())}; {time.time() - t0:.0f}s")
    front_sets, FT, FA, FB = make_fronts(H, majors, TB, lg)
    Dh = build_days(DAY, front_sets["headline"], majors, None)
    binstack_check(Dh, res, TB, a.part, lg)
    lg(f"done {time.time() - t0:.0f}s")
    if not test:
        open(f"{Q.R7}/FRONT_{C.REGION}_{a.part}_DONE", "w").write(Q.pt_now() + "\n")


def make_fronts(H, majors, TB, lg):
    fams = set(int(x) for x in H.fam.unique())
    FA = fronts_s3(majors, fams, lg); FB = fronts_own(majors, fams, lg)
    front_sets = {"headline": {**FB, **FA}, "S3only": FA, "ownonly": FB}
    ftab = []
    for nm, F in front_sets.items():
        for (f, e), t in F.items():
            ftab.append(dict(set=nm, fam=f, major=e, front_utc=str(t), start=str(Q.U.start[e]),
                             rel=float((pd.Timestamp(t).normalize() - pd.Timestamp(Q.U.start[e])).days)))
    FT = pd.DataFrame(ftab); FT.to_csv(f"{TB}/r7_front_times.csv", index=False)
    both = [k for k in FA if k in FB]
    if both:
        dfd = np.array([(pd.Timestamp(FA[k]).normalize() - pd.Timestamp(FB[k]).normalize()).days for k in both])
        lg(f"S3 vs own front dates on {len(both)} common pairs: median {np.median(dfd):+.1f} d, |Δ| ≤ 1 d for {np.mean(np.abs(dfd) <= 1):.2f}")
    return front_sets, FT, FA, FB


def build_days(DAY, F, majors, sub, kind="data", bser=None):
    mi = {e: i for i, e in enumerate(majors)}
    X = DAY[(DAY.kind == kind) & DAY.episode.isin(majors)].copy()
    if sub is not None:
        X = X[sub(X)]
    fr = pd.to_datetime(pd.Series([F.get((int(f), e), pd.NaT) for f, e in zip(X.fam, X.episode)], dtype="object"), utc=False)
    fr = pd.DatetimeIndex(fr).astype("datetime64[ns]").normalize()
    dd = pd.DatetimeIndex(pd.to_datetime(X.date)).astype("datetime64[ns]")
    X["lag"] = ((dd - fr) / pd.Timedelta(days=1)).to_numpy(float)
    X = X[np.isfinite(X.lag)]
    X["bin"] = bin_of(X.lag.to_numpy()); X["major"] = X.episode; X["mi"] = X.major.map(mi)
    X["base"] = X.series.map(bser) if bser is not None else np.nan; X["grp"] = X.series + "#" + X.major
    return X.reset_index(drop=True)


def analyse_part(a, test, TB, lg, t0, H, inj, majors, ad_ser, far_ser):
    DAY = pd.concat([pd.read_parquet(f"{TB}/r7_front_days_{p}.parquet") for p in ("near", "far") if os.path.exists(f"{TB}/r7_front_days_{p}.parquet")],
                    ignore_index=True)
    lg(f"day table: {len(DAY)} rows ({(DAY.kind == 'data').sum()} data, {(DAY.kind == 'ctrl').sum()} control; far series "
       f"{DAY[DAY.far].series.nunique()}); railed/missing {int(DAY.dvv.isna().sum())}")
    front_sets, FT, FA, FB = make_fronts(H, majors, TB, lg)
    base = DAY[(DAY.kind == "data") & ~DAY.episode.isin(majors) & DAY.dvv.notna() & (DAY.w > 0)]
    bser = (base.w * base.dvv).groupby(base.series).sum() / base.w.groupby(base.series).sum()
    mi = {e: i for i, e in enumerate(majors)}
    summ, curves = [], []
    scr_store = {}

    def build(F, kind="data", sub=None):
        return build_days(DAY, F, majors, sub, kind, bser)

    ftH = FT[FT.set == "headline"][["fam", "major", "rel"]]
    plan = [("headline", "headline", "data", None, True), ("ctrl_timeofday", "headline", "ctrl", None, False),
            ("S3only", "S3only", "data", None, True), ("ownonly", "ownonly", "data", None, False),
            ("near", "headline", "data", lambda X: ~X.far, False), ("far", "headline", "data", lambda X: X.far, False)]
    med_ad = float(np.nanmedian(ad_ser.to_numpy()))
    plan += [("a_d_high", "headline", "data", lambda X: X.a_d >= med_ad, False), ("a_d_low", "headline", "data", lambda X: X.a_d < med_ad, False)]
    Dh = None
    for label, fs, kind, sub, do_scr in plan:
        Dx = build(front_sets[fs], kind, sub)
        if label == "headline":
            Dh = Dx
        ft = FT[FT.set == fs][["fam", "major", "rel"]] if do_scr else None
        s, cr = analyse(Dx, majors, lg, label, ft, scr=do_scr)
        if s is None:
            lg(f"[{label}] too few days"); continue
        if "_scr_C" in s:
            scr_store[label] = s.pop("_scr_C")
        summ.append(s); curves += cr
    S = pd.DataFrame(summ); S.to_csv(f"{TB}/r7_front_summary.csv", index=False)
    pd.DataFrame(curves).to_csv(f"{TB}/r7_front_curves.csv", index=False)
    if scr_store:
        np.savez(f"{TB}/r7_front_scrambled.npz", **{k: v for k, v in scr_store.items()})
    # front-scrambled curves (for the figure): 20 draws of the headline shape curve
    rng = np.random.default_rng(Q.SEED + 7); sc_rows = []
    Dv = Dh[(Dh.w > 0) & np.isfinite(Dh.dvv)].reset_index(drop=True); shp = shape_vals(Dv)
    key_fm = list(zip(Dv.fam.to_numpy(), Dv.major.to_numpy()))
    uniq = sorted(set(key_fm)); ui = {k: i for i, k in enumerate(uniq)}; kidx = np.array([ui[k] for k in key_fm])
    r0, alts = scr_tables(ftH, uniq)
    for dr in range(20):
        lag_s = Dv.lag.to_numpy() + scr_shift(r0, alts, rng)[kidx]
        p = pooled(Dv.assign(lag=lag_s, bin=bin_of(lag_s)), shp)
        sc_rows += [dict(draw=dr, bin=nm, shape=p[nm]) for nm, _, _ in BINS]
    pd.DataFrame(sc_rows).to_csv(f"{TB}/r7_front_scr_curves.csv", index=False)
    CK = pd.concat([pd.read_csv(f"{TB}/r7_front_binstack_check_{p}.csv") for p in ("near", "far")
                    if os.path.exists(f"{TB}/r7_front_binstack_check_{p}.csv")], ignore_index=True)
    CK.to_csv(f"{TB}/r7_front_binstack_check.csv", index=False)
    ok = CK.stack_dvv.notna() & CK.daypool_dvv.notna() if len(CK) else pd.Series(dtype=bool)
    if ok.sum() > 5:
        lg(f"bin-stack check: {int(ok.sum())} (series, bin) cells, corr {np.corrcoef(CK.stack_dvv[ok], CK.daypool_dvv[ok])[0, 1]:.2f}, "
           f"median (stack − day pool) {np.median(CK.stack_dvv[ok] - CK.daypool_dvv[ok]):+.4f} %")
    # bias on ΔD
    src = Q.e4_reference(H, inj, Q.LSEL)
    hd = S[S.label == "headline"].iloc[0]
    b_up = abs(hd.b) + 2 * hd.b_sig
    Dv2 = Dh[(Dh.w > 0) & np.isfinite(Dh.dvv)]
    xl = np.clip(Dv2.lag, -1, 10)
    Lall = float(np.sum(Dv2.w * xl) / Dv2.w.sum())
    brow = []
    gc = []
    for e in majors:
        m = Dv2.major == e
        if not m.any() or e not in src["g"]:
            continue
        Le = float(np.sum(Dv2.w[m] * xl[m]) / Dv2.w[m].sum())
        bias = b_up * (Le - Lall)
        brow.append(dict(major=e, Lbar=Le, Lall=Lall, bias_up=bias, sig_e=src["sig"][e], ratio=abs(bias) / src["sig"][e],
                         dD=src["dD"][e], n_days=int(m.sum()))); gc.append(src["g"][e])
    B = pd.DataFrame(brow)
    gcm = np.mean(gc, 0); sc_ = float(np.sqrt(gcm @ src["F"]["C"] @ gcm))
    comp_bias = float(B.bias_up.mean()) if len(B) else np.nan
    B.to_csv(f"{TB}/r7_front_bias.csv", index=False)
    ctl = S[S.label == "ctrl_timeofday"]
    ctl_ok = bool(len(ctl) and abs(ctl.iloc[0].C) < 2 * ctl.iloc[0].C_sig)
    det = bool(abs(hd.C) >= 3 * hd.C_sig and hd.scr_p_C < 0.01 and ctl_ok)
    V = dict(region=C.REGION, verdict="FRONT-LOCKED CHANGE" if det else "NONE SEEN", C=hd.C, C_sig=hd.C_sig, z=hd.C / hd.C_sig,
             scr_p_C=hd.scr_p_C, b=hd.b, b_sig=hd.b_sig, scr_p_b=hd.scr_p_b, bound2_C=abs(hd.C) + 2 * hd.C_sig, bound2_b=abs(hd.b) + 2 * hd.b_sig,
             abs_C=hd.abs_C, abs_C_sig=hd.abs_C_sig, ctrl_C=float(ctl.iloc[0].C) if len(ctl) else np.nan,
             ctrl_C_sig=float(ctl.iloc[0].C_sig) if len(ctl) else np.nan, ctrl_ok=ctl_ok,
             bias_max_ratio=float(B.ratio.max()) if len(B) else np.nan, bias_comp=comp_bias, sig_comp=sc_,
             bias_comp_ratio=abs(comp_bias) / sc_ if sc_ > 0 else np.nan,
             front_biases_dD=bool(len(B) and (B.ratio.max() >= 0.5 or abs(comp_bias) >= 0.5 * sc_)),
             n_fronts_S3=len(FA), n_fronts_own=len(FB), written=Q.pt_now())
    Q.dump(V, f"{TB}/r7_front_verdict.json")
    lg("VERDICT " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in V.items()))
    lg(f"done {time.time() - t0:.0f}s")
    if not test:
        open(f"{Q.R7}/FRONT_{C.REGION}_DONE", "w").write(Q.pt_now() + "\n")


if __name__ == "__main__":
    main()
