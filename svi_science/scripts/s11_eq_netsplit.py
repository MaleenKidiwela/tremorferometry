#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 5): where does the WA LFE κ excess come from?
(a) κ per station class (permanent UW/US/PB vs temporary XU/TA/YW) + class-label scrambled control + earthquake column;
(b) per-station residuals of tS - tP vs the SVI-calibrated model (κ 1.773, CVM v1.7 rays);
(c) B013 instrument-step check;
(d) S-pick scatter correction for near-station fits (IV headline; E3, Deming, SIMEX columns; synthetic closure;
    random-subset scrambled control).

  OMP_NUM_THREADS=1 python s11_eq_netsplit.py [--workers 3] [--parts a,b,c,d]
writes tomo/eq/tables/netsplit_*.csv
"""
import argparse
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_eq_diff as D  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

ROOT = "/home/jovyan/tremorferometry"
TT = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
EQ = f"{ROOT}/svi_science/results/s11_depth/tomo/eq"
OUT = f"{EQ}/tables"
XY = PR.coords()
AC = sys.modules["a4f_common"]
NET = AC.NET
KCAL, KCAL2 = 1.773, 1.775
EQK = {"WA_NW": 1.775, "WA_SW": 1.735, "SVI": 1.755}
PROF = {"WA_NW": "royer_NW_WA", "WA_SW": "royer_SW_WA", "SVI": "south_VI"}
PICKS = {"WA_NW": "BV3", "WA_SW": "CV3", "SVI": "EallV3"}
PERM_NETS, TEMP_NETS = {"UW", "US", "PB", "CN", "PO", "C8"}, {"XU", "TA", "YW", "XG"}
PO.EST = "inv"


def cls(sta):
    n = sta.split(".")[0] if "." in sta else NET.get(sta, "?")
    return "permanent" if n in PERM_NETS else ("temporary" if n in TEMP_NETS else "?")


def epi(lat, lon, sta):
    if sta not in XY:
        return np.nan
    la, lo = XY[sta][0], XY[sta][1]
    return float(np.hypot((lon - lo) * 111.19 * np.cos(np.radians(lat)), (lat - la) * 111.19))


def ray_tp(reg, depth, dist):
    """CVM v1.7 direct P time; depth clipped to 25-50 km; vectorised per rounded depth."""
    pf = D.load_profile(PROF[reg])
    depth = np.clip(np.asarray(depth, float), 25, 50)
    dist = np.asarray(dist, float)
    out = np.full(len(depth), np.nan)
    zr = np.round(depth * 2) / 2
    for z in np.unique(zr[np.isfinite(zr)]):
        j = np.where((zr == z) & np.isfinite(dist))[0]
        if len(j):
            out[j] = D.ray_partials(pf, float(z), dist[j]).sum(1)
    return out


def load_rows():
    L = pd.read_csv(f"{OUT}/swidth2_lfe_rows.csv")
    L["fam"] = L.fam.astype(str)
    L = L[(L.rP >= 1.2) & (L.rS >= 1.2)].dropna(subset=["tP", "tS", "epi_km"])
    L = L[L.tS > L.tP].copy()
    L["net"] = [s.split(".")[0] if "." in s else NET.get(s, "?") for s in L.sta]
    L["cls"] = L.sta.map(cls)
    L["z"] = np.nan
    for reg, g in L.groupby("reg"):
        L.loc[g.index, "z"] = ray_tp(reg, g.depth_km.values, g.epi_km.values)
    return L.reset_index(drop=True)


def q_of(d):
    return pd.DataFrame(dict(fam=d.fam.values, sta=d.sta.values, tP=d.tP.values, tS=d.tS.values,
                             z=d.z.values if "z" in d else np.nan))


def full_zone(args):
    """PO.zone with gates (bootstrap 1000, scrambled control 300)."""
    label, q, nmin, seed = args
    PO.NB, PO.NS, PO.EST, PO.NSTA_MIN = 1000, 300, "inv", nmin
    q = q.groupby("fam").filter(lambda g: len(g) >= nmin)
    row = dict(label=label, nmin=nmin, n_fam_in=q.fam.nunique(), n_rows_in=len(q), n_sta_in=q.sta.nunique())
    if q.fam.nunique() >= 5:
        o, _ = PO.zone(q, np.random.default_rng(seed), label, True)
        row.update({k: o.get(k) for k in ("n_fam", "n_rows", "kappa", "sig", "sig_fam", "sig_sta", "rms", "z_vs_scr",
                                          "loso_max", "loso_sta", "d_site", "pass")})
        row["kappa_profile"] = prof_kappa(q, nmin)
    return row


def prof_kappa(q, minrows=3):
    """start-free inlier-profile κ (Amendment 7 score, 0.25 s), argmax only."""
    q = q.groupby("fam").filter(lambda g: len(g) >= minrows)
    if q.fam.nunique() < 5:
        return np.nan
    sc = sum(PR.fam_counts(g.tP.values, g.tS.values)[0] for _, g in q.groupby("fam"))
    return float(PR.KG[int(np.argmax(sc))])


def clip_keep(q, nmin):
    PO.NSTA_MIN = nmin
    x, y = q.tS.values, q.tP.values
    b, keep, rms, nf = PO.clip_fit(q.fam.values, x, y)
    return (1 / b if np.isfinite(b) and b != 0 else np.nan), keep


# ------------------------------------------------------------------ estimators on kept rows
def demean(f, *cols):
    df = pd.DataFrame({"f": f, **{f"c{i}": c for i, c in enumerate(cols)}})
    return [(df[f"c{i}"] - df.groupby("f")[f"c{i}"].transform("mean")).values for i in range(len(cols))]


def iv(f, tp, ts, z):
    tPd, tSd, zd = demean(f, tp, ts, z)
    den = (zd * tPd).sum()
    return (zd * tSd).sum() / den if den != 0 else np.nan


def e1(f, tp, ts):
    tPd, tSd = demean(f, tp, ts)
    return (tSd ** 2).sum() / (tSd * tPd).sum()


def boot(Q, fn, rng, nb=1000, by="fam"):
    keys = Q[by].unique()
    idx = {k: np.where(Q[by].values == k)[0] for k in keys}
    out = []
    for _ in range(nb):
        pick = rng.choice(keys, len(keys))
        ii = np.concatenate([idx[k] for k in pick])
        if by == "fam":
            lab = np.concatenate([[j] * len(idx[k]) for j, k in enumerate(pick)])
        else:
            lab = Q.fam.values[ii]
        out.append(fn(lab, Q.tP.values[ii], Q.tS.values[ii], Q.z.values[ii]))
    return np.array(out, float)


def estimators(Q, rng, sigP, sigS, simex=True):
    """Q = kept rows (fam, sta, tP, tS, z)."""
    f, tp, ts, z = Q.fam.values, Q.tP.values, Q.tS.values, Q.z.values
    o = dict(n_fam=Q.fam.nunique(), n_rows=len(Q), n_sta=Q.sta.nunique())
    o["kappa_E1"] = e1(f, tp, ts)
    tPd, tSd, zd = demean(f, tp, ts, z)
    o["kappa_IV"] = iv(f, tp, ts, z)
    o["r_z_tP"] = float(np.corrcoef(zd, tPd)[0, 1]); o["r_z_tS"] = float(np.corrcoef(zd, tSd)[0, 1])
    o["slope_tP_z"] = float((zd * tPd).sum() / (zd ** 2).sum()); o["slope_tS_z"] = float((zd * tSd).sum() / (zd ** 2).sum())
    o["var_w_tS"] = float(np.mean(tSd ** 2)); o["span_tS_med"] = float(Q.groupby("fam").tS.agg(lambda v: v.max() - v.min()).median())
    bf = boot(Q, iv, rng); bs = boot(Q, iv, rng, by="sta")
    o["sig_IV_fam"], o["sig_IV_sta"] = float(np.nanstd(bf)), float(np.nanstd(bs))
    o["sig_IV"] = max(o["sig_IV_fam"], o["sig_IV_sta"])
    be = boot(Q, lambda a, b, c, d: e1(a, b, c), rng, nb=500)
    o["sig_E1_fam"] = float(np.nanstd(be))
    lo = {}
    for st in Q.sta.unique():
        m = Q.sta.values != st
        q2 = Q[m].groupby("fam").filter(lambda g: len(g) >= 2)
        if q2.fam.nunique() >= 3:
            lo[st] = iv(q2.fam.values, q2.tP.values, q2.tS.values, q2.z.values) - o["kappa_IV"]
    if lo:
        smax = max(lo, key=lambda s: abs(lo[s]))
        o.update(loso_IV=float(abs(lo[smax])), loso_IV_sta=smax)
    for tag, s_s in (("", sigS), ("_2s", 2 * sigS)):
        RS = 1 - s_s ** 2 / o["var_w_tS"]
        o[f"R_S{tag}"] = RS
        o[f"kappa_E3{tag}"] = o["kappa_E1"] * RS
    o["kappa_E2"] = 1 / PR.deming(tSd, tPd, sigP ** 2 / sigS ** 2)
    o["sigP"], o["sigS"] = sigP, sigS
    if simex:
        lam = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
        bl = [1 / e1(f, tp, ts)]
        for la in lam[1:]:
            bl.append(np.mean([1 / e1(f, tp, ts + rng.normal(0, np.sqrt(la) * sigS, len(ts))) for _ in range(200)]))
        c = np.polyfit(lam, bl, 2)
        o["kappa_SIMEX"] = float(1 / np.polyval(c, -1.0))
    return o


# ------------------------------------------------------------------ (a) network split
def part_a(L, workers):
    tasks, meta = [], []
    for reg in ("WA_NW", "WA_SW"):
        d = L[L.reg == reg]
        subsets = [("all", d), ("permanent", d[d.cls == "permanent"]), ("temporary", d[d.cls == "temporary"]),
                   ("PB", d[d.net == "PB"]), ("UW_US", d[d.net.isin(["UW", "US"])]), ("XU", d[d.net == "XU"]),
                   ("TA_YW", d[d.net.isin(["TA", "YW"])])]
        # common families: >= 3 rows in both classes
        cp = d[d.cls == "permanent"].fam.value_counts(); ct = d[d.cls == "temporary"].fam.value_counts()
        com = sorted(set(cp.index[cp >= 3]) & set(ct.index[ct >= 3]))
        subsets += [("permanent_common", d[(d.cls == "permanent") & d.fam.isin(com)]),
                    ("temporary_common", d[(d.cls == "temporary") & d.fam.isin(com)]),
                    ("permanent_le50", d[(d.cls == "permanent") & (d.epi_km <= 50)]),
                    ("temporary_le50", d[(d.cls == "temporary") & (d.epi_km <= 50)])]
        for lab, s in subsets:
            for nmin in ((4, 3) if not lab.endswith("le50") else (3,)):
                tasks.append((f"{reg}|{lab}", q_of(s), nmin, 1000 + len(tasks)))
                meta.append(dict(reg=reg, subset=lab, n_common_fam=len(com),
                                 med_epi=float(s.epi_km.median()) if len(s) else np.nan,
                                 stations=",".join(sorted(s.sta.unique()))))
    with Pool(workers) as pool:
        res = pool.map(full_zone, tasks, chunksize=1)
    A = pd.DataFrame([dict(m, **r) for m, r in zip(meta, res)])
    A.to_csv(f"{OUT}/netsplit_class_kappa.csv", index=False)
    print(A[["reg", "subset", "nmin", "n_fam", "n_rows", "kappa", "sig", "z_vs_scr", "loso_max", "d_site", "pass"]]
          .round(3).to_string(), flush=True)
    return A


def scr_task(args):
    reg, d, real, ndraw, seed, strat = args
    rng = np.random.default_rng(seed)
    stas = d.groupby("sta").agg(c=("cls", "first"), e=("epi_km", "median")).reset_index()
    stas = stas[stas.c.isin(["permanent", "temporary"])]
    if strat:
        stas["terc"] = pd.qcut(stas.e, 3, labels=False)
    out = []
    for _ in range(ndraw):
        lab = stas.c.values.copy()
        if strat:
            for t in range(3):
                j = np.where(stas.terc.values == t)[0]; lab[j] = lab[rng.permutation(j)]
        else:
            lab = lab[rng.permutation(len(lab))]
        m = dict(zip(stas.sta, lab))
        kk = []
        for c in ("permanent", "temporary"):
            q = q_of(d[d.sta.map(m) == c]).groupby("fam").filter(lambda g: len(g) >= 4)
            kk.append(clip_keep(q, 4)[0] if q.fam.nunique() >= 5 else np.nan)
        out.append(kk)
    return reg, strat, np.array(out)


def part_a_scr(L, A, workers):
    tasks = []
    for reg in ("WA_NW", "WA_SW"):
        d = L[L.reg == reg]
        for strat in (False, True):
            for k in range(3):
                tasks.append((reg, d, None, 100, 77 + 10 * k + (5 if strat else 0) + (0 if reg == "WA_NW" else 300), strat))
    with Pool(workers) as pool:
        res = pool.map(scr_task, tasks, chunksize=1)
    rows = []
    for reg in ("WA_NW", "WA_SW"):
        g = A[(A.reg == reg) & (A.nmin == 4)].set_index("subset")
        kp, kt = g.loc["permanent", "kappa"], g.loc["temporary", "kappa"]
        for strat in (False, True):
            arr = np.vstack([r[2] for r in res if r[0] == reg and r[1] == strat])
            dk = arr[:, 0] - arr[:, 1]
            dk = dk[np.isfinite(dk)]
            real = kp - kt
            rows.append(dict(reg=reg, scramble="within distance terciles" if strat else "all stations",
                             kappa_perm=kp, kappa_temp=kt, d_real=real, n_draws=len(dk),
                             scr_med=float(np.median(dk)), scr_p025=float(np.percentile(dk, 2.5)),
                             scr_p975=float(np.percentile(dk, 97.5)), real_pct=float(100 * np.mean(dk < real)),
                             scr_perm_med=float(np.nanmedian(arr[:, 0])), scr_temp_med=float(np.nanmedian(arr[:, 1]))))
    S = pd.DataFrame(rows)
    S.to_csv(f"{OUT}/netsplit_class_scrambled.csv", index=False)
    print(S.round(3).to_string(), flush=True)


def part_a_eq(workers):
    """picker κ on WA earthquake records per class (event intercepts)."""
    R = pd.read_csv(f"{OUT}/pickertest_rows.csv")
    R = R[(R.variant == "raw") & (R.epi_km <= 100)].copy()
    R["cls"] = R.sta.map(cls)
    tasks, meta = [], []
    for reg in ("WA_NW", "WA_SW"):
        for pk, tp, rp in (("V3", "tPv3", "rPv3"), ("frozen", "tP0", "rP0"), ("catalogue", "cat_tP", None)):
            d = R[R.reg == reg]
            if rp:
                d = d[(d[rp] >= 1.2) & (d.rS >= 1.2)]
                ts = d.tS
            else:
                ts = d.cat_tS
            q0 = pd.DataFrame(dict(fam=d.cc_id.values, sta=d.sta.values, tP=d[tp].values, tS=ts.values,
                                   cls=d.cls.values)).dropna(subset=["tP", "tS"])
            q0 = q0[q0.tS > q0.tP]
            for c in ("all", "permanent", "temporary"):
                q = q0 if c == "all" else q0[q0.cls == c]
                tasks.append((f"eq|{reg}|{pk}|{c}", q.drop(columns="cls"), 4, 5000 + len(tasks)))
                meta.append(dict(reg=reg, picks=pk, subset=c))
    with Pool(workers) as pool:
        res = pool.map(full_zone, tasks, chunksize=1)
    E = pd.DataFrame([dict(m, **r) for m, r in zip(meta, res)])
    E.to_csv(f"{OUT}/netsplit_eq_class_kappa.csv", index=False)
    print(E[["reg", "picks", "subset", "n_fam", "n_rows", "kappa", "sig", "pass"]].round(3).to_string(), flush=True)


# ------------------------------------------------------------------ (b) residual maps
def part_b(L):
    rows = []
    out = []
    for reg in ("WA_NW", "WA_SW", "SVI"):
        d = L[L.reg == reg].copy()
        _, keep = clip_keep(q_of(d), 4)
        d["kept_allfit"] = keep
        for kc, tag in ((KCAL, ""), (KCAL2, "_k1775")):
            d[f"r_abs{tag}"] = (d.tS - d.tP) - (kc - 1) * d.z
            d[f"r_rel{tag}"] = d[f"r_abs{tag}"] - d.groupby("fam")[f"r_abs{tag}"].transform("median")
        out.append(d)
        for which, dd in (("all_rows", d), ("kept_rows", d[d.kept_allfit])):
            # per-class within-family slope of r_rel vs z  (κ_obs - κ_cal on model times)
            for c in ("all", "permanent", "temporary"):
                q = dd if c == "all" else dd[dd.cls == c]
                q = q.groupby("fam").filter(lambda g: len(g) >= 3)
                if q.fam.nunique() >= 5:
                    zd, rd = demean(q.fam.values, q.z.values, q.r_abs.values)
                    rows.append(dict(reg=reg, rows=which, cls=c, n=len(q), slope_r_vs_z=float((zd * rd).sum() / (zd ** 2).sum()),
                                     med_r_abs=float(q.r_abs.median()), med_r_rel=float(q.r_rel.median())))
    R = pd.concat(out)
    R.to_csv(f"{OUT}/netsplit_residual_rows.csv", index=False)
    st = []
    for (reg, sta), g in R.groupby(["reg", "sta"]):
        k = g[g.kept_allfit]
        st.append(dict(reg=reg, sta=sta, net=NET.get(sta, "?"), cls=cls(sta), lat=XY.get(sta, (np.nan,))[0],
                       lon=XY.get(sta, (np.nan, np.nan))[1], n=len(g), n_kept=len(k),
                       med_epi=float(g.epi_km.median()), med_r_abs=float(g.r_abs.median()),
                       med_r_rel=float(g.r_rel.median()), med_r_rel_kept=float(k.r_rel.median()) if len(k) else np.nan,
                       med_r_abs_kept=float(k.r_abs.median()) if len(k) else np.nan,
                       robsd_r_rel=float(1.4826 * np.median(np.abs(g.r_rel - g.r_rel.median()))),
                       med_r_rel_k1775=float(g.r_rel_k1775.median())))
    S = pd.DataFrame(st)
    S.to_csv(f"{OUT}/netsplit_station_residuals.csv", index=False)
    C = pd.DataFrame(rows)
    C.to_csv(f"{OUT}/netsplit_residual_class.csv", index=False)
    print(C.round(3).to_string(), flush=True)
    print(S[S.reg != "SVI"].sort_values(["reg", "med_epi"]).round(3).to_string(), flush=True)


# ------------------------------------------------------------------ (c) B013
def part_c(L):
    P = pd.read_csv(f"{TT}/picks_BV3.csv", low_memory=False)
    rows = [dict(item="era of the WA NW stack data", value=f"{P.era_start.min()} .. {P.era_end.max()}")]
    C = pd.read_csv(f"{TT}/picks_CV3.csv", low_memory=False)
    rows.append(dict(item="era of the WA SW stack data", value=f"{C.era_start.min()} .. {C.era_end.max()}"))
    d = L[L.reg == "WA_NW"]
    for dmax, nmin in ((1000, 4), (50, 3)):
        base = q_of(d[d.epi_km <= dmax])
        b0 = base.groupby("fam").filter(lambda g: len(g) >= nmin)
        k0, kp0 = clip_keep(b0, nmin)
        rows.append(dict(item=f"κ E1 epi<={dmax} (nmin {nmin}) baseline", value=round(k0, 4), n_kept=int(kp0.sum()),
                         kappa_profile=prof_kappa(b0, nmin),
                         site_share_kept=float(b0[kp0].sta.isin(["B013", "B943", "W030"]).mean())))
        for lab, fn in (("B013 tS -5 ms", lambda q: q.assign(tS=np.where(q.sta == "B013", q.tS - 0.005, q.tS))),
                        ("B013 tS -50 ms", lambda q: q.assign(tS=np.where(q.sta == "B013", q.tS - 0.050, q.tS))),
                        ("B013 dropped", lambda q: q[q.sta != "B013"]),
                        ("B013/B943/W030 site dropped", lambda q: q[~q.sta.isin(["B013", "B943", "W030"])])):
            q = fn(base).groupby("fam").filter(lambda g: len(g) >= nmin)
            k, kp = clip_keep(q, nmin)
            x, y = q.tS.values[kp], q.tP.values[kp]
            rms = float(np.sqrt(np.mean(PO.pooled_fit(q.fam.values[kp], x, y)[1] ** 2))) if kp.sum() > 5 else np.nan
            rows.append(dict(item=f"κ E1 epi<={dmax} (nmin {nmin}) {lab}", value=round(k, 4), delta=round(k - k0, 4),
                             n_kept=int(kp.sum()), rms=rms, kappa_profile=prof_kappa(q, nmin)))
    W = d.pivot_table(index="fam", columns="sta", values=["tP", "tS"])
    for a, b in (("B013", "B943"), ("B013", "W030"), ("B943", "W030")):
        for ph in ("tP", "tS"):
            if (ph, a) in W and (ph, b) in W:
                x = (W[(ph, a)] - W[(ph, b)]).dropna()
                rows.append(dict(item=f"{a}-{b} {ph} (s), same families", value=round(float(x.median()), 4),
                                 robsd=round(float(1.4826 * np.median(np.abs(x - x.median()))), 4), n=len(x)))
    B = pd.DataFrame(rows)
    B.to_csv(f"{OUT}/netsplit_b013.csv", index=False)
    print(B.to_string(), flush=True)


# ------------------------------------------------------------------ (d) scatter correction
DSETS = (("all", 1000, 4), ("le50", 50, 3), ("le50", 50, 4), ("le40", 40, 3), ("le40", 40, 4))


def d_task(args):
    reg, lab, dmax, nmin, d, sigP, sigS, seed = args
    rng = np.random.default_rng(seed)
    q = q_of(d[d.epi_km <= dmax]).groupby("fam").filter(lambda g: len(g) >= nmin)
    row = dict(reg=reg, set=lab, dmax=dmax, nmin=nmin, n_fam_in=q.fam.nunique())
    if q.fam.nunique() < 5:
        return row, None
    k, keep = clip_keep(q, nmin)
    row["kappa_E1_clip"] = k
    row["kappa_profile"] = prof_kappa(q, nmin)
    Q = q[keep].reset_index(drop=True)
    row.update(estimators(Q, rng, sigP, sigS))
    return row, Q


def part_d(L, workers):
    sig = {}
    for reg, ps in PICKS.items():
        P = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False)
        P["fam"] = P.fam.astype(str)
        sig[reg] = PR.colocated_sigma(P, XY)[:2]
    print("sigma P/S", sig, flush=True)
    tasks = []
    for reg in ("WA_NW", "WA_SW", "SVI"):
        d = L[L.reg == reg]
        for lab, dmax, nmin in DSETS:
            tasks.append((reg, lab, dmax, nmin, d, *sig[reg], 300 + len(tasks)))
    with Pool(workers) as pool:
        res = pool.map(d_task, tasks, chunksize=1)
    K = pd.DataFrame([r for r, _ in res])
    kept = {(r["reg"], r["set"], r["nmin"]): Q for r, Q in res if Q is not None}
    # joint family bootstrap of the near - all difference (IV and E1)
    rng = np.random.default_rng(91)
    for i, r in K.iterrows():
        if r.set == "all" or (r.reg, r.set, r.nmin) not in kept:
            continue
        Qn, Qa = kept[(r.reg, r.set, r.nmin)], kept[(r.reg, "all", 4)]
        fams = np.array(sorted(set(Qn.fam) | set(Qa.fam)))
        idn = {f: np.where(Qn.fam.values == f)[0] for f in fams}; ida = {f: np.where(Qa.fam.values == f)[0] for f in fams}
        div, de1 = [], []
        for _ in range(1000):
            pick = rng.choice(fams, len(fams))
            out = []
            for Q, idx in ((Qn, idn), (Qa, ida)):
                ii = np.concatenate([idx[f] for f in pick]); lab = np.concatenate([[j] * len(idx[f]) for j, f in enumerate(pick)])
                if len(np.unique(lab)) < 3:
                    out.append((np.nan, np.nan)); continue
                out.append((iv(lab, Q.tP.values[ii], Q.tS.values[ii], Q.z.values[ii]), e1(lab, Q.tP.values[ii], Q.tS.values[ii])))
            div.append(out[0][0] - out[1][0]); de1.append(out[0][1] - out[1][1])
        a = K.loc[(K.reg == r.reg) & (K.set == "all")].iloc[0]
        K.loc[i, "dIV_vs_all"] = r.kappa_IV - a.kappa_IV
        K.loc[i, "sig_dIV"] = float(np.nanstd(div))
        K.loc[i, "dE1_vs_all"] = r.kappa_E1 - a.kappa_E1
        K.loc[i, "sig_dE1"] = float(np.nanstd(de1))
    # scrambled control: same per-family row counts drawn at random from all that family's rows
    for i, r in K.iterrows():
        if r.set == "all" or r.nmin != 3 or not np.isfinite(r.get("kappa_IV", np.nan)):
            continue
        d = L[L.reg == r.reg]
        q = q_of(d[d.epi_km <= r.dmax]).groupby("fam").filter(lambda g: len(g) >= 3)
        cnt = q.fam.value_counts()
        pool_ = q_of(d[d.fam.isin(cnt.index)]).reset_index(drop=True)
        gi = {f: np.where(pool_.fam.values == f)[0] for f in cnt.index}
        rs = np.random.default_rng(500 + i)
        ks, ke = [], []
        for _ in range(200):
            ii = np.concatenate([rs.choice(gi[f], cnt[f], replace=False) for f in cnt.index])
            s = pool_.iloc[ii].reset_index(drop=True)
            k, keep = clip_keep(s, 3)
            if not np.isfinite(k):
                continue
            Q = s[keep]
            ks.append(iv(Q.fam.values, Q.tP.values, Q.tS.values, Q.z.values)); ke.append(k)
        ks = np.array(ks)
        K.loc[i, "scr_IV_med"] = float(np.median(ks)); K.loc[i, "scr_IV_p16"] = float(np.percentile(ks, 16))
        K.loc[i, "scr_IV_p84"] = float(np.percentile(ks, 84)); K.loc[i, "real_IV_pct_in_scr"] = float(100 * np.mean(ks < r.kappa_IV))
        K.loc[i, "scr_E1_med"] = float(np.median(ke))
    K["eq_kappa"] = K.reg.map(EQK)
    K.to_csv(f"{OUT}/netsplit_scatter_kappa.csv", index=False)
    cols = ["reg", "set", "nmin", "n_fam", "n_rows", "kappa_E1", "kappa_IV", "sig_IV", "loso_IV", "kappa_E3", "kappa_E3_2s",
            "kappa_E2", "kappa_SIMEX", "dIV_vs_all", "sig_dIV", "dE1_vs_all", "sig_dE1", "scr_IV_med", "real_IV_pct_in_scr",
            "r_z_tP", "span_tS_med"]
    print(K[[c for c in cols if c in K]].round(3).to_string(), flush=True)
    # class split with IV (column for part a)
    rows = []
    for reg in ("WA_NW", "WA_SW"):
        d = L[L.reg == reg]
        for c in ("permanent", "temporary"):
            q = q_of(d[d.cls == c]).groupby("fam").filter(lambda g: len(g) >= 4)
            k, keep = clip_keep(q, 4)
            o = dict(reg=reg, cls=c, kappa_E1=k)
            if np.isfinite(k):
                o.update(estimators(q[keep].reset_index(drop=True), np.random.default_rng(7), *sig[reg], simex=False))
            rows.append(o)
    pd.DataFrame(rows).to_csv(f"{OUT}/netsplit_class_iv.csv", index=False)
    print(pd.DataFrame(rows)[["reg", "cls", "kappa_E1", "kappa_IV", "sig_IV", "kappa_E3"]].round(3).to_string(), flush=True)


def synth_task(args):
    reg, kap, draw, g, seed = args
    rows = []
    for lab, dmax, nmin in DSETS:
        q = g[g.epi_km <= dmax][["fam", "sta", "tP", "tS", "z"]].groupby("fam").filter(lambda x: len(x) >= nmin)
        row = dict(reg=reg, kappa_true=kap, draw=draw, set=lab, nmin=nmin, n_fam_in=q.fam.nunique())
        if q.fam.nunique() >= 5:
            k, keep = clip_keep(q, nmin)
            Q = q[keep]
            row.update(kappa_E1=k, n_rows=len(Q), n_fam=Q.fam.nunique(),
                       kappa_IV=iv(Q.fam.values, Q.tP.values, Q.tS.values, Q.z.values),
                       kappa_IV_true_z=iv(Q.fam.values, Q.tP.values, Q.tS.values, Q.z_true.values) if "z_true" in Q else np.nan)
        rows.append(row)
    return rows


def part_d_synth(workers):
    S = pd.read_csv(f"{OUT}/synth_rows.csv")
    S = S[(S.leak == 0.3) & (S.draw >= 1) & (S.rPv3 >= 1.2) & (S.rS >= 1.2)].copy()
    S["fam"] = S.fam.astype(str)
    rng = np.random.default_rng(13)
    parts = []
    for reg, ps in (("WA_NW", "BV3"), ("WA_SW", "CV3")):
        P = pd.read_csv(f"{TT}/picks_{ps}.csv", low_memory=False)[["fam", "lat", "lon", "depth_km"]].drop_duplicates("fam")
        P["fam"] = P.fam.astype(str)
        P["dx"] = rng.normal(0, 5, len(P)); P["dy"] = rng.normal(0, 5, len(P)); P["dz"] = rng.normal(0, 5, len(P))
        g = S[S.reg == reg].merge(P, on="fam")
        g["epi_km"] = [epi(r.lat, r.lon, r.sta) for r in g.itertuples()]
        g["lat_p"] = g.lat + g.dy / 111.19
        g["lon_p"] = g.lon + g.dx / (111.19 * np.cos(np.radians(g.lat)))
        g["epi_p"] = [epi(r.lat_p, r.lon_p, r.sta) for r in g.itertuples()]
        g["z"] = ray_tp(reg, g.depth_km + g.dz, g.epi_p)
        g["z_true"] = g.tp_true
        g["tP"] = g.tPv3
        parts.append(g)
    G = pd.concat(parts)
    G = G[G.tS > G.tP]
    tasks = [(reg, kap, dr, g[["fam", "sta", "tP", "tS", "z", "z_true", "epi_km"]].reset_index(drop=True), 0)
             for (reg, kap, dr), g in G.groupby(["reg", "kappa_true", "draw"])]
    with Pool(workers) as pool:
        res = pool.map(synth_task, tasks, chunksize=1)
    R = pd.DataFrame([r for rr in res for r in rr])
    R.to_csv(f"{OUT}/netsplit_synth_rows_kappa.csv", index=False)
    M = R.groupby(["reg", "kappa_true", "set", "nmin"]).agg(
        kappa_E1=("kappa_E1", "mean"), kappa_E1_sd=("kappa_E1", "std"), kappa_IV=("kappa_IV", "mean"),
        kappa_IV_sd=("kappa_IV", "std"), kappa_IV_true_z=("kappa_IV_true_z", "mean"), n_rows=("n_rows", "mean"),
        n_fam=("n_fam", "mean")).reset_index()
    M["bias_E1"] = M.kappa_E1 - M.kappa_true; M["bias_IV"] = M.kappa_IV - M.kappa_true
    M.to_csv(f"{OUT}/netsplit_synth_kappa.csv", index=False)
    print(M.round(3).to_string(), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--parts", default="a,b,c,d,s")
    a = ap.parse_args()
    parts = a.parts.split(",")
    L = load_rows()
    print("rows", L.groupby(["reg", "cls"]).size().to_dict(), flush=True)
    print("unknown-class stations", sorted(L[L.cls == "?"].sta.unique()), flush=True)
    if "c" in parts:
        part_c(L)
    if "b" in parts:
        part_b(L)
    if "d" in parts:
        part_d(L, a.workers)
    if "s" in parts:
        part_d_synth(a.workers)
    if "a" in parts:
        A = part_a(L, a.workers)
        part_a_scr(L, A, a.workers)
        part_a_eq(a.workers)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
