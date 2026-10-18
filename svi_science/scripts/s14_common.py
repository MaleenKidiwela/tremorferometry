"""s14_common -- S14 (stretching CC as its own observable): loaders, episode classes, flags, Fisher z, REG model.
PLAN: svi_science/results/s14_cc/PLAN.md (written before any CC-vs-time number).
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
RES = f"{ROOT}/svi_science/results"
OUT = f"{RES}/s14_cc"
TB, FIG, LOG = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
S11 = f"{RES}/s11_depth"
EP_FILE = {"SVI": f"{S11}/tables/episode_raw.parquet", "WA_NW": f"{S11}/wa_nw/tables/episode_raw.parquet"}
UNITS = {"SVI": [f"{RES}/s10_dvv/multi/tables/episode_units.csv", f"{S11}/tables/ext/episode_units_new.csv"],
         "WA_NW": [f"{S11}/wa_nw/tables/inv/episode_units.csv"]}
WINDOWS = ["25", "2035", "3560", "36"]           # 3C windows past S (s): 2-5, 2-3.5, 3.5-6, 3-6
PILOT_FAMS = [10, 87, 17]
PILOT_STAS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
COLOC = {"SVI": [["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]],
         "WA_NW": [["PB.B005", "PB.B006", "PB.B007", "UW.LRIV"], ["PB.B013", "PB.B943"]]}
GLITCH = [("CN.PGC", "2019-11"), ("CN.PGC", "2021-01"), ("PB.B010", "2009-03"), ("PB.B010", "2022-03"), ("PB.B010", "2024-08")]
B013_STEP = ("PB.B013", "2015-01-01")


def units(region):
    u = pd.concat([pd.read_csv(f) for f in UNITS[region] if os.path.exists(f)], ignore_index=True)
    u["start"], u["end"] = pd.to_datetime(u.start), pd.to_datetime(u.end)
    u["mid"] = u.start + (u.end - u.start) / 2
    return u.drop_duplicates("episode").set_index("episode")


def fz(cc):
    return np.arctanh(np.clip(np.asarray(cc, float), -0.99, 0.99))


def load_episode(region, pilot=False):
    E = pd.read_parquet(EP_FILE[region]).copy()
    if pilot:
        E = E[E.fam.isin(PILOT_FAMS) & E.netsta.isin(PILOT_STAS)]
    U = units(region)
    E["date"] = pd.to_datetime(E.date_mid).astype("datetime64[ns]")
    E["cls"] = E.episode.map(U.cls).fillna("burst")
    E["series"] = E.fam.astype(str) + "|" + E.netsta + "|" + E.epoch
    E["doy"] = E.date.dt.dayofyear
    E["year"] = E.date.dt.year + (E.doy - 1) / 365.25
    E["dvv25"] = -100 * E.eps3c_25
    E["sig25"] = 100 * E.sig3c_25
    for w in WINDOWS:
        E[f"z_{w}"] = fz(E[f"cc3c_{w}"])
        E[f"cz_{w}"] = fz(E[f"c2_cc3c_{w}"])
    E["z_Z25"], E["z_H25"] = fz(E.ccZ_25), fz(E.ccH_25)
    E["lnn"] = np.log(E.n.clip(lower=1))
    E["lnr"] = np.log(E.pre_rms_ratio.clip(lower=1e-3))
    E["lnd"] = np.log(E.n_days.clip(lower=1))
    E["dcc_shift"] = E.cc2p_25 - E.cc3c_25
    ym = E.date.dt.strftime("%Y-%m")
    g = np.zeros(len(E), bool)
    for ns, m in GLITCH:
        g |= (E.netsta == ns).to_numpy() & (ym == m).to_numpy()
    g |= ((E.netsta == B013_STEP[0]) & (E.date >= pd.Timestamp(B013_STEP[1]))).to_numpy()
    # 2026-10-14: instrument veto windows (svi_science/instrument_vetoes.csv: B013 2012-07-13..2013-07-24, B009
    # 2007-01-08..2009-10-01) are flagged too; the S14 tables on disk were written before this line.
    vf = f"{ROOT}/svi_science/instrument_vetoes.csv"
    if os.path.exists(vf):
        Vt = pd.read_csv(vf); Vt = Vt[Vt.decision.str.strip() == "veto"]
        for r in Vt.itertuples():
            g |= ((E.netsta == r.netsta) & (E.date >= pd.Timestamp(r.start)) & (E.date <= pd.Timestamp(r.end))).to_numpy()
    E["glitch"] = g
    E["shift_flag"] = E.shift2p_25.abs() >= 0.02
    # position relative to majors (for C2): days after the nearest earlier major end / before the next major start
    M = U[U.cls == "major"].sort_values("start")
    ends, starts = M.end.to_numpy(), M.start.to_numpy()
    d = E.date.to_numpy()
    i_prev = np.searchsorted(ends, d, side="right") - 1
    i_next = np.searchsorted(starts, d, side="left")
    E["d_after_major"] = np.where(i_prev >= 0, (d - ends[np.clip(i_prev, 0, None)]) / np.timedelta64(1, "D"), np.nan)
    E["d_before_major"] = np.where(i_next < len(starts), (starts[np.clip(i_next, None, len(starts) - 1)] - d) / np.timedelta64(1, "D"), np.nan)
    return E.reset_index(drop=True)


def demean_fe(y, groups):
    """subtract group means (series fixed effects) from columns of y"""
    y = np.asarray(y, float)
    df = pd.DataFrame(y if y.ndim == 2 else y[:, None])
    return (df - df.groupby(np.asarray(groups)).transform("mean")).to_numpy()


def reg_fit(D, ycol, covs, terms=None, fe="series"):
    """ANCOVA with fixed effects: y = a_fe + b.covs + g.terms. Returns dict(coef, se, resid, names). OLS, FE by demeaning."""
    cols = list(covs) + list(terms or [])
    ok = np.isfinite(D[ycol].to_numpy())
    for c in cols:
        ok &= np.isfinite(D[c].to_numpy(float))
    Dk = D[ok]
    if len(Dk) < len(cols) + 5:
        return None
    X = demean_fe(Dk[cols].to_numpy(float), Dk[fe]) if cols else np.zeros((len(Dk), 0))
    y = demean_fe(Dk[ycol].to_numpy(float), Dk[fe])[:, 0]
    keep = np.ones(X.shape[1], bool)
    keep &= X.std(0) > 1e-12
    Xk = X[:, keep]
    b, *_ = np.linalg.lstsq(Xk, y, rcond=None)
    r = y - Xk @ b
    nfe = Dk[fe].nunique()
    dof = max(len(y) - Xk.shape[1] - nfe, 1)
    s2 = r @ r / dof
    cov = s2 * np.linalg.pinv(Xk.T @ Xk)
    names = [c for c, k in zip(cols, keep) if k]
    resid = pd.Series(np.nan, index=D.index); resid[Dk.index] = r
    return dict(coef=dict(zip(names, b)), se=dict(zip(names, np.sqrt(np.diag(cov)))), resid=resid, n=len(y), names=names,
                Xk=Xk, y=y, idx=Dk.index)


def gamma_fast(Xc, y, t):
    """coefficient of one demeaned term t given demeaned covariates Xc (Frisch-Waugh)"""
    if Xc.shape[1]:
        bt, *_ = np.linalg.lstsq(Xc, t, rcond=None); tr = t - Xc @ bt
        by, *_ = np.linalg.lstsq(Xc, y, rcond=None); yr = y - Xc @ by
    else:
        tr, yr = t, y
    den = tr @ tr
    return np.nan if den <= 0 else float(tr @ yr / den)


def season_block_perm(D, mask_cols, rng, half=45):
    """shuffle the (boolean) label columns among cells of the same series within +-half days of year.
    Implementation: within each series, sort by doy, shuffle labels inside consecutive 2*half-day doy blocks with a random
    phase offset per draw."""
    out = {c: D[c].to_numpy().copy() for c in mask_cols}
    off = rng.integers(0, 2 * half)
    blk = ((D.doy.to_numpy() + off) // (2 * half)).astype(int)
    key = pd.factorize(D.series.astype(str) + "#" + pd.Series(blk, index=D.index).astype(str))[0]
    order = np.argsort(key, kind="stable")
    ks = key[order]
    bounds = np.flatnonzero(np.diff(ks)) + 1
    for seg in np.split(order, bounds):
        if len(seg) < 2:
            continue
        p = rng.permutation(len(seg))
        for c in mask_cols:
            out[c][seg] = out[c][seg[p]]
    return out


def pt_now():
    import datetime as dt, zoneinfo
    return dt.datetime.now(zoneinfo.ZoneInfo("America/Los_Angeles")).strftime("%Y-%m-%d %-I:%M %p PT")
