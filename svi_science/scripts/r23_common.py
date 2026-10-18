"""r23_common -- shared loaders for the R23 studies (rules: svi_science/results/r23/PLAN.md)."""
from __future__ import annotations
import json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
R23 = f"{ROOT}/svi_science/results/r23"
TB = f"{R23}/tables"
FIG = f"{R23}/figures"
S0 = f"{ROOT}/svi_science/results/s0"
PNSN_T0 = "2010-01-06"; P_END = "2026-05-27"
DEEP_GROUPS = (1, 4, 10, 12)
SEG_ORDER = ["CN4f", "ERA", "B2f", "B3f", "S15f", "FWD"]
PERIODS = {"2003-13": ("2003-01-01", "2013-12-31"), "2014-26": ("2014-01-01", P_END), "2003-26": ("2003-01-01", P_END)}


class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, m):
        s = f"[{time.strftime('%H:%M:%S')}] {m}"
        print(s, flush=True); self.f.write(s + "\n"); self.f.flush()


def to_pt(t):
    return pd.to_datetime(t).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")


def families():
    F = pd.read_csv(f"{S0}/families.csv")
    F["deep"] = F.group.isin(DEEP_GROUPS)
    return F.set_index("fid", drop=False)


def episodes():
    E = pd.read_csv(f"{TB}/episodes_2003_2026.csv")
    E["start_ts"] = pd.to_datetime(E.start); E["end_ts"] = pd.to_datetime(E.end)
    return E.sort_values("start").reset_index(drop=True)


def pad_mask(dates, E, pad=2):
    d = pd.to_datetime(pd.Series(np.asarray(dates))).values
    m = np.zeros(len(d), bool)
    for s, e in zip(E.start_ts.values, E.end_ts.values):
        m |= (d >= s - np.timedelta64(pad, "D")) & (d <= e + np.timedelta64(pad, "D"))
    return m


def tremor_famday():
    """(fam, date) -> n30, n30_100, n15 from the source of the date (env < 2010-01-06 <= PNSN <= 2026-05-27)."""
    X = pd.read_parquet(f"{TB}/famday_tremor.parquet")
    main = X[((X.src == "env") & (X.date < PNSN_T0)) | ((X.src == "pnsn") & (X.date >= PNSN_T0))]
    envcol = X[X.src == "env"][["fam", "date", "n30"]].rename(columns={"n30": "n30_env"})
    return main[["fam", "date", "n30", "n30_100", "n15"]], envcol


def attach(T, F, X, Xenv):
    T = T.merge(X, on=["fam", "date"], how="left").merge(Xenv, on=["fam", "date"], how="left")
    for c in ("n30", "n30_100", "n15", "n30_env"):
        T[c] = T[c].fillna(0).astype(int)
    T["dep"] = F.dep.reindex(T.fam).to_numpy(); T["group"] = F.group.reindex(T.fam).to_numpy()
    T["along_km"] = F.along_strike_km.reindex(T.fam).to_numpy(); T["downdip_km"] = F.downdip_km.reindex(T.fam).to_numpy()
    T["deep"] = F.deep.reindex(T.fam).to_numpy()
    return T


def classify(T, E, tremor_end=P_END):
    T = T.copy()
    T["in_pad"] = pad_mask(T.date, E, 2)
    T["in_ep_nopad"] = pad_mask(T.date, E, 0)
    ok = T.date <= tremor_end
    T["BE"] = ~T.in_pad & (T.n30 == 0) & ok
    T["BE_broad"] = ~T.in_pad & ok
    T["BE_envsrc"] = ~T.in_pad & (T.n30_env == 0) & ok        # column: env tremor source for 2010-13 too
    return T


def load_events(F, E):
    V = pd.read_parquet(f"{TB}/events_2003_2026.parquet")
    X, Xe = tremor_famday()
    V = attach(V, F, X, Xe)
    V = classify(V, E)
    V["time_PT"] = to_pt(V.time_utc)
    V["burst_day_f"] = V.burst_day.fillna(False).astype(bool)      # forward rows have no burst-day flag (NA -> False)
    return V


def load_famdays(F, E):
    D = pd.read_parquet(f"{TB}/famdays_observable.parquet")
    D = D[D.date <= P_END]
    X, Xe = tremor_famday()
    D = attach(D, F, X, Xe)
    return classify(D, E)


def tau_clock(dates, starts, ends):
    """days to next start (tau_next) and since previous end (tau_prev) for UTC date strings; index of next start."""
    d = pd.to_datetime(pd.Series(np.asarray(dates))).values.astype("datetime64[D]")
    ms = np.asarray(starts).astype("datetime64[D]"); me = np.asarray(ends).astype("datetime64[D]")
    i = np.searchsorted(ms, d, side="right")
    tn = np.full(len(d), np.nan); ok = i < ms.size
    tn[ok] = (ms[i[ok]] - d[ok]).astype(int)
    j = np.searchsorted(me, d, side="left") - 1
    tp = np.full(len(d), np.nan); okp = j >= 0
    tp[okp] = (d[okp] - me[j[okp]]).astype(int)
    return tn, tp, np.where(ok, i, -1)


def pools(D, mask, by_seg=True):
    """family (and segment) -> array of eligible dates."""
    X = D[mask]
    key = ["fam", "seg"] if by_seg else ["fam"]
    return {k: g.date.to_numpy() for k, g in X.groupby(key)}


def shuffle_dates(keys, pool, rng, n):
    out = np.empty((n, len(keys)), dtype=object)
    for j, k in enumerate(keys):
        p = pool[k]; out[:, j] = p[rng.integers(0, p.size, n)]
    return out


def chain_bursts(t_utc, gap_h=24):
    t = pd.to_datetime(pd.Series(t_utc)).astype("int64").to_numpy() / 3.6e12
    o = np.argsort(t)
    b = np.empty(len(t), int)
    b[o] = np.r_[0, np.cumsum(np.diff(t[o]) >= gap_h)]
    return b


def pct_hi(obs, ctrl):
    ctrl = np.asarray(ctrl, float); return float((np.sum(ctrl >= obs) + 1) / (ctrl.size + 1))


def pct_lo(obs, ctrl):
    ctrl = np.asarray(ctrl, float); return float((np.sum(ctrl <= obs) + 1) / (ctrl.size + 1))


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
