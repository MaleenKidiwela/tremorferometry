#!/usr/bin/env python
"""p3_quiet_month (2026-10-07) -- P3 step 3: the post-ETS quiet month with three episode-end definitions
(D1 regional end, D2 end + 7 d, D3 the family's own last local-tremor day in that major), on SVI 2014-2026 (S7 tables)
and Washington-NW (item A tables). Bostock 2003-2013 is not testable (ETS windows only; PLAN_P3.md).
    python svi_science/scripts/p3_quiet_month.py
"""
import os, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/p3_bostock"
S0 = f"{ROOT}/svi_science/results/s0"
A_DIR = f"{ROOT}/svi_science/results/a_wa_between"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
NSH = 1000
RNG = np.random.default_rng(20261007)
SVI_P = ("2014-01-01", "2026-05-27")
WA_ERAS = {"E3": ("2015-01-01", "2026-05-27"), "E2": ("2012-01-01", "2014-12-31"), "E1": ("2009-08-06", "2011-12-31")}
LOG = open(f"{OUT}/p3_quiet_month.log", "w")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def pad_mask(dates, E, pad=2):
    d = pd.to_datetime(pd.Series(dates)).values
    m = np.zeros(len(d), bool)
    for s, e in zip(E.start.values, E.end.values):
        m |= (d >= s - np.timedelta64(pad, "D")) & (d <= e + np.timedelta64(pad, "D"))
    return m


def local_ends(maj, X, fids):
    """X: fid, date (str), n30. -> dict fid -> array(len(maj)) of datetime64[D] local ends (regional end if none)."""
    Xe = X[X.n30 > 0]
    out = {}
    me = maj.end.values.astype("datetime64[D]")
    for f in fids:
        d = np.sort(pd.to_datetime(Xe.date[Xe.fid == f]).values.astype("datetime64[D]"))
        le = me.copy()
        for j, (s, e) in enumerate(zip(maj.start.values.astype("datetime64[D]"), me)):
            k = d[(d >= s) & (d <= e)]
            if k.size:
                le[j] = k[-1]
        out[f] = le
    return out


def taus(dates, fids, maj, LE):
    d = pd.to_datetime(pd.Series(dates)).values.astype("datetime64[D]")
    me = maj.end.values.astype("datetime64[D]")
    j = np.searchsorted(me, d, side="left") - 1
    t1 = np.where(j >= 0, (d - me[np.maximum(j, 0)]).astype(int), -1).astype(float); t1[j < 0] = np.nan
    t3 = np.full(d.size, np.nan)
    fa = np.asarray(fids)
    for i in np.flatnonzero(j >= 0):
        t3[i] = (d[i] - LE[fa[i]][j[i]]).astype(int)
    t2 = t1 - 7.0
    t2[t1 < 7] = np.nan
    return t1, t2, t3


def f30_test(ev, EL, col):
    e = ev[np.isfinite(ev[col])]
    pool = {f: g[col].to_numpy() for f, g in EL.groupby("fid")}
    pool = {f: v[np.isfinite(v)] for f, v in pool.items()}
    e = e[e.fid.map(lambda f: f in pool and pool[f].size > 0)]
    if len(e) == 0:
        return dict(n=0)
    obs = float(np.mean(e[col] < 30))
    sh = np.empty((NSH, len(e)))
    for j, f in enumerate(e.fid.to_numpy()):
        p = pool[f]; sh[:, j] = p[RNG.integers(0, p.size, NSH)]
    f30s = np.mean(sh < 30, axis=1)
    return dict(n=int(len(e)), F30=obs, shuffle_F30_median=float(np.median(f30s)), shuffle_F30_p01=float(np.percentile(f30s, 1)),
                p_low=float((np.sum(f30s <= obs) + 1) / (NSH + 1)), n_lt30=int((e[col] < 30).sum()),
                expected_lt30=float(np.median(f30s) * len(e)),
                n_lt14=int((e[col] < 14).sum()), expected_lt14=float(np.median(np.mean(sh < 14, axis=1)) * len(e)))


def verdict(r1, r2, r3):
    ok = lambda r: r.get("n", 0) >= 20 and r["F30"] < r["shuffle_F30_median"] and r["p_low"] < 0.01
    if not ok(r1):
        return "does not hold (fails at the regional end)"
    if ok(r2) and ok(r3):
        return "holds under all three ends"
    return "holds only at the regional end"


def run_set(name, ev, EL, maj, LE, rows):
    for x in (ev, EL):
        x["tau1"], x["tau2"], x["tau3"] = taus(x.date, x.fid, maj, LE)
    rr = {}
    for col, lab in (("tau1", "D1 regional end"), ("tau2", "D2 end + 7 d"), ("tau3", "D3 family's own end")):
        r = f30_test(ev, EL, col); r.update(set=name, end_definition=lab); rows.append(r); rr[col] = r
        log(f"{name} {lab}: {r}")
    v = verdict(rr["tau1"], rr["tau2"], rr["tau3"])
    log(f"{name}: VERDICT {v}")
    d3_minus_d1 = (ev.tau3 - ev.tau1).dropna()
    return v, dict(median_days_family_end_precedes_regional_end=float(d3_minus_d1.median()) if len(d3_minus_d1) else None)


def svi():
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    maj = E[E.major].sort_values("start").reset_index(drop=True)
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8", "n_trem30"])
    T = T[(T.date >= SVI_P[0]) & (T.date <= SVI_P[1])]
    T = T[(T.series == "data") & (T.P_same >= 0.75) & T.keep_xfam8].rename(columns={"fam": "fid"}).copy()
    T["in_pad"] = pad_mask(T.date, E, 2)
    ev = T[~T.in_pad & (T.n_trem30 == 0)][["fid", "date"]].copy()
    evb = T[~T.in_pad][["fid", "date"]].copy()
    D = pd.read_parquet(f"{S0}/famday_detectability.parquet", columns=["run", "fid", "date", "observable"])
    D = D[(D.run == "prod") & D.observable & (D.date >= SVI_P[0]) & (D.date <= SVI_P[1])][["fid", "date"]]
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    D = D.merge(X, on=["fid", "date"], how="left"); D["n30"] = D.n30.fillna(0)
    D["in_pad"] = pad_mask(D.date, E, 2)
    EL = D[~D.in_pad & (D.n30 == 0)][["fid", "date"]].copy(); ELb = D[~D.in_pad][["fid", "date"]].copy()
    log(f"SVI: BE {len(ev)} (S7: 76), BE_broad {len(evb)} (S7: 688), eligible {len(EL)} / {len(ELb)}; majors {len(maj)}")
    LE = local_ends(maj, X, sorted(set(D.fid)))
    return ev, evb, EL, ELb, maj, LE


def wa(era):
    E = pd.read_csv(f"{EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    maj = E[E.major].sort_values("start").reset_index(drop=True)
    a, b = WA_ERAS[era]
    T = pd.read_parquet(f"{A_DIR}/events_all.parquet", columns=["date", "series", "fid", "P_own_WAVE_PNET", "once_D", "day_class"])
    T = T[(T.date >= a) & (T.date <= b) & (T.series == "data") & (T.P_own_WAVE_PNET >= 0.75) & T.once_D].copy()
    T["in_pad"] = pad_mask(T.date, E, 2)
    ev = T[~T.in_pad & (T.day_class == "quiet")][["fid", "date"]].copy()
    evb = T[~T.in_pad & (T.day_class != "unknown")][["fid", "date"]].copy()
    X = pd.read_parquet(f"{A_DIR}/exposure_famday.parquet")
    Xa = X[(X.date >= a) & (X.date <= b)].copy()
    Xa["in_pad"] = pad_mask(Xa.date, E, 2)
    EL = Xa[~Xa.in_pad & (Xa.day_class == "quiet")][["fid", "date"]].copy()
    ELb = Xa[~Xa.in_pad & (Xa.day_class != "unknown")][["fid", "date"]].copy()
    log(f"WA {era}: BE {len(ev)}, BE_broad {len(evb)}, eligible {len(EL)} / {len(ELb)}; majors {len(maj)}")
    LE = local_ends(maj, X[["fid", "date", "n30"]].fillna(0), sorted(set(X.fid)))
    return ev, evb, EL, ELb, maj, LE


def main():
    rows = []; R = {}
    ev, evb, EL, ELb, maj, LE = svi()
    R["SVI_BE"] = run_set("SVI BE (P>=0.75 once, 2014-2026)", ev, EL, maj, LE, rows)
    R["SVI_BE_broad"] = run_set("SVI BE_broad (column)", evb, ELb, maj, LE, rows)
    for era in ("E3", "E2", "E1"):
        ev, evb, EL, ELb, maj, LE = wa(era)
        R[f"WA_{era}_BE"] = run_set(f"WA {era} BE", ev, EL, maj, LE, rows)
        R[f"WA_{era}_BE_broad"] = run_set(f"WA {era} BE_broad (column)", evb, ELb, maj, LE, rows)
    Q = pd.DataFrame(rows)
    Q = Q[["set", "end_definition", "n", "n_lt30", "expected_lt30", "F30", "shuffle_F30_median", "shuffle_F30_p01", "p_low",
           "n_lt14", "expected_lt14"]]
    Q.to_csv(f"{OUT}/quiet_month_three_ends.csv", index=False)
    json.dump(R, open(f"{OUT}/quiet_month_verdicts.json", "w"), indent=1, default=str)
    log("\n" + Q.round(3).to_string(index=False))
    log(json.dumps(R, indent=1, default=str))


if __name__ == "__main__":
    main()
