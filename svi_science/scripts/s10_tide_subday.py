#!/usr/bin/env python
"""s10_tide_subday -- S10 follow-up part B (rules: results/s10_dvv/tide_subday/PLAN.md, B1-B4).
6-h (headline) and 12-h stacks inside the nine biggest episodes, per (family, station, epoch); each bin with >= 20
repeats is measured against the leave-bin-out epoch grand stack (2-4 and 3-5 s). Front arrival: S3 tremor front at
the family where it exists, else the family's own onset (10th-percentile time of its high-P repeats). Change tests:
after-front minus before-front and last third minus first third, with a within-(series, episode) time-shuffled control.
    python s10_tide_subday.py
"""
from __future__ import annotations
import os, sys, time, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_measure as M
import s10_tide_common as T

EPISODES = ["2004-07", "2005-09", "2003-02", "2013-09", "2008-05", "2011-08", "2024-08", "2023-04", "2019-08"]
RES = {"6h": 6, "12h": 12}
NMIN = 20
NSHUF = 2000
SEED = 20261012


def fronts(R, EPI):
    S3 = pd.read_csv(f"{S.ROOT}/svi_science/results/s5/step3_front_arrivals.csv")
    PF = pd.read_csv(f"{T.PILOT_TB}/front_dates.csv")
    rows = []
    for fam in T.FAMS:
        for ep in EPISODES:
            e = EPI[EPI.episode == ep].iloc[0]
            a = pd.Timestamp(e.start) - pd.Timedelta("2D"); b = pd.Timestamp(e.end) + pd.Timedelta("3D")
            m = (R.fam == fam) & (R.series == "data") & R.hiP & (R.tt >= a) & (R.tt < b)
            onset = R.tt[m].quantile(0.10) if m.sum() else pd.NaT
            s3 = S3[(S3.fam == fam) & (S3.episode == ep)]
            pf = PF[(PF.fam == fam) & (PF.episode == ep)]
            tf = pd.Timestamp(s3.t_front.iloc[0]) if len(s3) else onset
            rows.append(dict(fam=fam, episode=ep, win_start=a, win_end=b, n_hiP=int(m.sum()), onset_p10=onset,
                             s3_front=pd.Timestamp(s3.t_front.iloc[0]) if len(s3) else pd.NaT,
                             front=tf, front_src="S3 tremor front" if len(s3) else "family onset (10th pct high-P)",
                             pilot_front_date=pf.front.iloc[0] if len(pf) else ""))
    return pd.DataFrame(rows)


def measure_bins(R, FR):
    pos = R.set_index("row_id")
    out = []
    for f in sorted(os.listdir(T.CACHE)):
        z = np.load(f"{T.CACHE}/{f}")
        fam, ns, ep = int(z["fam"]), str(z["netsta"]), str(z["epoch"])
        t = pd.to_datetime(pos.loc[z["row_id"], "t"].to_numpy(), unit="s")
        t2 = pd.to_datetime(pos.loc[z["row_id2"], "t"].to_numpy(), unit="s") if len(z["row_id2"]) else pd.DatetimeIndex([])
        Z = z["Z"].astype(np.float64); P = z["P"].astype(np.float64); tS = float(z["tS"])
        Gs = P @ Z; Pt = P.sum()
        Z2 = z["Z2"].astype(np.float64)
        for _, fr in FR[FR.fam == fam].iterrows():
            me = (t >= fr.win_start) & (t < fr.win_end)
            if me.sum() < NMIN:
                continue
            for rn, hh in RES.items():
                bins = t[me].floor(f"{hh}h")
                idx_e = np.flatnonzero(me)
                for bstart in np.unique(bins):
                    ib = idx_e[bins == bstart]
                    if len(ib) < NMIN:
                        continue
                    ref = (Gs - P[ib] @ Z[ib]) / (Pt - P[ib].sum())
                    cur = (P[ib] @ Z[ib]) / P[ib].sum()
                    r = dict(fam=fam, netsta=ns, epoch=ep, episode=fr.episode, res=rn, bin_start=pd.Timestamp(bstart),
                             bin_mid=pd.Timestamp(bstart) + pd.Timedelta(hours=hh / 2), n=len(ib), sumP=P[ib].sum(),
                             t_mean=t[ib].mean(), front=fr.front, front_src=fr.front_src)
                    for k, win in T.WINS.items():
                        Rf = M.Ref(ref, tS, win)
                        e, c = Rf.measure(cur)
                        r[f"eps{k}_raw"] = float(e[0]); r[f"cc{k}"] = float(c[0])
                        r[f"eps{k}"] = float(e[0]) if abs(e[0]) < T.RAIL else np.nan
                        if len(t2):
                            m2 = np.flatnonzero((t2 >= pd.Timestamp(bstart)) & (t2 < pd.Timestamp(bstart) + pd.Timedelta(hours=hh)))
                            r["n_c2"] = len(m2)
                            if len(m2) >= 3:
                                e2, c2 = Rf.measure(Z2[m2].mean(0))
                                r[f"c2_eps{k}"] = float(e2[0]) if abs(e2[0]) < T.RAIL else np.nan
                                r[f"c2_cc{k}"] = float(c2[0])
                    out.append(r)
    return pd.DataFrame(out)


def _wm(e, w):
    return np.sum(w * e) / np.sum(w) if np.sum(w) > 0 else np.nan


def tests(B, rng):
    rows = []
    for rn in RES:
        for famset in [(10,), (87,), (17,), (10, 87, 17)]:
            Bs = B[(B.res == rn) & B.fam.isin(famset)].sort_values("bin_mid")
            for k in T.WINS:
                units = []
                for key, g in Bs.groupby(["fam", "netsta", "epoch", "episode"]):
                    g = g.sort_values("bin_mid")
                    e = g[f"eps{k}"].to_numpy(); c = g[f"cc{k}"].to_numpy(); ok = np.isfinite(e)
                    if ok.sum() < 3:
                        continue
                    aft = (g.bin_mid > g.front).to_numpy()
                    units.append(dict(key=key, e=e[ok], w=c[ok] ** 2, aft=aft[ok]))
                for stat in ("front", "third", "front1_posthoc"):
                    def calc(perm):
                        num = den = 0.0; per_sta = {}
                        for u, pp in zip(units, perm):
                            e = u["e"][pp] if pp is not None else u["e"]; w = u["w"][pp] if pp is not None else u["w"]
                            aft = u["aft"]
                            if stat.startswith("front"):
                                nmin = 1 if stat == "front1_posthoc" else 2
                                if aft.sum() < nmin or (~aft).sum() < nmin:
                                    continue
                                a, b = w[aft].sum(), w[~aft].sum(); d = _wm(e[aft], w[aft]) - _wm(e[~aft], w[~aft])
                            else:
                                m = len(e) // 3
                                a, b = w[-m:].sum(), w[:m].sum(); d = _wm(e[-m:], w[-m:]) - _wm(e[:m], w[:m])
                            ww = 2 * a * b / (a + b)
                            num += ww * d; den += ww
                            s = T.site(u["key"][1]); ps = per_sta.setdefault(u["key"][1], [0.0, 0.0, s]); ps[0] += ww * d; ps[1] += ww
                        return (num / den if den > 0 else np.nan), per_sta, den
                    obs, per_sta, den = calc([None] * len(units))
                    if not np.isfinite(obs):
                        rows.append(dict(res=rn, fams="+".join(map(str, famset)), win=k, stat=stat, n_units=0)); continue
                    null = np.array([calc([rng.permutation(len(u["e"])) for u in units])[0] for _ in range(NSHUF)])
                    nn = null[np.isfinite(null)]
                    sta = {s: v[0] / v[1] for s, v in per_sta.items() if v[1] > 0}
                    sites_same = sorted({per_sta[s][2] for s, v in sta.items() if np.sign(v) == np.sign(obs)})
                    nm_ = 1 if stat == "front1_posthoc" else 2
                    nunits = sum(1 for u in units if (stat == "third") or (u["aft"].sum() >= nm_ and (~u["aft"]).sum() >= nm_))
                    rows.append(dict(res=rn, fams="+".join(map(str, famset)), win=k, stat=stat, n_units=nunits,
                                     delta_eps_pct=100 * obs, delta_dvv_pct=-100 * obs, shuffle_sd_pct=100 * nn.std(),
                                     p_two=(1 + np.sum(np.abs(nn - nn.mean()) >= abs(obs - nn.mean()))) / (1 + len(nn)),
                                     n_stations=len(sta), stations_same_sign=sum(np.sign(v) == np.sign(obs) for v in sta.values()),
                                     sites_same_sign=len(sites_same), sites_same_list=";".join(sites_same),
                                     bound99_dvv_pct=100 * (abs(obs) + 2.58 * nn.std()),
                                     per_station=json.dumps({s: round(-100 * v, 4) for s, v in sta.items()})))
    D = pd.DataFrame(rows)
    # pass rule (B4): p < 0.01 (2-4 s), same sign 3-5 s, >= 2 sites same sign
    out = []
    for (rn, fs, st), g in D.groupby(["res", "fams", "stat"]):
        a = g[g.win == "24"]; b = g[g.win == "35"]
        if len(a) and "p_two" in a and np.isfinite(a.p_two.iloc[0]):
            ok = (a.p_two.iloc[0] < 0.01) and len(b) and np.sign(a.delta_eps_pct.iloc[0]) == np.sign(b.delta_eps_pct.iloc[0]) \
                and a.sites_same_sign.iloc[0] >= 2
            out.append(dict(res=rn, fams=fs, stat=st, pass_rule=bool(ok)))
    return D.merge(pd.DataFrame(out), on=["res", "fams", "stat"], how="left") if out else D


def main():
    T.ensure(); log = S.Log(f"{T.LOGD}/subday.log")
    if "--tests-only" in sys.argv:
        B = pd.read_csv(f"{T.TB}/subday_bins.csv", parse_dates=["bin_start", "bin_mid", "t_mean", "front"])
        D = tests(B, np.random.default_rng(SEED)); D.to_csv(f"{T.TB}/subday_tests.csv", index=False); log("tests rewritten"); return
    R = pd.read_parquet(f"{T.TB}/rows_tide.parquet")
    R["tt"] = pd.to_datetime(R.t, unit="s")
    EPI = pd.read_csv(f"{S.ROOT}/svi_science/results/r23/tables/episodes_2003_2026.csv")
    FR = fronts(R, EPI)
    FR.to_csv(f"{T.TB}/subday_fronts.csv", index=False)
    log(f"fronts: {len(FR)}")
    B = measure_bins(R, FR)
    B.to_csv(f"{T.TB}/subday_bins.csv", index=False)
    log(f"bins: {len(B)} ({(B.res == '6h').sum()} 6-h, {(B.res == '12h').sum()} 12-h)")
    D = tests(B, np.random.default_rng(SEED))
    D.to_csv(f"{T.TB}/subday_tests.csv", index=False)
    log("tests written")
    open(f"{T.OUT}/SUBDAY_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
