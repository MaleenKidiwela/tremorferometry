#!/usr/bin/env python
"""r23_front (2026-10-11) -- R23 study 1: front arrest, SVI LFE check with 2003-2013 added (finding 8; Q6 rules).
Reuses q6_lfe.episode_tables / run_check / perm unchanged (outputs redirected). Rules: svi_science/results/r23/PLAN.md §2.
    python svi_science/scripts/r23_front.py
"""
from __future__ import annotations
import os, sys, json, importlib.util
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C

OUT = f"{C.TB}/study1"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/study1.log")
BIN = 5.0


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


ql = load("q6lfe", f"{C.ROOT}/svi_science/scripts/q6_lfe.py")
ql.OUT = OUT                      # walk_rows_*.csv written here, never into q6_front_arrest
b3, qt = ql.b3, ql.qt
SEGFAM = {"CN4f": "ERA", "ERA": "ERA", "B2f": "B", "B3f": "B", "S15f": "S15f", "FWD": "FWD"}


def along(lat, lon):
    fr = json.load(open(f"{C.S0}/frame_groups.json"))
    x = (np.asarray(lon) - fr["origin_lon"]) * 111.32 * np.cos(np.deg2rad(fr["origin_lat"]))
    y = (np.asarray(lat) - fr["origin_lat"]) * 110.574
    s = np.deg2rad(fr["strike_az"]); return x * np.sin(s) + y * np.cos(s)


def run_variant(name, V, E, targets, s3start, period, rng, trem):
    a, b = C.PERIODS[period]
    Ep = E[(E.start <= b) & (E.end >= a)].copy()
    Ew = b3.padded_windows(Ep[["episode", "start", "end"]].copy())
    T = V[(V.date >= a) & (V.date <= b)][["time_utc", "fam"]].rename(columns={"time_utc": "time"}).copy()
    T["ep"] = b3.assign(T.time, Ew); T = T[T.ep.notna()]
    F = C.families()
    T["along_strike_km"] = F.along_strike_km.reindex(T.fam).to_numpy(); T["bin"] = np.floor(T.along_strike_km / BIN).astype(int)
    fam_bins = set(np.floor(F.along_strike_km / BIN).astype(int))

    def st(ep, g):
        if ep in s3start.index:
            return float(s3start[ep])
        g = g.sort_values("time"); k = max(1, int(np.ceil(0.1 * len(g))))
        return float(np.median(g.along_strike_km.values[:k]))
    tg = set(targets) & set(T.ep.unique())
    Cc, be, W, S = ql.episode_tables(T[["time", "ep", "bin", "along_strike_km"]].copy(), Ew, fam_bins, st, tg)
    segmap = dict(zip(E.episode, E.seg))
    groups = pd.Series([SEGFAM[segmap[e]] for e in Cc.ep], index=Cc.index)
    res = ql.run_check(name, Cc, be, W, S, rng, groups=groups)
    res.update(n_episodes=int(len(Cc)), n_targets=int(Cc.big.sum()), n_walk=int(len(W)), n_stops=int((W.y == 0).sum()) if len(W) else 0,
               n_enter=int(len(S)), targets=sorted(tg))
    if len(W):
        ww = pd.read_csv(f"{OUT}/walk_rows_{name}.csv")
        ww["Tcls"] = pd.cut(ww.T_any, [0, 30, 90, 1e6], labels=["<30 d", "30-90 d", ">90 d"])
        tab = ww.groupby("Tcls", observed=False).y.agg(["size", lambda x: int((x == 0).sum())])
        tab.columns = ["steps", "stops"]; tab["stop_share"] = tab.stops / tab.steps.clip(lower=1)
        res["stop_share_by_Tany"] = tab.reset_index().to_dict("records")
        res["stops_list"] = ww[ww.y == 0][["ep", "bin", "dirn", "T_any", "T_big"]].to_dict("records")
        # tremor-history column with the same order shuffle
        res["trem"] = trem_column(name, Cc, W, Ew, trem, rng, groups)
    log(f"[{name}] " + json.dumps({k: v for k, v in res.items() if k not in ("stops_list",)}, default=str))
    return res


def trem_column(name, Cc, W, Ew, trem, rng, groups, ndraw=2000):
    """T_trem = days from the episode's first event (slot) back to the last tremor-active bin-day before its padded window."""
    t0 = pd.Timestamp("2000-01-01")
    w0 = Ew.set_index("episode").w0
    pre = {r.comp: r.t_first - (pd.Timestamp(w0[r.ep]) - t0) / pd.Timedelta(days=1) for r in Cc.itertuples()}
    keys = np.sort(trem.bin.to_numpy() * 1e5 + trem.td.to_numpy())

    def T_of(slot):
        qt_ = np.array([slot[c] - pre[c] for c in W.comp])
        q = W.bin.to_numpy() * 1e5 + qt_
        i = np.searchsorted(keys, q, side="left") - 1
        out = np.full(len(q), np.nan); ok = i >= 0
        v = keys[np.clip(i, 0, None)]; ok &= np.floor(v / 1e5 + 1e-9) == W.bin.to_numpy()
        out[ok] = (slot[W.comp.to_numpy()] - (v - W.bin.to_numpy() * 1e5))[ok]
        return out
    real = np.zeros(Cc.comp.max() + 1); real[Cc.comp.values] = Cc.t_first.values
    T = T_of(real)
    beta = qt.beta_T(T, W.y.values, W.dist.values, W.lsize.values, W.prop.values, np.ones(len(W), bool))
    nl = []
    for _ in range(ndraw):
        nl.append(qt.beta_T(T_of(ql.perm(Cc, rng, groups)), W.y.values, W.dist.values, W.lsize.values, W.prop.values, np.ones(len(W), bool)))
    nl = np.array(nl)
    m = np.isfinite(T)
    cls = pd.cut(pd.Series(T[m]), [0, 30, 90, 1e6], labels=["<30 d", "30-90 d", ">90 d"])
    tab = pd.DataFrame(dict(c=cls, y=W.y.values[m])).groupby("c", observed=False).y.agg(["size", lambda x: int((x == 0).sum())])
    tab.columns = ["steps", "stops"]
    return dict(beta=beta, p_pos=qt.pval(beta, nl, "+"), p_neg=qt.pval(beta, nl, "-"), null_mean=float(np.nanmean(nl)),
                n_rows=int(m.sum()), n_stops=int((W.y.values[m] == 0).sum()), by_T=tab.reset_index().to_dict("records"))


def main():
    rng = np.random.default_rng(0)
    F = C.families(); E = C.episodes()
    V = pd.read_parquet(f"{C.TB}/events_2003_2026.parquet", columns=["time_utc", "date", "fam", "series"])
    V = V[V.series == "data"]
    k3 = pd.read_csv(f"{C.ROOT}/svi_science/results/s3/episodes_kinematics.csv")
    s3start = k3[k3.region == "SVI"].set_index("episode").lfe_start10
    s3big = json.load(open(f"{C.ROOT}/svi_science/results/s3/build_info.json"))["big"]["SVI"]
    # tremor bin-days (>= 3 windows/epicentres within 30 km of a family, same UTC day, same 5-km bin)
    TR = pd.read_parquet(f"{C.TB}/tremor_points.parquet")
    TR = TR[(TR.dmin_fam_km <= 30) & (((TR.src == "env") & (TR.date < C.PNSN_T0)) | ((TR.src == "pnsn") & (TR.date >= C.PNSN_T0)))].copy()
    TR["bin"] = np.floor(along(TR.lat, TR.lon) / BIN).astype(int)
    tb = TR.groupby(["bin", "date"]).size().reset_index(name="n"); tb = tb[tb.n >= 3]
    tb["td"] = (pd.to_datetime(tb.date) - pd.Timestamp("2000-01-01")) / pd.Timedelta(days=1) + 0.5
    log(f"tremor-active bin-days: {len(tb)}")
    big = E.episode[E.big].tolist(); majors = E.episode[E.major].tolist()
    R = {}
    for period in ("2003-26", "2003-13", "2014-26"):
        R[f"big|{period}"] = run_variant(f"big_{period}", V, E, big, s3start, period, rng, tb)
        R[f"majors|{period}"] = run_variant(f"majors_{period}", V, E, majors, s3start, period, rng, tb)
    R["s3big|2014-26 (reproduction)"] = run_variant("s3big_2014-26", V, E, s3big, s3start, "2014-26", rng, tb)
    # verdict (joint, walk, T_any)
    j = R["big|2003-26"]; w = j.get("walk_Tany", {})
    ns = j["n_stops"]
    if ns < 20:
        verdict = "undecided (too few stops)"
    elif w["beta"] > 0 and w["p_pos"] < 0.01:
        verdict = "supported"
    elif w["beta"] > 0 and w["p_pos"] < 0.05:
        verdict = "suggestive"
    elif w["beta"] < 0 and w["p_neg"] < 0.01:
        verdict = "contradicted"
    else:
        verdict = "not supported"
    R["verdict"] = dict(n_stops=ns, beta=w.get("beta"), p_pos=w.get("p_pos"), p_neg=w.get("p_neg"), verdict=verdict)
    log(f"VERDICT: {R['verdict']}")
    C.dump(R, f"{OUT}/study1_summary.json")
    rows = []
    for k, r in R.items():
        if k == "verdict":
            continue
        for m in ("walk_Tany", "walk_Tbig", "enter_Tany", "enter_Tbig"):
            if m in r:
                rows.append(dict(variant=k, model=m, beta=r[m]["beta"], p_pos=r[m]["p_pos"], p_neg=r[m]["p_neg"], ctrl_mean=r[m]["null_mean"],
                                 rows=r[m]["n_rows"], stops_or_zeros=r[m]["n_y0"], n_targets=r["n_targets"], n_episodes=r["n_episodes"]))
        if "trem" in r:
            t = r["trem"]; rows.append(dict(variant=k, model="walk_Ttrem", beta=t["beta"], p_pos=t["p_pos"], p_neg=t["p_neg"], ctrl_mean=t["null_mean"],
                                            rows=t["n_rows"], stops_or_zeros=t["n_stops"], n_targets=r["n_targets"], n_episodes=r["n_episodes"]))
    pd.DataFrame(rows).to_csv(f"{OUT}/front_arrest_table.csv", index=False)
    log("table:\n" + pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
