"""Q8 RTR analysis: thresholds, Royer match (step 1), controls, questions a-e (PLAN.md §3-4).

Reads tables/streaks_all.parquet, members_data.parquet. Writes tables/*.csv and results.json.
"""
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L

ROOT = L.ROOT
OUT = f"{L.OUT}/tables"
RES = {}
RNG = np.random.default_rng(20261008)
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
SGRID = np.round(np.arange(2.0, 60.01, 0.25), 2)


def curve(S, cat, var, mode_data, mode_c1, mode_c2, cls="RTR"):
    d = S[(S.cat == cat) & (S["var"] == var) & (S.cls == cls)]
    nd = np.array([(d[d["mode"] == mode_data].S >= s).sum() for s in SGRID], float)
    reps1 = d[d["mode"] == mode_c1].rep.unique(); reps2 = d[d["mode"] == mode_c2].rep.unique()
    run = pd.read_csv(f"{OUT}/runs{SUF}.csv")
    n1 = run[(run.cat == cat) & (run["var"] == var) & (run["mode"] == mode_c1)].rep.nunique()
    n2 = run[(run.cat == cat) & (run["var"] == var) & (run["mode"] == mode_c2)].rep.nunique()
    nreps_d = run[(run.cat == cat) & (run["var"] == var) & (run["mode"] == mode_data)].rep.nunique()
    nd = nd / max(nreps_d, 1)
    c1 = np.array([(d[d["mode"] == mode_c1].S >= s).sum() for s in SGRID], float) / max(n1, 1)
    c2 = np.array([(d[d["mode"] == mode_c2].S >= s).sum() for s in SGRID], float) / max(n2, 1)
    return pd.DataFrame(dict(S=SGRID, data=nd, c1=c1, c2=c2, cat=cat, var=var, mode=mode_data))


def sstar(cv, frac=0.10):
    worst = np.maximum(cv.c1.values, cv.c2.values)
    ok = worst <= frac * cv.data.values
    # smallest s such that ok for all larger s
    bad = np.nonzero(~ok)[0]
    i = 0 if len(bad) == 0 else bad[-1] + 1
    return float(cv.S.values[i]) if i < len(cv) else np.inf


def match(our_t, ref_t, tol_h=6.0):
    """greedy one-to-one by smallest |dt|; returns list of (i_our, j_ref, dt_h)."""
    pairs = []
    for i, a in enumerate(our_t):
        for j, b in enumerate(ref_t):
            d = (a - b) / np.timedelta64(1, "h")
            if abs(d) <= tol_h:
                pairs.append((abs(d), i, j, d))
    pairs.sort()
    ui, uj, out = set(), set(), []
    for _, i, j, d in pairs:
        if i in ui or j in uj:
            continue
        ui.add(i); uj.add(j); out.append((i, j, d))
    return out


def main():
    S = pd.read_parquet(f"{OUT}/streaks_all{SUF}.parquet")
    M = pd.read_parquet(f"{OUT}/members_data{SUF}.parquet")
    loaders = {("bostock", "all"): qd.bostock(), ("svi", "head"): qd.ours_svi("head"), ("svi", "p050"): qd.ours_svi("p050"),
               ("svi", "raw075"): qd.ours_svi("raw075"), ("svi", "scr075"): qd.ours_svi("scr075"),
               ("nw", "head"): qd.ours_nw("head"), ("nw", "p050"): qd.ours_nw("p050"), ("nw", "raw075"): qd.ours_nw("raw075")}
    # absolute times
    def absify(D):
        D = D.copy()
        w0 = pd.Series({(c, v, e): r.w0 for (c, v), (_, E) in loaders.items() for r in E.itertuples() for e in [r.ep]})
        on = pd.Series({(c, v, e): r.onset for (c, v), (_, E) in loaders.items() for r in E.itertuples() for e in [r.ep]})
        key = list(zip(D.cat, D["var"], D.ep))
        D["w0"] = w0.reindex(key).values; D["onset"] = on.reindex(key).values
        D["start"] = D.w0 + pd.to_timedelta(D.t_first, unit="h"); D["end"] = D.w0 + pd.to_timedelta(D.t_last, unit="h")
        D["dur_h"] = D.t_last - D.t_first
        D["days_after_onset"] = (D.start - D.onset).dt.total_seconds() / 86400
        return D
    S = absify(S)

    # ---------------------------------------------------------------- thresholds
    curves = [curve(S, "bostock", "all", "data", "c1", "c2"), curve(S, "bostock", "all", "thin", "thin_c1", "thin_c2")]
    for c, v in [("svi", "head"), ("svi", "p050"), ("svi", "raw075"), ("nw", "head"), ("nw", "p050"), ("nw", "raw075")]:
        if not ((S.cat == c) & (S["var"] == v)).any():
            continue
        curves.append(curve(S, c, v, "data", "c1", "c2"))
    CV = pd.concat(curves); CV.to_csv(f"{OUT}/threshold_curves{SUF}.csv", index=False)
    SST = float(sstar(curves[0]))
    RES["S_star_rule_value"] = SST
    if not np.isfinite(SST):          # D3: 10 % rule cannot be met (position-permuted control keeps surges) -> all streaks
        SST = 2.0
        RES["S_star_rule_failed_use_all"] = True
    RES["S_star_frozen"] = SST
    own = {}
    for cv in curves:
        k = f"{cv.cat.iloc[0]}|{cv['var'].iloc[0]}|{cv['mode'].iloc[0]}"
        own[k] = sstar(cv)
    RES["S_star_own"] = own
    print("S* frozen", SST, own)

    # control summary at S* and own
    rows = []
    for cv in curves:
        k = f"{cv.cat.iloc[0]}|{cv['var'].iloc[0]}|{cv['mode'].iloc[0]}"
        for lab, s in [("frozen", SST), ("own", own[k])]:
            if not np.isfinite(s):
                rows.append(dict(set=k, threshold=lab, S=s, data=0, c1=np.nan, c2=np.nan)); continue
            i = np.searchsorted(cv.S.values, s)
            rows.append(dict(set=k, threshold=lab, S=s, data=cv.data.values[i], c1=cv.c1.values[i], c2=cv.c2.values[i],
                             false_frac=max(cv.c1.values[i], cv.c2.values[i]) / max(cv.data.values[i], 1e-9)))
    CS = pd.DataFrame(rows); CS.to_csv(f"{OUT}/controls_summary{SUF}.csv", index=False)
    print(CS.to_string())

    # ---------------------------------------------------------------- step 1: Royer match
    Ro = pd.read_csv(f"{L.OUT}/tables/royer2015_rtr_list.csv", parse_dates=["start_utc", "end_utc"])
    ep_of_year = {2003: "2003-02", 2004: "2004-07", 2005: "2005-09", 2008: "2008-05", 2010: "2010-08", 2011: "2011-08", 2012: "2012-08"}
    Ro["ep"] = Ro.year.map(ep_of_year)
    B = S[(S.cat == "bostock") & (S["mode"] == "data")].copy()
    RT = B[(B.cls == "RTR") & (B.S >= SST)].copy()
    RT7 = RT[RT.ep.isin(Ro.ep.unique())]
    step1 = {}
    step1["n_rtr_7eps"] = int(len(RT7)); step1["n_rtr_all10"] = int(len(RT))
    step1["per_episode"] = RT.groupby("ep").size().to_dict()
    step1["royer_per_episode"] = Ro.groupby("ep").size().to_dict()
    step1["median_days_after_onset"] = float(RT7.days_after_onset.median())
    step1["royer_median_days_after_onset"] = float(np.median([(r.start_utc - loaders[("bostock", "all")][1].set_index("ep").onset[r.ep]).total_seconds() / 86400 for r in Ro.itertuples()]))

    def recall(our, ref, shift_h=0.0, tol=6.0):
        n = 0; npair = 0; dts = []
        for ep in ref.ep.unique():
            a = our[our.ep == ep].start.values + np.timedelta64(int(shift_h * 3600), "s")
            b = ref[ref.ep == ep].start_utc.values
            m = match(a, b, tol); n += len(m); dts += [d for _, _, d in m]
        return n, dts
    n, dts = recall(RT7, Ro)
    step1["matched"] = n; step1["recall"] = n / len(Ro); step1["precision"] = n / max(len(RT7), 1)
    step1["median_abs_dt_h"] = float(np.median(np.abs(dts))) if dts else np.nan
    ch = [recall(RT7, Ro, s)[0] / len(Ro) for s in (-72, -48, -24, 24, 48, 72)]
    step1["chance_recall_dayshift_mean"] = float(np.mean(ch)); step1["chance_recall_dayshift"] = ch
    ch2 = [recall(RT7, Ro, s)[0] / len(Ro) for s in (-60, -36, -12, 12, 36, 60)]
    step1["chance_recall_halfdayshift"] = ch2
    step1["tz_diag_recall"] = {f"{h:+d}h": recall(RT7, Ro, -h)[0] / len(Ro) for h in (-8, -7, 0, 7, 8)}
    # recall using every streak class (diagnostic)
    ALL7 = B[(B.S >= SST) & B.ep.isin(Ro.ep.unique())]
    step1["recall_any_class"] = recall(ALL7, Ro)[0] / len(Ro)
    step1["n_any_class_7eps"] = int(len(ALL7))
    for c in ["RTR", "forward", "other"]:
        step1[f"recall_{c}_only_S2"] = recall(B[(B.cls == c) & B.ep.isin(Ro.ep.unique())], Ro)[0] / len(Ro)
    # recall vs threshold
    rv = []
    for s in [2, 4, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50]:
        q = B[(B.cls == "RTR") & (B.S >= s) & B.ep.isin(Ro.ep.unique())]
        rv.append(dict(S=s, n=len(q), recall=recall(q, Ro)[0] / len(Ro),
                       chance=np.mean([recall(q, Ro, sh)[0] for sh in (-48, -24, 24, 48)]) / len(Ro)))
    pd.DataFrame(rv).to_csv(f"{OUT}/royer_recall_vs_S{SUF}.csv", index=False)
    step1["PASS_count"] = bool(43 <= len(RT7) <= 96)
    step1["PASS_timing"] = bool(step1["median_days_after_onset"] >= 2)
    step1["PASS_match"] = bool(step1["recall"] >= 0.70)
    step1["PASS"] = bool(step1["PASS_count"] and step1["PASS_timing"] and step1["PASS_match"])
    RES["step1"] = step1
    print(json.dumps(step1, indent=1, default=str))
    # matched table
    rows = []
    for ep in Ro.ep.unique():
        a = RT7[RT7.ep == ep].reset_index(drop=True); b = Ro[Ro.ep == ep].reset_index(drop=True)
        m = match(a.start.values, b.start_utc.values, 6.0)
        mj = {j: (i, d) for i, j, d in m}
        for j, r in b.iterrows():
            if j in mj:
                i, d = mj[j]; q = a.iloc[i]
                rows.append(dict(rtr=r.rtr, ep=ep, royer_start=r.start_utc, our_start=q.start, dt_h=d, royer_v=r.v_kmh, our_v=abs(q.v_ev),
                                 royer_span=r.dist_strike_km, our_span=q.span_ev, royer_dur=r.dur_h, our_dur=q.dur_h, our_S=q.S, our_uid=q.uid))
            else:
                rows.append(dict(rtr=r.rtr, ep=ep, royer_start=r.start_utc))
    pd.DataFrame(rows).to_csv(f"{OUT}/royer_match{SUF}.csv", index=False)

    RES["SST"] = SST
    json.dump(RES, open(f"{L.OUT}/results_step1{SUF}.json", "w"), indent=1, default=str)
    S.to_parquet(f"{OUT}/streaks_all_abs{SUF}.parquet", index=False)


if __name__ == "__main__":
    main()
