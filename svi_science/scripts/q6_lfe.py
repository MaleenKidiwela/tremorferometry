"""Q6 LFE family-resolved checks.

Q6.1 (front arrest, walk + enter models, episode-order shuffle) on:
  SVI ours 2014+ (79 episodes MAJOR+minor history; S3 big 10 targets), Bostock 2003-2013 (10 ETS), joint SVI 2003-2026
  (Bostock 10 + our big 10, big-only history, shuffle within era), NW E3 column (S3 big 14; events file covers major
  windows only -> history = majors).
Q6.2 (Strait starts) on joint SVI+NW big LFE episodes on the PNSN margin axis.
Rules: svi_science/results/q6_front_arrest/PLAN.md.
"""
import importlib.util
import json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q6_front_arrest"
ER = f"{ROOT}/lfe_extension_wa/results/episode_rates"
S3 = f"{ROOT}/svi_science/results/s3"
BIN = 5.0
NDRAW = 2000


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


b3 = load("s3b", f"{ROOT}/svi_science/scripts/s3_build.py")
qb = load("q6b", f"{ROOT}/svi_science/scripts/q6_build.py")
qt = load("q6t", f"{ROOT}/svi_science/scripts/q6_tests.py")


def episode_tables(ev, E, fam_bins, start_of, big_set):
    """ev: time (UTC ns), ep, bin. E: episodes with w0, w1, episode. Returns comp table, pairs, rows (walk, enter)."""
    ev = ev.sort_values("time")
    t0 = pd.Timestamp("2000-01-01")
    ev["td"] = (ev.time - t0) / pd.Timedelta(days=1)
    comps = []
    be = ev.groupby(["ep", "bin"]).agg(n=("td", "size"), t_med=("td", "median")).reset_index()
    be["active"] = be.n >= 3
    fb = np.array(sorted(fam_bins)); bmin, bmax = fb.min(), fb.max()
    W, S = [], []
    for i, (ep, g) in enumerate(ev.groupby("ep")):
        comps.append(dict(ep=ep, comp=i, t_first=float(g.td.iloc[0]), t_last=float(g.td.iloc[-1]),
                          dur=float(g.td.iloc[-1] - g.td.iloc[0]), n=len(g), big=ep in big_set))
    C = pd.DataFrame(comps)
    be = be.merge(C[["ep", "comp"]], on="ep")
    act = {c: set(g.bin.values) for c, g in be[be.active].groupby("comp")}
    for _, r in C[C.big].iterrows():
        A = act.get(r.comp, set())
        if not A:
            continue
        s = start_of(r.ep, ev[ev.ep == r.ep])
        sb = int(np.floor(s / BIN))
        lo, hi, sb2 = qb.front_run(np.array(sorted(A)), sb, bmin - 1, bmax + 1)
        base = dict(comp=int(r.comp), ep=r.ep, lsize=np.log10(r.n), start_s=s)
        for step, end in ((+1, hi), (-1, lo)):
            bb = sb2 + step
            while (step > 0 and bb <= end) or (step < 0 and bb >= end):
                if bb in fam_bins:
                    W.append(dict(base, bin=bb, y=1, dirn=step))
                bb += step
            if bb in fam_bins:
                W.append(dict(base, bin=bb, y=0, dirn=step))       # stop inside the family-covered span
        for bb in fam_bins:
            S.append(dict(base, bin=bb, y=int(bb in A), dirn=0))
    W = pd.DataFrame(W); S = pd.DataFrame(S)
    cnt = ev.groupby("bin").size()
    for R in (W, S):
        R["dist"] = np.abs((R.bin + 0.5) * BIN - R.start_s)
        own = be.set_index(["comp", "bin"]).n
        ownn = np.array([own.get((c, b_), 0) for c, b_ in zip(R.comp, R.bin)])
        R["prop"] = np.log10(cnt.reindex(R.bin).fillna(0).values - ownn + 1.0)
    return C, be, W, S


class H2(qt.History):
    def __init__(self, be, C):
        a = be[be.active].merge(C[["comp", "t_first", "dur", "big"]], on="comp")
        self.pc = a.comp.values; self.pb = a.bin.values.astype(float) + 1000.0   # keep keys positive
        self.po = (a.t_med - a.t_first).values; self.pdur = a.dur.values; self.pbig = a.big.values


def perm(C, rng, groups=None):
    s = np.zeros(C.comp.max() + 1)
    g = C.big.astype(str) if groups is None else C.big.astype(str) + "_" + groups
    for k in g.unique():
        m = (g == k).values
        s[C.comp.values[m]] = C.t_first.values[m][rng.permutation(m.sum())]
    return s


def run_check(name, C, be, W, S, rng, groups=None):
    H = H2(be, C)
    real = np.zeros(C.comp.max() + 1); real[C.comp.values] = C.t_first.values
    out = {}

    def stats(slot):
        r = {}
        for nm, R in (("walk", W), ("enter", S)):
            if len(R) == 0:
                continue
            Ta, Tb, T7 = H.T_all(slot, R.bin.values.astype(float) + 1000.0, R.comp.values)
            for k, T in (("Tany", Ta), ("Tbig", Tb)):
                r[f"{nm}_{k}"] = qt.beta_T(T, R.y.values, R.dist.values, R.lsize.values, R.prop.values, np.ones(len(R), bool), f"{name}_{nm}_{k}")
        return r
    rr = stats(real)
    nn = {k: [] for k in rr}
    for _ in range(NDRAW):
        x = stats(perm(C, rng, groups))
        for k in x:
            nn[k].append(x[k])
    for k in rr:
        nl = np.array(nn[k])
        out[k] = dict(beta=rr[k], p_pos=qt.pval(rr[k], nl, "+"), p_neg=qt.pval(rr[k], nl, "-"),
                      null_mean=float(np.nanmean(nl)) if np.isfinite(nl).any() else np.nan,
                      n_rows=qt.NROWS.get(f"{name}_{k}", (0, 0))[0], n_y0=qt.NROWS.get(f"{name}_{k}", (0, 0))[1])
    # descriptive medians
    Ta, Tb, T7 = H.T_all(real, W.bin.values.astype(float) + 1000.0, W.comp.values) if len(W) else (np.array([]),) * 3
    if len(W):
        out["walk_median_Tany_continue"] = float(np.nanmedian(Ta[W.y.values == 1]))
        out["walk_median_Tany_stop"] = float(np.nanmedian(Ta[W.y.values == 0])) if (W.y == 0).any() else np.nan
        W = W.assign(T_any=Ta, T_big=Tb)
        W.to_csv(f"{OUT}/walk_rows_{name}.csv", index=False)
    return out


def main():
    rng = np.random.default_rng(0)
    res = {}
    k3 = pd.read_csv(f"{S3}/episodes_kinematics.csv")
    big = json.load(open(f"{S3}/build_info.json"))["big"]
    # ---------------- SVI ours 2014+
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    F["bin"] = np.floor(F.along_strike_km / BIN).astype(int)
    Esvi = b3.padded_windows(pd.read_csv(f"{ER}/episodes_svi_final.csv"))
    T = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet",
                        columns=["run", "series", "fam", "P_same", "keep_xfam8", "time_utc"], filters=[("P_same", ">=", 0.75)])
    T = T[(T.run.astype(str) == "prod") & (T.series.astype(str) == "data") & T.keep_xfam8].rename(columns={"time_utc": "time"})
    T["ep"] = b3.assign(T.time, Esvi); T = T[T.ep.notna()]
    T = T.merge(F[["fid", "bin", "along_strike_km"]].rename(columns={"fid": "fam"}), on="fam")
    s3start = k3[k3.region == "SVI"].set_index("episode").lfe_start10

    def st_svi(ep, g):
        return float(s3start[ep]) if ep in s3start.index else float(np.median(g.sort_values("time").along_strike_km.values[:max(1, int(np.ceil(0.1 * len(g))))]))
    fam_bins = set(F.bin.unique())
    C, be, W, S = episode_tables(T[["time", "ep", "bin", "along_strike_km"]].copy(), Esvi, fam_bins, st_svi, set(big["SVI"]))
    res["svi_ours"] = run_check("svi_ours", C, be, W, S, rng)
    res["svi_ours"]["n_episodes"] = int(len(C)); res["svi_ours"]["n_big"] = int(C.big.sum())
    print(json.dumps(res["svi_ours"], indent=1, default=float), flush=True)
    del T
    # ---------------- Bostock 2003-2013
    B = pd.read_parquet(f"{ROOT}/svi_science/results/p3_bostock/bostock_events_mapped.parquet", columns=["fid", "time_utc", "along_km"])
    B = B.rename(columns={"time_utc": "time"})
    Eb = pd.read_csv(f"{ROOT}/svi_science/results/p3_bostock/episodes_bostock_final.csv")
    Eb["w0"] = pd.to_datetime(Eb.start); Eb["w1"] = pd.to_datetime(Eb.end) + pd.Timedelta(days=1)
    B["ep"] = b3.assign(B.time, Eb); B = B[B.ep.notna()]
    B["bin"] = np.floor(B.along_km / BIN).astype(int)

    def st_med(ep, g):
        g = g.sort_values("time"); col = "along_km" if "along_km" in g else "along_strike_km"
        return float(np.median(g[col].values[:max(1, int(np.ceil(0.1 * len(g))))]))
    bb = set(np.floor(pd.read_csv(f"{S3}/families_svi_frame.csv").along_strike_km / BIN).astype(int))
    Cb, beb, Wb, Sb = episode_tables(B[["time", "ep", "bin", "along_km"]].copy(), Eb, bb, st_med, set(Eb.episode))
    res["bostock"] = run_check("bostock", Cb, beb, Wb, Sb, rng)
    print(json.dumps(res["bostock"], indent=1, default=float), flush=True)
    # ---------------- joint SVI 2003-2026 (big only)
    T = pd.read_parquet(f"{ROOT}/svi_science/results/s0/analysis_table.parquet",
                        columns=["run", "series", "fam", "P_same", "keep_xfam8", "time_utc"], filters=[("P_same", ">=", 0.75)])
    T = T[(T.run.astype(str) == "prod") & (T.series.astype(str) == "data") & T.keep_xfam8].rename(columns={"time_utc": "time"})
    Eo = Esvi[Esvi.episode.isin(big["SVI"])]
    T["ep"] = b3.assign(T.time, Eo); T = T[T.ep.notna()]
    T = T.merge(F[["fid", "bin", "along_strike_km"]].rename(columns={"fid": "fam"}), on="fam").rename(columns={"along_strike_km": "along_km"})
    T["ep"] = "o" + T.ep
    J = pd.concat([B[["time", "ep", "bin", "along_km"]], T[["time", "ep", "bin", "along_km"]]], ignore_index=True)

    def st_joint(ep, g):
        return float(s3start[ep[1:]]) if ep.startswith("o") else st_med(ep, g)
    Cj, bej, Wj, Sj = episode_tables(J, None, fam_bins | bb, st_joint, set(J.ep.unique()))
    era = pd.Series(np.where(Cj.ep.str.startswith("o"), "ours", "bostock"), index=Cj.index)
    res["svi_joint"] = run_check("svi_joint", Cj, bej, Wj, Sj, rng, groups=era)
    print(json.dumps(res["svi_joint"], indent=1, default=float), flush=True)
    del T, J, B
    # ---------------- NW E3 column (headline; major-window history only)
    N = pd.read_parquet(f"{S3}/events_nw.parquet")
    N = N[(N.series == "data") & N.selected_v2 & N.counted_once & N.quality.isin(["A", "B"])].copy()
    N["bin"] = np.floor(N.as_km / BIN).astype(int)
    s3n = k3[k3.region == "NW"].set_index("episode").lfe_start10

    def st_nw(ep, g):
        return float(s3n[ep]) if ep in s3n.index else st_med(ep, g.rename(columns={"as_km": "along_km"}))
    srcb = set(np.floor(pd.read_csv(f"{S3}/sources_nw.csv").query("quality in ['A','B']").as_km / BIN).astype(int))
    Cn, ben, Wn, Sn = episode_tables(N[["time", "ep", "bin", "as_km"]].copy(), None, srcb, st_nw, set(big["NW"]))
    res["nw_e3"] = run_check("nw_e3", Cn, ben, Wn, Sn, rng)
    print(json.dumps(res["nw_e3"], indent=1, default=float), flush=True)

    # ---------------- Q6.2 LFE: joint SVI+NW big episodes on the margin axis
    nodes = pd.read_csv(f"{OUT}/axis_nodes.csv")
    nodes["s_raw"] = nodes.s_raw
    ss0 = float((nodes.s_raw - nodes.s_km).iloc[0])
    Sv = pd.read_parquet(f"{S3}/events_svi.parquet")
    Sv = Sv[(Sv.series == "data") & (Sv.P_same >= 0.75) & Sv.keep_xfam8 & Sv.counted_once & Sv.ep.isin(big["SVI"])].copy()
    Nv = pd.read_parquet(f"{S3}/events_nw.parquet")
    Nv = Nv[(Nv.series == "data") & Nv.selected_v2 & Nv.counted_once & Nv.ep.isin(big["NW"])].copy()
    for D in (Sv, Nv):
        s, _ = qb.project(nodes, D.lat.values, D.lon.values); D["s"] = s - ss0
    ov = pd.read_csv(f"{S3}/cross_region_overlaps.csv")
    joint = []
    used_s, used_n = set(), set()
    for _, r in ov.iterrows():
        joint.append((f"SVI{r.svi_episode}+NW{r.nw_episode}", pd.concat([Sv[Sv.ep == r.svi_episode], Nv[Nv.ep == r.nw_episode]])))
        used_s.add(r.svi_episode); used_n.add(r.nw_episode)
    for e in big["SVI"]:
        if e not in used_s:
            joint.append((f"SVI{e}", Sv[Sv.ep == e]))
    for e in big["NW"]:
        if e not in used_n:
            joint.append((f"NW{e}", Nv[Nv.ep == e]))
    unmerged = [(f"SVI{e}", Sv[Sv.ep == e]) for e in big["SVI"]] + [(f"NW{e}", Nv[Nv.ep == e]) for e in big["NW"]]
    fs = pd.read_csv(f"{S3}/families_svi_frame.csv"); sn = pd.read_csv(f"{S3}/sources_nw.csv").query("quality in ['A','B']")
    fpos = np.r_[qb.project(nodes, fs.lat.values, fs.lon.values)[0], qb.project(nodes, sn.lat.values, sn.lon.values)[0]] - ss0
    lo, hi = fpos.min(), fpos.max()
    for lab, L in (("merged", joint), ("unmerged", unmerged)):
        starts, pe = [], []
        rows = []
        for nm, g in L:
            g = g.sort_values("time"); k = max(1, int(np.ceil(0.1 * len(g))))
            st = float(np.median(g.s.values[:k])); starts.append(st)
            meds = np.array([np.median(rng.choice(g.s.values, k, replace=False)) for _ in range(500)])
            pe.append(float((np.abs(meds) <= 15).mean()))
            rows.append(dict(episode=nm, n=len(g), start_s=st, p05=float(np.quantile(g.s, .05)), p95=float(np.quantile(g.s, .95))))
        starts = np.array(starts); kk = int((np.abs(starts) <= 15).sum()); n = len(starts)
        pu = 30 / (hi - lo); pm = float((np.abs(fpos) <= 15).mean())
        r = dict(n=n, k_strait=kk, k_strait_hw25=int((np.abs(starts) <= 25).sum()), uniform_expect=n * pu, matched_expect=n * pm,
                 own_expect=float(np.sum(pe)),
                 p_uniform=float(1 - __import__("scipy").stats.binom.cdf(kk - 1, n, pu)),
                 p_matched=float(1 - __import__("scipy").stats.binom.cdf(kk - 1, n, pm)),
                 p_own=qt.poisson_binomial_sf(np.array(pe), kk), family_span=[float(lo), float(hi)], episodes=rows)
        res[f"q62_lfe_{lab}"] = r
        print(lab, {k: v for k, v in r.items() if k != "episodes"}, flush=True)
    json.dump(res, open(f"{OUT}/lfe_checks.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
