"""Q2: is the active band behind the ETS slip front the same in every episode, 2003-2026 (SVI + Washington-NW)?
Rules: svi_science/results/q2_band/PLAN.md (frozen before running).

Writes svi_science/results/q2_band/: episodes.csv, units.csv, controls.csv, thinning.csv, s13.csv, tremor.csv,
result.json, q2.log
"""
import glob
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
R = f"{ROOT}/svi_science/results"
OUT = f"{R}/q2_band"
sys.path.insert(0, f"{ROOT}/svi_science/scripts")
import s3_build as S3B  # noqa: E402  (functions only)

rng = np.random.default_rng(20261007)
LOG = open(f"{OUT}/q2.log", "w")
MAXN = 1500
NSHUF = 500
NBOOT = 200
NCTRL = 200
NPERM = 5000


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def to_pt(t):
    return pd.Series(pd.to_datetime(t).values).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.tz_localize(None).values


# ================================================================== loading
K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv")
EA = pd.read_csv(f"{R}/s3/episodes_all.csv")
FAM = pd.read_csv(f"{R}/s0/families.csv")
SRC = pd.read_csv(f"{R}/s3/sources_nw.csv")
CUT = {"SVI": tuple(FAM.downdip_km.quantile([1 / 3, 2 / 3])), "NW": tuple(SRC.dd_km.quantile([1 / 3, 2 / 3]))}


def dclass(dd, reg):
    a, b = CUT[reg]
    return np.where(dd >= b, "deep", np.where(dd <= a, "shallow", "mid"))


def std_cols(E, region, era):
    E = E.copy()
    E["time"] = pd.to_datetime(E.time).astype("datetime64[ns]")
    E["region"] = region; E["era"] = era
    E["dcl"] = dclass(E.dd_km.values, region)
    if region == "SVI":
        E["patch"] = "G" + E.group.astype(int).astype(str)
    else:
        E["patch"] = "B" + (np.floor(E.as_km / 10) * 10).astype(int).astype(str) + "_" + E.dcl
    if "mag" not in E:
        E["mag"] = np.nan
    return E[["time", "ep", "unit", "as_km", "dd_km", "dcl", "patch", "mag", "region", "era"]].dropna(subset=["as_km"]).sort_values("time").reset_index(drop=True)


def load_svi_o():
    eps = K[K.region == "SVI"].episode.tolist()
    E = pd.read_parquet(f"{R}/s3/events_svi.parquet")
    E = E[(E.series == "data") & (E.P_same >= 0.75) & E.counted_once & E.ep.isin(eps)].copy()
    E["unit"] = "F" + E.unit.astype(str)
    return std_cols(E, "SVI", "SVI-O 2014-26")


def load_svi_b():
    B = pd.read_parquet(f"{R}/p3_bostock/bostock_events_mapped.parquet")
    EP = pd.read_csv(f"{R}/p3_bostock/episodes_bostock_final.csv")
    B = B.rename(columns={"time_utc": "time", "along_km": "as_km", "downdip_km": "dd_km"})
    B["time"] = pd.to_datetime(B.time).astype("datetime64[ns]")
    ep = np.full(len(B), None, dtype=object); tv = B.time.values
    for r in EP.itertuples():
        m = (tv >= np.datetime64(pd.Timestamp(r.start) - pd.Timedelta(days=2))) & (tv < np.datetime64(pd.Timestamp(r.end) + pd.Timedelta(days=3)))
        ep[m] = r.episode
    B["ep"] = ep; B = B[B.ep.notna()].copy()
    B["unit"] = "F" + B.fid.astype(int).astype(str)
    log(f"SVI-B events in episodes {len(B)} of mapped {len(tv)}")
    return std_cols(B, "SVI", "SVI-B 2003-13")


def load_nw_o():
    eps = K[K.region == "NW"].episode.tolist()
    E = pd.read_parquet(f"{R}/s3/events_nw.parquet")
    E = E[(E.series == "data") & (E.P_same >= 0.75) & E.selected_v2 & E.counted_once & E.ep.isin(eps)].copy()
    E["unit"] = E.source
    return std_cols(E, "NW", "NW-O 2015-26")


def load_nw_e():
    Bep = pd.read_csv(f"{R}/b_ets_record/episodes_b.csv")
    Bep = Bep[(Bep.region == "NW") & Bep.big & (pd.to_datetime(Bep.start) < "2015-01-01")].copy()
    C = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet",
                        columns=["det_id", "time", "region", "family", "P_own_pnet", "selected_v2"],
                        filters=[("region", "==", "NW"), ("selected_v2", "==", True)])
    C = C[(C.time < "2015-01-01") & (C.P_own_pnet >= 0.75)].copy()
    fl = []
    for f in sorted(glob.glob(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_20*.parquet")):
        if f.split("_")[-1][:7] >= "2015-01":
            continue
        q = pd.read_parquet(f, columns=["det_id", "series", "dup_D"])
        fl.append(q[(q.series == "data") & q.det_id.isin(C.det_id)][["det_id", "dup_D"]])
    fl = pd.concat(fl).drop_duplicates("det_id")
    C = C.merge(fl, on="det_id", how="left")
    log(f"NW-E selected P>=0.75 rows {len(C)}; without flag {int(C.dup_D.isna().sum())}")
    C = C[~C.dup_D.fillna(False).astype(bool)].copy()
    ep = np.full(len(C), None, dtype=object); tv = C.time.values.astype("datetime64[ns]")
    for r in Bep.itertuples():
        m = (tv >= np.datetime64(pd.Timestamp(r.w0))) & (tv < np.datetime64(pd.Timestamp(r.w1)))
        ep[m] = r.episode
    C["ep"] = ep; C = C[C.ep.notna()].copy()
    SG = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3.csv")
    AC = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3_array_confirmed.csv").rename(columns={"group": "source"})
    fid2src = SG.merge(AC, on="fid").set_index("fid").source
    C["unit"] = C.family.map(fid2src)
    log(f"NW-E counted-once rows in big windows {len(C)}; no source {int(C.unit.isna().sum())}")
    C = C.merge(SRC[["source", "as_km", "dd_km"]].rename(columns={"source": "unit"}), on="unit", how="left")
    cnt = C.groupby("ep").size()
    keep = cnt[cnt >= 300].index.tolist()  # S3 "big" rule
    log(f"NW-E episode counts {cnt.to_dict()} -> kept {keep}")
    return std_cols(C[C.ep.isin(keep)], "NW", "NW-E 2007-14")


def load_svi_s13():
    eps = K[K.region == "SVI"].episode.tolist()
    T = pd.read_parquet(f"{R}/s0/analysis_table_S13.parquet", columns=["series", "fam", "P_same", "keep_xfam8", "time_utc"])
    T = T[(T.series == "data") & (T.P_same >= 0.75) & T.keep_xfam8.astype(bool)].copy()
    ep = np.full(len(T), None, dtype=object); tv = T.time_utc.values
    for _, r in EA[(EA.region == "SVI") & EA.episode.isin(eps)].iterrows():
        m = (tv >= np.datetime64(pd.Timestamp(r.w0))) & (tv < np.datetime64(pd.Timestamp(r.w1)))
        ep[m] = r.episode
    T["ep"] = ep; T = T[T.ep.notna()]
    T = T.merge(FAM[["fid", "along_strike_km", "downdip_km", "group"]].rename(columns={"fid": "fam"}), on="fam")
    T = T.rename(columns={"time_utc": "time", "along_strike_km": "as_km", "downdip_km": "dd_km"})
    T["unit"] = "F" + T.fam.astype(str)
    return std_cols(T, "SVI", "SVI-O S13")


def load_tremor(EV):
    """PNSN epicentres within 10 km of the sources, per episode window (first..last LFE ± 2 d), as s3_build."""
    P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
    info = json.load(open(f"{R}/s3/build_info.json"))
    frs = json.load(open(f"{R}/s0/frame_groups.json"))
    fr = {"SVI": dict(origin_lat=frs["origin_lat"], origin_lon=frs["origin_lon"], downdip_az=frs["mean_downdip_az"], strike_az=frs["strike_az"]),
          "NW": info["nw_frame"]}
    good = SRC.quality.isin(["A", "B"])
    pos = {"SVI": (FAM.lat.values, FAM.lon.values), "NW": (SRC.lat[good].values, SRC.lon[good].values)}
    out = []
    for reg in ["SVI", "NW"]:
        m = S3B.near_mask(P.lat.values, P.lon.values, pos[reg][0], pos[reg][1], 10.0)
        Pr = P[m].copy()
        a, d, _, _ = S3B.project(fr[reg], Pr.lat.values, Pr.lon.values); Pr["as_km"] = a; Pr["dd_km"] = d
        for (era, ep), g in EV[EV.region == reg].groupby(["era", "ep"]):
            t0 = g.time.min() - pd.Timedelta(days=2); t1 = g.time.max() + pd.Timedelta(days=2)
            if t1 < P.time.min() + pd.Timedelta(days=5):
                continue
            x = Pr[(Pr.time >= t0) & (Pr.time <= t1)].copy()
            if len(x) < 50:
                continue
            x["ep"] = ep; x["era"] = era; x["region"] = reg
            out.append(x)
    T = pd.concat(out, ignore_index=True)
    T["unit"] = "T" + (np.floor(T.as_km / 5) * 5).astype(int).astype(str)
    return T


# ================================================================== measurements
def ts_slope(t, a):
    t = np.asarray(t, float); a = np.asarray(a, float)
    i, j = np.triu_indices(len(t), 1)
    dt = t[j] - t[i]; da = a[j] - a[i]
    m = (np.abs(dt) > 1e-9) & (np.abs(da) > 1e-6)
    s = da[m] / dt[m]
    return float(np.median(s)) if len(s) else np.nan


def sub_idx(n, k):
    return np.sort(rng.choice(n, k, replace=False)) if n > k else np.arange(n)


def speed(td, a, full=True):
    idx = sub_idx(len(td), MAXN); ts, as_ = td[idx], a[idx]
    v = ts_slope(ts, as_)
    if not full:
        return dict(v=v)
    ii, jj = np.triu_indices(len(ts), 1); dt = ts[jj] - ts[ii]; mk = np.abs(dt) > 1e-9; ii, jj, dt = ii[mk], jj[mk], dt[mk]
    vs = np.empty(NSHUF)
    for k in range(NSHUF):
        ap = as_[rng.permutation(len(ts))]; da = ap[jj] - ap[ii]; mm = np.abs(da) > 1e-6
        vs[k] = np.median(da[mm] / dt[mm])
    p = float((np.sum(np.abs(vs) >= abs(v)) + 1) / (NSHUF + 1))
    day = np.floor(ts).astype(int); ud = np.unique(day); grp = {d: np.where(day == d)[0] for d in ud}
    bs = []
    for _ in range(NBOOT):
        ii2 = np.concatenate([grp[d] for d in rng.choice(ud, len(ud), replace=True)])
        bs.append(ts_slope(ts[ii2] + rng.normal(0, 1e-6, len(ii2)), as_[ii2]))
    bs = np.array(bs)
    return dict(v=v, p_speed=p, v_ci05=float(np.quantile(bs, 0.05)), v_ci95=float(np.quantile(bs, 0.95)),
                sig_logv=float(np.std(np.log(np.clip(np.abs(bs), 1e-3, None)))))


def arrival_speed(td, a, nmin=10):
    b = np.floor(a / 5) * 5 + 2.5
    rows = [(bc, np.quantile(td[b == bc], 0.1)) for bc in np.unique(b) if np.sum(b == bc) >= nmin]
    if len(rows) < 3:
        return np.nan
    x = np.array(rows)
    return float(stats.theilslopes(x[:, 0], x[:, 1])[0])


def d80_sub(t, k, nsub):
    """median over nsub random k-subsamples of (q90 - q10); vectorised."""
    n = len(t)
    if n < k:
        return np.nan
    if n == k:
        return float(np.quantile(t, 0.9) - np.quantile(t, 0.1))
    idx = np.argsort(rng.random((nsub, n)), axis=1)[:, :k]
    s = t[idx]
    return float(np.median(np.quantile(s, 0.9, axis=1) - np.quantile(s, 0.1, axis=1)))


def unit_durations(td, units, k=30, nsub=200):
    out = {}
    for u in np.unique(units):
        t = td[units == u]
        if len(t) >= k:
            out[u] = d80_sub(t, k, nsub)
    return out


def w_snap(tpt_d, a, k=15, nsub=100):
    """12-h PT windows with >= k events; extent = median over k-subsamples of q90-q10 position; median over windows."""
    w = np.floor(tpt_d * 2).astype(int)
    ex = []
    for wi in np.unique(w):
        x = a[w == wi]
        if len(x) >= k:
            ex.append(d80_sub(x, k, nsub))
    return (float(np.median(ex)) if ex else np.nan), len(ex)


def episode_measure(g, full=True, nsub_d=200, nsub_s=100):
    """g: events of one episode (one region/era). Returns dict + unit table."""
    t0 = g.time.min()
    td = (g.time - t0).dt.total_seconds().values / 86400
    tpt = to_pt(g.time.values); tpt_d = (pd.to_datetime(tpt) - pd.Timestamp("2000-01-01")).total_seconds().values / 86400
    a = g.as_km.values; units = g.unit.values
    r = dict(n=len(g), year=float(pd.Timestamp(g.time.quantile(0.1)).year + (pd.Timestamp(g.time.quantile(0.1)).dayofyear - 1) / 365.25),
             extent_km=float(np.quantile(a, 0.95) - np.quantile(a, 0.05)), dur_ep_d=float(np.quantile(td, 0.9) - np.quantile(td, 0.1)),
             n_units=int(len(np.unique(units))))
    r.update(speed(td, a, full=full))
    r["v_arrival"] = arrival_speed(td, a)
    D = unit_durations(td, units, 30, nsub_d)
    D15 = unit_durations(td, units, 15, nsub_d)
    P = unit_durations(td, g.patch.values, 30, nsub_d)
    r["n_units30"] = len(D); r["n_units15"] = len(D15)
    r["medD80"] = float(np.median(list(D.values()))) if D else np.nan
    r["medD80_n15"] = float(np.median(list(D15.values()))) if D15 else np.nan
    r["medD80_patch"] = float(np.median(list(P.values()))) if P else np.nan
    av = abs(r["v"])
    r["W_dur"] = av * r["medD80"]; r["W_dur15"] = av * r["medD80_n15"]; r["W_dur_patch"] = av * r["medD80_patch"]
    r["W_dur_arr"] = abs(r["v_arrival"]) * r["medD80"] if np.isfinite(r["v_arrival"]) else np.nan
    r["W_snap"], r["n_snapwin"] = w_snap(tpt_d, a, 15, nsub_s)
    # front-frame width
    res = a - r["v"] * td
    r["W_front"] = float(np.quantile(res, 0.9) - np.quantile(res, 0.1))
    # sigma of log median D80 by resampling units
    if len(D) >= 2 and full:
        vals = np.array(list(D.values()))
        bm = np.median(vals[rng.integers(0, len(vals), (500, len(vals)))], axis=1)
        r["sig_logD"] = float(np.std(np.log(bm)))
    else:
        r["sig_logD"] = np.nan
    if full:
        r["sig_logW"] = float(np.hypot(r.get("sig_logv", np.nan), r["sig_logD"]))
    # per depth third
    dcl = g.groupby("unit").dcl.first(); ddu = g.groupby("unit").dd_km.mean()
    for c in ["shallow", "mid", "deep"]:
        vals = [D[u] for u in D if dcl[u] == c]
        r[f"medD80_{c}"] = float(np.median(vals)) if vals else np.nan
        r[f"nU_{c}"] = len(vals)
        r[f"W_dur_{c}"] = av * r[f"medD80_{c}"]
    U = pd.DataFrame(dict(unit=list(D.keys()), D80_n30=list(D.values())))
    U["dd"] = U.unit.map(ddu); U["dcl"] = U.unit.map(dcl); U["as_km"] = U.unit.map(g.groupby("unit").as_km.mean())
    U["n"] = U.unit.map(g.groupby("unit").size())
    return r, U


def control(g, v, nset=NCTRL):
    """positions (unit labels with their positions) permuted among events within the episode; front speed kept."""
    t0 = g.time.min(); td = (g.time - t0).dt.total_seconds().values / 86400
    tpt = to_pt(g.time.values); tpt_d = (pd.to_datetime(tpt) - pd.Timestamp("2000-01-01")).total_seconds().values / 86400
    units = g.unit.values; a = g.as_km.values
    W, S = [], []
    for _ in range(nset):
        p = rng.permutation(len(g))
        D = unit_durations(td, units[p], 30, 20)
        W.append(abs(v) * np.median(list(D.values())) if D else np.nan)
        S.append(w_snap(tpt_d, a[p], 15, 10)[0])
    return np.array(W), np.array(S)


# ================================================================== stats helpers
def cochran(logw, sig):
    m = np.isfinite(logw) & np.isfinite(sig) & (sig > 0)
    x, s = logw[m], sig[m]
    if len(x) < 3:
        return dict(n=int(len(x)))
    w = 1 / s ** 2; mu = np.sum(w * x) / np.sum(w)
    Q = float(np.sum(w * (x - mu) ** 2)); df = len(x) - 1
    tau2 = max(0.0, (Q - df) / (np.sum(w) - np.sum(w ** 2) / np.sum(w)))
    return dict(n=int(len(x)), Q=Q, df=df, p=float(stats.chi2.sf(Q, df)), I2=float(max(0, (Q - df) / Q)) if Q > 0 else 0.0,
                tau_pct=float(100 * np.sqrt(tau2)), median_sig_pct=float(100 * np.median(s)), mean_W=float(np.exp(mu)))


def trend(year, w, strata=None, nperm=10000):
    m = np.isfinite(w) & np.isfinite(year)
    y, x = np.asarray(year)[m], np.asarray(w)[m]
    if len(x) < 4:
        return dict(n=int(len(x)))
    sl, ic, lo, hi = stats.theilslopes(x, y, 0.95)
    st = None if strata is None else np.asarray(strata)[m]
    obs = stats.spearmanr(y, x)[0]
    nul = []
    for _ in range(nperm):
        if st is None:
            xp = rng.permutation(x)
        else:
            xp = x.copy()
            for s in np.unique(st):
                k = st == s; xp[k] = rng.permutation(x[k])
        nul.append(stats.spearmanr(y, xp)[0])
    nul = np.array(nul)
    span = y.max() - y.min()
    return dict(n=int(len(x)), slope_km_per_yr=float(sl), slope_lo=float(lo), slope_hi=float(hi), span_yr=float(span),
                change_over_span_lo=float(lo * span), change_over_span_hi=float(hi * span), rho=float(obs),
                p=float((np.sum(np.abs(nul) >= abs(obs)) + 1) / (nperm + 1)))


def mw(a, b, nboot=5000):
    a = np.asarray(a)[np.isfinite(a)]; b = np.asarray(b)[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return dict(na=len(a), nb=len(b))
    p = float(stats.mannwhitneyu(a, b).pvalue)
    d = np.median(a) - np.median(b)
    bs = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(nboot)]
    return dict(na=int(len(a)), nb=int(len(b)), med_a=float(np.median(a)), med_b=float(np.median(b)), diff=float(d),
                ci_lo=float(np.quantile(bs, 0.025)), ci_hi=float(np.quantile(bs, 0.975)), p=p)


def strat_perm_corr(x, y, strata, nperm=10000):
    m = np.isfinite(x) & np.isfinite(y)
    x, y, st = np.asarray(x)[m], np.asarray(y)[m], np.asarray(strata)[m]
    if len(x) < 5:
        return dict(n=int(len(x)))
    obs = stats.spearmanr(x, y)[0]; nul = []
    for _ in range(nperm):
        yp = y.copy()
        for s in np.unique(st):
            k = st == s; yp[k] = rng.permutation(y[k])
        nul.append(stats.spearmanr(x, yp)[0])
    return dict(n=int(len(x)), rho=float(obs), p=float((np.sum(np.abs(np.array(nul)) >= abs(obs)) + 1) / (nperm + 1)))


def perm_within(U, ycol="D80_n30", xcol="dd", nperm=NPERM, min_cells=3):
    groups = [g for _, g in U.dropna(subset=[ycol, xcol]).groupby("ep") if len(g) >= min_cells]
    if not groups:
        return dict(n_ep=0)
    data = [(g[xcol].values, g[ycol].values) for g in groups]
    f = lambda ds: np.nanmean([stats.spearmanr(x, y)[0] for x, y in ds])
    obs = f(data)
    nul = np.array([f([(rng.permutation(x), y) for x, y in data]) for _ in range(nperm)])
    return dict(n_ep=len(groups), n_units=int(sum(len(g) for g in groups)), rho=float(obs),
                p=float((np.sum(np.abs(nul) >= abs(obs)) + 1) / (nperm + 1)))


# ================================================================== main
def main():
    res = {}
    EV = pd.concat([load_svi_b(), load_svi_o(), load_nw_e(), load_nw_o()], ignore_index=True)
    log("events per era: " + EV.groupby("era").size().to_string().replace("\n", "; "))
    rows, UU = [], []
    for (era, ep), g in EV.groupby(["era", "ep"]):
        r, U = episode_measure(g.reset_index(drop=True))
        r.update(region=g.region.iloc[0], era=era, ep=ep)
        rows.append(r); UU.append(U.assign(region=r["region"], era=era, ep=ep))
        log(f"{era} {ep}: n {r['n']} v {r['v']:+.2f} (p {r['p_speed']:.3f}) medD80 {r['medD80']:.2f} W_dur {r['W_dur']:.1f} W_snap {r['W_snap']:.1f} units30 {r['n_units30']}")
    EP = pd.DataFrame(rows)
    EP["single_front"] = (EP.p_speed < 0.01) & (EP.v.abs() >= 2)
    EP["direction"] = np.where(EP.v < 0, "SE-ward", "NW-ward")
    U = pd.concat(UU, ignore_index=True)
    U["W_unit"] = U.ep.map(dict(zip(zip(EP.era, EP.ep), EP.v.abs()))) if False else np.nan
    vmap = {(a, b): abs(v) for a, b, v in zip(EP.era, EP.ep, EP.v)}
    U["W_unit"] = [vmap[(a, b)] * d for a, b, d in zip(U.era, U.ep, U.D80_n30)]
    U.to_csv(f"{OUT}/units.csv", index=False)

    # ---------------- T1 control
    C = []
    for (era, ep), g in EV.groupby(["era", "ep"]):
        v = EP[(EP.era == era) & (EP.ep == ep)].v.iloc[0]
        cw, cs = control(g.reset_index(drop=True), v, NCTRL)
        C.append(pd.DataFrame(dict(era=era, ep=ep, set=np.arange(NCTRL), W_dur_ctrl=cw, W_snap_ctrl=cs)))
        log(f"control {era} {ep}: W_dur ctrl median {np.nanmedian(cw):.1f}, W_snap ctrl {np.nanmedian(cs):.1f}")
    C = pd.concat(C); C.to_csv(f"{OUT}/controls.csv", index=False)
    res["T1_control"] = {}
    for era, g in EP.groupby("era"):
        sf = g[g.single_front]
        cg = C[C.era == era]
        out = {}
        for col, sub in (("W_dur", sf), ("W_snap", g)):
            obs = float(np.nanmedian(sub[col]))
            cc = cg[cg.ep.isin(sub.ep)].pivot(index="set", columns="ep", values=f"{col}_ctrl")
            nul = np.nanmedian(cc.values, axis=1)
            out[col] = dict(obs=obs, ctrl_median=float(np.nanmedian(nul)), ctrl_q01=float(np.nanquantile(nul, 0.01)),
                            p=float((np.sum(nul <= obs) + 1) / (len(nul) + 1)), n_ep=int(len(sub)),
                            pass_=bool(obs < np.nanquantile(nul, 0.01)))
        res["T1_control"][era] = out
    for col in ("W_dur", "W_snap"):
        EP[f"{col}_ctrl_med"] = [np.nanmedian(C[(C.era == a) & (C.ep == b)][f"{col}_ctrl"]) for a, b in zip(EP.era, EP.ep)]

    # ---------------- T6 thinning (Bostock)
    target = int(EP[EP.era == "SVI-O 2014-26"].n.median())
    res["T6_target_count"] = target
    TH = []
    for ep, g in EV[EV.era == "SVI-B 2003-13"].groupby("ep"):
        g = g.reset_index(drop=True)
        for k in range(50):
            gi = g.iloc[np.sort(rng.choice(len(g), min(target, len(g)), replace=False))].reset_index(drop=True)
            r, _ = episode_measure(gi, full=False, nsub_d=50, nsub_s=30)
            TH.append(dict(ep=ep, kind="random", draw=k, **{c: r[c] for c in ("n", "v", "medD80", "W_dur", "W_snap", "n_units30", "n_snapwin")}))
        gm = g.sort_values("mag", ascending=False).iloc[:target].sort_values("time").reset_index(drop=True)
        r, _ = episode_measure(gm, full=False, nsub_d=200, nsub_s=100)
        TH.append(dict(ep=ep, kind="top_mag", draw=0, **{c: r[c] for c in ("n", "v", "medD80", "W_dur", "W_snap", "n_units30", "n_snapwin")}))
        log(f"thinning {ep}: random W_dur {np.nanmedian([x['W_dur'] for x in TH if x['ep'] == ep and x['kind'] == 'random']):.1f} "
            f"top_mag W_dur {TH[-1]['W_dur']:.1f} W_snap {TH[-1]['W_snap']:.1f}")
    TH = pd.DataFrame(TH); TH.to_csv(f"{OUT}/thinning.csv", index=False)
    thr = TH[TH.kind == "random"].groupby("ep")[["v", "medD80", "W_dur", "W_snap", "n_units30"]].median().add_suffix("_thinR")
    thm = TH[TH.kind == "top_mag"].set_index("ep")[["v", "medD80", "W_dur", "W_snap", "n_units30"]].add_suffix("_thinM")
    EP = EP.merge(thr, left_on="ep", right_index=True, how="left").merge(thm, left_on="ep", right_index=True, how="left")
    b = EP[EP.era == "SVI-B 2003-13"]
    res["T6_ratio"] = {k: float(np.nanmedian(b[f"{c}_{s}"] / b[c])) for c in ("v", "medD80", "W_dur", "W_snap") for s, k in
                       (("thinR", f"{c}_random"), ("thinM", f"{c}_topmag")) for k in [k]}

    # ---------------- T2 heterogeneity
    res["T2_heterogeneity"] = {}
    for key, g in list(EP.groupby("era")) + list(EP.groupby("region")):
        sf = g[g.single_front]
        res["T2_heterogeneity"][key] = dict(**cochran(np.log(sf.W_dur.values), sf.sig_logW.values),
                                            CV_W_dur=float(sf.W_dur.std() / sf.W_dur.mean()),
                                            CV_W_snap=float(g.W_snap.std() / g.W_snap.mean()),
                                            median_W_dur=float(sf.W_dur.median()), iqr_W_dur=[float(sf.W_dur.quantile(.25)), float(sf.W_dur.quantile(.75))],
                                            median_W_snap=float(g.W_snap.median()), iqr_W_snap=[float(g.W_snap.quantile(.25)), float(g.W_snap.quantile(.75))],
                                            median_v=float(sf.v.abs().median()), median_D80=float(sf.medD80.median()))

    # ---------------- T3 trend, T4 eras, T5 regions
    res["T3_trend"] = {}; res["T4_eras"] = {}; res["T5_regions"] = {}
    for reg, g in EP.groupby("region"):
        sf = g[g.single_front]
        res["T3_trend"][reg] = dict(W_dur=trend(sf.year.values, sf.W_dur.values), W_snap=trend(g.year.values, g.W_snap.values),
                                    medD80=trend(sf.year.values, sf.medD80.values), v=trend(sf.year.values, sf.v.abs().values))
        if reg == "SVI":
            # thinned-Bostock version of the trend
            wt = np.where(sf.era == "SVI-B 2003-13", sf.W_dur_thinM, sf.W_dur)
            st = np.where(g.era == "SVI-B 2003-13", g.W_snap_thinM, g.W_snap)
            res["T3_trend"]["SVI_bostock_topmag"] = dict(W_dur=trend(sf.year.values, wt), W_snap=trend(g.year.values, st))
            wt = np.where(sf.era == "SVI-B 2003-13", sf.W_dur_thinR, sf.W_dur)
            st = np.where(g.era == "SVI-B 2003-13", g.W_snap_thinR, g.W_snap)
            res["T3_trend"]["SVI_bostock_random"] = dict(W_dur=trend(sf.year.values, wt), W_snap=trend(g.year.values, st))
    pairs = [("SVI-B 2003-13", "SVI-O 2014-26"), ("NW-E 2007-14", "NW-O 2015-26")]
    for ea, eb in pairs:
        A, B_ = EP[EP.era == ea], EP[EP.era == eb]
        Asf, Bsf = A[A.single_front], B_[B_.single_front]
        d = dict(W_dur=mw(Asf.W_dur, Bsf.W_dur), W_snap=mw(A.W_snap, B_.W_snap), medD80=mw(Asf.medD80, Bsf.medD80),
                 v=mw(Asf.v.abs(), Bsf.v.abs()))
        if ea.startswith("SVI-B"):
            for s in ("thinR", "thinM"):
                d[f"W_dur_{s}"] = mw(Asf[f"W_dur_{s}"], Bsf.W_dur); d[f"W_snap_{s}"] = mw(A[f"W_snap_{s}"], B_.W_snap)
                d[f"medD80_{s}"] = mw(Asf[f"medD80_{s}"], Bsf.medD80)
        res["T4_eras"][f"{ea} vs {eb}"] = d
    sv, nw = EP[EP.region == "SVI"], EP[EP.region == "NW"]
    res["T5_regions"]["all"] = dict(W_dur=mw(sv[sv.single_front].W_dur, nw[nw.single_front].W_dur), W_snap=mw(sv.W_snap, nw.W_snap),
                                    medD80=mw(sv[sv.single_front].medD80, nw[nw.single_front].medD80))
    so, no = EP[EP.era == "SVI-O 2014-26"], EP[EP.era == "NW-O 2015-26"]
    res["T5_regions"]["2014-26"] = dict(W_dur=mw(so[so.single_front].W_dur, no[no.single_front].W_dur), W_snap=mw(so.W_snap, no.W_snap))

    # ---------------- T7 depth
    res["T7_depth"] = {}
    for era, g in U.groupby("era"):
        res["T7_depth"][era] = perm_within(g)
        e = EP[(EP.era == era) & EP.single_front]
        res["T7_depth"][era]["W_dur_by_third_median"] = {c: float(e[f"W_dur_{c}"].median()) for c in ("shallow", "mid", "deep")}
        res["T7_depth"][era]["medD80_by_third_median"] = {c: float(e[f"medD80_{c}"].median()) for c in ("shallow", "mid", "deep")}
        log(f"T7 {era}: {res['T7_depth'][era]}")

    # ---------------- T11 tremor (also gives size iii)
    TR = load_tremor(EV)
    trow = []
    for (era, ep), g in TR.groupby(["era", "ep"]):
        r, _ = episode_measure(g.sort_values("time").reset_index(drop=True).assign(patch="x", dcl="x", dd_km=0.0), full=True, nsub_d=100, nsub_s=50)
        trow.append(dict(era=era, ep=ep, n_trem=r["n"], v_trem=r["v"], p_speed_trem=r["p_speed"], medD80_trem=r["medD80"],
                         W_dur_trem=r["W_dur"], W_snap_trem=r["W_snap"], n_bins30=r["n_units30"]))
    TT = pd.DataFrame(trow); TT.to_csv(f"{OUT}/tremor.csv", index=False)
    EP = EP.merge(TT, on=["era", "ep"], how="left")
    m = EP.single_front & EP.W_dur_trem.notna() & (EP.p_speed_trem < 0.01)
    res["T11_tremor"] = dict(n_ep=int(m.sum()),
                             rho_W_dur=float(stats.spearmanr(EP[m].W_dur, EP[m].W_dur_trem)[0]),
                             p_W_dur=float(stats.spearmanr(EP[m].W_dur, EP[m].W_dur_trem)[1]),
                             median_ratio_W_dur=float(np.median(EP[m].W_dur_trem / EP[m].W_dur)),
                             rho_v=float(stats.spearmanr(EP[m].v, EP[m].v_trem)[0]),
                             median_ratio_D80=float(np.median(EP[m].medD80_trem / EP[m].medD80)),
                             n_snap=int((EP.W_snap_trem.notna() & EP.W_snap.notna()).sum()),
                             rho_W_snap=float(stats.spearmanr(EP.W_snap, EP.W_snap_trem, nan_policy="omit")[0]),
                             median_ratio_W_snap=float(np.nanmedian(EP.W_snap_trem / EP.W_snap)),
                             by_region_W_dur_trem={reg: float(EP[m & (EP.region == reg)].W_dur_trem.median()) for reg in ("SVI", "NW")})
    for reg in ("SVI", "NW"):
        mm = m & (EP.region == reg)
        res["T11_tremor"][f"trend_W_dur_trem_{reg}"] = trend(EP[mm].year.values, EP[mm].W_dur_trem.values)

    # ---------------- T8 size, T9 direction
    EP["lfe_rank_in_era"] = EP.groupby("era").n.rank(pct=True)
    res["T8_size"] = {}; res["T9_direction"] = {}
    for reg, g in EP.groupby("region"):
        sf = g[g.single_front]
        res["T8_size"][reg] = {f"{w}_vs_{s}": strat_perm_corr(gg[s].values, gg[w].values, gg.era.values)
                               for w, gg in (("W_dur", sf), ("W_snap", g)) for s in ("extent_km", "dur_ep_d", "n_trem", "lfe_rank_in_era")}
        res["T8_size"][reg]["medD80_vs_v"] = strat_perm_corr(sf.v.abs().values, sf.medD80.values, sf.era.values)
        res["T8_size"][reg]["medD80_vs_v_by_era"] = {e: dict(rho=float(stats.spearmanr(x.v.abs(), x.medD80)[0]), p=float(stats.spearmanr(x.v.abs(), x.medD80)[1]), n=int(len(x)))
                                                     for e, x in sf.groupby("era") if len(x) >= 4}
        se, nwd = sf[sf.v < 0], sf[sf.v > 0]
        obs = np.median(se.W_dur) - np.median(nwd.W_dur) if len(se) and len(nwd) else np.nan
        lab = (sf.v < 0).values; nul = []
        for _ in range(10000):
            lp = lab.copy()
            for e in sf.era.unique():
                k = (sf.era == e).values; lp[k] = rng.permutation(lab[k])
            nul.append(np.median(sf.W_dur.values[lp]) - np.median(sf.W_dur.values[~lp]) if lp.any() and (~lp).any() else np.nan)
        nul = np.array(nul)
        res["T9_direction"][reg] = dict(n_SE=int(len(se)), n_NW=int(len(nwd)), med_SE=float(np.median(se.W_dur)) if len(se) else None,
                                        med_NW=float(np.median(nwd.W_dur)) if len(nwd) else None, diff=float(obs),
                                        p=float((np.sum(np.abs(nul[np.isfinite(nul)]) >= abs(obs)) + 1) / (np.isfinite(nul).sum() + 1)),
                                        W_snap=mw(g[g.v < 0].W_snap, g[g.v > 0].W_snap))

    # ---------------- T10 fixed stations
    S13 = load_svi_s13(); srow = []
    for ep, g in S13.groupby("ep"):
        r, _ = episode_measure(g.reset_index(drop=True), full=True)
        srow.append(dict(ep=ep, n_s13=r["n"], v_s13=r["v"], medD80_s13=r["medD80"], W_dur_s13=r["W_dur"], W_snap_s13=r["W_snap"], sig_logW_s13=r["sig_logW"]))
    S = pd.DataFrame(srow).merge(EP[EP.era == "SVI-O 2014-26"][["ep", "W_dur", "W_snap", "sig_logW", "single_front"]], on="ep")
    S["agree"] = np.abs(np.log(S.W_dur_s13 / S.W_dur)) <= np.hypot(S.sig_logW, S.sig_logW_s13)
    S.to_csv(f"{OUT}/s13.csv", index=False)
    res["T10_s13"] = dict(n=int(len(S)), n_agree=int(S.agree.sum()), median_ratio_W_dur=float(np.median(S.W_dur_s13 / S.W_dur)),
                          median_ratio_W_snap=float(np.nanmedian(S.W_snap_s13 / S.W_snap)))

    EP.to_csv(f"{OUT}/episodes.csv", index=False)
    json.dump(res, open(f"{OUT}/result.json", "w"), indent=1, default=float)
    log(json.dumps(res, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main()
