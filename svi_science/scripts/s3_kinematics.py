"""S3 kinematics: per big ETS episode start / end / along-strike migration / secondary bursts, LFE vs PNSN tremor,
time-shuffle control, detectability (thinning) check, repeatability tests. Pre-registration: results/s3/RESULTS_S3.md §0.

Writes results/s3/episodes_kinematics.csv, bursts.csv, family_onsets.csv, repeatability.json.
"""
import json
import zlib
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
TZ = "America/Los_Angeles"
NSHUF = 500
NBOOT = 200
MAXN = 1500
BIN = 5.0
rng0 = np.random.default_rng(20261007)


def ts_slope(t, a):
    """Theil-Sen slope (pairwise median) — fast numpy version, pairs with dt > 0 only."""
    t = np.asarray(t, float); a = np.asarray(a, float)
    i, j = np.triu_indices(len(t), 1)
    dt = t[j] - t[i]; da = a[j] - a[i]
    m = (np.abs(dt) > 1e-9) & (np.abs(da) > 1e-6)     # pairs at the same position (same family/source) carry no migration information
    s = da[m] / dt[m]
    return float(np.median(s)) if len(s) else np.nan


def sub(n, rng):
    return np.sort(rng.choice(n, MAXN, replace=False)) if n > MAXN else np.arange(n)


def load():
    ev = {}
    S = pd.read_parquet(f"{OUT}/events_svi.parquet")
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")[["fid", "slab2_depth_km"]].rename(columns={"fid": "unit"})
    S = S.merge(F, on="unit", how="left")
    S["selected"] = True
    ev["SVI"] = S
    N = pd.read_parquet(f"{OUT}/events_nw.parquet")
    src = pd.read_csv(f"{OUT}/sources_nw.csv")
    import importlib.util
    spec = importlib.util.spec_from_file_location("b", f"{ROOT}/svi_science/scripts/s3_build.py"); b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
    lons, lats, Z = b.slab_grid()
    src["slab2_depth_km"] = [-b._bilin(lons, lats, Z, lo, la) for lo, la in zip(src.lon, src.lat)]
    src.to_csv(f"{OUT}/sources_nw.csv", index=False)
    N = N.merge(src[["source", "slab2_depth_km"]].rename(columns={"source": "unit"}), on="unit", how="left")
    N["selected"] = N.selected_v2
    ev["NW"] = N
    tr = {"SVI": pd.read_parquet(f"{OUT}/tremor_svi.parquet"), "NW": pd.read_parquet(f"{OUT}/tremor_nw.parquet")}
    EP = pd.read_csv(f"{OUT}/episodes_all.csv", parse_dates=["w0", "w1"])
    return ev, tr, EP


def headline(g, variant="once075"):
    d = g[g.series == "data"]
    if variant == "once075":
        return d[(d.P_same >= 0.75) & d.selected & d.counted_once]
    if variant == "raw075":
        return d[(d.P_same >= 0.75) & d.selected]
    if variant == "once050":
        return d[d.counted_once]
    raise ValueError


def days(t, w0):
    return (pd.to_datetime(t).values - np.datetime64(w0)) / np.timedelta64(1, "D")


def basic(t, a):
    """start10 / end10 / first-2-day start / extent from times (days) and positions."""
    o = np.argsort(t); t = t[o]; a = a[o]; n = len(t)
    k = max(1, int(np.ceil(0.1 * n)))
    t1 = np.quantile(t, 0.01)
    m2 = (t >= t1) & (t <= t1 + 2)
    return dict(start10=float(np.median(a[:k])), end10=float(np.median(a[-k:])), start2d=float(np.median(a[m2])),
                n_first2d=int(m2.sum()), p05=float(np.quantile(a, 0.05)), p95=float(np.quantile(a, 0.95)),
                t_first=float(t1), t_last=float(np.quantile(t, 0.99)), t10=float(t[k - 1]), t90=float(t[-k]))


def kin(t, a, rng, nshuf=NSHUF, boot=True):
    """Theil-Sen speed + time-shuffle p + day-block bootstrap CI + start-offset shuffle p."""
    idx = sub(len(t), rng); ts, as_ = t[idx], a[idx]
    v = ts_slope(ts, as_)
    b = basic(t, a)
    so = b["start10"] - b["end10"]
    vs, ss = np.empty(nshuf), np.empty(nshuf)
    ii, jj = np.triu_indices(len(ts), 1); dt = ts[jj] - ts[ii]; mk = np.abs(dt) > 1e-9; ii, jj, dt = ii[mk], jj[mk], dt[mk]
    for i in range(nshuf):
        p = rng.permutation(len(ts)); ap = as_[p]; da = ap[jj] - ap[ii]; mm = np.abs(da) > 1e-6
        vs[i] = float(np.median(da[mm] / dt[mm]))
        bb = basic(ts, as_[p]); ss[i] = bb["start10"] - bb["end10"]
    p_v = float((np.sum(np.abs(vs) >= abs(v)) + 1) / (nshuf + 1))
    p_s = float((np.sum(np.abs(ss) >= abs(so)) + 1) / (nshuf + 1))
    out = dict(speed=v, speed_shuf_q005=float(np.quantile(vs, 0.005)), speed_shuf_q995=float(np.quantile(vs, 0.995)), p_speed=p_v,
               start_minus_end=so, p_start_end=p_s, **b)
    if boot:
        day = np.floor(ts).astype(int); ud = np.unique(day); groups = {d: np.where(day == d)[0] for d in ud}
        bs = []
        for _ in range(NBOOT):
            pick = rng.choice(ud, len(ud), replace=True)
            ii = np.concatenate([groups[d] for d in pick])
            # break exact ties in t from duplicated days by tiny jitter
            bs.append(ts_slope(ts[ii] + rng.normal(0, 1e-6, len(ii)), as_[ii]))
        out["speed_ci05"], out["speed_ci95"] = [float(x) for x in np.quantile(bs, [0.05, 0.95])]
    return out


def bursts(day, a, lo, hi, thr=3, quiet=3):
    """Returns: list of (bin_center, return_day, n_on_return_day, first_active_day)."""
    edges = np.arange(lo, hi + BIN, BIN)
    bi = np.clip(np.digitize(a, edges) - 1, 0, len(edges) - 2)
    D = np.arange(day.min(), day.max() + 1)
    M = np.zeros((len(edges) - 1, len(D)), int)
    np.add.at(M, (bi, day - D[0]), 1)
    res = []
    for b_ in range(M.shape[0]):
        row = M[b_]; act = np.where(row >= thr)[0]
        if len(act) == 0:
            continue
        first = act[0]; zeros = 0; seen_gap = False
        for d in range(first + 1, len(row)):
            if row[d] == 0:
                zeros += 1
            else:
                if zeros >= quiet and row[d: d + 2].max() >= thr:
                    res.append((float(edges[b_] + BIN / 2), int(D[d]), int(row[d]), int(D[first]), zeros))
                zeros = 0 if row[d] > 0 else zeros
    return res, M, edges, D


def pt_day(t, w0):
    """integer PT calendar day index relative to the window start's PT date."""
    tp = pd.to_datetime(t).dt.tz_localize("UTC").dt.tz_convert(TZ).dt.tz_localize(None).dt.normalize()
    w0p = pd.Timestamp(w0).tz_localize("UTC").tz_convert(TZ).tz_localize(None).normalize()
    return ((tp - w0p) / pd.Timedelta(days=1)).astype(int).to_numpy()


def family_onsets(h, w0, t_ep_first):
    rows = []
    for u, g in h.groupby("unit"):
        t = np.sort(days(g.time, w0))
        if len(t) < 5:
            continue
        rows.append(dict(unit=u, n=len(t), onset=float(np.quantile(t, 0.1)) - t_ep_first, as_km=float(g.as_km.iloc[0])))
    return pd.DataFrame(rows)


def thin(h, w0, rng, reps=50):
    cnt = h.groupby("unit").size(); cap = int(np.median(cnt))
    st, sp = [], []
    for _ in range(reps):
        parts = [g if len(g) <= cap else g.iloc[rng.choice(len(g), cap, replace=False)] for _, g in h.groupby("unit")]
        x = pd.concat(parts); t = days(x.time, w0); a = x.as_km.to_numpy()
        o = np.argsort(t); t, a = t[o], a[o]
        st.append(basic(t, a)["start10"]); idx = sub(len(t), rng); sp.append(ts_slope(t[idx], a[idx]))
    return cap, float(np.median(st)), float(np.median(sp)), float(np.mean(np.sign(sp) == np.sign(np.median(sp))))


def main():
    ev, tr, EP = load()
    rows, brows, orows = [], [], []
    tracks = {}
    for reg in ["SVI", "NW"]:
        E = ev[reg]; T = tr[reg]
        big = EP[(EP.region == reg) & EP.big].sort_values("w0")
        hall = headline(E[E.ep.isin(big.episode)])
        lo, hi = np.floor(np.quantile(hall.as_km, 0.0) / BIN) * BIN, np.ceil(np.quantile(hall.as_km, 1.0) / BIN) * BIN
        for _, r in big.iterrows():
            rng = np.random.default_rng(zlib.crc32(f"{reg}{r.episode}".encode()))
            g = E[E.ep == r.episode]; h = headline(g).sort_values("time")
            t = days(h.time, r.w0); a = h.as_km.to_numpy(); dd = h.dd_km.to_numpy()
            K = kin(t, a, rng)
            Kd = kin(t, dd, rng, nshuf=200, boot=False)
            row = dict(region=reg, episode=r.episode, window_start_PT=pd.Timestamp(r.w0).tz_localize("UTC").tz_convert(TZ).strftime("%Y-%m-%d"),
                       window_end_PT=pd.Timestamp(r.w1).tz_localize("UTC").tz_convert(TZ).strftime("%Y-%m-%d"),
                       first_event_PT=(pd.Timestamp(r.w0) + pd.Timedelta(days=K["t_first"])).tz_localize("UTC").tz_convert(TZ).strftime("%Y-%m-%d %H:%M"),
                       n_headline=len(h), n_raw075=int(len(headline(g, "raw075"))), n_once050=int(len(headline(g, "once050"))),
                       n_c1_075=int(((g.series == "c1") & (g.P_same >= 0.75)).sum()), n_units=int(h.unit.nunique()))
            row.update({f"lfe_{k}": v for k, v in K.items()})
            row.update(lfe_dd_speed=Kd["speed"], lfe_dd_p=Kd["p_speed"], lfe_start10_dd=Kd["start10"], lfe_end10_dd=Kd["end10"],
                       lfe_duration_d=K["t_last"] - K["t_first"])
            # variants
            for var in ["raw075", "once050"]:
                hv = headline(g, var); tv = days(hv.time, r.w0); av = hv.as_km.to_numpy()
                bv = basic(tv, av); idx = sub(len(tv), rng)
                row[f"{var}_start10"] = bv["start10"]; row[f"{var}_end10"] = bv["end10"]; row[f"{var}_speed"] = ts_slope(tv[idx], av[idx])
            # front speed from family onsets
            fo = family_onsets(h, r.w0, K["t_first"])
            fo["region"] = reg; fo["episode"] = r.episode; orows.append(fo)
            if len(fo) >= 5:
                row["front_speed_onsets"] = ts_slope(fo.onset.to_numpy() + rng.normal(0, 1e-6, len(fo)), fo.as_km.to_numpy())
                rho, p = stats.spearmanr(fo.n, fo.onset); row["onset_vs_count_rho"] = rho; row["onset_vs_count_p"] = p
            cap, st_thin, sp_thin, sgn = thin(h, r.w0, rng)
            row.update(thin_cap=cap, thin_start10=st_thin, thin_speed=sp_thin, thin_sign_agree=sgn)
            # bursts (LFE)
            pdays = pt_day(h.time, r.w0)
            bl, M, edges, D = bursts(pdays, a, lo, hi)
            row["n_bursts"] = len(bl)
            for b_ in bl:
                brows.append(dict(region=reg, episode=r.episode, series="LFE", bin_center_km=b_[0], return_day_PT=(pd.Timestamp(r.w0).tz_localize("UTC").tz_convert(TZ).normalize() + pd.Timedelta(days=b_[1])).strftime("%Y-%m-%d"),
                                  n_on_return_day=b_[2], quiet_days_before=b_[4], first_active_day_PT=(pd.Timestamp(r.w0).tz_localize("UTC").tz_convert(TZ).normalize() + pd.Timedelta(days=b_[3])).strftime("%Y-%m-%d"),
                                  return_day_idx=b_[1]))
            nb = []
            for _ in range(200):
                nb.append(len(bursts(pdays[rng.permutation(len(pdays))], a, lo, hi)[0]))
            row["n_bursts_shuffle_mean"] = float(np.mean(nb)); row["p_bursts_ge"] = float((np.sum(np.array(nb) >= len(bl)) + 1) / 201)
            # tremor
            tg = T[T.ep == r.episode].sort_values("time")
            row["n_tremor"] = len(tg)
            if len(tg) >= 50:
                tt = days(tg.time, r.w0); ta = tg.as_km.to_numpy()
                KT = kin(tt, ta, rng, nshuf=200, boot=False)
                row.update({f"trem_{k}": v for k, v in KT.items() if k in ("speed", "p_speed", "start10", "end10", "start2d", "p05", "p95", "t_first", "t_last", "start_minus_end", "p_start_end")})
                tb, _, _, _ = bursts(pt_day(tg.time, r.w0), ta, lo, hi)
                row["trem_n_bursts"] = len(tb)
                for b_ in tb:
                    brows.append(dict(region=reg, episode=r.episode, series="tremor", bin_center_km=b_[0], n_on_return_day=b_[2], quiet_days_before=b_[4], return_day_idx=b_[1]))
            rows.append(row)
            print(reg, r.episode, len(h), f"start {K['start10']:.0f} end {K['end10']:.0f} v {K['speed']:.2f} p {K['p_speed']:.3f}",
                  f"| trem start {row.get('trem_start10', np.nan):.0f} v {row.get('trem_speed', np.nan):.2f} p {row.get('trem_p_speed', np.nan):.3f}",
                  f"| bursts {len(bl)} (shuf {row['n_bursts_shuffle_mean']:.1f}) | thin start {st_thin:.0f} v {sp_thin:.2f}", flush=True)
    K = pd.DataFrame(rows)
    K["migrates"] = K.lfe_p_speed < 0.01
    K["direction"] = np.where(K.lfe_speed > 0, "NW-ward", "SE-ward")
    K["trem_migrates"] = K.trem_p_speed < 0.01
    K["agree_tremor"] = (np.sign(K.lfe_speed) == np.sign(K.trem_speed)) & ((K.lfe_start10 - K.trem_start10).abs() <= 15)
    K["thin_ok"] = ((K.thin_start10 - K.lfe_start10).abs() < 10) & (np.sign(K.thin_speed) == np.sign(K.lfe_speed))
    K.to_csv(f"{OUT}/episodes_kinematics.csv", index=False)
    pd.DataFrame(brows).to_csv(f"{OUT}/bursts.csv", index=False)
    pd.concat(orows).to_csv(f"{OUT}/family_onsets.csv", index=False)
    repeatability(K, ev, EP)


def repeatability(K, ev, EP):
    rng = np.random.default_rng(7)
    out = {}
    for reg in ["SVI", "NW"]:
        k = K[K.region == reg].reset_index(drop=True)
        E = ev[reg]; hall = headline(E[E.ep.isin(k.episode)])
        lo, hi = np.quantile(hall.as_km, [0.05, 0.95])
        sd = float(np.std(k.lfe_start10, ddof=1))
        U = rng.uniform(lo, hi, (10000, len(k))); sdu = U.std(axis=1, ddof=1)
        # stricter control: start replaced by median of a random 10% of each episode's own events
        sdr = []
        hs = {e: headline(E[E.ep == e]).as_km.to_numpy() for e in k.episode}
        for _ in range(2000):
            sdr.append(np.std([np.median(rng.choice(hs[e], max(1, int(np.ceil(0.1 * len(hs[e])))), replace=False)) for e in k.episode], ddof=1))
        sdr = np.array(sdr)
        mig = k[k.migrates]
        nnw = int((mig.lfe_speed > 0).sum()); nse = int((mig.lfe_speed < 0).sum())
        binom_p = float(stats.binomtest(nnw, nnw + nse, 0.5).pvalue) if nnw + nse else np.nan
        # start vs previous end
        s = k.lfe_start10.to_numpy()[1:]; e = k.lfe_end10.to_numpy()[:-1]
        rho = stats.spearmanr(s, e)[0]
        perm = [stats.spearmanr(s, rng.permutation(e))[0] for _ in range(5000)]
        p_rho = float((np.sum(np.abs(perm) >= abs(rho)) + 1) / 5001)
        # tremor start spread (for comparison)
        sd_t = float(np.std(k.trem_start10.dropna(), ddof=1))
        out[reg] = dict(n_big=len(k), start10_sd_km=sd, uniform_span_km=[float(lo), float(hi)], uniform_sd_q05=float(np.quantile(sdu, 0.05)),
                        p_uniform=float((np.sum(sdu <= sd) + 1) / 10001), random10_sd_q05=float(np.quantile(sdr, 0.05)),
                        random10_sd_median=float(np.median(sdr)), p_random10=float((np.sum(sdr <= sd) + 1) / 2001),
                        clusters_rule4=bool(sd < np.quantile(sdu, 0.05)), clusters_beyond_random10=bool(sd < np.quantile(sdr, 0.05)),
                        n_migrating=int(k.migrates.sum()), n_NWward=nnw, n_SEward=nse, binom_p=binom_p,
                        start_vs_prev_end_rho=float(rho), start_vs_prev_end_p=p_rho,
                        tremor_start10_sd_km=sd_t,
                        n_agree_tremor_among_both_migrating=int(k[k.migrates & k.trem_migrates].agree_tremor.sum()),
                        n_both_migrating=int((k.migrates & k.trem_migrates).sum()),
                        n_thin_ok=int(k.thin_ok.sum()), n_onset_count_detection_driven=int(((k.onset_vs_count_rho < -0.3) & (k.onset_vs_count_p < 0.01)).sum()))
    json.dump(out, open(f"{OUT}/repeatability.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
