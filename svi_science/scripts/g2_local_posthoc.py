#!/usr/bin/env python
"""g2_local_posthoc (2026-10-11, POST-HOC, written after the pre-registered G2 tests) -- chance checks for the
individual cases the pre-registered tests flagged.

(1) Named-rule rate on random pseudo-earthquakes: for every ETS-status / between-status trial, 20 pseudo-trials at a
    random time on a random day of the same status and segment (ETS: any day inside an episode; between: the trial's
    random-window pool rule), the same H6 window, the same local-tile / random-window background and the same rule
    (a)-(d). Expected number of 'named' trials by chance = pseudo-trial rate x number of real trials.
(2) LFE burst starts after local earthquakes: fraction of between-episode bursts (tables/q2_bursts_<R>.csv, data)
    whose first LFE falls 0-6 h after a region-set GK mainshock M >= 1.5, vs the same for 10,000 sets of random
    burst-start times on between-episode observable days (same time of day).
(3) M >= 3.5 intraslab / near-interface earthquakes with epicentre <= 15 km of a family, between episodes: how many
    fall inside a burst window (first LFE -> last LFE + 24 h) vs expected = sum over events of the fraction of
    between-episode time covered by burst windows (binomial p).

    python svi_science/scripts/g2_local_posthoc.py -> tables/posthoc_*.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from scipy import stats
import g2_local_tests as LT
import g_trigger_tests as GT

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/g2_local/tables"
NS = 10 ** 9; H = 3600 * NS; DAY = 86400 * NS
rng = np.random.default_rng(7)


def named_rate(reg):
    K, W, days, dP, ep = LT.load_lfe(reg)
    TT = pd.read_parquet(f"{T}/q1_trials_{reg}.parquet")
    dns = days.dns.to_numpy(); grp = days.grp.to_numpy()
    out = []
    for status in ("ets", "between"):
        S = TT[TT.status == status]
        if status == "ets":
            pool_mask = days.in_ep.to_numpy()
        else:
            pool_mask = (~days.in_pad).to_numpy()
        hits = 0; n = 0
        for _, r in S.iterrows():
            cand = dns[pool_mask & (grp == r.grp)]
            if not len(cand):
                continue
            for d in rng.choice(cand, 20):
                t0 = int(d + rng.integers(0, DAY // NS) * NS)
                N = GT.count(K["data"], t0, t0 + 6 * H)
                if status == "ets":
                    L = 6 * H
                    ts = np.concatenate([t0 - L * np.arange(1, 5), t0 + L * np.arange(1, 5)])
                    bg = GT.count(K["data"], ts, ts + L); bc1 = GT.count(K["c1"], ts, ts + L)
                    bcd = GT.count(K["cand_c1"], ts, ts + L)
                else:
                    ok = cand[np.abs(cand - d) > 3 * DAY]
                    st = rng.choice(ok, 2000) + (t0 - (t0 // DAY) * DAY)
                    bg = GT.count(K["data"], st, st + 6 * H); bc1 = GT.count(K["c1"], st, st + 6 * H)
                    bcd = GT.count(K["cand_c1"], st, st + 6 * H)
                mz = GT.count(K["mz"], t0, t0 + 6 * H); c1 = GT.count(K["c1"], t0, t0 + 6 * H)
                cc1 = GT.count(K["cand_c1"], t0, t0 + 6 * H)
                mu = max(bg.mean(), 1 / 2000)
                pp = stats.poisson.sf(N - 1, mu)
                named = (N > bg.max()) and (pp < 0.01 / len(TT[TT.status.isin(["ets", "between"])])) and (mz >= 2) \
                    and (c1 <= np.quantile(bc1, 0.99)) and (cc1 <= np.quantile(bcd, 0.99))
                hits += int(named); n += 1
        real = int(S.named_abcd.sum())
        exp = hits / max(n, 1) * len(S)
        out.append(dict(region=reg, status=status, n_trials=len(S), named_real=real, n_pseudo=n, named_pseudo=hits,
                        expected_named=exp, p_poisson=stats.poisson.sf(real - 1, exp) if exp > 0 else np.nan))
        print(out[-1], flush=True)
    return out


def burst_after_eq(reg, A, Mc=1.5):
    K, W, days, dP, ep = LT.load_lfe(reg)
    B = pd.read_csv(f"{T}/q2_bursts_{reg}.csv"); B = B[B.series == "data"]
    t_start = pd.to_datetime(B.start_PT.str.replace(r" P[DS]T$", "", regex=True)).dt.tz_localize(
        "America/Los_Angeles", ambiguous="NaT", nonexistent="NaT").dt.tz_convert("UTC").dt.tz_localize(None)
    ts = t_start.dropna().astype("datetime64[ns]").astype("int64").to_numpy()
    X = A[A.in_analysis & A.gk_main & (A[f"{reg}_dist_km"] <= 150) & (A.mag >= Mc)]
    te = np.sort(X.t.astype("datetime64[ns]").astype("int64").to_numpy())

    def frac(ts):
        i = np.searchsorted(te, ts, "right")
        prev = np.where(i > 0, te[np.maximum(i - 1, 0)], -(2 ** 62))
        return np.mean(ts - prev <= 6 * H)
    obs = frac(ts)
    btw = days.dns[(~days.in_pad)].to_numpy()
    tod = ts - (ts // DAY) * DAY
    R = np.array([frac(rng.choice(btw, len(ts)) + tod) for _ in range(10000)])
    return dict(region=reg, n_bursts=len(ts), frac_with_eq_0_6h_before=obs, random_mean=R.mean(),
                p=(1 + (R >= obs - 1e-12).sum()) / 10001)


def big_in_burst(reg):
    K, W, days, dP, ep = LT.load_lfe(reg)
    B = pd.read_csv(f"{T}/q2_bursts_{reg}.csv"); B = B[B.series == "data"]
    TT = pd.read_parquet(f"{T}/q1_trials_{reg}.parquet")
    S = TT[(TT.status == "between") & TT.cls.isin(["intraslab", "near_interface"]) & (TT.mag >= 3.5) & (TT.epi_km <= 15)]
    tl = lambda s: pd.to_datetime(s.str.replace(r" P[DS]T$", "", regex=True)).dt.tz_localize(
        "America/Los_Angeles", ambiguous="NaT", nonexistent="NaT").dt.tz_convert("UTC").dt.tz_localize(None).astype("datetime64[ns]").astype("int64")
    bs = tl(B.start_PT).to_numpy(); be = tl(B.end_PT).to_numpy() + DAY
    btw_time = (~days.in_pad).sum() * DAY
    cover = np.sum(be - bs) / btw_time
    inside = [bool(((t >= bs) & (t < be)).any()) for t in S.t0]
    k = int(np.sum(inside))
    return dict(region=reg, n_events=len(S), n_inside_burst=k, burst_time_fraction=cover,
                p_binom=stats.binom.sf(k - 1, len(S), cover),
                events=";".join(f"{a} M{m:.1f}{' IN' if i else ''}" for a, m, i in zip(S.t_PT, S.mag, inside)))


def main():
    A = pd.read_parquet(f"{T}/eq_local.parquet"); A["t"] = A.t.astype("datetime64[ns]")
    if not os.path.exists(f"{T}/posthoc_named_rate.csv"):
        n = []
        for reg in ("SVI", "WA"):
            n += named_rate(reg)
        pd.DataFrame(n).to_csv(f"{T}/posthoc_named_rate.csv", index=False)
    b = [burst_after_eq(r, A) for r in ("SVI", "WA")]
    pd.DataFrame(b).to_csv(f"{T}/posthoc_burst_after_eq.csv", index=False); print(b)
    c = [big_in_burst(r) for r in ("SVI", "WA")]
    pd.DataFrame(c).to_csv(f"{T}/posthoc_big_in_burst.csv", index=False); print(c)


if __name__ == "__main__":
    main()
