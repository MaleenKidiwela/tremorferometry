"""Q6 tests on PNSN tremor: front arrest vs time since last activity (walk / stop hazard and enter / fail-to-enter),
episode-order shuffle referee; Strait of Juan de Fuca start / stop / crossing tests with uniform, matched and
own-footprint / rigid-shift controls.

Rules: svi_science/results/q6_front_arrest/PLAN.md. Reads q6_build.py outputs.
"""
import json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q6_front_arrest"
BIN = 5.0
import os
NDRAW = int(os.environ.get('Q6_NDRAW', 2000))
KEYS = 1e5
NORTH = (-200.0, 200.0)
END = 25.0


# ----------------------------------------------------------------------------- helpers
def logit_fit(X, y, iters=25):
    """Logistic regression by IRLS; X includes intercept column. Returns coefficients."""
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta))
        w = p * (1 - p) + 1e-9
        H = X.T @ (X * w[:, None]) + 1e-6 * np.eye(X.shape[1])
        step = np.linalg.solve(H, X.T @ (y - p))
        b += step
        if np.max(np.abs(step)) < 1e-8:
            break
    return b


def zs(v):
    sd = v.std()
    return (v - v.mean()) / (sd if sd > 0 else 1.0)


NROWS = {}


def beta_T(T, y, dist, lsize, prop, mask, name=None):
    m = mask & np.isfinite(T) & (T > 0)
    if name is not None:
        NROWS[name] = (int(m.sum()), int((y[m] == 0).sum()))
    if m.sum() < 20 or len(np.unique(y[m])) < 2:
        return np.nan
    X = np.column_stack([np.ones(m.sum()), zs(np.log10(T[m])), dist[m] / 100.0, zs(lsize[m]), zs(prop[m])])
    return float(logit_fit(X, y[m].astype(float))[1])


class History:
    """Activity pairs (component, bin) with offsets from the component start; query last activity before a time."""

    def __init__(self, be, E):
        a = be[be.active].merge(E[["comp", "t_first", "dur", "big"]], on="comp")
        self.pc = a.comp.values; self.pb = a.bin.values.astype(float)
        self.po = (a.t_med - a.t_first).values; self.pdur = a.dur.values; self.pbig = a.big.values
        self.slot = E.set_index("comp").t_first

    def times(self, slot_of_comp):
        s = slot_of_comp[self.pc]
        return s + self.po, s + self.pdur

    @staticmethod
    def last_before(pb, pt, qb, qt):
        key = np.sort(pb * KEYS + pt)
        q = qb * KEYS + qt
        i = np.searchsorted(key, q, side="left") - 1
        out = np.full(len(qb), np.nan)
        ok = i >= 0
        v = key[np.clip(i, 0, None)]
        ok &= np.floor(v / KEYS) == qb
        out[ok] = v[ok] - qb[ok] * KEYS
        return out

    @staticmethod
    def last_ended_before(pb, pt, pend, qb, qcut):
        o = np.lexsort((pend, pb))
        key = pb[o] * KEYS + pend[o]
        cm = np.maximum.accumulate(pb[o] * KEYS + pt[o])
        i = np.searchsorted(key, qb * KEYS + qcut, side="right") - 1
        out = np.full(len(qb), np.nan)
        ok = i >= 0
        v = cm[np.clip(i, 0, None)]
        ok &= v >= qb * KEYS
        out[ok] = v[ok] - qb[ok] * KEYS
        return out

    def T_all(self, slot_of_comp, qb, qcomp):
        qt = slot_of_comp[qcomp]
        pt, pend = self.times(slot_of_comp)
        last = self.last_before(self.pb, pt, qb, qt)
        lastbig = self.last_before(self.pb[self.pbig], pt[self.pbig], qb, qt)
        last7 = self.last_ended_before(self.pb, pt, pend, qb, qt - 7.0)
        return qt - last, qt - lastbig, qt - last7


def perm_slots(E, rng):
    """Return array slot_of_comp (indexed by comp id): big comps get permuted big slots, others other slots."""
    s = np.zeros(E.comp.max() + 1)
    for flag in (True, False):
        c = E.comp.values[E.big.values == flag]; t = E.t_first.values[E.big.values == flag]
        s[c] = t[rng.permutation(len(t))]
    return s


def pval(real, null, side):
    null = np.asarray(null)[np.isfinite(null)]
    if side == "+":
        return float((1 + (null >= real).sum()) / (1 + len(null)))
    if side == "-":
        return float((1 + (null <= real).sum()) / (1 + len(null)))
    return float((1 + (np.abs(null) >= abs(real)).sum()) / (1 + len(null)))


# ----------------------------------------------------------------------------- Q6.1
def build_rows(E, be, info, P_bin_counts):
    edges0 = info["edges0"]; nb = info["n_bins"]
    smin, smax = info["domain_s"]
    centre = edges0 + (np.arange(nb) + 0.5) * BIN
    endzone = (centre < smin + END) | (centre > smax - END)
    own = be.set_index(["comp", "bin"]).n
    W, S = [], []
    B = E[E.big & E.start_bin.notna()]
    act = {c: set(g.bin.values) for c, g in be[be.active].groupby("comp")}
    for _, r in B.iterrows():
        c = int(r.comp); sb = int(r.start_bin); lo = int(r.run_lo_bin); hi = int(r.run_hi_bin)
        base = dict(comp=c, lsize=np.log10(r.n), era=("E12" if r.t_first < (pd.Timestamp("2015-01-01") - pd.Timestamp("2010-01-01")).days else "E3"))
        for step, end in ((+1, hi), (-1, lo)):
            b = sb + step
            while (step > 0 and b <= end) or (step < 0 and b >= end):
                if not endzone[b]:
                    W.append(dict(base, bin=b, y=1, dirn=step))
                b += step
            if 0 <= b < nb and not endzone[b]:
                W.append(dict(base, bin=b, y=0, dirn=step))
        for b in range(max(0, sb - 30), min(nb, sb + 31)):
            if endzone[b]:
                continue
            S.append(dict(base, bin=b, y=int(b in act.get(c, ())), dirn=0))
    out = []
    for rows in (W, S):
        R = pd.DataFrame(rows)
        R = R.merge(E[["comp", "start_s"]], on="comp")
        R["s"] = centre[R.bin.values]
        R["dist"] = np.abs(R.s - R.start_s)
        ownn = np.array([own.get((c, b), 0) for c, b in zip(R.comp, R.bin)])
        R["prop"] = np.log10(P_bin_counts[R.bin.values] - ownn + 1.0)
        R["north"] = R.s.between(*NORTH)
        out.append(R)
    return out


def q61(E, be, info, Pcounts, rng):
    H = History(be, E)
    W, S = build_rows(E, be, info, Pcounts)
    real_slot = np.zeros(E.comp.max() + 1); real_slot[E.comp.values] = E.t_first.values
    t0_cat = 5.0  # catalogue starts 2010-01-06 PT (day ~5.7)
    specs = []   # (name, rows, Tkind, mask, censor_kept)
    for nm, R in (("walk", W), ("enter", S)):
        allm = np.ones(len(R), bool)
        specs += [(f"{nm}_Tany_full", R, 0, allm, False), (f"{nm}_Tbig_full", R, 1, allm, False),
                  (f"{nm}_Tany_north", R, 0, R.north.values, False), (f"{nm}_Tbig_north", R, 1, R.north.values, False)]
        if nm == "walk":
            specs += [("walk_Tany7_full", R, 2, allm, False), ("walk_Tany_full_censkept", R, 0, allm, True),
                      ("walk_Tany_E12", R, 0, (R.era == "E12").values, False), ("walk_Tany_E3", R, 0, (R.era == "E3").values, False),
                      ("walk_Tany7_north", R, 2, R.north.values, False)]

    def run(slot, record=False):
        res = {}
        cache = {}
        for nm, R, k, m, ck in specs:
            key = id(R)
            if key not in cache:
                cache[key] = H.T_all(slot, R.bin.values.astype(float), R.comp.values)
            T = cache[key][k].copy()
            if ck:
                T[~np.isfinite(T)] = slot[R.comp.values][~np.isfinite(T)] - t0_cat
            res[nm] = beta_T(T, R.y.values, R.dist.values, R.lsize.values, R.prop.values, m, nm if record else None)
        return res, cache

    real, cache = run(real_slot, record=True)
    W["T_any"], W["T_big"], W["T_any7"] = cache[id(W)]
    S["T_any"], S["T_big"], S["T_any7"] = cache[id(S)]
    W.to_csv(f"{OUT}/walk_rows_pnsn.csv", index=False); S.to_csv(f"{OUT}/enter_rows_pnsn.csv", index=False)
    nulls = {k: [] for k in real}
    # stop-fraction by T_any quintile (walk, full) — real and shuffle
    qs = np.nanquantile(np.log10(W.T_any[W.T_any > 0]), [0.2, 0.4, 0.6, 0.8])

    def stopfrac(T):
        lt = np.log10(np.where(T > 0, T, np.nan)); k = np.digitize(lt, qs); ok = np.isfinite(lt)
        return [float((W.y.values[ok & (k == j)] == 0).mean()) for j in range(5)]
    sf_real = stopfrac(W.T_any.values); sf_null = []
    for i in range(NDRAW):
        slot = perm_slots(E, rng)
        r, c = run(slot)
        for k in r:
            nulls[k].append(r[k])
        sf_null.append(stopfrac(c[id(W)][0]))
        if i % 100 == 0:
            print("draw", i, pd.Timestamp.now(), flush=True)
    tests = {}
    for k in real:
        nl = np.array(nulls[k])
        tests[k] = dict(beta=real[k], p_pos=pval(real[k], nl, "+"), p_neg=pval(real[k], nl, "-"), p_two=pval(real[k] - np.nanmean(nl), nl - np.nanmean(nl), "2"),
                        null_mean=float(np.nanmean(nl)), null_q05=float(np.nanquantile(nl, 0.05)), null_q95=float(np.nanquantile(nl, 0.95)),
                        n_rows=NROWS[k][0], n_y0=NROWS[k][1])
    sf_null = np.array(sf_null)
    sfq = dict(edges_log10T=qs.tolist(), real=sf_real, null_q05=np.quantile(sf_null, 0.05, 0).tolist(),
               null_q95=np.quantile(sf_null, 0.95, 0).tolist(), null_mean=sf_null.mean(0).tolist())
    # descriptive: stop trials vs continue trials median T
    desc = dict(n_big=int(E.big.sum()), walk_rows=int(len(W)), walk_stops=int((W.y == 0).sum()),
                walk_median_T_continue=float(np.nanmedian(W.T_any[W.y == 1])), walk_median_T_stop=float(np.nanmedian(W.T_any[W.y == 0])),
                walk_censored_frac=float((~np.isfinite(W.T_any)).mean()), enter_rows=int(len(S)))
    return tests, sfq, desc


# ----------------------------------------------------------------------------- Q6.2
def poisson_binomial_sf(p, k):
    """P(sum >= k) for independent Bernoulli(p_i)."""
    dist = np.zeros(len(p) + 1); dist[0] = 1.0
    for q in p:
        dist[1:] = dist[1:] * (1 - q) + dist[:-1] * q; dist[0] *= (1 - q)
    return float(dist[k:].sum())


def q62(E, ep, info, rng):
    smin, smax = info["domain_s"]
    lo_ok, hi_ok = smin + END, smax - END
    sall = ep.s_km.values
    res = {}
    B = E[E.big & E.start_bin.notna()].copy()
    A30 = E[(E.n >= 30) & E.start_bin.notna()].copy()
    comp_s = {c: g.s_km.values for c, g in ep[ep.comp.isin(A30.comp)].groupby("comp")}
    comp_t = {c: g.tday.values for c, g in ep[ep.comp.isin(A30.comp)].groupby("comp")}
    for lab, D in (("big", B), ("all30", A30)):
        for dom, (a, b) in (("north", NORTH), ("full", (lo_ok, hi_ok))):
            st = D[(D.start_s >= max(a, lo_ok)) & (D.start_s <= min(b, hi_ok))]
            dens = sall[(sall >= max(a, lo_ok)) & (sall <= min(b, hi_ok))]
            for hw in (15, 10, 25):
                k = int((np.abs(st.start_s) <= hw).sum()); n = len(st)
                pu = (2 * hw) / (min(b, hi_ok) - max(a, lo_ok))
                pm = float((np.abs(dens) <= hw).mean())
                un = rng.binomial(n, pu, 10000); mn = rng.binomial(n, pm, 10000)
                r = dict(n_starts=n, k_strait=k, frac=k / max(n, 1), uniform_expect=n * pu, matched_expect=n * pm,
                         p_uniform=float((1 + (un >= k).sum()) / 10001), p_matched=float((1 + (mn >= k).sum()) / 10001))
                if lab == "big" and hw == 15:
                    pe = []
                    for _, e in st.iterrows():
                        s = comp_s[e.comp]; m10 = max(1, int(np.ceil(0.1 * len(s))))
                        meds = np.array([np.median(rng.choice(s, m10, replace=False)) for _ in range(300)])
                        pe.append(float((np.abs(meds) <= hw).mean()))
                    r["own_expect"] = float(np.sum(pe)); r["p_own"] = poisson_binomial_sf(np.array(pe), k)
                    r["episodes_at_strait"] = st[np.abs(st.start_s) <= hw][["start_date_PT", "n", "start_s", "stop_se_s", "stop_nw_s"]].to_dict("records")
                res[f"starts_{lab}_{dom}_hw{hw}"] = r
            # start2d column
            st2 = D[(D.start2d_s >= max(a, lo_ok)) & (D.start2d_s <= min(b, hi_ok))]
            res[f"starts2d_{lab}_{dom}_hw15"] = dict(n=len(st2), k=int((np.abs(st2.start2d_s) <= 15).sum()))
            # stops
            stops = []
            for _, e in D[(D.start_s >= max(a, lo_ok)) & (D.start_s <= min(b, hi_ok))].iterrows():
                for s in (e.stop_se_s, e.stop_nw_s):
                    if lo_ok <= s <= hi_ok and a <= s <= b:
                        stops.append(s)
            stops = np.array(stops)
            k = int((np.abs(stops) <= 15).sum()); n = len(stops)
            pu = 30 / (min(b, hi_ok) - max(a, lo_ok)); pm = float((np.abs(dens) <= 15).mean())
            un = rng.binomial(n, pu, 10000); mn = rng.binomial(n, pm, 10000)
            # rigid shift: runs of these episodes moved by U(-100, 100) km, kept inside [lo_ok, hi_ok]
            runs = D[(D.start_s >= max(a, lo_ok)) & (D.start_s <= min(b, hi_ok))][["stop_se_s", "stop_nw_s"]].values
            sh = []
            for _ in range(10000 if lab == "big" else 2000):
                off = rng.uniform(-100, 100, len(runs))
                lo2 = runs[:, 0] + off; hi2 = runs[:, 1] + off
                ok = (lo2 >= lo_ok) & (hi2 <= hi_ok)
                off = np.where(ok, off, 0.0); lo2 = runs[:, 0] + off; hi2 = runs[:, 1] + off
                ends = np.r_[lo2, hi2]
                ends = ends[(ends >= a) & (ends <= b) & (ends >= lo_ok) & (ends <= hi_ok)]
                sh.append(int((np.abs(ends) <= 15).sum()))
            sh = np.array(sh)
            res[f"stops_{lab}_{dom}_hw15"] = dict(n_stops=n, k_strait=k, frac=k / max(n, 1), uniform_expect=n * pu, matched_expect=n * pm,
                                                  shift_mean=float(sh.mean()),
                                                  p_uniform=float((1 + (un >= k).sum()) / 10001), p_matched=float((1 + (mn >= k).sum()) / 10001),
                                                  p_shift=float((1 + (sh >= k).sum()) / (1 + len(sh))))
    # crossing (big, north)
    cz = []
    for z0 in np.arange(NORTH[0], NORTH[1] - 30 + 1, 5):
        z1 = z0 + 30; arr = cr = 0
        for _, e in B.iterrows():
            if e.start_s < z0 and e.stop_nw_s >= z0:
                arr += 1; cr += int(e.stop_nw_s >= z1)
            elif e.start_s > z1 and e.stop_se_s <= z1:
                arr += 1; cr += int(e.stop_se_s <= z0)
        cz.append(dict(z0=float(z0), z1=float(z1), arrived=arr, crossed=cr, rate=cr / arr if arr else np.nan))
    CZ = pd.DataFrame(cz); CZ.to_csv(f"{OUT}/crossing_zones.csv", index=False)
    st = CZ[CZ.z0 == -15].iloc[0]
    others = CZ[(CZ.z1 <= -15) | (CZ.z0 >= 15)]
    from scipy.stats import fisher_exact
    # non-overlapping comparison: zones tiling [-200,200] by 30 km excluding the Strait zone
    tiles = CZ[((CZ.z0 + 15) % 30 == 0) & (CZ.z0 != -15)]
    pa, pc_ = int(tiles.arrived.sum()), int(tiles.crossed.sum())
    orr, pf = fisher_exact([[int(st.crossed), int(st.arrived - st.crossed)], [pc_, pa - pc_]])
    res["crossing"] = dict(strait_arrived=int(st.arrived), strait_crossed=int(st.crossed), strait_rate=float(st.rate),
                           tiles_arrived=pa, tiles_crossed=pc_, tiles_rate=pc_ / pa,
                           strait_rank_low_among_all_zones=int((CZ.rate.dropna() < st.rate).sum()) + 1, n_zones=int(CZ.rate.notna().sum()),
                           fisher_odds=float(orr), fisher_p=float(pf),
                           tiles=tiles[["z0", "z1", "arrived", "crossed", "rate"]].to_dict("records"))
    # detectability diagnostic
    zden = []
    for z0 in np.arange(-195, 166, 30):
        zden.append(dict(z0=float(z0), z1=float(z0 + 30), n_epi=int(((sall >= z0) & (sall < z0 + 30)).sum())))
    res["density_tiles_north"] = zden
    return res


def main():
    rng = np.random.default_rng(0)
    info = json.load(open(f"{OUT}/build_info.json"))
    E = pd.read_csv(f"{OUT}/episodes_pnsn.csv")
    be = pd.read_parquet(f"{OUT}/bin_episode.parquet")
    ep = pd.read_parquet(f"{OUT}/epicentres.parquet", columns=["s_km", "bin", "tday", "comp"])
    Pcounts = np.bincount(ep.bin.values, minlength=info["n_bins"]).astype(float)
    out = {}
    out["q62"] = q62(E, ep, info, rng)
    json.dump(out, open(f"{OUT}/tests.json", "w"), indent=1, default=float)
    print(json.dumps(out["q62"], indent=1, default=float)[:6000])
    tests, sfq, desc = q61(E, be, info, Pcounts, rng)
    out["q61"] = tests; out["q61_stopfrac_quintiles"] = sfq; out["q61_desc"] = desc
    json.dump(out, open(f"{OUT}/tests.json", "w"), indent=1, default=float)
    rows = [dict(test=k, **v) for k, v in tests.items()]
    pd.DataFrame(rows).to_csv(f"{OUT}/tests_table.csv", index=False)
    print(pd.DataFrame(rows).to_string()); print(json.dumps(desc, indent=1)); print(json.dumps(sfq, indent=1))


if __name__ == "__main__":
    main()
