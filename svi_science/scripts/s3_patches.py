"""S4 (inside S3): patch strength per big episode and persistence vs catch-up tests. Pre-registration: RESULTS_S3.md §0 rule 7.

Units: SVI = the 12 S0 location groups; NW = 10-km along-strike bins (source positions; sources outside the A/B span
go to the end bins). Headline counts: counted once, P >= 0.75 (NW selected_v2).
Writes results/s3/patch_strength.csv, patch_tests.json.
"""
import json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
rng = np.random.default_rng(11)
NPERM = 2000


def headline(g, variant="once075"):
    d = g[g.series == "data"]
    sel = d.selected_v2 if "selected_v2" in d else True
    if variant == "once075":
        return d[(d.P_same >= 0.75) & sel & d.counted_once]
    if variant == "raw075":
        return d[(d.P_same >= 0.75) & sel]
    return d[d.counted_once]


def counts(reg, variant="once075"):
    EP = pd.read_csv(f"{OUT}/episodes_all.csv", parse_dates=["w0", "w1"])
    big = EP[(EP.region == reg) & EP.big].sort_values("w0")
    E = pd.read_parquet(f"{OUT}/events_{reg.lower()}.parquet")
    h = headline(E[E.ep.isin(big.episode)], variant)
    if reg == "SVI":
        h = h.assign(patch=h.group.astype(int))
    else:
        src = pd.read_csv(f"{OUT}/sources_nw.csv"); ab = src[src.quality.isin(["A", "B"])]
        lo = np.floor(np.quantile(ab.as_km, 0.03) / 10) * 10; hi = np.ceil(np.quantile(ab.as_km, 0.97) / 10) * 10
        edges = np.arange(lo, hi + 10, 10)
        b = np.clip(np.digitize(h.as_km, edges) - 1, 0, len(edges) - 2)
        h = h.assign(patch=[f"{edges[i]:+.0f}..{edges[i + 1]:+.0f}" for i in b])
    C = h.pivot_table(index="patch", columns="ep", values="time", aggfunc="size", fill_value=0)
    C = C[big.episode.tolist()]
    if reg == "NW":
        C = C.loc[sorted(C.index, key=lambda s: float(s.split("..")[0]))]
    mid = (big.w0 + (big.w1 - big.w0) / 2).to_numpy()
    return C, mid


def resid(C, i, j):
    """Pearson residuals of episodes i and j against the pooled share of all OTHER episodes."""
    others = [k for k in range(C.shape[1]) if k not in (i, j)]
    p = C[:, others].sum(1) / C[:, others].sum()
    p = np.clip(p, 1e-4, None); p = p / p.sum()
    r = []
    for k in (i, j):
        Ek = C[:, k].sum() * p
        r.append((C[:, k] - Ek) / np.sqrt(Ek))
    return r


def pair_stats(C, mid, min_gap=90):
    n = C.shape[1]; allc, cons, seps = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            gap = abs((mid[j] - mid[i]) / np.timedelta64(1, "D"))
            if gap < min_gap:
                continue
            ri, rj = resid(C, i, j)
            c = np.corrcoef(ri, rj)[0, 1]
            allc.append(c); seps.append(gap)
            if j == i + 1:
                cons.append(c)
    return np.array(allc), np.array(cons), np.array(seps)


def excursion(C):
    """sum over groups of squared cumulative Pearson residuals (leave-one-out shares) in time order."""
    n = C.shape[1]; R = np.zeros(C.shape)
    for k in range(n):
        others = [m for m in range(n) if m != k]
        p = C[:, others].sum(1) / C[:, others].sum(); p = np.clip(p, 1e-4, None); p /= p.sum()
        Ek = C[:, k].sum() * p; R[:, k] = (C[:, k] - Ek) / np.sqrt(Ek)
    return float((np.cumsum(R, axis=1) ** 2).sum()), R


def multinomial_ctrl(C, nrep=500):
    p = C.sum(1) / C.sum()
    out = []
    for _ in range(nrep):
        S = np.stack([rng.multinomial(int(C[:, k].sum()), p) for k in range(C.shape[1])], 1)
        out.append(S)
    return out


def main():
    rows, tests = [], {}
    for reg in ["SVI", "NW"]:
        for variant in ["once075", "raw075"]:
            Cdf, mid = counts(reg, variant)
            C = Cdf.to_numpy().astype(float)
            allc, cons, seps = pair_stats(C, mid)
            exc, R = excursion(C)
            # control 1 (pre-registered): patch labels shuffled within each episode
            pm = []
            for _ in range(NPERM):
                Cp = np.stack([rng.permutation(C[:, k]) for k in range(C.shape[1])], 1)
                pm.append(pair_stats(Cp, mid)[0].mean())
            pm = np.array(pm)
            # control 2: independent multinomial episodes with the pooled shares (no memory, same sizes)
            mm, mc = [], []
            for S in multinomial_ctrl(C, 500):
                a, c, _ = pair_stats(S.astype(float), mid); mm.append(a.mean()); mc.append(c.mean() if len(c) else np.nan)
            mm, mc = np.array(mm), np.array(mc)
            # control 3: episode order shuffled (consecutive-pair mean and cumulative excursion)
            oc, oe = [], []
            for _ in range(NPERM):
                o = rng.permutation(C.shape[1])
                _, c, _ = pair_stats(C[:, o], mid)
                oc.append(c.mean() if len(c) else np.nan)
                oe.append(excursion(C[:, o])[0])
            oc, oe = np.array(oc), np.array(oe)
            # raw-share rank agreement (descriptive: are the same patches always the big ones?)
            sh = C / C.sum(0)
            rk = [stats.spearmanr(sh[:, i], sh[:, j])[0] for i in range(C.shape[1]) for j in range(i + 1, C.shape[1])]
            # correlation vs time separation
            rs = stats.spearmanr(seps, allc)
            t = dict(n_patches=int(C.shape[0]), n_episodes=int(C.shape[1]), n_pairs_ge90d=int(len(allc)),
                     mean_resid_corr=float(allc.mean()),
                     labelshuffle_mean=float(pm.mean()), labelshuffle_q99=float(np.quantile(pm, 0.99)),
                     p_persist_labelshuffle=float((np.sum(pm >= allc.mean()) + 1) / (NPERM + 1)),
                     multinomial_mean=float(mm.mean()), multinomial_q01=float(np.quantile(mm, 0.01)), multinomial_q99=float(np.quantile(mm, 0.99)),
                     p_persist_multinomial=float((np.sum(mm >= allc.mean()) + 1) / (len(mm) + 1)),
                     p_lower_multinomial=float((np.sum(mm <= allc.mean()) + 1) / (len(mm) + 1)),
                     n_consecutive_pairs=int(len(cons)), mean_consecutive_corr=float(cons.mean()) if len(cons) else None,
                     ordershuffle_consecutive_mean=float(np.nanmean(oc)),
                     p_catchup_ordershuffle=float((np.sum(oc <= cons.mean()) + 1) / (np.isfinite(oc).sum() + 1)) if len(cons) else None,
                     multinomial_consecutive_mean=float(np.nanmean(mc)),
                     excursion=exc, ordershuffle_excursion_q05=float(np.quantile(oe, 0.05)), ordershuffle_excursion_q95=float(np.quantile(oe, 0.95)),
                     p_excursion_low=float((np.sum(oe <= exc) + 1) / (NPERM + 1)), p_excursion_high=float((np.sum(oe >= exc) + 1) / (NPERM + 1)),
                     corr_vs_separation_rho=float(rs[0]), corr_vs_separation_p=float(rs[1]),
                     raw_share_rank_agreement_mean_spearman=float(np.mean(rk)),
                     chi2_overdispersion_per_df=float((R ** 2).sum() / ((C.shape[0] - 1) * (C.shape[1] - 1))))
            t["persistent"] = bool(t["mean_resid_corr"] > 0 and t["p_persist_labelshuffle"] < 0.01)
            t["catch_up"] = bool(len(cons) and t["mean_consecutive_corr"] < 0 and t["p_catchup_ordershuffle"] < 0.01 and exc < t["ordershuffle_excursion_q05"])
            tests[f"{reg}|{variant}"] = t
            print(reg, variant, json.dumps(t, indent=0))
            if variant == "once075":
                p_all = C.sum(1) / C.sum()
                for gi, g in enumerate(Cdf.index):
                    for k, e in enumerate(Cdf.columns):
                        Ek = C[:, k].sum() * p_all[gi]
                        rows.append(dict(region=reg, patch=g, episode=e, n=int(C[gi, k]), episode_total=int(C[:, k].sum()),
                                         share=C[gi, k] / C[:, k].sum(), longterm_share=p_all[gi],
                                         ratio_to_longterm=C[gi, k] / Ek if Ek > 0 else np.nan,
                                         pearson_resid_loo=R[gi, k]))
    pd.DataFrame(rows).to_csv(f"{OUT}/patch_strength.csv", index=False)
    json.dump(tests, open(f"{OUT}/patch_tests.json", "w"), indent=1)


if __name__ == "__main__":
    main()
