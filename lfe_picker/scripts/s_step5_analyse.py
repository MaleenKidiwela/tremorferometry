#!/usr/bin/env python
"""P3D STEP 5 analysis -- matched-FA tiers per network epoch, marginal detections, the a3b
withheld-site certification, and the certified recall gain against the Step-0 prediction.

Inputs: results/p3d/trackS/detections_<YYYY-MM>_<cfg>.parquet for cfg in
        old13 (frozen bank, withheld sites out), dense (frozen + transfer bank, withheld out),
        denseprod (frozen + transfer bank, every site).

Matched FA (Merlin Q2(i)): the two networks are compared at the SAME network-wide quiet-day false
events/day, not at the same MAD.  The reference rate is old13's at its own epoch tier; each other
configuration's threshold is the quantile of ITS quiet-month detections that reproduces that rate.

Certification (a3b): a marginal detection set is certified at a WITHHELD site by the stack excess
      excess_dB = 10 log10(1 + p^2 N (a/sigma)^2)
with (a/sigma) calibrated on the SAME family's own hp events at the SAME site (p = the mean p_real
of that control set).  Solving the same expression at the marginal set's N gives its purity p.

certified gain = (sum_f p_f N_marginal,f - sum_f p'_f N_reverse,f) / sum_f N_old,f
"""
from __future__ import annotations
import os, sys, json, glob, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_transfer as T                                                     # noqa: E402

TRK = f"{S.P3D}/trackS"
CFGS = ["old13", "dense", "denseprod"]
YEARS = list(range(2019, 2027))


def load_cfg(cfg, months=None):
    fs = sorted(glob.glob(f"{TRK}/detections_*_{cfg}.parquet"))
    D = []
    for f in fs:
        ym = os.path.basename(f).split("_")[1]
        if months and ym not in months:
            continue
        d = pd.read_parquet(f)
        d["ym"] = ym
        D.append(d)
    if not D:
        return pd.DataFrame()
    d = pd.concat(D, ignore_index=True)
    d["dt"] = pd.to_datetime(d.date) if "date" in d else pd.to_datetime(d.t, unit="s")
    d["episode"] = S.episode_of(d.dt)
    return d


def quiet_rate(d, quiet_months, thr):
    q = d[d.ym.isin(quiet_months)]
    if not len(q):
        return np.nan, 0
    nd = q.date.nunique() if "date" in q else q.dt.dt.date.nunique()
    return float((q.net_mad >= thr).sum()) / max(nd, 1), int(nd)


def thr_for_rate(d, quiet_months, target_rate):
    q = d[d.ym.isin(quiet_months)]
    if not len(q) or not np.isfinite(target_rate):
        return np.nan
    nd = q.date.nunique()
    v = np.sort(q.net_mad.to_numpy())[::-1]
    k = int(round(target_rate * nd))
    if k <= 0:
        return float(v[0]) + 1e-6
    if k >= len(v):
        return float(v[-1])
    return float(v[k - 1])


def match_one_to_one(ta, tb, tol=2.0):
    """greedy nearest one-to-one match of two sorted time arrays; returns boolean masks of matched."""
    ia = np.argsort(ta)
    ib = np.argsort(tb)
    a, b = ta[ia], tb[ib]
    ma = np.zeros(len(a), bool)
    mb = np.zeros(len(b), bool)
    j = 0
    for i in range(len(a)):
        while j < len(b) and b[j] < a[i] - tol:
            j += 1
        k = j
        best, bk = tol + 1, -1
        while k < len(b) and b[k] <= a[i] + tol:
            if not mb[k] and abs(b[k] - a[i]) < best:
                best, bk = abs(b[k] - a[i]), k
            k += 1
        if bk >= 0:
            ma[i] = True
            mb[bk] = True
    out_a = np.zeros(len(ta), bool)
    out_b = np.zeros(len(tb), bool)
    out_a[ia] = ma
    out_b[ib] = mb
    return out_a, out_b


def a3b_purity(dc, sta, bkey, t_marg, t_ctrl, p_ctrl, ttS, ST, rng):
    """a3b stack test at a WITHHELD site.  excess = 10 log10(1 + p^2 N (a/sigma)^2): calibrate
    (a/sigma) on the family's own hp control at the same site, then solve for p on the marginal set."""
    out = dict(sta=sta, bkey=int(bkey), n_ctrl=len(t_ctrl), n_marg=len(t_marg))
    if len(t_ctrl) < 15 or len(t_marg) < 8:
        out["status"] = "too_few"
        return out
    Lc = T.leg(dc, sta, bkey, t_ctrl, ttS, ST, rng, want_nulls=False, keep_stack=False)
    Lm = T.leg(dc, sta, bkey, t_marg, ttS, ST, rng, want_nulls=False, keep_stack=False)
    if Lc.get("status") == "no_data" or Lm.get("status") == "no_data":
        out["status"] = "no_data"
        return out
    ec, em = Lc["excess_db"], Lm["excess_db"]
    Nc, Nm = Lc["n_stacked"], Lm["n_stacked"]
    out.update(excess_ctrl_db=ec, excess_marg_db=em, N_ctrl=Nc, N_marg=Nm, p_ctrl=p_ctrl)
    if not (np.isfinite(ec) and np.isfinite(em)) or ec <= 0:
        out["status"] = "no_control"
        return out
    a2 = (10 ** (ec / 10.0) - 1.0) / max(p_ctrl ** 2 * Nc, 1e-9)     # (a/sigma)^2
    if a2 <= 0:
        out["status"] = "no_control"
        return out
    p2 = (10 ** (em / 10.0) - 1.0) / max(a2 * Nm, 1e-9)
    out["a_over_sigma"] = float(np.sqrt(a2))
    out["p_marg"] = float(np.clip(np.sqrt(max(p2, 0.0)), 0.0, 1.0))
    out["status"] = "ok"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet-months", default="")
    ap.add_argument("--test-episodes", default=",".join(S.TEST_EPISODES))
    ap.add_argument("--withheld-stations", default="B926,CLRS")
    ap.add_argument("--tol", type=float, default=2.0)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    ST = S.station_table()
    PR = pd.read_csv(f"{S.P3D}/premise_family_sites.csv").set_index("bkey")
    PRED = json.load(open(f"{S.P3D}/freeze.json"))["prediction_table"]
    eps = a.test_episodes.split(",")
    D = {c: load_cfg(c) for c in CFGS}
    have = {c: sorted(D[c].ym.unique()) for c in CFGS if len(D[c])}
    print(f"[step5] months per configuration: { {c: len(v) for c, v in have.items()} }")
    qm = [m for m in (a.quiet_months.split(",") if a.quiet_months else S.QUIET_MONTHS)
          if all(m in have.get(c, []) for c in CFGS)]
    print(f"[step5] quiet months usable for the FA-only rule: {qm}")

    # ---------------- matched-FA thresholds per network epoch
    rows = []
    for eid, s0, s1 in S.NET_EPOCHS:
        sel = {c: D[c][(D[c].dt >= s0) & (D[c].dt <= s1)] for c in CFGS if len(D[c])}
        if "old13" not in sel or not len(sel["old13"]):
            continue
        fz = json.load(open(f"{S.EXT}/results/frozen_operating_point.json"))
        Tc0 = float(fz["T_core"])
        tb = pd.read_csv(f"{S.EXT}/results/trackC/tiers_by_epoch.csv") \
            if os.path.exists(f"{S.EXT}/results/trackC/tiers_by_epoch.csv") else None
        if tb is not None:
            m = tb[tb.epoch == ("E2" if eid == "D1" else "E3")]
            if len(m):
                Tc0 = float(m.T_core.iloc[0])
        r0, nd = quiet_rate(sel["old13"], qm, Tc0)
        row = dict(epoch=eid, T_old13=Tc0, quiet_rate_per_day=r0, n_quiet_days=nd)
        for c in CFGS:
            if c == "old13" or c not in sel:
                continue
            row[f"T_{c}"] = thr_for_rate(sel[c], qm, r0)
        rows.append(row)
    TI = pd.DataFrame(rows)
    TI.to_csv(f"{S.P3D}/step5_tiers.csv", index=False)
    print(TI.to_string(index=False))

    # ---------------- marginal / reverse-marginal on the TEST episodes
    thr = {}
    for c in CFGS:
        col = "T_old13" if c == "old13" else f"T_{c}"
        thr[c] = {r.epoch: getattr(r, col, np.nan) for r in TI.itertuples()}

    def ep_of(ts):
        for eid, s0, s1 in S.NET_EPOCHS:
            if pd.Timestamp(s0) <= ts <= pd.Timestamp(s1):
                return eid
        return "D?"

    for c in CFGS:
        if len(D[c]):
            D[c] = D[c][D[c].episode.isin(eps)].copy()
            D[c]["netepoch"] = [ep_of(t) for t in D[c].dt]
            D[c]["keep"] = [r.net_mad >= thr[c].get(r.netepoch, np.inf) for r in D[c].itertuples()]
    A = D["old13"][D["old13"].keep]
    B = D["dense"][D["dense"].keep]
    fam_rows, marg, rev = [], {}, {}
    for f in sorted(set(A.fam) | set(B.fam)):
        ta = A[A.fam == f].t.to_numpy()
        tb = B[B.fam == f].t.to_numpy()
        ma, mb = match_one_to_one(ta, tb, a.tol)
        marg[int(f)] = tb[~mb]
        rev[int(f)] = ta[~ma]
        fam_rows.append(dict(bkey=int(f), n_old=len(ta), n_dense=len(tb), n_matched=int(ma.sum()),
                             n_marginal=int((~mb).sum()), n_reverse=int((~ma).sum()),
                             d45=int(PR.d45.get(int(f), 0)), d55=int(PR.d55.get(int(f), 0))))
    FM = pd.DataFrame(fam_rows)
    FM.to_csv(f"{S.P3D}/step5_family_counts.csv", index=False)

    # ---------------- a3b certification at the withheld site(s)
    ws = [s for s in a.withheld_stations.split(",") if s]
    idx = S.cache_index(years=YEARS, stas=ws)
    dc = S.DayCache(maxn=600, years=YEARS, stas=ws)
    tts = S.fam_tt(ws, ST)
    H = S.hp_events(0.8, episodes=eps)
    rng = np.random.default_rng(T.RNG_SEED)
    cert = []
    for f in sorted(marg):
        h = H[H.bkey == f]
        if not len(h):
            continue
        pc = float(h.p_real.mean())
        for w in ws:
            for kind, tt in (("marginal", marg[f]), ("reverse", rev[f])):
                if len(tt) < 8:
                    continue
                r = a3b_purity(dc, w, f, np.sort(tt), h.t.to_numpy(), pc, tts[(f, w)], ST, rng)
                r.update(kind=kind, bkey=int(f))
                cert.append(r)
        print(f"  [a3b] family {f}: marginal {len(marg[f])} reverse {len(rev[f])}", flush=True)
    CE = pd.DataFrame(cert)
    CE.to_csv(f"{S.P3D}/step5_a3b_certification.csv", index=False)

    # ---------------- the certified gain
    P = {}
    if len(CE):
        ok = CE[CE.status == "ok"]
        for (f, k), g in ok.groupby(["bkey", "kind"]):
            P[(int(f), k)] = float(g.p_marg.median())
    num = den = 0.0
    grows = []
    for r in FM.itertuples():
        pm = P.get((r.bkey, "marginal"), np.nan)
        pr = P.get((r.bkey, "reverse"), np.nan)
        g = (np.nan if not np.isfinite(pm) else
             (pm * r.n_marginal - (0.0 if not np.isfinite(pr) else pr * r.n_reverse)) / max(r.n_old, 1))
        grows.append(dict(bkey=r.bkey, d45=r.d45, d55=r.d55, n_old=r.n_old, n_dense=r.n_dense,
                          n_marginal=r.n_marginal, n_reverse=r.n_reverse,
                          p_marginal=pm, p_reverse=pr, certified_gain=g,
                          raw_gain=(r.n_dense - r.n_old) / max(r.n_old, 1)))
        if np.isfinite(g) and r.d45 >= 1:
            num += pm * r.n_marginal - (0.0 if not np.isfinite(pr) else pr * r.n_reverse)
            den += r.n_old
    G = pd.DataFrame(grows)
    G.to_csv(f"{S.P3D}/step5_gain_by_family.csv", index=False)
    pooled = num / den if den else np.nan
    sel = G[(G.d45 >= 1) & G.certified_gain.notna()]
    boot = []
    if len(sel) >= 5:
        rr = np.random.default_rng(11)
        for _ in range(2000):
            k = rr.integers(0, len(sel), len(sel))
            s = sel.iloc[k]
            boot.append(float((s.p_marginal * s.n_marginal
                               - s.p_reverse.fillna(0) * s.n_reverse).sum() / max(s.n_old.sum(), 1)))
    lb = float(np.percentile(boot, 2.5)) if boot else np.nan
    from scipy.stats import spearmanr
    sp = spearmanr(G.d45[G.certified_gain.notna()], G.certified_gain.dropna())
    raw_pooled = float((G[G.d45 >= 1].n_dense.sum() - G[G.d45 >= 1].n_old.sum())
                       / max(G[G.d45 >= 1].n_old.sum(), 1))
    pred = PRED.get("PREDICTED_RECALL_GAIN_d45ge1_by_epoch", {})
    out = dict(pooled_certified_gain=pooled, bootstrap_lb_2p5=lb,
               spearman_gain_vs_d45=float(sp.statistic), spearman_p=float(sp.pvalue),
               raw_pooled_gain_d45ge1=raw_pooled,
               n_families_scored=int(len(sel)), predicted_gain_by_epoch=pred,
               bar="certified gain >= +25 % pooled, bootstrap LB > 0, Spearman(gain, d45) > 0.3",
               verdict=("MET" if (np.isfinite(pooled) and pooled >= 0.25 and lb > 0
                                  and sp.statistic > 0.3) else "NOT MET"),
               quiet_months_used=qm, tiers=TI.to_dict("records"))
    json.dump(out, open(f"{S.P3D}/step5_verdict.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
