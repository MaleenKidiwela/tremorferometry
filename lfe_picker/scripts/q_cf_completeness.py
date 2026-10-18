#!/usr/bin/env python
"""P3 STEP 2 (CF kill test B) — COMPLETENESS: does the learned CF know which template detections are
real?

Merlin's residual use for the P2 characteristic function: a per-detection REALNESS feature that would
let the frozen template detector run below its own operating point without drowning in false alarms
(that is the only route from 12-18 % completeness towards 30-40 %).  The kill test:

  * candidates = every stored 2012 ETS-day detection of the frozen template detector
    (../lfe_extension/results/e5r_detections_a1_mergeB_L12.parquet, L = 12, net_mad >= 7);
  * feature   = the site-merged learned CF (log-odds, LOWHALF-normalised per station-day, averaged
    within the PGC and LRIV clusters, summed over sites) read at each station's PREDICTED S for that
    detection's own family;
  * truth     = the detection's own-family one-to-one Bostock match (the `match_own` column, Track D's
    +-2 s greedy matcher);
  * readout   = AUC inside net_mad bins 7-8, 8-9 and >= 9.

This is an UPPER BOUND: the CF was trained on 2012, so anything it can do here it can only do worse
on an unseen year.  Rule (pre-registered by Merlin): AUC < 0.65 in the 8-9 bin -> the completeness
thread is CLOSED; >= 0.75 -> recommend a retrain on 2008/09/12 only and a Track E feature; in
between -> report and decide.

    python scripts/q_cf_completeness.py --cf results/p2/cf_p2nnpu_2012 --year 2012
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
FS = P.FS


def logodds(p, lo=1e-4):
    p = np.clip(p.astype(np.float32), lo, 1 - lo)
    return np.log(p / (1 - p))


def lowhalf_norm(v, valid):
    """z with a NOISE-SIDE scale: median(|x-md| | x < md) == md - q25, 0.6745 sigma for a Gaussian.
    Same estimator as p2_backproj.BackProjector(norm='lowhalf'), applied per station-day."""
    x = v[valid]
    if x.size < 1000:
        return None
    md = float(np.median(x))
    q25 = float(np.quantile(x, 0.25))
    sd = (md - q25) * 1.4826
    if not np.isfinite(sd) or sd <= 0:
        return None
    return (v - md) / sd


def auc(score, y):
    score = np.asarray(score, float)
    y = np.asarray(y, bool)
    m = np.isfinite(score)
    score, y = score[m], y[m]
    if y.sum() == 0 or (~y).sum() == 0:
        return np.nan, int(y.sum()), int((~y).sum())
    r = pd.Series(score).rank().to_numpy()
    n1, n0 = y.sum(), (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), int(n1), int(n0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", default=f"{P.RES}/p2/cf_p2nnpu_2012")
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--min-mad", type=float, default=7.0)
    ap.add_argument("--out", default=f"{P3}/cf_completeness")
    a = ap.parse_args()

    D = pd.read_parquet(f"{P.EXT}/results/e5r_detections_a1_mergeB_L12.parquet")
    D = D[(D.L == 12.0) & (D.net_mad >= a.min_mad)].copy()
    D["t"] = D.t.astype(float)
    ts = pd.to_datetime(D.t, unit="s")
    D = D[ts.dt.year == a.year].copy()
    D["doy"] = ts[ts.dt.year == a.year].dt.dayofyear.to_numpy()
    ed = json.load(open(f"{P.RES}/p2/score_p2nnpu_{a.year}.json"))["ets_days"]
    D = D[D.doy.isin(ed)].reset_index(drop=True)
    print(f"[scope] {len(D)} template detections, {D.doy.nunique()} ETS days, "
          f"own-family match rate {D.match_own.mean():.3f}")

    ftt = P.fam_tt(P.S15)
    store = P.DayStore()
    sites = sorted({P.SITE_OF[s] for s in P.S15})
    net = np.zeros(len(D)); nsite = np.zeros(len(D), int)
    for doy, g in D.groupby("doy"):
        acc = np.zeros((len(sites), len(g))); cnt = np.zeros((len(sites), len(g)))
        for s in P.S15:
            f = f"{a.cf}/{s}_{a.year}_{doy}.npy"
            if not os.path.exists(f) or not store.has(s, a.year, int(doy)):
                continue
            i = store.idx(s, a.year, int(doy))
            cf = np.load(f).astype(np.float32)
            valid = cf > 0
            z = lowhalf_norm(logodds(cf), valid)
            if z is None:
                continue
            z = np.where(valid, z, np.nan)
            t0 = float(store.t0[i])
            j = np.rint((g.t.to_numpy() + np.array([ftt[(int(b), s)] for b in g.fam])
                         - t0) * FS).astype(np.int64)
            ok = (j >= 0) & (j < len(z))
            v = np.full(len(g), np.nan)
            v[ok] = z[j[ok]]
            k = sites.index(P.SITE_OF[s])
            good = np.isfinite(v)
            acc[k, good] += v[good]; cnt[k, good] += 1
        with np.errstate(invalid="ignore"):
            sm = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
        net[g.index.to_numpy()] = np.nansum(sm, 0)
        nsite[g.index.to_numpy()] = (cnt > 0).sum(0)
        print(f"  [{a.year}.{doy}] {len(g)} detections, sites {np.median((cnt>0).sum(0)):.0f}",
              flush=True)
    D["cf_net"] = net
    D["cf_nsite"] = nsite
    D = D[D.cf_nsite >= 5].reset_index(drop=True)
    D["madbin"] = pd.cut(D.net_mad, [7, 8, 9, 99], labels=["7-8", "8-9", ">=9"])
    rows = []
    for b, g in D.groupby("madbin", observed=True):
        A, n1, n0 = auc(g.cf_net, g.match_own.astype(bool))
        rows.append(dict(net_mad_bin=str(b), n=len(g), n_match=n1, n_nonmatch=n0, auc=A,
                         cf_net_median_match=float(np.nanmedian(g.cf_net[g.match_own.astype(bool)])),
                         cf_net_median_non=float(np.nanmedian(g.cf_net[~g.match_own.astype(bool)]))))
    A, n1, n0 = auc(D.cf_net, D.match_own.astype(bool))
    rows.append(dict(net_mad_bin="all", n=len(D), n_match=n1, n_nonmatch=n0, auc=A,
                     cf_net_median_match=float(np.nanmedian(D.cf_net[D.match_own.astype(bool)])),
                     cf_net_median_non=float(np.nanmedian(D.cf_net[~D.match_own.astype(bool)]))))
    R = pd.DataFrame(rows)
    R.to_csv(a.out + ".csv", index=False)
    D[["doy", "fam", "t", "net_mad", "match_own", "cf_net", "cf_nsite"]].to_parquet(
        a.out + "_detections.parquet", index=False)
    a89 = float(R[R.net_mad_bin == "8-9"].auc.iloc[0]) if (R.net_mad_bin == "8-9").any() else np.nan
    verdict = ("completeness thread CLOSED" if a89 < 0.65 else
               "recommend the retrain (do NOT run it here)" if a89 >= 0.75 else
               "in between - report")
    json.dump(dict(auc_8_9=a89, verdict=verdict, rows=R.to_dict("records"),
                   note="upper bound: the CF was trained on 2012"),
              open(a.out + ".json", "w"), indent=1, default=str)
    print("\n" + R.to_string(index=False))
    print(f"\n[verdict] AUC(8-9) = {a89:.3f} -> {verdict}")


if __name__ == "__main__":
    main()
