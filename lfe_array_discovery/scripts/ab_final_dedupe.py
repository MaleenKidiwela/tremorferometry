"""ONE final duplicate removal after a no-merge growth run (Merlin 09-28 step 2, user-approved).

Pair (a, b): for every detection of a, the b-detections within +-15 s give lags (0.5 s bins); modal lag = the
most populated bin; shared = detections with a b-partner in the modal bin +-1 bin (+-0.75 s); fraction =
shared / min(n_a, n_b).  Control: the same with every b time shifted by +60 s (lags measured in 45..75 s).
Duplicate when frac_real - frac_shift > 0.10.  Grouping: greedy star in descending (site-count score, n_det);
a family joins the representative with the largest excess > 0.10, else becomes a representative; so the kept
family of any duplicate set is the higher site-count score (tie: more detections).  Writes
  <out>/dedupe_pairs.csv   every pair with excess > 0.05 (a, b, n_a, n_b, lag_s, f_real, f_shift, excess)
  <out>/dedupe_families.csv grp, n_det, score, is_rep, rep, partners (all duplicate partners, ';'-joined)
  <out>/det_dedup.parquet   detections of the representatives only (unchanged rows)
usage: python ab_final_dedupe.py --det det_itN.parquet --score score.csv --score-col n_sites_3km --out DIR
"""
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_grow as G

WIN, SHIFT, BIN, NB, THR = 15.0, 60.0, 0.5, 1, 0.10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True); ap.add_argument("--score", required=True)
    ap.add_argument("--score-col", default="n_sites_3km"); ap.add_argument("--out", required=True)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    det = pd.read_parquet(a.det)
    sc = pd.read_csv(a.score).set_index("grp")[a.score_col]
    gids = np.sort(det.grp.unique()).astype(np.int64); Gn = gids.size
    gi = np.searchsorted(gids, det.grp.to_numpy(np.int64))
    n = np.bincount(gi, minlength=Gn).astype(float)
    pairs, _ = G._numba()
    nbin = int(round(WIN / BIN))
    acc = {}                                   # (ga, gb) -> [hist_real(2nbin+1), hist_shift]
    H = []
    for date, idx in det.groupby("date").indices.items():
        t = det.t.to_numpy(np.float64)[idx]; g = gi[idx]
        o = np.argsort(t, kind="stable"); t, g = t[o], g[o]
        ia, gb, d = pairs(t, g, SHIFT + WIN + 1e-6)
        if ia.size == 0:
            continue
        ga = g[ia]
        for lab, dd in (("r", d), ("s", d - SHIFT)):      # d = t_b - t_a ; shift: b at t_b - 60 => lag d - 60
            m = np.abs(dd) <= WIN
            k = np.round(dd[m] / BIN).astype(np.int64) + nbin
            H.append(pd.DataFrame(dict(a=ga[m], b=gb[m], k=k, c=lab, n=1)))
        if sum(len(h) for h in H) > 2e7:
            H = [pd.concat(H).groupby(["a", "b", "c", "k"], as_index=False).n.sum()]
    if not H:
        print("no pairs"); return
    Hc = pd.concat(H)
    Hc = Hc.groupby(["a", "b", "c", "k"], as_index=False).n.sum()
    rows = []
    for c in ("r", "s"):
        h = Hc[Hc.c == c]
        # modal bin per (a, b) -> count in modal +-NB bins
        piv = h.pivot_table(index=["a", "b"], columns="k", values="n", fill_value=0, aggfunc="sum")
        piv = piv.reindex(columns=range(2 * nbin + 1), fill_value=0).to_numpy()
        ker = np.ones(2 * NB + 1)
        sm = np.apply_along_axis(lambda r: np.convolve(r, ker, "same"), 1, piv) if len(piv) else piv
        kk = sm.argmax(1)
        idx = h.groupby(["a", "b"]).size().index
        rows.append(pd.DataFrame(dict(a=idx.get_level_values(0), b=idx.get_level_values(1),
                                      **{f"cnt_{c}": sm[np.arange(len(kk)), kk], f"lag_{c}": (kk - nbin) * BIN})))
    P = rows[0].merge(rows[1], on=["a", "b"], how="left").fillna({"cnt_s": 0})
    # symmetrise: use the direction a<b and average the two directions' counts
    P["lo"] = np.minimum(P.a, P.b); P["hi"] = np.maximum(P.a, P.b)
    Q = P.groupby(["lo", "hi"], as_index=False).agg(cnt_r=("cnt_r", "max"), cnt_s=("cnt_s", "max"),
                                                     lag_r=("lag_r", "first"))
    mn = np.minimum(n[Q.lo], n[Q.hi])
    Q["f_real"] = Q.cnt_r / mn; Q["f_shift"] = Q.cnt_s / mn; Q["excess"] = Q.f_real - Q.f_shift
    Q["a"] = gids[Q.lo]; Q["b"] = gids[Q.hi]; Q["n_a"] = n[Q.lo].astype(int); Q["n_b"] = n[Q.hi].astype(int)
    Q[Q.excess > 0.05][["a", "b", "n_a", "n_b", "lag_r", "f_real", "f_shift", "excess"]] \
        .rename(columns={"lag_r": "lag_s"}).to_csv(f"{a.out}/dedupe_pairs.csv", index=False)
    D = Q[Q.excess > THR]
    nb = {i: {} for i in range(Gn)}
    for lo, hi, e in zip(D.lo, D.hi, D.excess):
        nb[lo][hi] = e; nb[hi][lo] = e
    s = np.array([sc.get(g, -1) for g in gids], float)
    order = np.lexsort((np.arange(Gn), -n, -s))
    is_rep = np.zeros(Gn, bool); rep = np.arange(Gn)
    for i in order:
        best, be = -1, THR
        for j, e in nb[i].items():
            if is_rep[j] and e > be:
                best, be = j, e
        if best < 0:
            is_rep[i] = True
        else:
            rep[i] = best
    F = pd.DataFrame(dict(grp=gids, n_det=n.astype(int), score=s, is_rep=is_rep, rep=gids[rep],
                          partners=[";".join(str(gids[j]) for j in sorted(nb[i])) for i in range(Gn)],
                          absorbed=[";".join(str(gids[j]) for j in np.flatnonzero((rep == i) & ~is_rep))
                                    for i in range(Gn)]))
    F.to_csv(f"{a.out}/dedupe_families.csv", index=False)
    keep = set(gids[is_rep])
    det[det.grp.isin(keep)].to_parquet(f"{a.out}/det_dedup.parquet", index=False)
    print(f"[dedupe] {Gn} families, {len(D)} duplicate pairs -> {is_rep.sum()} representatives; "
          f"reps with score>=3: {(s[is_rep] >= 3).sum()} (all families >=3: {(s >= 3).sum()})")


if __name__ == "__main__":
    main()
