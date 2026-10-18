#!/usr/bin/env python
"""P3 REDISCOVERY SCORER — RELATIVE-CLAUSE version of scripts/g_rediscover.py (Merlin flaw #5).

g_rediscover.py judged a discovered cluster against ABSOLUTE clauses (>= 30 % one-to-one recall of
the hidden family's events, >= 50 % member purity).  Both clauses sit at or above the frozen TEMPLATE
detector's own per-family level — its per-family recall is 0.18-0.20 and its per-family purity 0.395
on 2012 — so the old scorer was asking a template-free detector to beat the templates before it could
be called a rediscovery, and it scored the templates' own clusters at only 9/30.

This version makes clauses (b) and (c) RELATIVE to the TRUTH TEMPLATE's own per-family values,
measured on the SAME year, the SAME window and with the SAME stations withheld:

    (a) cluster centroid within `r_km` (5 km) of the hidden family            [ABSOLUTE]
    (b) one-to-one recall (+-2 s) >= 0.5 x recall_ref(family)                 [RELATIVE]
    (c) member purity           >= 0.8 x purity_ref(family)                   [RELATIVE]
    (d) the cluster is not > 50 % explained by a KEPT family's events         [ABSOLUTE]

`recall_ref` / `purity_ref` come from `truth_reference()` below, applied to a frozen-template run
produced by scripts/q_run_month.py with `--drop-stations SNB,B926` (the frozen operating-point json
and the frozen template npz are read unchanged; only the loaded (family, station) tables are
filtered, in-process).  Families the truth template cannot detect at all on that network
(recall_ref == 0) are reported separately as NO-REFERENCE: a relative clause against zero is
meaningless, and such a family must not be counted either as a rediscovery or as a failure without
saying so.

    python scripts/q_rediscover.py --reference-only \
        --det results/p3/trackQ/detections_noSNBB926_2012-*.parquet --year 2012 \
        --windows 2012-08-15,2012-10-05
    python scripts/q_rediscover.py --clusters C.parquet --reference results/p3/truth_template_reference.csv ...
"""
from __future__ import annotations
import argparse, glob, json, os, sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
from g_rediscover import one2one, in_windows                     # noqa: E402  (same recipe)

RES = f"{G.RESULTS}"
P3 = f"{RES}/p3"


# ------------------------------------------------------------------ the truth-template reference
def truth_reference(det, events, windows=None, tol=2.0, tier_col="cat_core"):
    """Per-family recall and purity of the FROZEN TEMPLATE detector itself.

    det:    frozen-detector detections with columns fam, t (origin axis) and a tier flag;
    events: {bkey: array of Bostock t_origin}.
    Returns a DataFrame indexed by bkey with n_events, n_det, recall_ref, purity_ref."""
    rows = []
    for f, ev in events.items():
        ev = np.sort(np.asarray(ev, float)[in_windows(ev, windows)])
        d = det[det.fam.astype(int) == int(f)]
        if tier_col and tier_col in d.columns:
            d = d[d[tier_col].astype(bool)]
        td = np.sort(d.t.to_numpy(float)[in_windows(d.t.to_numpy(float), windows)])
        mr, _ = one2one(ev, td, tol)          # family side -> recall
        mp, _ = one2one(td, ev, tol)          # detector side -> purity
        rows.append(dict(bkey=int(f), n_events=int(len(ev)), n_det=int(len(td)),
                         recall_ref=float(mr.mean()) if len(ev) else np.nan,
                         purity_ref=float(mp.mean()) if len(td) else np.nan))
    return pd.DataFrame(rows).set_index("bkey")


# ------------------------------------------------------------------ the scorer
def score(clusters, events, families, hidden, reference, windows=None, tol=2.0, r_km=5.0,
          recall_frac=0.5, purity_frac=0.8, max_explained=0.50, split_frac=0.20,
          min_ref_recall=1e-9):
    hidden = set(hidden)
    kept = [f for f in families if f not in hidden]
    ev = {f: np.sort(np.asarray(t, float)[in_windows(t, windows)]) for f, t in events.items()}
    cl = {c: dict(times=np.sort(np.asarray(v["times"], float)[in_windows(v["times"], windows)]),
                  lat=v["lat"], lon=v["lon"]) for c, v in clusters.items()}
    expl = {}
    for c, v in cl.items():
        expl[c] = {f: int(one2one(v["times"], t, tol)[0].sum()) for f, t in ev.items()}
    rows = []
    for h in sorted(hidden):
        if h not in families:
            continue
        rr = float(reference.recall_ref.get(h, np.nan))
        pr = float(reference.purity_ref.get(h, np.nan))
        need_r = recall_frac * rr
        need_p = purity_frac * pr
        no_ref = not np.isfinite(rr) or rr <= min_ref_recall
        best = None
        for c, v in cl.items():
            n = len(v["times"])
            if n == 0:
                continue
            d = float(G.epi_km(v["lat"], v["lon"], *families[h]))
            mh, _ = one2one(ev[h], v["times"], tol)
            recall = float(mh.mean()) if len(ev[h]) else 0.0
            purity = expl[c][h] / n
            kbest = max((expl[c][k] / n for k in kept), default=0.0)
            kwho = max(kept, key=lambda k: expl[c][k], default=None) if kept else None
            ok = (not no_ref) and (d <= r_km and recall >= need_r and purity >= need_p
                                   and kbest <= max_explained)
            r = dict(hidden=h, cluster=c, dist_km=d, n_cluster=n, n_family=len(ev[h]),
                     recall=recall, purity=purity, recall_ref=rr, purity_ref=pr,
                     need_recall=need_r, need_purity=need_p,
                     kept_explained=kbest, kept_family=kwho,
                     no_reference=bool(no_ref), rediscovered=bool(ok))
            if best is None or (ok, recall) > (best["rediscovered"], best["recall"]):
                best = r
        rows.append(best or dict(hidden=h, cluster=None, rediscovered=False,
                                 no_reference=bool(no_ref), recall_ref=rr, purity_ref=pr))
    R = pd.DataFrame(rows)
    srows = []
    for k in kept:
        n = len(ev[k])
        if n == 0:
            continue
        parts = []
        for c, v in cl.items():
            nc = len(v["times"])
            if nc == 0:
                continue
            frac = float(one2one(ev[k], v["times"], tol)[0].mean())
            pur = expl[c][k] / nc
            pk = float(reference.purity_ref.get(k, np.nan))
            if frac >= split_frac and np.isfinite(pk) and pur >= purity_frac * pk:
                parts.append(c)
        srows.append(dict(family=k, n_events=n, n_parts=len(parts), split=bool(len(parts) >= 2),
                          parts=",".join(map(str, parts))))
    S = pd.DataFrame(srows)
    testable = R[~R.no_reference] if len(R) else R
    summ = dict(
        n_hidden=int(len(hidden)), n_kept=int(len(kept)), n_clusters=int(len(cl)),
        n_no_reference=int(R.no_reference.sum()) if len(R) else 0,
        n_testable=int(len(testable)),
        n_rediscovered=int(R.rediscovered.sum()) if len(R) else 0,
        rediscovery_rate_of_testable=(float(testable.rediscovered.mean())
                                      if len(testable) else np.nan),
        n_kept_split=int(S.split.sum()) if len(S) else 0,
        spurious_split_rate=float(S.split.mean()) if len(S) else np.nan,
        criterion=dict(r_km=r_km, tol_s=tol, recall_frac=recall_frac, purity_frac=purity_frac,
                       max_explained=max_explained, split_frac=split_frac,
                       relative_to="frozen template, same year/window, SNB+B926 withheld"))
    return R, S, summ


# ------------------------------------------------------------------ CLI
def load_events(year):
    sys.path.insert(0, HERE)
    import p_common as P
    B = P.bostock()
    B = B[B.year == year]
    return {int(b): g.t_origin.to_numpy() for b, g in B.groupby("bkey")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-only", action="store_true")
    ap.add_argument("--det", default=f"{P3}/trackQ/detections_noSNBB926_2012-*.parquet")
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--tier", default="cat_core", choices=["cat_core", "cat_ext"])
    ap.add_argument("--windows", default="2012-08-15,2012-10-05")
    ap.add_argument("--reference", default=f"{P3}/truth_template_reference.csv")
    ap.add_argument("--clusters", default="")
    ap.add_argument("--out", default=f"{P3}/q_rediscovery")
    a = ap.parse_args()
    windows = []
    for w in a.windows.split(";"):
        if w:
            x, y = w.split(",")
            windows.append((pd.Timestamp(x).timestamp(), pd.Timestamp(y).timestamp()))
    events = load_events(a.year)
    if a.reference_only or not os.path.exists(a.reference):
        files = sorted(glob.glob(a.det))
        if not files:
            raise SystemExit(f"no detection files match {a.det}")
        det = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
        ref = truth_reference(det, events, windows, tier_col=a.tier)
        os.makedirs(os.path.dirname(a.reference), exist_ok=True)
        ref.to_csv(a.reference)
        hid = {int(h["bkey"]): h["stratum"]
               for h in json.load(open(f"{RES}/holdouts.json"))["hidden_families"]}
        ref["stratum"] = pd.Series(hid).reindex(ref.index).fillna("kept")
        print(f"[reference] {len(files)} detection files, tier {a.tier}, window {a.windows}")
        print(ref.groupby("stratum").agg(n=("n_events", "size"), events=("n_events", "sum"),
                                         det=("n_det", "sum"),
                                         recall=("recall_ref", "median"),
                                         purity=("purity_ref", "median")).to_string())
        print("\nhidden families:")
        print(ref[ref.stratum != "kept"].sort_values(["stratum", "recall_ref"]).to_string())
        print(f"\n[done] -> {a.reference}")
        if a.reference_only:
            return
    ref = pd.read_csv(a.reference).set_index("bkey")
    C = (pd.read_parquet(a.clusters) if a.clusters.endswith("parquet")
         else pd.read_csv(a.clusters))
    sys.path.insert(0, HERE)
    import p_common as P
    F = P.families()
    families = {int(b): (float(F.lat[b]), float(F.lon[b])) for b in F.index}
    hidden = [int(h["bkey"]) for h in json.load(open(f"{RES}/holdouts.json"))["hidden_families"]]
    clusters = {c: dict(times=g.t.to_numpy(float), lat=float(g.lat.iloc[0]),
                        lon=float(g.lon.iloc[0])) for c, g in C.groupby("cid")}
    R, S, summ = score(clusters, events, families, hidden, ref, windows)
    R.to_csv(a.out + ".csv", index=False)
    S.to_csv(a.out + "_splits.csv", index=False)
    json.dump(summ, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
