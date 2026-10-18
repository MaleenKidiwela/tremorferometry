#!/usr/bin/env python
"""Score every locator variant against Bostock's PUBLISHED family coordinates.  Prefix `uu_`.

The calibration set: the 58 northern families that the frozen recipe rediscovered blind on the
2005 POLARIS bed (`results/repro/classify15.csv`, column `kept_bkey`), whose published
coordinates are `lfe_catalog_repair/results/e5r_families_2012.csv` (verified in
`uu_common.verify_truth_source` to be exactly the distinct (lat, lon, dep) triples carried by
Bostock et al. (2015)'s own catalogue export).  The 18 located SEALED families are scored once at
the end as a blind confirmation and no variant is chosen on them.

One representative template per family, FIXED before any variant is run (the template with the
highest one-to-one match fraction `frac` to that family), so every variant is scored on exactly
the same 58 objects.  Nothing is filtered: a family whose representative fails to locate in a
variant is reported as not located, not dropped.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import traceback

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                       # noqa: E402


def representatives():
    C = pd.read_csv(f"{U.RES}/repro/classify15.csv")
    F = U.bostock_families()
    out = []
    for col, tag in (("kept_bkey", "dev"), ("hid_bkey", "blind")):
        S = C[C[col] >= 0].copy()
        S["bkey"] = S[col]
        S = S[S.bkey.isin(F.index)]
        S = S.sort_values(["frac", "n_det"], ascending=False)
        R = S.groupby("bkey").head(1)[["grp", "bkey", "n_det", "n_days", "frac"]].copy()
        R["set"] = tag
        out.append(R)
    R = pd.concat(out, ignore_index=True)
    R["t_lat"] = F.lat[R.bkey].values
    R["t_lon"] = F.lon[R.bkey].values
    R["t_dep"] = F.dep[R.bkey].values
    return R


def noise_population():
    C = pd.read_csv(f"{U.RES}/repro/classify15.csv")
    return C[(C.kept_bkey < 0) & (C.hid_bkey < 0)][["grp", "n_det", "n_days"]].copy()


REQ_COLS = ("grp", "lat", "lon", "dep")


def why_not_a_location_table(L):
    """Return a short reason string if L cannot be scored, else "" (2026-09-23: uu_dd wrote its
    differential-time observations as loc_ZZdd_obs.csv, the glob below picked it up as a variant
    and the merge on `grp` killed the whole queue).  Nothing is dropped silently: the caller
    records the reason as a row in the ladder."""
    if L is None or not len(L):
        return "empty table"
    miss = [c for c in REQ_COLS if c not in L.columns]
    if miss:
        return "not a location table (missing " + ",".join(miss) + ")"
    return ""


def score(L, R):
    """L = a variant's location table; R = representatives with truth.  Returns per-family rows."""
    L = L.copy()
    if "located" not in L.columns:
        L["located"] = L.lat.notna() & L.lon.notna() & L.dep.notna()
    M = R.merge(L, on="grp", how="left", suffixes=("", "_loc"))
    M["located"] = M.located.fillna(False).astype(bool)
    M["located"] &= M.lat.notna() & M.lon.notna() & M.dep.notna()
    ok = M.located
    M["err_epi"] = np.nan
    M["err_dep"] = np.nan
    M["err_3d"] = np.nan
    if ok.any():
        e = U.epi_km(M.loc[ok, "lat"], M.loc[ok, "lon"], M.loc[ok, "t_lat"], M.loc[ok, "t_lon"])
        dz = M.loc[ok, "dep"].to_numpy(float) - M.loc[ok, "t_dep"].to_numpy(float)
        M.loc[ok, "err_epi"] = e
        M.loc[ok, "err_dep"] = dz
        M.loc[ok, "err_3d"] = np.sqrt(e ** 2 + dz ** 2)
    return M


def shuffled_null(M, n=200, seed=5):
    """Distance from each solution to a RANDOM family's published coordinate (the null the
    Chestler check reports): the same solutions, the truth labels permuted."""
    rng = np.random.default_rng(seed)
    ok = M[M.located].copy()
    if len(ok) < 3:
        return np.nan, np.nan
    med = []
    for _ in range(n):
        p = rng.permutation(len(ok))
        e = U.epi_km(ok.lat.to_numpy(), ok.lon.to_numpy(),
                     ok.t_lat.to_numpy()[p], ok.t_lon.to_numpy()[p])
        dz = ok.dep.to_numpy() - ok.t_dep.to_numpy()[p]
        med.append(np.median(np.sqrt(e ** 2 + dz ** 2)))
    return float(np.median(med)), float(np.quantile(med, 0.05))


def summarise(name, M, N=None):
    ok = M[M.located]
    d = dict(variant=name, n_families=len(M), n_located=int(len(ok)))
    if len(ok):
        d.update(epi_med=float(np.median(ok.err_epi)), epi_p90=float(np.quantile(ok.err_epi, .9)),
                 dep_med=float(np.median(np.abs(ok.err_dep))),
                 dep_p90=float(np.quantile(np.abs(ok.err_dep), .9)),
                 dep_bias=float(np.median(ok.err_dep)),
                 d3_med=float(np.median(ok.err_3d)),
                 rail=float(ok.edge.mean()) if "edge" in ok else np.nan,
                 within5=int((ok.err_epi <= 5).sum()), within10=int((ok.err_epi <= 10).sum()),
                 within20=int((ok.err_epi <= 20).sum()))
        nm, np05 = shuffled_null(M)
        d["null_d3_med"] = nm
    if N is not None and len(N):
        nl = N[N.located == True]                                           # noqa: E712
        d.update(noise_n=len(N), noise_located=int(len(nl)),
                 noise_rail=float(nl.edge.mean()) if len(nl) and "edge" in nl else np.nan)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--baseline", default=f"{U.RES}/repro/loc15_it0.csv")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    R = representatives()
    dev = R[R.set == "dev"]
    blind = R[R.set == "blind"]
    NOI = noise_population()
    R.to_csv(f"{a.out}/representatives.csv", index=False)
    json.dump(U.verify_truth_source(), open(f"{a.out}/truth_source.json", "w"), indent=1)

    rows, per = [], {}
    # --- the on-disk baseline, exactly as the chain produced it
    B = pd.read_csv(a.baseline)
    if "edge" not in B.columns:
        lats = np.arange(47.9, 49.31, 0.03)
        lons = np.arange(-124.9, -122.79, 0.04)
        deps = np.arange(18.0, 56.1, 2.0)
        B["edge"] = ((B.lat <= lats[0] + 1e-6) | (B.lat >= lats[-1] - 1e-6)
                     | (B.lon <= lons[0] + 1e-6) | (B.lon >= lons[-1] - 1e-6)
                     | (B.dep <= deps[0] + 1e-6) | (B.dep >= deps[-1] - 1e-6))
    bad = why_not_a_location_table(B)
    if bad:
        rows.append(dict(variant="BASELINE_on_disk(k_locate, tpl_it0)", n_families=len(dev),
                         n_located=0, status=bad))
        U.log(f"  BASELINE {a.baseline}: {bad}")
    else:
        M = score(B, dev)
        d = summarise("BASELINE_on_disk(k_locate, tpl_it0)", M,
                      NOI.merge(B, on="grp", how="left").assign(
                          located=lambda x: x.located.fillna(False)))
        d["status"] = "ok"
        rows.append(d)
        per["BASELINE_on_disk"] = M
        db = summarise("BASELINE_on_disk[BLIND sealed 18]", score(B, blind))
        db["status"] = "ok"
        rows.append(db)

    for p in sorted(glob.glob(f"{a.locdir}/loc_*.csv")):
        name = os.path.basename(p)[4:-4]
        try:
            L = pd.read_csv(p)
        except Exception as e:                                        # unreadable / truncated
            rows.append(dict(variant=name, n_families=len(dev), n_located=0,
                             status=f"unreadable: {type(e).__name__}: {e}"))
            U.log(f"  SKIP {name}: unreadable ({e})")
            continue
        bad = why_not_a_location_table(L)
        if bad:
            # A variant that produced NO solutions, or a file that is not a location table at
            # all, is REPORTED as such -- it never kills the queue.
            rows.append(dict(variant=name, n_families=len(dev), n_located=0, status=bad))
            U.log(f"  SKIP {name}: {bad} ({len(L) if L is not None else 0} rows, "
                  f"columns {list(L.columns)[:8] if L is not None else []})")
            continue
        try:
            M = score(L, dev)
            NL = NOI.merge(L, on="grp", how="left")
            NL["located"] = NL.located.fillna(False).astype(bool)
            d = summarise(name, M, NL)
            d["status"] = "ok" if d.get("n_located", 0) else "no family located"
            rows.append(d)
            per[name] = M
        except Exception as e:
            traceback.print_exc()
            rows.append(dict(variant=name, n_families=len(dev), n_located=0,
                             status=f"eval failed: {type(e).__name__}: {e}"))
            U.log(f"  SKIP {name}: eval failed ({type(e).__name__}: {e})")
    S = pd.DataFrame(rows)
    if "status" in S.columns:
        nb = S[S.status.fillna("") != "ok"]
        if len(nb):
            U.log(f"{len(nb)} of {len(S)} rows are not scorable and are reported, not dropped:\n"
                  + nb[["variant", "status"]].to_string(index=False))
    S.to_csv(f"{a.out}/ladder.csv", index=False)
    if per:
        PF = pd.concat([v.assign(variant=k) for k, v in per.items()], ignore_index=True)
        PF.to_csv(f"{a.out}/per_family.csv", index=False)
    else:
        U.log("no variant produced a scorable table; per_family.csv not written")
    with pd.option_context("display.width", 220, "display.max_columns", 40):
        U.log("\n" + S.round(2).to_string(index=False))
    U.log(f"wrote {a.out}/ladder.csv, per_family.csv")


if __name__ == "__main__":
    main()
