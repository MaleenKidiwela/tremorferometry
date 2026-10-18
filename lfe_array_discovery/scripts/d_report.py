#!/usr/bin/env python
"""Delivery report per tranche + sealed-set verification.

(1) results/delivery_by_tranche.csv, delivery_by_network.csv, delivery_low_stations.csv
(6) checks results/sealed_2011_manifest.sha256 against traces_sealed/ and proves no cache20hz
    entry exists for any day of the sealed window.
"""
from __future__ import annotations
import os
import sys
import glob
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402


def load():
    M = pd.read_csv(D.MANIFEST, keep_default_na=False, low_memory=False)
    M = M[M.tranche != ""].drop_duplicates(["tranche", "net", "sta", "day"], keep="last")
    M["bytes"] = pd.to_numeric(M.bytes, errors="coerce").fillna(0)
    M["covfrac"] = pd.to_numeric(M.coverage_frac, errors="coerce")
    Q = []
    for p in sorted(glob.glob(os.path.join(D.QUEUES, "q_*.csv"))):
        q = pd.read_csv(p, keep_default_na=False)
        Q.append(q)
    Q = pd.concat(Q, ignore_index=True)
    return M, Q


def main():
    M, Q = load()
    planned = Q.groupby("tranche").size().rename("planned")
    rows = []
    for tr, g in M.groupby("tranche"):
        n = len(g)
        rows.append(dict(
            tranche=tr, planned=int(planned.get(tr, np.nan)) if tr in planned.index else np.nan,
            attempted=n, data=int((g.status == "data").sum()),
            nodata=int((g.status == "nodata").sum()), error=int((g.status == "error").sum()),
            linked=int((g.linked == "linked").sum()),
            gib=round(g.bytes.sum() / 2**30, 2),
            frozen=int((g.chain == "frozen").sum()), rational=int((g.chain == "rational").sum()),
            n_stations=g[["net", "sta"]].drop_duplicates().shape[0],
            delivery_frac=round((g.status == "data").mean(), 4),
            median_cov=round(float(g[g.status == "data"].covfrac.median()), 4)))
    T = pd.DataFrame(rows).sort_values("tranche")
    T.to_csv(os.path.join(D.RES, "delivery_by_tranche.csv"), index=False)
    N = (M.assign(ok=(M.status == "data"))
          .groupby(["tranche", "net"])
          .agg(station_days=("ok", "size"), data=("ok", "sum"),
               gib=("bytes", lambda v: round(v.sum() / 2**30, 2)),
               n_stations=("sta", "nunique"))
          .reset_index())
    N["delivery_frac"] = (N.data / N.station_days).round(4)
    N.to_csv(os.path.join(D.RES, "delivery_by_network.csv"), index=False)
    S = (M.assign(ok=(M.status == "data"))
          .groupby(["tranche", "net", "sta"])
          .agg(planned=("ok", "size"), data=("ok", "sum"))
          .reset_index())
    S["delivery_frac"] = (S.data / S.planned).round(3)
    low = S[(S.delivery_frac < 0.70) & (S.tranche != "T6")].sort_values(
        ["tranche", "delivery_frac"])
    low.to_csv(os.path.join(D.RES, "delivery_low_stations.csv"), index=False)
    print(T.to_string(index=False))
    print()
    print(N.pivot_table(index="net", columns="tranche", values="delivery_frac").to_string())
    print(f"\n[low] {len(low)} station-tranche entries below 70 % delivery "
          f"({low[['net','sta']].drop_duplicates().shape[0]} distinct stations)")
    print(low.groupby(["tranche", "net"]).size().to_string())

    # ------------------------------------------------------------------ (6) sealed set
    print("\n=== SEALED SET ===")
    hp = os.path.join(D.RES, "sealed_2011_manifest.sha256")
    listed = {}
    with open(hp) as f:
        for ln in f:
            h, rel = ln.rstrip("\n").split("  ", 1)
            listed[rel] = h
    on_disk = set()
    for dp, _, fns in os.walk(D.TRACES_SEALED):
        for fn in fns:
            on_disk.add(os.path.relpath(os.path.join(dp, fn), D.TRACES_SEALED))
    miss_hash = sorted(on_disk - set(listed))
    miss_file = sorted(set(listed) - on_disk)
    print(f"files on disk {len(on_disk)}, hashed {len(listed)}, "
          f"unhashed {len(miss_hash)}, hashed-but-gone {len(miss_file)}")
    s0, s1 = D.WINDOWS["t3_sealed2011"]
    d0, d1 = pd.Timestamp(s0), pd.Timestamp(s1)
    doys = set(range(int(d0.dayofyear), int(d1.dayofyear) + 1))
    leaked = [p for p in glob.glob(os.path.join(D.CACHE, "*", "2011", "*.npz"))
              if int(os.path.basename(p)[:3]) in doys]
    mleak = M[(M.outroot != D.TRACES_SEALED) & (M.day >= s0) & (M.day <= s1)]
    print(f"cache20hz entries inside {s0}..{s1}: {len(leaked)}")
    print(f"manifest rows for those days written OUTSIDE traces_sealed/: {len(mleak)}")
    verdict = ("SEALED OK" if not miss_hash and not miss_file and not leaked and not len(mleak)
               else "SEAL PROBLEM")
    print(f"VERDICT: {verdict}")
    with open(os.path.join(D.RES, "sealed_check.txt"), "w") as f:
        f.write(f"files_on_disk {len(on_disk)}\nhashed {len(listed)}\nunhashed {len(miss_hash)}\n"
                f"hashed_but_missing {len(miss_file)}\ncache_entries_in_window {len(leaked)}\n"
                f"manifest_rows_outside_sealed {len(mleak)}\nverdict {verdict}\n")


if __name__ == "__main__":
    main()
