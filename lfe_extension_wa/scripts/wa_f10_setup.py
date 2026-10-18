#!/usr/bin/env python
"""WA FIXED-10 re-run (2026-10-07): build the run's own setup under OUT = $WA_PROD_DIR (results/wa_fixed10/) from the
production setup (results/wa_prod/), which is only READ.  Nothing in results/wa_prod/ is written.

  setup/   copies of the production setup tables (same days, exclusions, families, truth, sealed truth, scan days, day
           roles, raw inventories, node map), with three files restricted:
             stations.csv      only the 10 fixed stations (the matched filter loads only these)
             tpl_library.csv   only the W4a templates of the 10 stations (same frozen templates, same order)
             core_sets.json    E1 = E2 = E3 = the 9 PB boreholes (network 'core' = fixed 9, one sensor and one frozen
                               template set each, 2007-2026); network 'all' = every loaded station = the literal fixed 10
           + fixed_sets.json   F9 / F10 lists, GNW template epochs, per-family fixed template stations (strict days)
           + COPIED_FROM.csv   source path, sha256 of source and copy
  cache/   copies of the production cache manifests (the matched filter reads them to tell 'nodata' from 'missing')
  post/    copy of post/tremor_family_day.parquet (PNSN tremor <= 30 / 60 km per family-day; station independent)

  python scripts/wa_f10_setup.py
"""
from __future__ import annotations
import os, sys, json, glob, shutil
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P

SRC = f"{P.RES}/wa_prod"
OUT = P.OUT
F9 = ["PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007", "PB.B013", "PB.B017", "PB.B018"]
F10 = F9 + ["UW.GNW"]
GNW_EPOCHS = [("2000-06-21", "2010-11-30", "UW.GNW..BH.50@2000-06-21"),      # own 'build' templates (260 families)
              ("2010-12-01", "2019-05-06", "UW.GNW..BH.40@2010-12-01"),      # t2012new / step4b re-templates (52 families)
              ("2019-05-07", "2099-12-31", "UW.GNW..BH.40@2010-12-01")]      # HH.100 node carries BH.40 (W2-approved)
COPY = ["day_table.csv", "days_prod.json", "days_truth.json", "exclusions.csv", "families.csv", "raw_inventory_dl2.parquet",
        "raw_inventory_nodl2.parquet", "scan_days.csv", "truth.parquet", "truth_array_frac.csv", "tpl_nodemap.csv",
        "core_set_table.csv", "b014_check.json", "sealed/truth_sealed.parquet", "sealed/README.txt"]


def main():
    assert os.path.realpath(OUT) != os.path.realpath(SRC), "OUT must not be the production directory"
    P.ensure(f"{OUT}/setup/sealed", f"{OUT}/cache", f"{OUT}/post")
    rows = []

    def cp(src, dst):
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
        rows.append(dict(src=src, dst=dst, sha_src=P.sha256(src), sha_dst=P.sha256(dst)))
    for f in COPY:
        if os.path.exists(f"{SRC}/setup/{f}"):
            cp(f"{SRC}/setup/{f}", f"{OUT}/setup/{f}")
    for m in sorted(glob.glob(f"{SRC}/cache/manifest_*.csv")):
        cp(m, f"{OUT}/cache/{os.path.basename(m)}")
    cp(f"{SRC}/post/tremor_family_day.parquet", f"{OUT}/post/tremor_family_day.parquet")
    # restricted tables
    S = pd.read_csv(f"{SRC}/setup/stations.csv", keep_default_na=False)
    S = S[S.netsta.isin(F10)]
    assert sorted(S.netsta) == sorted(F10), S.netsta.tolist()
    S.to_csv(f"{OUT}/setup/stations.csv", index=False)
    L = pd.read_csv(f"{SRC}/setup/tpl_library.csv", dtype={"live": str})
    L10 = L[L.netsta.isin(F10)]
    L10.to_csv(f"{OUT}/setup/tpl_library.csv", index=False)
    CS = json.load(open(f"{SRC}/setup/core_sets.json"))
    for e in ("E1", "E2", "E3"):
        missing = [s for s in F10 if s not in CS[e]]
        assert not missing, (e, missing)               # the 10 are in every production era set
    J = dict(E1=F9, E2=F9, E3=F9, rule="FIXED-10 re-run (2026-10-07): network 'core' = the 9 PB boreholes in every era "
             "(one sensor + one frozen template set each); network 'all' = every loaded station = the literal fixed 10 "
             "(F9 + UW.GNW, whose templates change with its sensor in 2010-12). Stations = production E1 n E2 n E3.",
             raw_inventory=CS.get("raw_inventory"), production_sets={e: CS[e] for e in ("E1", "E2", "E3")})
    json.dump(J, open(f"{OUT}/setup/core_sets.json", "w"), indent=1)
    # per-family fixed template stations (strict days)
    fam = {}
    for f, g in L10.groupby("fid"):
        pb = sorted(set(g[g.netsta.isin(F9)].netsta))
        gn = g[g.netsta == "UW.GNW"]
        fam[int(f)] = dict(F9=pb, GNW_nodes=sorted(set(gn.node)))
    json.dump(dict(F9=F9, F10=F10, GNW_epochs=GNW_EPOCHS, families=fam), open(f"{OUT}/setup/fixed_sets.json", "w"), indent=1)
    pd.DataFrame(rows).to_csv(f"{OUT}/setup/COPIED_FROM.csv", index=False)
    bad = [r for r in rows if r["sha_src"] != r["sha_dst"]]
    assert not bad, bad
    P.log(f"[f10 setup] {len(rows)} files copied (sha256 equal); stations {len(S)}; templates {len(L10)} of {len(L)} "
          f"({L10.groupby('netsta').size().to_dict()}); core = F9 every era")


if __name__ == "__main__":
    main()
