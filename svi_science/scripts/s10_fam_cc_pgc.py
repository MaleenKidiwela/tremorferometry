"""Per-detection template correlation at PGC (cc_PGC, z_PGC) for one family, all detections (no P cut).
Joins catalogue rows (run, date, k) to the matched-filter feature files. User request 2026-10-11."""
import sys, os, numpy as np, pandas as pd
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
R = "/home/jovyan/tremorferometry"
b = pd.read_parquet(f"{R}/data/catalogue/svi_catalogue_v1_4_backward.parquet",
                    columns=["series", "fam", "date", "P", "primary", "counted_once", "detector", "stations_day", "run", "k"])
b = b[(b.fam == FAM) & (b.series == "data") & b.primary & b.stations_day.astype(str).str.contains("PGC")]
f = pd.read_parquet(f"{R}/lfe_extension/results/refine/conf/P_by_candidate.parquet",
                    columns=["series", "fam", "date", "P_same", "keep_xfam8", "run", "k"])
f = f[(f.fam == FAM) & (f.series == "data") & (f.run == "prod")].rename(columns={"P_same": "P", "keep_xfam8": "counted_once"})
rows = pd.concat([b.assign(src="backward"), f.assign(src="forward")], ignore_index=True)
out = []
for (run, date, src), g in rows.groupby(["run", "date", "src"]):
    base = f"{R}/lfe_extension/results/refine/feat/{run}" if src == "forward" else f"{R}/svi_science/results/backward/mf/feat/{run}"
    p = f"{base}/{date}.parquet"
    if not os.path.exists(p):
        continue
    F = pd.read_parquet(p)
    cols = [c for c in ["k", "cc_PGC", "z_PGC", "net_mad", "net_cc", "n_sta", "site_cc_max"] if c in F.columns]
    out.append(g.merge(F[cols], on="k", how="left"))
X = pd.concat(out, ignore_index=True); X["counted_once"] = X.counted_once.astype(float).fillna(0).astype(bool); X["primary"] = X.primary.astype(float).fillna(1).astype(bool); X["stations_day"] = X.stations_day.astype(str)
X.to_parquet(f"{R}/svi_science/results/s10_dvv/pilot/tables/fam{FAM}_cc_PGC.parquet", index=False)
print("rows", len(rows), "joined", len(X), "with cc_PGC", int(X.cc_PGC.notna().sum()))
