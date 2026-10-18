"""smoke: F10 matched-filter rows vs production rows on days where the networks are identical."""
import sys, json, numpy as np, pandas as pd
S, PR = sys.argv[1], "results/wa_prod/mf/prod"
pairs = [("2016-03-01", "core", "core"), ("2020-06-01", "core", "core")]
out = []
for d, nf, npd in pairs:
    A = pd.read_parquet(f"{S}/mf/prod/det_{d}.parquet"); B = pd.read_parquet(f"{PR}/det_{d}.parquet")
    A = A[A.net == nf]; B = B[B.net == npd]
    ja, jb = json.load(open(f"{S}/mf/prod/run_{d}.json")), json.load(open(f"{PR}/run_{d}.json"))
    same_mem = sorted(ja["members"][nf]) == sorted(jb["members"][npd])
    for ser in ("data", "c1", "rev"):
        a = A[A.series == ser].sort_values(["fid", "i"]).reset_index(drop=True); b = B[B.series == ser].sort_values(["fid", "i"]).reset_index(drop=True)
        r = dict(date=d, f10_net=nf, prod_net=npd, series=ser, same_members=same_mem, n_f10=len(a), n_prod=len(b))
        if len(a) == len(b) and (a.fid.values == b.fid.values).all() and (a.i.values == b.i.values).all():
            r["same_positions"] = True
            cols = [c for c in a.columns if c in b.columns and a[c].dtype.kind == "f"]
            r["max_abs_diff_float_cols"] = float(max(np.nanmax(np.abs(a[c].values - b[c].values)) if np.isfinite(a[c].values).any() else 0 for c in cols))
            r["max_abs_d_net_mad"] = float(np.abs(a.net_mad.values - b.net_mad.values).max()) if len(a) else 0
            r["n_float_cols"] = len(cols)
        else:
            r["same_positions"] = False
        out.append(r)
R = pd.DataFrame(out); print(R.to_string()); R.to_csv(f"{S}/../compare_prod_{S.rstrip('/').split('/')[-1]}.csv", index=False)
