"""Q8d: attach the frozen-model P (already applied in production to every series) to the WA scrambled-control rows.
Writes a NEW file; edits nothing in lfe_extension_wa."""
import glob
import numpy as np, pandas as pd
ROOT = "/home/jovyan/tremorferometry"
WA = f"{ROOT}/lfe_extension_wa/results/wa_prod"
setup = pd.read_csv(f"{WA}/setup/families.csv"); fids = set(setup[setup.reg.isin(["NW", "SW"])].fid)
out, chk = [], []
for f in sorted(glob.glob(f"{WA}/post_v2/dedup/flags/flags_20*.parquet")):
    ym = f.split("_")[-1][:7]
    if ym < "2015-01":
        continue
    fl = pd.read_parquet(f, columns=["det_id", "date", "series", "fid", "k", "t", "dup_D"])
    fl = fl[fl.series.isin(["c1", "data"]) & fl.fid.isin(fids)]
    p = pd.read_parquet(f"{WA}/conf/P/prod/P_{ym}.parquet", columns=["date", "series", "fid", "k", "t", "P_own_WAVE_PNET", "reg"])
    p = p[p.series.isin(["c1", "data"]) & p.fid.isin(fids)]
    m = fl.merge(p, on=["date", "series", "fid", "k"], how="left", suffixes=("", "_p"))
    chk.append(dict(month=ym, n_c1=int((m.series == "c1").sum()), miss_c1=int(m[m.series == "c1"].P_own_WAVE_PNET.isna().sum()),
                    max_tdiff=float((m.t - m.t_p).abs().max())))
    out.append(m[m.series == "c1"][["det_id", "date", "series", "fid", "k", "t", "dup_D", "reg", "P_own_WAVE_PNET"]])
O = pd.concat(out, ignore_index=True)
O.to_parquet(f"{ROOT}/svi_science/results/q8_rtr/tables/wa_scrambled_P_2015_2026.parquet", index=False)
K = pd.DataFrame(chk); K.to_csv(f"{ROOT}/svi_science/results/q8_rtr/tables/wa_scrambled_P_join_check.csv", index=False)
print(len(O), K[["n_c1", "miss_c1"]].sum().to_dict(), "max tdiff", K.max_tdiff.max())
print(O.groupby("reg").P_own_WAVE_PNET.describe())
