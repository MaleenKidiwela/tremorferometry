#!/usr/bin/env python
"""SH_KMRFIX task 5: per-window KMR flags for the cut windows (results/cut/sh_ext/index_YYYY-MM.parquet), without
touching the HDF5 files.  bad_for_dvv = NC.KMR HHE/HHN from 2025-07-25 00:00 (DC drift of the failing STS-2 starts
2025-07-25; HHE stuck at a constant from 2025-08-02/03) and every NC.KMR channel from 2025-09-17 00:00 (sensor swaps;
no usable ground motion, results/SH_CHECKS_RESULTS.md).  kmr_z_caution = NC.KMR HHZ 2025-07-25 .. 2025-09-16 (mass
drift and weakened LFE signal, results/sh_kmrfix/t1_kmr_channel_stacks_summary.csv).
Output results/cut/sh_ext/window_flags_v2.parquet (one row per index row: det_id, seed_id, h5_file, h5_start, t0, flags)."""
import glob
import numpy as np, pandas as pd
CUT = "/home/jovyan/tremorferometry/lfe_south/results/cut/sh_ext"
H0 = pd.Timestamp("2025-07-25").timestamp(); SW = pd.Timestamp("2025-09-17").timestamp()
parts = []
for f in sorted(glob.glob(f"{CUT}/index_*.parquet")):
    I = pd.read_parquet(f, columns=["det_id", "netsta", "cha", "seed_id", "t0", "h5_file", "h5_start", "n_samples"])
    k = I.netsta.eq("NC.KMR").to_numpy(); hz = I.cha.isin(["HHE", "HHN"]).to_numpy(); t = I.t0.to_numpy()
    bad_h = k & hz & (t >= H0) & (t < SW); bad_all = k & (t >= SW)
    I["bad_for_dvv"] = bad_h | bad_all
    I["flag_reason"] = np.where(bad_all, "KMR_after_2025-09-17_swap_no_ground_motion",
                                np.where(bad_h, "KMR_STS2_horizontal_failing_from_2025-07-25", ""))
    I["kmr_z_caution"] = k & I.cha.eq("HHZ").to_numpy() & (t >= H0) & (t < SW)
    parts.append(I.drop(columns=["netsta", "cha"]))
W = pd.concat(parts, ignore_index=True)
W.to_parquet(f"{CUT}/window_flags_v2.parquet", index=False)
print(len(W), "windows;", int(W.bad_for_dvv.sum()), "bad_for_dvv;", W.flag_reason.value_counts().to_dict(), "z caution", int(W.kmr_z_caution.sum()))
