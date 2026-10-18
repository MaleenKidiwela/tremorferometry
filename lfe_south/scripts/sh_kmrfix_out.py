#!/usr/bin/env python
"""SH_KMRFIX: joinable outputs (new files only): results/cut/sh_ext/detections_v2.parquet (every >= 8x v2 candidate:
det_id_v2 + the cut-window det_id to read it from) and results/cut/sh_ext/conf_P_v2_by_det_id.parquet (P(same)/P(any) v2)."""
import os, json
import numpy as np, pandas as pd
RES = "/home/jovyan/tremorferometry/lfe_south/results"; OUT = f"{RES}/sh_kmrfix"; CUT = f"{RES}/cut/sh_ext"
for p in [f"{CUT}/detections_v2.parquet", f"{CUT}/conf_P_v2_by_det_id.parquet"]:
    assert not os.path.exists(p), p
R = pd.read_parquet(f"{OUT}/candidates8x_v2.parquet"); R["fam"] = R.fam.astype(str)
R["key"] = np.round(R.t_epoch * 20).astype(np.int64)
F = pd.read_parquet(f"{OUT}/conf_v2_rows.parquet"); F = F[(F.run == "fwd") & (F.series == "data")].copy()
F["key"] = np.round(F.t_epoch * 20).astype(np.int64)
FN = F[F.seg == "NK"].set_index("det_id_v2"); FO = F[F.seg != "NK"]
cols = ["P_same_v2", "P_any_v2", "P_model_v2", "feat_src", "P_same_FULL", "P_any_FULL", "role", "m_same", "m_any"]
a = R[R.seg != "NK"].merge(FO[["fam", "key"] + cols], on=["fam", "key"], how="left")
b = R[R.seg == "NK"].join(FN[cols], on="det_id_v2")
J = pd.concat([a, b], ignore_index=True).sort_values("det_id_v2").reset_index(drop=True)
assert J.det_id_v2.is_unique and len(J) == len(R) and J.P_same_v2.notna().all(), (len(J), len(R), J.P_same_v2.isna().sum())
J = J.rename(columns={"P_same_FULL": "P_same_v1model", "P_any_FULL": "P_any_v1model", "role": "conf_role",
                      "m_same": "his_same_fam_2s_nearest", "m_any": "his_any_fam_2s_nearest"})
J["P_same_v2a"] = np.where(J.t < pd.Timestamp("2025-07-25"), J.P_same_v1model, J.P_same_v2)
J["P_any_v2a"] = np.where(J.t < pd.Timestamp("2025-07-25"), J.P_any_v1model, J.P_any_v2)
J["in_sample"] = J.conf_role.eq("insample")
J["det_id"] = J.det_id.where(J.map_type != "none", -1).astype(np.int64)
J["window_unmapped"] = J.map_type.eq("none")
det = J[["det_id_v2", "t", "t_epoch", "fam", "seg", "net_mad", "ccsum", "meancc", "nlive", "nchan", "nsta", "in_S1_9x_v2",
         "det_id", "map_type", "map_dt_s", "map_fam", "window_unmapped"]].rename(columns={"net_mad": "xmad"})
det.to_parquet(f"{CUT}/detections_v2.parquet", index=False)
conf = J[["det_id_v2", "det_id", "fam", "t_epoch", "seg", "in_S1_9x_v2", "map_type", "P_same_v2", "P_any_v2", "P_model_v2", "feat_src",
          "P_same_v2a", "P_any_v2a", "P_same_v1model", "P_any_v1model", "conf_role", "in_sample", "his_same_fam_2s_nearest", "his_any_fam_2s_nearest"]]
conf.to_parquet(f"{CUT}/conf_P_v2_by_det_id.parquet", index=False)
# summaries
nk = J[J.seg == "NK"]
s = dict(n_v2_candidates=len(J), n_S1=int(J.in_S1_9x_v2.sum()),
         NK_rows=len(nk), NK_P_model=nk.P_model_v2.value_counts().to_dict(), NK_S1_P_model=nk[nk.in_S1_9x_v2].P_model_v2.value_counts().to_dict(),
         NK_sumP_same_share_fallback=float(nk[nk.P_model_v2 == "FALLBACK_v2"].P_same_v2.sum() / nk.P_same_v2.sum()),
         NK_unmapped_windows=int(nk.window_unmapped.sum()), NK_unmapped_S1=int((nk.window_unmapped & nk.in_S1_9x_v2).sum()),
         NK_unmapped_sumP_same=float(nk[nk.window_unmapped].P_same_v2.sum()),
         meanP_same_by_maptype=nk.groupby("map_type").P_same_v2.mean().round(3).to_dict(),
         meanP_same_by_model=nk.groupby("P_model_v2").P_same_v2.mean().round(3).to_dict())
json.dump(s, open(f"{OUT}/t5_outputs_summary.json", "w"), indent=1, default=float); print(json.dumps(s, indent=1, default=float))
