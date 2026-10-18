#!/usr/bin/env python
"""bw_prod_kelb_templates -- era template set with KELB added (user 2026-10-09): KELB = the KLNB site before its
2004-06-03 rename (same coordinates to 0.0005 deg, same Guralp CMG3ESP_100sec / Nanometrics Trident).  KELB gets the KLNB
era templates, placements (sp; the detector recomputes travel times from KELB's own coordinates) and weights unchanged.
Metadata difference: overall gain 6.17e7 (2002-09-12 -> 2003-07-31) and 3.07e8 (-> 2004-06-03) vs KLNB 3.97e8.  The cc
is amplitude-normalised, so only the rel_amp feature sees the gain; KELB's amplitude reference is KLNB's x (3.07/3.97)^2
(later KELB epoch; before 2003-07-31 KELB rel_amp is therefore ~1/25 too small -- noted, not corrected).
Out: svi_science/results/backward/templates/kelb/g3_era_templates.npz, g3_era_weights.csv, g3_era_ampref.csv (+ SHA256SUMS)."""
import numpy as np, pandas as pd
T = "/home/jovyan/tremorferometry/svi_science/results/backward/templates"
z = np.load(f"{T}/g3_era_templates.npz", allow_pickle=True)
tpl, sp = z["tpl"].item(), z["sp"].item()
add = {(f, "KELB"): v for (f, s), v in tpl.items() if s == "KLNB"}
tpl.update(add); sp.update({(f, "KELB"): sp[(f, "KLNB")] for (f, _) in add})
np.savez_compressed(f"{T}/kelb/g3_era_templates.npz", tpl=np.array(tpl, dtype=object), sp=np.array(sp, dtype=object))
W = pd.read_csv(f"{T}/g3_era_weights.csv"); W = pd.concat([W, W[W.sta == "KLNB"].assign(sta="KELB")])
W.to_csv(f"{T}/kelb/g3_era_weights.csv", index=False)
A = pd.read_csv(f"{T}/g3_era_ampref.csv")
A = pd.concat([A, A[A.sta == "KLNB"].assign(sta="KELB", e3_ref=lambda d: d.e3_ref * (3.07305e8 / 3.96523e8) ** 2)])
A.to_csv(f"{T}/kelb/g3_era_ampref.csv", index=False)
print(f"KELB templates added: {len(add)}; total era templates {len(tpl)}")
