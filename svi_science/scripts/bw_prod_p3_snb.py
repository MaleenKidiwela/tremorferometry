#!/usr/bin/env python
"""bw_prod_p3_snb -- P3: SNB era templates for the CMG40T epoch (before 2006-04-20), the one marginal G2 case.
Built exactly like the G3 era templates (e5r_templates.do_station) from Bostock's 2004-07 BUILD-block events re-timed
onto the frozen axis by the OTHER CN stations (PGC, VGZ, NLLB; bw_gate_retime.py).  Weights on the frozen 2012 scale.
Out: svi_science/results/backward/templates/snb40t_2004.npz (+ _weights.csv, _qc.csv)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bw_gate_g2_epoch as E, bw_gate_retime as RT
G = E.G
OD = "/home/jovyan/tremorferometry/svi_science/results/backward/templates"
S = E.setup(); T = S["T"]
ev = RT.retime(E.events("2004-07", build_blocks=True), 2004, [G.CACHE_BW], ["PGC", "VGZ", "NLLB"])
ev = ev.assign(bkey=ev.fam.astype(int))
F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv"); cnt = ev.groupby("bkey").size()
F = F.assign(n_build=F.bkey.map(cnt).fillna(0).astype(int))
fam_tbl = F[F.n_build >= 30][["bkey", "sb_fam_id", "lat", "lon", "dep", "n_build"]].to_dict("records")
T.YEAR = 2004; T.CACHES = [G.CACHE_BW]; T.QUIET_DOYS = range(156, 170); T.LENS = (8.0, 12.0)
st = S["stab"]
sta, rows, TPL, SP, STK, ref = T.do_station(("SNB", "CN.SNB", float(st.loc["SNB", "lat"]), float(st.loc["SNB", "lon"]), "surface_bb",
                                             fam_tbl, ev[["bkey", "t_origin"]], {}, 777))
Q = pd.DataFrame(rows); Q.to_csv(f"{OD}/snb40t_2004_qc.csv", index=False)
use = Q[Q.usable & Q.template]
TPLu = {k: v for k, v in TPL.items() if ((use.bkey == int(k[0])) & (use.sta == k[1])).any()}
np.savez_compressed(f"{OD}/snb40t_2004.npz", tpl=np.array(TPLu, dtype=object), sp=np.array({k: SP[k] for k in TPLu}, dtype=object))
W12 = pd.read_csv(f"{G.REP}/results/e5r_weights_snr_2012.csv"); M = float(np.median(10 ** (W12.snr_db / 20.0)))
pd.DataFrame(dict(fam=use.bkey.astype(float), sta="SNB", w=(np.minimum(10 ** (use.snr_db / 20.0), 3 * M) / M).to_numpy())).to_csv(
    f"{OD}/snb40t_2004_weights.csv", index=False)
G.log(f"[p3 snb] {len(TPLu)} usable SNB-40T era templates of {len(Q)} stacks")
