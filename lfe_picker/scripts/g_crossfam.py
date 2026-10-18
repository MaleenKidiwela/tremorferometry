#!/usr/bin/env python
"""Per-station CROSS-FAMILY cc p95 (the Gate-0 chance level) and the LFE-core HF-ratio distribution.

Track A-2 step 4 measured the cross-family chance level only for the six stations it tested
(PGC 0.67, NLLB 0.52, SNB 0.50, VGZ 0.42, LRIV 0.42, B926 0.32).  Gate 0 needs it at EVERY station,
so it is recomputed here from the frozen 2012 from-scratch stacks: for every station, correlate each
family's stack against every OTHER family's 12-s template over +-1 s and take the p95 of |cc|.

Writes results/gates/g_crossfam_cc.csv.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402


def main():
    Z = np.load(f"{G.REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    TPL = Z["tpl"].item()
    qc = pd.read_csv(f"{G.REP}/results/e5r_template_qc_2012.csv")
    qc = qc[qc.usable & qc.template]
    rows = []
    for sta, g in qc.groupby("sta"):
        fams = sorted(int(b) for b in g.bkey.unique())
        stacks, tpls = {}, {}
        for b in fams:
            k = f"{sta}|{b}"
            if k in Z.files and (float(b), sta) in TPL and 12.0 in TPL[(float(b), sta)]:
                stacks[b] = np.asarray(Z[k], float)
                tpls[b] = np.asarray(TPL[(float(b), sta)][12.0], float)
        fams = sorted(stacks)
        cc = []
        for i in fams:
            for j in fams:
                if i == j:
                    continue
                cc.append(abs(G.cc_template(stacks[i], tpls[j])))
        cc = np.asarray(cc, float)
        if cc.size == 0:
            continue
        rows.append(dict(sta=sta, n_families=len(fams), n_pairs=cc.size,
                         cc_med=float(np.median(cc)), cc_p90=float(np.quantile(cc, .90)),
                         cc_p95=float(np.quantile(cc, .95)), cc_p99=float(np.quantile(cc, .99)),
                         frac_ge_0p6=float((cc >= 0.6).mean())))
    R = pd.DataFrame(rows).sort_values("cc_p95", ascending=False)
    os.makedirs(G.GATES, exist_ok=True)
    R.to_csv(f"{G.GATES}/g_crossfam_cc.csv", index=False)
    print(R.to_string(index=False))


if __name__ == "__main__":
    main()
