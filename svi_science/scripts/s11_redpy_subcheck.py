#!/usr/bin/env python
"""s11_redpy_subcheck -- PLAN §1 follow-ups: station consistency (adjusted Rand index of sub-family labels on shared events,
vs 500 label shuffles), stack-vs-stack cc between sub-families (same waveform or not), and whether membership tracks
single-event SNR (cc with the series grand stack)."""
import os, sys, itertools
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R
from sklearn.metrics import adjusted_rand_score as ari
fam = int(sys.argv[1]) if len(sys.argv) > 1 else 10
L = pd.read_parquet(f"{R.TB}/subfam_labels_fam{fam}.parquet")
rng = np.random.default_rng(5); out = []
E1 = L[L.epoch == "E1"]
for s1, s2 in itertools.combinations(["CN.PGC", "PB.B009", "PB.B010", "PB.B011"], 2):
    a = E1[E1.netsta == s1].set_index("row_id"); b = E1[E1.netsta == s2].set_index("row_id")
    sh = a.index.intersection(b.index)
    for col in ("labA", "labB"):
        x, y = a.loc[sh, col].to_numpy(), b.loc[sh, col].to_numpy()
        if col == "labA":
            x, y = (x == 1).astype(int), (y == 1).astype(int)      # biggest sub-family vs rest
        r0 = ari(x, y); perm = np.array([ari(rng.permutation(x), y) for _ in range(500)])
        out.append(dict(sta1=s1, sta2=s2, labels=col if col == "labB" else "labA(sub1 vs rest)", n_shared=len(sh), ari=r0,
                        perm95=float(np.percentile(perm, 95)), p=float((np.sum(perm >= r0) + 1) / 501),
                        same_frac=float((x == y).mean())))
A = pd.DataFrame(out); A.to_csv(f"{R.TB}/subfam_station_ari_fam{fam}.csv", index=False); print(A.round(3).to_string())
# stack-vs-stack cc and SNR association
rows = []
for ns in ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]:
    D = R.load_series(fam, ns)
    for ep, s in D.items():
        q = L[(L.netsta == ns) & (L.epoch == ep)].set_index("row_id").loc[s["g"].row_id]
        F = R.feat(s["W"], s["tS"]); Gs = F.sum(0); ccg = (F @ Gs - (F * F).sum(1)) / np.linalg.norm(Gs[None] - F, axis=1)
        for col in ("labA", "labB"):
            lab = q[col].to_numpy(); ks = [k for k in np.unique(lab) if (lab == k).sum() >= 20]
            St = {k: R.wstack(s["W"][lab == k], np.ones((lab == k).sum())) for k in ks}
            Fs = {k: R.feat(St[k][None], s["tS"])[0] for k in ks}
            for k in ks:
                others = [j for j in ks if j != k]
                rows.append(dict(netsta=ns, epoch=ep, labels=col, sub=int(k), size=int((lab == k).sum()),
                                 ccgrand_med=float(np.median(ccg[lab == k])), ccgrand_med_rest=float(np.median(ccg[lab != k])),
                                 stack_cc_max_other=float(max(Fs[k] @ Fs[j] / np.linalg.norm(Fs[k]) / np.linalg.norm(Fs[j]) for j in others)) if others else np.nan,
                                 stack_cc_min_other=float(min(Fs[k] @ Fs[j] / np.linalg.norm(Fs[k]) / np.linalg.norm(Fs[j]) for j in others)) if others else np.nan,
                                 P_med=float(q.P[lab == k].median())))
B = pd.DataFrame(rows); B.to_csv(f"{R.TB}/subfam_stackcc_fam{fam}.csv", index=False); print(B.round(3).to_string())
