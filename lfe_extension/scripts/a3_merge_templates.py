#!/usr/bin/env python
"""Track A step 3 -- merge the two iterated template builds (non-PGC-site stations from the gated selection, PGC-site
stations from the leave-one-site-out selection) over the frozen 2012 set: a (family, station) pair takes the ITERATED
template where the iterated build made a usable one, and keeps the FROZEN template otherwise.  Weights are recomputed
over the merged usable set exactly as e5r_templates.py does (w = 10^(S-SNR/20), capped at 3x the median, median-normalised),
using the iterated S-SNR where the iterated template is used and the frozen S-SNR otherwise."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="a3_it1_nonPGC,a3_it1_PGC")
    ap.add_argument("--out-tag", default="a3_it1")
    a = ap.parse_args()
    zf = np.load(f"{IN}/e5r_templates_2012.npz", allow_pickle=True)
    tpl = dict(zf["tpl"].item()); sp = dict(zf["sp"].item())
    stk = {k: zf[k] for k in zf.files if "|" in k}
    qf = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    qf = qf[qf.usable & qf.template][["sta", "bkey", "snr_db", "snr_null_db", "s_pick_s", "n_events", "cc_vs_seisbench"]].copy()
    qf["source"] = "frozen"
    parts = []
    n_new = 0
    for p in a.parts.split(","):
        z = np.load(f"{OUT}/e5r_templates_{p}.npz", allow_pickle=True)
        t = z["tpl"].item(); s = z["sp"].item()
        tpl.update(t); sp.update(s); n_new += len(t)
        for k in z.files:
            if "|" in k:
                stk[k] = z[k]
        q = pd.read_csv(f"{OUT}/e5r_template_qc_{p}.csv")
        q = q[q.usable & q.template][["sta", "bkey", "snr_db", "snr_null_db", "s_pick_s", "n_events"] +
                                     (["cc_vs_seisbench"] if "cc_vs_seisbench" in q.columns else [])].copy()
        q["source"] = "iterated"
        parts.append(q)
    QN = pd.concat(parts, ignore_index=True)
    keys_new = set(zip(QN.bkey.astype(int), QN.sta))
    Q = pd.concat([QN, qf[~qf.set_index(["bkey", "sta"]).index.isin(keys_new)]], ignore_index=True)
    Q["usable"] = True; Q["template"] = True
    Q = Q.sort_values(["sta", "bkey"]).reset_index(drop=True)
    Q.to_csv(f"{OUT}/e5r_template_qc_{a.out_tag}.csv", index=False)
    keep = {(float(b), s) for b, s in zip(Q.bkey.astype(int), Q.sta)}
    TPL = {k: v for k, v in tpl.items() if k in keep}
    SP = {k: sp[k] for k in TPL}
    np.savez_compressed(f"{OUT}/e5r_templates_{a.out_tag}.npz", tpl=np.array(TPL, dtype=object),
                        sp=np.array(SP, dtype=object), **{k: v for k, v in stk.items()
                                                          if (float(k.split("|")[1]), k.split("|")[0]) in keep})
    w = 10 ** (Q.snr_db / 20.0); w = np.minimum(w, 3 * w.median()); w = w / w.median()
    pd.DataFrame(dict(fam=Q.bkey.astype(float), sta=Q.sta, w=w.to_numpy(), snr_db=Q.snr_db,
                      n_events=Q.n_events)).to_csv(f"{OUT}/e5r_weights_snr_{a.out_tag}.csv", index=False)
    print(f"[merge] {n_new} iterated templates over the frozen 1,144 -> {len(TPL)} pairs, {Q.bkey.nunique()} families; "
          f"iterated {int((Q.source=='iterated').sum())}, frozen kept {int((Q.source=='frozen').sum())}")
    print(Q.groupby(["sta", "source"]).size().unstack(fill_value=0).to_string())
    print(f"[write] {OUT}/e5r_templates_{a.out_tag}.npz, e5r_template_qc_{a.out_tag}.csv, e5r_weights_snr_{a.out_tag}.csv")


if __name__ == "__main__":
    main()
