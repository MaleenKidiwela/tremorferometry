#!/usr/bin/env python
"""Track C — per-station REGISTRATION of the T_core detections: how well does each station's stack of
the month's T_core detections match its template?

For every (family, station) it sums the daily stacks of `stacks_<tag>.npz` weighted by their event
counts, cuts the 12-s template window out of the 16-s stored window (the template starts `pre_s` into
it) and forms exactly the matched-filter statistic: the mean over components of the per-component
L2-normalised dot product with the template.  A null is formed at a random offset in the same window.

  python scripts/c_registration.py --tags 2018-06 2012-09ets
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C


def reg(tag, tpl, L, out_dir=None):
    p = f"{out_dir or C.TRC}/stacks_{tag}.npz"
    if not os.path.exists(p):
        raise SystemExit(f"[reg] {p} missing")
    z = np.load(p, allow_pickle=True)
    fam = z["fam"]; sta = z["sta"]; n = z["n"]; S = z["stack"]
    nsp = int(round(float(z["pre_s"]) * C.FS))
    nL = int(round(L * C.FS))
    rng = np.random.default_rng(7)
    acc = {}
    for i in range(len(fam)):
        acc.setdefault((float(fam[i]), str(sta[i])), [0.0, np.zeros_like(S[i])])
        acc[(float(fam[i]), str(sta[i]))][0] += float(n[i])
        acc[(float(fam[i]), str(sta[i]))][1] += float(n[i]) * S[i]
    rows = []
    for (f, s), (w, acc_s) in acc.items():
        if (f, s) not in tpl or w <= 0:
            continue
        T = tpl[(f, s)][L]                              # (3, nL), per-component L2-normalised
        st = acc_s / w
        seg = st[:, nsp:nsp + nL].astype(np.float64)
        nn = np.linalg.norm(seg, axis=1)
        if np.any(nn <= 0):
            continue
        cc = float(((seg / nn[:, None]) * T).sum(1).mean())
        j = int(rng.integers(0, max(st.shape[1] - nL, 1)))
        sg2 = st[:, j:j + nL].astype(np.float64); n2 = np.linalg.norm(sg2, axis=1)
        cc0 = float(((sg2 / np.maximum(n2, 1e-30)[:, None]) * T).sum(1).mean()) if sg2.shape[1] == nL else np.nan
        rows.append(dict(tag=tag, fam=int(f), sta=s, n_events=int(w), cc=cc, cc_null=cc0))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", nargs="+", required=True)
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--out", default=f"{C.TRC}/registration.csv")
    a = ap.parse_args()
    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, tpl_tag = C.resolve_templates(fz)
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    D = pd.concat([reg(t, tpl, L, a.out_dir or None) for t in a.tags], ignore_index=True)
    D.to_csv(a.out, index=False)
    pd.set_option("display.width", 250)
    piv = D.groupby(["sta", "tag"]).agg(n_fam=("fam", "nunique"), n_events=("n_events", "sum"),
                                        cc_median=("cc", "median"), cc_null_median=("cc_null", "median")).reset_index()
    print(piv.round(4).to_string(index=False))
    print(f"[reg] -> {a.out}")


if __name__ == "__main__":
    main()
