#!/usr/bin/env python
"""r23b_f1_posthoc (2026-10-11, POST-HOC, PLAN Revision 2): is the joint 2003-26 'persistence' (S3 rule 7a) a property
of lasting patches or of mixing full-box and half-box ETS in the pooled share? Rule 7a rerun on full-box big ETS only
(tremor in-reach >= 0.9 of patches), and on pair types. Writes tables/f1/posthoc_persistence.json."""
import importlib.util, json, sys
import numpy as np, pandas as pd
sys.argv = [sys.argv[0]]
ROOT = "/home/jovyan/tremorferometry"
spec = importlib.util.spec_from_file_location("f1", f"{ROOT}/svi_science/scripts/r23b_f1.py")
f1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f1)
W, PE, FE, D = f1.load()
out = {}
for per in ["2003-26", "2014-26"]:
    des = f1.design(D, W, "ceq", per)
    cov = f1.Q.covariates(des, D, "SVI", "tremor"); inr, _, _, _ = f1.Q.kin_rows(des, cov)
    frac = inr.mean(0); eps = [des.keys[j][1] for j in des.tslots]
    full = [e for e, f in zip(eps, frac) if f >= 0.9]
    d2 = f1.design(D, W, "ceq", per, targets=full)
    mid = (d2.H.w0 + (d2.H.w1 - d2.H.w0) / 2).to_numpy()[d2.tslots]
    out[f"{per}_fullbox"] = dict(targets=full, **f1.s3_tests(d2.C[:, d2.tslots].astype(float), mid))
    out[f"{per}_in_reach_frac"] = dict(zip(eps, np.round(frac, 2).tolist()))
    print(per, {k: out[f"{per}_fullbox"][k] for k in ["n_episodes", "mean_resid_corr", "labelshuffle_mean", "p_persist_labelshuffle",
                                                       "p_persist_multinomial", "mean_consecutive_corr", "p_catchup_ordershuffle"]}, flush=True)
json.dump(out, open(f"{ROOT}/svi_science/results/r23b/tables/f1/posthoc_persistence.json", "w"), indent=1, default=float)

# second post-hoc check: drop the era-detector ETS (they dominate the pooled share through their large counted-once N)
for name, sel in [("2003-26_noERA", lambda W_: W_.big & (W_.seggroup != "ERA")),
                  ("2003-26_fullbox_noERA", lambda W_: W_.big & (W_.seggroup != "ERA") & W_.episode.isin(out["2003-26_fullbox"]["targets"])),
                  ("2003-26_equalN", None)]:
    if sel is None:
        des = f1.design(D, W, "ceq", "2003-26")
        C = des.C[:, des.tslots].astype(float); C = C / C.sum(0) * 1000.0      # every ETS weighted equally in the pool
        mid = (des.H.w0 + (des.H.w1 - des.H.w0) / 2).to_numpy()[des.tslots]
        out[name] = f1.s3_tests(C, mid)
    else:
        t = W[sel(W)].episode.tolist(); d2 = f1.design(D, W, "ceq", "2003-26", targets=t)
        mid = (d2.H.w0 + (d2.H.w1 - d2.H.w0) / 2).to_numpy()[d2.tslots]
        out[name] = dict(targets=t, **f1.s3_tests(d2.C[:, d2.tslots].astype(float), mid))
    print(name, {k: out[name][k] for k in ["n_episodes", "mean_resid_corr", "labelshuffle_mean", "p_persist_labelshuffle",
                                          "p_persist_multinomial", "mean_consecutive_corr", "p_catchup_ordershuffle"]}, flush=True)
json.dump(out, open(f"{ROOT}/svi_science/results/r23b/tables/f1/posthoc_persistence.json", "w"), indent=1, default=float)
