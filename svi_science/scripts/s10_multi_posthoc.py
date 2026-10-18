#!/usr/bin/env python
"""s10_multi_posthoc -- post-hoc checks of the annual cycle (multi/PLAN.md Revision 1 item 5; added after it was seen):
per-station annual fits (episode stacks and day means), split halves 2003-2014 vs 2015-2026, residual chi2 after the
annual fit, largest residual episodes. Writes tables/posthoc_seasonal.csv and tables/posthoc.json."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_multi_common as M
import s10_multi_analysis as A


def fit(t_doy, x, sg):
    ph = 2 * np.pi * (t_doy + 10) / 365.25
    Am = np.vstack([np.cos(ph), np.sin(ph), np.ones(len(ph))]).T; w = 1 / sg ** 2
    b, *_ = np.linalg.lstsq(Am * np.sqrt(w)[:, None], x * np.sqrt(w), rcond=None)
    cov = np.linalg.pinv(Am.T @ (Am * w[:, None]))
    chi1 = np.sum(w * (x - Am @ b) ** 2); wm = np.sum(w * x) / w.sum(); chi0 = np.sum(w * (x - wm) ** 2)
    n = len(x); amp = float(np.hypot(b[0], b[1]))
    J = np.array([b[0], b[1]]) / max(amp, 1e-12); s_amp = float(np.sqrt(J @ cov[:2, :2] @ J * max(chi1 / max(n - 3, 1), 1)))
    return dict(n=n, amp=amp, amp_sig=s_amp, peak_doy=float((np.degrees(np.arctan2(b[1], b[0])) % 360) / 360 * 365.25 - 10) % 365.25,
                p_F=float(stats.f.sf(((chi0 - chi1) / 2) / (chi1 / max(n - 3, 1)), 2, max(n - 3, 1))),
                chi2red_after=float(chi1 / max(n - 3, 1)), chi2_after_p=float(stats.chi2.sf(chi1, max(n - 3, 1))))


def main():
    D, U = A.load()
    ES = A.es_frame(U)
    rows = []
    for unit, DD in (("es", ES), ("dm", D)):
        for grp, stas in [("core5", M.CORE5), ("all7", M.STATIONS)] + [(s, [s]) for s in M.STATIONS]:
            A.GROUPS[grp] = stas
            P, X = A.pooled(DD, "d24", grp, ["major"])
            if len(P) < 6:
                continue
            P["sig"] = A.pooled_boot(X, "d24", M.SEED).reindex(P.index)
            P = P.join(U.set_index("episode")[["start"]])
            P = P[(P.n_fam >= 3) & (P.sig > 0)]
            for half, m in (("all", np.ones(len(P), bool)), ("2003-2014", P.start < "2015"), ("2015-2026", P.start >= "2015")):
                Q = P[m]
                if len(Q) < 6:
                    continue
                r = fit(pd.to_datetime(Q.start).dt.dayofyear.to_numpy(), Q.value.to_numpy(), Q.sig.to_numpy())
                rows.append(dict(unit=unit, stations=grp, years=half, **r))
    T = pd.DataFrame(rows); T.to_csv(f"{M.TB}/posthoc_seasonal.csv", index=False)
    print(T.round(4).to_string())
    # largest residual episodes after the annual fit (core5, episode stacks)
    R = json.load(open(f"{M.TB}/result.json"))
    out = {}
    for tag in ("es", "dm"):
        S = pd.read_csv(f"{M.TB}/pooled_series_{tag}.csv", index_col=0)
        S = S[(S.cls == "major") & (S.core5_d24_n_fam >= 3)]
        b0, b1 = R[f"annual_model_{tag}"]; ph = 2 * np.pi * (pd.to_datetime(S.start).dt.dayofyear + 10) / 365.25
        S["resid"] = S.core5_d24_value - (b0 * np.cos(ph) + b1 * np.sin(ph)); S["z"] = S.resid / S.core5_d24_sigma_famboot
        out[tag] = S.sort_values("z", key=np.abs, ascending=False)[["start", "core5_d24_value", "resid", "core5_d24_sigma_famboot", "z",
                                                                    "core5_d24_n_fam"]].head(6).round(4).reset_index().to_dict("records")
    json.dump(out, open(f"{M.TB}/posthoc.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
