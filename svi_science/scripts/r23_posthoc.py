#!/usr/bin/env python
"""r23_posthoc (2026-10-11) -- POST-HOC column for R23 study 3 (PLAN Revision 3): S7 quiet-month F30 restricted to the
deep groups (1, 4, 10, 12), fallback and union episode lists. Not a verdict input.   python svi_science/scripts/r23_posthoc.py
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C
OUT = f"{C.TB}/posthoc"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/posthoc.log")
RNG = np.random.default_rng(20261011)


def main():
    F = C.families(); E0 = C.episodes()
    sys.argv.append("--union")
    import r23_quiet as Q
    rows = []
    for lab, E in (("fallback (pre-registered)", E0), ("union (diagnostic)", Q.union_majors(E0))):
        V = C.load_events(F, E); D = C.load_famdays(F, E)
        for setname, flag in (("deep BE", "BE"), ("deep BE_broad", "BE_broad")):
            pool = C.pools(D, (D[flag] & D.deep).to_numpy())
            for p, (a, b) in C.PERIODS.items():
                ev = V[(V.date >= a) & (V.date <= b) & V[flag] & V.deep & (V.series == "data")]
                sh = C.shuffle_dates(list(zip(ev.fam, ev.seg)), pool, RNG, 1000)
                for cname, ref in (("all majors", E[E.major]), ("big ETS", E[E.big])):
                    tp = C.tau_clock(ev.date, ref.start_ts.values, ref.end_ts.values)[1]
                    ts = C.tau_clock(sh.ravel(), ref.start_ts.values, ref.end_ts.values)[1].reshape(1000, -1)
                    m = np.isfinite(tp); f = float(np.mean(tp[m] < 30)); fs = np.array([np.mean(r[np.isfinite(r)] < 30) for r in ts])
                    f60 = float(np.mean(tp[m] < 60)); fs60 = np.array([np.mean(r[np.isfinite(r)] < 60) for r in ts])
                    rows.append(dict(episodes=lab, set=setname, period=p, clock=cname, n=int(m.sum()), F30=f, shuffle_F30_median=float(np.median(fs)),
                                     p_low=C.pct_lo(f, fs), F60=f60, shuffle_F60_median=float(np.median(fs60)), p60_low=C.pct_lo(f60, fs60)))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/deep_quiet_month_F30.csv", index=False)
    log("POST-HOC deep-only quiet month:\n" + R.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
