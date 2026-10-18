#!/usr/bin/env python
"""r23_ramp (2026-10-11) -- R23 study 5: pre-ETS ramp-up of deep families, 2003-2026 (S7 test 2 rules + scrambled-date control).
Rules: svi_science/results/r23/PLAN.md §6.    python svi_science/scripts/r23_ramp.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C

OUT = f"{C.TB}/study5"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/study5.log")
RNG = np.random.default_rng(20261011)
NSH = 1000
BINS = [0, 7, 14, 30, 60, 120, 240, 5000]; BLAB = ["0-7", "7-14", "14-30", "30-60", "60-120", "120-240", ">240"]


def main():
    F = C.families(); E = C.episodes()
    V = C.load_events(F, E); D = C.load_famdays(F, E)
    maj = E[E.major].reset_index(drop=True); big = E[E.big].reset_index(drop=True)
    clocks = {"all majors": maj, "big ETS": big}
    rows, epo = [], []
    sets = {"deep BE_broad": ("BE_broad", True), "all-family BE (S7)": ("BE", False), "deep BE": ("BE", True)}
    for p, (a, b) in C.PERIODS.items():
        for sname, (flag, deep_only) in sets.items():
            dm = D[flag] & (D.deep if deep_only else True)
            for cname, ref in clocks.items():
                st, en = ref.start_ts.values, ref.end_ts.values
                Dp = D[dm].copy()
                Dp["tn"], Dp["tp"], _ = C.tau_clock(Dp.date, st, en)
                for cond in ("all", "prev>=30"):
                    Dq = Dp if cond == "all" else Dp[Dp.tp >= 30]
                    pool = {k: g.date.to_numpy() for k, g in Dq.groupby(["fam", "seg"])}
                    for ser in ("data", "c1"):
                        ev = V[(V.date >= a) & (V.date <= b) & V[flag] & (V.series == ser) & ((V.deep) if deep_only else True)].copy()
                        ev["tn"], ev["tp"], _ = C.tau_clock(ev.date, st, en)
                        if cond == "prev>=30":
                            ev = ev[ev.tp >= 30]
                        ev = ev[np.isfinite(ev.tn) & np.array([k in pool for k in zip(ev.fam, ev.seg)], bool)]
                        if len(ev) == 0:
                            continue
                        f30 = float(np.mean(ev.tn < 30)); f14 = float(np.mean(ev.tn < 14))
                        row = dict(period=p, set=sname, clock=cname, cond=cond, series=ser, n=int(len(ev)), F30=f30, F14=f14)
                        if ser == "data":
                            sh = C.shuffle_dates(list(zip(ev.fam, ev.seg)), pool, RNG, NSH)
                            tn = C.tau_clock(sh.ravel(), st, en)[0].reshape(NSH, -1)
                            f30s = np.array([np.mean(r[np.isfinite(r)] < 30) for r in tn]); f14s = np.array([np.mean(r[np.isfinite(r)] < 14) for r in tn])
                            row.update(shuffle_F30_median=float(np.median(f30s)), shuffle_F30_p99=float(np.percentile(f30s, 99)), p_hi=C.pct_hi(f30, f30s),
                                       shuffle_F14_median=float(np.median(f14s)), p14_hi=C.pct_hi(f14, f14s))
                            # burst level (each burst once at its first event)
                            evs = ev.sort_values("time_utc").copy(); evs["burst"] = C.chain_bursts(evs.time_utc)
                            fb = evs.groupby("burst").first()
                            shb = C.shuffle_dates(list(zip(fb.fam, fb.seg)), pool, RNG, NSH)
                            tnb = C.tau_clock(shb.ravel(), st, en)[0].reshape(NSH, -1)
                            fbo = float(np.mean(fb.tn < 30)); fbs = np.array([np.mean(r[np.isfinite(r)] < 30) for r in tnb])
                            row.update(n_bursts=int(len(fb)), F30_bursts=fbo, shuffle_F30_bursts_median=float(np.median(fbs)), p_hi_bursts=C.pct_hi(fbo, fbs))
                            if cond == "all":
                                obs = np.histogram(ev.tn, BINS)[0]; expo = np.histogram(Dq[(Dq.date >= a) & (Dq.date <= b)].tn.dropna(), BINS)[0]
                                shh = np.mean([np.histogram(r[np.isfinite(r)], BINS)[0] for r in tn], 0)
                                for lab, o, e_, s_ in zip(BLAB, obs, expo, shh):
                                    epo.append(dict(period=p, set=sname, clock=cname, bin=lab, events=int(o), elig_famdays=int(e_),
                                                    per_1000=1000 * o / max(e_, 1), shuffle_mean=float(s_)))
                        rows.append(row)
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/ramp_F30.csv", index=False)
    pd.DataFrame(epo).to_csv(f"{OUT}/superposed_epoch.csv", index=False)
    log("F30 before next onset:\n" + A.round(4).to_string(index=False))
    # ---------------- control 2: scrambled pre-window position inside the same inter-major gap
    TR = pd.read_parquet(f"{C.TB}/tremor_points.parquet", columns=["date", "src", "dmin_deep_km"])
    TR = TR[(TR.dmin_deep_km <= 15) & (((TR.src == "env") & (TR.date < C.PNSN_T0)) | ((TR.src == "pnsn") & (TR.date >= C.PNSN_T0)))]
    trd = TR.groupby("date").size()
    c2 = []
    for p, (a, b) in C.PERIODS.items():
        for cname, ref in clocks.items():
            ev = V[V.BE_broad & V.deep & (V.series == "data")]
            evd = ev.groupby("date").size()
            dd = D[D.BE_broad & D.deep].groupby("date").size()
            gaps = []
            for i in range(1, len(ref)):
                g0 = ref.end_ts.iloc[i - 1] + pd.Timedelta(days=3); g1 = ref.start_ts.iloc[i] - pd.Timedelta(days=1)
                if not (a <= str(ref.start_ts.iloc[i].date()) <= b) or ref.start_ts.iloc[i] > pd.Timestamp(C.P_END):
                    continue
                days = pd.date_range(g0, g1).strftime("%Y-%m-%d")
                if len(days) < 1:
                    continue
                y = evd.reindex(days).fillna(0).to_numpy(); n = dd.reindex(days).fillna(0).to_numpy(); tt = trd.reindex(days).fillna(0).to_numpy()
                gaps.append((ref.episode.iloc[i], y, n, tt))
            if not gaps:
                continue
            L = 30

            def stat(offsets):
                Y = N = TT = 0.0
                for (nm, y, n, tt), o in zip(gaps, offsets):
                    w = min(L, len(y)); Y += y[o:o + w].sum(); N += n[o:o + w].sum(); TT += tt[o:o + w].sum()
                return Y, N, TT
            obs = stat([max(len(y) - L, 0) for _, y, _, _ in gaps])
            sims = np.array([stat([RNG.integers(0, max(len(y) - L, 0) + 1) for _, y, _, _ in gaps]) for _ in range(NSH)])
            ro = obs[0] / max(obs[1], 1); rs = sims[:, 0] / np.maximum(sims[:, 1], 1)
            c2.append(dict(period=p, clock=cname, n_gaps=len(gaps), events_prewindow=obs[0], famdays_prewindow=obs[1], rate_per_1000=1000 * ro,
                           scrambled_rate_median=1000 * float(np.median(rs)), p_hi=C.pct_hi(ro, rs),
                           tremor_prewindow=obs[2], scrambled_tremor_median=float(np.median(sims[:, 2])), p_hi_tremor=C.pct_hi(obs[2], sims[:, 2])))
    C2 = pd.DataFrame(c2); C2.to_csv(f"{OUT}/ramp_scrambled_dates.csv", index=False)
    log("control 2 (scrambled pre-window position in the same gap):\n" + C2.round(4).to_string(index=False))
    # ---------------- verdict per period (deep BE_broad, all majors clock)
    ver = {}
    for p in C.PERIODS:
        r1 = A[(A.period == p) & (A.set == "deep BE_broad") & (A.clock == "all majors") & (A.cond == "all") & (A.series == "data")].iloc[0]
        r2 = A[(A.period == p) & (A.set == "deep BE_broad") & (A.clock == "all majors") & (A.cond == "prev>=30") & (A.series == "data")].iloc[0]
        c = C2[(C2.period == p) & (C2.clock == "all majors")].iloc[0]
        sup = r1.F30 > r1.shuffle_F30_p99 and r2.F30 > r2.shuffle_F30_p99 and c.p_hi < 0.01
        notsup = r1.F30 <= r1.shuffle_F30_median
        ver[p] = dict(F30=r1.F30, shuffle_median=r1.shuffle_F30_median, p=r1.p_hi, F30_prev30=r2.F30, p_prev30=r2.p_hi, c2_p=c.p_hi,
                      verdict="supported" if sup else ("not supported" if notsup else "undecided"))
    log(f"VERDICT: {ver}")
    C.dump(dict(verdict=ver), f"{OUT}/study5_summary.json")
    figure(pd.DataFrame(epo))
    log("done")


def figure(EP):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for k, p in enumerate(("2003-13", "2014-26")):
        e = EP[(EP.period == p) & (EP.set == "deep BE_broad") & (EP.clock == "all majors")]
        x = np.arange(len(e))
        ax[k].bar(x - 0.2, e.events, 0.4, color="#08519c", label="observed")
        ax[k].bar(x + 0.2, e.shuffle_mean, 0.4, color="#bdbdbd", label="shuffled days (mean)")
        ax[k].set_xticks(x); ax[k].set_xticklabels(e.bin, fontsize=8); ax[k].set_xlabel("days before the next major ETS starts")
        ax[k].set_ylabel("deep-family LFEs between episodes"); ax[k].set_title(f"{p}"); ax[k].legend(fontsize=8)
    fig.suptitle("Pre-ETS ramp-up test: deep families between episodes vs time to the next major")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/r23_preets_ramp.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
