#!/usr/bin/env python
"""p3_figures (2026-10-07) -- figures for P3: Bostock-era timeline, recharge curves by era, quiet month (three ends)."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/p3_bostock"
FIG = f"{ROOT}/svi_science/figures"
S0 = f"{ROOT}/svi_science/results/s0"
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK2 = "#52514e"; GRID = "#e6e5e0"
DEEP_GROUPS = [1, 4, 10, 12]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": "#0b0b0b", "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})


def timeline():
    F = pd.read_csv(f"{S0}/families.csv")
    deep = set(F.fid[F.group.isin(DEEP_GROUPS)])
    M = pd.read_parquet(f"{OUT}/bostock_events_mapped.parquet")
    M["deep"] = M.fid.isin(deep)
    T = pd.read_parquet(f"{S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8"])
    T = T[(T.series == "data") & (T.P_same >= 0.75) & T.keep_xfam8]
    T["deep"] = T.fam.isin(deep)
    Eb = pd.read_csv(f"{OUT}/episodes_bostock_final.csv")
    fig = plt.figure(figsize=(11, 7.2))
    ax = fig.add_axes([0.07, 0.62, 0.9, 0.32])
    for src, lab0 in ((M, "Bostock"), (T, "ours")):
        mm = src.assign(m=pd.to_datetime(src.date).dt.to_period("M")).groupby(["m", "deep"]).size().unstack(fill_value=0)
        x = mm.index.to_timestamp()
        ax.bar(x, mm.get(False, 0), width=25, color=C1, label=f"{lab0}: other families" if lab0 == "ours" else None)
        ax.bar(x, mm.get(True, 0), bottom=mm.get(False, 0), width=25, color=C2, label=f"{lab0}: deep groups" if lab0 == "ours" else None)
    ax.set_yscale("log"); ax.set_ylim(0.8, 3e5)
    ax.axvspan(pd.Timestamp("2003-01-01"), pd.Timestamp("2014-01-01"), color=GRID, alpha=0.6, zorder=0)
    ax.text(pd.Timestamp("2003-06-01"), 9e4, "2003–2013 Bostock catalogue: only the 10 ETS windows were scanned\n(grey = not scanned, NOT quiet)",
            fontsize=8, color=INK2, va="top")
    ax.text(pd.Timestamp("2014-03-01"), 9e4, "2014–2026 ours: every day scanned (P ≥ 0.75, counted once;\nrecall ~0.2 of Bostock's)",
            fontsize=8, color=INK2, va="top")
    ax.set_ylabel("LFEs per month (log)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C1, label="other families"), Patch(color=C2, label="deep groups 1/4/10/12")],
              loc="lower left", bbox_to_anchor=(0.55, 1.0), frameon=False, ncol=2)
    ax.set_title("Bostock-era LFE timeline (2003–2013) and ours (2014–2026)", loc="left", fontsize=10)
    eps = Eb[Eb.major].reset_index(drop=True)
    for i, r in eps.iterrows():
        a = fig.add_axes([0.07 + (i % 5) * 0.185, 0.31 - (i // 5) * 0.26, 0.16, 0.19])
        d = M[(M.date >= r.start) & (M.date <= r.end)].groupby(["date", "deep"]).size().unstack(fill_value=0)
        x = pd.to_datetime(d.index)
        a.bar(x, d.get(False, 0), width=0.85, color=C1)
        a.bar(x, d.get(True, 0), bottom=d.get(False, 0), width=0.85, color=C2)
        sh = d.get(True, 0).sum() / d.values.sum()
        a.set_title(f"{r.episode}: {int(r.n_events):,} LFEs, deep {100*sh:.0f} %", fontsize=8, loc="left")
        a.tick_params(labelsize=6.5)
        a.xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=7))
        a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%m-%d"))
        if i % 5 == 0:
            a.set_ylabel("LFEs / day", fontsize=7)
    fig.savefig(f"{FIG}/p3_bostock_timeline.png", dpi=150); plt.close(fig)


def model(t, r):
    return 1 - r["a"] * np.exp(-t / r["tau"])


def recharge():
    C = pd.read_csv(f"{OUT}/recharge_curves.csv")
    Fi = pd.read_csv(f"{OUT}/recharge_fits.csv").set_index("name")
    mids = {"0-7": 4, "7-14": 11, "14-30": 22, "30-60": 45, "60-90": 75, "90-150": 120, "150-240": 195, "240-400": 320}
    tt = np.linspace(1, 400, 400)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    panels = [
        ("A. Ours 2014–2026, LFEs (deep vs up-dip)", [("LFE deep groups 1/4/10/12 | ref 39 majors | D1", C1, "deep groups"),
                                                      ("LFE up-dip groups 3/6/8/11 | ref 39 majors | D1", C3, "up-dip groups (no fit)"),
                                                      ("LFE deep groups 1/4/10/12 | ref big list | D1", C4, "deep, clock = big ETS only")]),
        ("B. Ours, same detector: first vs second half", [("LFE deep groups | ref 39 majors | D1 | 2014-2019", C1, "2014–2019"),
                                                          ("LFE deep groups | ref 39 majors | D1 | 2020-2026", C2, "2020–2026")]),
        ("C. PNSN deep tremor: 2010–2013 vs 2014–2026", [("PNSN deep tremor | PNSN majors | D1 | 2010-2013", C2, "2010–2013"),
                                                         ("PNSN deep tremor | PNSN majors | D1 | 2014-2026", C1, "2014–2026")]),
    ]
    for ax, (ttl, items) in zip(axs, panels):
        for nm, col, lab in items:
            if nm not in Fi.index:
                continue
            c = C[C.set == nm]; r = Fi.loc[nm]
            x = [mids[b] for b in c.bin]
            mean = r.n_events / r.exposure
            ax.plot(x, c.rate / mean, "o", color=col, ms=6, mec="white", mew=1,
                    label=f"{lab}: {int(r.n_events)} events, τ {r.tau:.0f} d" + (" (at fit limit)" if r.tau >= 399 else ""))
            if r.n_events >= 100:
                ax.plot(tt, r.r_inf * model(tt, r) / mean, "-", color=col, lw=2)
        ax.axhline(1, color=INK2, lw=0.8, ls=":")
        ax.axvspan(0, 30, color=GRID, alpha=0.5, zorder=0)
        ax.set_xscale("log"); ax.set_xlim(2, 420)
        ax.set_xlabel("days since the last major ETS ended")
        ax.set_title(ttl, loc="left", fontsize=9)
        ax.legend(frameon=False, fontsize=7, loc="upper left")
    axs[0].set_ylabel("rate / that curve's mean rate")
    axs[0].set_ylim(0, 3.2)
    fig.suptitle("Recharge curves: deep activity after a major ETS (each curve divided by its own mean rate = shape only; grey = first 30 days; lines = fitted recovery curve)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/p3_recharge_curves.png", dpi=150); plt.close(fig)


def quiet():
    Q = pd.read_csv(f"{OUT}/quiet_month_three_ends.csv")
    Q = Q[~Q.set.str.contains("broad")]
    sets = Q.set.unique()
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ends = ["D1 regional end", "D2 end + 7 d", "D3 family's own end"]
    cols = [C1, C2, C3]
    for i, s in enumerate(sets):
        for k, (e, col) in enumerate(zip(ends, cols)):
            r = Q[(Q.set == s) & (Q.end_definition == e)].iloc[0]
            y = i + (k - 1) * 0.22
            ax.plot([r.shuffle_F30_p01, r.shuffle_F30_median], [y, y], color=GRID, lw=5, solid_capstyle="round")
            ax.plot(r.shuffle_F30_median, y, "|", color=INK2, ms=10)
            ax.plot(r.F30, y, "o", color=col, ms=7, mec="white", mew=1, label=e if i == 0 else None)
            ax.text(0.56, y, f"p={r.p_low:.3f}, n={int(r.n)}", fontsize=7, va="center", color=INK2)
    lab = {s: s.replace(" (P>=0.75 once, 2014-2026)", " 2014–26 (discovery)").replace("WA E3 BE", "WA 2015–26 (confirmatory)")
              .replace("WA E2 BE", "WA 2012–14").replace("WA E1 BE", "WA 2009–11") for s in sets}
    ax.set_yticks(range(len(sets))); ax.set_yticklabels([lab[s] for s in sets]); ax.invert_yaxis()
    ax.set_xlim(0, 0.75)
    ax.set_xlabel("share of between-episode LFEs in the first 30 days after a major ETS ends\n(dot = observed; tick = random-day expectation; grey bar down to its 1st percentile)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Quiet month after ETS, three end definitions\n(Bostock 2003–2013: not testable, ETS windows only)", loc="left", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/p3_quiet_month.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    timeline(); quiet(); recharge()
