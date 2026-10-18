#!/usr/bin/env python
"""g_trigger_figs (2026-10-11) -- Study G summary figures (descriptive; the tests are in g_trigger_tests.py).
  g_superposed.png   counts in 30-min bins, t_P - 24 h .. t_P + 24 h, strong between-episode trials, per region:
                     counted-once LFEs (data / scrambled control) and all candidates (data / scrambled), vs the same
                     clock on 200 random same-segment between-episode days per trial (mean and 1-99 % band)
  g_dose.png         per-trial LFE response vs peak dynamic stress (instantaneous; delayed +1..+10 d)
  g_ets_onsets.png   major ETS onsets and earthquakes (stress at the region) through time
  g_sigma_calib.png  estimated vs measured PGC stress
    python svi_science/scripts/g_trigger_figs.py
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import g_trigger_tests as GT

ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger"
T = f"{G}/tables"; F = f"{G}/figures"
NS = 10 ** 9; H = 3600 * NS; DAY = 86400 * NS
RNG = np.random.default_rng(1)
BINS = np.arange(-24 * H, 24 * H + 1, H // 2)


def superposed():
    fig, axs = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    for j, reg in enumerate(("SVI", "WA")):
        K, days, ep = GT.load_region(reg)
        TT = pd.read_csv(f"{T}/inst_trials_{reg}.csv")
        S = TT[(TT.status == "between") & (TT.sigma_kpa >= 1)]
        pool = days[~days.in_pad]
        obs = {k: np.zeros(len(BINS) - 1) for k in ("data", "c1", "cand_data", "cand_c1")}
        rnd = {k: np.zeros((200, len(BINS) - 1)) for k in obs}
        for _, r in S.iterrows():
            for k in obs:
                obs[k] += np.histogram(K[k] - r.inst_s, BINS)[0]
            cand = pool.dns[pool.grp == r.grp].to_numpy()
            tod = r.inst_s % DAY
            pick = RNG.choice(cand, 200) + tod
            for k in obs:
                a = K[k]
                lo = np.searchsorted(a, pick[:, None] + BINS[None, :])
                rnd[k] += np.diff(lo, axis=1)
        x = (BINS[:-1] + H / 4) / H
        for i, (ka, kb, lab) in enumerate((("data", "c1", "counted-once LFEs (P rule)"),
                                           ("cand_data", "cand_c1", "all candidates, any P (artefact gauge)"))):
            ax = axs[i, j]
            for k, c, nm in ((ka, "C0", "data"), (kb, "C3", "scrambled control")):
                m = rnd[k].mean(0); lo, hi = np.percentile(rnd[k], [1, 99], axis=0)
                ax.fill_between(x, lo, hi, color=c, alpha=0.15, step="mid")
                ax.plot(x, m, color=c, lw=0.8, ls="--")
                ax.step(x, obs[k], where="mid", color=c, lw=1.2, label=f"{nm} (EQ-aligned)")
            ax.axvline(0, color="k", lw=0.8)
            ax.set_title(f"{reg}: {lab}\n{len(S)} strong (≥1 kPa) between-episode earthquakes", fontsize=8)
            ax.set_ylabel("count per 30 min (stacked)")
            if i == 1:
                ax.set_xlabel("hours after P arrival at the region")
            ax.legend(fontsize=7)
        axs[0, j].text(0.01, 1.0, "dashed = random same-clock days, band = 1–99 %", transform=axs[0, j].transAxes,
                       fontsize=7, va="bottom")
    fig.tight_layout(); fig.savefig(f"{F}/g_superposed.png", dpi=120); plt.close(fig)


def dose():
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    for reg, c in (("SVI", "C0"), ("WA", "C1")):
        TT = pd.read_csv(f"{T}/inst_trials_{reg}.csv")
        B = TT[TT.status == "between"]
        y = (B.inst_data + 0.5) / (B.mu_data + 0.5)
        axs[0].scatter(B.sigma_kpa, y, s=10, color=c, alpha=0.6, label=f"{reg} between ({len(B)})")
        E_ = TT[TT.status == "ets"]
        axs[1].scatter(E_.sigma_kpa, (E_.inst_data + 0.5) / (E_.mu_data + 0.5), s=10, color=c, alpha=0.6,
                       label=f"{reg} during ETS ({len(E_)})")
        D = pd.read_csv(f"{T}/delayed_trials_{reg}_all_10.csv")
        axs[2].scatter(D.sigma_kpa, D.after_exc * 1000, s=10, color=c, alpha=0.6, label=f"{reg} ({len(D)})")
    for a in axs:
        a.set_xscale("log"); a.set_xlabel("peak dynamic stress at the region (kPa)"); a.legend(fontsize=7)
    axs[0].set_yscale("log"); axs[1].set_yscale("log")
    axs[0].set_ylabel("(N + 0.5) / (random-window mean + 0.5)"); axs[0].set_title("instantaneous, between episodes", fontsize=9)
    axs[1].set_ylabel("(N + 0.5) / (local ±24 h mean + 0.5)"); axs[1].set_title("instantaneous, during episodes", fontsize=9)
    axs[2].set_ylabel("LFEs per 1000 family-days, days +1..+10\n(data − scrambled control)")
    axs[2].set_title("delayed, between episodes", fontsize=9)
    for a in axs[:2]:
        a.axhline(1, color="k", lw=0.6)
    fig.tight_layout(); fig.savefig(f"{F}/g_dose.png", dpi=120); plt.close(fig)


def ets_fig():
    E = pd.read_csv(f"{T}/eq_list.csv", parse_dates=["t0_utc"])
    fig, axs = plt.subplots(2, 1, figsize=(13, 5.5), sharex=True)
    for ax, reg, f in ((axs[0], "SVI", "svi_episodes.csv"), (axs[1], "WA", "wa_episodes.csv")):
        ep = pd.read_csv(f"{T}/{f}", parse_dates=["start"])
        for s in ep.start[ep.major.astype(bool)]:
            ax.axvline(s, color="0.6", lw=0.8)
        ax.scatter(E.t0_utc, E[f"{reg}_sigma_kpa"], s=8, c=np.where(E[f"{reg}_sigma_kpa"] >= 1, "C3", "C0"))
        ax.set_yscale("log"); ax.set_ylabel(f"{reg} stress (kPa)")
        ax.set_title(f"{reg}: grey lines = major ETS onsets; dots = earthquakes (red ≥ 1 kPa)", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{F}/g_ets_onsets.png", dpi=120); plt.close(fig)


def calib():
    E = pd.read_csv(f"{T}/eq_list.csv")
    ok = E.pgc_ok.fillna(False).astype(bool)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(E.SVI_sig_est_kpa[ok], E.sig_pgc_kpa[ok], s=8, c=np.where(E.SVI_ddeg[ok] < 20, "C1", "C0"))
    ax.plot([1e-3, 1e2], [1e-3, 1e2], "k", lw=0.6)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("estimated (Ms relation) kPa"); ax.set_ylabel("measured at PGC kPa")
    r = np.log10(E.sig_pgc_kpa[ok] / E.SVI_sig_est_kpa[ok])
    ax.set_title(f"n={ok.sum()}, median log10 ratio {r.median():+.2f}, sd {r.std():.2f}\norange = Δ < 20°", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{F}/g_sigma_calib.png", dpi=120); plt.close(fig)
    return dict(n=int(ok.sum()), med=float(r.median()), sd=float(r.std()),
                med_tele=float(r[E.SVI_ddeg[ok] >= 20].median()), med_reg=float(r[E.SVI_ddeg[ok] < 20].median()))


if __name__ == "__main__":
    print(calib())
    dose(); ets_fig(); superposed()
