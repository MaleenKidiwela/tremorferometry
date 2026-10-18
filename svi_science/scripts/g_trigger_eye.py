#!/usr/bin/env python
"""g_trigger_eye (2026-10-11) -- Study G eye check (PLAN.md §6.3): 2-8 Hz record sections + LFE timeline per trial.
20 Hz caches (bw_p4env_detect.cache_path) for the envelopes, raw CN.PGC for the 0.02-2 Hz broadband trace.
Single process, light.

    python svi_science/scripts/g_trigger_eye.py REGION:TRIAL[:label] ...   -> figures/g_eye_<region>_<trial>.png
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
from obspy import UTCDateTime, read_inventory
import bw_p4env_detect as E
import g_trigger_pgv as PG

ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger"
FS = 20.0
SOS = butter(4, [2 / (FS / 2), 8 / (FS / 2)], "band", output="sos")
NS = 10 ** 9
C = E.coords()
REF = {"SVI": (48.4837, -123.7485), "WA": (47.8432, -123.0282)}


def env_day(sta, d):
    p = E.cache_path(sta, d)
    if p is None:
        return None
    z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"])
    day0 = pd.Timestamp(d).value / 1e9
    n = int(86400 * FS)
    X = np.zeros((3 if x.ndim == 2 else 1, n)); x = np.atleast_2d(x)
    i0 = int(round((t0 - day0) * FS)); a, b = max(i0, 0), min(i0 + x.shape[1], n)
    if b <= a:
        return None
    X[:x.shape[0], a:b] = x[:, a - i0:b - i0]
    X[~np.isfinite(X)] = 0
    gap = ~np.any(X != 0, axis=0)
    Y = sosfiltfilt(SOS, X, axis=1)
    e = np.sqrt((Y ** 2).sum(0))
    e[gap] = np.nan
    k = int(10 * FS)
    m = len(e) // k
    return np.sqrt(np.nanmean(e[:m * k].reshape(m, k) ** 2, 1))       # 10-s RMS, 8640 values


def env_span(sta, t1, t2):
    out = []
    d = pd.Timestamp(t1 // NS * NS).normalize()
    while d.value < t2:
        e = env_day(sta, d.date())
        tt = d.value + np.arange(8640) * 10 * NS
        out.append((tt, e if e is not None else np.full(8640, np.nan)))
        d += pd.Timedelta(days=1)
    tt = np.concatenate([o[0] for o in out]); ee = np.concatenate([o[1] for o in out])
    m = (tt >= t1) & (tt < t2)
    return tt[m], ee[m]


def main(specs):
    PG.INV = read_inventory(f"{ROOT}/codex-LFE/results/full_history_metadata_v1/CN.xml").select(station="PGC")
    fams = {"SVI": pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv").set_index("fid").along_strike_km,
            "WA": pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv").set_index("fid").along_km}
    ev = {"SVI": pd.read_parquet(f"{G}/tables/svi_events.parquet").rename(columns={"fam": "fid"}),
          "WA": pd.read_parquet(f"{G}/tables/wa_events.parquet")}
    trem = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    trem["tns"] = pd.to_datetime(trem.time).astype("datetime64[ns]").astype("int64")
    for spec in specs:
        reg, trial = spec.split(":")[:2]; lab = spec.split(":")[2] if spec.count(":") >= 2 else ""
        r = "WA" if reg.startswith("WA") else "SVI"
        TT = pd.read_csv(f"{G}/tables/inst_trials_{reg}.csv").set_index("trial")
        row = TT.loc[int(trial)]
        t1 = int(row.inst_s - 45 * 60 * NS); t2 = int(min(max(row.inst_e, row.inst_s + 3 * 3600 * NS) + 60 * 60 * NS,
                                                         row.inst_s + 8 * 3600 * NS))
        tm = lambda t: (np.asarray(t) - row.t0_ns) / 60e9          # minutes after origin
        sta = sorted([s for s in C if E.cache_path(s, pd.Timestamp(row.inst_s).date()) is not None],
                     key=lambda s: -C[s][0])
        fig, ax = plt.subplots(3, 1, figsize=(11, 3 + 0.42 * len(sta) + 4), sharex=True,
                               gridspec_kw=dict(height_ratios=[1.3, max(3, 0.42 * len(sta)), 2.2]))
        # broadband
        st = PG.load(UTCDateTime(t1 / 1e9), UTCDateTime(t2 / 1e9))
        for tr in st:
            tr.stats.sampling_rate = round(tr.stats.sampling_rate, 3)
        try:
            st.merge(method=1, fill_value=0)
        except Exception:
            st = st.select(channel="?H?")[:0]
        for tr in st:
            if tr.stats.channel[:2] not in ("BH", "HH"):
                continue
            tr.detrend("demean"); tr.taper(0.02)
            if tr.stats.sampling_rate > 20:
                tr.filter("lowpass", freq=8.0, zerophase=True); tr.decimate(int(tr.stats.sampling_rate // 20), no_filter=True)
            tr.data = tr.data / PG.INV.get_response(tr.id, tr.stats.starttime).instrument_sensitivity.value
            tr.filter("bandpass", freqmin=0.02, freqmax=2, zerophase=True)
            tt = (tr.stats.starttime.ns + np.arange(tr.stats.npts) * tr.stats.delta * NS)
            ax[0].plot(tm(tt), tr.data * 1e6, lw=0.4, label=tr.stats.channel)
        ax[0].set_ylabel("PGC µm/s\n0.02–2 Hz"); ax[0].legend(loc="upper right", fontsize=7, ncol=3)
        # record section
        for i, s in enumerate(sta):
            tt, ee = env_span(s, t1, t2)
            pre = ee[tt < row.inst_s]
            ref = np.nanmedian(pre) if np.isfinite(pre).any() else np.nanmedian(ee)
            y = np.log10(ee / ref)
            ax[1].plot(tm(tt), i + 0.6 * np.clip(y, -0.5, 1.5), lw=0.6, color="k")
        ax[1].set_yticks(range(len(sta))); ax[1].set_yticklabels(sta, fontsize=7)
        ax[1].set_ylabel("2–8 Hz envelope\nlog10(/pre-median), N→S")
        # LFEs
        X = ev[r]; X = X[(X.time_utc.astype("int64") >= t1) & (X.time_utc.astype("int64") < t2)]
        for ser, mk, c in (("data", "o", "C0"), ("c1", "x", "C3")):
            q = X[X.series == ser]
            y = q.fid.map(fams[r]).to_numpy()
            ax[2].scatter(tm(q.time_utc.astype("int64")), y, s=14 if ser == "data" else 30, marker=mk, color=c,
                          label=f"{ser} ({len(q)})" if ser == "data" else f"scrambled control ({len(q)})")
            if ser == "data":
                mz = q.n_z2 >= 2
                ax[2].scatter(tm(q.time_utc.astype("int64"))[mz.to_numpy()], y[mz.to_numpy()], s=40, facecolors="none",
                              edgecolors="C2", label=f"n_z2 ≥ 2 ({int(mz.sum())})")
        la, lo = REF[r]
        dk = 111.195 * np.hypot(trem.lat - la, (trem.lon - lo) * np.cos(np.radians(la)))
        tq = trem.tns[(dk <= 100) & (trem.tns >= t1) & (trem.tns < t2)]
        if len(tq):
            ax[2].plot(tm(tq), np.full(len(tq), ax[2].get_ylim()[0] if len(X) else -60), "|", color="0.5",
                       label=f"PNSN tremor ≤100 km ({len(tq)})")
        ax[2].set_ylabel("family along-strike km"); ax[2].legend(fontsize=7, loc="upper right")
        for a in ax:
            for t, ls, c in ((row.inst_s, "-", "C1"), (row.surf_s, "--", "C4"), (row.inst_e, "-", "C1")):
                a.axvline(tm(t), ls=ls, color=c, lw=0.8)
            a.axvline(tm(row.t0_ns + row.dkm / 2.5 * NS * 1.0), ls=":", color="C4", lw=0.8)
        ax[2].set_xlabel("minutes after origin time  (orange = W_inst start/end, purple dashed = 4.5 km/s, dotted = 2.5 km/s)")
        fig.suptitle(f"{reg} trial {trial}: {row.place}, M{row.mag}, {row.t0_PT}; Δ {row.dkm:.0f} km, σ {row.sigma_kpa:.2f} kPa "
                     f"({row.sigma_src}); status {row.status}\nW_inst LFEs data {row.inst_data}, scrambled {row.inst_c1}, "
                     f"random-window mean {row.get('mu_data', np.nan):.3f}  {lab}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{G}/figures/g_eye_{reg}_{trial}.png", dpi=110); plt.close(fig)
        print("wrote", reg, trial, len(sta), "stations", len(X), "LFE rows")


if __name__ == "__main__":
    main(sys.argv[1:])
