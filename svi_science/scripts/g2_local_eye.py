#!/usr/bin/env python
"""g2_local_eye (2026-10-11) -- Study G2 eye check (PLAN.md §5): 2-8 Hz record sections around a LOCAL earthquake,
PGC broadband, LFE rows (data / scrambled / n_z2 >= 2), PNSN tremor.  Reuses g_trigger_eye helpers.

    python svi_science/scripts/g2_local_eye.py REGION:TRIAL[:label] ...  -> results/g2_local/figures/g2_eye_<region>_<trial>.png
    (REGION = SVI | WA | WA_pre; TRIAL = row 'trial' of tables/q1_trials_<REGION>.parquet)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import UTCDateTime, read_inventory
import g_trigger_eye as GE

ROOT = "/home/jovyan/tremorferometry"
G2 = f"{ROOT}/svi_science/results/g2_local"
GTT = f"{ROOT}/svi_science/results/g_trigger/tables"
NS = 10 ** 9


def main(specs):
    PG = GE.PG
    PG.INV = read_inventory(f"{ROOT}/codex-LFE/results/full_history_metadata_v1/CN.xml").select(station="PGC")
    fams = {"SVI": pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv").set_index("fid").along_strike_km,
            "WA": pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv").set_index("fid").along_km}
    ev = {"SVI": pd.read_parquet(f"{GTT}/svi_events.parquet").rename(columns={"fam": "fid"}),
          "WA": pd.read_parquet(f"{GTT}/wa_events.parquet")}
    trem = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    trem["tns"] = pd.to_datetime(trem.time).astype("datetime64[ns]").astype("int64")
    EQ = pd.read_parquet(f"{G2}/tables/eq_local.parquet")
    EQ = EQ[EQ.in_analysis].assign(tns=lambda d: d.t.astype("datetime64[ns]").astype("int64"))
    for spec in specs:
        parts = spec.split(":"); reg, trial = parts[0], int(parts[1]); lab = parts[2] if len(parts) > 2 else ""
        r = "WA" if reg.startswith("WA") else "SVI"
        TT = pd.read_parquet(f"{G2}/tables/q1_trials_{reg}.parquet").set_index("trial")
        row = TT.loc[trial]
        t0 = int(row.t0)
        t1 = t0 - 60 * 60 * NS; t2 = t0 + 7 * 3600 * NS
        tm = lambda t: (np.asarray(t, dtype=np.int64) - t0) / 60e9
        sta = sorted([s for s in GE.C if GE.E.cache_path(s, pd.Timestamp(t0).date()) is not None], key=lambda s: -GE.C[s][0])
        fig, ax = plt.subplots(4, 1, figsize=(11, 5 + 0.42 * len(sta) + 5), sharex=True,
                               gridspec_kw=dict(height_ratios=[1.2, max(3, 0.42 * len(sta)), 2.2, 1.0]))
        try:
            st = PG.load(UTCDateTime(t1 / 1e9), UTCDateTime(t2 / 1e9))
            for tr in st:
                tr.stats.sampling_rate = round(tr.stats.sampling_rate, 3)
            st.merge(method=1, fill_value=0)
            for tr in st:
                if tr.stats.channel[:2] not in ("BH", "HH"):
                    continue
                tr.detrend("demean"); tr.taper(0.02)
                if tr.stats.sampling_rate > 20:
                    tr.filter("lowpass", freq=8.0, zerophase=True); tr.decimate(int(tr.stats.sampling_rate // 20), no_filter=True)
                tr.data = tr.data / PG.INV.get_response(tr.id, tr.stats.starttime).instrument_sensitivity.value
                tr.filter("bandpass", freqmin=2, freqmax=8, zerophase=True)
                tt = tr.stats.starttime.ns + np.arange(tr.stats.npts) * tr.stats.delta * NS
                ax[0].plot(tm(tt), tr.data * 1e6, lw=0.3, label=tr.stats.channel)
            ax[0].legend(loc="upper right", fontsize=7, ncol=3)
        except Exception as e:
            ax[0].text(0.5, 0.5, f"PGC broadband not available ({type(e).__name__})", transform=ax[0].transAxes, ha="center")
        ax[0].set_ylabel("PGC µm/s\n2–8 Hz")
        for i, s in enumerate(sta):
            tt, ee = GE.env_span(s, t1, t2)
            pre = ee[tt < t0]
            ref = np.nanmedian(pre) if np.isfinite(pre).any() else np.nanmedian(ee)
            ax[1].plot(tm(tt), i + 0.6 * np.clip(np.log10(ee / ref), -0.5, 1.5), lw=0.6, color="k")
        ax[1].set_yticks(range(len(sta))); ax[1].set_yticklabels(sta, fontsize=7)
        ax[1].set_ylabel("2–8 Hz envelope (10-s RMS)\nlog10(/pre-median), N→S")
        X = ev[r]; X = X[(X.time_utc.astype("int64") >= t1) & (X.time_utc.astype("int64") < t2)]
        for ser, mk, c in (("data", "o", "C0"), ("c1", "x", "C3")):
            q = X[X.series == ser]
            y = q.fid.map(fams[r]).to_numpy()
            ax[2].scatter(tm(q.time_utc.astype("int64")), y, s=14 if ser == "data" else 30, marker=mk, color=c,
                          label=f"LFE data ({len(q)})" if ser == "data" else f"scrambled control ({len(q)})")
            if ser == "data":
                mz = (q.n_z2 >= 2).fillna(False).to_numpy()
                ax[2].scatter(tm(q.time_utc.astype("int64"))[mz], y[mz], s=40, facecolors="none", edgecolors="C2",
                              label=f"n_z2 ≥ 2 ({int(mz.sum())})")
        la, lo = GE.REF[r]
        dk = 111.195 * np.hypot(trem.lat - la, (trem.lon - lo) * np.cos(np.radians(la)))
        tq = trem.tns[(dk <= 100) & (trem.tns >= t1) & (trem.tns < t2)]
        if len(tq):
            ax[2].plot(tm(tq), np.full(len(tq), np.nanmin(list(fams[r].values)) - 5), "|", color="0.5",
                       label=f"PNSN tremor ≤100 km ({len(tq)})")
        ax[2].set_ylabel("family along-strike km"); ax[2].legend(fontsize=7, loc="upper right")
        # local earthquakes (catalogue) in the span, within 150 km of the region
        q = EQ[(EQ.tns >= t1) & (EQ.tns < t2) & (EQ[f"{r}_dist_km"] <= 150)]
        ax[3].scatter(tm(q.tns), q.mag, s=12, c="C1")
        ax[3].set_ylabel("local EQs\n≤150 km (M)")
        for a in ax:
            a.axvline(0, color="C1", lw=0.8); a.axvline(360, color="C1", lw=0.8, ls="--")
        ax[3].set_xlabel("minutes after origin (orange: origin, dashed: end of the 6-h window)")
        fig.suptitle(f"{reg} trial {trial}: M{row.mag:.1f} {row.cls}, depth {row.depth:.0f} km, {row.place}, {row.t_PT}; "
                     f"{row.dist_km:.0f} km from nearest family; status {row.status}\n"
                     f"6-h window: LFEs data {row.H6_data}, scrambled {row.H6_c1}, random mean {row.get('mu_data', np.nan):.3f}; "
                     f"Δσ_d ~{row.dyn_kpa:.2g} kPa  {lab}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{G2}/figures/g2_eye_{reg}_{trial}.png", dpi=105); plt.close(fig)
        print("wrote", reg, trial, len(sta), "stations", len(X), "LFE rows", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
