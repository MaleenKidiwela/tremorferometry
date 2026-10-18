#!/usr/bin/env python
"""Track A-2 step 5(b) — STATION CORRECTIONS from intraslab EARTHQUAKES, not from Bostock's families.

Merlin F10(ii): the far boreholes and LRIV carry systematic S residuals of +1.7 to +3.25 s against the
Crosson 1-D model; taking those corrections from Bostock's own hypocentres and then "recovering" his
hypocentres is circular.  So they are measured here on catalogued intraslab earthquakes (USGS, M >= 2,
25-50 km depth, inside the family box 48.0-49.0 N / -124.5..-123.0 E) picked on OUR OWN stations, with
the USGS hypocentre and origin time as the reference.

Picks: P on Z, S on sqrt(H1^2+H2^2), 2-15 Hz at the NATIVE sample rate (the 20 Hz / 1-8 Hz LFE cache is
too narrow for a M2 P onset), envelope-onset + AIC refinement in a +-3 s (P) / +-4 s (S) window around
the Crosson prediction, accepted at peak/noise-rms >= --min-snr.  Residual = observed - predicted; the
station correction is the MEDIAN residual over events (reported with its MAD and n).

Outputs: results/a5_station_corrections.csv, results/a5_eq_picks.csv, figures/a5_eq_picks.png
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d3_vel1d as V                                          # noqa: E402
import a5_locate as A                                         # noqa: E402
import e4_chanorder as chanorder                              # noqa: E402

EXT = os.path.dirname(HERE)
REP = A.REP
TRACE_ROOTS = [f"{EXT}/traces", f"{REP}/traces"]
STAS = ["PGC", "VGZ", "SNB", "NLLB", "LRIV", "B926", "B001", "B003", "B004",
        "B005", "B006", "B007", "B009", "B010", "B011"]


def smooth(x, n):
    return np.convolve(x, np.ones(n) / n, mode="same") if n >= 2 else x


def aic_pick(x, i0, i1):
    seg = np.asarray(x[i0:i1], float)
    n = seg.size
    if n < 12:
        return (i0 + i1) // 2
    aic = np.full(n, np.inf)
    for k in range(5, n - 5):
        va, vb = seg[:k].var(), seg[k:].var()
        if va > 0 and vb > 0:
            aic[k] = k * np.log(va) + (n - k) * np.log(vb)
    return i0 + (int(np.argmin(aic)) if np.isfinite(aic).any() else n // 2)


def env_onset(x, lo, hi, nnoise, fs):
    """Envelope-onset pick inside [lo, hi) -> (index, peak/noise-rms).  d3_stage4_locate recipe."""
    from scipy.signal import hilbert
    x = np.asarray(x, float)
    env = smooth(np.abs(hilbert(x - x.mean())), max(2, int(round(0.25 * fs))))
    noise = float(np.median(env[:nnoise])) if nnoise > 2 else 0.0
    nrms = float(np.sqrt(np.mean(x[:nnoise] ** 2))) if nnoise > 2 else 0.0
    lo, hi = max(lo, nnoise), min(hi, len(x) - 2)
    if hi <= lo or nrms <= 0:
        return None, 0.0
    j = lo + int(np.argmax(env[lo:hi]))
    snr = float(env[j] / nrms)
    lvl = noise + 0.25 * (env[j] - noise)
    k = j
    while k > max(lo - int(2.0 * fs), nnoise) and env[k] > lvl:
        k -= 1
    w = int(round(0.75 * fs))
    return aic_pick(x, max(k - w, nnoise), min(k + w, len(x) - 1)), snr


def read_event(netsta, t0, pre, post):
    """3-C (Z, N|1, E|2) native-rate window [t0-pre, t0+post] or None."""
    from obspy import read, UTCDateTime
    t = UTCDateTime(t0)
    out = []
    for tt in (t - pre, t + post):
        pass
    for root in TRACE_ROOTS:
        p = f"{root}/{netsta}/{t.year}/{t.julday:03d}.mseed"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                st = read(p, starttime=t - pre - 30, endtime=t + post + 30)
            except Exception:
                continue
            if len(st) == 0:
                continue
            for tr in st:
                r = float(tr.stats.sampling_rate)
                if r > 0 and abs(r - round(r)) / r < 1e-4 and r != round(r):
                    tr.stats.sampling_rate = float(round(r))
            st.merge(method=1, fill_value=0)
            try:
                order = chanorder.resolve_order({tr.stats.channel for tr in st}, station=netsta)
            except Exception:
                continue
            st.trim(t - pre, t + post, pad=True, fill_value=0)
            try:
                X = np.vstack([st.select(channel=c)[0].data.astype(float) for c in order])
            except Exception:
                continue
            return X, float(st[0].stats.sampling_rate), float(st.select(channel=order[0])[0].stats.starttime.timestamp)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cat", default=f"{EXT}/results/a5_eq_catalog_raw.csv")
    ap.add_argument("--min-snr", type=float, default=4.0)
    ap.add_argument("--min-events", type=int, default=3)
    ap.add_argument("--fig-event", default="")
    a = ap.parse_args()
    from scipy.signal import butter, sosfiltfilt

    stab = pd.read_csv(f"{REP}/results/e5r_stations.csv").set_index("sta")
    E = pd.read_csv(a.cat)
    E["t"] = pd.to_datetime(E.time)
    E = E[E.n.fillna(0) >= 5].reset_index(drop=True) if "n" in E.columns else E
    print(f"[eq] {len(E)} catalogued intraslab earthquakes with >= 5 stations on disk", flush=True)

    rows = []
    for ie, e in enumerate(E.itertuples()):
        tor = e.t.timestamp()
        for s in STAS:
            ns = str(stab.loc[s, "netsta"])
            d = A.epi_km(e.lat, e.lon, float(stab.loc[s, "lat"]), float(stab.loc[s, "lon"]))
            ttP, ttS = float(V.tt(d, e.dep, "P")), float(V.tt(d, e.dep, "S"))
            got = read_event(ns, tor + ttP, 20.0, (ttS - ttP) + 40.0)
            if got is None:
                continue
            X, fs, tstart = got
            f2 = min(15.0, 0.45 * fs)
            sos = butter(4, [2.0 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
            X = sosfiltfilt(sos, X, axis=-1)
            iP0 = int(round((tor + ttP - tstart) * fs))
            iS0 = int(round((tor + ttS - tstart) * fs))
            nnoise = max(int(10 * fs), 50)
            iP, snrP = env_onset(X[0], iP0 - int(3 * fs), iP0 + int(3 * fs), nnoise, fs)
            H = np.sqrt(X[1] ** 2 + X[2] ** 2)
            iS, snrS = env_onset(H, iS0 - int(4 * fs), iS0 + int(4 * fs), nnoise, fs)
            rP = (tstart + iP / fs) - (tor + ttP) if iP is not None else np.nan
            rS = (tstart + iS / fs) - (tor + ttS) if iS is not None else np.nan
            rows.append(dict(ie=ie, time=str(e.t), sta=s, dist_km=d, dep=e.dep, mag=e.mag, fs=fs,
                             ttP=ttP, ttS=ttS, snrP=snrP, snrS=snrS,
                             resP=rP if snrP >= a.min_snr else np.nan,
                             resS=rS if snrS >= a.min_snr else np.nan,
                             resP_raw=rP, resS_raw=rS))
        print(f"  [eq {ie}] {e.time} M{e.mag} {e.dep:.0f} km -> "
              f"{sum(1 for r in rows if r['ie'] == ie)} stations read", flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(f"{EXT}/results/a5_eq_picks.csv", index=False)

    out = []
    for s in STAS:
        g = P[P.sta == s]
        out.append(dict(sta=s, n_P=int(g.resP.notna().sum()), n_S=int(g.resS.notna().sum()),
                        corr_P=float(g.resP.median()) if g.resP.notna().sum() else np.nan,
                        corr_S=float(g.resS.median()) if g.resS.notna().sum() else np.nan,
                        mad_P=float((g.resP - g.resP.median()).abs().median()) if g.resP.notna().sum() else np.nan,
                        mad_S=float((g.resS - g.resS.median()).abs().median()) if g.resS.notna().sum() else np.nan))
    C = pd.DataFrame(out)
    C["usable"] = (C.n_S >= a.min_events)
    C.loc[~C.usable, ["corr_P", "corr_S"]] = np.nan
    C.to_csv(f"{EXT}/results/a5_station_corrections.csv", index=False)
    print(C.to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for i, ph in enumerate(("P", "S")):
        for s in STAS:
            g = P[(P.sta == s)]
            ax[i].scatter(g.dist_km, g[f"res{ph}"], s=18, label=s)
        ax[i].axhline(0, color="k", lw=.8)
        ax[i].set_xlabel("epicentral distance (km)"); ax[i].set_ylabel(f"{ph} residual vs Crosson (s)")
        ax[i].set_title(f"{ph}: {int(P[f'res{ph}'].notna().sum())} accepted picks "
                        f"(SNR >= {a.min_snr}) on {len(E)} earthquakes")
        ax[i].grid(alpha=.3)
    ax[1].legend(fontsize=6, ncol=3)
    plt.tight_layout(); plt.savefig(f"{EXT}/figures/a5_eq_picks.png", dpi=130)
    print("[write] results/a5_station_corrections.csv, results/a5_eq_picks.csv, figures/a5_eq_picks.png")


if __name__ == "__main__":
    main()
