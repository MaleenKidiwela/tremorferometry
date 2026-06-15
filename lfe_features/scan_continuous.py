"""Deploy the per-window picker v0 as a CONTINUOUS SCANNER.

Slides 40 s / 3-comp windows (step --step s) across a day, computes features, and applies the
saved 3-class model -> P(LFE) and P(EQ) traces for the whole day. Auto-picks a busy-LFE
(ETS) day and a quiet day for contrast; overlays Lin LFE times and ANSS EQ times; quantifies
correlation of the P(LFE) trace with the hourly Lin-LFE rate.

Usage: PYTHONPATH=src python lfe_features/scan_continuous.py --net PB --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import joblib
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all
HERE = os.path.dirname(__file__); PRE, POST = 10.0, 30.0


def scan_day(net, sta, year, julday, step, cols):
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    st = read(path)

    def get(*ch):
        for c in ch:
            s = st.select(channel=f"*{c}")
            if len(s):
                s = s.merge(method=1, fill_value=0); return s[0]
        return None
    trZ, trH1, trH2 = get("Z"), get("1", "N"), get("2", "E")
    fs = float(trZ.stats.sampling_rate)
    Z, H1, H2 = trZ.data.astype(float), trH1.data.astype(float), trH2.data.astype(float)
    L = min(len(Z), len(H1), len(H2)); W = int((PRE + POST) * fs); st0 = int(step * fs)
    day0 = UTCDateTime(year=year, julday=julday).timestamp
    feats = []; times = []
    for i in range(0, L - W, st0):
        z = Z[i:i + W]
        if np.std(z) == 0 or not np.all(np.isfinite(z)):
            continue
        f = compute_all(z, fs, PRE, H1[i:i + W], H2[i:i + W])
        if any(c not in f or not np.isfinite(f[c]) for c in cols):
            continue
        feats.append([f[c] for c in cols]); times.append(day0 + (i + W * (PRE / (PRE + POST))) / fs)
    return np.array(feats), np.array(times)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--step", type=float, default=5.0)
    a = ap.parse_args(); s = a.sta.lower()
    B = joblib.load(f"{HERE}/models/picker_v0_{s}.joblib")
    model, scaler, cols, classes = B["model"], B["scaler"], B["cols"], B["classes"]
    iL = classes.index("LFE"); iE = classes.index("EQ")

    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0]))
               for f in glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed"))
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dk = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dk < 35)].copy()
    cur["t"] = pd.to_datetime(cur.OT, utc=True); cur["y"] = cur.t.dt.year; cur["j"] = cur.t.dt.dayofyear
    cur = cur[[ (r.y, r.j) in have for r in cur.itertuples(index=False)]]
    daycount = cur.groupby(["y", "j"]).size().sort_values(ascending=False)
    busy = daycount.index[0]
    quiet = [k for k in have if k not in set(daycount.index)]
    quiet = sorted(quiet)[len(quiet) // 2] if quiet else daycount.index[-1]
    ev = pd.read_csv(f"data/eq_catalog_{s}.csv"); ev["t"] = pd.to_datetime(ev.ot, utc=True)

    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=False)
    for ax, (yy, jj), tag in [(axes[0], busy, "BUSY (ETS)"), (axes[1], quiet, "QUIET")]:
        print(f"[{a.sta}] scanning {tag} day {yy}-{jj:03d}...", flush=True)
        X, T = scan_day(a.net, a.sta, yy, jj, a.step, cols)
        if len(X) == 0:
            continue
        P = model.predict_proba(scaler.transform(X))
        hrs = (T - T.min()) / 3600.0
        ax.plot(hrs, P[:, iL], color="#1a9850", lw=0.8, label="P(LFE)")
        ax.plot(hrs, P[:, iE], color="#7570b3", lw=0.6, alpha=0.7, label="P(EQ)")
        from obspy import UTCDateTime
        d0 = UTCDateTime(year=yy, julday=jj).timestamp
        lin_d = cur[(cur.y == yy) & (cur.j == jj)]
        for tt in lin_d.t:
            ax.axvline((tt.timestamp() - d0) / 3600.0, color="green", alpha=0.25, lw=0.8)
        eqd = ev[(ev.t.dt.year == yy) & (ev.t.dt.dayofyear == jj)]
        for tt in eqd.t:
            ax.axvline((tt.timestamp() - d0) / 3600.0, color="purple", alpha=0.4, lw=0.8, ls=":")
        # hourly correlation P(LFE) vs Lin-LFE rate
        hb = np.arange(0, 25)
        lin_rate, _ = np.histogram((lin_d.t.values.astype("datetime64[ns]").astype("int64") / 1e9 - d0) / 3600, bins=hb)
        pl_h = np.array([P[(hrs >= h) & (hrs < h + 1), iL].mean() if ((hrs >= h) & (hrs < h + 1)).any() else 0 for h in hb[:-1]])
        r = np.corrcoef(lin_rate, pl_h)[0, 1] if lin_rate.std() > 0 else np.nan
        ax.set_title(f"{a.sta} {tag} {yy}-{jj:03d}: {len(lin_d)} Lin LFEs (green), {len(eqd)} EQ (purple) | "
                     f"hourly corr P(LFE)-vs-LFErate r={r:.2f}", fontsize=9)
        ax.set_ylabel("probability"); ax.set_ylim(0, 1); ax.legend(loc="upper right", fontsize=8)
        print(f"   {tag}: {len(X)} windows, mean P(LFE)={P[:,iL].mean():.3f}, hourly corr r={r:.2f}", flush=True)
    axes[1].set_xlabel("hour of day (UTC)")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/scan_{s}.png", dpi=120)
    print(f"saved figures/scan_{s}.png")
    print("SCAN DONE")


if __name__ == "__main__":
    main()
