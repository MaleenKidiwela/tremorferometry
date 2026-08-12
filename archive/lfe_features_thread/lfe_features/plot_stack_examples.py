"""Visual verdict: plot coherent stacks (waveform + amplitude spectrum) for example
B011 families (GOLD + FAIL) next to LIN-positive and RAND-null stacks. Do GOLD stacks
show emergent, low-frequency, S-dominant LFE character?

Usage: PYTHONPATH=src python lfe_features/plot_stack_examples.py --net PB --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import signal
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

PRE, POST, ANCHOR = 10.0, 30.0, 11.0
HERE = os.path.dirname(__file__)
POL = ["hv_ratio", "rectilinearity", "planarity"]


def cut_day(payload):
    net, sta, year, julday, items = payload
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return {}
    try:
        st = read(path)
    except Exception:
        return {}
    s = st.select(channel="*Z")
    if len(s) == 0:
        return {}
    try: s = s.merge(method=1, fill_value=0)
    except Exception: pass
    trZ = s[0]; fs = float(trZ.stats.sampling_rate)
    exp = int(round((PRE + POST) * fs))
    acc = {}
    for key, anch in items:
        try:
            z = trZ.slice(UTCDateTime(anch) - PRE, UTCDateTime(anch) + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        z = z[:exp] / (np.std(z[:exp]) + 1e-20)
        if key not in acc:
            acc[key] = [np.zeros(exp), 0]
        acc[key][0] += z; acc[key][1] += 1
    return acc


def gather(net, sta, jobs, workers=12):
    payloads = [(net, sta, k[0], k[1], v) for k, v in jobs.items()]
    tot = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            for key, val in f.result().items():
                if key not in tot:
                    tot[key] = val
                else:
                    tot[key][0] += val[0]; tot[key][1] += val[1]
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--per", type=int, default=300); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)
    fs = 100.0 if a.net == "PB" else 40.0

    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files)

    A = pd.read_csv("results/family_verdicts/ALL_families_labeled.csv")
    A = A[A.station == a.sta]
    gold = A[A.grade == "GOLD"].nlargest(4, "coda_sigma").fam.tolist()
    fail = A[A.grade == "FAIL"].fam.tolist()[:2]
    pick = [("GOLD", f) for f in gold] + [("FAIL", f) for f in fail]

    m = pd.read_csv(f"data/mf_{s}_all.csv", usecols=["template", "time"])
    m["t"] = pd.to_datetime(m.time, utc=True); m["year"] = m.t.dt.year; m["julday"] = m.t.dt.dayofyear
    m = m[m.template.isin([f for _, f in pick])]
    m = m[[ (r.year, r.julday) in have for r in m.itertuples(index=False)]]
    jobs = defaultdict(list)
    for fam, g in m.groupby("template"):
        for r in g.sample(min(a.per, len(g)), random_state=a.seed).itertuples(index=False):
            jobs[(int(r.year), int(r.julday))].append(((("FAM", fam)), r.t.value / 1e9))

    # LIN + RAND control windows
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dkm = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dkm < 35)].copy()
    cur["t"] = pd.to_datetime(cur.OT, utc=True); cur["year"] = cur.t.dt.year; cur["julday"] = cur.t.dt.dayofyear
    cur = cur[[ (r.year, r.julday) in have for r in cur.itertuples(index=False)]].sample(min(a.per, len(cur)), random_state=1)
    for r in cur.itertuples(index=False):
        jobs[(int(r.year), int(r.julday))].append((("CTRL", "LIN"), r.t.value / 1e9 + ANCHOR))
    daylist = sorted(have)
    from datetime import datetime
    for j in range(a.per):
        y, jd = daylist[rng.randint(len(daylist))]
        base = datetime(y, 1, 1).timestamp() + (jd - 1) * 86400 + rng.uniform(60, 86340)
        jobs[(y, jd)].append((("CTRL", "RAND"), base + ANCHOR))

    print(f"[{a.sta}] gathering stacks for {len(pick)} families + LIN/RAND over {len(jobs)} days...", flush=True)
    tot = gather(a.net, a.sta, jobs)

    # train P(LFE) model
    ref = pd.read_parquet(f"{HERE}/data/feat_ref_{s}.parquet")
    vert = [c for c in FEATURE_ORDER if c in ref.columns and c not in POL]
    cols = vert
    ref = ref.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)
    sub = ref[ref.label.isin(["LIN", "RAND"])]
    yb = (sub.label == "LIN").astype(int).values
    sc = StandardScaler().fit(sub[cols]); rf = RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0).fit(sc.transform(sub[cols]), yb)

    def score(z):
        f = compute_all(z, fs, PRE)
        row = {c: f.get(c, np.nan) for c in cols}
        if any(not np.isfinite(v) for v in row.values()):
            return np.nan
        return float(rf.predict_proba(sc.transform(pd.DataFrame([row])[cols]))[:, 1][0])

    panels = [(("CTRL", "LIN"), "LIN (real LFEs)", "green")] + \
             [((("FAM", f)), f"GOLD {f.split('__')[-1]}", "#1a9850") for _, f in pick if _ == "GOLD"] + \
             [((("FAM", f)), f"FAIL {f.split('__')[-1]}", "#d73027") for _, f in pick if _ == "FAIL"] + \
             [(("CTRL", "RAND"), "RAND (noise)", "gray")]
    t = np.arange(int((PRE + POST) * fs)) / fs - PRE
    fig, axes = plt.subplots(len(panels), 2, figsize=(12, 2.0 * len(panels)))
    for i, (key, title, c) in enumerate(panels):
        if key not in tot or tot[key][1] < 10:
            continue
        z = tot[key][0] / tot[key][1]
        p = score(z)
        ax = axes[i, 0]
        ax.plot(t, z, c=c, lw=0.6)
        ax.set_ylabel(f"{title}\nN={tot[key][1]}", fontsize=8)
        ax.axvline(0, color="k", ls=":", lw=0.5)
        ax.text(0.98, 0.9, f"P(LFE)={p:.2f}", transform=ax.transAxes, ha="right", fontsize=9,
                bbox=dict(fc="white", ec=c))
        if i == 0: ax.set_title("coherent stack — vertical waveform (S~0 s)")
        zb = signal.detrend(z)
        f_, P = signal.welch(zb, fs=fs, nperseg=min(len(zb), 256))
        ax2 = axes[i, 1]
        ax2.semilogy(f_, P, c=c, lw=0.8); ax2.set_xlim(0, min(40, fs / 2))
        if i == 0: ax2.set_title("amplitude spectrum")
        ax2.axvspan(2, 8, color="yellow", alpha=0.15)
    axes[-1, 0].set_xlabel("time rel. anchor (s)"); axes[-1, 1].set_xlabel("Hz")
    fig.suptitle(f"{a.sta}: family stacks vs Lin LFEs vs noise (yellow = 2-8 Hz LFE band)")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = f"{HERE}/figures/stack_examples_{s}.png"
    fig.savefig(out, dpi=120)
    print(f"saved {out}")
    print("PLOT DONE")


if __name__ == "__main__":
    main()
