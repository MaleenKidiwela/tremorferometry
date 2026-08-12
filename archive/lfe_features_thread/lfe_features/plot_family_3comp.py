"""Plot one LFE family's coherently-stacked signal on all 3 components (Z, EH1, EH2).

Stacks high-cc matched-filter detections (already template-aligned) at the family's detection
times; each detection's 3 channels are normalized by a SINGLE factor (max|amp| over the 3) so
the inter-component amplitude ratio is preserved (LFEs are S-dominant -> bigger on horizontals).

Usage: PYTHONPATH=src python lfe_features/plot_family_3comp.py --net PB --sta B011 \
         --fam 48.900_-124.750__c86 --mf data/mf_b011p70_all.csv --cc-min 0.88 --n 2000
"""
import argparse, os, sys, glob
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(__file__)
PRE, POST = 3.0, 10.0    # 13 s window


def cut_day(payload):
    net, sta, year, julday, times = payload
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return None
    try:
        st = read(path)
    except Exception:
        return None
    def get(*c):
        for ch in c:
            s = st.select(channel=f"*{ch}")
            if len(s):
                try: s = s.merge(method=1, fill_value=0)
                except Exception: pass
                return s[0]
        return None
    trZ, trH1, trH2 = get("Z"), get("1", "N"), get("2", "E")
    if trZ is None or trH1 is None or trH2 is None:
        return None
    fs = float(trZ.stats.sampling_rate)
    sos = signal.butter(4, [2, 8], "band", fs=fs, output="sos")
    n = int(round((PRE + POST) * fs))
    acc = np.zeros((3, n)); cnt = 0
    for t in times:
        a = UTCDateTime(t)
        try:
            z = trZ.slice(a - PRE, a + POST).data.astype(float)
            h1 = trH1.slice(a - PRE, a + POST).data.astype(float)
            h2 = trH2.slice(a - PRE, a + POST).data.astype(float)
        except Exception:
            continue
        if min(len(z), len(h1), len(h2)) < n:
            continue
        z, h1, h2 = z[:n], h1[:n], h2[:n]
        if not (np.all(np.isfinite(z)) and np.all(np.isfinite(h1)) and np.all(np.isfinite(h2))):
            continue
        W = np.vstack([signal.sosfiltfilt(sos, signal.detrend(x)) for x in (z, h1, h2)])
        m = np.abs(W).max()
        if m == 0:
            continue
        acc += W / m; cnt += 1
    return (acc, cnt, fs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--fam", required=True); ap.add_argument("--mf", required=True)
    ap.add_argument("--cc-min", type=float, default=0.88); ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args()
    mf = pd.read_csv(a.mf, usecols=["template", "time", "cc"])
    d = mf[(mf.template == a.fam) & (mf.cc >= a.cc_min)].copy()
    if len(d) > a.n:
        d = d.nlargest(a.n, "cc")
    d["t"] = pd.to_datetime(d.time, utc=True)
    print(f"[{a.fam}] {len(d)} detections (cc>={a.cc_min}) to stack", flush=True)
    jobs = defaultdict(list)
    for r in d.itertuples(index=False):
        jobs[(r.t.year, r.t.dayofyear)].append(r.t.value / 1e9)
    payloads = [(a.net, a.sta, k[0], k[1], v) for k, v in jobs.items()]
    acc = None; total = 0; fs = 100.0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            r = f.result()
            if r is None: continue
            A, c, fs = r
            if c == 0: continue
            acc = A if acc is None else acc + A
            total += c
    stack = acc / total
    t = np.arange(stack.shape[1]) / fs - PRE
    np.save(f"{HERE}/data/stack3c_{a.sta.lower()}_{a.fam.replace('.','p')}.npy", stack)

    labels = ["Z (vertical, EHZ)", "EH1 (horizontal)", "EH2 (horizontal)"]
    colors = ["#222222", "#1f77b4", "#d62728"]
    ymax = np.abs(stack).max() * 1.15
    fig, ax = plt.subplots(3, 1, figsize=(11, 7), sharex=True, sharey=True)
    for i in range(3):
        ax[i].plot(t, stack[i], color=colors[i], lw=0.8)
        ax[i].axvline(0, color="green", ls=":", lw=1, label="detection / S anchor")
        ax[i].set_ylabel(labels[i], fontsize=9)
        ax[i].set_ylim(-ymax, ymax); ax[i].grid(alpha=0.25)
        rms = np.sqrt(np.mean(stack[i] ** 2))
        ax[i].text(0.99, 0.9, f"RMS={rms:.3f}", transform=ax[i].transAxes, ha="right", fontsize=8)
    ax[0].legend(fontsize=8, loc="upper left")
    hv = np.sqrt((stack[1] ** 2 + stack[2] ** 2).mean()) / (np.sqrt((stack[0] ** 2).mean()) + 1e-9)
    ax[0].set_title(f"{a.sta} LFE family {a.fam}  —  {total} detections stacked (cc>={a.cc_min}), "
                    f"3-comp coherent stack, 2-8 Hz   |  H/V = {hv:.1f}", fontsize=10)
    ax[2].set_xlabel("time relative to detection (s)")
    fig.tight_layout()
    out = f"{HERE}/figures/family3c_{a.sta.lower()}_{a.fam.replace('.','p')}.png"
    fig.savefig(out, dpi=130)
    print(f"stacked {total} detections; H/V={hv:.2f}; saved {out}")
    print("FAMILY 3COMP DONE")


if __name__ == "__main__":
    main()
