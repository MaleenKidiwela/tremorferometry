"""Categorize each discovered family by running the PICKER on its coherent 3-comp stack.

For each family: stack its member detections (3-comp, matched-filter-aligned), extract the
hand-crafted features, and apply the saved tremor-picker (LFE / EQ / BLAST / TNOISE). This is
the uncompressed family-level categorization (not biased by the candidate P-prefilter).

Usage: PYTHONPATH=src python lfe_features/score_family_stacks_picker.py --net PB --sta B011 \
         --members data/b011_disc_p50.members.parquet --summary data/b011_disc_p50.summary.csv --out data/family_picker_p50.csv
"""
import argparse, os, sys
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import joblib
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all
PRE, POST = 10.0, 30.0
HERE = os.path.dirname(__file__)


def cut_day(payload):
    net, sta, year, julday, items, target_fs, zonly = payload   # items: (family_id, epoch)
    from obspy import read, UTCDateTime
    path = f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"
    if not os.path.exists(path):
        return {}
    try:
        st = read(path)
    except Exception:
        return {}
    # resample the day to target_fs BEFORE stacking so the stack features match the picker's training fs
    # (broadband picker = 40 Hz; native 100 Hz EHZ/HHZ stacks would give wrong-Nyquist features).
    if target_fs:
        from fractions import Fraction
        from scipy.signal import resample_poly
        for tr in st:
            nat = float(tr.stats.sampling_rate)
            if abs(nat - target_fs) > 1e-6:
                fr = Fraction(int(round(target_fs)), int(round(nat))).limit_denominator(1000)
                try:
                    tr.data = resample_poly(tr.data.astype("float64"), fr.numerator, fr.denominator).astype("float32")
                    tr.stats.sampling_rate = target_fs
                except Exception:
                    pass
    def get(*c):
        for ch in c:
            s = st.select(channel=f"*{ch}")
            if len(s):
                try: s = s.merge(method=1, fill_value=0)
                except Exception: pass
                return s[0]
        return None
    trZ = get("Z")
    if trZ is None:
        return {}
    trH1, trH2 = (None, None) if zonly else (get("1", "N"), get("2", "E"))
    nr = 3 if (trH1 is not None and trH2 is not None) else 1     # 3-comp (borehole) or Z-only (broadband/EHZ)
    fs = float(trZ.stats.sampling_rate); exp = int(round((PRE + POST) * fs))
    acc = {}
    for fam, ep in items:
        a = UTCDateTime(ep)
        try:
            z = trZ.slice(a - PRE, a + POST).data.astype(float)[:exp]
        except Exception:
            continue
        if len(z) < exp or np.std(z) == 0:
            continue
        if nr == 3:
            try:
                h1 = trH1.slice(a - PRE, a + POST).data.astype(float)[:exp]
                h2 = trH2.slice(a - PRE, a + POST).data.astype(float)[:exp]
            except Exception:
                continue
            if min(len(h1), len(h2)) < exp:
                continue
            m = max(np.abs(z).max(), np.abs(h1).max(), np.abs(h2).max()) + 1e-20
            W = np.vstack([z, h1, h2]) / m
        else:
            m = np.abs(z).max() + 1e-20
            W = (z / m).reshape(1, -1)
        if fam not in acc:
            acc[fam] = [np.zeros((W.shape[0], exp)), 0]
        if acc[fam][0].shape[0] != W.shape[0]:
            continue
        acc[fam][0] += W; acc[fam][1] += 1
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--members", required=True); ap.add_argument("--summary", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--max-per-fam", type=int, default=400)
    a = ap.parse_args()
    B = joblib.load(f"{HERE}/models/tremor_picker_{a.sta.lower()}.joblib")
    model, scaler, cols, classes = B["model"], B["scaler"], B["cols"], B["classes"]
    fs = 100.0 if a.net == "PB" else 40.0            # feature fs (matches the picker's training fs)
    target_fs = None if a.net == "PB" else 40.0      # non-PB (broadband/land) days -> resample to 40 before stacking
    POL = {"hv_ratio", "rectilinearity", "planarity"}
    zonly = not (set(cols) & POL)                    # vertical-only picker (broadband) -> stack Z only (EHZ has no H)

    mem = pd.read_parquet(a.members); mem["t"] = pd.to_datetime(mem.time, utc=True)
    parts = [g.sample(min(a.max_per_fam, len(g)), random_state=1) for _, g in mem.groupby("family_id")]
    mem = pd.concat(parts, ignore_index=True)
    mem["y"] = mem.t.dt.year; mem["j"] = mem.t.dt.dayofyear
    jobs = defaultdict(list)
    for r in mem.itertuples(index=False):
        jobs[(int(r.y), int(r.j))].append((r.family_id, r.t.value / 1e9))
    payloads = [(a.net, a.sta, k[0], k[1], v, target_fs, zonly) for k, v in jobs.items()]
    print(f"[{a.sta}] stacking {mem.family_id.nunique()} families over {len(payloads)} days...", flush=True)
    tot = {}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(cut_day, p) for p in payloads]):
            for fam, val in f.result().items():
                if fam not in tot: tot[fam] = val
                else: tot[fam][0] += val[0]; tot[fam][1] += val[1]
    rows = []
    for fam, (acc, n) in tot.items():
        if n < 3:                       # seed clusters are 3-9 members; stack what we have
            continue
        stk = acc / n
        z, h1, h2 = (stk[0], None, None) if stk.shape[0] == 1 else tuple(stk)   # Z-only or 3-comp
        feats = compute_all(z, fs, PRE, h1, h2)
        row = {c: feats.get(c, np.nan) for c in cols}
        if any(not np.isfinite(v) for v in row.values()):
            continue
        proba = model.predict_proba(scaler.transform(pd.DataFrame([row])[cols]))[0]
        pr = dict(zip(classes, proba))
        rows.append(dict(fam=fam, nstack=n, pred=classes[int(np.argmax(proba))],
                         p_lfe=pr.get("LFE", np.nan), p_eq=pr.get("EQ", np.nan),
                         p_blast=pr.get("BLAST", np.nan), p_tnoise=pr.get("TNOISE", np.nan)))
    d = pd.DataFrame(rows)
    s = pd.read_csv(a.summary).set_index("family_id")
    d = d.join(s[["lat", "lon"]], on="fam")
    d.to_csv(a.out, index=False)
    print(f"\n=== PICKER categorization of {len(d)} family stacks ===")
    print(d.pred.value_counts().to_string())
    print(f"\nfamily-stack P(LFE): median {d.p_lfe.median():.2f}")
    for thr in [0.5, 0.7, 0.9]:
        print(f"  P(LFE)>={thr}: {(d.p_lfe>=thr).sum()} families")
    print(f"saved {a.out}")
    print("PICKER FAMILY SCORING DONE")


if __name__ == "__main__":
    main()
