"""DECISIVE test: do family STACKS look like LFEs in feature space?

Per-window scoring is swamped by noise-matches (~95% of densified detections). Coherent
stacking (aligned at t_d) sums the genuine repeating LFE and averages noise down -- the
same logic as the trust battery, but scored by the Lin-trained LFE feature classifier.

For station <sta>:
  - build a coherent mean-window stack per family (Z + horizontals)
  - controls: stacks of equal-N random subsets of curated LIN windows (positive) and RAND (null)
  - train LFE-vs-RAND on single ref windows; score every stack
Outputs: data/stack_plfe_<sta>.csv + figures/stack_plfe_<sta>.png + printed verdict.

Usage: PYTHONPATH=src python lfe_features/family_stack_test.py --net PB --sta B011 --slat 48.65 --slon -123.448
"""
import argparse, os, sys, glob, math
import numpy as np, pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import compute_all, FEATURE_ORDER
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier

PRE, POST, ANCHOR = 10.0, 30.0, 11.0
HERE = os.path.dirname(__file__)
POL = ["hv_ratio", "rectilinearity", "planarity"]


def day_file(net, sta, year, julday):
    return f"data/waveforms/{net}.{sta}/{year}/{julday:03d}.mseed"


def cut_day(payload):
    """Return key -> (z_sum, h1_sum, h2_sum, n, exp) partial coherent sums for this day."""
    net, sta, year, julday, items = payload   # items: (key, anchor_epoch)
    from obspy import read, UTCDateTime
    path = day_file(net, sta, year, julday)
    if not os.path.exists(path):
        return {}
    try:
        st = read(path)
    except Exception:
        return {}

    def get(*chars):
        for ch in chars:
            s = st.select(channel=f"*{ch}")
            if len(s):
                try: s = s.merge(method=1, fill_value=0)
                except Exception: pass
                return s[0]
        return None
    trZ = get("Z")
    if trZ is None:
        return {}
    fs = float(trZ.stats.sampling_rate)
    trH1 = get("1", "N"); trH2 = get("2", "E")
    exp = int(round((PRE + POST) * fs))
    acc = {}
    for key, anch in items:
        a = UTCDateTime(anch)
        try:
            z = trZ.slice(a - PRE, a + POST).data.astype(float)
        except Exception:
            continue
        if len(z) < exp or not np.all(np.isfinite(z)) or np.std(z) == 0:
            continue
        z = z[:exp]
        z = z / (np.std(z) + 1e-20)                 # normalize per window before stacking
        h1 = h2 = None
        if trH1 is not None and trH2 is not None:
            try:
                aa = trH1.slice(a - PRE, a + POST).data.astype(float)[:exp]
                bb = trH2.slice(a - PRE, a + POST).data.astype(float)[:exp]
                if len(aa) == exp and len(bb) == exp:
                    s = np.std(z) if False else 1.0
                    h1 = aa / (np.std(aa) + 1e-20); h2 = bb / (np.std(bb) + 1e-20)
            except Exception:
                h1 = h2 = None
        if key not in acc:
            acc[key] = [np.zeros(exp), np.zeros(exp), np.zeros(exp), 0, exp]
        acc[key][0] += z
        if h1 is not None:
            acc[key][1] += h1; acc[key][2] += h2
        acc[key][3] += 1
    return acc


def gather(net, sta, jobs, workers):
    payloads = [(net, sta, k[0], k[1], v) for k, v in jobs.items()]
    tot = defaultdict(lambda: None)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(cut_day, p) for p in payloads]
        for f in as_completed(futs):
            for key, val in f.result().items():
                if tot[key] is None:
                    tot[key] = val
                else:
                    tot[key][0] += val[0]; tot[key][1] += val[1]
                    tot[key][2] += val[2]; tot[key][3] += val[3]
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True); ap.add_argument("--sta", required=True)
    ap.add_argument("--slat", type=float, required=True); ap.add_argument("--slon", type=float, required=True)
    ap.add_argument("--per-fam", type=int, default=400); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--n-ctrl", type=int, default=4000); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args(); s = a.sta.lower(); rng = np.random.RandomState(a.seed)

    files = glob.glob(f"data/waveforms/{a.net}.{a.sta}/*/*.mseed")
    have = set((int(f.split("/")[-2]), int(os.path.basename(f).split(".")[0])) for f in files)
    yrs = sorted(set(y for y, _ in have))
    print(f"[{a.sta}] waveform years {yrs[0]}-{yrs[-1]} ({len(have)} days)", flush=True)

    # ---- family detection windows (only days we have) ----
    m = pd.read_csv(f"data/mf_{s}_all.csv", usecols=["template", "time", "cc"])
    m["t"] = pd.to_datetime(m.time, utc=True)
    m["year"] = m.t.dt.year; m["julday"] = m.t.dt.dayofyear
    m = m[[ (r.year, r.julday) in have for r in m.itertuples(index=False)]]
    A = pd.read_csv("results/family_verdicts/ALL_families_labeled.csv")
    grade = dict(zip(A[A.station == a.sta].fam, A[A.station == a.sta].grade))
    fam_jobs = defaultdict(list)
    for fam, g in m.groupby("template"):
        gg = g.sample(min(a.per_fam, len(g)), random_state=a.seed)
        for r in gg.itertuples(index=False):
            fam_jobs[(int(r.year), int(r.julday))].append((("FAM", fam), r.t.value / 1e9))
    print(f"[{a.sta}] stacking {m.template.nunique()} families over {len(fam_jobs)} days...", flush=True)
    fam_acc = gather(a.net, a.sta, fam_jobs, a.workers)

    # ---- control windows: LIN (positive) and RAND (null), grouped into equal-N stacks ----
    lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv")
    dkm = np.sqrt(((lin.lat - a.slat) * 111) ** 2 + ((lin.lon - a.slon) * 111 * math.cos(math.radians(a.slat))) ** 2)
    cur = lin[(lin.N >= 6) & (lin.residual < 1.2) & (dkm < 35)].copy()
    cur["t"] = pd.to_datetime(cur.OT, utc=True); cur["year"] = cur.t.dt.year; cur["julday"] = cur.t.dt.dayofyear
    cur = cur[[ (r.year, r.julday) in have for r in cur.itertuples(index=False)]]
    cur = cur.sample(min(a.n_ctrl, len(cur)), random_state=a.seed)
    n_per = int(np.median([v[3] for v in fam_acc.values() if v]))   # match family stack depth
    n_per = max(20, n_per)
    ctrl_jobs = defaultdict(list)
    # LIN groups
    cur = cur.reset_index(drop=True)
    for i, r in enumerate(cur.itertuples(index=False)):
        ctrl_jobs[(int(r.year), int(r.julday))].append((("LIN", i // n_per), r.t.value / 1e9 + ANCHOR))
    # RAND groups
    daylist = sorted(have)
    from datetime import datetime
    for j in range(a.n_ctrl):
        y, jd = daylist[rng.randint(len(daylist))]
        base = datetime(y, 1, 1).timestamp() + (jd - 1) * 86400 + rng.uniform(60, 86340)
        ctrl_jobs[(y, jd)].append((("RAND", j // n_per), base + ANCHOR))
    print(f"[{a.sta}] control stacks at depth ~{n_per}; gathering...", flush=True)
    ctrl_acc = gather(a.net, a.sta, ctrl_jobs, a.workers)

    # ---- features on each stack, score with Lin-trained model ----
    ref = pd.read_parquet(f"{HERE}/data/feat_ref_{s}.parquet")
    # ENV features are alignment-sensitive (family stacks are matched-filter-aligned, LIN
    # stacks are loosely OT-anchored -> unfair). ROBUST = spectral + polarization only.
    ENV = ["env_crest", "env_kurt", "env_skew", "rise_time", "duration", "zcr"]
    vert = [c for c in FEATURE_ORDER if c in ref.columns and c not in POL]
    pol = [c for c in POL if c in ref.columns and ref[c].notna().mean() > 0.9]
    cols_full = vert + pol
    cols_rob = [c for c in cols_full if c not in ENV]
    ref = ref.replace([np.inf, -np.inf], np.nan).dropna(subset=cols_full)
    sub = ref[ref.label.isin(["LIN", "RAND"])]
    yb = (sub.label == "LIN").astype(int).values
    models = {}
    for tag, cc in [("full", cols_full), ("robust", cols_rob)]:
        scl = StandardScaler().fit(sub[cc])
        rf = RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0).fit(scl.transform(sub[cc]), yb)
        models[tag] = (scl, rf, cc)
    fs_native = 100.0 if a.net == "PB" else 40.0

    stacks = {}   # save mean stacks

    def score_stack(key, val):
        z, h1, h2, n, exp = val
        if n < 10:
            return None
        z = z / n
        hh1 = h1 / n if np.any(h1) else None
        hh2 = h2 / n if np.any(h2) else None
        stacks["__".join(map(str, key))] = z
        feats = compute_all(z, fs_native, PRE, hh1, hh2)
        out = {}
        for tag, (scl, rf, cc) in models.items():
            row = {c: feats.get(c, np.nan) for c in cc}
            if any(not np.isfinite(v) for v in row.values()):
                return None
            out[tag] = float(rf.predict_proba(scl.transform(pd.DataFrame([row])[cc]))[:, 1][0])
        out["n"] = n
        return out

    rows = []
    for key, val in fam_acc.items():
        r = score_stack(key, val)
        if r:
            rows.append(dict(kind="FAM", fam=key[1], grade=grade.get(key[1], "NA"),
                             p_full=r["full"], p_robust=r["robust"], n=r["n"]))
    for key, val in ctrl_acc.items():
        r = score_stack(key, val)
        if r:
            rows.append(dict(kind=key[0], fam=str(key[1]), grade=key[0],
                             p_full=r["full"], p_robust=r["robust"], n=r["n"]))
    df = pd.DataFrame(rows)
    df.to_csv(f"{HERE}/data/stack_plfe_{s}.csv", index=False)
    np.savez_compressed(f"{HERE}/data/stacks_{s}.npz", **stacks)

    print("\n=== STACK P(LFE):  full model  vs  alignment-ROBUST model (spectral+polarization only) ===")
    print(f"{'grade':14s}{'n':>5s}{'P_full':>9s}{'P_robust':>10s}{'depth':>7s}")
    for g in ["LIN", "RAND", "GOLD", "TRUSTED", "UNDETERMINED", "FAIL", "NA"]:
        d = df[df.grade == g]
        if len(d):
            print(f"{g:14s}{len(d):>5d}{d.p_full.median():>9.3f}{d.p_robust.median():>10.3f}{int(d.n.median()):>7d}")
    df["p_lfe"] = df["p_robust"]

    fig, ax = plt.subplots(figsize=(9, 5))
    order = [g for g in ["LIN", "GOLD", "TRUSTED", "UNDETERMINED", "FAIL", "RAND"] if g in df.grade.values]
    col = {"LIN": "#1a9850", "GOLD": "#1a9850", "TRUSTED": "#91cf60", "UNDETERMINED": "#fee08b", "FAIL": "#d73027", "RAND": "#999999"}
    bp = ax.boxplot([df[df.grade == g].p_lfe.values for g in order], labels=order, patch_artist=True)
    for p, g in zip(bp["boxes"], order):
        p.set_facecolor(col[g])
    ax.set_ylabel("stack P(LFE)"); ax.set_title(f"{a.sta}: FAMILY STACKS vs Lin/RAND controls")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/stack_plfe_{s}.png", dpi=120)
    print(f"saved figures/stack_plfe_{s}.png + data/stack_plfe_{s}.csv")
    print("STACK TEST DONE")


if __name__ == "__main__":
    main()
