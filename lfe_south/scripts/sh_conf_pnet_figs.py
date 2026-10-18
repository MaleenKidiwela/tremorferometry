#!/usr/bin/env python
"""PNET fine-tune figures for SH_CONF (2026-10-01). Reads saved files only (no waveform data needed):
results/sh_conf/pnet/models/train_ft1.csv, pnet/train/{X,Y,meta}, models, eval_*.csv, candidates_conf.parquet."""
import sys, numpy as np, pandas as pd, torch, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts"); import p_scan
R = "/home/jovyan/tremorferometry/lfe_south/results/sh_conf"; FS = 20.0
W05 = "/home/jovyan/tremorferometry/lfe_picker/results/p1/models/pnet_w05.pt"; FT = f"{R}/pnet/models/ft1.pt"
HO = ["BK.JEWT", "NC.KBU", "PB.B933"]

h = pd.read_csv(f"{R}/pnet/models/train_ft1.csv")
Mt = pd.read_parquet(f"{R}/pnet/train/meta.parquet")
dates = sorted(Mt.date.unique()); vd = set(dates[5::10])            # same validation dates as the fine-tune
vidx = np.nonzero(Mt.date.isin(vd).to_numpy())[0]
rng = np.random.default_rng(1); vs = np.sort(rng.choice(vidx, min(20000, len(vidx)), replace=False))
X = np.load(f"{R}/pnet/train/X.npy", mmap_mode="r"); Y = np.load(f"{R}/pnet/train/Y.npy", mmap_mode="r")
x = np.asarray(X[vs], dtype=np.float32); y = np.asarray(Y[vs], dtype=np.float32) / 255.0
x /= np.maximum(np.abs(x).max(axis=(1, 2), keepdims=True), 1e-12)
kind = Mt.kind.to_numpy()[vs]; sta = Mt.netsta.to_numpy()[vs]

def run(path):
    m, _ = p_scan.load_model(path); out = []
    with torch.no_grad():
        for i in range(0, len(x), 1024):
            out.append(m(torch.from_numpy(x[i:i + 1024]).cuda()).float().cpu().numpy())
    return np.concatenate(out)
Pw, Pf = run(W05), run(FT)
from scipy.ndimage import maximum_filter1d
def peak(P):  # value at the labelled S (positives) / window centre (negatives), +-0.5 s max
    mp = maximum_filter1d(P, 21, axis=1, mode="constant")
    pos = kind == 0; at = np.where(pos, y.argmax(1), 600)
    return mp[np.arange(len(P)), at]
kw, kf = peak(Pw), peak(Pf); pos = kind == 0
def auc(a, b):
    from scipy.stats import rankdata
    r = rankdata(np.r_[a, b]); return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))

st = pd.read_csv(f"{R}/eval_pnet_station.csv")
cand = pd.read_parquet(f"{R}/candidates_conf.parquet", columns=["series", "role", "m_any", "pnet_max_f", "pnet_max_w"])
ind = cand[cand.role.isin(["ho2024", "fwd_ovl"])]
ea = pd.read_csv(f"{R}/eval_auc.csv")

fig = plt.figure(figsize=(15, 13)); gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.55)
# A learning curve
a = fig.add_subplot(gs[0, 0]); a.plot(h.step, h.train_loss, "o-", c="C0", label="train loss")
a.plot(h.step, h.val_loss, "s-", c="C1", label="validation loss"); a.set_xlabel("step"); a.set_ylabel("loss (masked BCE)")
b = a.twinx(); b.plot(h.step, h.val_auc, "^--", c="C2", label="validation AUC"); b.set_ylabel("validation AUC", color="C2")
b.axvline(5000, c="0.6", ls=":"); b.text(5050, 0.70, "kept\n(step 5000)", fontsize=8, color="0.4")
a.legend(loc="upper center", fontsize=8); a.set_title("A  Learning curve (validation = every 10th build date)", fontsize=10, loc="left")
# B peak distributions on validation windows
a = fig.add_subplot(gs[0, 1]); bins = np.linspace(0, 1, 41)
for k, c, lab in [(kw, "C7", "Vancouver Is. weights"), (kf, "C3", "retrained")]:
    a.hist(k[pos], bins, histtype="step", lw=2, color=c, density=True, label=f"{lab}: LFE arrival")
    a.hist(k[~pos], bins, histtype="step", lw=1.2, ls="--", color=c, density=True, label=f"{lab}: no LFE")
a.set_yscale("log"); a.set_xlabel("PNET peak (±0.5 s)"); a.set_ylabel("density"); a.legend(fontsize=7)
a.set_title(f"B  Validation windows: AUC {auc(kw[pos], kw[~pos]):.2f} → {auc(kf[pos], kf[~pos]):.2f}", fontsize=10, loc="left")
# C per-station AUC on validation windows (training stations only exist here)
a = fig.add_subplot(gs[0, 2]); rows = []
for s in np.unique(sta):
    m = sta == s
    if (m & pos).sum() >= 50 and (m & ~pos).sum() >= 50:
        rows.append((s, auc(kw[m & pos], kw[m & ~pos]), auc(kf[m & pos], kf[m & ~pos])))
rows = pd.DataFrame(rows, columns=["sta", "w", "f"])
a.scatter(rows.w, rows.f, c="C3"); [a.annotate(r.sta.split(".")[1], (r.w, r.f), fontsize=7) for r in rows.itertuples()]
a.plot([0.4, 1], [0.4, 1], c="0.6", ls=":"); a.set_xlabel("AUC, Vancouver Is. weights"); a.set_ylabel("AUC, retrained")
a.set_title("C  Per station, validation windows", fontsize=10, loc="left")
# D per-station AUC at candidate times (independent days), any family
a = fig.add_subplot(gs[1, 0]); p = st.pivot(index="netsta", columns="model", values="auc_any"); ho = p.index.isin(HO)
a.scatter(p["off-the-shelf"][~ho], p["fine-tuned"][~ho], c="C0", label="used in retraining")
a.scatter(p["off-the-shelf"][ho], p["fine-tuned"][ho], c="C3", marker="D", label="held out of retraining")
for s, r in p.iterrows(): a.annotate(s.split(".")[1], (r["off-the-shelf"], r["fine-tuned"]), fontsize=7)
a.plot([0.4, 0.8], [0.4, 0.8], c="0.6", ls=":"); a.set_xlabel("AUC, Vancouver Is. weights"); a.set_ylabel("AUC, retrained")
a.legend(fontsize=7); a.set_title("D  Per station at candidate times (matches any family)", fontsize=10, loc="left")
# E same, same family
a = fig.add_subplot(gs[1, 1]); p = st.pivot(index="netsta", columns="model", values="auc_same")
a.scatter(p["off-the-shelf"][~ho], p["fine-tuned"][~ho], c="C0"); a.scatter(p["off-the-shelf"][ho], p["fine-tuned"][ho], c="C3", marker="D")
for s, r in p.iterrows(): a.annotate(s.split(".")[1], (r["off-the-shelf"], r["fine-tuned"]), fontsize=7)
a.plot([0.4, 0.7], [0.4, 0.7], c="0.6", ls=":"); a.set_xlabel("AUC, Vancouver Is. weights"); a.set_ylabel("AUC, retrained")
a.set_title("E  Per station at candidate times (same family)", fontsize=10, loc="left")
# F network max at candidates: matched / unmatched / scrambled control
a = fig.add_subplot(gs[1, 2])
for sel, c, lab in [((ind.series == "data") & (ind.m_any == 1), "C2", "data, matches his event"),
                    ((ind.series == "data") & (ind.m_any == 0), "C1", "data, no match"),
                    (ind.series == "c1", "k", "scrambled control")]:
    a.hist(ind.pnet_max_f[sel].dropna(), bins, histtype="step", lw=1.8, color=c, density=True, label=f"{lab} (n={sel.sum()})")
a.set_xlabel("best-station PNET (retrained)"); a.set_ylabel("density"); a.legend(fontsize=7)
a.set_title("F  Candidates on independent days", fontsize=10, loc="left")
# G confidence AUC with / without PNET
a = fig.add_subplot(gs[2, 0]); e = ea[ea.eval_set.isin(["ho2024", "fwd_ovl"])]
cols = [("auc_MF_S1", "detector + rules"), ("auc_NO_PNET", "+ site/density (no PNET)"), ("auc_PNET_W_ONLY", "+ PNET old weights"), ("auc_FULL", "+ PNET retrained")]
xx = np.arange(len(cols)); wdt = 0.2
for j, (es, lb) in enumerate([(s, l) for s in ["ho2024", "fwd_ovl"] for l in ["same", "any"]]):
    r = e[(e.eval_set == es) & (e.label == lb)].iloc[0]
    a.bar(xx + (j - 1.5) * wdt, [r[c] for c, _ in cols], wdt, label=f"{es} · {lb}")
a.set_xticks(xx); a.set_xticklabels([l for _, l in cols], rotation=15, fontsize=8); a.set_ylim(0.75, 0.9); a.set_ylabel("confidence AUC")
a.legend(fontsize=7, ncol=2); a.set_title("G  What PNET adds to the per-repeat confidence", fontsize=10, loc="left")
# H examples
t = np.arange(1200) / FS - 30
mid = (y.argmax(1) > 300) & (y.argmax(1) < 900)
ex = [np.nonzero(pos & mid & (kf > 0.6) & (kw < 0.3))[0], np.nonzero(pos & mid & (kf > 0.5) & (kw > 0.5))[0], np.nonzero(~pos & (kw > 0.5) & (kf < 0.2))[0]]
titles = ["LFE: missed by old weights, found after retraining", "LFE: found by both", "no LFE: old weights false alarm, removed"]
sub = gs[2, 1:].subgridspec(1, 3, wspace=0.25)
for j, (ii, tt) in enumerate(zip(ex, titles)):
    a = fig.add_subplot(sub[0, j])
    if len(ii):
        i = ii[len(ii) // 2]
        a.plot(t, x[i, 0] * 0.45 + 1.6, c="0.3", lw=0.5); a.plot(t, y[i] * 0 + 0, alpha=0)
        a.plot(t, Pw[i], c="C7", label="old weights"); a.plot(t, Pf[i], c="C3", label="retrained")
        if pos[i]: a.axvline(t[y[i].argmax()], c="C2", ls="--", lw=1, label="his S")
        a.set_title(f"{sta[i]}  ({tt})", fontsize=8)
        c0 = t[y[i].argmax()] if pos[i] else 0.0
    a.set_xlim(c0 - 12, c0 + 12); a.set_ylim(-0.05, 2.1); a.set_xlabel("s (window time)"); a.set_yticks([0, 0.5, 1])
    if j == 0: a.legend(fontsize=7, loc="upper left"); a.set_ylabel("PNET probability (trace above: Z)")
fig.text(0.67, 0.335, "H  Example validation windows", fontsize=10)
fig.suptitle("PNET retrained on Shelly's Mendocino events (build days, 14 stations; JEWT, KBU, B933 held out)", fontsize=12)
out = f"{R}/figures/pnet_finetune.png"; fig.savefig(out, dpi=110, bbox_inches="tight"); print(out)
print(rows.round(3).to_string()); print("val AUC w/f", auc(kw[pos], kw[~pos]), auc(kf[pos], kf[~pos]), "n", pos.sum(), (~pos).sum())
print("examples", [len(i) for i in ex])
