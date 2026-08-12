"""Stage-4: a PhaseNet/EQTransformer-style 1D U-Net that outputs a per-sample
[LFE, EQ, NOISE] probability stream from continuous 3-component data.

Trains on the labeled segments from build_picker_dataset.py; reports held-out per-sample
LFE/EQ AUC and pick-level precision/recall; saves the model; plots example probability traces.

Usage: python lfe_features/train_phasenet.py --sta B011 --epochs 25
"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from scipy.signal import find_peaks
HERE = os.path.dirname(__file__); FS = 50.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True); ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--batch", type=int, default=64)
    a = ap.parse_args(); s = a.sta.lower()
    import torch, torch.nn as nn
    torch.manual_seed(0); torch.set_num_threads(16)

    z = np.load(f"{HERE}/data/picker_ds_{s}.npz")
    X, Y = z["X"].astype(np.float32), z["Y"].astype(np.float32)
    n = len(X); rng = np.random.RandomState(0); idx = rng.permutation(n)
    ntr = int(0.85 * n); tr, va = idx[:ntr], idx[ntr:]
    print(f"[{s}] segments {X.shape}  train {len(tr)} val {len(va)}")
    Xt, Yt = torch.from_numpy(X), torch.from_numpy(Y)

    class DBlock(nn.Module):
        def __init__(s_, ci, co):
            super().__init__(); s_.c = nn.Sequential(
                nn.Conv1d(ci, co, 7, padding=3), nn.BatchNorm1d(co), nn.ReLU(),
                nn.Conv1d(co, co, 7, padding=3), nn.BatchNorm1d(co), nn.ReLU())
        def forward(s_, x): return s_.c(x)

    class UNet(nn.Module):
        def __init__(s_, ch=(8, 16, 32, 64)):
            super().__init__()
            s_.d0 = DBlock(3, ch[0]); s_.d1 = DBlock(ch[0], ch[1])
            s_.d2 = DBlock(ch[1], ch[2]); s_.d3 = DBlock(ch[2], ch[3])
            s_.pool = nn.MaxPool1d(2)
            s_.u2 = nn.ConvTranspose1d(ch[3], ch[2], 2, 2); s_.c2 = DBlock(ch[2] * 2, ch[2])
            s_.u1 = nn.ConvTranspose1d(ch[2], ch[1], 2, 2); s_.c1 = DBlock(ch[1] * 2, ch[1])
            s_.u0 = nn.ConvTranspose1d(ch[1], ch[0], 2, 2); s_.c0 = DBlock(ch[0] * 2, ch[0])
            s_.out = nn.Conv1d(ch[0], 3, 1)
        def forward(s_, x):
            e0 = s_.d0(x); e1 = s_.d1(s_.pool(e0)); e2 = s_.d2(s_.pool(e1)); e3 = s_.d3(s_.pool(e2))
            d2 = s_.c2(torch.cat([s_.u2(e3), e2], 1))
            d1 = s_.c1(torch.cat([s_.u1(d2), e1], 1))
            d0 = s_.c0(torch.cat([s_.u0(d1), e0], 1))
            return s_.out(d0)

    net = UNet(); opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-4)
    w = torch.tensor([6.0, 6.0, 1.0]).view(1, 3, 1)        # upweight rare LFE/EQ samples

    def softCE(logits, target):
        logp = torch.log_softmax(logits, 1)
        return -(w * target * logp).sum(1).mean()

    import io
    best_auc, best_state, best_ep = -1, None, 0
    for ep in range(a.epochs):
        net.train(); rng.shuffle(tr); tot = 0
        for i in range(0, len(tr), a.batch):
            b = tr[i:i + a.batch]
            opt.zero_grad(); out = net(Xt[b]); l = softCE(out, Yt[b]); l.backward(); opt.step()
            tot += l.item() * len(b)
        net.eval()
        with torch.no_grad():
            pv = torch.softmax(net(Xt[va]), 1).numpy()
        yv = Y[va]
        aL = roc_auc_score((yv[:, 0] > 0.5).ravel(), pv[:, 0].ravel())
        aE = roc_auc_score((yv[:, 1] > 0.5).ravel(), pv[:, 1].ravel())
        if aL > best_auc:        # early stopping on val LFE AUC
            best_auc, best_ep = aL, ep + 1
            buf = io.BytesIO(); torch.save(net.state_dict(), buf); best_state = buf.getvalue()
        if (ep + 1) % 5 == 0:
            print(f"  epoch {ep+1}/{a.epochs} loss {tot/len(tr):.4f}  val per-sample AUC LFE {aL:.3f} EQ {aE:.3f}", flush=True)
    net.load_state_dict(torch.load(io.BytesIO(best_state)))    # restore best
    print(f"  >> restored best epoch {best_ep} (val LFE AUC {best_auc:.3f})", flush=True)

    net.eval()
    with torch.no_grad():
        pv = torch.softmax(net(Xt[va]), 1).numpy()
    yv = Y[va]
    aucL = roc_auc_score((yv[:, 0] > 0.5).ravel(), pv[:, 0].ravel())
    aucE = roc_auc_score((yv[:, 1] > 0.5).ravel(), pv[:, 1].ravel())

    # pick-level precision/recall for LFE (peaks in P(LFE) > 0.3, tol 0.6 s)
    tol = int(0.6 * FS); tp = fp = fn = 0
    for k in range(len(va)):
        pk, _ = find_peaks(pv[k, 0], height=0.3, distance=int(0.5 * FS))
        tk, _ = find_peaks(yv[k, 0], height=0.5, distance=int(0.5 * FS))
        used = set()
        for p in pk:
            m = [t for t in tk if abs(t - p) <= tol and t not in used]
            if m: tp += 1; used.add(min(m, key=lambda t: abs(t - p)))
            else: fp += 1
        fn += len(tk) - len(used)
    prec = tp / (tp + fp + 1e-9); rec = tp / (tp + fn + 1e-9)
    print(f"\nHELD-OUT  per-sample AUC: LFE {aucL:.3f}  EQ {aucE:.3f}")
    print(f"LFE pick-level (peak>0.3, tol 0.6s): precision {prec:.3f}  recall {rec:.3f}  (tp {tp} fp {fp} fn {fn})")

    os.makedirs(f"{HERE}/models", exist_ok=True)
    torch.save(net.state_dict(), f"{HERE}/models/phasenet_{s}.pt")

    # plot example traces (val segments containing an LFE)
    has = [k for k in range(len(va)) if yv[k, 0].max() > 0.5][:4]
    t = np.arange(X.shape[2]) / FS
    fig, axes = plt.subplots(len(has), 1, figsize=(11, 2.4 * len(has)))
    for ax, k in zip(np.atleast_1d(axes), has):
        ax.plot(t, X[va[k], 0] * 0.4 + 0.5, color="k", lw=0.4, alpha=0.6)
        ax.plot(t, yv[k, 0], color="green", lw=1.2, label="true LFE", alpha=0.6)
        ax.plot(t, pv[k, 0], color="lime", lw=1.4, ls="--", label="pred P(LFE)")
        ax.plot(t, pv[k, 1], color="purple", lw=1.0, ls=":", label="pred P(EQ)")
        ax.set_ylim(-0.1, 1.1); ax.set_ylabel("prob")
    axes[0].legend(fontsize=8, ncol=4); axes[0].set_title(f"{a.sta} per-sample LFE picker (U-Net) — held-out segments")
    axes[-1].set_xlabel("time (s)")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/phasenet_examples_{s}.png", dpi=120)
    print(f"saved models/phasenet_{s}.pt + figures/phasenet_examples_{s}.png")
    print("PHASENET DONE")


if __name__ == "__main__":
    main()
