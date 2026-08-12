"""Stage-3: spectrogram convolutional AUTOENCODER -> learned fingerprints.

Trains a small conv AE (unsupervised) on log-spectrograms of LIN/EQ/RAND/FAM windows, takes
the latent vector as a learned fingerprint, embeds with UMAP, clusters (HDBSCAN), and asks:
do LFE/EQ/noise separate in LEARNED space, and where do the densified families land?

Usage: python lfe_features/train_autoencoder.py --sta B011 --epochs 40
"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(__file__)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True); ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--latent", type=int, default=32); ap.add_argument("--batch", type=int, default=256)
    a = ap.parse_args(); s = a.sta.lower()
    import torch, torch.nn as nn
    torch.manual_seed(0)

    z = np.load(f"{HERE}/data/spec_{s}.npz", allow_pickle=True)
    X = z["X"].astype(np.float32); y = z["y"].astype(str); fam = z["fam"].astype(str)
    NF, NT = X.shape[1], X.shape[2]
    print(f"[{s}] specs {X.shape}, classes {dict(zip(*np.unique(y, return_counts=True)))}")
    Xt = torch.from_numpy(X).unsqueeze(1)

    class AE(nn.Module):
        def __init__(self, latent):
            super().__init__()
            self.enc = nn.Sequential(
                nn.Conv2d(1, 16, 3, 2, 1), nn.ReLU(),
                nn.Conv2d(16, 32, 3, 2, 1), nn.ReLU(),
                nn.Conv2d(32, 64, 3, 2, 1), nn.ReLU())
            with torch.no_grad():
                self.shape = self.enc(torch.zeros(1, 1, NF, NT)).shape
            flat = int(np.prod(self.shape[1:]))
            self.fc1 = nn.Linear(flat, latent); self.fc2 = nn.Linear(latent, flat)
            self.dec = nn.Sequential(
                nn.ConvTranspose2d(64, 32, 3, 2, 1, output_padding=1), nn.ReLU(),
                nn.ConvTranspose2d(32, 16, 3, 2, 1, output_padding=1), nn.ReLU(),
                nn.ConvTranspose2d(16, 1, 3, 2, 1, output_padding=1))

        def encode(self, x):
            return self.fc1(self.enc(x).flatten(1))

        def forward(self, x):
            h = self.fc1(self.enc(x).flatten(1))
            d = self.fc2(h).view(-1, *self.shape[1:])
            return self.dec(d)[:, :, :NF, :NT]

    net = AE(a.latent); opt = torch.optim.Adam(net.parameters(), 1e-3); loss_fn = nn.MSELoss()
    n = len(Xt); idx = np.arange(n)
    for ep in range(a.epochs):
        np.random.shuffle(idx); tot = 0
        for i in range(0, n, a.batch):
            b = Xt[idx[i:i + a.batch]]
            opt.zero_grad(); out = net(b); l = loss_fn(out, b); l.backward(); opt.step()
            tot += l.item() * len(b)
        if (ep + 1) % 10 == 0:
            print(f"  epoch {ep+1}/{a.epochs}  recon MSE {tot/n:.4f}", flush=True)

    net.eval()
    with torch.no_grad():
        Z = np.concatenate([net.encode(Xt[i:i + 512]).numpy() for i in range(0, n, 512)])
    np.savez_compressed(f"{HERE}/data/ae_latent_{s}.npz", Z=Z, y=y, fam=fam)

    # supervised check: linear-probe LFE-vs-rest on learned latent
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    base = np.isin(y, ["LIN", "EQ", "RAND"])
    yb = (y[base] == "LIN").astype(int)
    auc = cross_val_score(LogisticRegression(max_iter=2000, class_weight="balanced"),
                          Z[base], yb, cv=4, scoring="roc_auc").mean()
    print(f"\nlearned-latent linear-probe LIN-vs-(EQ+RAND) AUC: {auc:.3f}")

    # embed + cluster
    try:
        import umap
        emb = umap.UMAP(n_neighbors=30, min_dist=0.1, random_state=0).fit_transform(Z)
        tag = "UMAP"
    except Exception:
        from sklearn.decomposition import PCA
        emb = PCA(2).fit_transform(Z); tag = "PCA"
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    col = {"LIN": "#1a9850", "EQ": "#7570b3", "RAND": "#bbbbbb",
           "FAM_GOLD": "#d95f02", "FAM_TRUSTED": "#e6ab02", "FAM_UNDETERMINED": "#a6761d",
           "FAM_FAIL": "#e7298a", "FAM_NA": "#666666"}
    for g in ["RAND", "EQ", "LIN"]:
        m = y == g
        ax[0].scatter(emb[m, 0], emb[m, 1], s=4, alpha=0.3, c=col[g], label=f"{g} ({m.sum()})")
    ax[0].legend(markerscale=4); ax[0].set_title(f"{a.sta} learned latent ({tag}): LFE/EQ/noise")
    for g in [k for k in col if k.startswith("FAM")]:
        m = y == g
        if m.sum():
            ax[1].scatter(emb[m, 0], emb[m, 1], s=5, alpha=0.4, c=col[g], label=f"{g} ({m.sum()})")
    # overlay LIN centroid region
    mlin = y == "LIN"
    ax[1].scatter(emb[mlin, 0], emb[mlin, 1], s=3, alpha=0.08, c="green")
    ax[1].legend(markerscale=4); ax[1].set_title(f"{a.sta}: densified FAMILIES in learned space (green haze=LIN)")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/ae_umap_{s}.png", dpi=120)
    print(f"saved figures/ae_umap_{s}.png + data/ae_latent_{s}.npz")
    print("AUTOENCODER DONE")


if __name__ == "__main__":
    main()
