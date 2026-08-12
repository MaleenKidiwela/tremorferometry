"""Stage-2: structure in feature space.

- standardize features, PCA (+ UMAP if available)
- separation: how distinct are FAIL vs GOLD families? (per-detection logistic AUC + per-family)
- feature variability & redundancy (correlation, within- vs between-family variance)
- 2-D scatter colored by grade

Usage: python lfe_features/cluster_explore.py --sta B927 [--sta2 ...]  (pools given stations)
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
HERE = os.path.dirname(__file__)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", nargs="+", required=True)
    a = ap.parse_args()
    dfs = []
    for s in a.sta:
        p = f"{HERE}/data/feat_{s.lower()}.parquet"
        if os.path.exists(p):
            d = pd.read_parquet(p); d["station"] = s.upper(); dfs.append(d)
    df = pd.concat(dfs, ignore_index=True)
    POL = ["hv_ratio", "rectilinearity", "planarity"]
    vert = [c for c in FEATURE_ORDER if c in df.columns and c not in POL]
    df = df.replace([np.inf, -np.inf], np.nan)
    # require vertical features present; only add polarization if non-NaN for (nearly) all rows
    df = df.dropna(subset=vert)
    pol = [c for c in POL if c in df.columns and df[c].notna().mean() > 0.9]
    if pol:
        df = df.dropna(subset=pol)
    feat = vert + pol
    print(f"pooled windows: {len(df):,}  features: {len(feat)} (vert {len(vert)}+pol {len(pol)})  stations: {a.sta}")
    print("grade counts:\n", df.grade.value_counts().to_string())

    X = StandardScaler().fit_transform(df[feat].values)

    # --- GOLD vs FAIL separability (per-window) ---
    mask = df.grade.isin(["GOLD", "FAIL"])
    if mask.sum() > 50 and df[mask].grade.nunique() == 2:
        y = (df[mask].grade == "GOLD").astype(int).values
        clf = LogisticRegression(max_iter=2000)
        auc = cross_val_score(clf, X[mask.values], y, cv=4, scoring="roc_auc").mean()
        clf.fit(X[mask.values], y)
        imp = sorted(zip(feat, clf.coef_[0]), key=lambda t: -abs(t[1]))
        print(f"\nGOLD-vs-FAIL per-window logistic AUC (4-fold): {auc:.3f}")
        print("top discriminative features (|coef|):")
        for name, c in imp[:10]:
            print(f"   {name:14s} {c:+.2f}")

    # --- PCA scatter ---
    pca = PCA(n_components=2).fit(X)
    Z = pca.transform(X)
    print(f"\nPCA explained var: PC1 {pca.explained_variance_ratio_[0]:.2f}  PC2 {pca.explained_variance_ratio_[1]:.2f}")
    fig, ax = plt.subplots(figsize=(8, 7))
    colors = {"GOLD": "#1a9850", "TRUSTED": "#91cf60", "UNDETERMINED": "#fee08b",
              "FAIL": "#d73027", "NA": "#999999"}
    for g, c in colors.items():
        m = df.grade == g
        if m.sum():
            ax.scatter(Z[m.values, 0], Z[m.values, 1], s=3, alpha=0.25, c=c, label=f"{g} ({m.sum()})")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title(f"LFE feature space (PCA) — {'+'.join(a.sta)}")
    ax.legend(markerscale=4, fontsize=8)
    os.makedirs(f"{HERE}/figures", exist_ok=True)
    fig_path = f"{HERE}/figures/featspace_pca_{'_'.join(s.lower() for s in a.sta)}.png"
    fig.tight_layout(); fig.savefig(fig_path, dpi=120)
    print(f"saved {fig_path}")

    # --- feature redundancy ---
    corr = pd.DataFrame(X, columns=feat).corr().abs()
    hi = [(feat[i], feat[j], corr.iloc[i, j]) for i in range(len(feat))
          for j in range(i + 1, len(feat)) if corr.iloc[i, j] > 0.85]
    print(f"\nhighly-correlated feature pairs (|r|>0.85): {len(hi)}")
    for x, yv, r in sorted(hi, key=lambda t: -t[2])[:8]:
        print(f"   {x} ~ {yv}: {r:.2f}")
    print("CLUSTER EXPLORE DONE")


if __name__ == "__main__":
    main()
