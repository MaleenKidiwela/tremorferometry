"""Multi-class separability: LFE (LIN) vs EARTHQUAKE (EQ) vs NOISE (RAND).

The picker must reject earthquakes, not just silence. This quantifies LIN-vs-EQ separability
(the hard test), the per-feature LFE-vs-EQ contrast, and a 3-class embedding.

Usage: PYTHONPATH=src python lfe_features/analyze_multiclass.py --sta B011
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.decomposition import PCA
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
HERE = os.path.dirname(__file__); POL = ["hv_ratio", "rectilinearity", "planarity"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--sta", required=True); a = ap.parse_args(); s = a.sta.lower()
    ref = pd.read_parquet(f"{HERE}/data/feat_ref_{s}.parquet")
    eq = pd.read_parquet(f"{HERE}/data/feat_eq_{s}.parquet")
    df = pd.concat([ref, eq], ignore_index=True)
    vert = [c for c in FEATURE_ORDER if c in df.columns and c not in POL]
    pol = [c for c in POL if c in df.columns and df[c].notna().mean() > 0.9]
    cols = vert + pol
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)
    print(f"=== {a.sta} multiclass: " + ", ".join(f"{k}={v}" for k, v in df.label.value_counts().items()) +
          f" | {len(cols)} features ===")
    cv = StratifiedKFold(5, shuffle=True, random_state=0)

    def auc(c1, c2, drop_amp=False):
        use = [c for c in cols if not (drop_amp and c in ("snr", "env_crest"))]
        m = df.label.isin([c1, c2]); y = (df[m].label == c1).astype(int).values
        X = StandardScaler().fit_transform(df[m][use].values)
        return cross_val_score(RandomForestClassifier(300, class_weight="balanced", n_jobs=-1, random_state=0),
                               X, y, cv=cv, scoring="roc_auc").mean()

    print("\nPairwise RandomForest CV-AUC:")
    print(f"  LIN vs RAND : {auc('LIN','RAND'):.3f}")
    print(f"  LIN vs EQ   : {auc('LIN','EQ'):.3f}   (HARD test: reject earthquakes)")
    print(f"  LIN vs EQ   : {auc('LIN','EQ',drop_amp=True):.3f}   (without amplitude features)")
    print(f"  EQ  vs RAND : {auc('EQ','RAND'):.3f}")

    # per-feature LIN vs EQ effect size + class means
    def d(c1, c2, f):
        a1 = df[df.label == c1][f].values; a0 = df[df.label == c2][f].values
        sp = np.sqrt((np.nanvar(a1) + np.nanvar(a0)) / 2) + 1e-9
        return (np.nanmean(a1) - np.nanmean(a0)) / sp
    print("\nLIN-vs-EQ top features (Cohen's d):")
    for c, v in sorted(((c, d("LIN", "EQ", c)) for c in cols), key=lambda t: -abs(t[1]))[:10]:
        print(f"   {c:14s} d={v:+.2f}")
    print("\nClass means on key features:")
    keyf = ["sp_centroid", "sp_slope", "hf_ratio", "bf_2_4", "bf_16up", "hv_ratio", "rectilinearity"]
    keyf = [k for k in keyf if k in cols]
    print(df.groupby("label")[keyf].mean().round(3).to_string())

    # 3-class PCA embedding
    X = StandardScaler().fit_transform(df[cols].values); Z = PCA(2).fit_transform(X)
    fig, ax = plt.subplots(figsize=(8, 7))
    col = {"LIN": "#1a9850", "EQ": "#7570b3", "RAND": "#999999"}
    for g in ["RAND", "EQ", "LIN"]:
        m = (df.label == g).values
        ax.scatter(Z[m, 0], Z[m, 1], s=4, alpha=0.3, c=col[g], label=f"{g} ({m.sum()})")
    ax.legend(markerscale=4); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title(f"{a.sta}: LFE vs EARTHQUAKE vs NOISE feature space")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/multiclass_pca_{s}.png", dpi=120)
    print(f"\nsaved figures/multiclass_pca_{s}.png")
    print("MULTICLASS DONE")


if __name__ == "__main__":
    main()
