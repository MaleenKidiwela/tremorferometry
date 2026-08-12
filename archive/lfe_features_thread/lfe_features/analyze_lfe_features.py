"""Core analysis: are the features credible LFE discriminators, and where do the
densified trust-graded families land relative to ground-truth Lin LFEs?

Q1 (foundational): do features separate curated Lin LFEs (LIN) from random noise (RAND)?
    -> logistic + RandomForest CV-AUC; per-feature effect sizes (Cohen's d).
Q2 (the payoff): train LFE-vs-RAND on the ref set, score every densified family window,
    and ask: do GOLD families score like LFEs and FAIL like noise?

Usage: PYTHONPATH=src python lfe_features/analyze_lfe_features.py --sta B011 --densified
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
HERE = os.path.dirname(__file__)
POL = ["hv_ratio", "rectilinearity", "planarity"]


def feat_cols(df):
    vert = [c for c in FEATURE_ORDER if c in df.columns and c not in POL]
    pol = [c for c in POL if c in df.columns and df[c].notna().mean() > 0.9]
    return vert, pol


def clean(df, cols):
    return df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True)
    ap.add_argument("--densified", action="store_true")
    a = ap.parse_args(); s = a.sta.lower()

    ref = pd.read_parquet(f"{HERE}/data/feat_ref_{s}.parquet")
    vert, pol = feat_cols(ref); cols = vert + pol
    ref = clean(ref, cols)
    print(f"=== {a.sta} REF: " + ", ".join(f"{k}={v}" for k, v in ref.label.value_counts().items()) +
          f" | features {len(cols)} (vert {len(vert)}+pol {len(pol)}) ===")

    sub = ref[ref.label.isin(["LIN", "RAND"])]
    y = (sub.label == "LIN").astype(int).values
    sc = StandardScaler().fit(sub[cols].values)
    X = sc.transform(sub[cols].values)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)

    # ---- Q1: LIN vs RAND separability ----
    auc_lr = cross_val_score(LogisticRegression(max_iter=3000, class_weight="balanced"),
                             X, y, cv=cv, scoring="roc_auc").mean()
    rf = RandomForestClassifier(300, class_weight="balanced", n_jobs=-1, random_state=0)
    auc_rf = cross_val_score(rf, X, y, cv=cv, scoring="roc_auc").mean()
    print(f"\nQ1  LIN-vs-RAND CV-AUC:  logistic {auc_lr:.3f}   RandomForest {auc_rf:.3f}   (n_LIN={y.sum()}, n_RAND={(1-y).sum()})")

    # SKEPTIC CONTROL: is it just amplitude/SNR? drop SNR + envelope-amplitude features, redo
    shape_cols = [c for c in cols if c not in ("snr", "env_crest")]
    Xs = StandardScaler().fit_transform(sub[shape_cols].values)
    auc_shape = cross_val_score(RandomForestClassifier(300, class_weight="balanced", n_jobs=-1, random_state=0),
                                Xs, y, cv=cv, scoring="roc_auc").mean()
    isnr = cols.index("snr") if "snr" in cols else None
    dsnr = ((X[y == 1, isnr].mean() - X[y == 0, isnr].mean()) /
            (np.sqrt((X[y == 1, isnr].var() + X[y == 0, isnr].var()) / 2) + 1e-9)) if isnr is not None else np.nan
    print(f"    skeptic control: AUC without SNR/amplitude features = {auc_shape:.3f}  (snr Cohen's d = {dsnr:+.2f})")

    # per-feature effect size (Cohen's d) LIN vs RAND
    d = {}
    for i, c in enumerate(cols):
        a1 = X[y == 1, i]; a0 = X[y == 0, i]
        sp = np.sqrt((a1.var() + a0.var()) / 2) + 1e-9
        d[c] = (a1.mean() - a0.mean()) / sp
    print("   top discriminating features (Cohen's d, LIN-RAND):")
    for c, v in sorted(d.items(), key=lambda t: -abs(t[1]))[:12]:
        print(f"     {c:14s} d={v:+.2f}")

    # ---- Q2: place densified families ----
    if a.densified and os.path.exists(f"{HERE}/data/feat_{s}.parquet"):
        rf.fit(X, y)                                   # LFE-vs-RAND model on ground truth
        den = pd.read_parquet(f"{HERE}/data/feat_{s}.parquet")
        den = clean(den, cols)
        den["p_lfe"] = rf.predict_proba(sc.transform(den[cols].values))[:, 1]
        # reference scores (sanity): LIN windows should score high, RAND low
        sub2 = sub.copy(); sub2["p_lfe"] = rf.predict_proba(X)[:, 1]
        print(f"\nQ2  ref window P(LFE): LIN median {sub2[sub2.label=='LIN'].p_lfe.median():.2f}  "
              f"RAND median {sub2[sub2.label=='RAND'].p_lfe.median():.2f}")
        fam = den.groupby(["fam", "grade"]).p_lfe.mean().reset_index()
        print("   densified FAMILY mean P(LFE) by grade (higher = more LFE-like):")
        order = ["GOLD", "TRUSTED", "UNDETERMINED", "FAIL", "NA"]
        summ = fam.groupby("grade").p_lfe.agg(["mean", "median", "count"]).reindex([g for g in order if g in fam.grade.values])
        print(summ.to_string())

        # figure
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = {"GOLD": "#1a9850", "TRUSTED": "#91cf60", "UNDETERMINED": "#fee08b", "FAIL": "#d73027", "NA": "#999999"}
        data = [fam[fam.grade == g].p_lfe.values for g in order if g in fam.grade.values]
        labs = [g for g in order if g in fam.grade.values]
        bp = ax.boxplot(data, labels=labs, patch_artist=True)
        for patch, g in zip(bp["boxes"], labs):
            patch.set_facecolor(colors[g])
        ax.axhline(sub2[sub2.label == "LIN"].p_lfe.median(), color="green", ls="--", label="LIN median")
        ax.axhline(sub2[sub2.label == "RAND"].p_lfe.median(), color="red", ls="--", label="RAND median")
        ax.set_ylabel("family mean P(LFE)"); ax.set_title(f"{a.sta}: densified families vs Lin ground truth")
        ax.legend()
        os.makedirs(f"{HERE}/figures", exist_ok=True)
        fig.tight_layout(); fig.savefig(f"{HERE}/figures/plfe_by_grade_{s}.png", dpi=120)
        fam.to_csv(f"{HERE}/data/plfe_family_{s}.csv", index=False)
        print(f"   saved figures/plfe_by_grade_{s}.png + data/plfe_family_{s}.csv")
    print("ANALYZE DONE")


if __name__ == "__main__":
    main()
