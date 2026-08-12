"""Picker v0 — a calibrated per-window LFE-probability model (LFE / EQ / NOISE).

This is the feature-based precursor to the EQTransformer-style per-sample picker: it takes a
40 s window's features and returns calibrated P(LFE), P(EQ), P(NOISE). Trained on Lin LFEs +
ANSS earthquakes + random noise at a station. Reports CV multi-class metrics, per-class AUC,
confusion matrix, feature importances; saves the model bundle for reuse / scanning.

Usage: PYTHONPATH=src python lfe_features/build_picker_v0.py --sta B011
"""
import argparse, os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report
import joblib
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
HERE = os.path.dirname(__file__); POL = ["hv_ratio", "rectilinearity", "planarity"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--sta", required=True); a = ap.parse_args(); s = a.sta.lower()
    ref = pd.read_parquet(f"{HERE}/data/feat_ref_{s}.parquet")
    eq = pd.read_parquet(f"{HERE}/data/feat_eq_{s}.parquet")
    df = pd.concat([ref, eq], ignore_index=True)
    df["cls"] = df.label.map({"LIN": "LFE", "EQ": "EQ", "RAND": "NOISE"})
    vert = [c for c in FEATURE_ORDER if c in df.columns and c not in POL]
    pol = [c for c in POL if c in df.columns and df[c].notna().mean() > 0.9]
    cols = vert + pol
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)
    X = np.asarray(df[cols].to_numpy(), dtype=float)      # de-arrow
    yc = np.asarray(df.cls.tolist())
    print(f"=== picker v0 [{a.sta}] : " + ", ".join(f"{k}={v}" for k, v in df.cls.value_counts().items()) +
          f" | {len(cols)} features ===")

    scaler = StandardScaler().fit(X); Xs = scaler.transform(X)
    base = RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)

    # cross-validated probabilities for honest metrics
    proba = cross_val_predict(base, Xs, yc, cv=cv, method="predict_proba", n_jobs=-1)
    classes = base.fit(Xs, yc).classes_
    pred = classes[proba.argmax(1)]
    print("\nper-class one-vs-rest CV-AUC:")
    for i, c in enumerate(classes):
        print(f"   {c:6s} {roc_auc_score((yc == c).astype(int), proba[:, i]):.3f}")
    print("\nconfusion matrix (rows=true, cols=pred) order=%s:" % list(classes))
    print(confusion_matrix(yc, pred, labels=classes))
    print("\n" + classification_report(yc, pred, labels=classes, digits=3))

    # feature importances
    imp = sorted(zip(cols, base.feature_importances_), key=lambda t: -t[1])
    print("top feature importances:")
    for c, v in imp[:12]:
        print(f"   {c:14s} {v:.3f}")

    # final calibrated model on all data, saved
    cal = CalibratedClassifierCV(RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0),
                                 method="isotonic", cv=4).fit(Xs, yc)
    os.makedirs(f"{HERE}/models", exist_ok=True)
    bundle = dict(model=cal, scaler=scaler, cols=cols, classes=list(cal.classes_), station=a.sta)
    joblib.dump(bundle, f"{HERE}/models/picker_v0_{s}.joblib")
    meta = dict(station=a.sta, n=len(df), features=cols, classes=list(cal.classes_),
                cv_auc={c: float(roc_auc_score((yc == c).astype(int), proba[:, i])) for i, c in enumerate(classes)})
    json.dump(meta, open(f"{HERE}/models/picker_v0_{s}.json", "w"), indent=2)
    print(f"\nsaved models/picker_v0_{s}.joblib (+ .json)")

    # P(LFE) histogram by true class
    ip = list(cal.classes_).index("LFE")
    pl = cal.predict_proba(Xs)[:, ip]
    fig, ax = plt.subplots(figsize=(8, 4))
    for c, col in [("LFE", "#1a9850"), ("EQ", "#7570b3"), ("NOISE", "#999999")]:
        ax.hist(pl[df.cls.values == c], bins=40, range=(0, 1), alpha=0.5, density=True, label=c, color=col)
    ax.set_xlabel("calibrated P(LFE)"); ax.set_ylabel("density"); ax.legend()
    ax.set_title(f"{a.sta} picker v0 — P(LFE) by true class")
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/picker_v0_plfe_{s}.png", dpi=120)
    print(f"saved figures/picker_v0_plfe_{s}.png")
    print("PICKER V0 DONE")


if __name__ == "__main__":
    main()
