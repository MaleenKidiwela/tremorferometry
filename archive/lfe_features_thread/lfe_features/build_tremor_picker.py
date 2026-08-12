"""THE GOAL: pick LFEs inside PNSN tremor windows while rejecting earthquakes / blasts /
cultural+ambient (in-tremor) noise.

Classes: LFE (curated Lin) vs EQ (ANSS earthquake) vs BLAST (ANSS explosion) vs TNOISE
(random in-tremor-window times, not near an LFE = the ambient/cultural the discovery detector grabs).
Trains a calibrated classifier, reports LFE-vs-each-contaminant AUC, the pooled LFE-vs-all AUC, a
precision/recall curve, and a base-rate-aware operating point (how pure is the kept-LFE set if the
candidate pool is only k% LFE). Saves the model for use as a discovery filter.

Usage: PYTHONPATH=src python lfe_features/build_tremor_picker.py --sta B011
"""
import argparse, os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, confusion_matrix, precision_recall_curve
import joblib
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
HERE = os.path.dirname(__file__); POL = ["hv_ratio", "rectilinearity", "planarity"]
SRC = {"LFE": ("feat_ref", "LIN"), "EQ": ("feat_eq", "EQ"), "BLAST": ("feat_blast", "BLAST"),
       "TNOISE": ("feat_tnoise", "TNOISE")}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--sta", required=True); a = ap.parse_args(); s = a.sta.lower()
    frames = []
    for cls, (pre, lab) in SRC.items():
        p = f"{HERE}/data/{pre}_{s}.parquet"
        d = pd.read_parquet(p)
        d = d[d.label == lab] if "label" in d.columns else d
        d = d.copy(); d["cls"] = cls; frames.append(d)
    df = pd.concat(frames, ignore_index=True)
    vert = [c for c in FEATURE_ORDER if c in df.columns and c not in POL]
    pol = [c for c in POL if c in df.columns and df[c].notna().mean() > 0.9]
    cols = vert + pol
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols)
    print(f"=== {a.sta} tremor-window picker: " + ", ".join(f"{k}={v}" for k, v in df.cls.value_counts().items()) +
          f" | {len(cols)} feats ===")
    X = np.asarray(df[cols].to_numpy(), float); yc = np.asarray(df.cls.tolist())
    Xs = StandardScaler().fit(X).transform(X)
    cv = StratifiedKFold(5, shuffle=True, random_state=0)
    proba = cross_val_predict(RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0),
                              Xs, yc, cv=cv, method="predict_proba", n_jobs=-1)
    classes = sorted(set(yc))                       # cross_val_predict columns follow sorted classes
    iL = classes.index("LFE"); pL = proba[:, iL]

    print("\nLFE-vs-each-contaminant AUC (can we tell LFEs from each confuser?):")
    for c in ["EQ", "BLAST", "TNOISE"]:
        m = np.isin(yc, ["LFE", c])
        print(f"   LFE vs {c:7s}: {roc_auc_score((yc[m]=='LFE').astype(int), pL[m]):.3f}")
    aucALL = roc_auc_score((yc == "LFE").astype(int), pL)
    print(f"   LFE vs ALL contaminants pooled: {aucALL:.3f}")

    print("\nconfusion (rows true, cols pred), order=%s:" % classes)
    pred = np.array(classes)[proba.argmax(1)]
    print(confusion_matrix(yc, pred, labels=classes))

    # precision/recall for LFE vs pooled contaminants
    yb = (yc == "LFE").astype(int)
    prec, rec, thr = precision_recall_curve(yb, pL)
    def prec_at(r):
        i = np.argmin(np.abs(rec - r)); return prec[i], (thr[min(i, len(thr)-1)])
    print("\nLFE detection precision/recall (negatives = EQ+BLAST+TNOISE, equal-sampled):")
    for r in (0.9, 0.8, 0.7, 0.6):
        p, t = prec_at(r); print(f"   recall {r:.1f} -> precision {p:.3f}  (P(LFE) thr {t:.2f})")

    # base-rate-aware: if candidate pool is k% LFE, purity of kept set at recall 0.8
    p80, t80 = prec_at(0.8)
    # at threshold t80, contaminant pass-rate = fraction of contaminants above thr
    cont_pass = (pL[yc != "LFE"] >= t80).mean()
    print(f"\nOperating point (P(LFE)>={t80:.2f}, LFE recall 0.8, contaminant pass-rate {cont_pass:.3f}):")
    for k in (0.02, 0.05, 0.10, 0.20):
        purity = (k * 0.8) / (k * 0.8 + (1 - k) * cont_pass + 1e-9)
        print(f"   candidate pool {int(k*100)}% LFE -> kept-set PURITY {purity:.2f}  (enrichment x{purity/k:.1f})")

    # save calibrated model
    cal = CalibratedClassifierCV(RandomForestClassifier(400, class_weight="balanced", n_jobs=-1, random_state=0),
                                 method="isotonic", cv=4).fit(Xs, yc)
    sc = StandardScaler().fit(X)
    joblib.dump(dict(model=cal, scaler=sc, cols=cols, classes=list(cal.classes_), station=a.sta),
                f"{HERE}/models/tremor_picker_{s}.joblib")
    json.dump(dict(station=a.sta, classes=list(cal.classes_), n=len(df), auc_LFE_vs_all=float(aucALL),
                   cont_pass_at_rec80=float(cont_pass)), open(f"{HERE}/models/tremor_picker_{s}.json", "w"), indent=2)

    # figure: P(LFE) by class + PR curve
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    col = {"LFE": "#1a9850", "EQ": "#7570b3", "BLAST": "#e6550d", "TNOISE": "#999999"}
    for c in classes:
        ax[0].hist(pL[yc == c], bins=40, range=(0, 1), alpha=0.5, density=True, label=c, color=col[c])
    ax[0].set_xlabel("P(LFE)"); ax[0].set_ylabel("density"); ax[0].legend(); ax[0].set_title(f"{a.sta}: P(LFE) by true class")
    ax[1].plot(rec, prec, lw=2, color="#1a9850"); ax[1].set_xlabel("LFE recall"); ax[1].set_ylabel("LFE precision")
    ax[1].set_ylim(0, 1.02); ax[1].set_title(f"LFE vs EQ+BLAST+TNOISE (AUC {aucALL:.3f})"); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/tremor_picker_{s}.png", dpi=120)
    print(f"\nsaved models/tremor_picker_{s}.joblib + figures/tremor_picker_{s}.png")
    print("TREMOR PICKER DONE")


if __name__ == "__main__":
    main()
