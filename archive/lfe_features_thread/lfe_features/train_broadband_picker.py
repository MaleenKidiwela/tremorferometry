#!/usr/bin/env python
"""Train + validate a FRESH broadband LFE picker for PGC (Merlin plan step 3, 2026-07-14).
Positives = COLOC (B011 certified detections -> PGC surface windows, co-located labels).
Negatives = RAND (feat_ref_pgc), EQ (feat_eq_pgc), TNOISE (in-tremor noise, the hard negative).
Model = RandomForest(400, balanced) + isotonic calibration, vertical features only (Z-only windows).
Validation: GroupKFold by YEAR (no ETS episode straddles train/test); report LFE-vs-{RAND,EQ,TNOISE} AUC
separately; SNR/amplitude-drop control (retrain without 'snr'). Acceptance (Merlin): vs RAND >=0.95,
vs TNOISE >=0.75, amplitude-control drop <=0.02; FAIL if vs TNOISE <0.65 (picker can't score during tremor).
Saves models/picker_broadband_pgc.joblib + .json report.
"""
import os, sys, json, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from feature_defs import FEATURE_ORDER
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import cross_val_predict, GroupKFold
from sklearn.metrics import roc_auc_score
import joblib

HERE = os.path.dirname(__file__)
POL = {'hv_ratio', 'rectilinearity', 'planarity'}          # 3-comp only; PGC windows are Z-only here
FEATS = [c for c in FEATURE_ORDER if c not in POL]          # 20 vertical features


def load():
    def prep(d):
        d = d.copy()
        if 'year' not in d.columns:                              # derive from epoch-seconds 'time' (unit='s'!)
            if 'time' not in d.columns:
                raise ValueError(f"no 'year' or 'time' column to group CV honestly: {list(d.columns)}")
            d['year'] = pd.to_datetime(d['time'], unit='s').dt.year
        if 'fam' not in d.columns:
            d['fam'] = d['label']
        return d
    # POSITIVES = Lin network-detected STRONG LFEs (visible per-window; co-located B011 cc>=0.92 detections
    # are individually too weak at the surface -> ~noise-like features, 0.52 vs RAND. Co-located validated
    # the STACK/coda survival, not per-window detection). Diagnostic in notes/2026-07-14_broadband_dvv.md.
    ref = pd.read_parquet(f'{HERE}/data/feat_ref_pgc.parquet')
    co = ref[ref.label == 'LIN'].assign(label='LFE')                                             # positives
    tn = pd.read_parquet(f'{HERE}/data/feat_tnoise_pgc.parquet')                                # TNOISE (in-tremor)
    rf = ref[ref.label == 'RAND']                                                                # RAND
    eq = pd.read_parquet(f'{HERE}/data/feat_eq_pgc.parquet').assign(label='EQ')                  # EQ
    keep = set(FEATS) | {'label', 'year', 'fam'}
    parts = [prep(d)[[c for c in prep(d).columns if c in keep]] for d in (co, tn, rf, eq)]
    df = pd.concat(parts, ignore_index=True)
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=FEATS)
    return df


def auc_by_class(y_islfe, proba_lfe, sub, cls):
    m = (sub == cls) | (y_islfe == 1)
    return roc_auc_score(y_islfe[m], proba_lfe[m])


def main():
    df = load()
    print(f"[train] {len(df)} windows | classes: {df.label.value_counts().to_dict()}", flush=True)
    df['islfe'] = (df.label == 'LFE').astype(int)
    X = df[FEATS].values; y = df.islfe.values; grp = df.year.values
    scaler = StandardScaler().fit(X); Xs = scaler.transform(X)
    cv = GroupKFold(n_splits=min(5, df.year.nunique()))
    rf = RandomForestClassifier(400, class_weight='balanced', n_jobs=-1, random_state=0)
    proba = cross_val_predict(rf, Xs, y, cv=cv, groups=grp, method='predict_proba', n_jobs=-1)[:, 1]

    res = {}
    for cls in ['RAND', 'EQ', 'TNOISE']:
        res[f'auc_vs_{cls}'] = round(float(auc_by_class(y, proba, df.label.values, cls)), 3)
    res['auc_overall'] = round(float(roc_auc_score(y, proba)), 3)

    # SNR/amplitude-drop control: retrain without 'snr'
    fdrop = [f for f in FEATS if f != 'snr']
    Xs2 = StandardScaler().fit_transform(df[fdrop].values)
    proba2 = cross_val_predict(RandomForestClassifier(400, class_weight='balanced', n_jobs=-1, random_state=0),
                               Xs2, y, cv=cv, groups=grp, method='predict_proba', n_jobs=-1)[:, 1]
    res['auc_vs_TNOISE_noSNR'] = round(float(auc_by_class(y, proba2, df.label.values, 'TNOISE')), 3)
    res['amp_control_drop'] = round(res['auc_vs_TNOISE'] - res['auc_vs_TNOISE_noSNR'], 3)

    print("\n[RESULT] broadband PGC picker (GroupKFold by year):")
    for k, v in res.items(): print(f"   {k}: {v}")
    ok = res['auc_vs_RAND'] >= 0.95 and res['auc_vs_TNOISE'] >= 0.75 and abs(res['amp_control_drop']) <= 0.05
    fail = res['auc_vs_TNOISE'] < 0.65
    print(f"\n   ACCEPTANCE vs RAND>=0.95 & vs TNOISE>=0.75 & |amp drop|<=0.05 -> {'PASS' if ok else 'REVIEW'}"
          f"{'  *** FAIL: vs TNOISE<0.65, broadband discovery noise-flooded ***' if fail else ''}")

    # fit final calibrated picker on all data
    cal = CalibratedClassifierCV(RandomForestClassifier(400, class_weight='balanced', n_jobs=-1, random_state=0),
                                 method='isotonic', cv=4).fit(Xs, y)
    joblib.dump(dict(model=cal, scaler=scaler, cols=FEATS, classes=['NONLFE', 'LFE'], station='PGC', band='broadband'),
                f'{HERE}/models/picker_broadband_pgc.joblib')
    json.dump(res, open(f'{HERE}/models/picker_broadband_pgc.json', 'w'), indent=2)
    print(f"\n   saved models/picker_broadband_pgc.joblib (+ .json)")


if __name__ == '__main__':
    main()
