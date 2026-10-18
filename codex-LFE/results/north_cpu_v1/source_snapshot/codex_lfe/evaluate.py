"""Evaluate the FINAL rows, using one-to-one time/family matches.

Reference agreement is not precision against an exhaustive truth set. Missing reference
events remain unverified. All input times must share a declared physical time basis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.stats import beta


def match_events(detections, reference, tolerance_s=1.0, same_family=True):
    """Maximum-cardinality bipartite match, with stable sorted input.

    Returned indices are positional, not DataFrame labels. Within ambiguous clusters
    this does not claim a unique/minimum-residual physical assignment. An event can
    contribute at most once on either side. NaN family IDs are never family matches.
    """
    if not np.isfinite(tolerance_s) or tolerance_s < 0:
        raise ValueError("tolerance_s must be finite and nonnegative")
    for frame in [detections, reference]:
        if not {"time_s", "family_id", "time_basis"}.issubset(frame):
            raise ValueError("need time_s, family_id and time_basis")
        if not np.isfinite(frame.time_s.to_numpy(float)).all():
            raise ValueError("nonfinite event time")
    bases = set(detections.time_basis) | set(reference.time_basis)
    if len(bases) > 1 or any(not isinstance(b, str) or not b for b in bases):
        raise ValueError("incompatible or missing time bases")
    cols = ["detection_index", "reference_index", "residual_s"]
    if not len(detections) or not len(reference):
        return pd.DataFrame(columns=cols)
    dt = detections.time_s.to_numpy(float)
    rt = reference.time_s.to_numpy(float)
    order = np.argsort(rt, kind="stable")
    sorted_rt = rt[order]
    df = detections.family_id.to_numpy()
    rf = reference.family_id.to_numpy()
    indptr, indices = [0], []
    for i, time in enumerate(dt):
        lo, hi = np.searchsorted(sorted_rt, [time - tolerance_s, time + tolerance_s], side="left")
        hi = np.searchsorted(sorted_rt, time + tolerance_s, side="right")
        js = order[lo:hi]
        if same_family:
            js = js[(rf[js] == df[i]) & pd.notna(rf[js]) & pd.notna(df[i])]
        # Prefer closest edges; cardinality is guaranteed, minimum distance is not.
        js = js[np.argsort(np.abs(rt[js] - time), kind="stable")]
        indices.extend(js.tolist())
        indptr.append(len(indices))
    graph = csr_matrix((np.ones(len(indices), dtype=np.int8), indices, indptr),
                       shape=(len(dt), len(rt)))
    matched = maximum_bipartite_matching(graph, perm_type="column")
    ii = np.flatnonzero(matched >= 0)
    jj = matched[ii]
    return pd.DataFrame(dict(detection_index=ii, reference_index=jj, residual_s=dt[ii] - rt[jj]))


def binomial_interval(k, n, alpha=0.05):
    """Exact interval; describes binomial uncertainty, not event dependence."""
    if n == 0:
        return [None, None]
    return [0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1)),
            1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))]


def score_catalog(detections, reference, tolerance_s=1.0):
    own = match_events(detections, reference, tolerance_s, True)
    any_family = match_events(detections, reference, tolerance_s, False)
    return {
        "n_detections": len(detections), "n_reference": len(reference),
        "n_matched_same_family": len(own), "n_matched_any_family": len(any_family),
        "same_family_recall": len(own) / len(reference) if len(reference) else None,
        "any_family_recall": len(any_family) / len(reference) if len(reference) else None,
        "reference_agreement_fraction": len(any_family) / len(detections) if len(detections) else None,
        "same_family_recall_binomial_ci95": binomial_interval(len(own), len(reference)),
        "reference_is_exhaustive": False,
        "warning": "Reference agreement is not measured precision; CI ignores clustering.",
    }


def null_summary(observed_count, observed_exposure_s, null_counts, null_exposures_s):
    """Exposure-adjusted control/observed count ratios; not a calibrated FDR claim."""
    counts = np.asarray(null_counts, float)
    exp = np.asarray(null_exposures_s, float)
    if observed_exposure_s <= 0 or counts.shape != exp.shape or np.any(exp <= 0):
        raise ValueError("positive exposure and matching control arrays required")
    ratios = counts / exp * observed_exposure_s / observed_count if observed_count else None
    return {"observed_count": int(observed_count), "observed_exposure_s": float(observed_exposure_s),
            "n_controls": len(counts), "control_counts": counts.astype(int).tolist(),
            "control_exposure_s": exp.tolist(),
            "control_to_observed_mean": float(ratios.mean()) if ratios is not None and len(ratios) else None,
            "control_to_observed_range": [float(ratios.min()), float(ratios.max())]
            if ratios is not None and len(ratios) else [None, None],
            "interpretation": "Empirical control ratio, not an event probability or proven FDR."}
