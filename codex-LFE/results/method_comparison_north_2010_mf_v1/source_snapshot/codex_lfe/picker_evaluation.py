"""Score arrival picks only inside explicitly reviewed, usable windows.

An unreviewed/ambiguous window is never a negative example. Scores on enriched
review panels are conditional on that sampling scheme, not operational precision.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .evaluate import match_events
from .provenance import sha256, write_json

BASIS = 'station_phase_arrival_utc'
WINDOW_COLUMNS = ['window_id', 'station', 'start_s', 'end_s', 'split', 'sampling', 'review_status', 'reviewer']
LABEL_COLUMNS = ['label_id', 'window_id', 'phase', 'time_s', 'uncertainty_s', 'event_class']


def validate_review(windows, labels):
    if not set(WINDOW_COLUMNS) <= set(windows) or not set(LABEL_COLUMNS) <= set(labels):
        raise ValueError('missing review schema columns')
    if windows.window_id.duplicated().any() or labels.label_id.duplicated().any():
        raise ValueError('duplicate review IDs')
    if not np.isfinite(windows[['start_s', 'end_s']].to_numpy(float)).all() or (windows.end_s <= windows.start_s).any():
        raise ValueError('invalid review intervals')
    if not windows.review_status.isin(['pending', 'complete', 'ambiguous', 'unusable']).all():
        raise ValueError('unknown review status')
    if not windows.sampling.isin(['random', 'reference_enriched', 'picker_enriched']).all():
        raise ValueError('unknown sampling scheme')
    if not windows.split.isin(['development', 'test']).all():
        raise ValueError('invalid split')
    if not windows.station.astype(str).str.fullmatch(r'[A-Z0-9]+\.[A-Z0-9]+').all():
        raise ValueError('full network.station required')
    complete = windows.review_status == 'complete'
    if windows.loc[complete, 'reviewer'].fillna('').astype(str).str.strip().eq('').any():
        raise ValueError('complete review needs reviewer identity')
    # Overlap would double count station exposure and leak observations across splits.
    for _, group in windows.groupby('station'):
        g = group.sort_values('start_s')
        if np.any(g.start_s.to_numpy()[1:] < g.end_s.to_numpy()[:-1]):
            raise ValueError('overlapping station review windows')
    if not labels.phase.isin(['P', 'S']).all() or not labels.event_class.isin(['lfe', 'ordinary_earthquake']).all():
        raise ValueError('invalid arrival phase or event class')
    if not np.isfinite(labels[['time_s', 'uncertainty_s']].to_numpy(float)).all() or (labels.uncertainty_s < 0).any():
        raise ValueError('invalid label time/uncertainty')
    by_id = windows.set_index('window_id')
    for row in labels.itertuples():
        if row.window_id not in by_id.index:
            raise ValueError('label has unknown window')
        w = by_id.loc[row.window_id]
        if not w.start_s <= row.time_s < w.end_s:
            raise ValueError('label outside review window')


def score_picks(picks, windows, labels, threshold, tolerance_s, split, coverage):
    validate_review(windows, labels)
    if not np.isfinite([threshold, tolerance_s]).all() or tolerance_s < 0:
        raise ValueError('finite threshold and nonnegative tolerance required')
    if split not in ('development', 'test'):
        raise ValueError('invalid split')
    if not {'station', 'phase', 'time_s', 'model_score', 'time_basis'} <= set(picks):
        raise ValueError('missing prediction schema columns')
    if not picks.phase.isin(['P', 'S']).all() or not picks.time_basis.eq(BASIS).all():
        raise ValueError('invalid predicted phase or time basis')
    if not np.isfinite(picks[['time_s', 'model_score']].to_numpy(float)).all():
        raise ValueError('nonfinite prediction')
    if not picks.station.astype(str).str.fullmatch(r'[A-Z0-9]+\.[A-Z0-9]+').all():
        raise ValueError('full network.station required for predictions')
    selected = windows[windows.split == split]
    if not {'station', 'phase', 'start_s', 'end_s'} <= set(coverage):
        raise ValueError('explicit prediction coverage intervals required')
    if not coverage.phase.isin(['P', 'S']).all() or not np.isfinite(coverage[['start_s', 'end_s']].to_numpy(float)).all() or (coverage.end_s <= coverage.start_s).any():
        raise ValueError('invalid prediction coverage')
    covered = []
    for w in selected.itertuples():
        both = True
        for phase in ['P', 'S']:
            frontier = w.start_s
            intervals = coverage[(coverage.station == w.station) & (coverage.phase == phase)].sort_values('start_s')
            for interval in intervals.itertuples():
                if interval.start_s > frontier:
                    break
                frontier = max(frontier, interval.end_s)
            both &= frontier >= w.end_s
        if both:
            covered.append(w.window_id)
    # Uncertain arrival labels invalidate the window for this tolerance. Otherwise
    # excluding just those labels would turn nearby predictions into false picks.
    wide = set(labels.loc[labels.uncertainty_s > tolerance_s, 'window_id'])
    usable = selected[(selected.review_status == 'complete') & ~selected.window_id.isin(wide) & selected.window_id.isin(covered)]
    if usable.empty:
        raise ValueError('no completely reviewed windows with adequate label certainty')
    p = picks[picks.model_score >= threshold].copy().reset_index(drop=True)
    p['prediction_index'] = np.arange(len(p))
    details = []
    for w in usable.itertuples():
        sub = p[(p.station == w.station) & (p.time_s >= w.start_s) & (p.time_s < w.end_s)]
        truth = labels[(labels.window_id == w.window_id) & (labels.event_class == 'lfe')]
        for phase in ['P', 'S']:
            pred = sub[sub.phase == phase].copy().reset_index(drop=True)
            ref = truth[truth.phase == phase].copy().reset_index(drop=True)
            # Group by station AND phase AND window; no cross-station credit.
            for frame in [pred, ref]:
                frame['family_id'] = phase
                frame['time_basis'] = BASIS
            matches = match_events(pred, ref, tolerance_s=tolerance_s, same_family=True)
            details.append({'window_id': w.window_id, 'station': w.station, 'sampling': w.sampling,
                'phase': phase, 'reviewed_s': w.end_s - w.start_s,
                'n_predictions': len(pred), 'n_reference_arrivals': len(ref), 'n_matched': len(matches),
                'residuals_s': matches.residual_s.tolist()})
    frame = pd.DataFrame(details)
    reports = []
    # Keep enriched and random panels separate; no unweighted operational rates.
    for (sampling, phase), g in frame.groupby(['sampling', 'phase']):
        npred, nref, matched = map(int, [g.n_predictions.sum(), g.n_reference_arrivals.sum(), g.n_matched.sum()])
        residuals = np.asarray([v for rs in g.residuals_s for v in rs])
        reports.append({'sampling': sampling, 'phase': phase, 'n_predictions': npred,
            'n_reviewed_lfe_arrivals': nref, 'n_matched': matched,
            'recall': matched / nref if nref else None,
            'review_panel_precision': matched / npred if npred else None,
            'unmatched_predictions': npred - matched,
            'reviewed_station_hours': float(g.reviewed_s.sum() / 3600),
            'unmatched_picks_per_reviewed_station_hour': (npred - matched) / (g.reviewed_s.sum() / 3600),
            'median_signed_residual_s': float(np.median(residuals)) if len(residuals) else None,
            'median_absolute_residual_s': float(np.median(abs(residuals))) if len(residuals) else None,
            'p90_absolute_residual_s': float(np.quantile(abs(residuals), .9)) if len(residuals) else None})
    return {'split': split, 'threshold': threshold, 'tolerance_s': tolerance_s,
        'n_windows_selected': len(selected), 'n_windows_scored': len(usable),
        'n_windows_missing_complete_prediction_coverage': int((~selected.window_id.isin(covered)).sum()),
        'n_windows_excluded_review_status': int((selected.review_status != 'complete').sum()),
        'n_complete_windows_excluded_label_uncertainty': int(((selected.review_status == 'complete') & selected.window_id.isin(wide)).sum()),
        'results': reports,
        'limitations': ['Scores apply to reviewed identifiable LFE arrivals, not all natural LFEs.',
            'Enriched-panel precision/rates are not operational prevalence estimates.',
            'Prediction timing residuals are relative to reviewed nominal times, which retain uncertainty.',
            'Matching maximizes cardinality; crowded arrivals can have nonunique timing assignments.',
            'No confidence intervals: station/event clustering must be handled in the final benchmark.',
            'No phase-confusion estimate yet; an incorrect phase contributes an unmatched pick and missed arrival.']}, details


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--picks', required=True)
    ap.add_argument('--review', required=True, help='exported JSON containing windows and labels')
    ap.add_argument('--coverage', required=True, help='Parquet: station,phase,start_s,end_s of valid prediction coverage')
    ap.add_argument('--threshold', required=True, type=float)
    ap.add_argument('--tolerance-s', required=True, type=float)
    ap.add_argument('--split', required=True, choices=['development', 'test'])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    review = json.loads(Path(args.review).read_text())
    windows = pd.DataFrame(review['windows'], columns=WINDOW_COLUMNS)
    labels = pd.DataFrame(review['labels'], columns=LABEL_COLUMNS)
    report, rows = score_picks(pd.read_parquet(args.picks), windows, labels, args.threshold, args.tolerance_s, args.split,
                               pd.read_parquet(args.coverage))
    out.mkdir(parents=True, exist_ok=True)
    report['provenance'] = {'picks_sha256': sha256(args.picks), 'coverage_sha256': sha256(args.coverage), 'review_sha256': sha256(args.review), 'source_sha256': sha256(__file__)}
    write_json(out / 'report.json', report)
    write_json(out / 'window_scores.json', rows)

if __name__ == '__main__':
    main()
