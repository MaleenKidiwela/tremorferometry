"""Validate experimental detections at a station never used by the detector.

Training learns a held-out-station template and its delay to the reference anchor.
Evaluation uses signed CC at the exact predicted lag (no test-window maximization).
Repeated within-day shifts retain event clustering; their empirical tail fraction
is a diagnostic, not a rigorously calibrated event-level false-alarm probability.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .matched_filter import align_to_anchor, normalized_cc
from .pilot import ROOT, date_start
from .provenance import sha256, write_json
from .reference import read_bostock
from .templates import build_station_template, cut
from .waveforms import load_day


def sample_series(series, times, start_s, fs):
    idx = (np.asarray(times) - start_s) * fs
    lo, hi = np.floor(idx).astype(int), np.ceil(idx).astype(int)
    valid = (lo >= 0) & (hi < len(series))
    out = np.full(len(idx), np.nan)
    ii = np.flatnonzero(valid)
    fraction = idx[ii] - lo[ii]
    out[ii] = series[lo[ii]] * (1 - fraction) + series[hi[ii]] * fraction
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--config', default='configs/heldout_station.yaml')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.output)
    cfg = json.loads((run / 'manifest.json').read_text())['config']
    spec = yaml.safe_load(Path(args.config).read_text())
    if spec['station'] in cfg['stations'] or spec['site_id'] in {s['site_id'] for s in cfg['stations'].values()}:
        raise ValueError('validation station/site was used by the detector')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    training = json.loads((run / 'training.json').read_text())
    anchors = training['anchors_from_published_s']
    reference = read_bostock((ROOT / cfg['reference_csv']).resolve())
    manifests = {}

    def load(date):
        d = pd.Timestamp(date)
        path = ROOT / cfg['raw_root'] / spec['station'] / str(d.year) / f'{d.dayofyear:03d}.mseed'
        wave, manifest = load_day(path, ROOT / 'cache' / 'waveforms', spec['channels'], spec['site_id'],
                                 date_start(date), fs=cfg['fs'])
        manifests[date] = manifest
        return wave

    windows = {family: [] for family in anchors}
    for date in cfg['train_dates']:
        print(f'[heldout training] {date}', flush=True)
        wave = load(date)
        r = reference[(reference.date == date) & reference.family_id.isin(anchors)]
        for family, group in r.groupby('family_id'):
            for time in group.time_s:
                seg = cut(wave, time - 10, int(70 * cfg['fs']))
                if seg is not None:
                    windows[family].append(seg)
    templates, qc = {}, []
    for family, segments in windows.items():
        template, row = build_station_template(segments, cfg['fs'], template_s=cfg['template_s'],
            min_events=spec['min_training_events'], min_split_cc=spec['min_training_split_cc'])
        row['family_id'] = family
        qc.append(row)
        if template is not None:
            templates[family] = (template, row['template_start_from_published_s'] - anchors[family])
    pd.DataFrame(qc).to_csv(out / 'training_qc.csv', index=False)
    # Freeze the template choice and validation rule before reading test candidates.
    write_json(out / 'validation_design.json', {'config': spec, 'config_sha256': sha256(args.config),
        'detector_manifest_sha256': sha256(run / 'manifest.json'), 'families_with_usable_template': list(templates),
        'source_sha256': sha256(__file__), 'detector_site_ids': [s['site_id'] for s in cfg['stations'].values()]})
    if not templates:
        write_json(out / 'report.json', {'status': 'inconclusive_no_coherent_heldout_template',
                                        'interpretation': 'This is not evidence against the detections.'})
        write_json(out / 'input_manifest.json', manifests)
        return
    candidates = pd.read_parquet(run / 'test_experimental_candidates.parquet')
    reference_test = pd.read_parquet(run / 'test_eligible_reference.parquet')
    values = {(family, kind): [] for family in templates for kind in ['candidates', 'reference_positive_control']}
    controls = {(family, kind): [[] for _ in range(spec['null_replicates'])] for family, kind in values}
    sample_rows = []
    for date in cfg['test_dates']:
        print(f'[heldout evaluation] {date}', flush=True)
        wave = load(date)
        rng = np.random.default_rng(np.random.SeedSequence([spec['seed'], int(date.replace('-', ''))]))
        shifts = rng.uniform(spec['shift_min_s'], spec['shift_max_s'], spec['null_replicates']) * rng.choice([-1, 1], spec['null_replicates'])
        for family, (template, offset) in templates.items():
            cc, valid = normalized_cc(wave.x, template, wave.valid)
            series, _ = align_to_anchor(cc, valid, offset, wave.fs, wave.x.shape[1])
            for kind, frame in [('candidates', candidates), ('reference_positive_control', reference_test)]:
                sub = frame[(frame.family_id == family) & (frame.date == date)]
                t = sub.time_s.to_numpy()
                obs = sample_series(series, t, wave.start_s, wave.fs)
                values[family, kind].extend(obs[np.isfinite(obs)].tolist())
                for row_time, val in zip(t, obs):
                    sample_rows.append(dict(family_id=family, kind=kind, date=date, time_s=row_time,
                                            signed_cc=val, valid=bool(np.isfinite(val))))
                for k, shift in enumerate(shifts):
                    null = sample_series(series, t + shift, wave.start_s, wave.fs)
                    good = np.isfinite(obs) & np.isfinite(null)
                    # Store paired effects, so data gaps cannot change the observed
                    # event subset without making the same change to the control.
                    controls[family, kind][k].extend((obs[good] - null[good]).tolist())
    # Preserve the same shift index across families when pooling. This checks for
    # aggregate unused-site signal; it does not certify every pooled candidate.
    for kind in ['candidates', 'reference_positive_control']:
        values['ALL_USABLE_GROUPS', kind] = [v for family in templates for v in values[family, kind]]
        controls['ALL_USABLE_GROUPS', kind] = [
            [v for family in templates for v in controls[family, kind][k]]
            for k in range(spec['null_replicates'])]
    results = []
    for key, obs in values.items():
        family, kind = key
        differences = [float(np.mean(x)) for x in controls[key] if len(x) >= spec['min_validation_events']]
        enough = len(obs) >= spec['min_validation_events'] and len(differences) >= spec['null_replicates'] * 0.8
        results.append({'family_id': family, 'kind': kind, 'n_valid_observations': len(obs),
                        'mean_signed_cc': float(np.mean(obs)) if obs else None,
                        'median_signed_cc': float(np.median(obs)) if obs else None,
                        'n_usable_shift_realizations': len(differences),
                        'paired_mean_cc_excess_median': float(np.median(differences)) if differences else None,
                        'paired_mean_cc_excess_range': [float(np.min(differences)), float(np.max(differences))] if differences else None,
                        'surrogate_nonpositive_fraction': (1 + sum(x <= 0 for x in differences)) / (1 + len(differences)) if differences else None,
                        'status': 'diagnostic_only' if enough else 'inconclusive_low_exposure'})
    pd.DataFrame(sample_rows).to_parquet(out / 'signed_cc_observations.parquet', index=False)
    write_json(out / 'report.json', {'status': 'heldout_station_diagnostic_not_precision_estimate', 'results': results,
        'limitations': ['Single unused site', 'Reference source-group identity and anchor time inherited from training only',
                       'No test-lag maximization', 'Shift tail fraction is not a calibrated event-level probability',
                       'Positive aggregate signal does not validate every event or every family assignment']})
    write_json(out / 'input_manifest.json', manifests)


if __name__ == '__main__':
    main()
