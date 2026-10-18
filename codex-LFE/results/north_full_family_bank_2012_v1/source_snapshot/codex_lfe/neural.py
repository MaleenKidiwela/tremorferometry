"""Local LFEDetect inference for new-source discovery, never an event catalog.

The published model runs with its own batch normalization. Our explicit 20-Hz
preprocessing and central-window stitching are recorded methodological choices,
not a claim to reproduce the authors' full processing chain. All weights are local;
no model downloads or waveform fetches occur in this runner.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from .provenance import environment, sha256, write_json
from .waveforms import Waveform, load_day, runs

ROOT = Path(__file__).resolve().parents[1]


def orientation_matrix(orientations):
    """Rows project geographic [up, north, east] onto each sensor's positive axis."""
    az = np.deg2rad([o['azimuth'] for o in orientations])
    dip = np.deg2rad([o['dip'] for o in orientations])
    matrix = np.column_stack([-np.sin(dip), np.cos(dip) * np.cos(az), np.cos(dip) * np.sin(az)])
    if np.linalg.cond(matrix) > 10:
        raise ValueError('ill-conditioned sensor orientation')
    return matrix


def verified_zne(wave, inventory):
    """Verify metadata epochs, component response compatibility, then rotate.

    Full frequency-dependent response removal is NOT performed. This permits
    within-instrument gain equalization only when normalized component responses
    agree within 5% throughout 1-8 Hz. Other instruments must be response-corrected
    first rather than silently fed to the model in distorted coordinates.
    """
    from obspy import UTCDateTime
    begin = UTCDateTime(wave.start_s)
    end = UTCDateTime(wave.start_s + (wave.x.shape[1] - 1) / wave.fs)
    orientations, gains, responses = [], [], []
    for channel in wave.channels:
        a = inventory.get_orientation(channel, begin)
        b = inventory.get_orientation(channel, end)
        if a != b:
            raise ValueError(f'orientation changes within requested epoch: {channel}')
        ra, rb = inventory.get_response(channel, begin), inventory.get_response(channel, end)
        if str(ra) != str(rb):
            raise ValueError(f'response epoch changes: {channel}')
        gain = float(ra.instrument_sensitivity.value)
        if gain <= 0 or not np.isfinite(gain):
            raise ValueError(f'invalid instrument sensitivity: {channel}')
        orientations.append(a)
        gains.append(gain)
        responses.append(ra.get_evalresp_response_for_frequencies(np.linspace(1, 8, 16), output='VEL') / gain)
    responses = np.asarray(responses)
    if np.max(np.abs(responses / responses[0] - 1)) > 0.05:
        raise ValueError('component responses differ: full response correction required before rotation')
    matrix = orientation_matrix(orientations)
    corrected = wave.x * (np.median(gains) / np.asarray(gains))[:, None]
    zne = np.linalg.solve(matrix, corrected).astype(np.float32)
    valid = np.broadcast_to(wave.valid.all(axis=0), zne.shape).copy()
    zne[~valid] = 0
    return Waveform(zne, valid, wave.start_s, wave.fs, ('Z', 'N', 'E'), wave.site_id), {
        'native_channels': list(wave.channels), 'orientations': orientations, 'sensitivities': gains,
        'rotation_matrix': matrix.tolist(), 'full_response_removed': False,
        'component_response_agreement_tolerance': 0.05, 'orientation_verified': True}


def infer_probabilities(wave, model, device='cpu', batch_size=64):
    """60-s windows, 30-s step, central 30-s output; windows touching gaps rejected."""
    import torch
    if wave.channels != ('Z', 'N', 'E') or wave.fs != 20:
        raise ValueError('verified ZNE at 20 Hz required')
    length = int(model.in_samples)
    if length != 1200 or float(model.sampling_rate) != wave.fs:
        raise ValueError('runner expects published 60-second/20-Hz LFEDetect model')
    step, edge = length // 2, length // 4
    starts = [i for i in range(0, wave.x.shape[1] - length + 1, step)
              if wave.valid[:, i:i + length].all() and np.max(np.abs(wave.x[:, i:i + length])) > 0]
    probabilities = np.full((2, wave.x.shape[1]), np.nan, np.float32)
    model = model.to(device).eval()
    with torch.inference_mode():
        for first in range(0, len(starts), batch_size):
            ids = starts[first:first + batch_size]
            batch = torch.as_tensor(np.stack([wave.x[:, i:i + length] for i in ids]), device=device)
            batch = model.annotate_batch_pre(batch, {})
            result = model(batch).cpu().numpy()
            if result.shape != (len(ids), 2, length):
                raise ValueError(f'unexpected model output shape: {result.shape}')
            for i, values in zip(ids, result):
                probabilities[:, i + edge:i + length - edge] = values[:, edge:length - edge]
    return probabilities, len(starts)


def picks_from_probabilities(probabilities, start_s, fs, station, threshold=0.3):
    rows = []
    for phase, values in zip(['P', 'S'], probabilities):
        finite = np.isfinite(values)
        peaks, _ = find_peaks(np.where(finite, values, -np.inf), height=threshold, distance=int(fs))
        for i in peaks:
            # Reject peaks abutting a missing prediction interval (edge artifacts).
            lo, hi = max(0, i - int(fs)), min(len(values), i + int(fs) + 1)
            if hi - lo != 2 * int(fs) + 1 or not finite[lo:hi].all():
                continue
            rows.append({'station': station, 'phase': phase, 'time_s': start_s + i / fs,
                         'model_score': float(values[i]), 'time_basis': 'station_phase_arrival_utc',
                         'quality_status': 'unassociated_model_candidate'})
    return pd.DataFrame(rows, columns=['station', 'phase', 'time_s', 'model_score', 'time_basis', 'quality_status'])


def prediction_coverage(probabilities, start_s, fs, station):
    """Intervals where the picker's complete +/-1 second edge check is valid."""
    rows = []
    edge = int(fs)
    for phase, values in zip(['P', 'S'], probabilities):
        for lo, hi in runs(np.isfinite(values)):
            lo, hi = lo + edge, hi - edge
            if hi > lo:
                rows.append({'station': station, 'phase': phase,
                             'start_s': start_s + lo / fs, 'end_s': start_s + hi / fs})
    return pd.DataFrame(rows, columns=['station', 'phase', 'start_s', 'end_s'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', required=True)
    ap.add_argument('--inventory', required=True)
    ap.add_argument('--channels', required=True, help='comma-separated FULL SEED IDs')
    ap.add_argument('--site', required=True)
    ap.add_argument('--date', required=True)
    ap.add_argument('--model', default='../lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--batch-size', type=int, default=64)
    ap.add_argument('--duration-s', type=int, default=86400)
    ap.add_argument('--threshold', type=float, default=0.3)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('SEISBENCH_CACHE_ROOT', str(ROOT / 'cache' / 'seisbench'))
    import torch
    import seisbench.models as sbm
    from obspy import read_inventory
    torch.set_num_threads(4)
    if args.device != 'cpu' and not torch.cuda.is_available():
        raise RuntimeError('CUDA requested but unavailable')
    started = time.monotonic()
    channels = args.channels.split(',')
    wave, manifest = load_day(args.raw, ROOT / 'cache' / 'waveforms', channels, args.site,
                              pd.Timestamp(args.date, tz='UTC').timestamp(), duration_s=args.duration_s)
    wave, orientation = verified_zne(wave, read_inventory(args.inventory))
    model = sbm.LFEDetect.load(args.model, map_location='cpu', weights_only=True)
    scores, nw = infer_probabilities(wave, model, args.device, args.batch_size)
    picks = picks_from_probabilities(scores, wave.start_s, wave.fs, '.'.join(channels[0].split('.')[:2]), args.threshold)
    picks.to_parquet(out / 'unassociated_picks.parquet', index=False)
    prediction_coverage(scores, wave.start_s, wave.fs, '.'.join(channels[0].split('.')[:2])).to_parquet(
        out / 'prediction_coverage.parquet', index=False)
    np.savez(out / 'probabilities.npz', probabilities=scores, start_s=wave.start_s, fs=wave.fs)
    write_json(out / 'manifest.json', {'waveform': manifest, 'orientation': orientation,
        'inventory_sha256': sha256(args.inventory), 'model_config_sha256': sha256(args.model + '.json'),
        'weights_sha256': sha256(args.model + '.pt'), 'environment': environment(), 'device': args.device,
        'source_sha256': sha256(__file__), 'threshold': args.threshold, 'windows_processed': nw,
        'valid_prediction_seconds': float(np.isfinite(scores[0]).sum() / wave.fs), 'n_picks': len(picks),
        'elapsed_s': time.monotonic() - started, 'status': 'picks_only_not_events',
        'model_training': model.weights_docstring})
    print(f'{len(picks)} unassociated picks; {nw} windows; {time.monotonic() - started:.1f}s', flush=True)


if __name__ == '__main__':
    main()
