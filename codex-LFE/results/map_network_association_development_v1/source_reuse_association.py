"""Development association of saved E2 picks, restricted to the explicit roster.

Reuses picks and velocity table, not previous event assignments or event labels.
No GPU or waveform download. Outputs are candidates, not validated LFEs.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from .provenance import digest, environment, sha256, write_json
from .roster import read_roster


def prepare_picks(frame, allowed, threshold, start, end, boundary_guard=1.0):
    """Keep network identity, valid P/S peaks, and fixed UTC analysis interval."""
    frame = frame.copy()
    frame['station'] = frame.net + '.' + frame.sta
    frame['source_pick_id'] = np.arange(len(frame))
    frame['epoch'] = pd.to_datetime(frame.time, utc=True).dt.as_unit('ns').astype('int64') / 1e9
    member = frame.station.isin(allowed)
    selected = member & (frame.thr_first <= threshold) & (frame.prob >= threshold)
    selected &= (frame.epoch >= start) & (frame.epoch < end)
    frame = frame[selected].copy()
    if not frame.phase.isin(['P', 'S']).all():
        raise ValueError('invalid phase')
    if not np.isfinite(frame[['epoch', 'prob']].to_numpy()).all():
        raise ValueError('nonfinite time/probability')
    if not frame.prob.between(0, 1).all():
        raise ValueError('invalid probability')
    if frame.duplicated(['station', 'phase', 'epoch']).any():
        raise ValueError('duplicate station/phase/time picks')
    before = len(frame)
    sod = frame.epoch % 86400
    frame = frame[(sod >= boundary_guard) & (sod < 86400-boundary_guard)]
    out = pd.DataFrame({'station': frame.station, 'phase': frame.phase,
                        'time': frame.epoch, 'probability': frame.prob,
                        'source_pick_id': frame.source_pick_id})
    return out.sort_values('time').reset_index(drop=True), before-len(frame)


def shift_picks(picks, mode, seed, start, end):
    """Shift each station or station-hour; retain pick identities and count.

    Use a fixed hour axis for reproducibility across thresholds. Event scoring
    excludes outer 3120 seconds, shared by observed and controls. This does not
    fix internal waveform gaps or guarantee exchangeable null data.
    """
    out = picks.copy()
    if mode == 'observed':
        return out, {}
    if mode not in ('station', 'station_hour'):
        raise ValueError('invalid control mode')
    rng = np.random.default_rng(seed)
    sites = sorted(out.station.unique())
    si = out.station.map({s: i for i, s in enumerate(sites)}).to_numpy()
    shape = (len(sites),) if mode == 'station' else (len(sites), int(np.ceil((end-start)/3600)))
    offsets = rng.uniform(300, 3000, shape) * rng.choice([-1., 1.], shape)
    applied = offsets[si] if mode == 'station' else offsets[si, ((out.time-start)//3600).astype(int)]
    out['time'] += applied
    out['shift_s'] = applied
    return out.sort_values('time').reset_index(drop=True), dict(
        mode=mode, seed=seed, stations=sites, offsets_s=offsets.tolist())


def finish_events(events, assignments, minimum, start, end):
    if events.empty:
        return events, assignments
    counts = assignments.groupby('event_idx').agg(
        n_stations=('station', 'nunique'), n_picks=('pick_idx', 'count'),
        n_s_stations=('station', lambda x: x[assignments.loc[x.index, 'phase']=='S'].nunique()),
        rms_residual_s=('residual', lambda x: float(np.sqrt(np.mean(x*x)))))
    events = events.join(counts, on='idx')
    events = events[(events.n_stations >= minimum) & (events.n_s_stations >= minimum)
                    & (events.time >= start) & (events.time < end)].copy()
    assignments = assignments[assignments.event_idx.isin(events.idx)].copy()
    events['time_utc'] = pd.to_datetime(events.time, unit='s', utc=True)
    events['near_depth_boundary'] = (events.depth < 21) | (events.depth > 49)
    return events, assignments


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--picks', default='../lfe_catalog_repair/results/e2_picks_vi2010_cascadia.parquet')
    parser.add_argument('--metadata', default='../lfe_catalog_repair/results/e2_stations.csv')
    parser.add_argument('--velocity', default='../lfe_catalog_repair/results/e2_pyocto_crosson')
    parser.add_argument('--roster', default='results/target_network_from_map/stations.csv')
    parser.add_argument('--output', required=True)
    parser.add_argument('--thresholds', default='0.3,0.5')
    parser.add_argument('--minimum-stations', default='4,5')
    parser.add_argument('--seeds', default='1,2,3')
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()
    import pyocto
    start = pd.Timestamp('2010-08-10', tz='UTC').timestamp()
    end = pd.Timestamp('2010-09-14', tz='UTC').timestamp()
    guard = 3120.0
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    sources = {str(Path(p).resolve()): sha256(p) for p in
               [args.picks, args.metadata, args.velocity, args.roster, __file__]}
    config = {'sources': sources, 'arguments': vars(args), 'start': start, 'end': end,
              'event_guard_s': guard, 'utc_day_pick_guard_s': 1.0,
              'lat': [47.7, 49.6], 'lon': [-125.4, -122.3], 'zlim': [20., 50.],
              'pick_match_tolerance': 1.2, 'edt_pick_std': 0.5, 'exponential_edt': False,
              'pyocto': importlib.metadata.version('pyocto')}
    fingerprint = digest(config)
    manifest = out/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text())['fingerprint'] != fingerprint:
        raise ValueError('output exists with different configuration or source hashes')
    write_json(manifest, {'fingerprint': fingerprint, 'config': config, 'environment': environment(),
                         'status': 'development; reused picks have no per-sample inference coverage mask'})
    allowed = read_roster(args.roster)
    raw = pd.read_parquet(args.picks)
    if raw.model.nunique() != 1 or raw['set'].nunique() != 1:
        raise ValueError('one model and source set required')
    meta = pd.read_csv(args.metadata)
    meta['id'] = meta.net + '.' + meta.sta
    if meta.id.duplicated().any():
        raise ValueError('duplicate station metadata')
    summaries = []
    for threshold in map(float, args.thresholds.split(',')):
        base, removed = prepare_picks(raw, allowed, threshold, start, end)
        used = set(base.station)
        if not used or not used <= set(meta.id):
            raise ValueError('missing station metadata or no picks')
        st = meta[meta.id.isin(used)].rename(columns={'lat': 'latitude', 'lon': 'longitude', 'elev_m': 'elevation'})
        st = st[['id', 'latitude', 'longitude', 'elevation']]
        st.to_csv(out/'stations.csv', index=False)
        base.to_parquet(out/f'input_t{threshold}.parquet', index=False)
        modes = [('observed', 0)] + [(mode, seed) for mode in ['station', 'station_hour']
                                   for seed in map(int, args.seeds.split(','))]
        for minimum in map(int, args.minimum_stations.split(',')):
            for mode, seed in modes:
                name = f't{threshold}_n{minimum}_{mode}_{seed}'
                result_path = out/f'{name}.json'
                if result_path.exists():
                    saved = json.loads(result_path.read_text())
                    for file, checksum in saved['output_sha256'].items():
                        if sha256(out/file) != checksum:
                            raise ValueError('cached association output changed')
                    summaries.append(saved)
                    continue
                picks, shifts = shift_picks(base, mode, seed, start, end)
                write_json(out/f'{name}_shifts.json', shifts)
                vm = pyocto.VelocityModel1D(args.velocity, tolerance=2., association_cutoff_distance=250.)
                assoc = pyocto.OctoAssociator.from_area(
                    lat=config['lat'], lon=config['lon'], zlim=config['zlim'], velocity_model=vm,
                    time_before=120., n_picks=minimum, n_s_picks=minimum, n_p_picks=0,
                    n_p_and_s_picks=0, pick_match_tolerance=1.2, edt_pick_std=.5,
                    exponential_edt=False, min_interevent_time=3., n_threads=args.threads)
                stations = st.copy()
                assoc.transform_stations(stations)
                tick = time.monotonic()
                print(f'RUN {name} {len(picks)} picks, {len(stations)} stations', flush=True)
                events, assignments = assoc.associate(picks, stations)
                if len(events):
                    assoc.transform_events(events)
                    assignments['source_pick_id'] = picks.iloc[assignments.pick_idx.to_numpy()].source_pick_id.to_numpy()
                events, assignments = finish_events(events, assignments, minimum, start+guard, end-guard)
                files = [f'{name}_events.parquet', f'{name}_assignments.parquet']
                events.to_parquet(out/files[0], index=False)
                assignments.to_parquet(out/files[1], index=False)
                summary = {'name': name, 'threshold': threshold, 'minimum_stations': minimum,
                           'mode': mode, 'seed': seed, 'n_events': len(events),
                           'n_input_picks': len(picks), 'boundary_picks_removed': removed,
                           'near_depth_boundary': int(events.near_depth_boundary.sum()) if len(events) else 0,
                           'seconds': time.monotonic()-tick,
                           'output_sha256': {f: sha256(out/f) for f in files}}
                write_json(result_path, summary)
                summaries.append(summary)
                print(f'DONE {name}: {len(events)} candidates in {summary["seconds"]:.1f}s', flush=True)
    frame = pd.DataFrame(summaries).drop(columns='output_sha256')
    frame.to_csv(out/'runs.csv', index=False)
    rows = []
    for (threshold, minimum), group in frame.groupby(['threshold', 'minimum_stations']):
        observed = int(group[group['mode']=='observed'].n_events.iloc[0])
        for mode in ['station', 'station_hour']:
            counts = group[group['mode']==mode].n_events
            rows.append({'threshold': threshold, 'minimum_stations': minimum, 'mode': mode,
                         'observed': observed, 'control_counts': counts.tolist(),
                         'mean_control_over_observed': float(counts.mean()/observed) if observed else None})
    write_json(out/'comparison.json', rows)
    print(json.dumps(rows, indent=2), flush=True)


if __name__ == '__main__':
    main()
