"""Small FDSN metadata acquisitions and explicit triplet-epoch inventories.

This never requests waveforms. Metadata coverage is an upper bound on recording
availability; it cannot establish station noise, uptime or LFE detectability.
"""
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

from .neural import orientation_matrix
from .provenance import sha256, write_json


def triplet_epochs(inventory, start_s, end_s, min_fs=20):
    rows = []
    for network in inventory:
        for station in network:
            groups = {}
            for channel in station:
                if channel.sample_rate < min_fs:
                    continue
                groups.setdefault((channel.location_code, channel.code[:2]), []).append(channel)
            for (location, prefix), channels in groups.items():
                for selected in combinations(channels, 3):
                    if len({c.code for c in selected}) != 3:
                        continue
                    lo = max([start_s] + [float(c.start_date) for c in selected if c.start_date is not None])
                    hi = min([end_s] + [float(c.end_date) for c in selected if c.end_date is not None])
                    if hi <= lo or any(c.azimuth is None or c.dip is None for c in selected):
                        continue
                    try:
                        matrix = orientation_matrix([{'azimuth': c.azimuth, 'dip': c.dip} for c in selected])
                    except ValueError:
                        continue
                    ids = sorted(f'{network.code}.{station.code}.{location}.{c.code}' for c in selected)
                    rows.append({'network': network.code, 'station': station.code, 'location': location,
                        'band_instrument': prefix, 'channels': ','.join(ids),
                        'latitude': float(station.latitude), 'longitude': float(station.longitude),
                        'elevation_m': float(station.elevation), 'start_utc': pd.Timestamp(lo, unit='s', tz='UTC').isoformat(),
                        'end_utc_exclusive': pd.Timestamp(hi, unit='s', tz='UTC').isoformat(),
                        'duration_days_metadata': (hi - lo) / 86400,
                        'min_sampling_rate_hz': min(float(c.sample_rate) for c in selected),
                        'orientation_condition': float(np.linalg.cond(matrix)),
                        'uncompressed_float32_bytes': int((hi - lo) * sum(float(c.sample_rate) for c in selected) * 4),
                        'waveform_availability_verified': False})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bounds', required=True, help='minlat,maxlat,minlon,maxlon')
    ap.add_argument('--start', required=True)
    ap.add_argument('--end', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--endpoint', default='https://service.earthscope.org/fdsnws/station/1/query')
    args = ap.parse_args()
    lo, hi, west, east = map(float, args.bounds.split(','))
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    params = {'minlatitude': lo, 'maxlatitude': hi, 'minlongitude': west, 'maxlongitude': east,
              'starttime': args.start, 'endtime': args.end, 'channel': '?H?', 'level': 'channel', 'format': 'xml'}
    url = args.endpoint + '?' + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read(20 * 1024 * 1024 + 1)
    if len(data) > 20 * 1024 * 1024:
        raise ValueError('metadata exceeds 20 MB request budget')
    p = out / 'stations.xml'
    p.write_bytes(data)
    from obspy import read_inventory
    table = triplet_epochs(read_inventory(str(p)), pd.Timestamp(args.start, tz='UTC').timestamp(),
                           pd.Timestamp(args.end, tz='UTC').timestamp())
    table.to_csv(out / 'triplet_epochs.csv', index=False)
    write_json(out / 'manifest.json', {'url': url, 'parameters': params, 'sha256': sha256(p),
        'source_code_sha256': sha256(__file__), 'n_triplet_epochs': len(table),
        'n_stations': int(table[['network', 'station']].drop_duplicates().shape[0]) if len(table) else 0,
        'interpretation': 'Metadata only; no waveform availability/noise check; instrument alternatives are not additive.'})
    print(f'{len(table)} triplet epochs; {len(data)} metadata bytes; no waveforms downloaded', flush=True)


if __name__ == '__main__':
    main()
