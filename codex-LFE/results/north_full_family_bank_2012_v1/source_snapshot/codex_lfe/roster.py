"""Explicit project-network membership; station codes alone are not identities."""
from pathlib import Path
import re
import pandas as pd
from .provenance import sha256


def read_roster(path):
    path = Path(path)
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not {'network', 'station'}.issubset(frame):
        raise ValueError('roster CSV requires network,station columns; no inferred network codes')
    identities = []
    for row in frame.itertuples():
        net, sta = row.network.strip(), row.station.strip()
        if not re.fullmatch(r'[A-Z0-9]{1,8}', net) or not re.fullmatch(r'[A-Z0-9]{1,8}', sta):
            raise ValueError('invalid explicit network/station identity')
        identities.append(f'{net}.{sta}')
    if not identities or len(set(identities)) != len(identities):
        raise ValueError('roster must be nonempty and contain unique network/station pairs')
    return set(identities)


def validate_membership(stations, path):
    allowed = read_roster(path)
    outside = set(stations) - allowed
    if outside:
        raise ValueError(f'stations outside target roster: {sorted(outside)}')
    return {'path': str(Path(path).resolve()), 'sha256': sha256(path),
            'n_allowed': len(allowed), 'n_used': len(set(stations)),
            'used': sorted(set(stations)),
            'limitation': 'Membership only; longevity, instrument epochs, and usable waveform coverage require independent verification.'}
