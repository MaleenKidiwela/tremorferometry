"""Import original source records without inheriting previous pipeline corrections."""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd


def read_bostock(path):
    d = pd.read_csv(path)
    if d.duplicated().any():
        raise ValueError("duplicate source rows; resolve provenance before importing")
    day = pd.to_datetime(dict(year=d.year, month=d.month, day=d.day), utc=True)
    time = day + pd.to_timedelta(d.hour * 3600 + d['min'] * 60 + d.sec, unit="s")
    # Explicit ns conversion: source precision is never inferred from int64 magnitude.
    epoch = time.astype("datetime64[ns, UTC]").astype("int64").to_numpy() / 1e9
    keys = [f"{la:.5f},{lo:.5f},{de:.5f}" for la, lo, de in zip(d.lat, d.lon, d.dep)]
    # The published table has no original family ID. Exact coordinate groups are
    # reproducible source groups, not a claim to recover the original family map.
    family = ["bostock_xyz_" + hashlib.sha256(k.encode()).hexdigest()[:12] for k in keys]
    return pd.DataFrame(dict(reference_id=[f"bostock_row_{i}" for i in range(len(d))],
                             time_s=epoch, time_basis="bostock_published_time_unresolved",
                             family_id=family, latitude=d.lat, longitude=d.lon, depth_km=d.dep,
                             published_magnitude=d.mag, date=time.dt.strftime("%Y-%m-%d")))


def choose_families(reference, dates, number, separation_km=10):
    """Training-only source counts and coordinates. Never inspect test performance."""
    d = reference[reference.date.isin(dates)]
    counts = d.groupby("family_id").size().sort_values(ascending=False)
    chosen, points = [], []
    for family in counts.index:
        row = d[d.family_id == family].iloc[0]
        point = np.array([row.latitude, row.longitude])
        if any(np.hypot((point[0] - p[0]) * 111.195,
                        (point[1] - p[1]) * 111.195 * np.cos(np.deg2rad(point[0]))) < separation_km for p in points):
            continue
        chosen.append(family)
        points.append(point)
        if len(chosen) == number:
            break
    return chosen
