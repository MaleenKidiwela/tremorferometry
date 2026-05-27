"""Download tremor events from the PNSN tremor API.

Endpoint:
    GET https://tremorapi.pnsn.org/api/v3.0/events?starttime=YYYY-MM-DD&endtime=YYYY-MM-DD

Response is GeoJSON-style:
    {"count": N, "features": [{"geometry": {"coordinates": [lon, lat], ...},
                                "properties": {"time": "RFC1123", "depth": ...,
                                                "duration": ..., "magnitude": ...,
                                                "num_stas": ..., "id": ...}}, ...]}

We chunk requests in time windows (default 14 days) to keep responses moderate,
parse into a DataFrame with `time` (UTC), `lat`, `lon`, `depth`, plus the rest
of the PNSN fields preserved.
"""

from __future__ import annotations

import logging
from datetime import datetime, timedelta
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests

log = logging.getLogger(__name__)

API_URL = "https://tremorapi.pnsn.org/api/v3.0/events"
DEFAULT_CHUNK_DAYS = 14
DEFAULT_TIMEOUT = 60


def _parse_pnsn_time(s: str) -> datetime:
    """PNSN returns RFC1123 (e.g. 'Mon, 01 Jan 2024 01:47:30 GMT'); convert to naive UTC."""
    dt = parsedate_to_datetime(s)
    if dt.tzinfo is not None:
        dt = dt.astimezone(tz=None).replace(tzinfo=None)
    return dt


def _date_chunks(t0: datetime, t1: datetime, days: int) -> Iterable[tuple[datetime, datetime]]:
    cur = t0
    while cur < t1:
        nxt = min(cur + timedelta(days=days), t1)
        yield cur, nxt
        cur = nxt


def fetch_events(
    t_start: datetime,
    t_end: datetime,
    chunk_days: int = DEFAULT_CHUNK_DAYS,
    timeout: int = DEFAULT_TIMEOUT,
    session: requests.Session | None = None,
) -> pd.DataFrame:
    """Fetch all tremor events in [t_start, t_end), chunked. Returns a DataFrame."""
    sess = session or requests.Session()
    rows: list[dict] = []
    for c0, c1 in _date_chunks(t_start, t_end, chunk_days):
        params = {"starttime": c0.date().isoformat(), "endtime": c1.date().isoformat()}
        log.info("PNSN tremor fetch %s..%s", params["starttime"], params["endtime"])
        r = sess.get(API_URL, params=params, timeout=timeout)
        if r.status_code == 404:
            # API returns 404 for date ranges with no events; treat as empty.
            log.info("  no events in this chunk (404)")
            continue
        r.raise_for_status()
        payload = r.json()
        for feat in payload.get("features", []):
            geom = feat.get("geometry", {})
            coords = geom.get("coordinates", [None, None])
            props = feat.get("properties", {})
            t_str = props.get("time")
            if not t_str:
                continue
            rows.append(
                {
                    "time": _parse_pnsn_time(t_str),
                    "lat": float(coords[1]) if coords[1] is not None else None,
                    "lon": float(coords[0]) if coords[0] is not None else None,
                    "depth": float(props["depth"]) if props.get("depth") is not None else None,
                    "magnitude": props.get("magnitude"),
                    "num_stas": props.get("num_stas"),
                    "duration": props.get("duration"),
                    "energy": props.get("energy"),
                    "id": props.get("id"),
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df = df.sort_values("time").reset_index(drop=True)
    return df


def filter_bbox(
    df: pd.DataFrame, bbox: tuple[float, float, float, float]
) -> pd.DataFrame:
    """Restrict to lat_min, lat_max, lon_min, lon_max."""
    lat_min, lat_max, lon_min, lon_max = bbox
    mask = (
        (df["lat"] >= lat_min)
        & (df["lat"] <= lat_max)
        & (df["lon"] >= lon_min)
        & (df["lon"] <= lon_max)
    )
    return df.loc[mask].reset_index(drop=True)


def write_csv(df: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
