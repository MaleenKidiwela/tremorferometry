"""FDSN waveform fetching with on-disk day-file cache.

Layout:  data/waveforms/{net}.{sta}/{year}/{jday}.mseed
One MSEED per station per UTC day; all channels concatenated.

Parallelism: ThreadPoolExecutor across (station, day) — FDSN is I/O bound.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

log = logging.getLogger(__name__)

DEFAULT_PROVIDERS = ("IRIS", "NCEDC")


def _path_for(root: Path, net: str, sta: str, day: datetime) -> Path:
    return root / f"{net}.{sta}" / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"


def daterange(t0: datetime, t1: datetime) -> Iterable[datetime]:
    cur = datetime(t0.year, t0.month, t0.day)
    while cur < t1:
        yield cur
        cur += timedelta(days=1)


def fetch_day(
    client: Client,
    root: Path,
    net: str,
    sta: str,
    day: datetime,
    location: str = "*",
    channel: str = "HH?,BH?,EH?",
    overwrite: bool = False,
) -> Path | None:
    out = _path_for(root, net, sta, day)
    if out.exists() and not overwrite:
        return out
    t0 = UTCDateTime(day)
    t1 = t0 + 86400
    try:
        st = client.get_waveforms(net, sta, location, channel, t0, t1)
    except Exception as e:  # noqa: BLE001
        log.warning("fetch failed %s.%s %s: %s", net, sta, day.date(), e)
        return None
    if len(st) == 0:
        return None
    out.parent.mkdir(parents=True, exist_ok=True)
    st.write(str(out), format="MSEED")
    return out


def fetch_many(
    root: Path,
    network: str,
    stations: Iterable[str],
    t_start: datetime,
    t_end: datetime,
    provider: str = "IRIS",
    max_workers: int = 16,
    channel: str = "HH?,BH?,EH?",
) -> list[Path]:
    """Parallel FDSN fetch across stations x days. Returns list of written paths."""
    client = Client(provider)
    tasks = [
        (net, sta, day)
        for net in network.split(",")
        for sta in stations
        for day in daterange(t_start, t_end)
    ]
    paths: list[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {
            ex.submit(fetch_day, client, root, net, sta, day, channel=channel): (net, sta, day)
            for net, sta, day in tasks
        }
        for f in as_completed(futures):
            p = f.result()
            if p is not None:
                paths.append(p)
    return paths
