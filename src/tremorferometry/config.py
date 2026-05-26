"""YAML-backed config for one ETS episode."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Episode:
    t_start: datetime
    t_end: datetime
    bbox: tuple[float, float, float, float]  # lat_min, lat_max, lon_min, lon_max


@dataclass(frozen=True)
class Stations:
    network: str
    list: tuple[str, ...]


@dataclass(frozen=True)
class Match:
    cc_threshold: float = 0.4
    mad_multiplier: float = 8.0


@dataclass(frozen=True)
class Stack:
    bin_days: int = 2


@dataclass(frozen=True)
class Dvv:
    coda_window: tuple[float, float] = (5.0, 25.0)
    freq_band: tuple[float, float] = (2.0, 8.0)
    stretch_range: float = 0.02
    stretch_steps: int = 201
    min_detections: int = 20
    min_cc: float = 0.6
    reference_window: tuple[int, int] = (-90, -30)


@dataclass(frozen=True)
class Config:
    event_id: str
    episode: Episode
    stations: Stations
    match: Match = field(default_factory=Match)
    stack: Stack = field(default_factory=Stack)
    dvv: Dvv = field(default_factory=Dvv)


def _to_datetime(v) -> datetime:
    if isinstance(v, datetime):
        return v
    if isinstance(v, date):
        return datetime(v.year, v.month, v.day)
    if isinstance(v, str):
        return datetime.fromisoformat(v)
    raise TypeError(f"cannot parse datetime from {v!r}")


def load_config(path: str | Path) -> Config:
    """Load a YAML config and validate into the frozen Config dataclass."""
    with open(path) as f:
        raw = yaml.safe_load(f)

    ep = raw["episode"]
    episode = Episode(
        t_start=_to_datetime(ep["t_start"]),
        t_end=_to_datetime(ep["t_end"]),
        bbox=tuple(ep["bbox"]),  # type: ignore[arg-type]
    )
    st = raw["stations"]
    stations = Stations(network=str(st["network"]), list=tuple(st["list"]))

    match = Match(**raw.get("match", {}))
    stack = Stack(**raw.get("stack", {}))

    dvv_raw = raw.get("dvv", {})
    dvv = Dvv(
        coda_window=tuple(dvv_raw.get("coda_window", (5.0, 25.0))),  # type: ignore[arg-type]
        freq_band=tuple(dvv_raw.get("freq_band", (2.0, 8.0))),  # type: ignore[arg-type]
        stretch_range=float(dvv_raw.get("stretch_range", 0.02)),
        stretch_steps=int(dvv_raw.get("stretch_steps", 201)),
        min_detections=int(dvv_raw.get("min_detections", 20)),
        min_cc=float(dvv_raw.get("min_cc", 0.6)),
        reference_window=tuple(dvv_raw.get("reference_window", (-90, -30))),  # type: ignore[arg-type]
    )

    return Config(
        event_id=str(raw["event_id"]),
        episode=episode,
        stations=stations,
        match=match,
        stack=stack,
        dvv=dvv,
    )
