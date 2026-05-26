"""Round-trip test for the YAML config loader."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from tremorferometry.config import Config, load_config

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = REPO_ROOT / "configs" / "ets_2010_vi.yaml"


def test_loads_example_config() -> None:
    cfg = load_config(EXAMPLE)
    assert isinstance(cfg, Config)
    assert cfg.event_id == "ets_2010_vi"
    assert cfg.episode.t_start == datetime(2010, 8, 1)
    assert cfg.episode.t_end == datetime(2010, 9, 1)
    assert cfg.episode.bbox == (48.0, 49.5, -124.5, -122.5)
    assert "SNB" in cfg.stations.list
    assert cfg.dvv.coda_window == (5.0, 25.0)
    assert cfg.dvv.freq_band == (2.0, 8.0)
    assert cfg.dvv.stretch_steps == 201


def test_defaults_applied(tmp_path: Path) -> None:
    p = tmp_path / "min.yaml"
    p.write_text(
        "event_id: minimal\n"
        "episode:\n"
        "  t_start: 2020-01-01\n"
        "  t_end:   2020-02-01\n"
        "  bbox: [40.0, 41.0, -125.0, -124.0]\n"
        "stations:\n"
        "  network: CN\n"
        "  list: [AAA]\n"
    )
    cfg = load_config(p)
    assert cfg.match.cc_threshold == 0.4
    assert cfg.stack.bin_days == 2
    assert cfg.dvv.stretch_range == 0.02
