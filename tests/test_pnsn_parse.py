"""Schema test for PNSN time parsing and bbox filter (no network)."""

from __future__ import annotations

from datetime import datetime

import pandas as pd

from tremorferometry.pnsn import _parse_pnsn_time, filter_bbox


def test_parse_pnsn_time() -> None:
    assert _parse_pnsn_time("Mon, 01 Jan 2024 01:47:30 GMT") == datetime(2024, 1, 1, 1, 47, 30)
    assert _parse_pnsn_time("Sat, 18 Oct 2025 12:00:00 GMT") == datetime(2025, 10, 18, 12, 0, 0)


def test_filter_bbox() -> None:
    df = pd.DataFrame({
        "time": pd.to_datetime(["2025-01-01", "2025-01-02", "2025-01-03"]),
        "lat": [40.0, 48.0, 50.0],
        "lon": [-124.0, -123.0, -122.0],
    })
    out = filter_bbox(df, (47.5, 49.5, -125.0, -122.5))
    assert len(out) == 1
    assert out["lat"].iloc[0] == 48.0
