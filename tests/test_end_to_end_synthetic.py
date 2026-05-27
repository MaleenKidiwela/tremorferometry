"""End-to-end integration: synthetic stacks -> measure -> assert recovery."""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

from tremorferometry.measure import measure_many
from tremorferometry.synthetic import (
    dvv_ets_pattern,
    make_bin_edges,
    master_template,
    write_synthetic_family,
)


def test_recovers_imposed_dvv_pattern(tmp_path: Path) -> None:
    fs = 100.0
    ets_start = datetime(2024, 1, 15)
    ets_end = ets_start + timedelta(days=10)
    t_start = ets_start - timedelta(days=60)
    t_end = ets_end + timedelta(days=30)

    edges = make_bin_edges(t_start, t_end, bin_days=2)
    template = master_template(fs=fs)

    family_ids = ["FAMa", "FAMb"]
    stations = ["S00", "S01"]
    for k, fid in enumerate(family_ids):
        write_synthetic_family(
            out_path=tmp_path / f"{fid}.h5",
            family_id=fid,
            stations=stations,
            bin_edges=edges,
            template=template,
            fs=fs,
            ets_start=ets_start,
            ets_end=ets_end,
            dvv_min=-0.005,
            noise_level=0.02,
            n_det_per_bin=50,
            seed=k,
        )

    df = measure_many(
        h5_paths=sorted(tmp_path.glob("*.h5")),
        coda_window=(5.0, 25.0),
        fs=fs,
        eps_max=0.02,
        n_eps=401,
        ref_start=ets_start + timedelta(days=-60),
        ref_end=ets_start + timedelta(days=-20),
        min_cc=0.5,
        n_workers=1,
    )
    assert not df.empty, "measurement returned no rows"

    df["t_center"] = pd.to_datetime(df["t_center"])
    df["dvv_true"] = df["t_center"].apply(
        lambda t: dvv_ets_pattern(t.to_pydatetime(), ets_start, ets_end, dvv_min=-0.005)
    )
    residual = df["dvv"] - df["dvv_true"]
    assert residual.std() < 5e-4, f"recovery rms too high: {residual.std():.5e}"
    assert abs(residual.mean()) < 2e-4, f"bias too high: {residual.mean():+.5e}"

    # the ETS-interior bins must be visibly negative (the imposed signal)
    inside_ets = df[
        (df["t_center"] >= pd.Timestamp(ets_start))
        & (df["t_center"] <= pd.Timestamp(ets_end))
    ]
    assert not inside_ets.empty
    assert inside_ets["dvv"].mean() < -2e-3, (
        f"failed to recover negative dv/v inside ETS: mean = {inside_ets['dvv'].mean():.4e}"
    )
