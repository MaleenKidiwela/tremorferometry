"""LFE template matching.

Default backend is EQcorrscan's `Tribe.detect` which can dispatch to
`fast-matched-filter` (GPU) when installed — controlled by EQcorrscan's
`correlate` config. We don't reimplement matched-filter here; we just glue
the Bostock templates + continuous data into EQcorrscan's API and write
detections to parquet.

For the L40S GPU on this host, install the optional `gpu` extra:

    pip install -e .[gpu]

EQcorrscan will then pick up `fast-matched-filter` automatically when
`xcorr_func='fmf'`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)


def detect_family(
    template_path: Path,
    waveform_paths: list[Path],
    cc_threshold: float,
    mad_multiplier: float,
    use_gpu: bool = True,
) -> pd.DataFrame:
    """Run matched-filter detection for one family across the supplied data.

    Returns a DataFrame with columns matching `io.DETECTIONS_COLUMNS`.

    This is intentionally a thin wrapper:  templates are read from MSEED via
    obspy, continuous data is opened with obspy.read, EQcorrscan's `Tribe` and
    `Tribe.detect` do the work. The GPU path is selected by setting
    `xcorr_func='fmf'` via `eqcorrscan.utils.correlate.set_xcorr`.
    """
    # Lazy imports keep the package importable without the heavy stack.
    from eqcorrscan.core.match_filter import Tribe
    from obspy import read

    if use_gpu:
        try:
            from eqcorrscan.utils.correlate import set_xcorr

            set_xcorr("fmf")
            log.info("matched filter using GPU (fast-matched-filter)")
        except Exception as e:  # noqa: BLE001
            log.warning("FMF unavailable, falling back to CPU: %s", e)

    template_stream = read(str(template_path))
    family_id = template_path.stem
    tribe = Tribe().construct(
        method="from_meta_file",
        meta_file=None,
        st=template_stream,
        process=False,
        lowcut=2.0,
        highcut=8.0,
        samp_rate=template_stream[0].stats.sampling_rate,
        length=template_stream[0].stats.endtime - template_stream[0].stats.starttime,
        prepick=0.5,
    )

    rows: list[dict] = []
    for wf_path in waveform_paths:
        try:
            st = read(str(wf_path))
        except Exception as e:  # noqa: BLE001
            log.warning("read failed %s: %s", wf_path, e)
            continue
        party = tribe.detect(
            stream=st,
            threshold=cc_threshold,
            threshold_type="MAD",
            trig_int=6.0,
            plot=False,
            parallel_process=True,
            xcorr_func="fmf" if use_gpu else None,
        )
        for fam in party:
            for det in fam.detections:
                rows.append(
                    {
                        "family_id": family_id,
                        "station": det.template_name,
                        "channel": "",
                        "time": det.detect_time.datetime,
                        "cc": float(det.detect_val),
                        "shift": 0.0,
                    }
                )
    return pd.DataFrame(rows)
