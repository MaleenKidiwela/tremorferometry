"""Review retained probes at native sample resolution without downloading again."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from obspy import read
from .provenance import sha256, write_json


def sample_grid_complete(intervals, start, end):
    """Allow a sub-sample leading edge, but no internal data gaps.

    Intervals use last-sample time plus delta as the exclusive upper support.
    The final sample must reach the requested end in that convention.
    """
    spans = sorted((lo, hi, dt) for lo, hi, dt in intervals if hi > start and lo < end)
    if not spans or spans[0][0]-start > spans[0][2]+1e-5:
        return False
    cursor = spans[0][1]
    for lo, hi, dt in spans[1:]:
        if lo > cursor+1e-5:
            return False
        cursor = max(cursor, hi)
    return bool(cursor >= end-1e-5)


def main():
    source = Path('results/temporal_waveform_probes_v1')
    if not (source/'report.json').exists():
        raise RuntimeError('Wait for the completed acquisition report before reviewing')
    out = Path('results/temporal_waveform_probe_review_v1')
    out.mkdir(exist_ok=True)
    records = json.loads((source/'receipts.json').read_text())
    details = []
    for record in records:
        complete = False
        channels = {}
        if record.get('sha256'):
            path = source/record['month']/(record['station_id']+'.mseed')
            if sha256(path) != record['sha256']:
                raise ValueError('Retained waveform checksum changed')
            stream = read(str(path))
            start = pd.Timestamp(record['start_utc']).timestamp()
            end = pd.Timestamp(record['end_utc_exclusive']).timestamp()
            for channel in record['requested_channels']:
                spans = [(float(tr.stats.starttime), float(tr.stats.endtime)+tr.stats.delta, tr.stats.delta)
                         for tr in stream.select(id=channel) if np.isfinite(tr.data).all()]
                channels[channel] = sample_grid_complete(spans, start, end)
            complete = all(channels.values())
        details.append(dict(month=record['month'], station_id=record['station_id'],
                            site_id=record['site_id'], region=record['region'],
                            request_status=record['status'], full_native_grid_coverage=complete,
                            strict_continuous_coverage=record['full_three_component_coverage'],
                            channel_native_grid_coverage=json.dumps(channels)))
    frame = pd.DataFrame(details)
    frame.to_csv(out/'station_results.csv', index=False)
    rows = []
    for month, monthly in frame.groupby('month'):
        for region in ['all', 'north', 'washington', 'oregon', 'southern_edge']:
            group = monthly if region == 'all' else monthly[monthly.region == region]
            good = group[group.full_native_grid_coverage]
            rows.append(dict(month=month, region=region, metadata_stations=len(group),
                             metadata_sites=group.site_id.nunique(), received_stations=int((group.request_status=='waveform_received').sum()),
                             complete_native_grid_stations=len(good), complete_native_grid_sites=good.site_id.nunique()))
    summary = pd.DataFrame(rows)
    summary.to_csv(out/'summary.csv', index=False)
    write_json(out/'report.json', dict(status='complete_native_sample_grid_review', summary=rows,
        source_report_sha256=sha256(source/'report.json'), source_sha256=sha256(__file__), raw_bytes_downloaded=0,
        criterion='All three requested channels finite, first sample within one native sample interval of requested start, continuous internal sample support, final sample plus delta at or beyond exclusive end.',
        limitations=['Thirty-second availability only, not monthly uptime or LFE sensitivity.',
                     'Original strict continuous-duration measurements are retained unchanged.',
                     'Data existence and finite samples do not establish nonzero signal, response quality or correct timing.']))
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
