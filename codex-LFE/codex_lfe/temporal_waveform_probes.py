"""Small waveform-availability snapshots, one month at a time, newest first."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import io
import json
import urllib.parse
import urllib.request
import numpy as np
import pandas as pd
from obspy import read
from .provenance import sha256, write_json


def covered_seconds(intervals, start, end):
    cursor = start
    total = 0.
    for lo, hi in sorted(intervals):
        lo, hi = max(lo, start, cursor), min(hi, end)
        if hi > lo:
            total += hi-lo
            cursor = hi
    return total


def main():
    out = Path('results/temporal_waveform_probes_v1')
    out.mkdir(exist_ok=True)
    source = Path('results/network_temporal_scope_v1/scoped_triplet_epochs.csv')
    epochs = pd.read_csv(source)
    epochs['lo'] = pd.to_datetime(epochs.start_utc, utc=True)
    epochs['hi'] = pd.to_datetime(epochs.end_utc_exclusive, utc=True)
    dates = ['2026-08-15T12:00:00Z', '2019-01-15T12:00:00Z', '2011-01-15T12:00:00Z']
    write_json(out/'design.json', dict(dates=dates, duration_s=30, max_workers=2,
        byte_limit_per_request=2000000, processing_order='One snapshot month at a time, newest first',
        purpose='Test advertised historical three-component availability, not detector performance or monthly completeness',
        metadata_sha256=sha256(source), source_sha256=sha256(__file__),
        retained_raw_policy='Small validation snippets retained with receipts and hashes; no continuous archive download'))
    records = []
    for date in dates:
        start = pd.Timestamp(date)
        end = start + pd.Timedelta(seconds=30)
        active = epochs[(epochs.lo <= start) & (epochs.hi >= end)]
        month = start.strftime('%Y-%m')
        stage = out/month
        stage.mkdir(exist_ok=True)
        def probe(row):
            path = stage/f'{row.id}.mseed'
            receipt = stage/f'{row.id}.json'
            if receipt.exists():
                previous = json.loads(receipt.read_text())
                if previous.get('sha256') and (not path.exists() or sha256(path) != previous['sha256']):
                    raise ValueError('Previously retained probe changed')
                return previous
            channels = row.channels.split(',')
            location = channels[0].split('.')[2] or '--'
            endpoint = ('https://service.ncedc.org' if row.network in ('BK', 'NC') else
                        'https://service.earthscope.org') + '/fdsnws/dataselect/1/query'
            params = dict(net=row.network, sta=row.station, loc=location,
                          cha=','.join(c.split('.')[-1] for c in channels),
                          starttime=start.strftime('%Y-%m-%dT%H:%M:%S'),
                          endtime=end.strftime('%Y-%m-%dT%H:%M:%S'))
            url = endpoint+'?'+urllib.parse.urlencode(params)
            record = dict(month=month, station_id=row.id, site_id=row.site_id, region=row.region,
                          start_utc=start.isoformat(), end_utc_exclusive=end.isoformat(),
                          requested_channels=channels, url=url, bytes=0, full_three_component_coverage=False)
            try:
                with urllib.request.urlopen(url, timeout=20) as response:
                    data = response.read(2000001)
                    record['http_status'] = response.status
                if not data:
                    record['status'] = 'no_data_at_probe_time'
                else:
                    if len(data) > 2000000:
                        raise ValueError('Probe exceeds 2 MB limit')
                    stream = read(io.BytesIO(data))
                    path.write_bytes(data)
                    coverage = {}
                    for channel in channels:
                        spans = [(float(tr.stats.starttime), float(tr.stats.endtime)+tr.stats.delta)
                                 for tr in stream.select(id=channel)
                                 if np.isfinite(tr.data).all()]
                        coverage[channel] = covered_seconds(spans, start.timestamp(), end.timestamp())
                    record.update(status='waveform_received', bytes=len(data),
                                  sha256=hashlib.sha256(data).hexdigest(), channel_coverage_s=coverage,
                                  received_channels=sorted({tr.id for tr in stream}),
                                  full_three_component_coverage=all(s >= 30-1e-5 for s in coverage.values()))
            except Exception as exc:
                record.update(status='request_failed_not_evidence_of_absence', error=f'{type(exc).__name__}: {exc}')
            write_json(receipt, record)
            print(month, row.id, record['status'], record['full_three_component_coverage'], flush=True)
            return record
        with ThreadPoolExecutor(max_workers=2) as pool:
            for record in pool.map(probe, active.itertuples()):
                records.append(record)
        write_json(stage/'complete.json', dict(month=month, stations=len(active),
            complete_three_component_stations=sum(r['full_three_component_coverage'] for r in records if r['month'] == month)))
    frame = pd.DataFrame(records)
    frame.to_json(out/'receipts.json', orient='records', indent=2)
    summary = []
    for (month, region), group in frame.groupby(['month', 'region']):
        good = group[group.full_three_component_coverage]
        summary.append(dict(month=month, region=region, metadata_stations=len(group),
                            full_coverage_stations=len(good), full_coverage_sites=good.site_id.nunique(),
                            metadata_sites=group.site_id.nunique(), downloaded_bytes=int(group.bytes.sum()),
                            failed_requests=int((group.status=='request_failed_not_evidence_of_absence').sum())))
    pd.DataFrame(summary).to_csv(out/'summary.csv', index=False)
    write_json(out/'report.json', dict(status='complete_bounded_temporal_availability_probe',
        summary=summary, raw_bytes_downloaded=int(frame.bytes.sum()), production_validated=False,
        limitations=['One 30-second window per station and epoch, not monthly uptime.',
                     'No-data or failed requests do not establish archive-wide absence.',
                     'Coverage does not establish sensor quality, response correctness, or LFE detectability.']))


if __name__ == '__main__':
    main()
