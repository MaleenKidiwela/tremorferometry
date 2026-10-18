"""Bounded recent waveform probes; not a measurement of monthly coverage."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import io
import urllib.error
import urllib.parse
import urllib.request
import pandas as pd
from obspy import read
from .provenance import sha256, write_json


def main():
    out = Path('results/recent_waveform_probes_v1')
    out.mkdir(exist_ok=True)
    epochs = pd.read_csv('results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv')
    start = pd.Timestamp('2026-09-13T12:00:00Z')
    end = start + pd.Timedelta(seconds=30)
    epochs = epochs[(pd.to_datetime(epochs.start_utc, utc=True) <= start) &
                    (pd.to_datetime(epochs.end_utc_exclusive, utc=True) >= end)]
    def probe(row):
        station = f'{row.network}.{row.station}'
        receipt = out / f'{station}.json'
        if receipt.exists():
            return __import__('json').loads(receipt.read_text())
        channels = row.channels.split(',')
        location = channels[0].split('.')[2] or '--'
        endpoint = ('https://service.ncedc.org' if row.network in ('BK', 'NC') else
                    'https://service.earthscope.org') + '/fdsnws/dataselect/1/query'
        params = dict(net=row.network, sta=row.station, loc=location,
                      cha=','.join(c.split('.')[-1] for c in channels),
                      starttime=start.strftime('%Y-%m-%dT%H:%M:%S'),
                      endtime=end.strftime('%Y-%m-%dT%H:%M:%S'))
        url = endpoint + '?' + urllib.parse.urlencode(params)
        record = dict(station=station, url=url, start_utc=start.isoformat(),
                      end_utc_exclusive=end.isoformat(), requested_channels=channels)
        try:
            with urllib.request.urlopen(url, timeout=20) as response:
                data = response.read(2_000_001)
                record['http_status'] = response.status
            if not data:
                record['status'] = 'no_data_at_probe_time'
            else:
                if len(data) > 2_000_000:
                    raise ValueError('Probe exceeded 2 MB limit')
                stream = read(io.BytesIO(data))
                path = out / f'{station}.mseed'
                path.write_bytes(data)
                record.update(status='waveform_received', bytes=len(data), sha256=sha256(path),
                              received_channels=sorted(set(tr.id for tr in stream)),
                              all_requested_channels_present=all(stream.select(id=c) for c in channels))
        except Exception as exc:
            record.update(status='request_failed_not_evidence_of_absence', error=f'{type(exc).__name__}: {exc}')
        write_json(receipt, record)
        print(station, record['status'], flush=True)
        return record
    with ThreadPoolExecutor(max_workers=2) as pool:
        records = list(pool.map(probe, epochs.itertuples()))
    write_json(out/'manifest.json', dict(records=records, source_sha256=sha256(__file__),
        interpretation='Thirty-second probes only. Success does not establish monthly coverage; failure is not evidence of archival absence. No detector or catalog validation implied.'))

if __name__ == '__main__':
    main()
