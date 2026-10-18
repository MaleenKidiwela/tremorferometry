"""Exact metadata exposure for the requested catalog scope; never data availability."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .site_selection import representative_sites
from .provenance import sha256, write_json

START = pd.Timestamp('2011-01-01', tz='UTC')
FALLBACK = pd.Timestamp('2019-01-01', tz='UTC')
END = pd.Timestamp('2026-09-01', tz='UTC')


def exposure_summary(epochs, lo, hi):
    """Integrate simultaneous station/site counts over exact half-open epochs."""
    active = epochs[(epochs.lo < hi) & (epochs.hi > lo)]
    edges = sorted({lo, hi} | {max(lo, t) for t in active.lo}
                   | {min(hi, t) for t in active.hi})
    records = []
    for a, b in zip(edges[:-1], edges[1:]):
        current = active[(active.lo <= a) & (active.hi >= b)]
        records.append(((b-a).total_seconds(), current.id.nunique(), current.site_id.nunique()))
    values = np.asarray(records, dtype=float)
    return dict(mean_stations=float(np.average(values[:, 1], weights=values[:, 0])),
                min_stations=int(values[:, 1].min()), max_stations=int(values[:, 1].max()),
                mean_sites=float(np.average(values[:, 2], weights=values[:, 0])),
                min_sites=int(values[:, 2].min()), max_sites=int(values[:, 2].max()))


def main():
    out = Path('results/network_temporal_scope_v1')
    out.mkdir(exist_ok=True)
    source = Path('results/full_history_work_plan_v1/nonoverlapping_triplet_epochs.csv')
    roster_path = Path('results/target_network_from_map/stations.csv')
    regions_path = Path('results/catalog_regional_preflight_v1/station_cached_waveform_inventory.csv')
    roster = pd.read_csv(roster_path)
    roster['id'] = roster.network + '.' + roster.station
    regions = pd.read_csv(regions_path)
    regions['id'] = regions.network + '.' + regions.station
    roster = roster.merge(regions[['id', 'region']], on='id', validate='one_to_one')
    assert len(roster) == 183
    representatives, groups = representative_sites(roster)
    roster['site_id'] = roster.id.map(groups)
    # A site's region stays fixed even if a different colocated sensor is active.
    site_regions = representatives.set_index('id').region
    roster['region'] = roster.site_id.map(site_regions)
    roster.to_csv(out/'fixed_site_groups.csv', index=False)
    epochs = pd.read_csv(source)
    epochs['id'] = epochs.network + '.' + epochs.station
    assert set(epochs.id) <= set(roster.id)
    epochs = epochs.merge(roster[['id', 'site_id', 'region']], on='id', validate='many_to_one')
    epochs['lo'] = pd.to_datetime(epochs.start_utc, utc=True).clip(lower=START)
    epochs['hi'] = pd.to_datetime(epochs.end_utc_exclusive, utc=True).clip(upper=END)
    epochs = epochs[epochs.lo < epochs.hi].copy()
    epochs['bytes_per_second'] = epochs.uncompressed_float32_bytes / epochs.metadata_seconds
    epochs['start_utc'] = epochs.lo.map(lambda x: x.isoformat())
    epochs['end_utc_exclusive'] = epochs.hi.map(lambda x: x.isoformat())
    epochs['metadata_seconds'] = (epochs.hi-epochs.lo).dt.total_seconds()
    epochs['uncompressed_float32_bytes'] = (epochs.metadata_seconds * epochs.bytes_per_second).astype('int64')
    epochs.drop(columns=['lo', 'hi', 'bytes_per_second']).to_csv(out/'scoped_triplet_epochs.csv', index=False)
    edges = pd.date_range(START, END, freq='MS')
    rows, jobs, queue = [], [], []
    region_order = ['north', 'washington', 'oregon', 'southern_edge']
    for lo, hi in zip(edges[:-1], edges[1:]):
        for region in ['all'] + region_order:
            subset = epochs if region == 'all' else epochs[epochs.region == region]
            rows.append(dict(month=lo.strftime('%Y-%m'), region=region,
                             **exposure_summary(subset, lo, hi)))
        active = epochs[(epochs.lo < hi) & (epochs.hi > lo)]
        month_bytes = 0
        for e in active.itertuples():
            a, b = max(lo, e.lo), min(hi, e.hi)
            nbytes = int((b-a).total_seconds()*e.bytes_per_second)
            month_bytes += nbytes
            jobs.append(dict(month=lo.strftime('%Y-%m'), station_id=e.id, site_id=e.site_id,
                             channels=e.channels, start_utc=a.isoformat(), end_utc_exclusive=b.isoformat(),
                             planned_uncompressed_float32_bytes=nbytes,
                             status='awaiting_regional_validation_and_waveform_availability'))
        queue.append(dict(month=lo.strftime('%Y-%m'), start_utc=lo.isoformat(),
                          end_utc_exclusive=hi.isoformat(), station_epoch_blocks=len(active),
                          planned_uncompressed_float32_bytes=month_bytes,
                          in_fallback_range=bool(lo >= FALLBACK),
                          status='awaiting_regional_validation_and_waveform_availability'))
    monthly = pd.DataFrame(rows)
    monthly.to_csv(out/'monthly_network_exposure.csv', index=False)
    pd.DataFrame(jobs).sort_values(['month', 'station_id'], ascending=[False, True]).to_csv(out/'station_month_work_plan.csv', index=False)
    pd.DataFrame(queue).sort_values('month', ascending=False).to_csv(out/'descending_month_queue.csv', index=False)
    comparison = monthly[monthly.month.isin(['2011-01', '2019-01', '2026-08'])]
    comparison.to_csv(out/'endpoint_comparison.csv', index=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True, layout='constrained')
    dates = pd.to_datetime(monthly[monthly.region == 'all'].month)
    whole = monthly[monthly.region == 'all']
    axes[0].plot(dates, whole.mean_stations, label='Station codes', color='0.5')
    axes[0].plot(dates, whole.mean_sites, label='Physical site groups (1 km)', color='black')
    axes[0].set_ylabel('Whole-roster count')
    axes[0].legend()
    for region in region_order:
        g = monthly[monthly.region == region]
        line, = axes[1].plot(dates, g.mean_sites, label=region.replace('_', ' ').title())
        axes[1].fill_between(dates, g.min_sites, g.max_sites, color=line.get_color(), alpha=.12)
    axes[1].set_ylabel('Regional physical sites')
    axes[1].legend(ncol=2)
    for ax in axes:
        ax.axvline(pd.Timestamp('2019-01-01'), color='gray', linestyle='--', alpha=.6)
        ax.grid(alpha=.2)
    fig.suptitle('183-station roster: January 2011–August 2026\nMonthly mean metadata coverage; shading = within-month minimum/maximum')
    fig.supxlabel('Metadata epochs indicate potential recording, not waveform uptime or LFE sensitivity')
    fig.savefig(out/'network_history.png', dpi=160)
    plt.close(fig)
    report = dict(status='metadata_scope_assessment_complete', preferred_start=START.isoformat(),
                  fallback_start=FALLBACK.isoformat(), end_utc_exclusive=END.isoformat(),
                  preferred_months=len(queue), fallback_months=sum(q['in_fallback_range'] for q in queue),
                  roster_stations=len(roster), fixed_physical_site_groups=roster.site_id.nunique(),
                  station_epoch_blocks=len(jobs), raw_bytes_downloaded=0, production_validated=False,
                  endpoint_comparison=comparison.to_dict('records'),
                  method='Exact duration-weighted simultaneous coverage of selected three-component instrument epochs; fixed connected 1 km site groups from the exact roster.',
                  limitations=['Metadata is not actual waveform availability, response suitability, or LFE sensitivity.',
                               'Regional bins describe station geography, not source-region detection coverage.',
                               'Physical proximity grouping does not prove independent noise between groups.',
                               'BH/HH/EH/SH/VH planning preference is inherited, not a validated detector instrument choice.',
                               'Station movements and instrument changes require epoch-specific waveform and response checks.',
                               'No additional reviewed P/S labels exist, per the user; previous input requests remain withdrawn.'],
                  input_sha256={str(p): sha256(p) for p in [source, roster_path, regions_path, Path(__file__)]})
    write_json(out/'report.json', report)
    print(comparison.to_string(index=False))
    print(json.dumps({k: report[k] for k in ['preferred_months', 'fallback_months', 'station_epoch_blocks']}))


if __name__ == '__main__':
    main()
