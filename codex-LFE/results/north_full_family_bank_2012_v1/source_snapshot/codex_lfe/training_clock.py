"""Development reference-clock calibration using training-date S picks only.

This estimates a model-dependent event clock. It neither creates manually
reviewed phase labels nor validates absolute origin times.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import yaml

from .provenance import sha256, write_json
from .reference import choose_families, read_bostock


def peak_offset(reference_times, pick_times, half_width=30., bin_s=.25):
    """Peak in all-pair residual histogram; no nearest-pick zero-lag bias."""
    picks = np.sort(np.asarray(pick_times))
    parts = []
    for t in reference_times:
        lo, hi = np.searchsorted(picks, [t-half_width, t+half_width])
        parts.append(picks[lo:hi]-t)
    delta = np.concatenate(parts) if parts else np.array([])
    edges = np.arange(-half_width, half_width+bin_s/2, bin_s)
    hist, _ = np.histogram(delta, edges)
    center = (edges[:-1]+edges[1:])/2
    k = int(np.argmax(hist))
    offset = float(center[k])
    background = float(np.mean(hist[np.abs(center-offset)>3]))
    contrast = float(hist[k]/max(background, 1.))
    return {'offset_s': offset, 'peak_count': int(hist[k]), 'background_per_bin': background,
            'peak_over_background': contrast, 'n_reference': len(reference_times),
            'passed_peak_check': bool(hist[k]>=5 and contrast>=3 and abs(offset)<half_width-1)}


def accepted_site_clocks(rows, minimum_sites=2, max_spread_s=1.):
    """Nearby instruments share one vote; discordant sites invalidate the clock."""
    passed = rows[rows.accepted_station]
    if passed.empty:
        return None, 0, None
    by_site = passed.groupby('site_id').offset_s.median()
    spread = float(by_site.max()-by_site.min())
    offset = float(by_site.median()) if len(by_site)>=minimum_sites and spread<=max_spread_s else None
    return offset, len(by_site), spread


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True)
    ap.add_argument('--picks', default='../lfe_catalog_repair/results/e2_picks_vi2010_cascadia.parquet')
    ap.add_argument('--travel-table', default='../lfe_catalog_repair/results/e2_ttable_crosson.npz')
    ap.add_argument('--metadata', default='../lfe_catalog_repair/results/e2_stations.csv')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path(args.config).read_text())
    dates = sorted(cfg['train_dates'])
    early, late = dates[:len(dates)//2], dates[len(dates)//2:]
    if not early or not late or set(dates)&set(cfg['calibration_dates']+cfg['test_dates']):
        raise ValueError('need separated training halves and disjoint evaluation dates')
    ref = read_bostock(cfg['reference_csv'])
    families = choose_families(ref, dates, cfg['family_count'], cfg['family_separation_km'])
    ref = ref[ref.date.isin(dates)&ref.family_id.isin(families)].copy()
    picks = pd.read_parquet(args.picks)
    picks['station'] = picks.net+'.'+picks.sta
    picks['time_s'] = pd.to_datetime(picks.time, utc=True).dt.as_unit('ns').astype('int64')/1e9
    picks['date'] = pd.to_datetime(picks.time, utc=True).dt.strftime('%Y-%m-%d')
    picks = picks[picks.date.isin(dates)&picks.station.isin(cfg['stations'])
                  &(picks.phase=='S')&(picks.prob>=.3)&(picks.thr_first<=.3)].copy()
    if picks.model.nunique()!=1:
        raise ValueError('one model required')
    sod = picks.time_s%86400
    picks = picks[(sod>=1)&(sod<86399)]
    metadata = pd.read_csv(args.metadata)
    metadata['station'] = metadata.net+'.'+metadata.sta
    metadata = metadata.set_index('station')
    with np.load(args.travel_table) as table:
        travel = RegularGridInterpolator((table['x'],table['z']),table['S'],bounds_error=True)
    rows = []
    for family, frame in ref.groupby('family_id'):
        source = frame.iloc[0]
        for station, info in cfg['stations'].items():
            m = metadata.loc[station]
            epi = np.hypot((source.latitude-m.lat)*111.195,
                           (source.longitude-m.lon)*111.195*np.cos(np.deg2rad((source.latitude+m.lat)/2)))
            tt = float(travel([[epi,source.depth_km]])[0])
            records = {}
            for label, half in [('early',early),('late',late),('all',dates)]:
                g = frame[frame.date.isin(half)]
                p = picks[(picks.station==station)&picks.date.isin(half)]
                records[label] = peak_offset(g.time_s.to_numpy()+tt,p.time_s.to_numpy())
            a,b = records['early'],records['late']
            drift = abs(a['offset_s']-b['offset_s'])
            accepted = a['passed_peak_check'] and b['passed_peak_check'] and drift<=1.
            rows.append({'family_id':family,'station':station,'site_id':info['site_id'],
                         'travel_s':tt,'offset_s':records['all']['offset_s'],
                         'accepted_station':bool(accepted),'early_late_difference_s':drift,
                         **{f'{label}_{key}':value for label,r in records.items() for key,value in r.items()}})
    station_rows = pd.DataFrame(rows)
    station_rows.to_csv(out/'station_clocks.csv',index=False)
    clocks = []
    for family, frame in station_rows.groupby('family_id'):
        offset,n,spread = accepted_site_clocks(frame)
        clocks.append({'family_id':family,'origin_minus_published_s':offset,
                       'n_passing_site_groups':n,'site_offset_range_s':spread,'accepted':offset is not None})
    write_json(out/'family_clocks.json',clocks)
    sources = [args.config,args.picks,args.travel_table,args.metadata,cfg['reference_csv'],__file__]
    write_json(out/'manifest.json',{'sources':{str(Path(p).resolve()):sha256(p) for p in sources},
        'train_dates':dates,'early_dates':early,'late_dates':late,'families_selected':families,
        'reference_ids_used':ref.reference_id.tolist(),'n_training_s_picks':len(picks),
        'criteria':{'min_peak_bin_count':5,'min_peak_over_background':3.,'max_early_late_drift_s':1.,
                    'minimum_independent_site_groups':2,'max_site_clock_range_s':1.},
        'origin_time_basis':'crosson_model_training_pick_calibrated_utc',
        'absolute_origin_time_validated':False,
        'limitations':['S picks not independently hand labeled','Travel-time model systematics',
                       'Pretrained exposure not excluded','QC thresholds are development choices',
                       'No complete inference coverage mask for legacy picks']})
    print(pd.DataFrame(clocks).to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
