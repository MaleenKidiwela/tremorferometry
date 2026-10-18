"""Post-fit loss diagnostics on existing development data; no detector retuning."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from .matched_filter import align_to_anchor, combine_sites, normalized_cc, robust_scale
from .pilot import date_start, load_waves
from .provenance import sha256, write_json


def phase_window(values, start, fs, expected, radius, pick_times, pick_scores):
    """Inclusive arrival window; reject missing prediction context, never wrap days.

    A threshold crossing and an exported peak are distinct diagnostics. All
    samples in the window and its one-second peak guard must have predictions.
    """
    lo = int(np.ceil((expected - radius - start) * fs - 1e-5))
    hi = int(np.floor((expected + radius - start) * fs + 1e-5)) + 1
    edge = int(round(fs))
    if lo < edge or hi + edge > len(values) or hi <= lo:
        return False, np.nan, np.nan, np.nan
    if not np.isfinite(values[lo-edge:hi+edge]).all():
        return False, np.nan, np.nan, np.nan
    maximum = float(np.max(values[lo:hi]))
    a = np.searchsorted(pick_times, expected-radius-1e-6, side='left')
    b = np.searchsorted(pick_times, expected+radius+1e-6, side='right')
    if a == b:
        return True, maximum, np.nan, np.nan
    k = a + int(np.argmax(pick_scores[a:b]))
    return True, maximum, float(pick_scores[k]), float(pick_times[k]-expected)


def support_counts(frame, threshold):
    """Count unique S stations and physical sites, never P/S or nearby duplicates."""
    s = frame[frame.phase == 'S']
    hits = s[s.max_pick_score >= threshold]
    return {'all_s_covered': bool(len(s) and s.covered.all()),
            'n_s_stations': hits.station.nunique(), 'n_s_sites': hits.site.nunique(),
            'n_s_crossing_sites': s[s.max_probability >= threshold].site.nunique()}


def matcher_window(z, support, good, index, radius, minimum_sites, threshold):
    ix = np.arange(max(0, index-radius), min(len(z), index+radius+1))
    ix = ix[good[ix]]
    if not len(ix):
        return {'status': 'no_valid_sample'}
    k = ix[np.argmax(z[ix])]
    supported = ix[support[ix] >= minimum_sites]
    max_supported = float(z[supported].max()) if len(supported) else np.nan
    if z[k] < threshold:
        status = 'below_network_threshold'
    elif not len(supported) or max_supported < threshold:
        status = 'network_threshold_passes_but_site_support_fails'
    else:
        status = 'joint_threshold_and_site_support_pass'
    return {'status': status, 'max_network_mad': float(z[k]),
            'max_supported_mad': max_supported, 'support_at_network_max': int(support[k]),
            'lag_at_network_max_samples': int(k-index)}


def matcher_date(task):
    cfg, run_path, date, threshold = task
    run = Path(run_path); inputs = {}; rows = []
    refs = pd.read_parquet(run/'calibration_eligible_reference.parquet')
    refs = refs[refs.date == date]
    waves = load_waves(cfg, date, inputs)
    fs = cfg['fs']; n = int(86400*fs)
    block_n = int(cfg['normalization_block_s']*fs)
    guard = int(cfg['block_guard_s']*fs)
    radius = int(round(cfg['reference_tolerance_s']*fs))
    with np.load(run/'templates.npz') as archive:
        for family, sub in refs.groupby('family_id'):
            aligned = {}
            for station in cfg['stations']:
                key = f'{family}__{station}'
                if key+'__template' not in archive:
                    continue
                cc, valid = normalized_cc(waves[station].x, archive[key+'__template'], waves[station].valid)
                aligned[station], _ = align_to_anchor(cc, valid, float(archive[key+'__offset_s']), fs, n)
            sites = combine_sites(aligned, {s:cfg['stations'][s]['site_id'] for s in aligned})
            data = np.asarray(list(sites.values()))
            blocks = ((sub.time_s-date_start(date))*fs//block_n).astype(int)
            for block, group in sub.groupby(blocks):
                first = block*block_n; x = data[:, first:first+block_n]
                good = np.isfinite(x).all(axis=0)
                good[:guard] = False; good[-guard:] = False
                stats = [robust_scale(v[good]) for v in x]
                scale = robust_scale(x[:, good].mean(axis=0))
                if scale is None or any(s is None for s in stats):
                    raise ValueError('eligible reference has no valid normalization')
                z = (x.mean(axis=0)-scale[0])/scale[1]
                support = np.sum([(v-s[0])/s[1] >= cfg['site_support_mad'] for v,s in zip(x,stats)],axis=0)
                for ref in group.itertuples():
                    index = int(round((ref.time_s-date_start(date))*fs))-first
                    row = matcher_window(z,support,good,index,radius,cfg['min_sites'],threshold)
                    row.update(reference_id=ref.reference_id,family_id=family,date=date,
                               reference_anchor_s=ref.time_s)
                    rows.append(row)
    print(f'[matcher diagnosis] {date}: {len(rows)} references', flush=True)
    return rows, inputs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matcher',required=True); ap.add_argument('--neural',required=True)
    ap.add_argument('--shared',required=True); ap.add_argument('--clocks',required=True)
    ap.add_argument('--association',required=True); ap.add_argument('--output',required=True)
    ap.add_argument('--workers',type=int,default=4)
    ap.add_argument('--travel-table',default='../lfe_catalog_repair/results/e2_ttable_crosson.npz')
    ap.add_argument('--station-metadata',default='../lfe_catalog_repair/results/e2_stations.csv')
    args = ap.parse_args(); out = Path(args.output)
    if not 1 <= args.workers <= 8: raise ValueError('workers must be 1..8')
    if out.exists() and any(out.iterdir()): raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True); started=time.monotonic(); sources={}
    def record(path): sources[str(Path(path).resolve())]=sha256(path)
    mf, neural, shared, clocks, assoc = map(Path,[args.matcher,args.neural,args.shared,args.clocks,args.association])
    cfg=json.loads((mf/'manifest.json').read_text())['config']
    accepted={r['family_id']:r['origin_minus_published_s'] for r in
              json.loads((clocks/'family_clocks.json').read_text()) if r['accepted']}
    transfers=pd.read_csv(clocks/'station_transfers.csv')
    transfers=transfers[transfers.accepted_waveform_transfer & transfers.family_id.isin(accepted)]
    offsets=transfers.set_index(['family_id','station','phase']).phase_arrival_after_published_s.to_dict()
    meta=pd.read_csv(args.station_metadata); meta.index=meta.net+'.'+meta.sta
    with np.load(args.travel_table) as table:
        travel={p:RegularGridInterpolator((table['x'],table['z']),table[p],bounds_error=True) for p in ['P','S']}
    picks=pd.read_parquet(neural/'picks.parquet')
    pairs=pd.read_csv(shared/'matched_pairs.csv')
    for p in [Path(__file__), mf/'manifest.json',mf/'templates.npz',mf/'training.json',
              mf/'calibration_eligible_reference.parquet',neural/'manifest.json',neural/'picks.parquet',
              clocks/'family_clocks.json',clocks/'station_transfers.csv',shared/'manifest.json',
              shared/'matched_pairs.csv',assoc/'manifest.json',args.station_metadata,args.travel_table]: record(p)
    rows=[]
    for split in ['calibration','test']:
        path=shared/f'{split}_reference.parquet'; record(path); refs=pd.read_parquet(path)
        for date, day in refs.groupby('date'):
            for station, station_cfg in cfg['stations'].items():
                path=neural/f'{date}_{station}_probabilities.npz';record(path)
                with np.load(path) as archive: scores=archive['probabilities']
                if scores.shape != (2,int(86400*cfg['fs'])): raise ValueError('unexpected probability dimensions')
                m=meta.loc[station]
                station_picks=picks[(picks.station==station)&(picks.date==date)]
                for phase_index,phase in enumerate(['P','S']):
                    pp=station_picks[station_picks.phase==phase].sort_values('time_s')
                    pt,ps=pp.time_s.to_numpy(),pp.model_score.to_numpy()
                    for ref in day.itertuples():
                        distance=np.hypot((ref.latitude-m.lat)*111.195,
                            (ref.longitude-m.lon)*111.195*np.cos(np.deg2rad((ref.latitude+m.lat)/2)))
                        modeled=ref.time_s+float(travel[phase]([[distance,ref.depth_km]])[0])
                        centers={'model':modeled}
                        key=(ref.family_id,station,phase)
                        if key in offsets:
                            centers['transferred_stack']=ref.time_s-accepted[ref.family_id]+offsets[key]
                        for basis,center in centers.items():
                            for shift in [0.,-300.,-120.,120.,300.]:
                                for radius in [1.2,2.,5.]:
                                    covered,maximum,peak,residual=phase_window(scores[phase_index],date_start(date),
                                        cfg['fs'],center+shift,radius,pt,ps)
                                    rows.append({'split':split,'date':date,'reference_id':ref.reference_id,
                                        'family_id':ref.family_id,'station':station,'site':station_cfg['site_id'],
                                        'phase':phase,'arrival_basis':basis,'shift_s':shift,'radius_s':radius,
                                        'expected_arrival_s':center+shift,'covered':covered,
                                        'max_probability':maximum,'max_pick_score':peak,'pick_residual_s':residual})
            print(f'[neural diagnosis] {split} {date}',flush=True)
    frame=pd.DataFrame(rows);frame.to_parquet(out/'phase_windows.parquet',index=False)
    summaries=[]; events=[]
    phase_keys=['split','station','phase','arrival_basis','shift_s','radius_s']
    for key,g in frame.groupby(phase_keys):
        for threshold in [.1,.3,.5]:
            covered=g[g.covered]
            summaries.append(dict(zip(phase_keys,key))|{'threshold':threshold,'n_reference_windows':len(g),
                'n_covered':len(covered),'n_with_probability_crossing':int((covered.max_probability>=threshold).sum()),
                'n_with_exported_peak':int((covered.max_pick_score>=threshold).sum())})
    pd.DataFrame(summaries).to_csv(out/'phase_summary.csv',index=False)
    keys=['split','reference_id','family_id','shift_s','radius_s']
    for key,g in frame[frame.arrival_basis=='model'].groupby(keys):
        for threshold in [.1,.3,.5]:
            row=dict(zip(keys,key))|{'threshold':threshold}|support_counts(g,threshold)
            for minimum in [4,5]:
                hit=pairs[(pairs.arm=='neural_association')&(pairs.setting==f'p{threshold:g}_n{minimum}')&
                          (pairs.tolerance_s==2)&(pairs.reference_id==row['reference_id'])]
                row[f'associated_n{minimum}_within_2s']=bool(len(hit)) if threshold in [.3,.5] and row['shift_s']==0 else None
            events.append(row)
    events=pd.DataFrame(events);events.to_parquet(out/'event_support.parquet',index=False)
    summary=[]
    for key,g in events.groupby(['split','shift_s','radius_s','threshold']):
        covered=g[g.all_s_covered]
        summary.append(dict(zip(['split','shift_s','radius_s','threshold'],key))|{
            'n_references':len(g),'n_all_s_covered':len(covered),
            'n_four_s_sites':int((covered.n_s_sites>=4).sum()),
            'n_four_s_sites_five_s_stations':int(((covered.n_s_sites>=4)&(covered.n_s_stations>=5)).sum()),
            'n_four_crossing_sites':int((covered.n_s_crossing_sites>=4).sum())})
    pd.DataFrame(summary).to_csv(out/'event_support_summary.csv',index=False)
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        completed=list(pool.map(matcher_date,[(cfg,str(mf),d,9.) for d in cfg['calibration_dates']]))
    mf_rows=[r for rows,_ in completed for r in rows]
    pd.DataFrame(mf_rows).to_parquet(out/'matcher_calibration_reference_scores.parquet',index=False)
    pd.DataFrame(mf_rows).groupby('status').size().rename('n_references').to_csv(out/'matcher_loss_summary.csv')
    write_json(out/'manifest.json',{'sources':sources,'matcher_waveform_inputs':[x for _,x in completed],
        'elapsed_s':time.monotonic()-started,'workers':args.workers,
        'status':'post_fit_development_miss_diagnosis_no_retuning',
        'limitations':['Modeled/transferred stack arrivals are not adjudicated single-event phases',
          'Reference-centered window evidence is not one-to-one pick recall; burst references may share peaks',
          'Shifted windows measure local background coincidence, not false-discovery rate or known noise',
          'Fixed-location S-support is necessary evidence at the assumed time/model, not a rerun of association',
          'Threshold 0.1 is diagnostic only; no new association or detector setting selected',
          'Matcher reference maxima precede peak suppression and family resolution; not detection counts',
          'Low model scores do not establish that a visible physical phase exists or is absent']})
    print(pd.DataFrame(summary).query('shift_s==0 and radius_s==2').to_string(index=False),flush=True)
    print(pd.DataFrame(mf_rows).groupby('status').size().to_string(),flush=True)


if __name__=='__main__': main()
