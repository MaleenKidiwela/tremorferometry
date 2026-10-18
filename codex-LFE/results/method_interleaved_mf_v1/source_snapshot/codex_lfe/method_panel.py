"""Collect same-date development diagnostics without pretending clocks are resolved."""
import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from .provenance import sha256, write_json
from .reference import read_bostock


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True)
    ap.add_argument('--matcher', required=True)
    ap.add_argument('--association', required=True)
    ap.add_argument('--clocks', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path(args.config).read_text())
    mf, neural, clocks = map(Path, [args.matcher, args.association, args.clocks])
    report = json.loads((mf/'report.json').read_text())
    reference = read_bostock(cfg['reference_csv'])
    training = json.loads((mf/'training.json').read_text())
    cohort = []
    for split in ['train','calibration','test']:
        g = reference[reference.date.isin(cfg[f'{split}_dates'])]
        cohort.append({'split':split,'all_reference_events':len(g),
                       'selected_group_events':int(g.family_id.isin(training['families_selected']).sum()),
                       'usable_template_group_events':int(g.family_id.isin(training['families_usable']).sum()),
                       'all_source_groups':g.family_id.nunique()})
    pd.DataFrame(cohort).to_csv(out/'source_cohort.csv',index=False)
    rows = []
    for split, settings in report['splits'].items():
        for threshold, value in settings.items():
            r,n = value['reference'],value['null']
            rows.append({'split':split,'threshold_mad':float(threshold),
                         'candidates':r['n_detections'],'reference_eligible':r['n_reference'],
                         'matched_same_source':r['n_matched_same_family'],
                         'reference_recovery':r['same_family_recall'],
                         'reference_agreement_not_precision':r['reference_agreement_fraction'],
                         'control_counts':json.dumps(n['control_counts']) if n else None,
                         'mean_control_ratio':n['control_to_observed_mean'] if n else None,
                         'max_control_ratio':n['control_to_observed_range'][1] if n else None})
    pd.DataFrame(rows).to_csv(out/'matched_filter_diagnostics.csv',index=False)
    runs = pd.read_csv(neural/'runs.csv')
    raw_rows = []
    for row in runs.itertuples():
        events = pd.read_parquet(neural/f'{row.name}_events.parquet')
        dates = (pd.to_datetime(events.time,unit='s',utc=True).dt.strftime('%Y-%m-%d')
                 if len(events) else pd.Series(dtype=str))
        for split in ['calibration','test']:
            subset = events[dates.isin(cfg[f'{split}_dates'])] if len(events) else events
            raw_rows.append({'split':split,'threshold':row.threshold,
                             'minimum_stations':row.minimum_stations,'mode':row.mode,'seed':row.seed,
                             'candidates':len(subset),
                             'near_depth_boundary':int(subset.near_depth_boundary.sum())
                             if len(subset) else 0})
    counts = pd.DataFrame(raw_rows)
    counts.to_csv(out/'association_counts_by_run.csv',index=False)
    nn = []
    for key,g in counts.groupby(['split','threshold','minimum_stations']):
        observed = g[g['mode']=='observed'].iloc[0]
        for mode in ['station','station_hour']:
            control = g[g['mode']==mode].candidates
            nn.append({'split':key[0],'threshold':key[1],'minimum_stations':int(key[2]),
                       'mode':mode,'candidates':int(observed.candidates),
                       'near_depth_boundary':int(observed.near_depth_boundary),
                       'control_counts':json.dumps(control.tolist()),
                       'mean_control_ratio':float(control.mean()/observed.candidates)
                       if observed.candidates else None})
    pd.DataFrame(nn).to_csv(out/'association_diagnostics.csv',index=False)
    family_clocks = json.loads((clocks/'family_clocks.json').read_text())
    files = [Path(args.config),mf/'report.json',mf/'training.json',neural/'manifest.json',
             clocks/'manifest.json',clocks/'family_clocks.json',Path(cfg['reference_csv']),Path(__file__)]
    files += sorted(neural.glob('*_events.parquet'))
    write_json(out/'report.json',{
        'status':'same_date_development_diagnostics_not_a_method_ranking',
        'sources':{str(p.resolve()):sha256(p) for p in files},
        'matcher_threshold_selection':report['selection'],
        'source_cohort':cohort,
        'accepted_multisite_family_clocks':sum(c['accepted'] for c in family_clocks),
        'neural_reference_recovery_computed':False,
        'why_no_winner':['Reference origin clocks did not pass training-only multisite checks',
                         'Different source inventories: known-family matching versus open association',
                         'Common station/date panel but no common valid-exposure mask yet',
                         'Different control shift lengths and realizations',
                         'No independently adjudicated event labels',
                         'Dates are development data already examined by earlier experiments'],
        'next_evidence':'Reviewed training-stack phase arrivals and candidate waveforms; reconcile clocks and exposure before ranking.'})
    print('Matched-filter diagnostics:')
    print(pd.DataFrame(rows).drop(columns='control_counts').to_string(index=False))
    print('Association on the same dates (different clocks/exposure/source inventory):')
    print(pd.DataFrame(nn).to_string(index=False))


if __name__ == '__main__':
    main()
